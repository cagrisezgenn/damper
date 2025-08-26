%% =====================================================================
%  10-Katlı Çerçeve — ODE-only viskoz damper modeli (CLEAN / AHENK + Switch)
%  (Sıkışabilirlik + Hat ataletı + Cd(Re) orifis + kavitasyon + 2-düğümlü ısı)
%  Laminer kayıp hidrolikte: Δp_lam(T) = R_lam(T)*Q
%  Kavitasyon: p2_eff = max(p2, cav_sf * p_vap(T)) hem akışta hem kuvvette
%  Solver: ode23tb (+ güvenli deval + ode15s fallback)
%  NOT: Antoine katsayılarını yağınıza göre kalibre edin.
% =====================================================================

clear; clc; close all;
% Global log anahtarları
global LOG; LOG = struct('verbose_decode', false);
global PARETO;
if isempty(whos('global','PARETO')) || isempty(PARETO)
    PARETO = struct('J1',[],'J2',[],'F',[],'Pen',[], ...
                    'set',[],'x',{{}},'feas',[]);
end

%% -------------------- Kullanıcı anahtarları ---------------------------
% (1) Kayıt/yön seçimi + görselleştirme
vis.make_plots          = true;
vis.sel_rec             = 1;         % 1..R
vis.sel_dir             = 'x';       % 'X' veya 'Y'
vis.dual_run_if_needed  = true;      % plot_source ≠ sim_source ise ikinci koşu yap

% (2) Örnekleme / yeniden örnekleme anahtarları
prep.target_dt      = 0.005;     % hedef dt (s)
prep.resample_mode  = 'auto';     % 'auto'|'off'|'force'
prep.regrid_method  = 'pchip';   % 'pchip'|'linear'
prep.tol_rel        = 1e-6;      % dt eşitlik toleransı (göreli)

% (3) Şiddet eşitleme / PSA / CMS anahtarları
pp.on.intensity      = true;       % band-ortalama Sa(+CMS) normalizasyonu (yalnız GA için)
pp.on.CMS            = false;      % cms_target.mat varsa ve kullanmak istersen true
pp.gammaCMS          = 0.50;       % hibrit ağırlık (0→yalnız band, 1→yalnız CMS)

pp.PSA.zeta          = 0.05;       % SDOF sönüm oranı
pp.PSA.band_fac      = [0.8 1.2];  % T1 bandı (T1±%20)
pp.PSA.Np_band       = 15;         % band içi periyot sayısı
pp.PSA.use_augmented = true;       % tek ODE ile çok periyot (hızlı)
pp.PSA.downsample_dt = 0.02;       % SA hesabı için isteğe bağlı downsample (<=0 kapalı)
pp.PSA.use_parfor    = false;      % Parallel Toolbox varsa denersin

% (4) Arias penceresi & kuyruk
pp.on.arias          = true;
pp.tail_sec          = 20;

% (5) Simülasyon ve grafik için veri kaynağı seçimi
%     'raw'    : ham ivmeler (genlik korunur)
%     'scaled' : şiddet eşitlenmiş ivmeler (GA/normalize amaçlı)
sel.sim_source  = 'scaled';
sel.plot_source = 'raw';

%% -------------------- Amaç fonksiyonu anahtarları ---------------------
% (A) Hedef katlar ve metrik
obj.idx_disp_story   = 10;        % d_10 → tepe yer değiştirme
obj.idx_acc_story    = 3;         % a_3  → ivme metriği
obj.acc_metric       = 'rms+p95';     % 'rms' | 'energy' | 'rms+p95' (hibrit)
obj.p95_penalty_w    = 0.20;      % hibritte pik cezası ağırlığı
%% -------------------- Kısıt anahtarları -------------------------------
% Kısıtları aç/kapa ve eşik/ceza ayarları
cons.on.spring_tau     = true;    % K1: yay kesme gerilmesi (τ_max ≤ τ_allow)
cons.on.spring_slender = true;    % K2: L_free/D_m ≤ λ_max
cons.on.stroke         = true;    % K3: max|drift| ≤ 0.9*L_gap
cons.on.force_cap      = true;   % K4: max|F_story| ≤ F_cap  (varsayılan: kapalı; F_cap tanımla)
cons.on.dp_quant       = true;    % K5: q≈0.99 Δp_orf ≤ dP_cap
cons.on.thermal_dT     = true;    % K6: ΔT_est ≤ ΔT_cap
cons.on.cav_frac       = false;    % K7: kavitasyon payı sınırı
cons.on.qsat_margin    = true;    % K8: Q 95p ≤ margin*Qcap_big
cons.on.fail_bigM      = true;    % K9: Simülasyon başarısızsa büyük ceza

% Eşikler / limitler
cons.spring.tau_allow   = 300e6;  % [Pa] yay çeliği için tipik, gerekirse güncelle
cons.spring.lambda_max  = 12.0;   % boy/çap sınırı
cons.spring.L_free_mode = 'auto'; % 'auto' veya 'fixed'
cons.spring.L_free_fix  = NaN;    % 'fixed' için metre cinsinden serbest boy
cons.spring.L_free_auto_fac = 2.2; % 'auto' modda L_free ≈ fac * L_gap

cons.stroke.util_factor = 0.90;   % izinli strok = 0.90*L_gap

cons.force.F_cap        = inf;    % [N] cihaz kuvvet sınırı; Inf → devre dışı
cons.dp.q               = 0.99;   % Δp_orf zaman-içi quantile
cons.dp.agg             = 'max';  % kayıtlar arası: 'max' | 'cvar'
cons.alpha_CVaR_cons    = 0.20;   % yalnız 'cvar' seçilirse kullanılır

cons.thermal.cap_C      = 30.0;   % [°C] yağ ısınma sınırı
cons.hyd.cav_frac_cap   = 0.15;   % kavitasyon zaman oranı sınırı (95p)
cons.hyd.Q_margin       = 0.90;   % Q 95p ≤ margin*Qcap_big

% Ceza ayarları
cons.pen.power  = 2;              % hinge^power
cons.pen.bigM   = 1e6;            % simülasyon başarısızlığında eklenecek ceza
cons.pen.lambda = struct( ...     % kısıt başına ağırlıklar
    'spring_tau',     0.3, ...%0.3
    'spring_slender', 0.2, ...%0.2
    'stroke',         1.2, ...%1.2
    'force_cap',      0, ...%0
    'dp_quant',       0.5, ...%0.5
    'thermal_dT',     0.2, ...%0.2
    'cav_frac',       1.5, ...%1.5
    'qsat_margin',    0.3, ...%0.3
    'fail_bigM',      1 );%

% Kısıt simülasyon kaynağı ve μ seti (amaçla tutarlı kalsın)
cons.src_for_constraints = 'raw';           % kısıtları 'raw' ile değerlendir
if isfield(obj,'mu_scenarios')
    cons.mu_scenarios = obj.mu_scenarios;
else
    cons.mu_scenarios = [0.75 1.00 1.25];  % varsayılan; aşağıda obj tanımlanınca senkronlayacağız
end

% (B) Arias penceresi içi ölçüm
obj.use_arias_window = true;      % true → [t5,t95] aralığında metrik
obj.window_source    = 'same';    % 'same' → kullanılan seri üzerinden hesapla

% (C) Yön zarfı
obj.dir_mode         = 'envelope'; % 'envelope' (=max(X,Y)) | 'Xonly' | 'Yonly'

% (D) μ-robustluk
obj.mu_scenarios = [0.9 1.0 1.1];
obj.mu_aggregate = 'weighted';
obj.mu_weights   = [0.25 0.50 0.25];
% Kısıt μ-senaryolarını amaçla hizala
cons.mu_scenarios = obj.mu_scenarios;

% (E) Risk agregasyonu (kayıtlar üzerinde)
obj.use_scaled_for_goal = true;   % amaç için 'scaled' setiyle çalış (önerilir)
obj.alpha_CVaR       = 0.20;      % tail payı
obj.weights_da       = [0.5 0.5]; % [w_disp, w_acc] toplam 1.0

% (F) Referans tanımı
% d_ref ve a_ref: aynı kayıt/yön için (damper yok, μ=1.0), sabit baz koşusundan.
% NOT: Bu referanslar bir kez hesaplanır ve tasarımla değişmez.
%% -------------------- Model anahtarları (blok-bazlı aç/kapa) ----------
% NOT: Bu blok dosyanın BAŞINDA olmalı. Aşağıdaki varsayılanları
% dilediğin gibi değiştir; ayrıca ensure_cfg_defaults() eksikleri tamamlar.
cfg.use_orifice = true;
cfg.use_thermal = true;
cfg.on.CdRe            = true;
cfg.on.Rlam            = true;
cfg.on.Rkv             = true;
cfg.on.Qsat            = true;
cfg.on.cavitation      = true;
cfg.on.dP_cap          = true;
cfg.on.hyd_inertia     = true;
cfg.on.leak            = true;
cfg.on.pressure_ode    = true;
cfg.on.pressure_force  = true;
cfg.on.mu_floor        = true;


% Basınç-kuvvet rampa/kazanç (PF)
cfg.PF.mode  = 'ramp';
cfg.PF.t_on  = nan;      % t5 sonrası atanacak  <-- BURASI NaN KALSIN
cfg.PF.tau   = 0.9;      % 4.0 → 1.5  (önerim)
cfg.PF.gain  = 1.7;      % 0.45 → 1.0 (önerim)
cfg.on.pf_resistive_only = true;

% Eksik/yanlış alanlar için güvenli tamamlayıcı (guard)
cfg = ensure_cfg_defaults(cfg);
%% -------------------- GA/Opt tasarım anahtarları ----------------------
ga.enable      = false;    % GA tasarım vektörü uygula? (false → baz set)
ga.design_set  = 1;        % 1|2|3 (aşağıdaki set tanımları)
ga.x           = [];       % Örn: set-1 için 9x1 vektör, boşsa uygulanmaz
% Örnek: otomatik orta-nokta denemesi
% [lb,ub,~,~] = ga_get_bounds(ga.design_set); ga.x = 0.5*(lb+ub); ga.enable=true;
%% -------------------- GA ÇALIŞTIRICI (opsiyonel) ----------------------
% --- Değişken sınırları ve IntCon hazırlığı ---
[lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
if isempty(int_idx)
    IntCon = [];           % tamsayı değişken yok
else
    IntCon = int_idx(:)';  % tamsayı indeks vektörü
end

% GA anahtarları
ga.opt.enable       = false;
ga.opt.seed         = 11;
ga.opt.popsize      = 12;
% Çok aşamalı akış: önce set-1, sonra set-3
ga.opt.multi_stage  = [1 1 2 2 3 3];
% Aşama başına nesil sayısı (opsiyonel):
ga.opt.maxgen_stage = [20 12 20 12 16 12];  % her blok için nesil
ga.opt.maxgen       = ga.opt.maxgen_stage(1);   % ilk aşama için başlangıç
ga.opt.use_parallel  = true;      % GA.UseParallel
ga.opt.lhs_refine    = true;      % 2-aşamalı LHS (coarse→refine)
ga.opt.cache_enable  = true;      % aynı x için memoize
ga.opt.fail_early_k  = 6;         % erken çıkış eşiği (başarısız sim sayısı)
ga.opt.save_log      = 'runLog_ga.mat';  % nüfus & en iyi çözüm logu
% daraltma parametreleri
ga.opt.keep_top     = 0.20;   % en iyi %20'yi tut
ga.opt.buffer_fac   = 0.10;   % p10–p90 etrafına %10 tampon

%% -------------------- Girdiler (7 kayıt; sütunlar: t, ax, (ops) ay) ---
Sall = load('acc_matrix.mat');
fn   = fieldnames(Sall);
fn   = fn(startsWith(fn,'acc_matrix'));
R    = numel(fn);
if R==0, error('acc_matrix.mat içinde acc_matrix* isimli dizi bulunamadı.'); end

%% -------------------- Yapı (T1 için gerekli) --------------------------
n  = 10;
m  = 2.2e6 * ones(n,1);
k  = 2.95e8 * ones(n,1);
c0 = 2.55e6 * ones(n,1);
[M,K,Cstr] = make_KCM(n,m,k,c0);
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);
% --- IDR için kat yükseklik(leri) ---
h_story_m = 3.0 * ones(n-1,1);   % tüm katlar 3.0 m ise
% h_story_m = [h1; h2; ...; h_{n-1}];   % kat kat farklı ise (alternatif)


%% -------------------- PSA fonksiyon seçimi ----------------------------
if pp.PSA.use_augmented
    f_band = @sdof_PSA_band_avg_aug;   % tek ODE, çok periyot
    f_vec  = @sdof_PSA_vec_aug_ode;
else
    f_band = @sdof_PSA_band_avg_ode;   % periyot başına ayrı ODE
    f_vec  = @sdof_PSA_vec_ode;
end

%% -------------------- Kayıtları oku → RAW (tekil & hedef dt kontrolü) -
t_rawX = cell(R,1); t_rawY = cell(R,1);
a_rawX = cell(R,1); a_rawY = cell(R,1);

for r=1:R
    A = Sall.(fn{r});
    if size(A,2)<2, error('%s: en az iki sütun (t, ax) olmalı.', fn{r}); end
    t0 = A(:,1); ax0 = A(:,2); ay0 = []; if size(A,2)>=3, ay0 = A(:,3); end
    [t0,iu] = unique(t0,'stable'); ax0 = ax0(iu); if ~isempty(ay0), ay0=ay0(iu); end

    [tX,ax] = regrid_to_target(t0, ax0, prep);
    if ~isempty(ay0), [~,ay] = regrid_to_target(t0, ay0, prep); else, ay=[]; end

    t_rawX{r} = tX; a_rawX{r} = ax;
    t_rawY{r} = tX; a_rawY{r} = ay;   % Y varsa aynı ızgara

    % Hedef dt kontrol uyarısı (resample_mode='off' iken dahi)
    dtX = median(diff(tX),'omitnan');
    tol = max(prep.tol_rel*max(prep.target_dt,eps), 1e-12);
    if abs(dtX - prep.target_dt) > tol
        warning('Kayıt #%d dt=%.6g s, hedef=%.6g s (resample kapalı).', r, dtX, prep.target_dt);
    end
end

%% -------------------- Şiddet eşitleme (scaled set) --------------------
t_sclX = t_rawX; t_sclY = t_rawY;
a_sclX = a_rawX; a_sclY = a_rawY;   % default: scaled = raw

if pp.on.intensity
    zeta_SA = pp.PSA.zeta; band_fac = pp.PSA.band_fac; Np_band = pp.PSA.Np_band;

    % Opsiyonel CMS hedefi
    useCMS = false; T_cms = []; Sa_cms = [];
    if pp.on.CMS && exist('cms_target.mat','file')
        Scms = load('cms_target.mat');
        if isfield(Scms,'T_cms') && isfield(Scms,'Sa_cms') && numel(Scms.T_cms)==numel(Scms.Sa_cms)
            T_cms  = Scms.T_cms(:); Sa_cms = Scms.Sa_cms(:); useCMS = true;
        end
    end

    % Band Sa (parfor opsiyonel)
    Sa_band = zeros(R,1);
    canPar  = pp.PSA.use_parfor && ~isempty(ver('parallel'));
    if canPar
        parfor r=1:R
            [tPSA,agPSA] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
            Sa_band(r) = f_band(tPSA, agPSA, T1, zeta_SA, band_fac, Np_band);
        end
    else
        for r=1:R
            [tPSA,agPSA] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
            Sa_band(r) = f_band(tPSA, agPSA, T1, zeta_SA, band_fac, Np_band);
        end
    end
    Sa_band_target = median(Sa_band);

    % Her kayıt için ölçek uygula (yalnız GA amaçlı; solver varsayılanı RAW)
    for r=1:R
        [tPSA,agPSA] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
        Sab_r  = f_band(tPSA, agPSA, T1, zeta_SA, band_fac, Np_band);
        s_band = Sa_band_target / max(Sab_r,eps);

        if useCMS
            Sa_rec = f_vec(tPSA, agPSA, T_cms, zeta_SA);
            nume   = max(sum(Sa_rec.*Sa_cms),eps);
            deno   = max(sum(Sa_rec.*Sa_rec),eps);
            s_cms  = nume/deno;
            s_hyb  = s_band^(1-pp.gammaCMS) * s_cms^(pp.gammaCMS);
        else
            s_hyb  = s_band;
        end

        a_sclX{r} = s_hyb * a_rawX{r};
        if ~isempty(a_rawY{r}), a_sclY{r} = s_hyb * a_rawY{r}; end
    end
end

%% -------------------- Arias pencereleri (raw & scaled) -----------------
[t5x_raw,t95x_raw,t5y_raw,t95y_raw] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));
[t5x_scl,t95x_scl,t5y_scl,t95y_scl] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));

for r=1:R
    if pp.on.arias
        [t5x_raw(r), t95x_raw(r)] = arias_win(t_rawX{r}, a_rawX{r}, 0.05, 0.95);
        if ~isempty(a_rawY{r}), [t5y_raw(r), t95y_raw(r)] = arias_win(t_rawY{r}, a_rawY{r}, 0.05, 0.95); end

        [t5x_scl(r), t95x_scl(r)] = arias_win(t_sclX{r}, a_sclX{r}, 0.05, 0.95);
        if ~isempty(a_sclY{r}), [t5y_scl(r), t95y_scl(r)] = arias_win(t_sclY{r}, a_sclY{r}, 0.05, 0.95); end
    else
        t5x_raw(r)=t_rawX{r}(1);  t95x_raw(r)=t_rawX{r}(end);
        t5x_scl(r)=t_sclX{r}(1);  t95x_scl(r)=t_sclX{r}(end);
        if ~isempty(a_rawY{r}), t5y_raw(r)=t5x_raw(r);  t95y_raw(r)=t95x_raw(r); end
        if ~isempty(a_sclY{r}), t5y_scl(r)=t5x_scl(r);  t95y_scl(r)=t95x_scl(r); end
    end
end
%% -------------------- Baz parametre setleri (parametrebulur uyumlu) ---
% Damper geometri + malzeme (BAZ)
geom.Dp    = 0.12;                 % piston çapı [m]
geom.Lgap  = 0.20;                 % etkin strok boşluğu [m]
geom.d_o   = 0.0022;               % orifis çapı [m]
geom.Lori  = 0.03;                 % orifis uzunluğu [m]
geom.Kd    = 1.8e9;                % yağın k_b (sıkışabilirlikten gelen) için ölçek [Pa]
geom.Ebody = 2.1e11;               % gövde elastisite modülü [Pa]
geom.Ap    = pi*geom.Dp^2/4;       % piston alanı [m^2] (türetilen)

% Yay (spiral) malzeme/geo (BAZ)
sh.G      = 79e9;                  % kayma modülü [Pa]
sh.d_w    = 0.018;                 % tel çapı [m]
sh.D_m    = 0.08;                  % yay ortalama çapı [m]
sh.n_turn = 26;                    % sarım sayısı [-]

% Orifis / akışkan (BAZ)  — parametrebulur ile aynı alan adları
if ~exist('orf','var')   || ~isstruct(orf),   orf   = struct(); end
if ~exist('hyd','var')   || ~isstruct(hyd),   hyd   = struct(); end
if ~exist('therm','var') || ~isstruct(therm), therm = struct(); end
if ~exist('cfg','var')   || ~isstruct(cfg),   cfg   = struct(); end
if ~exist('num','var')   || ~isstruct(num),   num   = struct(); end

% Orifis defaultları
orf.n_orf  = getfield_default(orf, 'n_orf',  2);
orf.Cd0    = getfield_default(orf, 'Cd0',    0.6);
orf.CdInf  = getfield_default(orf, 'CdInf',  0.9);
orf.Rec    = getfield_default(orf, 'Rec',    3800);
orf.p_exp  = getfield_default(orf, 'p_exp',  1.1);
orf.p_amb  = getfield_default(orf, 'p_amb',  1.0e5);
orf.cav_sf = getfield_default(orf, 'cav_sf', 1.05);

% Termal + T-bağımlı (BAZ) — parametrebulur alanları
therm.antoine_A = getfield_default(therm,'antoine_A',5.0);
therm.antoine_B = getfield_default(therm,'antoine_B',1700);
therm.antoine_C = getfield_default(therm,'antoine_C',-80);

therm.T0_C      = getfield_default(therm,'T0_C',25);
therm.Ts0_C     = getfield_default(therm,'Ts0_C',25);
therm.T_env_C   = getfield_default(therm,'T_env_C',25);

therm.hA_o_env  = getfield_default(therm,'hA_o_env',800);
therm.hA_s_env  = getfield_default(therm,'hA_s_env',600);
therm.hA_os     = getfield_default(therm,'hA_os',   800);

therm.resFactor = getfield_default(therm,'resFactor',22);  % emniyet alt sınırı
therm.cp_oil    = getfield_default(therm,'cp_oil',   1800);
therm.cp_steel  = getfield_default(therm,'cp_steel', 500);

therm.rho_ref   = getfield_default(therm,'rho_ref',  850);
therm.T_ref_C   = getfield_default(therm,'T_ref_C',   25);
therm.alpha_rho = getfield_default(therm,'alpha_rho',7e-4);
therm.beta0     = getfield_default(therm,'beta0',   1.6e9);
therm.b_beta    = getfield_default(therm,'b_beta', -0.0035);
therm.mu_ref    = getfield_default(therm,'mu_ref',   1.2);
therm.b_mu      = getfield_default(therm,'b_mu',   -0.011);

% Hidrolik / sayısı
hyd.n_parallel = getfield_default(hyd,'n_parallel', 2);
hyd.K_leak     = getfield_default(hyd,'K_leak',     0);
hyd.Lh         = getfield_default(hyd,'Lh',     0.0009);  % emniyet
hyd.Vmin_fac   = getfield_default(hyd,'Vmin_fac',0.98); % emniyet

% Numerik guardlar (parametrebulur ile aynı isimler)
num.dP_cap       = getfield_default(num,'dP_cap',       5e7);
num.mu_min_phys  = getfield_default(num,'mu_min_phys',  0.25);
num.softmin_eps  = getfield_default(num,'softmin_eps',  1e3);



%% -------------------- Türetilen parametreler --------------------------
% Alanlar
geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);

% Paralel damper etkileri
nd               = max(1, getfield_default(hyd,'n_parallel',1));
geom.Ap_eff      = nd * geom.Ap;
orf.Ao_eff       = nd * orf.Ao;
hyd.n_parallel   = nd;   % garanti

% Qsat/doygunluk için referans debi (ilk kodla uyumlu)
cd_ref        = max(orf.CdInf, orf.Cd0);
dp_for_qcap   = getfield_default(num,'dP_cap', 3e8);     % num.dP_cap varsa onu kullan
Ae_ref        = max(cd_ref * orf.Ao_eff, 1e-12);
num.Qcap_big  = getfield_default(num,'Qcap_big', ...
                    hyd.Vmin_fac * Ae_ref * sqrt( 2*dp_for_qcap / max(therm.rho_ref,100) ) );


% Yay rijitlikleri
k_p   = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);   % coil spring
k_h   = geom.Kd*geom.Ap^2/geom.Lgap;            % hidrolik (tek damper)
k_s   = geom.Ebody*geom.Ap/geom.Lgap;           % gövde (tek damper)
k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));    % seri birleşim
k_sd  = nd * (k_hyd + k_p);                     % kat başına efektif yay

% Yağ hacmi ve ısı kapasiteleri (parametrebulur → eq_demo_* ile uyumlu)
nStories         = n - 1;
steel_to_oil_mass_ratio = 1.5;

hyd.V0           = 0.5 * (geom.Ap * (2*geom.Lgap));        % tek hazne tahmini
V_oil_per        = therm.resFactor * (geom.Ap * (2*geom.Lgap));  % tek damper yağ [m^3]
nDtot            = nStories * nd;

m_oil_tot        = nDtot * (therm.rho_ref * V_oil_per);
m_steel_tot      = steel_to_oil_mass_ratio * m_oil_tot;
therm.C_oil      = m_oil_tot   * therm.cp_oil;    % [J/K]
therm.C_steel    = m_steel_tot * therm.cp_steel;  % [J/K]

% Bilgi: referans laminer direnç (T_ref’te)
R_lam0 = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / nd;


%% -------------------- Seçim yardımcıları ------------------------------
pickA = @(SRC,Rid,DIR) pick_series(SRC,Rid,DIR, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl);

%% ==================== GA Koşusu (isteğe bağlı) ========================
% Çok aşamalı GA: önce set-1 (geometri/rijitlik), sonra set-3 (akışkan/termal)
% Kullanım: ga.opt.multi_stage = [1 3]; // yoksa mevcut tek aşamalı akış sürer
if isfield(ga,'opt') && isfield(ga.opt,'enable') && ga.opt.enable

    % --- Aşama listesi (yoksa tek aşama: mevcut design_set) ---
    if isfield(ga.opt,'multi_stage') && ~isempty(ga.opt.multi_stage)
        stage_list = ga.opt.multi_stage(:).';
    else
        stage_list = ga.design_set;
    end

    last_best = []; last_set = NaN;

    for si = 1:numel(stage_list)
    ga.design_set = stage_list(si);

    % --- Değişken sınırları + override ---
    [lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
    if isfield(ga.opt,'lb_override') && ~isempty(ga.opt.lb_override)
        lb = max(lb, ga.opt.lb_override);
    end

    if isfield(ga.opt,'ub_override') && ~isempty(ga.opt.ub_override)
        ub = min(ub, ga.opt.ub_override);
    end
    if isempty(int_idx), IntCon = [];
    else, IntCon = int_idx(:)'; end

    % --- Bilgi satırı (banner) ---
    switch ga.design_set
      case 1
        fprintf('[GA] set=1 ... bounds(d_o)=%.1f..%.1f mm | Dp=%.0f..%.0f mm | Lgap=%.0f..%.0f mm\n', ...
    ga.opt.popsize, ga.opt.maxgen, 1e3*lb(1),1e3*ub(1), 1e3*lb(8),1e3*ub(8), 1e3*lb(9),1e3*ub(9));

      case 2
        fprintf(['[GA] set=2 | pop=%d | gen=%d | bounds(n_orf)=%g..%g | Cd0=%.2f..%.2f | CdInf=%.2f..%.2f | ' ...
                 'Rec=%g..%g | p_exp=%.1f..%.1f | cav_{sf}=%.2f..%.2f | Lh=%.1f..%.1f mm | K_{leak}=%.e..%.e | resFactor=%g..%g\n'], ...
            ga.opt.popsize, ga.opt.maxgen, ...
            lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), ...
            1e3*lb(7),1e3*ub(7), lb(8),ub(8), lb(9),ub(9));
      case 3
        fprintf(['[GA] set=3 | pop=%d | gen=%d | bounds(n_orf)=%g..%g | Cd0=%.2f..%.2f | ' ...
                 'mu_{ref}=%.2g..%.2g Pa·s | b_{\\mu}=%.3g..%.3g 1/°C | \\beta_0=%.1e..%.1e Pa | ' ...
                 'b_{\\beta}=%.3g..%.3g 1/°C | hA_{os}=%g..%g W/K | dP_{cap}=%.1e..%.1e Pa | ' ...
                 'Vmin_{fac}=%.2f..%.2f | resFactor=%g..%g\n'], ...
            ga.opt.popsize, ga.opt.maxgen, ...
            lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), ...
            lb(7),ub(7), lb(8),ub(8), lb(9),ub(9), lb(10),ub(10));
    end

    % --- Başlangıç nüfusu (LHS) ---
    P0 = lhs_population(lb,ub,ga.opt.popsize);
    if isfield(ga.opt,'lhs_refine') && ga.opt.lhs_refine
        P0(1,:) = 0.5*(lb+ub);
    end

    % --- Aşama-özel kısıt/PF ayarı ---
    cons_stage = cons;
    cfg_stage  = cfg;
   switch ga.design_set
    case 1  % geometri/rijitlik
        cons_stage.on.dp_quant    = true;
        cons_stage.on.thermal_dT  = true;
        cons_stage.on.cav_frac    = false;
        cons_stage.on.qsat_margin = true;

        cons_stage.pen.lambda.spring_tau     = 0.8;
        cons_stage.pen.lambda.spring_slender = 0.4;

        cfg_stage.PF.gain = 1.7;
        cfg_stage.PF.tau  = 0.9;

    case 2  % hidrolik (kavitasyonu önce bastır)
        cons_stage.on.dp_quant    = true;
        cons_stage.on.cav_frac    = true;
        cons_stage.on.qsat_margin = true;
        cons_stage.on.thermal_dT  = false;

        % cezaları biraz dengeli tut
        cons_stage.pen.lambda.cav_frac  = 1.5;
        cons_stage.pen.lambda.dp_quant  = 0.5;
        cons_stage.pen.lambda.qsat_margin = 0.3;

        % guard rails (decode sonrası yine clamp ediliyor)
        hyd.Lh          = max(hyd.Lh, 3.5e-3);
        hyd.Vmin_fac    = max(hyd.Vmin_fac, 0.93);
        therm.resFactor = max(therm.resFactor, 12);

        cfg_stage.PF.gain = 0.80;
        cfg_stage.PF.tau  = 0.9;

    case 3  % termal/akışkan tamamlayıcı
        cons_stage.on.dp_quant    = true;
        cons_stage.on.cav_frac    = true;
        cons_stage.on.qsat_margin = true;
        cons_stage.on.thermal_dT  = true;

        cons_stage.pen.lambda.thermal_dT = 0;   % aşırı ısınmayı izle ama kilitleme
        cons_stage.pen.lambda.cav_frac   = 0.2;

        % dP_cap geniş; solver güvenliği için yerinde kalsın
        num.dP_cap = max(num.dP_cap, 3e8);

        cfg_stage.PF.gain = 0.80;
        cfg_stage.PF.tau  = 0.90;
end

    % --- Fitness sarıcı ---
   inner_fitfun = @(xx) eval_fitness_for_x(xx, ga.design_set, ...
    obj, cons_stage, pp.tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_stage, ...
    ga.opt.cache_enable, ga.opt.fail_early_k);



    fhandle = @(xx) compact_log_wrapper(xx, inner_fitfun);

    % --- GA seçenekleri ---
    opts = optimoptions('ga', ...
        'UseParallel', ga.opt.use_parallel, ...
        'InitialPopulationMatrix', P0, ...
        'PopulationSize', ga.opt.popsize, ...
        'MaxGenerations', ga.opt.maxgen, ...
        'Display','off');

    % --- Güvenli ön-tanımlar (xbest henüz yokken kullanmak YASAK) ---
    xbest = []; fbest = []; pop = []; scores = []; exitflag = [];

    % --- Çalıştır ---
    [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts);
    

    % --- Fallback: GA bir şey döndürmediyse orta noktayı kullan ---
    if isempty(xbest)
        xbest = 0.5*(lb+ub);
        fbest = inf;
        pop   = []; scores = [];
    end

    % --- Log & kalıcılaştır ---
    runLog = struct('stage',ga.design_set,'xbest',xbest,'fbest',fbest,'output',output, ...
                    'pop',pop,'scores',scores,'exitflag',exitflag, ...
                    'ga_options',opts,'bounds',struct('lb',lb,'ub',ub),'seed',ga.opt.seed);
    if isfield(ga.opt,'save_log') && ~isempty(ga.opt.save_log)
        try
            [pth,base,ext] = fileparts(ga.opt.save_log); if isempty(pth), pth='.'; end
            save(fullfile(pth, sprintf('%s_stage%d%s',base,ga.design_set,ext)), 'runLog');
        catch ME
            warning('runLog kaydedilemedi: %s', ME.message);
        end
    end

    % --- Bu aşamanın en iyisini uygula → bir sonraki aşama bunun üstünde arar ---
    ga.enable = true; ga.x = xbest;
    [geom, sh, orf, hyd, therm, num, ga] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num);
    ga.enable = false; ga.x = [];

    last_best = xbest; last_set = ga.design_set;
    fprintf('\n=== GA Stage %d Bitti ===\nBest f = %.6g | xbest = [%s]\n', ...
        ga.design_set, fbest, join(string(xbest.'),', '));

    % --- Sınır daraltma: bir SONRAKİ aşama aynı set ise
    if si < numel(stage_list)
        next_set = stage_list(si+1);
        try
            [lb_sh, ub_sh] = shrink_bounds_from_pop(pop, scores, lb, ub, ga.opt.keep_top, ga.opt.buffer_fac);
            if ~isempty(IntCon)
                lb_sh(IntCon) = ceil(lb_sh(IntCon));
                ub_sh(IntCon) = floor(ub_sh(IntCon));
                lb_sh = max(lb, lb_sh);
                ub_sh = min(ub, ub_sh);
            end
            if next_set == ga.design_set
                ga.opt.lb_override = lb_sh;
                ga.opt.ub_override = ub_sh;
            else
                ga.opt.lb_override = [];
                ga.opt.ub_override = [];
            end
        catch
            ga.opt.lb_override = [];
            ga.opt.ub_override = [];
        end
    end

    % --- Aşama bazlı nesil sayısı (varsa) ---
    if isfield(ga.opt,'maxgen_stage') && numel(ga.opt.maxgen_stage)>=si+1
        ga.opt.maxgen = ga.opt.maxgen_stage(si+1);
    end
end


    % Son aşamanın x’i overlay/raporlar için kalsın
    ga.enable = true; ga.x = last_best; ga.best_x = last_best; ga.best_set = last_set;
    ga_dbg = struct('enable',true,'design_set',ga.best_set,'x',ga.best_x);
[geom_dbg, ~, orf_dbg, hyd_dbg, therm_dbg, ~, ~] = decode_design_apply(ga_dbg, geom, sh, orf, hyd, therm, num);

fprintf('DBG set-%d: n_orf=%d | d_o=%.3f mm | Ao=%.3e m^2 | Lgap=%.1f mm | Vmin_fac=%.2f | Lh=%.3f mm | resFactor=%.0f\n', ...
    ga.best_set, orf_dbg.n_orf, 1e3*geom_dbg.d_o, orf_dbg.n_orf*(pi*geom_dbg.d_o^2/4), ...
    1e3*geom_dbg.Lgap, hyd_dbg.Vmin_fac, 1e3*hyd_dbg.Lh, therm_dbg.resFactor);

end



%% -------------------- (1) Tek koşu: görselleştirme/diagnostic ---------
rec = min(max(1, vis.sel_rec), R);
dir = upper(string(vis.sel_dir));
[t_sim, ag_sim, t5_sim, t95_sim] = pickA(sel.sim_source, rec, dir);

% Kuyruk ekle (sim)
dt_sim  = median(diff(t_sim));
t_tail  = (t_sim(end)+dt_sim:dt_sim:t_sim(end)+pp.tail_sec).';
t_sim   = [t_sim; t_tail];
ag_sim  = [ag_sim; zeros(size(t_tail))];

% PF rampası (sim penceresine göre) — guard
cfg = set_pf_ton_if_nan(cfg, t5_sim, 0.5);

fprintf('SIM source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
        sel.sim_source, rec, dir, numel(t_sim), dt_sim, t5_sim, t95_sim);

[x0_sim,a0_sim] = lin_MCK_consistent(t_sim, ag_sim, M, Cstr, K);
[xD_sim, aD_sim, dlog_sim, vD_sim] = mck_with_damper_adv( ...
    t_sim, ag_sim, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);

% Grafik seti (opsiyonel ikinci koşu)
use_plot_run = false;
cfg_plot = cfg;   % tern(...) çağrısı için hazır dursun
if vis.dual_run_if_needed && ~strcmpi(sel.plot_source, sel.sim_source)
    [t_plot, ag_plot, t5_plot, t95_plot] = pickA(sel.plot_source, rec, dir);
    dt_plot = median(diff(t_plot));
    t_tail  = (t_plot(end)+dt_plot:dt_plot:t_plot(end)+pp.tail_sec).';
    t_plot  = [t_plot; t_tail];
    ag_plot = [ag_plot; zeros(size(t_tail))];

    cfg_plot = set_pf_ton_if_nan(cfg_plot, t5_plot, 0.5);

    [x0,a0] = lin_MCK_consistent(t_plot, ag_plot, M, Cstr, K);
    [xD,aD,dlog,vD] = mck_with_damper_adv( ...
        t_plot, ag_plot, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg_plot);

    use_plot_run = true;
    fprintf('PLOT source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
            sel.plot_source, rec, dir, numel(t_plot), dt_plot, t5_plot, t95_plot);
else
    % Plot, sim ile aynı
t_plot=t_sim; ag_plot=ag_sim; t5_plot=t5_sim; t95_plot=t95_sim;
x0=x0_sim; a0=a0_sim; xD=xD_sim; aD=aD_sim; dlog=dlog_sim; vD=vD_sim;
dt_plot = dt_sim;   % başlıktaki dt için

end

% ---- Enerji/diagnostik (sadece yazdırma) ----
active_cfg  = tern(use_plot_run, cfg_plot, cfg);
dlog_active = tern(use_plot_run, dlog, dlog_sim);

E_orif  = dlog_active.E_cum(end);                 % orifis enerjisi (J)
P_struc = sum( (vD * Cstr) .* vD, 2 );
E_struc = trapz(t_plot, P_struc);

fprintf('E_orifice = %.3e J | E_struct = %.3e J | oran = %.3f\n', ...
        E_orif, E_struc, E_orif/max(E_struc,eps));

P_mech = sum( (dlog_active.F_story .* (vD(:,2:end)-vD(:,1:end-1))), 2 );
fprintf('⟨P_mech⟩ = %.3e W (negatif ise net sönüm)\n', mean(P_mech,'omitnan'));

fprintf('CHECK → leak=%.2e, Lh=%.3e, Ao=%.2e m^2, PF=%s\n', ...
    hyd.K_leak, hyd.Lh, orf.n_orf*(pi*geom.d_o^2/4), active_cfg.PF.mode);

fprintf('dP_orf q95 = %.3e Pa | Q95 = %.3e m^3/s\n', ...
    prctile(dlog_active.dP_orf_time_max,95), dlog_active.Q_abs_p95);

if isfield(dlog_active,'cav_margin_min')
    fprintf('cav_margin_min = %.1f kPa\n', 1e-3*dlog_active.cav_margin_min);
end


% === Tek FIGÜR (ilk koddaki stil): üstte yer değiştirme, altta mutlak ivme ===
if vis.make_plots
    j_mon = min(max(1, obj.idx_disp_story), size(xD,2));   % izlenen kat

    % Mutlak ivme (dampersiz/damperli)
    a_abs0 = a0 + ag_plot * ones(1, size(a0,2));
    a_absD = aD + ag_plot * ones(1, size(aD,2));

    figure('Color','w','Visible','on','Name', ...
        sprintf('%d. kat: dampersiz vs damperli (TERMAL+Cd(Re)+Kaçak+Vmin+NUM+β(T)+P_{sat})', j_mon));

    % --- Üst: yer değiştirme ---
    subplot(2,1,1); hold on; grid on;
    plot(t_plot, x0(:, j_mon), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, xD(:, j_mon), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('x_{%d} [m]', j_mon)); legend('show','Location','best');
    title(sprintf('N=%d, dt=%.4fs | d_o=%.1f mm, gain=%.2f, \\tau=%.2f s', ...
        numel(t_plot), dt_plot, 1e3*geom.d_o, active_cfg.PF.gain, active_cfg.PF.tau));

    % --- Alt: mutlak ivme ---
    subplot(2,1,2); hold on; grid on;
    plot(t_plot, a_abs0(:, j_mon), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, a_absD(:, j_mon), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('a_{%d} [m/s^2]', j_mon)); xlabel('t [s]'); legend('show','Location','best');
end

%% -------------------- (2) Amaç fonksiyonu — kompakt (ilk koddaki gibi) ----
% Girdi olarak şunlar zaten mevcut olmalı:
% t_plot, dt_plot, vD, dlog_active, M, K, k_sd, obj, cons, orf, hyd, therm, num, t5_plot, t95_plot

% Kapasiteler
dPcap_eff = (isfield(num,'dP_cap') && isfinite(num.dP_cap)) * num.dP_cap + ...
            (~(isfield(num,'dP_cap') && isfinite(num.dP_cap))) * 3e8;

cd_ref   = max(orf.CdInf, orf.Cd0);
Ae_ref   = cd_ref * max(orf.Ao_eff, 1e-12);                 % paralel ve Cd dahil
rho_ref  = max(therm.rho_ref, 100);
Qcap_ref = hyd.Vmin_fac * Ae_ref * sqrt(2*dPcap_eff / rho_ref);

% Gözlenenler (dlog_active üzerinden)
Qp95_max = dlog_active.Q_abs_p95;                            % m^3/s
dp95_max = prctile(dlog_active.dP_orf_time_max, 95);         % Pa

% Uygunluk bayrakları
Qmargin = (isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin') && isfinite(cons.hyd.Q_margin)) ...
            * cons.hyd.Q_margin + (~(isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin'))) * 0.90;
okQ  = (Qp95_max <= Qmargin * Qcap_ref);
okdp = (dp95_max <= dPcap_eff);

% zeta_eq tahmini (1. mod)
[PHI,LAM] = eig(K,M);
WW = sqrt(builtin('diag', LAM));               % 'diag' fonksiyonu gölgelenmesin diye builtin

w1   = min(WW);
phi1 = PHI(:, WW==w1); phi1 = phi1(:,1);
phi1 = phi1 / max(abs(phi1));
m_eff = phi1.'*M*phi1;
k_eff = phi1.'*K*phi1;

% Arias penceresinde izlenen katın hız RMS'i
j_mon = min(max(1, obj.idx_disp_story), size(vD,2));
maskA = (t_plot >= t5_plot) & (t_plot <= t95_plot);
v_rms = sqrt(mean(vD(maskA, j_mon).^2, 'omitnan'));
Xeq   = v_rms / max(w1, 1e-6);
Est   = 0.5 * (k_eff + k_sd) * Xeq^2;

% Mekanik güçten efektif sönüm (yalnız sönümleyici kısım)
P_mech  = sum( dlog_active.F_story .* (vD(:,2:end) - vD(:,1:end-1)), 2 );
P_diss  = max(-P_mech, 0);                                  % negatif → sönüm kabulü
zeta_eq = sum(P_diss) * dt_plot / max(4*pi*Est, eps);

% <P_mech> işareti (pozitif ⇒ net sönüm)
Pmech_avg = -mean(P_mech, 'omitnan');
okP = (Pmech_avg > 0);

% Skor
score = zeta_eq ...
        - 0.7*double(~okP) ...
        - 0.5*double(~okQ) ...
        - 0.5*double(~okdp);

fprintf('\n=== AMAÇ (kompakt) ===\n');
fprintf('zeta_eq=%.3f | <P_mech>_diss=%.3e W | Q95=%.3e | dp95=%.3e | ok=[Q %d, dp %d, P %d] | score=%.3f\n', ...
    zeta_eq, Pmech_avg, Qp95_max, dp95_max, okQ, okdp, okP, score);

% GA minimizasyonu ile uyum için amaç değeri
J = -score;

%% -------------------- (3) Kısıtlar: Penalty hesap (R kayıt) ----------
[Penalty, cons_detail] = evaluate_constraints_over_records( ...
    cons, cons.src_for_constraints, obj, pp.tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, tern(ga.enable,ga.design_set,0), tern(ga.enable,ga.x,[]));

Fitness = J + Penalty;
fprintf('λpen = [tau=%.2f, stroke=%.2f, dT=%.2f, cav=%.2f]\n', ...
    cons.pen.lambda.spring_tau, cons.pen.lambda.stroke, ...
    cons.pen.lambda.thermal_dT, cons.pen.lambda.cav_frac);
fprintf('PF: mode=%s, t_on=%.2fs, tau=%.2f, gain=%.2f\n', ...
    cfg.PF.mode, cfg.PF.t_on, cfg.PF.tau, cfg.PF.gain);

fprintf('\n================ Kısıt & Ceza Sonuçları ================\n');
fprintf('Penalty = %.6g  |  Fitness = J + Penalty = %.6g\n', Penalty, Fitness);
% === Kısıt/cihaz özeti → CSV ===
try, mkdir('out'); end

ratio = cons_detail.ratios;
okflag = @(r) (r<=1+1e-12);

Names   = {'spring_tau','spring_slender','stroke','force_cap','dp_quant','thermal_dT','cav_frac','qsat_margin'};
Limits  = [cons.spring.tau_allow, cons.spring.lambda_max, cons.stroke.util_factor*geom.Lgap, ...
           cons.force.F_cap, num.dP_cap, cons.thermal.cap_C, cons.hyd.cav_frac_cap, cons.hyd.Q_margin * getfield_default(num,'Qcap_big', ...
           0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) )];

Ratios  = [ratio.spring_tau, ratio.spring_slender, ratio.stroke, ratio.force_cap, ...
           ratio.dp_quant, ratio.thermal_dT, ratio.cav_frac, ratio.qsat_margin];

Flags = strings(numel(Ratios),1);
okmask = Ratios <= (1 + 1e-12);
Flags(okmask)  = "OK";
Flags(~okmask) = "VIOL";


Vals_Fmax   = max(cons_detail.Fmax_records,   [], 'omitnan');
Vals_stroke = max(cons_detail.stroke_records, [], 'omitnan');
Vals_dpq    = max(cons_detail.dpq_records,    [], 'omitnan');
Vals_dT     = max(cons_detail.dT_records,     [], 'omitnan');
Vals_cav    = max(cons_detail.cav_records,    [], 'omitnan');
Vals_Qp95   = max(cons_detail.Qp95_records,   [], 'omitnan');

T = table(Names.', Ratios.', Flags, ...
    'VariableNames', {'Constraint','Ratio','Flag'});

% Ek cihaz metrikleri ayrı tablo: (zarf değerleri)
T_dev = table(Vals_Fmax, Vals_stroke, Vals_dpq, Vals_dT, Vals_cav, Vals_Qp95, ...
    'VariableNames', {'Fmax_N','stroke_m','dp_q95_Pa','dT_C','cav95','Q95_m3s'});

writetable(T,    fullfile('out','cons_summary.csv'));
writetable(T_dev,fullfile('out','device_summary.csv'));
fprintf('CSV yazıldı: out/cons_summary.csv, out/device_summary.csv\n');


% ---------- Yardımcılar ----------
hinge  = @(r) max(0, r - 1);
norm0  = @(x) (abs(x) < 1e-12) * 0 + (abs(x) >= 1e-12) .* x;  % -0.000 yerine 0.000
pwr    = cons.pen.power;
lam    = cons.pen.lambda;

pen_sum = 0;   % bileşenlerden yeniden hesaplanan toplam (kontrol amaçlı)

% ---------- Özet oranlar + bireysel ceza katkıları ----------
cdt = cons_detail;  % kısaltma

if cons.on.spring_tau
    r = norm0(cdt.ratios.spring_tau);
    pen_i = lam.spring_tau * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('τ_max/τ_allow         = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.spring_slender
    r = norm0(cdt.ratios.spring_slender);
    pen_i = lam.spring_slender * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('L_free/D_m / λ_max    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.stroke
    r = norm0(cdt.ratios.stroke);
    pen_i = lam.stroke * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('stroke/(0.9*L_gap)    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.force_cap && isfinite(cons.force.F_cap)
    r = norm0(cdt.ratios.force_cap);
    pen_i = lam.force_cap * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('F_max/F_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.dp_quant
    r = norm0(cdt.ratios.dp_quant);
    pen_i = lam.dp_quant * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('q_Δp/dP_cap           = %.3f (q=%.3f)   [pen=%.3g]\n', r, cons.dp.q, pen_i);
end

if cons.on.thermal_dT
    r = norm0(cdt.ratios.thermal_dT);
    pen_i = lam.thermal_dT * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('ΔT/ΔT_cap             = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.cav_frac
    r = norm0(cdt.ratios.cav_frac);
    pen_i = lam.cav_frac * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('cav95/γ_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.qsat_margin
    r = norm0(cdt.ratios.qsat_margin);
    pen_i = lam.qsat_margin * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('Q95/(margin*Qcap_big) = %.3f   [pen=%.3g]\n', r, pen_i);
end

% BigM (sim başarısızlığı) bilgilendirmesi ve katkısı
if cons.on.fail_bigM && cdt.any_fail
    pen_BigM = cons.pen.bigM * lam.fail_bigM;
    pen_sum  = pen_sum + pen_BigM;
    fprintf('Sim başarısız (≥1 koşu): BigM cezası eklendi. [pen=%.3g]\n', pen_BigM);
end

fprintf('--- Penalty (bileşenlerden) ≈ %.6g\n', pen_sum);

% ---------- İsteğe bağlı tanılama: zarf/metrik özetleri ----------
if isfield(cons_detail,'Fmax_records')
    % Kayıt zarfı / agregeler (dpq_all için agresyon kuralını uygula)
    Fmax_all   = max(cons_detail.Fmax_records,   [], 'omitnan');
    stroke_all = max(cons_detail.stroke_records, [], 'omitnan');
    if isfield(cons,'dp') && isfield(cons.dp,'agg') && strcmpi(cons.dp.agg,'cvar')
        dpq_all = cvar_from_samples(cons_detail.dpq_records(:), cons.alpha_CVaR_cons);
    else
        dpq_all = max(cons_detail.dpq_records, [], 'omitnan');
    end
    dT_all   = max(cons_detail.dT_records,   [], 'omitnan');
    cav_all  = max(cons_detail.cav_records,  [], 'omitnan');
    Qp95_all = max(cons_detail.Qp95_records, [], 'omitnan');

    fprintf('--- Tanılama (zarf değerleri) ---\n');
    fprintf('Fmax=%.3e N | stroke=%.3e m | dpq=%.3e Pa | ΔT=%.2f C | cav95=%.3f | Q95=%.3e m^3/s\n', ...
        Fmax_all, stroke_all, dpq_all, dT_all, cav_all, Qp95_all);
end



% ======================================================================
%                             FONKSİYONLAR
% ======================================================================
function [lb2,ub2] = shrink_bounds_from_pop(pop, scores, lb, ub, keep_top, buf)
    if isempty(pop) || isempty(scores)
        lb2 = lb; ub2 = ub; return;
    end
    [~,ix] = sort(scores(:),'ascend');                 % en iyi küçük
    K = max(1, ceil(keep_top * size(pop,1)));
    P = pop(ix(1:K), :);                               % elitler
    p10 = prctile(P,10,1);
    p90 = prctile(P,90,1);
    span = max(p90 - p10, 1e-12);
    lb2 = max(lb, p10 - buf.*span);
    ub2 = min(ub, p90 + buf.*span);
end

function y = tern(cond,a,b)
    if cond
        y = a;
    else
        y = b;
    end
end


function [J, out] = compute_objective_over_records( ...
src, obj, tail_sec, ...
t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end
    ...
    % local_design_ratios_one_dir içinde:
    % Simülasyon penceresi: aynı veriye kuyruk ekleyelim (PF ramp guard uyumu)
    % Kaynak seçici
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R = numel(tX);
    out(R) = struct('d_rel',NaN,'a_rel',NaN,'J_r',NaN);

    % Önce REFERANSLAR: (damper yok, μ=1.0) → her kayıt & (X,Y mevcutsa) yön için sabitlenir
    ref = struct('X',struct('d',nan(R,1),'a',nan(R,1)), ...
                 'Y',struct('d',nan(R,1),'a',nan(R,1)));

    for r=1:R
        % X yönü referans
        if ~isempty(aX{r})
            [dref, aref] = local_ref_metrics(tX{r}, aX{r}, t5x(r), t95x(r), obj, M,Cstr,K, n);
            ref.X.d(r)=dref; ref.X.a(r)=aref;
        end
        % Y yönü referans (varsa)
        if ~isempty(aY{r})
            [dref, aref] = local_ref_metrics(tY{r}, aY{r}, t5y(r), t95y(r), obj, M,Cstr,K, n);
            ref.Y.d(r)=dref; ref.Y.a(r)=aref;
        end
    end

    % Sonra TASARIM: her kayıt → yön zarfı → μ-agg → J_r
    for r=1:R
        [d_rel_X, a_rel_X] = local_design_ratios_one_dir('X', r, ...
            tX, aX, t5x, t95x, tY, aY, t5y, t95y, ref, tail_sec, obj, ...
            design_set, x_ga, M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg);
        [d_rel_Y, a_rel_Y] = local_design_ratios_one_dir('Y', r, ...
            tX, aX, t5x, t95x, tY, aY, t5y, t95y, ref, tail_sec, obj, ...
            design_set, x_ga, M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg);

        switch lower(obj.dir_mode)
            case 'xonly', d_rel = d_rel_X; a_rel = a_rel_X;
            case 'yonly', d_rel = d_rel_Y; a_rel = a_rel_Y;
            otherwise     % 'envelope'
                d_rel = max([d_rel_X, d_rel_Y], [], 'omitnan');
                a_rel = max([a_rel_X, a_rel_Y], [], 'omitnan');
        end

        % Ağırlıklı toplam
        J_r = obj.weights_da(1)*d_rel + obj.weights_da(2)*a_rel;

        out(r).d_rel = d_rel;
        out(r).a_rel = a_rel;
        out(r).J_r   = J_r;
    end

    % CVaR(α) hesap
    alpha = min(max(obj.alpha_CVaR,eps),0.99);
    Jlist = [out.J_r].';
    J = cvar_from_samples(Jlist, alpha);

end

function [d_ref, a_ref] = local_ref_metrics(t, ag, t5, t95, obj, M, C, K, n)
    % Kuyruk eklemeden baz referans (damper yok)
    [x0,a0] = lin_MCK_consistent(t, ag, M, C, K);
    d_ref = max(abs(x0(:,min(obj.idx_disp_story,n))));
    a_ref = acc_metric_from_series(t, a0(:,min(obj.idx_acc_story,n)), t5, t95, obj);
    d_ref = max(d_ref, eps);
    a_ref = max(a_ref, eps);
end

function [d_rel_dir, a_rel_dir] = local_design_ratios_one_dir(which, r, ...
    tX, aX, t5x, t95x, tY, aY, t5y, t95y, ref, tail_sec, obj, ...
    design_set, x_ga, M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg)

    switch upper(string(which))
        case "X"
            t=tX{r}; ag=aX{r}; t5=t5x(r); t95=t95x(r);
            d_ref=ref.X.d(r); a_ref=ref.X.a(r);
        case "Y"
            t=tY{r}; ag=aY{r}; t5=t5y(r); t95=t95y(r);
            d_ref=ref.Y.d(r); a_ref=ref.Y.a(r);
    end
    if isempty(ag) || ~isfinite(d_ref) || ~isfinite(a_ref)
        d_rel_dir = NaN; a_rel_dir = NaN; return;
    end

    % --- Simülasyon penceresi: aynı veriye kuyruk ekleyelim ---
    dt = median(diff(t));
    t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
    t_s    = [t;  t_tail];
    ag_s   = [ag; zeros(size(t_tail))];

    % PF ramp t_on: Arias t5 + 3
    cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);

    % μ senaryoları
    mus   = obj.mu_scenarios(:).';
    d_vals = nan(size(mus));  a_vals = nan(size(mus));

    for k = 1:numel(mus)
        resp = simulate( ...
            design_set, x_ga, mus(k), t_s, ag_s, ...
            M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir);

        if ~resp.ok
            fail_ratio = 5;              % istersen 3–10 arası dene
            d_vals(k) = fail_ratio * d_ref;
            a_vals(k) = fail_ratio * a_ref;
            continue;
        end

        x  = resp.y.x(1:numel(t), :);   % kuyruğu at
        aR = resp.y.a(1:numel(t), :);

        d  = max(abs(x(:,min(obj.idx_disp_story,n))));
        am = acc_metric_from_series(t, aR(:,min(obj.idx_acc_story,n)), t5, t95, obj);

        d_vals(k) = d;  a_vals(k) = am;
    end

    % μ-aggregation
    switch lower(obj.mu_aggregate)
        case 'weighted'
            w = obj.mu_weights(:); w = w/sum(w);
            d_agg = nansum(w.*d_vals(:));
            a_agg = nansum(w.*a_vals(:));
        otherwise
            d_agg = max(d_vals, [], 'omitnan');
            a_agg = max(a_vals, [], 'omitnan');
    end

    d_rel_dir = d_agg / max(d_ref, eps);
    a_rel_dir = a_agg / max(a_ref, eps);
end

function val = acc_metric_from_series(t, a, t5, t95, obj)
    % Arias penceresi içi metrikler
    if obj.use_arias_window
        w = (t>=t5 & t<=t95);
    else
        w = true(size(t));
    end
    aa = a(w); tt = t(w);
    if isempty(aa) || numel(aa)<2
        val = NaN; return;
    end
    switch lower(obj.acc_metric)
        case 'energy'
            val = trapz(tt, aa.^2);                 % enerji
        case 'rms+p95'
            rmsv = sqrt(mean(aa.^2,'omitnan'));
            p95  = prctile(abs(aa),95);
            val  = rmsv + obj.p95_penalty_w * p95;  % hibrit
        otherwise % 'rms'
            val = sqrt(mean(aa.^2,'omitnan'));
    end
end

function v = cvar_from_samples(x, alpha)
    % Örneklerden CVaR_α (Average Value-at-Risk)
    x = x(:); x = x(isfinite(x));
    if isempty(x), v = NaN; return; end
    q = quantile(x, 1-alpha);
    tail = x(x>=q);  % büyük-kötü kuyruk
    if isempty(tail), v = q; else, v = mean(tail); end
end
function [J1, out] = compute_J1_IDR_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    % Optional GA args
    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end

    % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R  = numel(tX);
    Ns = n - 1;
    if isscalar(h_story_m), h_story = repmat(h_story_m, Ns, 1);
    else,                   h_story = h_story_m(:); end

    out = struct('d_rel', nan(R,1));
    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü referans IDR
        dref_X = NaN; dagg_X = NaN;
        if ~isempty(aX{r})
            [dref_X, dagg_X] = local_idr_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), ...
                M, Cstr, K, n, h_story, tail_sec, cfg, mus, isWeighted, wmu, ...
                design_set, x_ga, geom, sh, orf, hyd, therm, num);
        end

        % --- Y yönü (varsa)
        dref_Y = NaN; dagg_Y = NaN;
        if ~isempty(aY{r})
            [dref_Y, dagg_Y] = local_idr_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), ...
                M, Cstr, K, n, h_story, tail_sec, cfg, mus, isWeighted, wmu, ...
                design_set, x_ga, geom, sh, orf, hyd, therm, num);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', d_ref = dref_X; d_agg = dagg_X;
            case 'yonly', d_ref = dref_Y; d_agg = dagg_Y;
            otherwise
                d_ref = max([dref_X, dref_Y], [], 'omitnan');
                d_agg = max([dagg_X, dagg_Y], [], 'omitnan');
        end

        out.d_rel(r) = d_agg / max(d_ref, eps);
    end

    J1 = cvar_from_samples(out.d_rel, obj.alpha_CVaR);

end

function [d_ref, d_agg] = local_idr_ref_and_agg(t, ag, t5, t95, ...
    M, Cstr, K, n, h_story, tail_sec, cfg, mus, isWeighted, wmu, ...
    design_set, x_ga, geom, sh, orf, hyd, therm, num)

    % Referans (damper yok, μ=1.0), pencere içinde IDR zarfı
    [x0,~] = lin_MCK_consistent(t, ag, M, Cstr, K);
    w = (t >= t5 & t <= t95);
    drift0 = x0(:,2:end) - x0(:,1:end-1);
    idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)');
    d_ref  = max(idr0,[],'all');

    % Tasarım: kuyruk ekle + PF ramp guard
    dt     = median(diff(t));
    t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
    t_s    = [t; t_tail];
    ag_s   = [ag; zeros(size(t_tail))];

    cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);

    vals = nan(size(mus));
    for k = 1:numel(mus)
        resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                        M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir);
        if ~resp.ok
            vals(k) = 5 * d_ref;
            continue;
        end
        xD   = resp.y.x(1:numel(t), :);
        driftD = xD(:,2:end) - xD(:,1:end-1);
        idrD   = abs(driftD(w,:)) ./ (ones(sum(w),1) * h_story(:)');
        vals(k)= max(idrD,[],'all');
    end

    if isWeighted
        d_agg = nansum(wmu .* vals(:));
    else
        d_agg = max(vals, [], 'omitnan');
    end
end


function [J2, out] = compute_J2_ACC_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end
     
 
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    w_p = obj.p95_penalty_w;
    out = struct('a_rel', nan(R,1));

    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü: A_env ref ve agg
        aref_X = NaN; aagg_X = NaN;
        if ~isempty(aX{r})
            [aref_X, aagg_X] = local_acc_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), ...
                w_p, M, Cstr, K, n, tail_sec, cfg, mus, isWeighted, wmu, ...
                design_set, x_ga, geom, sh, orf, hyd, therm, num);
        end

        % --- Y yönü (varsa)
        aref_Y = NaN; aagg_Y = NaN;
        if ~isempty(aY{r})
            [aref_Y, aagg_Y] = local_acc_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), ...
                w_p, M, Cstr, K, n, tail_sec, cfg, mus, isWeighted, wmu, ...
                design_set, x_ga, geom, sh, orf, hyd, therm, num);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', a_ref = aref_X; a_agg = aagg_X;
            case 'yonly', a_ref = aref_Y; a_agg = aagg_Y;
            otherwise
                a_ref = max([aref_X, aref_Y], [], 'omitnan');
                a_agg = max([aagg_X, aagg_Y], [], 'omitnan');
        end

        out.a_rel(r) = a_agg / max(a_ref, eps);
    end

    J2 = cvar_from_samples(out.a_rel, obj.alpha_CVaR);

end

function [A_ref, A_agg] = local_acc_ref_and_agg(t, ag, t5, t95, ...
    w_p, M, Cstr, K, n, tail_sec, cfg, mus, isWeighted, wmu, ...
    design_set, x_ga, geom, sh, orf, hyd, therm, num)

    % Referans (damper yok) mutlak ivme zarfı (RMS+p95)
    [~,a_rel0] = lin_MCK_consistent(t, ag, M, Cstr, K);  % relatif
    a_abs0 = a_rel0 + ag(:) * ones(1, size(a_rel0,2));  % her kat için mutlak
    w = (t >= t5 & t <= t95);
    A_i = zeros(1,size(a_abs0,2));
    for i=1:size(a_abs0,2)
        ai = a_abs0(w,i);
        rmsv = sqrt(mean(ai.^2,'omitnan'));
        p95  = prctile(abs(ai),95);
        A_i(i) = rmsv + w_p * p95;
    end
    A_ref = max(A_i);   % kat zarfı

    % Tasarım: kuyruk ekle + PF guard
    dt     = median(diff(t));
    t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
    t_s    = [t; t_tail];
    ag_s   = [ag; zeros(size(t_tail))];
    cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);

    vals = nan(size(mus));
    for k = 1:numel(mus)
        resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                        M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir);
        if ~resp.ok
            vals(k) = 5 * A_ref;  % fail durumunda güvenli büyük ceza
            continue;
        end
        a_relD = resp.y.a(1:numel(t), :);           % relatif, kuyruğu at
        a_absD = a_relD + ag(:) * ones(1,size(a_relD,2));

        Ai = zeros(1,size(a_absD,2));
        for i=1:size(a_absD,2)
            ai = a_absD(w,i);
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            Ai(i)= rmsv + w_p * p95;
        end
        vals(k) = max(Ai);   % kat zarfı
    end

    if isWeighted
        A_agg = nansum(wmu .* vals(:));
    else
        A_agg = max(vals, [], 'omitnan');
    end
end

function [Penalty, out] = evaluate_constraints_over_records( ...
    cons, src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end
    ...

      % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    mus = cons.mu_scenarios(:).';

    % ----- Erken çıkış sayaç ayarları -----
    fail_early_k = 0;
    if isfield(cons,'fail_early_k') && ~isempty(cons.fail_early_k) && isfinite(cons.fail_early_k)
        fail_early_k = max(0, round(cons.fail_early_k));
    end
    fail_count_global = 0;
    early_break_all   = false;

    % Toplayıcılar
    Fmax_records   = nan(R,1);
    stroke_records = nan(R,1);
    dpq_records    = nan(R,1);
    dT_records     = nan(R,1);
    cav_records    = nan(R,1);
    Qp95_records   = nan(R,1);
    any_fail_rec   = false(R,1);

    for r = 1:R
        % ---- Yön zarfı için X/Y ölçüleri
        rec_metrics = struct('Fmax',[],'stroke',[],'dpq',[],'dT',[],'cav95',[],'Qp95',[],'fail',[]);

        for dir = ["X","Y"]
            % Y yönü yoksa atla
            if dir=="Y" && (isempty(aY{r}) || all(isnan(aY{r})))
                continue;
            end

            % Kayıt/yön serileri
            if dir=="X"
                t = tX{r}; ag = aX{r}; t5=t5x(r); t95=t95x(r);
            else
                t = tY{r}; ag = aY{r}; t5=t5y(r); t95=t95y(r);
            end

            % ---- Tail ekle (simülasyon penceresi)
            dt    = median(diff(t));
            t_tail= (t(end)+dt : dt : t(end)+tail_sec).';
            t_s   = [t; t_tail];
            ag_s  = [ag; zeros(size(t_tail))];

            % ---- PF ramp koruması
cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);


            % ---- μ-senaryosu zarfı
            Fmax_mu   = -Inf; stroke_mu = -Inf; dpq_mu = -Inf; dT_mu = -Inf; cav_mu = -Inf; Qp95_mu = -Inf;
            fail_mu   = false;

            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ (BURADA) <<<<<<<<<<<<<<<<<
            for k = 1:numel(mus)
                resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);

                if ~resp.ok
                    fail_mu = true;
                    fail_count_global = fail_count_global + 1;

                    if fail_early_k > 0 && fail_count_global >= fail_early_k
                        early_break_all = true;   % tüm döngülerden çık
                        break;                    % μ-döngüsünü kır
                    end
                    continue;  % sıradaki μ
                end

                % ---- ölçümler (başarılı koşu)
                Fmax_mu   = max(Fmax_mu,   resp.F_max);
                stroke_mu = max(stroke_mu, resp.stroke_max);
                dpq_mu    = max(dpq_mu,    resp.dP_q_time(cons.dp.q));
                dT_mu     = max(dT_mu,     resp.dT_est);
                mask = (t_s >= t5) & (t_s <= t95);
    if any(mask)
        cav_p95_win = prctile(resp.diag.cav_frac_t(mask), 95);
    else
        cav_p95_win = prctile(resp.diag.cav_frac_t, 95);  % emniyetli geri dönüş
    end
    cav_mu = max(cav_mu, cav_p95_win);
                Qp95_mu   = max(Qp95_mu,   resp.metrics.Q_abs_p95);
            end
            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ BİTTİ <<<<<<<<<<<<<<<<<<<<

            % μ-döngüsü eşik nedeniyle kırıldıysa, yön döngüsünü de bırak
            if early_break_all
                break;
            end

            % Bu yönün metriklerini biriktir
            rec_metrics.Fmax   = [rec_metrics.Fmax,   Fmax_mu];
            rec_metrics.stroke = [rec_metrics.stroke, stroke_mu];
            rec_metrics.dpq    = [rec_metrics.dpq,    dpq_mu];
            rec_metrics.dT     = [rec_metrics.dT,     dT_mu];
            rec_metrics.cav95  = [rec_metrics.cav95,  cav_mu];
            rec_metrics.Qp95   = [rec_metrics.Qp95,   Qp95_mu];
            rec_metrics.fail   = [rec_metrics.fail,   fail_mu];
        end

        % Yön döngüsünden erken çıktıysak kalan kayıtları INF kabul et
        if early_break_all
    any_fail_rec(r:end) = true;
    break;
end

if early_break_all
    Penalty = cons.pen.bigM * cons.pen.lambda.fail_bigM;
    out = struct();
    out.ratios = struct('spring_tau',0,'spring_slender',0,'stroke',0, ...
                        'force_cap',0,'dp_quant',0,'thermal_dT',0, ...
                        'cav_frac',0,'qsat_margin',0);
    out.any_fail = true;
    out.dpq_records    = dpq_records;
    out.Fmax_records   = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records     = dT_records;
    out.cav_records    = cav_records;
    out.Qp95_records   = Qp95_records;
    return
end


        % ---- Yön zarfı (X,Y → max)
        Fmax_records(r)   = max(rec_metrics.Fmax,   [], 'omitnan');
        stroke_records(r) = max(rec_metrics.stroke, [], 'omitnan');
        dpq_records(r)    = max(rec_metrics.dpq,    [], 'omitnan');
        dT_records(r)     = max(rec_metrics.dT,     [], 'omitnan');
        cav_records(r)    = max(rec_metrics.cav95,  [], 'omitnan');
        Qp95_records(r)   = max(rec_metrics.Qp95,   [], 'omitnan');
        any_fail_rec(r)   = any(rec_metrics.fail);
    end


    % ---- Kayıt agregasyonları
    % Fmax ve stroke en-kötü kayıt
    Fmax_all   = max(Fmax_records,   [], 'omitnan');
    stroke_all = max(stroke_records, [], 'omitnan');
    dT_all     = max(dT_records,     [], 'omitnan');
    cav_all    = max(cav_records,    [], 'omitnan');
    Qp95_all   = max(Qp95_records,   [], 'omitnan');

    % Δp quantile: kayıtlar arası 'max' veya 'CVaR'
    switch lower(cons.dp.agg)
        case 'cvar'
            dpq_all = cvar_from_samples(dpq_records(:), cons.alpha_CVaR_cons);
        otherwise
            dpq_all = max(dpq_records, [], 'omitnan');
    end

    % ---- Oranlar (≥1 → ihlal)
    ratios = struct();

   % K1: yay kesme gerilmesi (yalnız yay kolu kuvveti)
if cons.on.spring_tau
    k_p_est = sh.G*sh.d_w^4 / (8*sh.n_turn*sh.D_m^3);           % [N/m]
    Fspring_max = k_p_est * stroke_all;                          % [N]
    ratios.spring_tau = spring_tau_ratio(Fspring_max, sh, cons.spring.tau_allow);
else
    ratios.spring_tau = 0;
end

    % K2: burkulma/serbest boy oranı
    if cons.on.spring_slender
        if strcmpi(cons.spring.L_free_mode,'fixed') && isfinite(cons.spring.L_free_fix)
            L_free = cons.spring.L_free_fix;
        else
            L_free = cons.spring.L_free_auto_fac * geom.Lgap; % yaklaşıklama
        end
        lambda = (L_free / max(sh.D_m,eps)) / max(cons.spring.lambda_max,eps);
        ratios.spring_slender = lambda;
    else
        ratios.spring_slender = 0;
    end

    % K3: strok
    if cons.on.stroke
        ratios.stroke = stroke_all / max(cons.stroke.util_factor*geom.Lgap, eps);
    else
        ratios.stroke = 0;
    end

    % K4: cihaz kuvvet limiti
    if cons.on.force_cap && isfinite(cons.force.F_cap)
        ratios.force_cap = Fmax_all / max(cons.force.F_cap, eps);
    else
        ratios.force_cap = 0;
    end

    % K5: Δp quantile
    if cons.on.dp_quant
        ratios.dp_quant = dpq_all / max(num.dP_cap, eps);
    else
        ratios.dp_quant = 0;
    end

    % K6: termal ΔT
    if cons.on.thermal_dT
        ratios.thermal_dT = dT_all / max(cons.thermal.cap_C, eps);
    else
        ratios.thermal_dT = 0;
    end

    % K7: kavitasyon
    if cons.on.cav_frac
        ratios.cav_frac = cav_all / max(cons.hyd.cav_frac_cap, eps);
    else
        ratios.cav_frac = 0;
    end

    % K8: Q satürasyon marjı
    if cons.on.qsat_margin
        Qcap = getfield_default(num,'Qcap_big', ...
    0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
ratios.qsat_margin = Qp95_all / max(cons.hyd.Q_margin * Qcap, eps);

    else
        ratios.qsat_margin = 0;
    end

    % ---- Ceza hesabı
    hinge = @(r) max(0, r - 1);
    pwr   = cons.pen.power;

    pen = 0;
    if cons.on.spring_tau,     pen = pen + cons.pen.lambda.spring_tau     * hinge(ratios.spring_tau    )^pwr; end
    if cons.on.spring_slender, pen = pen + cons.pen.lambda.spring_slender * hinge(ratios.spring_slender)^pwr; end
    if cons.on.stroke,         pen = pen + cons.pen.lambda.stroke         * hinge(ratios.stroke        )^pwr; end
    if cons.on.force_cap && isfinite(cons.force.F_cap)
                               pen = pen + cons.pen.lambda.force_cap      * hinge(ratios.force_cap     )^pwr; end
    if cons.on.dp_quant,       pen = pen + cons.pen.lambda.dp_quant       * hinge(ratios.dp_quant      )^pwr; end
    if cons.on.thermal_dT,     pen = pen + cons.pen.lambda.thermal_dT     * hinge(ratios.thermal_dT    )^pwr; end
    if cons.on.cav_frac,       pen = pen + cons.pen.lambda.cav_frac       * hinge(ratios.cav_frac      )^pwr; end
    if cons.on.qsat_margin,    pen = pen + cons.pen.lambda.qsat_margin    * hinge(ratios.qsat_margin   )^pwr; end

    any_fail = any(any_fail_rec);
    if cons.on.fail_bigM && any_fail
        pen = pen + cons.pen.bigM * cons.pen.lambda.fail_bigM;
    end

    Penalty = pen;

    % Çıkış detayları
    out = struct();
    out.ratios    = ratios;
    out.any_fail  = any_fail;
    out.dpq_records = dpq_records;
    out.Fmax_records = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records = dT_records;
    out.cav_records = cav_records;
    out.Qp95_records = Qp95_records;
end
function ratio = spring_tau_ratio(Fmax, sh, tau_allow)
    C  = max(sh.D_m, eps) / max(sh.d_w, eps);
    Kw = (4*C - 1)/(4*C - 4) + 0.615/C;
    tau_max = (8 * Fmax * max(sh.D_m,eps) / max(pi*sh.d_w^3, eps)) * Kw;
    ratio = tau_max / max(eps, tau_allow);
end

function P = lhs_population(lb,ub,N)
    D = numel(lb); P = zeros(N,D);
    for d=1:D
        edges = linspace(0,1,N+1);
        centers = (edges(1:end-1)+edges(2:end))/2;
        P(:,d) = lb(d) + centers(randperm(N))' .* (ub(d)-lb(d));
    end
end

function [f, aux] = eval_fitness_for_x(x, design_set, ...
    obj, cons, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, ...
    use_cache, fail_early_k)

    % ---- Önbellek anahtarı
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType','char','ValueType','any'); end
    key = sprintf('set%d|%s', design_set, mat2str(x,8));

    if use_cache && isKey(memo,key)
        data = memo(key); f = data.f; aux = data.aux; return;
    end

    % ---- Tasarıma uygula ve değerlendir
    try
        % GA decode’u simulate içinde değil, doğrudan burada yapmaya gerek yok;
        % compute/evaluate fonksiyonları simulate’i çağırırken set/x’ı geçiriyoruz.

        % Amaç: J
        goal_src = tern(obj.use_scaled_for_goal,'scaled','raw');
        [J, ~] = compute_objective_over_records( ...
    goal_src, obj, tail_sec, ...
            t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
            t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);

        % Kısıt: Penalty (erken çıkış desteği ile)
        cons_loc = cons; cons_loc.fail_early_k = fail_early_k;
        [Penalty, cons_detail] = evaluate_constraints_over_records( ...
            cons_loc, cons.src_for_constraints, obj, tail_sec, ...
            t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
            t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);
      
% === J1 & J2 (split) hesap — Pareto günlüğü için ===
[J1_split, J2_split] = compute_objectives_split( ...
    tern(obj.use_scaled_for_goal,'scaled','raw'), obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num, ...
    cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);

% aux içine koy
aux_J1 = J1_split; aux_J2 = J2_split;

% --- Pareto günlüğüne yaz ---
global PARETO;
PARETO.J1(end+1,1)  = aux_J1;
PARETO.J2(end+1,1)  = aux_J2;
PARETO.F(end+1,1)   = J + Penalty;
PARETO.Pen(end+1,1) = Penalty;
PARETO.set(end+1,1) = design_set;
PARETO.x{end+1,1}   = x(:).';
PARETO.feas(end+1,1)= (Penalty <= 1e-6);    % eşik: cezasız ≈ fizibıl

        f = J + Penalty;
        aux = struct('J',J,'Penalty',Penalty,'cons',cons_detail);

        if use_cache, memo(key) = struct('f',f,'aux',aux); end
    catch ME
        % Güvenli büyük ceza
        f = 1e9; aux = struct('err',ME.message);
    end
end

function cfg2 = cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg)
    % simulate() tasarımı içerde decode ediyor; burada yalnız cfg’yi aynen geçiriyoruz.
    %#ok<INUSD>
    cfg2 = ensure_cfg_defaults(cfg);
end
function [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts)
% GA için sürüm-uyumlu çağrı (6/4/2 çıktı destekler)
    nvars = numel(lb);
    try
        % Yeni sürümler (6 çıktı)
        [xbest, fbest, exitflag, output, pop, scores] = ...
            ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
    catch
        try
            % Orta sürümler (4 çıktı)
            [xbest, fbest, exitflag, output] = ...
                ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
            pop = []; scores = [];
        catch
            % Eski sürümler (2 çıktı)
            [xbest, fbest] = ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
            exitflag = []; output = struct(); pop = []; scores = [];
        end
    end
end
function v = aux_if(field)
    v = NaN;
    try, v = output.bestfval; catch, end %#ok<*CTCH>
end
function v = cons_detail_if(which)
    v = NaN;
    try
        % compute from last simulate if elinizde varsa; basit placeholder
        if strcmpi(which,'E_ratio'), v = NaN; end
        if strcmpi(which,'cav95'),   v = NaN; end
    catch
    end
end

% =============================== Alt Yapı ===============================
function f = compact_log_wrapper(x, inner_fitfun)
% Tek satır log: [idx  #feval  J  (J+Penalty)  nViol]
% inner_fitfun: [f, aux] döndürebilir (f=J+Penalty), aux.J ve aux.cons.ratios içerebilir.
    persistent ROW FEVAL
    if isempty(ROW),   ROW   = 0; end
    if isempty(FEVAL), FEVAL = 0; end
    ROW   = ROW + 1;
    FEVAL = FEVAL + 1;

    J = NaN; nviol = 0;

    try
        [f_val, aux] = inner_fitfun(x);   % iki çıktı destekli
        if isstruct(aux)
            if isfield(aux,'J'), J = aux.J; end
            if isfield(aux,'cons') && isfield(aux.cons,'ratios') && isstruct(aux.cons.ratios)
                fn = fieldnames(aux.cons.ratios);
                vals = zeros(numel(fn),1);
                for i=1:numel(fn), vals(i) = aux.cons.ratios.(fn{i}); end
                nviol = sum(vals > 1+1e-12);
            end
        end
        f = f_val;
    catch
        % inner_fitfun tek çıktı verirse
        try
            f = inner_fitfun(x);
        catch
            f = 1e9;   % güvenli büyük ceza
        end
    end

    fprintf('%6d %15d %13.3f %13.3f %8d\n', ROW, FEVAL, J, f, nviol);
end



% ---- PF t_on guard helper ---------------------------------------------
function cfg2 = set_pf_ton_if_nan(cfg2, t5, dt_on)
    if nargin<3 || isempty(dt_on), dt_on = 0.5; end
    cfg2 = ensure_cfg_defaults(cfg2);
    if ~isfield(cfg2,'PF') || ~isfield(cfg2.PF,'t_on') || isnan(cfg2.PF.t_on)
        cfg2.PF.t_on = t5 + dt_on;
    end
end

function [t2,a2] = regrid_to_target(t1,a1,prep)
    % Tekil zaman düzelt, hedef dt'ye göre (auto/off/force)
    [t1,iu]=unique(t1,'stable'); a1=a1(iu);
    dt1 = median(diff(t1),'omitnan');
    tgt = prep.target_dt; tol = max(prep.tol_rel*max(tgt,eps), 1e-12);
    switch lower(prep.resample_mode)
        case 'off'
            t2=t1; a2=a1;
        case 'force'
            t2 = (t1(1):tgt:t1(end)).';
            a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
        otherwise % 'auto'
            if abs(dt1 - tgt) <= tol
                t2=t1; a2=a1;
            else
                t2 = (t1(1):tgt:t1(end)).';
                a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
                warning('Resample: dt=%.6g→%.6g s | N=%d', dt1, tgt, numel(t2));
            end
    end
end

function [t,a,t5,t95] = pick_series(src, rid, dir, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl)

    src = lower(string(src)); dir = upper(string(dir));
    switch src
        case "raw"
            if dir=="X", t=t_rawX{rid}; a=a_rawX{rid}; t5=t5x_raw(rid); t95=t95x_raw(rid);
            else,         t=t_rawY{rid}; a=a_rawY{rid}; t5=t5y_raw(rid); t95=t95y_raw(rid); end
        otherwise % 'scaled'
            if dir=="X", t=t_sclX{rid}; a=a_sclX{rid}; t5=t5x_scl(rid); t95=t95x_scl(rid);
            else,         t=t_sclY{rid}; a=a_sclY{rid}; t5=t5y_scl(rid); t95=t95y_scl(rid); end
    end
    if isempty(a), error('Kayıt #%d için %s yönü mevcut değil.', rid, dir); end
end

function [M,K,C] = make_KCM(n,mv,kv,cv)
    M = diag(mv); K = zeros(n); C = zeros(n);
    for i=1:n
        kL=kv(i); cL=cv(i);
        if i<n, kU=kv(i+1); cU=cv(i+1); else, kU=0; cU=0; end
        K(i,i) = kL + (i<n)*kU;   C(i,i) = cL + (i<n)*cU;
        if i>1, K(i,i-1)=-kL; C(i,i-1)=-cL; end
        if i<n, K(i,i+1)=-kU; C(i,i+1)=-cU; end
    end
end

function [x,a_rel] = lin_MCK_consistent(t, ag, M, C, K)
    n  = size(M,1); r  = ones(n,1);
    dt = median(diff(t));
    agf = griddedInterpolant(t,ag,'linear','nearest');
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',2e-3,'AbsTol',1e-6,'MaxStep',max(dt*10,2e-3),'InitialStep',max(dt*0.25,1e-3));
    sol = ode23tb(odef,[t(1) t(end)],z0,opts);
    t_end = sol.x(end); idx = find(t <= t_end + 1e-12);
    if isempty(idx), x=nan(numel(t),n); a_rel=x; warning('lin_MCK_consistent: early stop'); return; end
    t_use = t(idx); Z = deval(sol,t_use).';
    x_use = Z(:,1:n); v_use = Z(:,n+1:end);
    a_use = ( -(M\(C*v_use.' + K*x_use.')).' - ag(1:numel(t_use)).*r.' );
    x=nan(numel(t),n); a_rel=x; x(1:numel(t_use),:)=x_use; a_rel(1:numel(t_use),:)=a_use;
end


% --- NOTE: Artık 4. opsiyonel çıktı "v" döndürülebilir (mevcut çağrılar çalışmaya devam eder)
function [x,a,diag,varargout] = mck_with_damper_adv(t,ag,M,C,K,k_sd,geom,orf,hyd,therm,num,cfg)
    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, round(hyd.n_parallel)); end
    Ao1 = max(orf.n_orf * (pi*geom.d_o^2/4), 1e-12);
    Ao  = nd * Ao1;              % toplam orifis alanı
    Ap1 = geom.Ap;
    Ap  = nd * Ap1;              % toplam piston alanı

    n = size(M,1); r = ones(n,1); Ns = n-1;
    z0 = zeros(2*n + 2*Ns + Ns + 2, 1);
    z0(2*n + (1:Ns)) = orf.p_amb;  z0(2*n + Ns + (1:Ns)) = orf.p_amb;
    z0(end-1)=therm.T0_C; z0(end)=therm.Ts0_C;

    % --- Ölçekli toleranslar (basınca ve akışa göre) ---
    rho0  = therm.rho_ref;
    dPcap = max(num.dP_cap, 1e5);
    Qcap_est = Ao * sqrt( 2*dPcap / max(rho0,100) );   % Ao = nd*Ao1

    AbsX = 1e-4;  AbsV = 1e-3;  AbsP = max(5e3, 1e-5 * dPcap);
    AbsQ = max(1e-5, 1e-3 * Qcap_est);  AbsT = 5e-3;
    AbsTol = [AbsX*ones(n,1); AbsV*ones(n,1); AbsP*ones(2*Ns,1); AbsQ*ones(Ns,1); AbsT*ones(2,1)];

    dt  = median(diff(t));
    opts = odeset('RelTol',7e-2,'AbsTol',AbsTol,'JPattern',local_JPattern(n), ...
                  'MaxStep',max(dt*15,3e-3),'InitialStep',max(dt*0.30,1.5e-3));

    agf = griddedInterpolant(t, ag, 'pchip', 'nearest');

    odef=@(tt,z) mck_rhs(tt,z,n,Ns,nd,hyd,geom,therm,num,cfg,Ap,Ao,k_sd,M,C,K,agf,r,orf);
    sol=ode23tb(odef,[t(1) t(end)],z0,opts);

    t_end_sol=sol.x(end); t_use=t(t<=t_end_sol+1e-12); Z_use=deval(sol,t_use).';
    if t_end_sol < t(end)-1e-9
        opts2 = odeset(opts, 'RelTol',8e-2, 'MaxOrder', 2);
        try
            sol2=ode15s(odef,[t_end_sol t(end)],sol.y(:,end),opts2);
            t_use2=t(t>t_end_sol); Z_use2=deval(sol2,t_use2).'; Z=[Z_use;Z_use2];
            warning('ode23tb early stop at t=%.3f → continued to t=%.3f.',t_end_sol,t(end));
        catch
            Z=[Z_use; nan(numel(t)-numel(t_use), size(Z_use,2))];
            warning('Both solvers stopped early at t=%.3f s.',t_end_sol);
        end
    else
        Z=Z_use;
    end

    x = Z(:,1:n); v = Z(:,n+1:2*n);
    Fdev = dev_force_from_story(t, x, Z, k_sd, geom, hyd, therm, num, orf, cfg);
    a = ( -(M\(C*v.' + K*x.' + Fdev)).' - ag.*r.' );

    [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
     cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
     cav_margin_t, cav_margin_min] = ...
        diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg);

    diag = struct('drift',drift,'F_story',F_story,'dP_orf_time',dP_orf_t, ...
                  'dP_orf_time_max',max(dP_orf_t,[],2),'T_oil',T_oil,'T_steel',T_steel, ...
                  'mu_t',mu_t,'E_cum',E_cum,'cav_frac_t',cav_frac_t,'dP_q50',dP_q50,'dP_q95',dP_q95, ...
                  'Q_abs_med',Q_abs_med,'Q_abs_p95',Q_abs_p95, ...
                  'cav_margin_t',cav_margin_t,'cav_margin_min',cav_margin_min);

    if nargout>=4, varargout{1} = v; end

end

function dz = mck_rhs(tt,z,n,Ns,nd,hyd,geom,therm,num,cfg,Ap,Ao,k_sd,M,C,K,agf,r,orf)
    x  = z(1:n); v = z(n+1:2*n);
    p1 = z(2*n + (1:Ns)); p2 = z(2*n + Ns + (1:Ns));
    Q  = z(2*n + 2*Ns + (1:Ns)); T_o = z(end-1); T_s = z(end);

    mu_raw = therm.mu_ref * exp(therm.b_mu*(T_o - therm.T_ref_C));
    if cfg.use_thermal && cfg.on.mu_floor
        mu = max(num.mu_min_phys,mu_raw);
    else
        mu = cfg.use_thermal*mu_raw + (~cfg.use_thermal)*therm.mu_ref;
    end
    if cfg.use_thermal
        rho = max(100,therm.rho_ref/(1+therm.alpha_rho*(T_o-therm.T_ref_C)));
        beta = max(1e8,therm.beta0*exp(therm.b_beta*(T_o-therm.T_ref_C)));
        p_vap = p_vap_Antoine(T_o,therm,orf);
    else
        rho = therm.rho_ref; beta = therm.beta0; p_vap = p_vap_Antoine(therm.T_ref_C,therm,orf);
    end

    drift = x(2:end) - x(1:end-1);
    dvel  = v(2:end) - v(1:end-1);

    V1 = nd*hyd.V0 + Ap*drift;
    V2 = nd*hyd.V0 - Ap*drift;
    Vmin = nd * hyd.Vmin_fac * hyd.V0; V1=max(V1,Vmin); V2=max(V2,Vmin);

    % --- Orifis hidrolik kayıpları ---
    if cfg.use_orifice
        if cfg.on.CdRe
            Re = rho .* abs(Q) .* geom.d_o ./ max(Ao * mu, 1e-9);
            Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            Cd = max(min(Cd, 1.2), 0.2);
        else
            Cd = orf.CdInf;
        end
        if cfg.on.Rkv
            RQ = rho ./ max(2 * (Cd .* Ao).^2, 1e-12);
        else
            RQ = 0 * Q;
        end

        Qcap = getfield_default(num,'Qcap_big', ...
            0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
        if cfg.on.Qsat
            Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9));
        else
            Q_h = Q;
        end

        dP_kv = RQ .* Q_h .* abs(Q_h);
        if cfg.on.Rlam
            R_lam  = (128 * mu * geom.Lori / (pi * geom.d_o^4)) / max(1, nd);
            dP_lam = R_lam .* Q;
        else
            dP_lam = 0 * Q;
        end
        dP_h = dP_lam + dP_kv;
    else
        dP_h = 0 * Q; Q_h = 0 * Q;
    end

    % --- Kavitasyon, dP_cap, atalet, leak ---
    if cfg.on.cavitation
        p2_eff = max(p2, orf.cav_sf * p_vap);
    else
        p2_eff = p2;
    end
    dP_raw = p1 - p2_eff;
    if cfg.on.dP_cap && isfinite(num.dP_cap)
        dP_eff = num.dP_cap * tanh(dP_raw ./ max(num.dP_cap, 1));
    else
        dP_eff = dP_raw;
    end
    if cfg.use_orifice && cfg.on.hyd_inertia
        dQ = (dP_eff - dP_h) ./ max(hyd.Lh, 1e-12);
    else
        dQ = 0 * Q;
    end
    if cfg.on.leak
        Q_leak = hyd.K_leak * (p1 - p2);
    else
        Q_leak = 0 * Q;
    end

    % --- Basınç ODE'leri ---
    if cfg.on.pressure_ode
        dp1 = (beta ./ V1) .* ( -Q - Q_leak - Ap * dvel );
        dp2 = (beta ./ V2) .* ( +Q + Q_leak + Ap * dvel );
    else
        dp1 = 0 * p1; dp2 = 0 * p2;
    end
    % Cavitation clamp
    if cfg.on.cavitation
        m1 = (p1 <= p_vap) & ((-Q - Q_leak - Ap*dvel) < 0);
        m2 = (p2 <= p_vap) & ((+Q + Q_leak + Ap*dvel) < 0);
        dp1(m1) = 0; dp2(m2) = 0;
    end

    % --- Damper kuvveti (yay + PF)
    w_pf    = pf_weight_scalar(tt, cfg) * cfg.PF.gain;
    F_story = k_sd * drift;
    if cfg.on.pressure_force
        dp_pf = (p1 - p2_eff);
        if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel);
            dp_pf = s .* max(0, s .* dp_pf);
        end
        F_story = F_story + Ap * (w_pf .* dp_pf);
    end

    % --- Kuvvet dağıtımı, yapı ODE ---
    F        = zeros(n,1);
    F(1)     = -F_story(1);
    if n > 2, F(2:n-1) = F_story(1:end-1) - F_story(2:end); end
    F(n)     =  F_story(end);

    dv = M \ ( -C*v - K*x - F - M*r*agf(tt) );

    % --- Isı ODE'leri ---
    if cfg.use_thermal
        if cfg.use_orifice && cfg.on.Rlam, P_lam = dP_lam .* Q; else, P_lam = 0; end
        if cfg.use_orifice && cfg.on.Rkv,  P_kv  = dP_kv  .* Q; else, P_kv  = 0; end
        P_sum = sum(P_lam + P_kv); P_sum = max(P_sum, 0);

        dT_o = ( P_sum ...
                 - therm.hA_os   * (T_o - T_s) ...
                 - therm.hA_o_env* (T_o - therm.T_env_C) ) / max(therm.C_oil,   eps);
        dT_s = ( + therm.hA_os   * (T_o - T_s) ...
                 - therm.hA_s_env* (T_s - therm.T_env_C) ) / max(therm.C_steel, eps);
    else
        dT_o = 0; dT_s = 0;
    end

    dz = [ v; dv; dp1; dp2; dQ; dT_o; dT_s ];
end

function Jp = local_JPattern(nn)
    Ns   = nn - 1;
    Ntot = 2*nn + 2*Ns + Ns + 2;
    Jp   = sparse(ones(Ntot, Ntot));
end

function [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
          cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
          cav_margin_t, cav_margin_min] = ...
          diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg)

    nd  = 1; if isfield(hyd,'n_parallel'), nd = hyd.n_parallel; end
    Ao1 = orf.n_orf * (pi*geom.d_o^2/4);
    Ao  = nd * Ao1;
    Ap  = nd * geom.Ap;

    X  = Z(:,1:n);
    V  = Z(:,n+1:2*n);
    p1 = Z(:,2*n + (1:Ns));
    p2 = Z(:,2*n + Ns + (1:Ns));
    Q  = Z(:,2*n + 2*Ns + (1:Ns));
    T_oil   = Z(:,end-1);
    T_steel = Z(:,end);

    drift = X(:,2:end) - X(:,1:end-1);
    dvel  = V(:,2:end) - V(:,1:end-1);      % Nt x Ns

    mu_raw  = therm.mu_ref * exp(therm.b_mu*(T_oil - therm.T_ref_C));
    if cfg.use_thermal && cfg.on.mu_floor, mu_t = max(num.mu_min_phys, mu_raw);
    else,                                   mu_t = cfg.use_thermal.*mu_raw + (~cfg.use_thermal).*therm.mu_ref; end

    rho_t   = max(100, therm.rho_ref ./ (1 + therm.alpha_rho*(T_oil - therm.T_ref_C)));
    p_vap_t = p_vap_Antoine(T_oil, therm, orf);

    cav_margin_t  = min(p2 - p_vap_t, [], 2);
    cav_margin_min = min(cav_margin_t, [], 'omitnan');

    if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap_t); else, p2_eff = p2; end

    % ---- PF kuvveti: RESISTIVE-ONLY burada uygulanıyor ----
    w_pf_vec = pf_weight_vec(t, cfg) * cfg.PF.gain;
    if cfg.on.pressure_force
        dp_pf = (p1 - p2_eff);                         % Nt x Ns
        if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuşatma
                            % Nt x Ns
            % İstersen: s = tanh(20*dvel);
            dp_pf = s .* max(0, s .* dp_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;  % Nt x Ns → (.*) yayımlı
    else
        F_story = k_sd*drift;
    end

    % ---- Hidrolik kayıplar (dP_h) + dP_orf_t ----
    if cfg.use_orifice
        if cfg.on.CdRe
            Re = rho_t .* abs(Q) .* geom.d_o ./ max(Ao .* mu_t, 1e-9);
            Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            Cd = max(min(Cd, 1.2), 0.2);
        else
            Cd = orf.CdInf;
        end
        if cfg.on.Rkv, RQ = rho_t ./ max(2 * (Cd .* Ao).^2, 1e-12); else, RQ = 0 * Q; end

        Qcap = getfield_default(num,'Qcap_big', ...
            0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
        if cfg.on.Qsat, Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9)); else, Q_h = Q; end

        dP_kv = RQ .* Q_h .* abs(Q_h);
        if cfg.on.Rlam
            R_lam_t = (128 * mu_t .* geom.Lori ./ (pi * geom.d_o^4)) / max(1, nd);
            dP_lam  = R_lam_t .* Q;
        else
            dP_lam  = 0 * Q;
        end
        dP_h = dP_lam + dP_kv;
    else
        Q_h = 0 * Q; dP_h = 0 * Q;
    end

    dP_raw = p1 - p2_eff;
    if cfg.on.dP_cap && isfinite(num.dP_cap)
        dP_eff = num.dP_cap * tanh( dP_raw ./ max(num.dP_cap, 1) );
    else
        dP_eff = dP_raw;
    end
    epsm     = max(1e3, double(num.softmin_eps));
    dP_orf_t = 0.5 * ( dP_eff + dP_h - sqrt( (dP_eff - dP_h).^2 + epsm^2 ) );

    cav_mask   = p2 < p_vap_t;
    cav_frac_t = mean(cav_mask, 2);

    if cfg.use_orifice
        P_lam = dP_lam .* Q; P_kv  = dP_kv  .* Q;
        P_sum = sum(P_lam + P_kv, 2); P_sum = max(P_sum, 0);
    else
        P_sum = zeros(size(t));
    end
    E_cum = cumtrapz(t, P_sum);

    dP_q50    = prctile(dP_orf_t, 50, 2);
    dP_q95    = prctile(dP_orf_t, 95, 2);
    qQ        = prctile(abs(Q(:)), [50 95]);
    Q_abs_med = qQ(1);
    Q_abs_p95 = qQ(2);
end

function F = dev_force_from_story(t, X, Z, k_sd, geom, hyd, therm, num, orf, cfg)
    Nt=size(X,1); n=size(X,2); Ns=n-1;
    drift = X(:,2:end) - X(:,1:end-1);
    V     = Z(:,n+1:2*n);
    dvel  = V(:,2:end) - V(:,1:end-1);    % Nt x Ns
    p1    = Z(:,2*n+(1:Ns)); 
    p2    = Z(:,2*n+Ns+(1:Ns)); 
    T_o   = Z(:,end-1);

    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, hyd.n_parallel); end
    Ap = nd * geom.Ap;

    if cfg.use_thermal, p_vap=p_vap_Antoine(T_o,therm,orf); else, p_vap=p_vap_Antoine(therm.T_ref_C,therm,orf); end
    if cfg.on.cavitation, p2_eff=max(p2,orf.cav_sf*p_vap); else, p2_eff=p2; end

    w_pf_vec = pf_weight_vec(t,cfg)*cfg.PF.gain;

    if cfg.on.pressure_force
        dp_pf = (p1 - p2_eff);                  % Nt x Ns
        if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuşatma
            % yumuşak istersen: s = tanh(20*dvel);
            dp_pf = s .* max(0, s .* dp_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;
    else
        F_story = k_sd*drift;
    end

    F = zeros(n,Nt);
    F(1,:) = -F_story(:,1).';
    if n>2, F(2:n-1,:) = (F_story(:,1:end-1) - F_story(:,2:end)).'; end
    F(n,:) =  F_story(:,end).';
end


% ---- PF ramp ağırlığı -------------------------------------------------
function w = pf_weight_scalar(tt,cfg)
    if ~cfg.on.pressure_force, w=0; return; end
    switch lower(cfg.PF.mode)
        case 'off', w=0;
        case 'on',  w=1;
        otherwise
            dt=max(tt-cfg.PF.t_on,0); w=1-exp(-dt/max(cfg.PF.tau,1e-6));
    end
    w=min(max(w,0),1);
end
function wv = pf_weight_vec(t,cfg)
    if ~cfg.on.pressure_force, wv=zeros(size(t)); return; end
    switch lower(cfg.PF.mode)
        case 'off', wv=zeros(size(t));
        case 'on',  wv=ones(size(t));
        otherwise
            dt=max(t-cfg.PF.t_on,0); wv=1-exp(-dt/max(cfg.PF.tau,1e-6));
    end
    wv=min(max(wv,0),1);
end

% ---- Buhar basıncı (Antoine) -----------------------------------------
function p_v = p_vap_Antoine(T_C, therm, ~)
    if isfield(therm,'antoine_A') && isfield(therm,'antoine_B') && isfield(therm,'antoine_C')
        A=therm.antoine_A; B=therm.antoine_B; C=therm.antoine_C;
    else, A=5.0; B=1700; C=-80; end
    T_C=double(T_C); p_v = 10.^(A - B./(C + T_C));
p_v = min(max(p_v, 5), 5e2);     % 5–500 Pa

end

% ---- Arias penceresi --------------------------------------------------
function [t5,t95] = arias_win(t,ag,p1,p2)
    if nargin<3||isempty(p1), p1=0.05; end
    if nargin<4||isempty(p2), p2=0.95; end
    t=t(:); ag=ag(:); mask=isfinite(t)&isfinite(ag); t=t(mask); ag=ag(mask);
    if numel(t)<2, t5=t(1); t95=t(end); return; end
    [t,ia]=unique(t,'stable'); ag=ag(ia); dt=diff(t); a2=ag.^2;
    Eincr=0.5*(a2(1:end-1)+a2(2:end)).*dt; E=[0; cumsum(Eincr)]; Eend=E(end);
    if Eend<=eps, t5=t(1); t95=t(end); return; end
    t5=arias_taf(0.05, E, Eend, t); t95=arias_taf(0.95, E, Eend, t); t5=max(min(t5,t(end)),t(1)); t95=max(min(t95,t(end)),t(1));
    if t95<=t5, t5=t(1); t95=t(end); end
end

function tv = arias_taf(frac, E, Eend, t)
    target=frac*Eend; k=find(E>=target,1,'first');
    if isempty(k), tv=t(end); elseif k==1, tv=t(1);
    else
        E0=E(k-1); E1=E(k);
        r=(target-E0)/max(E1-E0,eps); r=min(max(r,0),1);
        tv=t(k-1)+r*(t(k)-t(k-1));
    end
end

% --------- PSA (klasik) -----------------------------------------------
function Sab = sdof_PSA_band_avg_ode(t, ag, T1, zeta, band_fac, Np)
    Tvec = linspace(band_fac(1)*T1, band_fac(2)*T1, Np);
    Sa   = sdof_PSA_vec_ode(t, ag, Tvec, zeta);
    Sab  = mean(Sa);
end

function Sa_vec = sdof_PSA_vec_ode(t, ag, T_vec, zeta)
    T_vec = T_vec(:);
    Sa_vec = zeros(numel(T_vec),1);
    for i=1:numel(T_vec)
        Sa_vec(i) = sdof_PSA_ode(t, ag, T_vec(i), zeta);
    end
end

% --------- PSA (augmented — hızlı) ------------------------------------
function Sab = sdof_PSA_band_avg_aug(t, ag, T1, zeta, band_fac, Np)
    Tvec = linspace(band_fac(1)*T1, band_fac(2)*T1, Np);
    Sa   = sdof_PSA_vec_aug_ode(t, ag, Tvec, zeta);
    Sab  = mean(Sa);
end
function Sa_vec = sdof_PSA_vec_aug_ode(t, ag, T_vec, zeta)
    T_vec = T_vec(:); Np = numel(T_vec);
    wn = 2*pi./max(T_vec,eps);
    agf = griddedInterpolant(t, ag, 'linear','nearest');

    z0 = zeros(2*Np,1);
    dt = median(diff(t));
    opts = odeset('RelTol',1e-4,'AbsTol',1e-6,...
                  'MaxStep',max(10*dt,1e-3),'InitialStep',max(0.25*dt,1e-4));
    odef = @(tt,zz) aug_rhs(tt,zz,wn,zeta,agf);
    sol  = ode23tb(odef, [t(1) t(end)], z0, opts);

    Z = deval(sol, t).';                 % Nt x (2*Np)
    X = Z(:,1:2:end);                    % Nt x Np
    V = Z(:,2:2:end);                    % Nt x Np
    Nt = size(Z,1);
    Aabs = abs( -2*zeta*(ones(Nt,1)*wn.').*V ...
                - (ones(Nt,1)*(wn.'.^2)).*X ...
                - ag(:)*ones(1,Np) );    % Nt x Np
    Sa_vec = max(Aabs,[],1).';
end
function dz = aug_rhs(tt,zz,wn,zeta,agf)
    x = zz(1:2:2*numel(wn)); v = zz(2:2:2*numel(wn));
    agt = agf(tt);
    ax  = v;
    av  = -2*zeta*wn.*v - (wn.^2).*x - agt;
    dz  = zeros(2*length(wn),1);
    dz(1:2:end) = ax;
    dz(2:2:end) = av;
end

% --------- PSA yardımcı (yalnız PSA için downsample) -------------------
function [t2,a2] = psa_grid(t,a,dt_target)
    if isempty(dt_target) || dt_target<=0
        t2=t; a2=a; return;
    end
    dt = median(diff(t));
    if dt >= dt_target-1e-12
        t2=t; a2=a;
    else
        t2 = (t(1):dt_target:t(end)).';
        a2 = pchip(t,a,t2);
    end
end

% ===================== GA Setleri & Decode Yardımcıları =================
function [lb,ub,int_idx,names] = ga_get_bounds(set_id)
    switch set_id
        case 1
           %      d_o   Lori   mu_ref   Kd     d_w    D_m   n_turn    Dp    Lgap
            lb = [2e-3, 0.02, 1, 1.5e9, 3e-3, 0.06,   20, 0.10, 0.10];
            ub = [6e-3, 0.05, 1.50, 2.0e9,15e-3, 0.10,  40, 0.14, 0.24];
            int_idx = 7; % n_turn
            names = {'d_o','Lori','mu_ref','Kd','d_w','D_m','n_turn','Dp','Lgap'};

        case 2
            % [n_orf, Cd0, CdInf, Rec, p_exp, cav_sf, Lh, K_leak, resFactor]
            lb = [2, 0.55, 0.90, 2000, 1.0, 0.95, 1e-3, 1e-9, 10];
            ub = [2, 0.70, 1.10, 6000, 1.6, 1.00, 8e-3, 6e-9, 18];
            int_idx = 1; % n_orf
            names = {'n_orf','Cd0','CdInf','Rec','p_exp','cav_sf','Lh','K_leak','resFactor'};

                otherwise % 3  (n_orf SABİT=2)
            % [n_orf, Cd0, mu_ref, b_mu, beta0, b_beta, hA_os, dP_cap, Vmin_fac]
            lb = [2, 0.64, 0.6, -0.012, 1.2e9, -6e-3, 450, 3e8, 0.90,  10];
            ub = [2, 0.68, 1.2, -0.007, 2.2e9, -2e-3, 800, 6e8, 0.95, 18];
            int_idx = 1; % n_orf (sabit ama interface için sorun değil)
            names = {'n_orf','Cd0','mu_ref','b_mu','beta0','b_beta','hA_os','dP_cap','Vmin_fac'};

    end
end


function [geom, sh, orf, hyd, therm, num, ga] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num)
    % --- set_id'i güvene al (string/char gelirse sayıya çevir) ---
    set_id = ga.design_set;
    if isstring(set_id) || ischar(set_id)
        set_id = str2double(set_id);
    end
    if isempty(set_id) || isnan(set_id)
        set_id = 1;
    end
    ga.design_set = double(set_id);
    % GA kapalıysa hiçbir şeyi değiştirme
    if ~isfield(ga,'enable') || ~ga.enable || isempty(ga.x)
        [lb,ub,ii,nn] = ga_get_bounds(ga.design_set);
        ga.lb=lb; ga.ub=ub; ga.int_idx=ii; ga.names=nn; ga.x_use=[];
        return;
    end

    [lb,ub,int_idx,names] = ga_get_bounds(ga.design_set);
    x = ga.x(:);
    if numel(x)~=numel(lb)
        error('ga.x uzunluğu set-%d için %d olmalı.', ga.design_set, numel(lb));
    end
    % Kenetleme + tamsayı
x = x(:);
lb = lb(:);
ub = ub(:);
x = max(lb, min(ub, x));
    if numel(int_idx)>=1
        idx = int_idx;
        x(idx) = round(x(idx));
        x = max(lb, min(ub, x));
    end

    switch ga.design_set
        case 1
            % [d_o, Lori, mu_ref, Kd, d_w, D_m, n_turn, Dp, Lgap]
            geom.d_o  = x(1);
            geom.Lori = x(2);
            therm.mu_ref = x(3);
            geom.Kd   = x(4);
            sh.d_w    = x(5);
            sh.D_m    = x(6);
            sh.n_turn = x(7);
            geom.Dp   = x(8);
            geom.Lgap = x(9);
        case 2
            % [n_orf, Cd0, CdInf, Rec, p_exp, cav_sf, Lh, K_leak, resFactor]
            orf.n_orf   = x(1);
            orf.Cd0     = x(2);
            orf.CdInf   = max(x(3), x(2)+0.05); % tutarlılık
            orf.Rec     = x(4);
            orf.p_exp   = x(5);
            orf.cav_sf  = x(6);
            hyd.Lh      = x(7);
            hyd.K_leak  = x(8);
            therm.resFactor = x(9);
        otherwise % 3
            % [n_orf, Cd0, mu_ref, b_mu, beta0, b_beta, hA_os, dP_cap, Vmin_fac]
            orf.n_orf   = x(1);
            orf.Cd0     = x(2);
            therm.mu_ref= x(3);
            therm.b_mu  = x(4);
            therm.beta0 = x(5);
            therm.b_beta= x(6);
            therm.hA_os = x(7);
            num.dP_cap  = x(8);
            hyd.Vmin_fac= x(9);
therm.resFactor = x(10);
    end

    ga.lb=lb; ga.ub=ub; ga.int_idx=int_idx; ga.names=names; ga.x_use=x;
    % ... (x hazırlandıktan SONRA)
    xdisp = strjoin( compose('%.3g', x(:).'), ', ' );  % 9 sayı -> "a, b, c..."
 % (Sessiz mod) — yalnızca LOG.verbose_decode=true ise yaz
try
    if evalin('base','exist(''LOG'',''var'')') && evalin('base','isstruct(LOG)') ...
            && evalin('base','isfield(LOG,''verbose_decode'') && LOG.verbose_decode')
        fprintf('GA decode: set=%s | x_use = [%s]\n', num2str(ga.design_set), xdisp);
    end
catch
    % base workspace erişimi yoksa sessizce geç
end
end

% ===================== İNCE SARMA – simulate(...) ======================
% Standart, güvenli tek çağrı noktası (GA ve çoklu denemeler için)
function resp = simulate(design_set, x, mu_mult, t, ag, ...
    M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg)
design_set = double(design_set);
    if ~isfinite(design_set) || ~ismember(design_set,[1 2 3])
        design_set = 1;
    end
% Standart, güvenli tek çağrı noktası (GA ve çoklu denemeler için)

    resp = struct('ok',false,'msg','','y',struct(),'drift',[],'dvel',[], ...
                  'dP_orf_env',[],'c_lam',NaN,'dT_est',NaN,'metrics',struct());
    try
        % ---- guard + opsiyonel viskozite çarpanı
        cfg = ensure_cfg_defaults(cfg);
        if ~isempty(mu_mult) && isfinite(mu_mult) && mu_mult>0
            therm.mu_ref = therm.mu_ref * mu_mult;
        end

        % ---- GA decode (varsa)
        ga_local = struct('enable',~isempty(x),'design_set',design_set,'x',x);
        [geom, sh, orf, hyd, therm, num, ~] = decode_design_apply(ga_local, geom, sh, orf, hyd, therm, num);

        % ---- Türetilenler
        geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);
nd = max(1, getfield_default(hyd,'n_parallel',1));
geom.Ap_eff = nd * geom.Ap;         % toplam piston alanı
orf.Ao_eff  = nd * orf.Ao;          % toplam orifis alanı
hyd.n_parallel = nd;                % fonksiyonlara iletmek için
        k_p       = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
        k_h       = geom.Kd*geom.Ap^2/geom.Lgap;
        k_s       = geom.Ebody*geom.Ap/geom.Lgap;
        k_hyd     = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
k_sd = nd * (k_hyd + k_p);

        % hacim ve ısı kapasiteleri
    nStories = n-1;
steel_to_oil_mass_ratio = 1.5;

hyd.V0 = 0.5*(geom.Ap*(2*geom.Lgap));                 % tek damper/tek hazne
V_oil_per = therm.resFactor*(geom.Ap*(2*geom.Lgap));  % tek damper başına yağ
nDtot     = nStories * nd;

m_oil_tot    = nDtot*(therm.rho_ref*V_oil_per);
m_steel_tot  = steel_to_oil_mass_ratio*m_oil_tot;
therm.C_oil   = m_oil_tot*therm.cp_oil;
therm.C_steel = m_steel_tot*therm.cp_steel;

    
       % ---- PF guard: t5 + 0.5
[t5_sim,~] = arias_win(t, ag, 0.05, 0.95);
cfg = set_pf_ton_if_nan(cfg, t5_sim, 0.5);   % <-- BURASI: t5_plot → t5_sim

        % ---- Çözüm (hız isteğe bağlı çıktıyla)
        [xD, aD, diag, vD] = mck_with_damper_adv(t, ag, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);
% ---- NaN/Inf emniyeti (özellikle dP_orf_time_max için) ----
if ~isfield(diag,'dP_orf_time_max') || any(~isfinite(diag.dP_orf_time_max))
    if isfield(diag,'dP_orf_time') && ~isempty(diag.dP_orf_time)
        diag.dP_orf_time_max = max(diag.dP_orf_time, [], 2, 'omitnan');
    else
        diag.dP_orf_time_max = zeros(numel(t),1); % güvenli fallback
    end
end

% Ana dizilerde tekil NaN/Inf temizliği (fail yerine yumuşat)
fixnf = @(A) (A + 0.*A);          % sınıfı korumak için num trick
xD(~isfinite(xD)) = 0;  aD(~isfinite(aD)) = 0;

        % ---- Standart çıktı paketi
        resp.y = struct('x',xD,'v',vD,'a',aD,'t',t);
        resp.drift = diag.drift;
        resp.dvel  = vD(:,2:end) - vD(:,1:end-1);
        resp.dP_orf_env = diag.dP_orf_time_max;               % kat zarfı (max across stories)
        resp.c_lam = 128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4);
        resp.dT_est = diag.T_oil(end) - diag.T_oil(1);

        % pratik metrikler
        x10_max = max(abs(xD(:,min(10,size(xD,2)))));
        a3_rms  = sqrt(mean(aD(:,min(3,size(aD,2))).^2,'omitnan'));
        cav_p95 = prctile(diag.cav_frac_t,95);
        Toil_max= max(diag.T_oil);
        resp.metrics = struct('x10_max',x10_max,'a3_rms',a3_rms, ...
                      'cav_frac_p95',cav_p95,'Q_abs_p95',diag.Q_abs_p95, ...
                      'T_oil_max',Toil_max);

        resp.diag = diag;   % cav_frac_t dahil tüm zaman serisi diagnostiklerini dışarı ver

        % Zarf kuvvet (hikaye bazında max → zaman serisi)
        F_story_env = max(abs(diag.F_story),[],2);   % Nt x 1
        resp.F_env  = F_story_env;
        resp.F_max  = max(F_story_env);              % skaler

        % Zarf strok (hikayeler üzerinde max)
        drift_env = max(abs(resp.drift),[],2);       % Nt x 1
        resp.stroke_max = max(drift_env);            % skaler

        % Δp_orf için zaman-içi quantile hesap kolaylığı
        resp.dP_q_time = @(qq) prctile(resp.dP_orf_env, 100*qq);
% ---- Sağlamlık kontrolleri
resp.dP_orf_env = diag.dP_orf_time_max;  % (yukarıda garanti edildi)

isBad = @(A) isempty(A) || any(~isfinite(A(:)));  % NaN/Inf veya boş

if isBad(xD) || isBad(aD) || isBad(vD) || isBad(resp.dP_orf_env)
    resp.ok  = false;
    resp.msg = 'NaN/Inf/empty tespit edildi';
    return;
end

% dP_orf_env zaten üstte garanti edildi; x/a boyut kontrolü:
if numel(t) ~= size(xD,1) || numel(t) ~= size(aD,1)
    resp.ok  = false;
    resp.msg = 'Zaman uzunluğu uyuşmazlığı';
    return;
end


        if numel(t)~=size(xD,1) || numel(t)~=size(aD,1)
            resp.ok=false; resp.msg='Zaman uzunluğu uyuşmazlığı'; return;
        end

        resp.ok=true; resp.msg='ok';

    catch ME
        % ---- Hata halinde NaN dolgulu güvenli dönüş
        N = numel(t); ntry = n;
        resp.ok=false; resp.msg=sprintf('simulate: %s',ME.message);
        resp.y = struct('x',nan(N,ntry),'v',nan(N,ntry),'a',nan(N,ntry),'t',t(:));
        resp.drift = nan(N,max(ntry-1,1));
        resp.dvel  = nan(N,max(ntry-1,1));
        resp.dP_orf_env = nan(N,1);
        resp.c_lam = NaN; resp.dT_est = NaN;
        resp.metrics = struct('x10_max',NaN,'a3_rms',NaN,'cav_frac_p95',NaN,'Q_abs_p95',NaN,'T_oil_max',NaN);
    end
end
function v = getfield_default(S, fname, defaultVal)
% Güvenli alan okuma: alan yoksa/boşsa/NaN-Inf ise default döndür.
    if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
        v = defaultVal;
        return;
    end
    val = S.(fname);
    if isnumeric(val)
        if isempty(val) || any(~isfinite(val(:)))
            v = defaultVal;
        else
            v = val;
        end
    else
        % Sayısal olmayan türler için yalnızca boşluk kontrolü yaptık
        v = val;
    end
end



function Jp = local_JPattern(n)
% Güvenli, tutucu (full) JPattern — küçük boyutlarda performans yeterli
% İstersen kendi bant/blok yapını sonra geri koyarsın.
    Ns = n - 1;
    Ntot = 2*n + 2*Ns + Ns + 2;
    Jp = sparse(ones(Ntot, Ntot));  % tüm girişlerin potansiyel nonzero olduğunu varsay
end

%% ===================== EKSİK/YARDIMCI FONKSİYONLAR =====================

function val = getfield_default(S, name, defaultVal)
% S alanı yoksa default döndürür
    if ~isstruct(S) || ~isfield(S,name) || isempty(S.(name))
        val = defaultVal;
    else
        val = S.(name);
    end
end

function [tPSA, agPSA] = psa_grid(t, ag, down_dt)
% PSA için opsiyonel yeniden örnekleme (yalnız downsample; upsample etmez)
    t = t(:); ag = ag(:);
    if nargin<3 || isempty(down_dt) || down_dt<=0
        tPSA = t; agPSA = ag; return;
    end
    dt0 = median(diff(t),'omitnan');
    if down_dt <= dt0*(1+1e-12)
        tPSA = t; agPSA = ag; return;
    end
    tPSA  = (t(1):down_dt:t(end)).';
    agPSA = interp1(t, ag, tPSA, 'pchip', 'extrap');
end



function [J1, J2] = compute_objectives_split( ...
    src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga)

% J1: IDR zarf oranı CVaR; J2: mutlak ivme (RMS+p95) zarf oranı CVaR
    h_story_m = 3.0 * ones(n-1,1); % zaten üstte de var; buraya da koyduk
    [J1, ~] = compute_J1_IDR_over_records( ...
        src, obj, h_story_m, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga);

    [J2, ~] = compute_J2_ACC_over_records( ...
        src, obj, h_story_m, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga);
end

function S = decode_design_apply(ga, geom, sh, orf, hyd, therm, num)
% Şimdilik NO-OP (GA kapalıyken akışı bozmamak için).
% GA'yı kullanacaksan, burada ga.x'i set-id'e göre geom/sh/orf/hyd/therm/num'a uygula.
    S = deal(geom); %#ok<NASGU>
    % ( intentionally empty )
end

function resp = simulate(design_set, x_ga, mu_scale, t, ag, ...
    M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg)

    %#ok<*INUSD>
    ag_mu = mu_scale * ag;

    % k_sd aynı formülle türetilir (ana gövdeyle tutarlı)
    nd      = max(1, getfield_default(hyd,'n_parallel',1));
    k_p     = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
    k_h     = geom.Kd*geom.Ap^2/geom.Lgap;
    k_s     = geom.Ebody*geom.Ap/geom.Lgap;
    k_hyd   = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
    k_sd    = nd * (k_hyd + k_p);

    try
        [xD, aD, dlog, vD] = mck_with_damper_adv( ...
            t, ag_mu, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);

        ok = all(isfinite(xD(:))) && all(isfinite(aD(:)));

        Fmax   = max(abs(dlog.F_story), [], 'all', 'omitnan');
        stroke = max(abs(dlog.drift),   [], 'all', 'omitnan');

        dp_env = dlog.dP_orf_time_max;    % Nt x 1
        dP_q_time = @(q) prctile(dp_env, 100*q);

        % ΔT tahmini: başlangıç yağ sıcaklığından artış
        i0 = find(isfinite(dlog.T_oil),1,'first');
        if isempty(i0), dT_est = 0;
        else, dT_est = max(dlog.T_oil,[],'omitnan') - dlog.T_oil(i0);
        end

        resp = struct( ...
            'ok', ok, ...
            'y', struct('x', xD, 'a', aD, 'v', vD), ...
            'F_max', Fmax, ...
            'stroke_max', stroke, ...
            'dP_q_time', dP_q_time, ...
            'dT_est', dT_est, ...
            'diag', dlog, ...
            'metrics', struct('Q_abs_p95', dlog.Q_abs_p95) );
    catch
        resp = struct('ok', false);
    end
end


% ---------------- PSA (augmented varyantları, klasik fonksiyonlara delege) --
function Sab = sdof_PSA_band_avg_aug(t, ag, T1, zeta, band_fac, Np)
% Hızlı/tek ODE seçeneği yerine klasik fonksiyona delege – güvenli ve yeterli
    Sab = sdof_PSA_band_avg_ode(t, ag, T1, zeta, band_fac, Np);
end

function Sa_vec = sdof_PSA_vec_aug_ode(t, ag, T_vec, zeta)
    Sa_vec = sdof_PSA_vec_ode(t, ag, T_vec, zeta);
end

% ---------------- PSA (klasik tek periyot) – KESİLEN FONKSİYONUN TAMAMI ---
function Sa = sdof_PSA_ode(t, ag, T, zeta)
% Tek SDOF için PSA (mutlak ivme zarfı). ag: yer ivmesi (+yukarı pozitif).
    wn  = 2*pi/max(T,eps);
    agf = griddedInterpolant(t, ag, 'pchip', 'nearest');
    odef = @(tt, z) [ z(2);
                      -2*zeta*wn*z(2) - wn^2*z(1) - agf(tt) ];
    z0 = [0;0];
    % Naif ama sağlam seçenekler:
    dt  = median(diff(t), 'omitnan');
    opts = odeset('RelTol',2e-4,'AbsTol',1e-7,'MaxStep',max(10*dt,2e-3),'InitialStep',max(0.25*dt,1e-3));
    sol = ode23tb(odef, [t(1) t(end)], z0, opts);

    % yoğ. değerlendirme
    tt = t;                         % giriş ızgarasında num. değerlendirelim
    Z  = deval(sol, tt).';
    v  = Z(:,2);
    a_rel = -2*zeta*wn*v - wn^2*Z(:,1) - ag(:);  % ODE’nin 2. satırı + tekrar ag çıkardık → a_rel
    a_abs = a_rel + ag(:);                       % mutlak ivme
    Sa    = max(abs(a_abs));
end
