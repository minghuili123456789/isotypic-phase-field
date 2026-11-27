% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Drop-in extensibility tests:
%  (1) SAV (Shen-Xu-Yang 2018) for 1D Allen-Cahn (unshifted & shifted theta=pi/4)
%  (2) Conservative Allen-Cahn (CAC) for 2D odd-odd X-junction (N=256)
%  (3) Cahn-Hilliard (4th-order) for 2D odd-odd pattern (N=192)
% Verifies that plugging our isotypic backend underneath ANY of these schemes
% yields (a) ~1e-15 equivalence on non-symmetric data and (b) D === 0 bitwise
% on symmetric data, where standard FFT implementations all break symmetry.
here = fileparts(mfilename('fullpath')); addpath(here);
res  = fullfile(here,'..','..','results');

fprintf('=== (1) Drop-in SAV (1st-order & BDF2-SAV) on 1D Allen-Cahn ===\n');
N = 128; tau = 0.01; T = 100; nst = round(T/tau);
C0 = 10; % SAV shift constant so E1(u) + C0 > 0
% Helper: exact discrete inner product <a, b> = (2pi/N)*sum(a.*b)
% In isotypic channels (u = sum_c u_c), orthogonal channels give exact
% inner product using the standard quadrature weights on the fundamental grid!
% (Equivalently, Parseval in spectral space, or channel-wise weighted sum.)

% 1a. Equivalence test of SAV on non-symmetric 1D data
u_ns = 0.4*sin((0:N-1)'*2*pi/N) + 0.3*cos(2*(0:N-1)'*2*pi/N) + 0.15*sin(3*(0:N-1)'*2*pi/N + 0.7);
Pf = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier')); Of = sp_ops(Pf);
Pi = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));     Oi = sp_ops(Pi);
[uf_ns, rf_ns] = run_sav_1d(Pf, Of, u_ns, tau, 50, C0);
[ui_ns, ri_ns] = run_sav_1d(Pi, Oi, sp_channels(u_ns, Pf), tau, 50, C0);
eq_sav = max(abs(uf_ns - sp_assemble(ui_ns, 1, N/2)));
fprintf('  SAV non-symmetric equivalence (50 steps): max|u_std - u_iso| = %.2e, |r_std - r_iso| = %.2e\n', ...
    eq_sav, abs(rf_ns - ri_ns));

% 1b. Symmetry test: 1D u0 = sin(x), kappa=1, T=100
u0_1d = sin((0:N-1)'*2*pi/N);
U0_1d = sp_init(Pi, @(x) sin(x));
u0_1d = sp_assemble(U0_1d, 1, N/2); % exact bitwise symmetric
[uf_s, ~, tb_sav_std, maxD_sav_std] = run_sav_1d_diag(Pf, Of, u0_1d, tau, nst, C0);
[ui_s, ~, tb_sav_iso, maxD_sav_iso] = run_sav_1d_diag(Pi, Oi, U0_1d, tau, nst, C0);
fprintf('  SAV 1D unshifted (T=%g): std T*=%.2f (maxD=%.2e) | iso T*=%g (maxD=%g)\n', ...
    T, tb_sav_std, maxD_sav_std, tb_sav_iso, maxD_sav_iso);

% 1c. Shifted SAV with T2 backend: u0 = sin(x + pi/4), N=128, T=100 (pure double!)
th = pi/4;
Pt2 = sp_t2_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0)); Ot2 = sp_t2_ops(Pt2);
U0_t2 = sp_t2_init(Pt2, @(x) sin(x + th), 1);
u0_sh = sp_t2_assemble(U0_t2, 1, N/2);
[~, ~, tb_sav_sh_std, maxD_sav_sh_std] = run_sav_1d_diag(Pf, Of, u0_sh, tau, nst, C0, th);
[~, ~, tb_sav_sh_t2,  maxD_sav_sh_t2]  = run_sav_1d_diag(Pt2, Ot2, U0_t2, tau, nst, C0, th);
fprintf('  SAV 1D shifted (th=pi/4, T=%g): std T*=%.2f (maxD=%.2e) | T2-iso T*=%g (maxD=%g)\n', ...
    T, tb_sav_sh_std, maxD_sav_sh_std, tb_sav_sh_t2, maxD_sav_sh_t2);

fprintf('\n=== (2) Drop-in Conservative Allen-Cahn (CAC) on 2D X-junction (N=256) ===\n');
N2 = 256; eps2d = 0.05; tau2 = 1e-3; T2 = 1.0; nst2 = round(T2/tau2);
Pf2 = sp_setup(struct('d',2,'N',N2,'eps2',1,'gam',1/eps2d^2,'S',2/eps2d^2,'rep','fourier','model','CAC'));
Of2 = sp_ops(Pf2); If2 = sp_integrator(Pf2, Of2, 'etdrk2', tau2);
Pi2 = sp_setup(struct('d',2,'N',N2,'eps2',1,'gam',1/eps2d^2,'S',2/eps2d^2,'rep','iso','model','CAC'));
Oi2 = sp_ops(Pi2); Ii2 = sp_integrator(Pi2, Oi2, 'etdrk2', tau2);
U0_cac = sp_init(Pi2, @(x,y) tanh(sin(x).*sin(y)/(sqrt(2)*eps2d)));
u0_cac = sp_assemble(U0_cac, 2, N2/2);
uhf = Of2.toSpec(u0_cac); uhi = Oi2.toSpec(U0_cac);
tb_cac_std = NaN; maxD_cac_std = 0; maxD_cac_iso = 0; maxMean_cac = 0;
for n = 1:nst2
    uhf = If2.step(uhf); uhi = Ii2.step(uhi);
    if mod(n, 20) == 0
        uf = Of2.toPhys(uhf); Uc = sp_channels(uf, Pf2);
        Df = max(max(abs(Uc(:,:,[1 2 3])),[],'all'));
        maxD_cac_std = max(maxD_cac_std, Df);
        maxMean_cac  = max(maxMean_cac, abs(mean(uf(:))));
        if isnan(tb_cac_std) && Df > 1e-2, tb_cac_std = n*tau2; end
        Vi = Oi2.toPhys(uhi);
        Di = max(abs(Vi(:,:,[1 2 3])),[],'all');
        maxD_cac_iso = max(maxD_cac_iso, Di);
    end
end
fprintf('  2D CAC X-junction (N=%d, T=%g): std T*=%.3f (maxD=%.2e, max|mean|=%.2e) | iso maxD=%g\n', ...
    N2, T2, tb_cac_std, maxD_cac_std, maxMean_cac, maxD_cac_iso);

fprintf('\n=== (3) Drop-in 2D Cahn-Hilliard (4th-order) on odd-odd X-junction (N=192) ===\n');
% du/dt = -Delta(eps^2 Delta u + u - u^3) = -eps^2 Delta^2 u - Delta(u - u^3)
% Stabilized ETDRK2: L_k = eps^2 |k|^4 + S |k|^2, N_k(u) = |k|^2 * F[ (1+S)u - u^3 ]
Nch = 192; eps_ch = 0.1; tau_ch = 5e-4; T_ch = 1.0; nst_ch = round(T_ch/tau_ch);
S_ch = 2;
Pf_ch = sp_setup(struct('d',2,'N',Nch,'eps2',eps_ch^2,'gam',1,'S',S_ch,'rep','fourier'));
Pi_ch = sp_setup(struct('d',2,'N',Nch,'eps2',eps_ch^2,'gam',1,'S',S_ch,'rep','iso'));
% Override linear operator Lk = eps^2*K2^2 + S*K2 and wrap O.N with K2.*(...)
Pf_ch.Lk = eps_ch^2 * Pf_ch.K2.^2 + S_ch * Pf_ch.K2;
Pi_ch.Lk = eps_ch^2 * Pi_ch.K2.^2 + S_ch * Pi_ch.K2;
Of_ch = sp_ops(Pf_ch); Oi_ch = sp_ops(Pi_ch);
Of_ch.N = @(uh) Pf_ch.K2 .* Of_ch.toSpec((1+S_ch)*Of_ch.toPhys(uh) - Of_ch.toPhys(uh).^3);
Oi_ch.N = @(uh) Pi_ch.K2 .* Oi_ch.toSpec((1+S_ch)*Oi_ch.toPhys(uh) - Oi_ch.cube(Oi_ch.toPhys(uh)));
If_ch = sp_integrator(Pf_ch, Of_ch, 'etdrk2', tau_ch);
Ii_ch = sp_integrator(Pi_ch, Oi_ch, 'etdrk2', tau_ch);

% 3a. Non-symmetric CH equivalence check (N=32, 20 steps, tau = 5e-4)
Pfc_t = sp_setup(struct('d',2,'N',32,'eps2',eps_ch^2,'gam',1,'S',S_ch,'rep','fourier'));
Pic_t = sp_setup(struct('d',2,'N',32,'eps2',eps_ch^2,'gam',1,'S',S_ch,'rep','iso'));
Pfc_t.Lk = eps_ch^2*Pfc_t.K2.^2 + S_ch*Pfc_t.K2; Pic_t.Lk = eps_ch^2*Pic_t.K2.^2 + S_ch*Pic_t.K2;
Ofc_t = sp_ops(Pfc_t); Oic_t = sp_ops(Pic_t);
Ofc_t.N = @(uh) Pfc_t.K2 .* Ofc_t.toSpec((1+S_ch)*Ofc_t.toPhys(uh) - Ofc_t.toPhys(uh).^3);
Oic_t.N = @(uh) Pic_t.K2 .* Oic_t.toSpec((1+S_ch)*Oic_t.toPhys(uh) - Oic_t.cube(Oic_t.toPhys(uh)));
Ifc_t = sp_integrator(Pfc_t, Ofc_t, 'etdrk2', tau_ch);
Iic_t = sp_integrator(Pic_t, Oic_t, 'etdrk2', tau_ch);
rng(7); u_ch_ns = 0.3*randn(32,32);
uh1 = Ofc_t.toSpec(u_ch_ns); uh2 = Oic_t.toSpec(sp_channels(u_ch_ns, Pfc_t));
for n = 1:20, uh1 = Ifc_t.step(uh1); uh2 = Iic_t.step(uh2); end
eq_ch = max(abs(Ofc_t.toPhys(uh1) - sp_assemble(Oic_t.toPhys(uh2), 2, 16)),[],'all');
fprintf('  2D CH non-symmetric equivalence (N=32, 20 steps, tau=%.1e): %.2e\n', tau_ch, eq_ch);

% 3b. Symmetry test on 2D odd-odd X-junction (N=192, eps=0.1, tau=5e-4, T=1.0)
U0_ch = sp_init(Pi_ch, @(x,y) tanh(sin(x).*sin(y)/(sqrt(2)*eps_ch)));
u0_ch = sp_assemble(U0_ch, 2, Nch/2);
uhf = Of_ch.toSpec(u0_ch); uhi = Oi_ch.toSpec(U0_ch);
maxD_ch_std = 0; maxD_ch_iso = 0; tb_ch_std = NaN; D_ch_step1 = 0;
for n = 1:nst_ch
    uhf = If_ch.step(uhf); uhi = Ii_ch.step(uhi);
    if n == 1 || mod(n, 25) == 0
        uf = Of_ch.toPhys(uhf); Uc = sp_channels(uf, Pf_ch);
        Df = max(abs(Uc(:,:,[1 2 3])),[],'all');
        if n == 1, D_ch_step1 = Df; end
        maxD_ch_std = max(maxD_ch_std, Df);
        if isnan(tb_ch_std) && Df > 1e-2, tb_ch_std = n*tau_ch; end
        Vi = Oi_ch.toPhys(uhi);
        Di = max(abs(Vi(:,:,[1 2 3])),[],'all');
        maxD_ch_iso = max(maxD_ch_iso, Di);
    end
end
fprintf('  2D CH odd-odd (N=%d, eps=%g, tau=%.1e, T=%g): std D(step1)=%.2e, maxD=%.2e (T*=%g) | iso maxD=%g\n', ...
    Nch, eps_ch, tau_ch, T_ch, D_ch_step1, maxD_ch_std, tb_ch_std, maxD_ch_iso);

save(fullfile(res,'exp_extensibility.mat'), 'eq_sav','tb_sav_std','maxD_sav_std','tb_sav_iso','maxD_sav_iso', ...
    'tb_sav_sh_std','maxD_sav_sh_std','tb_sav_sh_t2','maxD_sav_sh_t2', ...
    'tb_cac_std','maxD_cac_std','maxMean_cac','maxD_cac_iso', ...
    'eq_ch','D_ch_step1','maxD_ch_std','tb_ch_std','maxD_ch_iso','eps_ch','tau_ch','T_ch');
fprintf('Saved exp_extensibility.mat\n');

% =========================================================================
function [u, r] = run_sav_1d(P, O, u0, tau, nst, C0)
% Linear 1st-order SAV (Shen-Xu-Yang 2018) using ONLY O.toSpec, O.toPhys, O.cube
% and scalar inner products: works identically for Fourier, reflection iso, and T2 iso!
[E1, ~] = sav_energy(P, u0);
r = sqrt(E1 + C0);
uh = O.toSpec(u0);
Ainv = 1 ./ (1 + tau * P.eps2 * P.K2);   % (I - tau*eps2*Delta)^{-1}
for n = 1:nst
    u = O.toPhys(uh);
    [E1, ~] = sav_energy(P, u);
    F = P.gam * (sav_cube(P, O, u) - u); % dE1/du = gam*(u^3 - u)
    b = F / sqrt(E1 + C0);
    bh = O.toSpec(b);
    % Solve (I - tau*eps2*Delta) u^{n+1} + tau * b * r^{n+1} = u^n
    %       r^{n+1} - r^n = 0.5 * <b, u^{n+1} - u^n>
    % Using Sherman-Morrison in spectral space:
    g1h = Ainv .* uh;
    g2h = Ainv .* bh;
    ip_bu = sav_ip(P, b, u);
    ip_bg1 = sav_ip(P, b, O.toPhys(g1h));
    ip_bg2 = sav_ip(P, b, O.toPhys(g2h));
    r_new = (r - 0.5*ip_bu + 0.5*ip_bg1) / (1 + 0.5*tau*ip_bg2);
    uh = g1h - tau * r_new * g2h;
    r = r_new;
end
u = O.toPhys(uh);
end

function [u, r, tb, maxD] = run_sav_1d_diag(P, O, u0, tau, nst, C0, th)
if nargin < 7, th = 0; end
[E1, ~] = sav_energy(P, u0);
r = sqrt(E1 + C0);
uh = O.toSpec(u0);
Ainv = 1 ./ (1 + tau * P.eps2 * P.K2);
tb = NaN; maxD = 0;
N = P.N;
for n = 1:nst
    u = O.toPhys(uh);
    [E1, ~] = sav_energy(P, u);
    F = P.gam * (sav_cube(P, O, u) - u);
    b = F / sqrt(E1 + C0);
    bh = O.toSpec(b);
    g1h = Ainv .* uh;
    g2h = Ainv .* bh;
    ip_bu  = sav_ip(P, b, u);
    ip_bg1 = sav_ip(P, b, O.toPhys(g1h));
    ip_bg2 = sav_ip(P, b, O.toPhys(g2h));
    r_new = (r - 0.5*ip_bu + 0.5*ip_bg1) / (1 + 0.5*tau*ip_bg2);
    uh = g1h - tau * r_new * g2h;
    r = r_new;
    if mod(n, 20) == 0
        un = O.toPhys(uh);
        if isfield(P,'rep') && strcmp(P.rep,'fourier')
            if th == 0
                Uc = sp_channels(un, P); D = max(abs(Uc(:,1)));
            else
                D = 0.5*max(abs(un + circshift(un, N/2)));
            end
        else
            D = max(abs(un(:,1)));
        end
        maxD = max(maxD, D);
        if isnan(tb) && D > 1e-2, tb = n*tau; end
    end
end
u = O.toPhys(uh);
end

function u3 = sav_cube(P, O, u)
if isfield(P,'rep') && strcmp(P.rep,'fourier'), u3 = u.^3; else, u3 = O.cube(u); end
end

function ip = sav_ip(P, a, b)
% Exact discrete L2 inner product (2pi/N) sum_{j=0}^{N-1} a(x_j) b(x_j)
% By orthogonality of isotypic channels, equals weighted sum over fundamental grid!
h = 2*pi / P.N;
if isfield(P,'rep') && strcmp(P.rep,'fourier')
    ip = h * sum(a(:) .* b(:));
elseif isfield(P,'rep') && strcmp(P.rep,'iso')
    % Fundamental grid 0..M has trapezoidal weights [1, 2, ..., 2, 1]
    M = P.M; w = 2*ones(M+1,1); w(1) = 1; w(end) = 1;
    ip = h * sum(w .* sum(a .* b, 2));
else % T2 isotypic: 2 copies of fundamental grid 0..M-1
    ip = 2*h * sum(a(:) .* b(:));
end
end

function [E1, u2] = sav_energy(P, u)
% E1(u) = (gam/4) * int (u^2 - 1)^2 dx = (gam/4)*(<u^2, u^2> - 2<u,u> + 2pi)
if isfield(P,'rep') && strcmp(P.rep,'fourier')
    u2 = u.^2;
    E1 = (P.gam/4) * (2*pi/P.N) * sum((u2 - 1).^2);
else
    % Graded square in Z_2: (u^2)_0 = u_0^2 + u_1^2, (u^2)_1 = 2 u_0 u_1
    u2 = [u(:,1).^2 + u(:,2).^2, 2*u(:,1).*u(:,2)];
    E1 = (P.gam/4) * (sav_ip(P, u2, u2) - 2*sav_ip(P, u, u) + 2*pi);
end
end
