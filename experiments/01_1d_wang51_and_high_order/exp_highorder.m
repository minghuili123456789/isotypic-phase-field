% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Comprehensive comparison across 7 classic HIGH-ORDER phase-field time integrators:
%  1. Strang2  (2nd-order operator splitting with exact odd nonlinear flow via orbit eval)
%  2. BDF2-SAV (2nd-order Scalar Auxiliary Variable, Shen-Xu-Yang 2018)
%  3. ARS443   (3rd-order 4-stage L-stable IMEX-RK, Ascher-Ruuth-Spiteri 1997)
%  4. SBDF3    (3rd-order IMEX-BDF, Ascher-Ruuth-Wetton 1995)
%  5. SBDF4    (4th-order IMEX-BDF, Ascher-Ruuth-Wetton 1995)
%  6. IFRK4    (4th-order Integrating Factor Runge-Kutta)
%  7. SDC4     (4th-order Gauss-Lobatto Spectral Deferred Correction, Feng-Tang-Yang 2015)
% (Supplementing IMEX1, ETDRK2, ETDRK4, Crank-Nicolson, ESP already in Table 1!)
here = fileparts(mfilename('fullpath')); addpath(here);
res  = fullfile(here,'..','..','results');

N = 128; T = 300;
Pf = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier')); Of = sp_ops(Pf);
Pi = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));     Oi = sp_ops(Pi);

x = (0:N-1)'*2*pi/N;
u_ns = 0.4*sin(x) + 0.25*cos(2*x) + 0.15*sin(3*x + 0.7);
U_ns = sp_channels(u_ns, Pf);

U0 = sp_init(Pi, @(xx) sin(xx));
u0 = sp_assemble(U0, 1, N/2);

methods = {'Strang2','BDF2-SAV','ARS443','SBDF3','SBDF4','IFRK4','SDC4'};
taus    = [0.02,      0.01,      0.02,    0.01,   0.01,   0.02,   0.05];
R = struct('name',{},'tau',{},'order',{},'equiv',{},'std_Tstar',{},'std_maxD',{},'iso_maxD',{});

% Reference solution at t=0.5 (ETDRK4 tau=1e-4)
Iref = sp_integrator(Pf, Of, 'etdrk4', 1e-4);
uh_ref = Of.toSpec(u_ns);
for k = 1:5000, uh_ref = Iref.step(uh_ref); end
u_ref = Of.toPhys(uh_ref);

fprintf('%-10s | %6s | %5s | %9s | %9s | %9s | %8s\n', ...
    'Method','tau','order','equiv','std T*','std maxD','iso maxD');
fprintf('%s\n', repmat('-',1,72));

for m = 1:numel(methods)
    nm = methods{m}; tau = taus(m);
    % 1. Convergence order check at t = 0.5
    t1 = 0.05; t2 = 0.025;
    if strcmp(nm,'SBDF4') || strcmp(nm,'SBDF3') || strcmp(nm,'BDF2-SAV'), t1 = 0.02; t2 = 0.01; end
    u_c1 = run_ho(nm, Pf, Of, u_ns, t1, round(0.5/t1), false);
    u_c2 = run_ho(nm, Pf, Of, u_ns, t2, round(0.5/t2), false);
    ord = log2(max(abs(u_c1 - u_ref)) / max(abs(u_c2 - u_ref)));

    % 2. Non-symmetric equivalence check (50 steps)
    uf_eq = run_ho(nm, Pf, Of, u_ns, tau, 50, false);
    ui_eq = run_ho(nm, Pi, Oi, U_ns, tau, 50, false);
    eq = max(abs(uf_eq - sp_assemble(ui_eq, 1, N/2)));

    % 3. Symmetry test up to T = 300
    nst = round(T/tau);
    [~, tb_f, maxD_f] = run_ho(nm, Pf, Of, u0, tau, nst, true);
    [~, ~,    maxD_i] = run_ho(nm, Pi, Oi, U0, tau, nst, true);

    R(m) = struct('name',nm,'tau',tau,'order',ord,'equiv',eq, ...
                  'std_Tstar',tb_f,'std_maxD',maxD_f,'iso_maxD',maxD_i);
    fprintf('%-10s | %6.3f | %5.2f | %9.2e | %9.2f | %9.2e | %8g\n', ...
        nm, tau, ord, eq, tb_f, maxD_f, maxD_i);
end

% Shifted irrational benchmark (theta = 0.37123456789, N = 512, T = 200)
fprintf('\n=== Shifted irrational (theta=0.3712, N=512, T=200) ===\n');
Nw = 512; th = 0.37123456789; Tw = 200;
xf = -pi + (0:Nw-1)'*2*pi/Nw;
u0w = sin(xf + th);
Pfw = sp_setup(struct('d',1,'N',Nw,'eps2',1,'gam',1,'S',0,'rep','fourier')); Ofw = sp_ops(Pfw);
Pt2 = sp_t2_setup(struct('d',1,'N',Nw,'eps2',1,'gam',1,'S',0)); Ot2 = sp_t2_ops(Pt2);
U0w = sp_t2_init(Pt2, @(xx) sin(xx + th), 1);

sh_names = {'Strang2','BDF2-SAV','SBDF4','IFRK4','SDC4'};
sh_taus  = [0.02,     0.01,      0.01,   0.02,   0.05];
Rsh = struct('name',{},'std_Tstar',{},'std_maxD',{},'t2_maxD',{});
for m = 1:numel(sh_names)
    nm = sh_names{m}; tau = sh_taus(m); nst = round(Tw/tau);
    [~, tb_s, mD_s] = run_ho(nm, Pfw, Ofw, u0w, tau, nst, true, true);
    [~, ~,    mD_t] = run_ho(nm, Pt2, Ot2, U0w, tau, nst, true, true);
    Rsh(m) = struct('name',nm,'std_Tstar',tb_s,'std_maxD',mD_s,'t2_maxD',mD_t);
    fprintf('  %-10s shifted: std T*=%6.2f (maxD=%.2e) | T2-iso maxD=%g\n', nm, tb_s, mD_s, mD_t);
end

save(fullfile(res,'exp_highorder.mat'), 'R', 'Rsh');
fprintf('Saved exp_highorder.mat\n');

% =========================================================================
function [u, tb, maxD] = run_ho(nm, P, O, u0, tau, nst, track, shifted)
if nargin < 8, shifted = false; end
tb = NaN; maxD = 0;
isF = isfield(P,'rep') && strcmp(P.rep,'fourier');
stride = max(1, round(0.1/tau));
    function check_sym(un, nstep)
        if ~track || mod(nstep, stride) ~= 0, return; end
        if isF
            if shifted
                D = 0.5*max(abs(un + circshift(un, P.N/2)));
            else
                Uc = sp_channels(un, P); D = max(abs(Uc(:,1)));
            end
        else
            D = max(abs(un(:,1)));
        end
        maxD = max(maxD, D);
        if isnan(tb) && D > 1e-2, tb = nstep*tau; end
    end

Lk = P.eps2 * P.K2;

switch nm
    case 'Strang2'
        E = exp(-0.5*tau*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            uh = E .* uh;
            uh = O.toSpec(exact_nl_flow(P, O.toPhys(uh), tau));
            uh = E .* uh;
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh);

    case 'IFRK4'
        % 4th-order Integrating Factor Runge-Kutta (strictly positive exp factors)
        E2 = exp(-0.5*tau*Lk);
        E  = exp(-tau*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            k1 = O.N(uh);
            u2 = E2 .* (uh + 0.5*tau*k1);      k2 = O.N(u2);
            u3 = E2 .* uh + 0.5*tau*k2;        k3 = O.N(u3);
            u4 = E  .* uh + tau*(E2 .* k3);    k4 = O.N(u4);
            uh = E.*uh + (tau/6)*(E.*k1 + 2*E2.*(k2 + k3) + k4);
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh);

    case 'SBDF3'
        I4 = sp_integrator(P, O, 'etdrk4', tau);
        uh0 = O.toSpec(u0); N0 = O.N(uh0);
        uh1 = I4.step(uh0); N1 = O.N(uh1);
        uh2 = I4.step(uh1); N2 = O.N(uh2);
        Dinv = 1 ./ (11/6 + tau*Lk);
        for n = 3:nst
            uh = Dinv .* (3*uh2 - 1.5*uh1 + (1/3)*uh0 + tau*(3*N2 - 3*N1 + N0));
            uh0 = uh1; N0 = N1; uh1 = uh2; N1 = N2; uh2 = uh; N2 = O.N(uh);
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh2);

    case 'SBDF4'
        I4 = sp_integrator(P, O, 'etdrk4', tau);
        uh0 = O.toSpec(u0); N0 = O.N(uh0);
        uh1 = I4.step(uh0); N1 = O.N(uh1);
        uh2 = I4.step(uh1); N2 = O.N(uh2);
        uh3 = I4.step(uh2); N3 = O.N(uh3);
        Dinv = 1 ./ (25/12 + tau*Lk);
        for n = 4:nst
            uh = Dinv .* (4*uh3 - 3*uh2 + (4/3)*uh1 - 0.25*uh0 + tau*(4*N3 - 6*N2 + 4*N1 - N0));
            uh0=uh1; N0=N1; uh1=uh2; N1=N2; uh2=uh3; N2=N3; uh3=uh; N3=O.N(uh);
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh3);

    case 'ARS443'
        a = zeros(5,5); ah = zeros(5,5);
        a(2,2) = 1/2;
        a(3,2) = 1/6;   a(3,3) = 1/2;
        a(4,2) = -1/2;  a(4,3) = 1/2;   a(4,4) = 1/2;
        a(5,2) = 3/2;   a(5,3) = -3/2;  a(5,4) = 1/2;   a(5,5) = 1/2;
        ah(2,1) = 1/2;
        ah(3,1) = 11/18; ah(3,2) = 1/18;
        ah(4,1) = 5/6;   ah(4,2) = -5/6; ah(4,3) = 1/2;
        ah(5,1) = 1/4;   ah(5,2) = 7/4;  ah(5,3) = 3/4;  ah(5,4) = -7/4;
        Dinv = 1 ./ (1 + 0.5*tau*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            L_st = cell(5,1); N_st = cell(5,1);
            N_st{1} = O.N(uh); L_st{1} = -Lk.*uh;
            for i = 2:5
                rhs = uh;
                for j = 1:i-1
                    if a(i,j)  ~= 0, rhs = rhs + (tau*a(i,j))*L_st{j}; end
                    if ah(i,j) ~= 0, rhs = rhs + (tau*ah(i,j))*N_st{j}; end
                end
                Ui = Dinv .* rhs;
                L_st{i} = -Lk .* Ui;
                if i < 5, N_st{i} = O.N(Ui); end
            end
            uh = Ui;
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh);

    case 'BDF2-SAV'
        C0 = 10;
        [E1, ~] = sav_energy(P, u0); r0 = sqrt(E1 + C0);
        uh0 = O.toSpec(u0);
        A1 = 1 ./ (1 + tau*Lk);
        b0 = P.gam*(sav_cube(P,O,u0) - u0) / sqrt(E1 + C0);
        g1h = A1.*uh0; g2h = A1.*O.toSpec(b0);
        r1 = (r0 - 0.5*sav_ip(P,b0,u0) + 0.5*sav_ip(P,b0,O.toPhys(g1h))) / (1 + 0.5*tau*sav_ip(P,b0,O.toPhys(g2h)));
        uh1 = g1h - tau*r1*g2h;
        u1 = O.toPhys(uh1);
        A2 = 1 ./ (1.5 + tau*Lk);
        for n = 2:nst
            ustar = 2*u1 - u0;
            [E1, ~] = sav_energy(P, ustar);
            b = P.gam*(sav_cube(P,O,ustar) - ustar) / sqrt(E1 + C0);
            bh = O.toSpec(b);
            rhs_u = 2*uh1 - 0.5*uh0;
            g1h = A2 .* rhs_u;
            g2h = A2 .* bh;
            g1 = O.toPhys(g1h); g2 = O.toPhys(g2h);
            rhs_r = (4*r1 - r0)/3 - (1/6)*sav_ip(P, b, 4*u1 - u0) + 0.5*sav_ip(P, b, g1);
            r2 = rhs_r / (1 + 0.5*tau*sav_ip(P, b, g2));
            uh2 = g1h - tau*r2*g2h;
            u2 = O.toPhys(uh2);
            u0 = u1; uh0 = uh1; r0 = r1;
            u1 = u2; uh1 = uh2; r1 = r2;
            check_sym(u2, n);
        end
        u = u1;

    case 'SDC4'
        % 4th-order IMEX-SDC on 3 Chebyshev-Gauss-Lobatto nodes {0, 1/2, 1}, 5 sweeps
        Sq = [5/24, 8/24, -1/24; -1/24, 8/24, 5/24];
        dt = tau/2;
        Dinv = 1 ./ (1 + dt*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            U = cell(3,1); L = cell(3,1); Nv = cell(3,1);
            U{1} = uh; L{1} = -Lk.*uh; Nv{1} = O.N(uh);
            for m = 1:2
                U{m+1} = Dinv .* (U{m} + dt*Nv{m});
                L{m+1} = -Lk.*U{m+1}; Nv{m+1} = O.N(U{m+1});
            end
            for sw = 1:5
                Un = cell(3,1); Ln = cell(3,1); Nn = cell(3,1);
                Un{1} = uh; Ln{1} = L{1}; Nn{1} = Nv{1};
                for m = 1:2
                    Iint = tau*(Sq(m,1)*(L{1}+Nv{1}) + Sq(m,2)*(L{2}+Nv{2}) + Sq(m,3)*(L{3}+Nv{3}));
                    rhs = Un{m} + dt*(Nn{m} - Nv{m}) - dt*L{m+1} + Iint;
                    Un{m+1} = Dinv .* rhs;
                    Ln{m+1} = -Lk.*Un{m+1}; Nn{m+1} = O.N(Un{m+1});
                end
                U = Un; L = Ln; Nv = Nn;
            end
            uh = U{3};
            check_sym(O.toPhys(uh), n);
        end
        u = O.toPhys(uh);
end
end

function v = exact_nl_flow(P, u, dt)
e2 = exp(-2*P.gam*dt);
f = @(w) w ./ sqrt(w.^2 + (1 - w.^2)*e2);
if isfield(P,'rep') && strcmp(P.rep,'fourier')
    v = f(u);
else
    p = u(:,1) + u(:,2); q = u(:,1) - u(:,2);
    fp = f(p); fq = f(q);
    v = [0.5*(fp + fq), 0.5*(fp - fq)];
end
end

function u3 = sav_cube(P, O, u)
if isfield(P,'rep') && strcmp(P.rep,'fourier'), u3 = u.^3; else, u3 = O.cube(u); end
end

function ip = sav_ip(P, a, b)
h = 2*pi / P.N;
if isfield(P,'rep') && strcmp(P.rep,'fourier')
    ip = h * sum(a(:) .* b(:));
elseif isfield(P,'rep') && strcmp(P.rep,'iso')
    M = P.M; w = 2*ones(M+1,1); w(1) = 1; w(end) = 1;
    ip = h * sum(w .* sum(a .* b, 2));
else
    ip = 2*h * sum(a(:) .* b(:));
end
end

function [E1, u2] = sav_energy(P, u)
if isfield(P,'rep') && strcmp(P.rep,'fourier')
    u2 = u.^2; E1 = (P.gam/4)*(2*pi/P.N)*sum((u2 - 1).^2);
else
    u2 = [u(:,1).^2 + u(:,2).^2, 2*u(:,1).*u(:,2)];
    E1 = (P.gam/4)*(sav_ip(P, u2, u2) - 2*sav_ip(P, u, u) + 2*pi);
end
end
