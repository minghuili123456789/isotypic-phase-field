% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Record full time-evolution trajectories (0 <= t <= 1000) and snapshots
% (t = 35, 44, 100, 1000) for all 7 high-order integrators in 1D
% (kappa = gamma = 1, S = 0, u0 = sin(x), N = 128) comparing standard FFT
% vs drop-in isotypic representation.
here = fileparts(mfilename('fullpath')); addpath(here);
res  = fullfile(here,'..','..','results');

N = 128; T_std = 100; T_long = 1000;
Pf = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier')); Of = sp_ops(Pf);
Pi = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));     Oi = sp_ops(Pi);

x_raw = (0:N-1)'*2*pi/N;
x = x_raw; x(x>=pi) = x(x>=pi) - 2*pi; [xs, ix] = sort(x); % sorted [-pi, pi)

U0 = sp_init(Pi, @(xx) sin(xx));
u0 = sp_assemble(U0, 1, N/2);

methods = {'BDF2-SAV','ARS443','SBDF3','SBDF4','IFRK4','SDC4','Strang2'};
taus    = [0.01,      0.02,    0.01,   0.01,   0.02,   0.05,  0.02];
tsnaps  = [35, 44, 100, 1000];
trec    = (0:0.5:100)';         % fine defect record on [0, 100]
tlong   = (0:5:1000)';          % long-time record on [0, 1000]

H = cell(numel(methods), 1);
for m = 1:numel(methods)
    nm = methods{m}; tau = taus(m);
    fprintf('Running %d/%d: %s (tau=%g) to T=%g...\n', m, numel(methods), nm, tau, T_long);
    [D_f, M_f, S_f, D_fl, M_fl] = run_track(nm, Pf, Of, u0, tau, T_long, trec, tlong, tsnaps, ix);
    [D_i, M_i, S_i, D_il, M_il] = run_track(nm, Pi, Oi, U0, tau, T_long, trec, tlong, tsnaps, ix);
    H{m} = struct('name',nm,'tau',tau,'D_f',D_f,'M_f',M_f,'S_f',{S_f}, ...
                  'D_fl',D_fl,'M_fl',M_fl,'D_i',D_i,'M_i',M_i,'S_i',{S_i}, ...
                  'D_il',D_il,'M_il',M_il);
    fprintf('  t=44: std mean=%+.3f D=%.3f | t=1000: std max|u|=%.3f | iso maxD(T=1000)=%g max|u(1000)|=%.4f\n', ...
        mean(S_f{2}), D_f(trec==44), max(abs(S_f{4})), max(D_il), max(abs(S_i{4})));
end

save(fullfile(res,'exp_highorder_evol.mat'), 'xs','tsnaps','trec','tlong','methods','taus','H');
fprintf('Saved exp_highorder_evol.mat\n');

% =========================================================================
function [Drec, Mrec, Snaps, Dlong, Mlong] = run_track(nm, P, O, u0, tau, Tmax, trec, tlong, tsnaps, ix)
isF = isfield(P,'rep') && strcmp(P.rep,'fourier');
nst = round(Tmax/tau);
Drec  = zeros(numel(trec),1);  Mrec  = Drec;
Dlong = zeros(numel(tlong),1); Mlong = Dlong;
Snaps = cell(numel(tsnaps),1);
r_step  = round(0.5/tau);
l_step  = round(5/tau);
s_steps = round(tsnaps/tau);

    function rec(un, nstep)
        if mod(nstep, r_step) == 0 && nstep*tau <= trec(end) + 1e-9
            idx = round(nstep/r_step) + 1;
            [D, mu] = calc_dm(un, P, isF);
            Drec(idx) = D; Mrec(idx) = mu;
        end
        if mod(nstep, l_step) == 0
            idx = round(nstep/l_step) + 1;
            [D, mu] = calc_dm(un, P, isF);
            Dlong(idx) = D; Mlong(idx) = mu;
        end
        k = find(s_steps == nstep, 1);
        if ~isempty(k)
            if isF, uf = un; else, uf = sp_assemble(un, 1, P.M); end
            Snaps{k} = uf(ix);
        end
    end

rec(u0, 0);
Lk = P.eps2 * P.K2;

switch nm
    case 'Strang2'
        E = exp(-0.5*tau*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            uh = E .* uh;
            uh = O.toSpec(exact_nl_flow(P, O.toPhys(uh), tau));
            uh = E .* uh;
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end

    case 'IFRK4'
        E2 = exp(-0.5*tau*Lk); E = exp(-tau*Lk);
        uh = O.toSpec(u0);
        for n = 1:nst
            k1 = O.N(uh);
            u2 = E2.*(uh + 0.5*tau*k1);      k2 = O.N(u2);
            u3 = E2.*uh + 0.5*tau*k2;        k3 = O.N(u3);
            u4 = E.*uh  + tau*(E2.*k3);      k4 = O.N(u4);
            uh = E.*uh + (tau/6)*(E.*k1 + 2*E2.*(k2 + k3) + k4);
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end

    case 'SBDF3'
        I4 = sp_integrator(P, O, 'etdrk4', tau);
        uh0 = O.toSpec(u0); N0 = O.N(uh0);
        uh1 = I4.step(uh0); N1 = O.N(uh1);
        uh2 = I4.step(uh1); N2 = O.N(uh2);
        Dinv = 1 ./ (11/6 + tau*Lk);
        for n = 3:nst
            uh = Dinv .* (3*uh2 - 1.5*uh1 + (1/3)*uh0 + tau*(3*N2 - 3*N1 + N0));
            uh0=uh1; N0=N1; uh1=uh2; N1=N2; uh2=uh; N2=O.N(uh);
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end

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
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end

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
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end

    case 'BDF2-SAV'
        C0 = 10;
        [E1, ~] = sav_energy(P, u0); r0 = sqrt(E1 + C0);
        uh0 = O.toSpec(u0);
        A1 = 1 ./ (1 + tau*Lk);
        b0 = P.gam*(sav_cube(P,O,u0) - u0) / sqrt(E1 + C0);
        g1h = A1.*uh0; g2h = A1.*O.toSpec(b0);
        r1 = (r0 - 0.5*sav_ip(P,b0,u0) + 0.5*sav_ip(P,b0,O.toPhys(g1h))) / (1 + 0.5*tau*sav_ip(P,b0,O.toPhys(g2h)));
        uh1 = g1h - tau*r1*g2h; u1 = O.toPhys(uh1);
        A2 = 1 ./ (1.5 + tau*Lk);
        for n = 2:nst
            ustar = 2*u1 - u0;
            [E1, ~] = sav_energy(P, ustar);
            b = P.gam*(sav_cube(P,O,ustar) - ustar) / sqrt(E1 + C0);
            bh = O.toSpec(b);
            g1h = A2 .* (2*uh1 - 0.5*uh0); g2h = A2 .* bh;
            g1 = O.toPhys(g1h); g2 = O.toPhys(g2h);
            rhs_r = (4*r1 - r0)/3 - (1/6)*sav_ip(P, b, 4*u1 - u0) + 0.5*sav_ip(P, b, g1);
            r2 = rhs_r / (1 + 0.5*tau*sav_ip(P, b, g2));
            uh2 = g1h - tau*r2*g2h; u2 = O.toPhys(uh2);
            u0 = u1; uh0 = uh1; r0 = r1;
            u1 = u2; uh1 = uh2; r1 = r2;
            if mod(n, r_step)==0 || any(s_steps==n), rec(u2, n); end
        end

    case 'SDC4'
        Sq = [5/24, 8/24, -1/24; -1/24, 8/24, 5/24];
        dt = tau/2; Dinv = 1 ./ (1 + dt*Lk);
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
            if mod(n, r_step)==0 || any(s_steps==n), rec(O.toPhys(uh), n); end
        end
end
end

function [D, mu] = calc_dm(un, P, isF)
if isF
    Uc = sp_channels(un, P); D = max(abs(Uc(:,1))); mu = abs(mean(un));
else
    D = max(abs(un(:,1)));
    % exact mean from even channel (0 when odd)
    M = P.M; w = 2*ones(M+1,1); w(1)=1; w(end)=1;
    mu = abs(sum(w.*un(:,1)) / P.N);
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
else
    M = P.M; w = 2*ones(M+1,1); w(1) = 1; w(end) = 1;
    ip = h * sum(w .* sum(a .* b, 2));
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
