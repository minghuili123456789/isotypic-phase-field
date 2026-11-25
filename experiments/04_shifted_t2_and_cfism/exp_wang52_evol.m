% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Long-time 1D shifted evolution comparison for Wang (2023) Ex. 5.2(a),
% theta = pi/4, kappa = 1, N = 512, T = 1000:
% Demonstrates visually that ALL existing methods fail while our
% non-interventional isotypic method (T2 ESP & CF-ISM) succeeds with D === 0.
here = fileparts(mfilename('fullpath')); addpath(here);
res  = fullfile(here,'..','..','results');

N = 512; th = pi/4; T = 1000;
x = -pi + (0:N-1)'*2*pi/N;
j_ref = mod(3*N/4 - (0:N-1)', N) + 1;          % exact index of -pi/2 - x_j
Dshift = @(u) 0.5*max(abs(u + u(j_ref)));      % shifted reflection defect
Danti  = @(u) 0.5*max(abs(u + circshift(u, N/2)));
u0 = sin(x + th);

tsnaps = [30, 38, 41, 42, 43, 45, 100, 1000];
tgrid  = (0:T)';                               % integer times 0..1000

% 1. Standard Fourier ETDRK4 (S=0, tau=0.02)
fprintf('Running 1/5: Standard Fourier ETDRK4 (S=0)...\n');
tau = 0.02; nst = round(T/tau);
P1 = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier'));
O1 = sp_ops(P1); I1 = sp_integrator(P1,O1,'etdrk4',tau);
uh = O1.toSpec(u0);
D1 = zeros(T+1,1); M1 = D1; S1 = cell(numel(tsnaps),1);
D1(1) = Dshift(u0); M1(1) = abs(mean(u0));
for n = 1:nst
    uh = I1.step(uh);
    if mod(n, round(1/tau)) == 0
        t = round(n*tau); u = O1.toPhys(uh);
        D1(t+1) = Dshift(u); M1(t+1) = abs(mean(u));
        k = find(tsnaps == t, 1); if ~isempty(k), S1{k} = u; end
    end
end

% 2. Wang (2023) ESP (S=1, Te=1, Ns=24, Niter=20) -- Wang Fig. 10(a) right!
fprintf('Running 2/5: Wang (2023) ESP (S=1, Te=1)...\n');
P2 = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',1,'rep','fourier'));
O2 = sp_ops(P2); I2 = sp_integrator(P2,O2,'esp',1,struct('Ns',24,'Niter',20));
uh = O2.toSpec(u0);
D2 = zeros(T+1,1); M2 = D2; S2 = cell(numel(tsnaps),1);
D2(1) = Dshift(u0); M2(1) = abs(mean(u0));
for t = 1:T
    uh = I2.step(uh); u = O2.toPhys(uh);
    D2(t+1) = Dshift(u); M2(t+1) = abs(mean(u));
    k = find(tsnaps == t, 1); if ~isempty(k), S2{k} = u; end
end

% 3. Origin-centered parity filter (LQTY 2022) on ETDRK4 (S=0, tau=0.02)
%    Filter applied on [-pi,pi) -> origin-odd Fourier modes are pure imaginary
%    on the [0,2pi) standard FFT ordering, so shift x = -pi..pi by exp(i*k*pi):
fprintf('Running 3/5: Origin-centered parity filter (LQTY 2022)...\n');
P3 = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier'));
O3 = sp_ops(P3); I3 = sp_integrator(P3,O3,'etdrk4',tau);
% On [0,2pi) grid u_std(j) = u0_shifted; origin-odd filter sets u -> 0.5*(u(x)-u(-x))
% Note: x_fft = (0:N-1)*2*pi/N; for our grid x = -pi + x_fft, -x corresponds to mod(-j, N)
% Wait: on x = -pi + j*2*pi/N, x=0 is at j=N/2, so -x_j = x_{N-j mod N} = x_j when j=0 or N/2!
j_orig = mod(N - (0:N-1)', N) + 1;             % index of -x_j on [-pi,pi) grid
uh = O3.toSpec(u0);
D3 = zeros(T+1,1); M3 = D3; S3 = cell(numel(tsnaps),1);
D3(1) = Dshift(u0); M3(1) = abs(mean(u0));
for n = 1:nst
    uh = I3.step(uh);
    u = O3.toPhys(uh);
    u = 0.5*(u - u(j_orig));                   % exact origin-centered odd projection
    uh = O3.toSpec(u);
    if mod(n, round(1/tau)) == 0
        t = round(n*tau);
        D3(t+1) = Dshift(u); M3(t+1) = abs(mean(u));
        k = find(tsnaps == t, 1); if ~isempty(k), S3{k} = u; end
    end
end

% 4. Proposed T2 isotypic ESP (EXACT same ESP S=1, Te=1, Ns=24, Niter=20; pure double!)
fprintf('Running 4/5: Proposed T2 isotypic ESP (S=1, Te=1, pure double)...\n');
P4 = sp_t2_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',1));
O4 = sp_t2_ops(P4); I4 = sp_integrator(P4,O4,'esp',1,struct('Ns',24,'Niter',20));
% sp_t2_init uses grid [0, 2*pi); since sin(x+pi/4) is 2pi-periodic, match x in [-pi,pi)
U0_t2 = sp_t2_init(P4, @(xx) sin(xx + th), 1);
vh = O4.toSpec(U0_t2);
D4_t2 = zeros(T+1,1); D4_sh = D4_t2; M4 = D4_t2; S4 = cell(numel(tsnaps),1);
u_init = circshift(sp_t2_assemble(U0_t2, 1, N/2), N/2);  % shift [0,2pi) -> [-pi,pi)
D4_t2(1) = max(abs(U0_t2(:,1))); D4_sh(1) = Dshift(u_init); M4(1) = abs(mean(u_init));
for t = 1:T
    vh = I4.step(vh); V = O4.toPhys(vh);
    D4_t2(t+1) = max(abs(V(:,1)));
    u = circshift(sp_t2_assemble(V, 1, N/2), N/2); % [-pi,pi) grid
    D4_sh(t+1) = Dshift(u); M4(t+1) = abs(mean(u));
    k = find(tsnaps == t, 1); if ~isempty(k), S4{k} = u; end
end

% 5. Proposed CF-ISM (d0=320 at t=0, pure double ETDRK4 S=0 for all t>0)
fprintf('Running 5/5: Proposed CF-ISM (d0=320 at t=0, pure double ETDRK4)...\n');
[c_cf, U0_cf] = cf_gauge_init(1, N, @(xx) sin(xx+th), @(xx) sin(xx+vpa(pi)/4), 320, 'refl');
P5 = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));
O5 = sp_ops(P5); I5 = sp_integrator(P5,O5,'etdrk4',tau);
vh = O5.toSpec(U0_cf);
D5 = zeros(T+1,1); S5 = cell(numel(tsnaps),1);
D5(1) = max(abs(U0_cf(:,1)));
kgrid = [0:N/2-1, -N/2:-1]';
shfac = exp(-1i*kgrid*(double(c_cf) - pi)); % map canonical [0,2pi) shifted by c_cf to [-pi,pi)
for n = 1:nst
    vh = I5.step(vh);
    if mod(n, round(1/tau)) == 0
        t = round(n*tau); V = O5.toPhys(vh);
        D5(t+1) = max(abs(V(:,1)));
        k = find(tsnaps == t, 1);
        if ~isempty(k)
            uc = sp_assemble(V, 1, N/2);
            S5{k} = real(ifft(fft(uc).*shfac));
        end
    end
end

fprintf('\n--- Snapshot diagnostics (theta = pi/4, Wang Ex. 5.2(a)) ---\n');
for k = 1:numel(tsnaps)
    t = tsnaps(k);
    fprintf('t=%4d | Std4 max|u|=%.3f D=%.2e | ESP max|u|=%.3f mean=%+.3f D=%.2e | LQTY err=%.3f | T2-ESP max|u|=%.4f D_t2=%g | CF D=%g diff_T2=%.2e\n', ...
        t, max(abs(S1{k})), D1(t+1), max(abs(S2{k})), mean(S2{k}), D2(t+1), ...
        max(abs(S3{k}-S4{k})), max(abs(S4{k})), D4_t2(t+1), D5(t+1), max(abs(S5{k}-S4{k})));
end

save(fullfile(res,'exp_wang52_evol.mat'), 'N','th','T','x','tsnaps','tgrid', ...
    'D1','M1','S1','D2','M2','S2','D3','M3','S3','D4_t2','D4_sh','M4','S4','D5','S5','c_cf');
fprintf('Saved exp_wang52_evol.mat\n');
