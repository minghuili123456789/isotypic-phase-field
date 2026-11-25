% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Exact reproduction of user's draft Section 5.7 (eps=1, N=4096, CN, T=100)
% plus shifted mixed-harmonic version via CF-ISM (pure double time stepping).
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');

N = 4096; tau = 0.01; T = 100;
u0 = @(x) 0.5*(sin(x) + sin(4*x));
snaps = [0 30 45 100];

% (1) Unshifted: standard CN, zero-mode CN, isotypic CN
P_f = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier'));
P_i = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));

fprintf('Running Draft 5.7 (eps=1, N=4096, CN, T=100)...\n');
R_std = sp_run(P_f,'cn',tau,T,u0,0,struct('tol',1e-10),struct('nout',20,'snap',snaps));
R_zm  = sp_run(P_f,'cn',tau,T,u0,0,struct('tol',1e-10),struct('nout',20,'snap',snaps,'post',@(uh) zm0(uh)));
R_iso = sp_run(P_i,'cn',tau,T,u0,0,struct('tol',1e-10),struct('nout',20,'snap',snaps));

fprintf('  Standard CN : maxD=%.3e, D(30)=%.2e, D(45)=%.2e, D(100)=%.2e, range=[%.4f,%.4f]\n', ...
    max(R_std.D), R_std.D(find(R_std.t>=30,1)), R_std.D(find(R_std.t>=45,1)), R_std.D(end), min(R_std.ufinal), max(R_std.ufinal));
fprintf('  Zero-mode CN: maxD=%.3e, D(30)=%.2e, D(45)=%.2e, D(100)=%.2e, range=[%.4f,%.4f]\n', ...
    max(R_zm.D), R_zm.D(find(R_zm.t>=30,1)), R_zm.D(find(R_zm.t>=45,1)), R_zm.D(end), min(R_zm.ufinal), max(R_zm.ufinal));
fprintf('  Isotypic CN : maxD=%g, range=[%.4f,%.4f]\n', max(R_iso.D), min(R_iso.ufinal), max(R_iso.ufinal));

% (2) Shifted Draft 5.7 (N=512, CN, T=100) for sh = pi/8, pi/4, 0.37123456789
Ns = 512; Ms = Ns/2;
shs = [pi/8, pi/4, 0.37123456789];
R_sh = cell(1, numel(shs));
for q = 1:numel(shs)
    sh = shs(q);
    f_d = @(x) 0.5*(sin(x + sh) + sin(4*(x + sh)));
    f_v = @(x) 0.5*(sin(x + vpa(sh)) + sin(4*(x + vpa(sh))));
    tic; [c_q, U0] = cf_gauge_init(1, Ns, f_d, f_v, 320, 'refl'); t_init = toc;
    % Evolve CF-ISM in pure double CN
    Pi = sp_setup(struct('d',1,'N',Ns,'eps2',1,'gam',1,'S',0,'rep','iso'));
    Oi = sp_ops(Pi); Ii = sp_integrator(Pi,Oi,'cn',tau,struct('tol',1e-10));
    vh = Oi.toSpec(U0); maxD_cf = max(abs(U0(:,1)));
    for n = 1:round(T/tau)
        vh = Ii.step(vh);
        if mod(n,100)==0, V = Oi.toPhys(vh); maxD_cf = max(maxD_cf, max(abs(V(:,1)))); end
    end
    V_cf = Oi.toPhys(vh);
    % Compare with fixed-origin iso & standard Fourier (both break!)
    Pf = sp_setup(struct('d',1,'N',Ns,'eps2',1,'gam',1,'S',0,'rep','fourier'));
    Of = sp_ops(Pf); If = sp_integrator(Pf,Of,'cn',tau,struct('tol',1e-10));
    xh = (0:Ns-1)'*2*pi/Ns; uh = Of.toSpec(f_d(xh));
    tb_f = NaN;
    for n = 1:round(T/tau)
        uh = If.step(uh);
        if isnan(tb_f) && abs(real(uh(1)))/Ns > 1e-2, tb_f = n*tau; end
    end
    u_f = Of.toPhys(uh);
    R_sh{q} = struct('sh',sh,'t_init',t_init,'init_ev',max(abs(U0(:,1))),'maxD_cf',maxD_cf, ...
        'cf_range',[min(V_cf(:,2)) max(V_cf(:,2))],'tb_std',tb_f,'std_range',[min(u_f) max(u_f)]);
    fprintf('  Shifted Draft5.7 sh=%.5f: t_init=%.2fs, init_ev=%g, CF-ISM maxD=%g | Std CN T*=%g, range=[%.3f,%.3f]\n', ...
        sh, t_init, max(abs(U0(:,1))), maxD_cf, tb_f, min(u_f), max(u_f));
end

save(fullfile(res,'exp_draft57.mat'),'R_std','R_zm','R_iso','R_sh','snaps');

function uh = zm0(uh)
uh(1) = 0;
end
