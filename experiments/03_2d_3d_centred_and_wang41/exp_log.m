% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% General nonlinearity: Flory--Huggins logarithmic potential (odd-exact form)
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
th = 0.8; thc = 1.6;
F = @(u) (th/2)*(log1p(u) - log1p(-u)) - thc*u;     % odd-exact evaluation
% (1) equivalence for non-symmetric data, 1D/2D
for d = 1:2
    N = [64 32]; N = N(d);
    base = struct('d',d,'N',N,'eps2',0.1,'gam',1,'S',2,'F',F);
    if d == 1, g0 = @(x) 0.5*sin(x)+0.3*cos(2*x)+0.1*sin(3*x+0.4);
    else, g0 = @(x,y) 0.5*sin(x).*cos(y)+0.2*cos(2*x+y)+0.1*sin(x-3*y); end
    Pf = base; Pf.rep='fourier'; Pf = sp_setup(Pf); Pi = base; Pi.rep='iso'; Pi = sp_setup(Pi);
    Of = sp_ops(Pf); Oi = sp_ops(Pi);
    If = sp_integrator(Pf,Of,'etdrk2',0.05); Ii = sp_integrator(Pi,Oi,'etdrk2',0.05);
    uf = Of.toSpec(sp_init(Pf,g0)); ui = Oi.toSpec(sp_init(Pi,g0));
    for n = 1:40, uf = If.step(uf); ui = Ii.step(ui); end
    a = sp_assemble(Oi.toPhys(ui), d, Pi.M); b = Of.toPhys(uf);
    fprintf('log-potential d=%d: |iso - fourier| = %.2e\n', d, max(abs(a(:)-b(:))));
end
% (2) symmetry: kappa=1, u0=0.9 sin x, S=0 (zero state unstable, rate 0.8)
out = struct;
for rep = {'fourier','iso'}
    P = sp_setup(struct('d',1,'N',128,'eps2',1,'gam',1,'S',0,'rep',rep{1},'F',F));
    R = sp_run(P,'etdrk2',0.01,400,@(x) 0.9*sin(x),0,struct,struct('nout',10));
    tb = R.t(find(R.D > 1e-2,1)); if isempty(tb), tb = NaN; end
    fprintf('log-potential 1D %-7s: T* = %g, max D = %g, final range [%g,%g]\n', rep{1}, tb, max(R.D), min(R.ufinal), max(R.ufinal));
    out.(rep{1}) = R;
end
% (3) 2D odd-odd log potential, eps=0.05 scaling, N=96 (non-equivariant FFT size)
ep = 0.05;
for rep = {'fourier','iso'}
    P = sp_setup(struct('d',2,'N',96,'eps2',1,'gam',1/ep^2,'S',8/ep^2,'rep',rep{1},'F',F));
    R = sp_run(P,'etdrk2',1e-3,5,@(x,y) 0.9*sin(x).*sin(y),[0 1 2],struct,struct('nout',50,'snap',[0 1 5]));
    tb = R.t(find(R.D > 1e-2,1)); if isempty(tb), tb = NaN; end
    fprintf('log-potential 2D %-7s: T* = %g, max D = %g\n', rep{1}, tb, max(R.D));
    out.(['d2_' rep{1}]) = R;
end
save(fullfile(res,'exp_log.mat'),'out');
