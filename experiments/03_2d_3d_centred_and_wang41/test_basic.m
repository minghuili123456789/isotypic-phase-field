% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% quick unit tests: equivalence of representations + exact invariance
addpath(fileparts(mfilename('fullpath')));
fprintf('sin(-x)==-sin(x): %d, tanh: %d\n', all(sin(-(0:1e-3:7))==-sin(0:1e-3:7)), all(tanh(-(0:1e-3:7))==-tanh(0:1e-3:7)));
for d = 1:3
    N = [64 32 16]; N = N(d);
    base = struct('d',d,'N',N,'eps2',0.1,'gam',1,'S',1);
    f = @(varargin) gen(varargin{:});
    Pf = base; Pf.rep='fourier'; Pf = sp_setup(Pf);
    Pi = base; Pi.rep='iso';     Pi = sp_setup(Pi);
    Of = sp_ops(Pf); Oi = sp_ops(Pi);
    for sch = {'imex1','etdrk2','etdrk4','cn'}
        If = sp_integrator(Pf,Of,sch{1},0.05); Ii = sp_integrator(Pi,Oi,sch{1},0.05);
        uf = Of.toSpec(sp_init(Pf,f)); ui = Oi.toSpec(sp_init(Pi,f));
        for n=1:20, uf = If.step(uf); ui = Ii.step(ui); end
        a = sp_assemble(Oi.toPhys(ui), d, Pi.M); b = Of.toPhys(uf);
        fprintf('d=%d %-7s  |iso-fourier| = %.2e\n', d, sch{1}, max(abs(a(:)-b(:))));
    end
end
% exact invariance test, 1D odd data, kappa=1
P = struct('d',1,'N',128,'eps2',1,'gam',1,'S',1,'rep','iso'); P = sp_setup(P);
R = sp_run(P,'etdrk4',0.05,60,@(x) sin(x),0,[],struct('nout',20));
fprintf('iso etdrk4: max defect over [0,60] = %g\n', max(R.D));
P.rep='fourier'; P = sp_setup(P);
R = sp_run(P,'etdrk4',0.05,60,@(x) sin(x),0,[],struct('nout',20));
fprintf('fourier etdrk4: defect at t=1,20,40,60: %s\n', mat2str(R.D([2 21 41 end])',3));

function v = gen(varargin)
switch nargin
  case 1, x=varargin{1}; v = 0.5*sin(x)+0.3*cos(2*x)+0.2*sin(3*x+0.4);
  case 2, [x,y]=varargin{:}; v = sin(x).*cos(y)+0.3*cos(2*x+y)+0.2*sin(x-3*y);
  case 3, [x,y,z]=varargin{:}; v = sin(x).*cos(y).*sin(z)+0.3*cos(2*x+y-z)+0.2*sin(x-3*y);
end
end
