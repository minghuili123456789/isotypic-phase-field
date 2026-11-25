% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Which grid sizes make the standard FFT implementation non-equivariant?
addpath(fileparts(mfilename('fullpath')));
ep = 0.05; g = 1/ep^2; tanhs = @(v) tanh(v/(sqrt(2)*ep));
tests = {struct('d',1,'u0',@(x) tanhs(0.5*(sin(x)+sin(4*x))),'f',0), ...
         struct('d',2,'u0',@(x,y) tanhs(sin(x).*sin(y)),'f',[0 1 2]), ...
         struct('d',2,'u0',@(x,y) tanhs(sin(x + 2.5*sin(3*y))),'f',[0 3]), ...
         struct('d',3,'u0',@(x,y,z) tanhs(sin(x).*sin(y).*sin(z)),'f',0:6)};
names = {'1D odd','2D odd-odd','2D C2-odd','3D odd^3'};
Ns = [32 48 64 96 128 160 192 256];
for q = 1:numel(tests)
    fprintf('%-11s', names{q});
    for N = Ns
        if tests{q}.d == 3 && N > 128, fprintf('   -    '); continue; end
        P = sp_setup(struct('d',tests{q}.d,'N',N,'eps2',1,'gam',g,'S',2*g,'rep','fourier'));
        R = sp_run(P,'etdrk2',1e-3,0.005,tests{q}.u0,tests{q}.f,struct,struct('nout',1));
        fprintf(' %8.1e', max(R.D));
    end
    fprintf('\n');
end
