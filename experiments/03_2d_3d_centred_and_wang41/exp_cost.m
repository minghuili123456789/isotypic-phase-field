% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% CPU cost per step (ETDRK4) and bitwise identity of the reduced isotypic run
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
maxNumCompThreads(4);
cases = {struct('d',1,'N',4096,'fun',@(x) 0.5*(sin(x)+sin(4*x)),'odd',1), ...
         struct('d',2,'N',256,'fun',@(x,y) tanh(sin(x).*sin(y)/0.07),'odd',3), ...
         struct('d',3,'N',64,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7), ...
         struct('d',3,'N',128,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7)};
T = table();
for q = 1:numel(cases)
    c = cases{q}; nst = 40; if c.d == 3 && c.N == 128, nst = 10; end
    tm = zeros(1,3); U = cell(1,3);
    for v = 1:3
        P = struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200);
        if v == 1, P.rep = 'fourier'; else, P.rep = 'iso'; end
        P = sp_setup(P);
        if v >= 2 && c.d >= 2, P.xform = 'mat'; end
        if v == 3, P.active = false(1,P.nc); P.active(c.odd+1) = true; end
        O = sp_ops(P); I = sp_integrator(P,O,'etdrk4',1e-3);
        uh = O.toSpec(sp_init(P,c.fun));
        uh = I.step(uh);                       % warm-up
        tic; for n = 1:nst, uh = I.step(uh); end; tm(v) = toc/nst;
        U{v} = O.toPhys(uh);
    end
    same = isequal(abs(U{2}), abs(U{3}));      % bitwise (up to signed zeros)
    T = [T; table(c.d, c.N, tm(1), tm(2), tm(3), same, ...
        'VariableNames',{'d','N','fourier','iso_full','iso_reduced','bitwise_identical'})]; %#ok<AGROW>
    disp(T(end,:));
end
save(fullfile(res,'exp_cost.mat'),'T');

