% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Paired CPU cost per ETDRK4 step (4 threads) on identical initial data
% for Standard Fourier, Reflection Isotypic (A_R: full & reduced), and
% Half-Period Translation Isotypic (T_2: full & reduced), verifying exact
% signed equality isequal(U_full, U_reduced) (zero sign-flip blind spots).
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
maxNumCompThreads(4);
cases = {struct('d',1,'N',4096,'fun',@(x) sin(x),'odd',1), ...
         struct('d',2,'N',256,'fun',@(x,y) tanh(sin(x).*sin(y)/0.07),'odd',3), ...
         struct('d',3,'N',64,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7), ...
         struct('d',3,'N',128,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7)};
T = table();
for q = 1:numel(cases)
    c = cases{q}; nst = 40; if c.d == 3 && c.N == 128, nst = 10; end
    tm = zeros(1,5);
    % 1. Standard Fourier
    Pf = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200,'rep','fourier'));
    Of = sp_ops(Pf); If = sp_integrator(Pf,Of,'etdrk4',1e-3);
    uh = Of.toSpec(sp_init(Pf,c.fun));
    uh = If.step(uh); % warm-up
    tic; for n = 1:nst, uh = If.step(uh); end; tm(1) = toc/nst;

    % 2 & 3. A_R full and reduced
    Uar = cell(1,2);
    for v = 1:2
        P = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200,'rep','iso'));
        if c.d >= 2, P.xform = 'mat'; end
        if v == 2, P.active = false(1,P.nc); P.active(c.odd+1) = true; end
        O = sp_ops(P); I = sp_integrator(P,O,'etdrk4',1e-3);
        uh = O.toSpec(sp_init(P,c.fun));
        uh = I.step(uh); % warm-up
        tic; for n = 1:nst, uh = I.step(uh); end; tm(1+v) = toc/nst;
        Uar{v} = O.toPhys(uh);
    end
    same_ar = isequal(Uar{1}, Uar{2}); % exact signed equality (up to +/-0)

    % 4 & 5. T_2 full and reduced (exact same initial datum & threads!)
    Ut2 = cell(1,2);
    for v = 1:2
        P = sp_t2_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200));
        if v == 2, P.active = false(1,P.nc); P.active(c.odd+1) = true; end
        O = sp_t2_ops(P); I = sp_integrator(P,O,'etdrk4',1e-3);
        U0 = sp_t2_init(P, c.fun, c.odd);
        uh = O.toSpec(U0);
        uh = I.step(uh); % warm-up
        tic; for n = 1:nst, uh = I.step(uh); end; tm(3+v) = toc/nst;
        Ut2{v} = O.toPhys(uh);
    end
    same_t2 = isequal(Ut2{1}, Ut2{2});

    T = [T; table(c.d, c.N, tm(1), tm(2), tm(3), tm(4), tm(5), same_ar, same_t2, ...
        'VariableNames',{'d','N','fourier','ar_full','ar_reduced','t2_full','t2_reduced','same_ar','same_t2'})]; %#ok<AGROW>
    disp(T(end,:));
end
save(fullfile(res,'exp_cost.mat'),'T');
