% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Paired CPU cost per ETDRK4 step (4 threads, median of 5 trials) on identical
% initial data for Standard Fourier, Reflection Isotypic (A_R: full & reduced),
% and Half-Period Translation Isotypic (T_2: full & reduced), verifying exact
% signed equality isequal(U_full, U_reduced).
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
maxNumCompThreads(4);
ntrials = 5;
cases = {struct('d',1,'N',4096,'fun',@(x) sin(x),'odd',1), ...
         struct('d',2,'N',256,'fun',@(x,y) tanh(sin(x).*sin(y)/0.07),'odd',3), ...
         struct('d',3,'N',64,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7), ...
         struct('d',3,'N',128,'fun',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/0.14),'odd',7)};
T = table();
tm_all = cell(1, numel(cases));
for q = 1:numel(cases)
    c = cases{q}; nst = 40; if c.d == 3 && c.N == 128, nst = 10; end
    trials = zeros(ntrials, 5);
    % Setup all 5 solvers once on identical initial datum
    Pf = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200,'rep','fourier'));
    Of = sp_ops(Pf); If = sp_integrator(Pf,Of,'etdrk4',1e-3);
    uh_f0 = Of.toSpec(sp_init(Pf,c.fun)); If.step(uh_f0); % warm-up

    Oar = cell(1,2); Iar = cell(1,2); uhar0 = cell(1,2);
    for v = 1:2
        P = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200,'rep','iso'));
        if c.d >= 2, P.xform = 'mat'; end
        if v == 2, P.active = false(1,P.nc); P.active(c.odd+1) = true; end
        Oar{v} = sp_ops(P); Iar{v} = sp_integrator(P,Oar{v},'etdrk4',1e-3);
        uhar0{v} = Oar{v}.toSpec(sp_init(P,c.fun)); Iar{v}.step(uhar0{v}); % warm-up
    end

    Ot2 = cell(1,2); It2 = cell(1,2); uht20 = cell(1,2);
    for v = 1:2
        P = sp_t2_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200));
        if v == 2, P.active = false(1,P.nc); P.active(c.odd+1) = true; end
        Ot2{v} = sp_t2_ops(P); It2{v} = sp_integrator(P,Ot2{v},'etdrk4',1e-3);
        uht20{v} = Ot2{v}.toSpec(sp_t2_init(P,c.fun,c.odd)); It2{v}.step(uht20{v}); % warm-up
    end

    Uar = cell(1,2); Ut2 = cell(1,2);
    for tr = 1:ntrials
        uh = uh_f0; tic; for n = 1:nst, uh = If.step(uh); end; trials(tr,1) = toc/nst;
        for v = 1:2
            uh = uhar0{v}; tic; for n = 1:nst, uh = Iar{v}.step(uh); end; trials(tr,1+v) = toc/nst;
            if tr == 1, Uar{v} = Oar{v}.toPhys(uh); end
        end
        for v = 1:2
            uh = uht20{v}; tic; for n = 1:nst, uh = It2{v}.step(uh); end; trials(tr,3+v) = toc/nst;
            if tr == 1, Ut2{v} = Ot2{v}.toPhys(uh); end
        end
    end
    same_ar = isequal(Uar{1}, Uar{2});
    same_t2 = isequal(Ut2{1}, Ut2{2});
    tm = median(trials, 1);
    tm_min = min(trials, [], 1);
    tm_max = max(trials, [], 1);
    tm_all{q} = trials;

    T = [T; table(c.d, c.N, tm(1), tm(2), tm(3), tm(4), tm(5), same_ar, same_t2, ...
        'VariableNames',{'d','N','fourier','ar_full','ar_reduced','t2_full','t2_reduced','same_ar','same_t2'})]; %#ok<AGROW>
    fprintf('d=%d N=%4d (med of %d): std=%.3e [%.2e,%.2e] | AR_f=%.3e (%.2fx) AR_r=%.3e (%.2fx) | T2_f=%.3e (%.2fx) T2_r=%.3e (%.2fx) | eq=(%d,%d)\n', ...
        c.d, c.N, ntrials, tm(1), tm_min(1), tm_max(1), ...
        tm(2), tm(2)/tm(1), tm(3), tm(3)/tm(1), tm(4), tm(4)/tm(1), tm(5), tm(5)/tm(1), same_ar, same_t2);
end
save(fullfile(res,'exp_cost.mat'),'T','tm_all','ntrials');
