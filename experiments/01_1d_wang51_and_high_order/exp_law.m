% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Breakdown-time law: kappa=1, u0 = sin x, S=0 (no stabiliser), N=128.
% (a) isotypic runs with a genuine constant perturbation delta -> T*(delta)
% (b) standard Fourier runs (no perturbation) in double and single precision
%     for several integrators and step sizes -> T*(u_mach)
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
base = struct('d',1,'N',128,'eps2',1,'gam',1,'S',0);
u0 = @(x) sin(x);
jobs = {};
for del = 10.^-(2:2:40), jobs{end+1} = struct('type','iso','sch','etdrk4','tau',0.01,'del',del,'prec','double','prm',struct); end %#ok<SAGROW>
schs = {'imex1','etdrk2','etdrk4','cn'};
for s = 1:numel(schs), for tau = [0.001 0.01 0.1], for pr = {'double','single'}
    jobs{end+1} = struct('type','fourier','sch',schs{s},'tau',tau,'del',0,'prec',pr{1},'prm',struct); %#ok<SAGROW>
end, end, end
for Te = [1 5 20], for pr = {'double','single'}
    jobs{end+1} = struct('type','fourier','sch','esp','tau',Te,'del',0,'prec',pr{1},'prm',struct('Ns',24,'Niter',20)); %#ok<SAGROW>
end, end
if isempty(gcp('nocreate')), cl=parcluster('Processes'); cl.NumWorkers=20; parpool(cl,8); end
out = cell(size(jobs));
parfor j = 1:numel(jobs)
    J = jobs{j};
    P = base; P.rep = J.type; P.prec = J.prec; P = sp_setup(P);
    prm = J.prm; if strcmp(J.sch,'cn'), prm.tol = 1e-13; if strcmp(J.prec,'single'), prm.tol = 1e-6; end, end
    del = J.del;
    pert = @(P) mkp(P, del);
    R = sp_run(P, J.sch, J.tau, 300, u0, 0, prm, ...
        struct('nout',max(1,round(0.05/J.tau)),'pert',pert,'stopD',0.5));
    tb = R.t(find(R.D > 1e-2,1)); if isempty(tb), tb = NaN; end
    J.Tstar = tb; J.t = R.t; J.D = R.D;
    out{j} = J;
end
save(fullfile(res,'exp_law.mat'),'out');
for j = 1:numel(out)
    J = out{j};
    fprintf('%-7s %-6s tau=%-6g prec=%-6s delta=%-8.1e T*=%.2f\n', J.type, J.sch, J.tau, J.prec, J.del, J.Tstar);
end

function p = mkp(P, del)
M = P.M;
if strcmp(P.rep,'iso'), p = zeros(M+1,2,P.prec); p(:,1) = del;
else, p = zeros(P.N,1,P.prec) + cast(del,P.prec); end
end
