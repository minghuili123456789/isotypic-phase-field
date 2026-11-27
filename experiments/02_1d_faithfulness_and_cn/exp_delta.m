% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Faithfulness test: odd data + genuine even perturbation of size delta.
% Exact solution breaks symmetry at T*(delta) ~ lambda^{-1} log(eta/delta).
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
cases = struct('name',{'k1','mk'}, ...
   'N',{128,1024},'eps2',{1,1},'gam',{1,400},'S',{1,800},'tau',{0.05,1e-3}, ...
   'T',{400,60},'u0',{@(x) sin(x), @(x) 0.5*(sin(x)+sin(4*x))}, ...
   'pfun',{@(x) ones(size(x)), @(x) cos(x)});
dl = [10.^-(2:2:30), 1e-50, 1e-100, 1e-200, 0];
variants = {'iso','fourier','single','filter'};
jobs = {};
for c = 1:numel(cases), for v = 1:numel(variants), for k = 1:numel(dl)
   jobs{end+1} = [c v k]; %#ok<SAGROW>
end, end, end
if isempty(gcp('nocreate')), cl=parcluster('Processes'); cl.NumWorkers=20; parpool(cl,6); end
out = cell(size(jobs));
parfor j = 1:numel(jobs)
   c = cases(jobs{j}(1)); v = variants{jobs{j}(2)}; del = dl(jobs{j}(3));
   rep = 'fourier'; prec = 'double'; post = [];
   if strcmp(v,'iso'), rep = 'iso'; end
   if strcmp(v,'single'), prec = 'single'; end
   if strcmp(v,'filter'), post = @(uh) 1i*imag(uh); end   % Li-Quan-Tang-Yang parity filter
   P = sp_setup(struct('d',1,'N',c.N,'eps2',c.eps2,'gam',c.gam,'S',c.S,'rep',rep,'prec',prec));
   pert = @(P) del*make_pert(P, c.pfun);
   R = sp_run(P,'etdrk4',c.tau,c.T,c.u0,0,struct, ...
       struct('nout',max(1,round(0.1/c.tau)),'pert',pert,'post',post,'stopD',0.5));
   out{j} = struct('case',c.name,'variant',v,'delta',del,'t',R.t,'D',R.D);
end
save(fullfile(res,'exp_delta.mat'),'out','dl','variants');
for j = 1:numel(out)
   o = out{j}; tb = o.t(find(o.D > 0.1,1)); if isempty(tb), tb = NaN; end
   fprintf('%s %-7s delta=%8.1e  T*=%8.2f\n', o.case, o.variant, o.delta, tb);
end

% Extended horizon T=500 for k1 isotypic at delta=1e-200 (needed by fig_faithful.m)
c1 = cases(1);
P_ext = sp_setup(struct('d',1,'N',c1.N,'eps2',c1.eps2,'gam',c1.gam,'S',c1.S,'rep','iso','prec','double'));
pert_ext = @(P) 1e-200*make_pert(P, c1.pfun);
R_ext = sp_run(P_ext,'etdrk4',c1.tau,500,c1.u0,0,struct, ...
    struct('nout',max(1,round(0.1/c1.tau)),'pert',pert_ext,'stopD',0.5));
tb200 = R_ext.t(find(R_ext.D > 0.1,1));
fprintf('k1 iso extended delta=1.0e-200 (T=500): T*=%8.2f\n', tb200);
save(fullfile(res,'exp_delta_extra.mat'),'tb200','R_ext');

function p = make_pert(P, pfun)
% even perturbation in physical representation P (exact even function)
M = P.M; xh = (0:M)'*P.h;
e = pfun(xh);
if strcmp(P.rep,'iso')
   p = zeros(M+1,2,P.prec); p(:,1) = e;
else
   p = cast(sp_assemble([e zeros(M+1,1)],1,M), P.prec);
end
end
