% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% 2D / 3D symmetry-breaking experiments: standard Fourier vs isotypic vs
% draft's zero-frequency suppression ('zmode', Fourier + etdrk2)
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
ep = 0.05; g = 1/ep^2;
tanhs = @(v) tanh(v/(sqrt(2)*ep));
cfg = {};
% (A) X-junction: odd in x and in y  -> only channel 3 (oo) allowed
%cfg{end+1} = struct('name','X2d','d',2,'N',256,'gam',g,'S',2*g,'tau',1e-3,'T',40, ...
%    'u0',@(x,y) tanhs(sin(x).*sin(y)),'forbid',[0 1 2],'snap',[0 5 10 20 40]);
% (B) draft long-time case: point-reflection odd (C2), channels ee,oo forbidden
cfg{end+1} = struct('name','C2d','d',2,'N',192,'gam',g,'S',2*g,'tau',2e-3,'T',200, ...
    'u0',@(x,y) tanhs(sin(x + 2.5*sin(3*y))),'forbid',[0 3],'snap',[0 1 5 20 100 200]);
% (C) 3D octant saddle: odd in x,y,z -> only channel 7 allowed
ep3 = 0.1; g3 = 1/ep3^2;
cfg{end+1} = struct('name','O3d','d',3,'N',96,'gam',g3,'S',2*g3,'tau',2e-3,'T',10, ...
    'u0',@(x,y,z) tanh(sin(x).*sin(y).*sin(z)/(sqrt(2)*ep3)),'forbid',0:6,'snap',[0 1 3 10]);
jobs = {};
for q = 1:numel(cfg)
    for v = {'fourier','iso','zmode'}
        jobs{end+1} = struct('q',q,'v',v{1}); %#ok<SAGROW>
    end
end
if isempty(gcp('nocreate')), cl=parcluster('Processes'); cl.NumWorkers=20; parpool(cl,6); end
out = cell(size(jobs));
parfor j = 1:numel(jobs)
    c = cfg{jobs{j}.q}; v = jobs{j}.v;
    rep = v; sch = 'etdrk2';
    if strcmp(v,'zmode'), rep = 'fourier'; sch = 'zmode'; end
    P = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',c.gam,'S',c.S,'rep',rep));
    R = sp_run(P, sch, c.tau, c.T, c.u0, c.forbid, struct, ...
        struct('nout',round(0.1/c.tau),'snap',c.snap));
    R.name = c.name; R.variant = v;
    out{j} = R;
end
save(fullfile(res,'exp2d3d_b.mat'),'out','-v7.3');
for j = 1:numel(out)
    R = out{j}; tb = R.t(find(R.D > 1e-2,1)); if isempty(tb), tb = NaN; end
    fprintf('%-5s %-8s Tbreak=%7.2f maxD=%.2e cpu=%.0fs\n', R.name, R.variant, tb, max(R.D), R.cpu);
end

