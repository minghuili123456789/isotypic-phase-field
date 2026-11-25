% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Exploration: which standard schemes break odd symmetry, and when?
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
cfg = {};
% Wang Ex.5.1: kappa=1, u0=sin x, N=128
for S = [0 1]
  for sch = {'imex1','etdrk2','etdrk4'}
    cfg{end+1} = struct('name',sprintf('k1_%s_S%d',sch{1},S),'N',128,'eps2',1,'gam',1,'S',S, ...
        'sch',sch{1},'tau',0.01,'T',3000,'u0',@(x) sin(x),'prm',struct); %#ok<SAGROW>
  end
end
cfg{end+1} = struct('name','k1_esp_Te1','N',128,'eps2',1,'gam',1,'S',1,'sch','esp','tau',1,'T',3000,'u0',@(x) sin(x),'prm',struct('Ns',24,'Niter',20));
cfg{end+1} = struct('name','k1_esp_Te20','N',128,'eps2',1,'gam',1,'S',1,'sch','esp','tau',20,'T',6000,'u0',@(x) sin(x),'prm',struct('Ns',24,'Niter',20));
% draft Sec 5.7: u_t = Lap u - f(u)/eps^2, eps=0.05, u0 = 0.5(sin x + sin 4x)
for sch = {'etdrk2','etdrk4','cn','imex1'}
  cfg{end+1} = struct('name',['mk_' sch{1}],'N',1024,'eps2',1,'gam',400,'S',800, ...
      'sch',sch{1},'tau',1e-3,'T',100,'u0',@(x) 0.5*(sin(x)+sin(4*x)),'prm',struct); %#ok<SAGROW>
end
% Wang Ex 5.2 (c): kappa^2=0.1, u0=(sin x + sin 4x)/2 ; (d) kappa^2=0.001, (sin x+sin 8x)/2
cfg{end+1} = struct('name','w52c_etdrk4','N',512,'eps2',0.1,'gam',1,'S',1,'sch','etdrk4','tau',0.05,'T',5000,'u0',@(x) 0.5*(sin(x)+sin(4*x)),'prm',struct);
cfg{end+1} = struct('name','w52d_etdrk4','N',512,'eps2',0.001,'gam',1,'S',1,'sch','etdrk4','tau',0.05,'T',5000,'u0',@(x) 0.5*(sin(x)+sin(8*x)),'prm',struct);
cfg{end+1} = struct('name','w52c_esp','N',512,'eps2',0.1,'gam',1,'S',1,'sch','esp','tau',2,'T',5000,'u0',@(x) 0.5*(sin(x)+sin(4*x)),'prm',struct('Ns',24,'Niter',20));
cfg{end+1} = struct('name','w52d_esp','N',512,'eps2',0.001,'gam',1,'S',1,'sch','esp','tau',10,'T',5000,'u0',@(x) 0.5*(sin(x)+sin(8*x)),'prm',struct('Ns',24,'Niter',20));
if isempty(gcp('nocreate')), cl=parcluster('Processes'); cl.NumWorkers=20; parpool(cl,12); end
out = cell(size(cfg));
parfor q = 1:numel(cfg)
  c = cfg{q};
  reps = {'fourier','iso'}; tmp = struct;
  for ir = 1:2
    P = sp_setup(struct('d',1,'N',c.N,'eps2',c.eps2,'gam',c.gam,'S',c.S,'rep',reps{ir}));
    nout = max(1,round(0.5/c.tau));
    R = sp_run(P, c.sch, c.tau, c.T, c.u0, 0, c.prm, struct('nout',nout));
    R = rmfield(R,'snap');
    tmp.(reps{ir}) = R;
  end
  tmp.cfg = rmfield(c,'u0'); out{q} = tmp;
end
save(fullfile(res,'explore1d.mat'),'out');
for q = 1:numel(out)
  F = out{q}.fourier; Iso = out{q}.iso;
  tb = F.t(find(F.D > 1e-2,1)); if isempty(tb), tb = NaN; end
  fprintf('%-14s fourier: Tbreak(D>1e-2)=%8.1f  maxD=%.2e | iso maxD=%g | cpu %.1f/%.1f s\n', ...
      out{q}.cfg.name, tb, max(F.D), max(Iso.D), F.cpu, Iso.cpu);
end

