% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Compare Standard FFT, Fixed-origin iso, T2 iso, and CF-ISM on:
% (1) 2D irrational-shift X-junction (N=128, eps=0.05, T=2):
%     Demonstrates BOTH: T2 prevents mean/even-harmonic collapse (D_T2 === 0),
%     while CF-ISM prevents EVEN reflection-breaking reconnection to diagonal stripes (D_refl === 0)!
% (2) 2D irrational-shift mixed-harmonic X-junction (N=128, eps=0.1, T=4)
%     with snapshots for figure!
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');

th1 = 0.37123456; th2 = 0.71828182;
ep = 0.1; g = 1/ep^2; S = 2*g; tau = 1e-3; Tm = 4; Nm = 128;
snaps = [0 0.5 1.5 4];
fm_d = @(x,y) tanh((sin(x+th1)+0.4*sin(4*(x+th1))).*(sin(y+th2)+0.4*sin(4*(y+th2)))/(sqrt(2)*ep));
fm_v = @(x,y) tanh((sin(x+vpa(th1))+vpa(0.4)*sin(4*(x+vpa(th1)))).*(sin(y+vpa(th2))+vpa(0.4)*sin(4*(y+vpa(th2))))/(sqrt(vpa(2))*vpa(ep)));

% CF-ISM (d0=320 at t=0, pure double ETDRK2 for t>0)
[c2_q, U0_cf] = cf_gauge_init(2, Nm, fm_d, fm_v, 320, 'refl');
Pm = sp_setup(struct('d',2,'N',Nm,'eps2',1,'gam',g,'S',S,'rep','iso')); Pm.xform = 'mat';
Om = sp_ops(Pm); Im = sp_integrator(Pm,Om,'etdrk2',tau);
vh = Om.toSpec(U0_cf);
nst = round(Tm/tau); nout = 50;
t = (0:nout:nst)'*tau; D_cf = zeros(size(t)); snap_cf = cell(size(snaps)); km = 0;
% shift back to original grid for plotting ONLY (never fed back!)
kf = [0:Nm/2, -Nm/2+1:-1]'; [KX,KY] = ndgrid(kf,kf);
c_d = double(c2_q); ph_back = exp(-1i*(KX*c_d(1) + KY*c_d(2)));
for n = 0:nst
    if mod(n,nout)==0
        km = km+1; V = Om.toPhys(vh);
        D_cf(km) = max(abs(reshape(V(:,:,1:3),[],1)));
    end
    s = find(abs(snaps - n*tau) < tau/2);
    if ~isempty(s)
        vf = sp_assemble(Om.toPhys(vh), 2, Pm.M);
        snap_cf{s} = real(ifftn(fftn(vf).*ph_back));
    end
    if n < nst, vh = Im.step(vh); end
end

% Standard Fourier on original grid (measure shifted reflection defect by shifting forward by +c_d)
Pmf = sp_setup(struct('d',2,'N',Nm,'eps2',1,'gam',g,'S',S,'rep','fourier'));
Omf = sp_ops(Pmf); Imf = sp_integrator(Pmf,Omf,'etdrk2',tau);
xm = (0:Nm-1)'*2*pi/Nm; [Xm,Ym] = ndgrid(xm,xm);
uh = Omf.toSpec(fm_d(Xm,Ym));
D_std = zeros(size(t)); mean_std = zeros(size(t)); snap_std = cell(size(snaps)); km = 0;
ph_fwd = conj(ph_back);
for n = 0:nst
    if mod(n,nout)==0
        km = km+1; u = Omf.toPhys(uh);
        mean_std(km) = abs(mean(u(:)));
        % Measure reflection defect about c_d via shifted spectral channels
        vs = real(ifftn(fftn(u).*ph_fwd));
        Uc = sp_channels(vs, Pmf);
        D_std(km) = max(abs(reshape(Uc(:,:,1:3),[],1)));
    end
    s = find(abs(snaps - n*tau) < tau/2);
    if ~isempty(s), snap_std{s} = Omf.toPhys(uh); end
    if n < nst, uh = Imf.step(uh); end
end
fprintf('2D shifted mixed-harm X-junction: Std maxD=%.3f (max|mean|=%.3e) | CF-ISM maxD=%g\n', ...
    max(D_std), max(mean_std), max(D_cf));
save(fullfile(res,'exp_cfism_2d.mat'),'t','D_cf','D_std','mean_std','snaps','snap_cf','snap_std','c_d');
