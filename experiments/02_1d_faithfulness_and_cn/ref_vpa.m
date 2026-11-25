% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Independent high-precision reference (variable-precision arithmetic, 60
% digits) for the faithfulness test, kappa=1, u0 = sin x + delta, N=32,
% ETDRK2 with S=1, tau=0.05.  Compares D(t)=max|even part| with the
% isotypic double-precision solver.
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
digits(120);
N = 32; tau = 0.05; T = 120; S = 1; dels = [1e-10 1e-30 1e-60];
nst = round(T/tau); nout = round(0.5/tau);
k = [0:N/2, -N/2+1:-1]';
j = (0:N-1)';
jv = vpa(j); Fv = exp(-1i*2*vpa(pi)*(jv*jv.')/N);   % DFT matrix in VPA
Fiv = conj(Fv)/N;
Lk = k.^2 + S; z = -tau*Lk;
E = exp(vpa(z)); Q1 = (exp(vpa(z))-1)./vpa(z)*tau;
Q2 = (exp(vpa(z))-1-vpa(z))./vpa(z).^2*tau;
x = vpa(2*pi)*j/N;
ref = struct;
for q = 1:numel(dels)
    u = sin(x) + vpa(dels(q));
    uh = Fv*u; Dt = []; tt = [];
    Nf = @(uh) Fv*real_(2*(Fiv*uh) - (Fiv*uh).^3);
    for n = 0:nst
        if mod(n,nout)==0
            uu = real_(Fiv*uh); e = (uu + uu([1 end:-1:2]))/2;
            Dt(end+1) = double(max(abs(e))); tt(end+1) = n*tau; %#ok<SAGROW>
        end
        if n < nst
            Nu = Nf(uh); a = E.*uh + Q1.*Nu;
            uh = a + Q2.*(Nf(a) - Nu);
        end
    end
    ref(q).delta = dels(q); ref(q).t = tt; ref(q).D = Dt;
    % isotypic double
    P = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',S,'rep','iso'));
    pert = @(P) [dels(q)*ones(P.M+1,1) zeros(P.M+1,1)];
    R = sp_run(P,'etdrk2',tau,T,@(x) sin(x),0,struct,struct('nout',nout,'pert',pert));
    ref(q).Diso = R.D';
    P = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',S,'rep','fourier'));
    pf = @(P) dels(q)*ones(P.N,1);
    R = sp_run(P,'etdrk2',tau,T,@(x) sin(x),0,struct,struct('nout',nout,'pert',pf));
    ref(q).Dfou = R.D';
    lin = ref(q).D < 1e-2;
    rel = max(abs(ref(q).Diso(lin) - ref(q).D(lin))./ref(q).D(lin));
    fprintf('delta=%g: max rel diff iso vs vpa = %.2e ; fourier vs vpa = %.2e\n', dels(q), rel, ...
        max(abs(ref(q).Dfou(lin) - ref(q).D(lin))./ref(q).D(lin)));
end
save(fullfile(res,'ref_vpa.mat'),'ref');

function y = real_(y)
y = real(y);
end

