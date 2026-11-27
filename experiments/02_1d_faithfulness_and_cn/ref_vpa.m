% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Independent 120-digit variable-precision arithmetic (VPA) reference for the
% faithfulness test: kappa=1, u0 = sin x + delta, N=32, ETDRK2 with S=1, tau=1/20.
% All constants (pi, tau=1/20, delta=10^-p, DFT matrix, ETD multipliers) are
% constructed symbolically in 120-digit arithmetic (zero double-precision input bias).
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');
digits(120);
N = 32; tau = 0.05; tau_v = vpa(1)/20; T = 120; S = 1;
exps = [10 30 60]; dels = 10.^(-exps);
nst = round(T/tau); nout = round(0.5/tau);
k = [0:N/2, -N/2+1:-1]';
j = (0:N-1)';
jv = vpa(j); kv = vpa(k);
Fv = exp(-1i*2*vpa(pi)*(jv*jv.')/N);   % exact 120-digit DFT matrix
Fiv = conj(Fv)/N;
Lk_v = kv.^2 + vpa(S); zv = -tau_v*Lk_v;
E = exp(zv); Q1 = (E - 1)./zv * tau_v;
Q2 = (E - 1 - zv)./zv.^2 * tau_v;
x = 2*vpa(pi)*jv/N;
ref = struct;
for q = 1:numel(exps)
    del_v = vpa(10)^(-exps(q));
    u = sin(x) + del_v;
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
    fprintf('delta=1e-%d: max rel diff iso vs 120-digit VPA = %.2e ; fourier vs VPA = %.2e\n', exps(q), rel, ...
        max(abs(ref(q).Dfou(lin) - ref(q).D(lin))./ref(q).D(lin)));
end
save(fullfile(res,'ref_vpa.mat'),'ref');

function y = real_(y)
y = real(y);
end
