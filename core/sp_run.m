function R = sp_run(P, scheme, tau, T, fun, forbid, prm, opts)
%SP_RUN  Integrate to time T and record symmetry diagnostics.
%  forbid : vector of isotypic channel indices (0-based) that must vanish
%           for the symmetry under study (e.g. 0 for odd data in 1D).
%  opts.nout   : monitor every nout steps
%  opts.snap   : times at which to store full-grid snapshots
%  opts.stopD  : stop when defect exceeds this value (default inf)
if nargin < 7 || isempty(prm), prm = struct; end
if nargin < 8, opts = struct; end
if ~isfield(opts,'nout'),  opts.nout = 1; end
if ~isfield(opts,'snap'),  opts.snap = []; end
if ~isfield(opts,'stopD'), opts.stopD = inf; end
if ~isfield(opts,'pert'),  opts.pert = []; end
O  = sp_ops(P);
I  = sp_integrator(P, O, scheme, tau, prm);
[u0, xh] = sp_init(P, fun);
if ~isempty(opts.pert)            % add a prescribed (asymmetric) perturbation
    u0 = u0 + opts.pert(P);
end
uh = O.toSpec(u0);
nsteps = round(T/tau);
nm = floor(nsteps/opts.nout)+1;
R.t = zeros(nm,1); R.D = zeros(nm,1); R.mean = zeros(nm,1);
R.u0c = zeros(nm,1); R.E = zeros(nm,1);
snapT = opts.snap(:)'; R.snapT = snapT; R.snap = cell(1,numel(snapT));
km = 0; tic;
for n = 0:nsteps
    t = n*tau;
    if mod(n,opts.nout) == 0
        u = O.toPhys(uh);
        U = gather(sp_channels(u, P));
        km = km+1;
        R.t(km) = t;
        R.D(km) = defect(U, forbid, P.d);
        [R.mean(km), R.u0c(km)] = means(U, P);
        if R.D(km) > opts.stopD, break; end
    end
    for s = find(abs(snapT - t) < tau/2)
        u = O.toPhys(uh); U = gather(sp_channels(u, P));
        R.snap{s} = sp_assemble(double(U), P.d, P.M);
    end
    if n < nsteps
        uh = I.step(uh);
        if isfield(opts,'post') && ~isempty(opts.post), uh = opts.post(uh); end
        if any(~isfinite(uh(:)))
            warning('blow-up at t=%g', t); break;
        end
    end
end
R.t = R.t(1:km); R.D = R.D(1:km); R.mean = R.mean(1:km); R.u0c = R.u0c(1:km);
R.E = R.E(1:km);
R.cpu = toc; R.scheme = scheme; R.tau = tau; R.rep = P.rep; R.prec = P.prec;
u = O.toPhys(uh); U = gather(sp_channels(u,P));
R.ufinal = sp_assemble(double(U), P.d, P.M);
R.xh = xh;
end

function D = defect(U, forbid, d)
colons = repmat({':'},1,d);
D = 0;
for c = forbid(:)'
    A = U(colons{:},c+1);
    D = max(D, double(max(abs(A(:)))));
end
end

function [mu, u0c] = means(U, P)
% spatial mean (trapezoid on the fundamental domain of the totally even
% channel) and the value at the symmetry centre x=0
d = P.d; M = P.M;
colons = repmat({':'},1,d);
A = double(U(colons{:},1));
w = ones(M+1,1); w([1 end]) = 1/2;
W = w; for i = 2:d, W = W.*reshape(w,[ones(1,i-1) M+1]); end
mu = sum(A(:).*W(:))/M^d;
u0c = 0;
for c = 1:2^d
    B = U(colons{:},c); u0c = u0c + double(B(1));
end
end
