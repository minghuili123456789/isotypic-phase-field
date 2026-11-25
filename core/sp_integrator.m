function I = sp_integrator(P, O, scheme, tau, prm)
%SP_INTEGRATOR  One-step maps uh -> uh_new in spectral space for several
% time integrators.  All integrators are built ONLY from
%   (i) diagonal Fourier multipliers m(|k|^2) applied channel-wise,
%  (ii) the nonlinear map O.N, and
% (iii) linear combinations,
% hence they inherit the exact floating-point invariance of the isotypic
% channels (Theorem in the paper) without any modification.
%
% scheme: 'imex1' | 'etdrk2' | 'etdrk4' | 'esp' | 'cn' | 'zmode' (draft's
%         zero-frequency suppression wrapped around etdrk2, for comparison)
if nargin < 5, prm = struct; end
Lk = P.Lk;                 % stabilised symbol  eps2|k|^2 + S
z  = -tau*Lk;              % z = tau*(linear operator)
I.tau = tau; I.scheme = scheme;
switch scheme
    case 'imex1'
        D = 1./(1 + tau*Lk);
        I.step = @(uh) D.*(uh + tau*O.N(uh));
    case {'etdrk2','zmode'}
        E  = exp(z);
        Q1 = tau*phik(z,1);
        Q2 = tau*phik(z,2);
        st = @(uh) etd2(uh, O, E, Q1, Q2);
        if strcmp(scheme,'zmode')
            I.step = @(uh) zm(st(uh));
        else
            I.step = st;
        end
    case 'etdrk4'
        E  = exp(z); E2 = exp(z/2);
        Q  = (tau/2)*phik(z/2,1);
        [f1,f2,f3] = etd4coef(z, tau);
        I.step = @(uh) etd4(uh, O, E, E2, Q, f1, f2, f3);
    case 'esp'
        Ns = prm.Ns; Nit = prm.Niter;     % Te = tau
        [Ej, W] = esp_weights(Lk, tau, Ns);
        I.step = @(uh) esp_step(uh, O, Ej, W, Ns, Nit);
    case 'cn'
        S  = P.S; A = P.eps2*P.K2;
        lhs = 1./(1 + tau/2*(A + S));
        rhs0 = (1 - tau/2*A);
        tol = 1e-12; if isfield(prm,'tol'), tol = prm.tol; end
        I.step = @(uh) cn_step(uh, O, P, lhs, rhs0, tau, tol);
    otherwise
        error('unknown scheme');
end
end

% ---------------------------------------------------------------------
function un = etd2(uh, O, E, Q1, Q2)
Nu = O.N(uh);
a  = E.*uh + Q1.*Nu;
un = a + Q2.*(O.N(a) - Nu);
end

function uh = zm(uh)
uh(1) = 0;          % draft's a-posteriori zero-frequency suppression
end

function un = etd4(uh, O, E, E2, Q, f1, f2, f3)
Nu = O.N(uh);
a  = E2.*uh + Q.*Nu;   Na = O.N(a);
b  = E2.*uh + Q.*Na;   Nb = O.N(b);
c  = E2.*a  + Q.*(2*Nb - Nu); Nc = O.N(c);
un = E.*uh + f1.*Nu + 2*f2.*(Na + Nb) + f3.*Nc;
end

function un = cn_step(uh, O, P, lhs, rhs0, tau, tol)
g = P.gam; S = P.S;
gfun = @(vh) O.N(vh) - S*vh;                 % gam*(u-u^3) in spectral space
g0 = gfun(uh);
base = rhs0.*uh + tau/2*g0;
v = uh;
for it = 1:200
    vn = lhs.*(base + tau/2*(gfun(v) + S*v));
    dv = O.toPhys(vn - v);
    v = vn;
    if max(abs(dv(:))) < tol, break; end
end
un = v;
end

function un = esp_step(uh, O, Ej, W, Ns, Nit)
% Wang (JCP 2023) Algorithm 1: constant initialisation + Gauss--Seidel sweeps
Nh = cell(1,Ns+1);
N0 = O.N(uh);
for j = 1:Ns+1, Nh{j} = N0; end
Uj = uh;
for it = 1:Nit
    for j = 2:Ns+1
        acc = Ej{j}.*uh;
        for i = 1:Ns+1
            acc = acc + W{i,j}.*Nh{i};
        end
        Uj = acc;
        Nh{j} = O.N(Uj);
    end
end
un = Uj;
end

% ---------------------------------------------------------------------
function p = phik(z, k)
% phi_k(z) evaluated by the Kassam--Trefethen contour integral (real z<=0)
Mq = 64; r = exp(1i*pi*((1:Mq)-.5)/Mq);
sz = size(z); zz = z(:);
if isa(zz,'gpuArray'), zz = gather(zz); end
cls = class(zz); zz = double(zz);
p = zeros(size(zz));
big = abs(zz) > 1;
% direct formulas where safe
zb = zz(big);
switch k
    case 1, p(big) = (exp(zb) - 1)./zb;
    case 2, p(big) = (exp(zb) - 1 - zb)./zb.^2;
    case 3, p(big) = (exp(zb) - 1 - zb - zb.^2/2)./zb.^3;
end
zs = zz(~big);
if ~isempty(zs)
    Zc = zs + r;   % circle of radius 1
    switch k
        case 1, F = (exp(Zc) - 1)./Zc;
        case 2, F = (exp(Zc) - 1 - Zc)./Zc.^2;
        case 3, F = (exp(Zc) - 1 - Zc - Zc.^2/2)./Zc.^3;
    end
    p(~big) = real(mean(F,2));
end
p = reshape(cast(p,cls), sz);
if isa(z,'gpuArray'), p = gpuArray(p); end
end

function [f1,f2,f3] = etd4coef(z, tau)
% Cox--Matthews ETDRK4 coefficients via contour integral (Kassam--Trefethen)
Mq = 64; r = exp(1i*pi*((1:Mq)-.5)/Mq);
sz = size(z); zz = z(:); isg = isa(zz,'gpuArray');
if isg, zz = gather(zz); end
cls = class(zz); zz = double(zz);
f1 = zeros(size(zz)); f2 = f1; f3 = f1;
big = abs(zz) > 2;
zb = zz(big);
f1(big) = (-4 - zb + exp(zb).*(4 - 3*zb + zb.^2))./zb.^3;
f2(big) = ( 2 + zb + exp(zb).*(-2 + zb))./zb.^3;
f3(big) = (-4 - 3*zb - zb.^2 + exp(zb).*(4 - zb))./zb.^3;
zs = zz(~big);
if ~isempty(zs)
    Zc = zs + 2*r;
    f1(~big) = real(mean((-4 - Zc + exp(Zc).*(4 - 3*Zc + Zc.^2))./Zc.^3, 2));
    f2(~big) = real(mean(( 2 + Zc + exp(Zc).*(-2 + Zc))./Zc.^3, 2));
    f3(~big) = real(mean((-4 - 3*Zc - Zc.^2 + exp(Zc).*(4 - Zc))./Zc.^3, 2));
end
f1 = reshape(cast(tau*f1,cls),sz); f2 = reshape(cast(tau*f2,cls),sz);
f3 = reshape(cast(tau*f3,cls),sz);
if isg, f1 = gpuArray(f1); f2 = gpuArray(f2); f3 = gpuArray(f3); end
end

function [Ej, W] = esp_weights(Lk, Te, Ns)
% ESP weights  W_{ij}(lam) = int_0^{t_j} exp((s-t_j) lam) ell_i(s) ds  on
% Chebyshev--Gauss--Lobatto nodes, computed by graded Gauss--Legendre
% quadrature in sigma = t_j - s (robust for very large lam).
cls = class(Lk); isg = isa(Lk,'gpuArray');
if isg, Lk = gather(Lk); end
lamAll = double(Lk(:));
[lam,~,ic] = unique(lamAll);
tn = Te*(1 - cos(pi*(0:Ns)/Ns))/2;          % CGL nodes on [0,Te]
wb = ones(1,Ns+1); wb([1 end]) = 1/2; wb(2:2:end) = -wb(2:2:end); % bary w.
[xg, wg] = gauleg(32);
Ej = cell(1,Ns+1); W = cell(Ns+1,Ns+1);
for j = 1:Ns+1
    Ej{j} = reshape(cast(exp(-tn(j)*lamAll),cls), size(Lk));
    if isg, Ej{j} = gpuArray(Ej{j}); end
end
for j = 2:Ns+1
    tj = tn(j);
    br = [0, tj*2.^(-(60:-1:1)), tj];       % graded breakpoints in sigma
    sig = []; ws = [];
    for p = 1:numel(br)-1
        a = br(p); b = br(p+1);
        sig = [sig; (a+b)/2 + (b-a)/2*xg]; %#ok<AGROW>
        ws  = [ws;  (b-a)/2*wg];           %#ok<AGROW>
    end
    s = tj - sig;                            % quadrature nodes in s
    Lmat = baryL(s, tn, wb);                 % ell_i(s): nq x (Ns+1)
    Ex = exp(-lam*sig');                     % nlam x nq
    Wl = Ex*(ws.*Lmat);                      % nlam x (Ns+1)
    for i = 1:Ns+1
        W{i,j} = reshape(cast(Wl(ic,i),cls), size(Lk));
        if isg, W{i,j} = gpuArray(W{i,j}); end
    end
end
for i = 1:Ns+1, W{i,1} = zeros(size(Lk),'like',Ej{1}); end
end

function Lm = baryL(s, tn, wb)
n = numel(tn); Lm = zeros(numel(s), n);
D = s(:) - tn;                         % nq x n
exact = abs(D) < 1e-300;
D(exact) = 1;
T = wb./D;
Lm = T./sum(T,2);
[r,c] = find(exact);
for q = 1:numel(r), Lm(r(q),:) = 0; Lm(r(q),c(q)) = 1; end
end

function [x,w] = gauleg(n)
b = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
[V,Dg] = eig(diag(b,1)+diag(b,-1));
x = diag(Dg); w = 2*V(1,:)'.^2;
[x,ix] = sort(x); w = w(ix);
end
