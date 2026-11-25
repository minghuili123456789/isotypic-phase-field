function O = sp_ops(P)
%SP_OPS  Return function handles for the spectral <-> physical transforms and
% the nonlinear map for representation P.rep.
%   O.toSpec(u), O.toPhys(uh), O.cube(u), O.N(uh) (stabilised nonlinear term
%   in spectral space:  N(u) = (gam+S) u - gam u^3 ).
d = P.d; M = P.M; nc = P.nc;
if strcmp(P.rep,'fourier')
    if d == 1
        O.toSpec = @(u) fft(u);
        O.toPhys = @(uh) real(ifft(uh));
    else
        O.toSpec = @(u) fftn(u);
        O.toPhys = @(uh) real(ifftn(uh));
    end
    O.cube = @(u) u.*u.*u;
else
    O.toSpec = @(u) iso_transform(u, P, +1);
    O.toPhys = @(uh) iso_transform(uh, P, -1);
    if isfield(P,'xform') && strcmp(P.xform,'mat')
        O.toSpec = @(u) iso_transform_mat(u, P, +1);
        O.toPhys = @(uh) iso_transform_mat(uh, P, -1);
    end
    if isfield(P,'active') && ~all(P.active)
        O.cube = @(u) iso_cube_red(u, d, nc, find(P.active)-1);
    else
        O.cube   = @(u) iso_cube(u, d, nc);
    end
end
g = P.gam; S = P.S;
if isfield(P,'F') && ~isempty(P.F)
    % general (odd, odd-exact) nonlinearity F:  N(u) = S u - gam F(u)
    if strcmp(P.rep,'fourier')
        O.Fval = @(u) P.F(u);
    else
        O.Fval = @(u) iso_orbit(u, d, nc, P.F);
    end
    O.N = @(uh) O.toSpec(S*O.toPhys(uh) - g*O.Fval(O.toPhys(uh)));
    return
end
if strcmp(P.model,'CAC')
    O.N = @(uh) cac_N(uh, O, g, S, P);
else
    O.N = @(uh) nlin(uh, O, g, S);
end
end

function Nh = nlin(uh, O, g, S)
u  = O.toPhys(uh);
Nh = O.toSpec((g+S)*u - g*O.cube(u));
end

function Nh = cac_N(uh, O, g, S, P)
% conservative Allen--Cahn: remove the spatial mean of f(u) (only the k=0
% coefficient of the totally even channel carries a mean).
u  = O.toPhys(uh);
F  = O.toSpec(g*O.cube(u) - g*u);       % gam*f(u)
F  = zero_mean(F, P);
Nh = S*uh - F;
end

function F = zero_mean(F, P)
if strcmp(P.rep,'fourier')
    F(1) = 0;
else
    F(1) = 0;   % linear index 1 = (k=0,...,k=0) of channel 0 (all even)
end
end

function w = iso_transform(v, P, dir)
d = P.d; M = P.M; nc = P.nc;
w = v;
colons = repmat({':'},1,d);
for c = 1:nc
    if isfield(P,'active') && ~P.active(c), continue; end   % structurally zero channel
    vc = v(colons{:},c);
    for i = 1:d
        vc = sp_r2r(vc, i, P.sig(c,i), M, dir);
    end
    w(colons{:},c) = vc;
end
end

function u3 = iso_cube(u, d, nc)
% Channel algebra of the cubic nonlinearity in the Z_2^d grading:
%   (u^2)_r = sum_a u_a u_{a xor r},   (u^3)_s = sum_r (u^2)_r u_{s xor r}.
% A product containing a vanishing channel is an exact (structural) zero.
sz = size(u);
U  = reshape(u, [], nc);
U2 = zeros(size(U),'like',U);
U3 = zeros(size(U),'like',U);
for r = 0:nc-1
    acc = U(:,1).*U(:,bitxor(0,r)+1);
    for a = 1:nc-1
        acc = acc + U(:,a+1).*U(:,bitxor(a,r)+1);
    end
    U2(:,r+1) = acc;
end
for s = 0:nc-1
    acc = U2(:,1).*U(:,bitxor(s,0)+1);
    for r = 1:nc-1
        acc = acc + U2(:,r+1).*U(:,bitxor(s,r)+1);
    end
    U3(:,s+1) = acc;
end
u3 = reshape(U3, sz);
end

function u3 = iso_cube_red(u, d, nc, act)
% reduced channel algebra: only the active (closed) channel set is computed;
% omitted products are exact zeros in the full algorithm (Theorem 4.2), so
% the result is bitwise identical whenever inactive channels vanish.
sz = size(u);
U  = reshape(u, [], nc);
U3 = zeros(size(U),'like',U);
isact = false(1,nc); isact(act+1) = true;
U2 = zeros(size(U),'like',U);
for r = 0:nc-1
    acc = []; 
    for a = 0:nc-1
        b = bitxor(a,r);
        if isact(a+1) && isact(b+1)
            t = U(:,a+1).*U(:,b+1);
            if isempty(acc), acc = t; else, acc = acc + t; end
        end
    end
    if ~isempty(acc), U2(:,r+1) = acc; end
end
for s = act
    acc = [];
    for r = 0:nc-1
        b = bitxor(s,r);
        if isact(b+1) && any(U2(:,r+1))
            t = U2(:,r+1).*U(:,b+1);
            if isempty(acc), acc = t; else, acc = acc + t; end
        end
    end
    if ~isempty(acc), U3(:,s+1) = acc; end
end
u3 = reshape(U3, sz);
end
