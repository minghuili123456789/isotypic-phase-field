function O = sp_t2_ops(P)
%SP_T2_OPS  Transforms and Z_2^d-graded nonlinear operator for T2 isotypic
% representation.  Same algebra as reflection isotypic, but transforms are
% standard M-point FFTs (with twiddle exp(-i*pi*j/M) on odd-parity axes).
d = P.d; M = P.M; nc = P.nc; g = P.gam; S = P.S;
O.toSpec = @(u) t2_xform(u, P, +1);
O.toPhys = @(uh) t2_xform(uh, P, -1);
if isfield(P,'active') && ~all(P.active)
    act = find(P.active) - 1;
    O.cube = @(u) cube_red(u, nc, act);
else
    O.cube = @(u) cube_full(u, nc);
end
if isfield(P,'F') && ~isempty(P.F)
    O.N = @(uh) O.toSpec(S*O.toPhys(uh) - g*t2_orbit(O.toPhys(uh), d, nc, P.F));
elseif strcmp(P.model,'CAC')
    O.N = @(uh) cac_nlin(uh, O, g, S);
else
    O.N = @(uh) O.toSpec((g+S)*O.toPhys(uh) - g*O.cube(O.toPhys(uh)));
end
end

function Nh = cac_nlin(uh, O, g, S)
u = O.toPhys(uh); F = O.toSpec(g*O.cube(u) - g*u);
F(1) = 0; Nh = S*uh - F;
end

function w = t2_xform(v, P, dir)
d = P.d; M = P.M; nc = P.nc; tw = P.tw;
if dir > 0, w = complex(zeros(size(v),'like',v)); else, w = zeros(size(v),'like',real(v)); end
colons = repmat({':'},1,d);
for c = 1:nc
    if isfield(P,'active') && ~P.active(c), continue; end
    vc = v(colons{:},c);
    for i = 1:d
        shp = ones(1,max(ndims(vc),i)); shp(i) = M;
        if P.par(c,i)
            if dir > 0
                vc = fft(vc .* reshape(tw, shp), [], i);
            else
                vc = ifft(vc, [], i) .* reshape(conj(tw), shp);
            end
        else
            if dir > 0, vc = fft(vc, [], i); else, vc = ifft(vc, [], i); end
        end
    end
    if dir < 0, vc = real(vc); end
    w(colons{:},c) = vc;
end
end

function u3 = cube_full(u, nc)
sz = size(u); U = reshape(u, [], nc);
U2 = zeros(size(U),'like',U); U3 = zeros(size(U),'like',U);
for r = 0:nc-1
    acc = U(:,1).*U(:,bitxor(0,r)+1);
    for a = 1:nc-1, acc = acc + U(:,a+1).*U(:,bitxor(a,r)+1); end
    U2(:,r+1) = acc;
end
for s = 0:nc-1
    acc = U2(:,1).*U(:,bitxor(s,0)+1);
    for r = 1:nc-1, acc = acc + U2(:,r+1).*U(:,bitxor(s,r)+1); end
    U3(:,s+1) = acc;
end
u3 = reshape(U3, sz);
end

function u3 = cube_red(u, nc, act)
sz = size(u); U = reshape(u, [], nc);
U2 = zeros(size(U),'like',U); U3 = zeros(size(U),'like',U);
isact = false(1,nc); isact(act+1) = true;
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

function Fc = t2_orbit(u, d, nc, F)
% Same Walsh-Hadamard orbit evaluation on the 2^d half-period shift orbit
sz = size(u); V = reshape(u, [], nc);
for i = 1:d
    W = V;
    for c = 0:nc-1
        if ~bitget(c,i)
            p = V(:,c+1); q = V(:,c+1+2^(i-1));
            W(:,c+1) = p + q; W(:,c+1+2^(i-1)) = p - q;
        end
    end
    V = W;
end
V = F(V);
for i = 1:d
    W = V;
    for c = 0:nc-1
        if ~bitget(c,i)
            p = V(:,c+1); q = V(:,c+1+2^(i-1));
            W(:,c+1) = (p + q)/2; W(:,c+1+2^(i-1)) = (p - q)/2;
        end
    end
    V = W;
end
Fc = reshape(V, sz);
end
