function [u, x] = sp_init(P, fun)
%SP_INIT  Initial data in representation P from a function handle
%   fun(X1,...,Xd) evaluated on the periodic grid x_j = j*h (j=0..N-1),
%   identified with [-L/2,L/2) by periodicity.
% The isotypic channels are obtained on the fundamental grid [0,L/2]^d by
% Walsh--Hadamard butterflies of fun evaluated at x and at -x (exact
% negation).  The full-grid (Fourier) data are *assembled* from the
% channels, so that data with an exact floating-point symmetry are handed
% bitwise-symmetric to BOTH representations (fair comparison).
d = P.d; M = P.M; h = P.h; nc = 2^d;
xh = (0:M)'*h;                      % fundamental grid [0,L/2]
switch d
    case 1, X = {xh};
    case 2, [a,b] = ndgrid(xh,xh); X = {a,b};
    case 3, [a,b,c] = ndgrid(xh,xh,xh); X = {a,b,c};
end
% values on the 2^d sign-orbits
vals = cell(1,nc);
for g = 0:nc-1
    Y = X;
    for i = 1:d
        if bitget(g,i), Y{i} = -X{i}; end
    end
    vals{g+1} = fun(Y{:});
end
% butterflies: dimension by dimension
for i = 1:d
    nv = vals;
    for g = 0:nc-1
        if ~bitget(g,i)
            p = vals{g+1}; q = vals{g+1+2^(i-1)};
            nv{g+1}          = (p + q)/2;     % even in x_i
            nv{g+1+2^(i-1)}  = (p - q)/2;     % odd  in x_i
        end
    end
    vals = nv;
end
% odd channels vanish on x_i in {0,L/2}: these are not degrees of freedom
for c = 0:nc-1
    for i = 1:d
        if bitget(c,i)
            vals{c+1} = zero_bnd(vals{c+1}, i, M);
        end
    end
end
U = cat(d+1, vals{:});
U = cast(U, P.prec);
if strcmp(P.rep,'iso')
    u = U;
else
    u = sp_assemble(U, d, M);
end
if P.gpu, u = gpuArray(u); end
x = xh;
end

function A = zero_bnd(A, i, M)
c = repmat({':'},1,ndims(A)); c{i} = [1 M+1];
A(c{:}) = 0;
end
