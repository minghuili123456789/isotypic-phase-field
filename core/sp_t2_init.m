function [U, uf] = sp_t2_init(P, fun, active_par)
%SP_T2_INIT  Initial channels on [0,L/2)^d for T2 isotypic representation.
% If active_par is omitted, evaluates fun on the 2^d half-period shift
% orbits x + g*(L/2) and splits via Walsh-Hadamard butterflies (general data).
% If numeric active_par (0-based channel index, e.g. 1 in 1D, 3 in 2D odd-odd,
% 7 in 3D odd^3) is given, initializes that single analytic parity channel on
% [0,L/2)^d with exact floating-point zeros in all transverse channels (avoiding
% O(u_mach) sampling error from floating-point evaluation of x + pi).
d = P.d; M = P.M; h = P.h; L = P.L; nc = P.nc;
xh = (0:M-1)'*h;
switch d
    case 1, X = {xh};
    case 2, [a,b] = ndgrid(xh,xh); X = {a,b};
    case 3, [a,b,c] = ndgrid(xh,xh,xh); X = {a,b,c};
end
if nargin >= 3 && isa(active_par, 'numeric')
    % Analytic anti-periodic benchmark: active_par is the 0-based channel index
    % (e.g. 1 in 1D, 3 in 2D odd-odd, 7 in 3D odd^3) — fun on [0,pi)^d defines
    % that channel and the other channels are identically 0 by periodicity
    % f(x + pi e_i) = (-1)^{par_i} f(x).  We also assemble full-grid uf so
    % standard Fourier gets the exact same bitwise anti-periodic vector!
    U = zeros([repmat(M,1,d) nc], P.prec);
    colons = repmat({':'},1,d);
    U(colons{:}, active_par+1) = cast(fun(X{:}), P.prec);
else
    % General data: evaluate on full 2^d shift orbit and split by butterflies
    vals = cell(1,nc);
    for g = 0:nc-1
        Y = X;
        for i = 1:d
            if bitget(g,i), Y{i} = X{i} + L/2; end
        end
        vals{g+1} = fun(Y{:});
    end
    for i = 1:d
        nv = vals;
        for g = 0:nc-1
            if ~bitget(g,i)
                p = vals{g+1}; q = vals{g+1+2^(i-1)};
                nv{g+1} = (p + q)/2; nv{g+1+2^(i-1)} = (p - q)/2;
            end
        end
        vals = nv;
    end
    U = cast(cat(d+1, vals{:}), P.prec);
end
if nargout >= 2
    uf = sp_t2_assemble(U, d, M);
end
end
