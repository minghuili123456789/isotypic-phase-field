function uf = sp_assemble(U, d, M)
%SP_ASSEMBLE  Full periodic grid values (N^d) from isotypic channels U on the
% fundamental grid ((M+1)^d x 2^d).  Index j in [0,N) maps to |j| = j (j<=M)
% or N-j (j>M), with sign (-1)^{sigma_i} for the reflected half.
N = 2*M; nc = 2^d;
map = [1:M+1, M:-1:2];                 % index of |x_j| on fundamental grid
sgn = [ones(1,M+1), -ones(1,M-1)];      % parity sign of the reflected half
colons = repmat({':'},1,d);
uf = [];
for c = 0:nc-1
    A = U(colons{:},c+1);
    for i = 1:d
        cc = repmat({':'},1,d); cc{i} = map;
        A = A(cc{:});
        if bitget(c,i)
            shp = ones(1,d); shp(i) = N;
            A = A .* reshape(cast(sgn,'like',A), [shp 1]);
        end
    end
    if isempty(uf), uf = A; else, uf = uf + A; end
end
end
