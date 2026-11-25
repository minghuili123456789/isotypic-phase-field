function uf = sp_t2_assemble(U, d, M)
%SP_T2_ASSEMBLE  Full N^d grid values from T2 isotypic channels (M^d x 2^d).
N = 2*M; nc = 2^d;
map = [1:M, 1:M];
sgn = [ones(1,M), -ones(1,M)];
colons = repmat({':'},1,d);
uf = [];
for c = 0:nc-1
    A = U(colons{:},c+1);
    for i = 1:d
        cc = repmat({':'},1,d); cc{i} = map; A = A(cc{:});
        if bitget(c,i)
            shp = ones(1,d); shp(i) = N;
            A = A .* reshape(cast(sgn,'like',A), [shp 1]);
        end
    end
    if isempty(uf), uf = A; else, uf = uf + A; end
end
end
