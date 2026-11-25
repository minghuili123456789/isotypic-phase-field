function U = sp_channels(u, P)
%SP_CHANNELS  Isotypic channels (fundamental grid, (M+1)^d x 2^d) of a
% physical state u in representation P.  For 'iso' this is the identity;
% for 'fourier' the full-grid field is split by reflection butterflies.
d = P.d; M = P.M; nc = 2^d; N = P.N;
if strcmp(P.rep,'iso'), U = u; return; end
vals = cell(1,nc); vals{1} = u;
refl = [1, N:-1:2];                    % j -> -j mod N
for i = 1:d
    nv = cell(1,nc);
    for g = 0:2^(i-1)-1
        p = vals{g+1};
        cc = repmat({':'},1,d); cc{i} = refl;
        q = p(cc{:});
        nv{g+1}         = (p + q)/2;
        nv{g+1+2^(i-1)} = (p - q)/2;
    end
    vals(1:2^i) = nv(1:2^i);
end
cc = repmat({1:M+1},1,d);
for c = 1:nc, vals{c} = vals{c}(cc{:}); end
U = cat(d+1, vals{:});
end
