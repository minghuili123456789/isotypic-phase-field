function w = iso_transform_mat(v, P, dir)
%ISO_TRANSFORM_MAT  Channel-wise DCT-I/DST-I by dense matrix products
% (one GEMM per direction and parity group).  Rounding errors remain
% channel-internal, so Theorems 4.2/4.3 apply unchanged.
persistent cache
d = P.d; M = P.M; nc = P.nc;
key = sprintf('%d_%s', M, class(v));
if isempty(cache) || ~isfield(cache, ['m' key])
    j = 0:M; [J,K] = ndgrid(j,j);
    Cs = cos(pi*mod(J.*K,2*M)/M);           % (k,j)
    wv = 2*ones(1,M+1); wv([1 end]) = 1;
    Cf = Cs.*wv;  Ci = (Cs.*wv)/(2*M);       % forward (k<-j), inverse (j<-k)
    Ss = sin(pi*mod(J.*K,2*M)/M); Ss([1 end],:) = 0; Ss(:,[1 end]) = 0;
    Sf = 2*Ss;    Si = Ss'/M;
    s.Cf = cast(Cf,class(v)); s.Ci = cast(Ci,class(v));
    s.Sf = cast(Sf,class(v)); s.Si = cast(Si,class(v));
    cache.(['m' key]) = s;
end
s = cache.(['m' key]);
if dir > 0, A0 = s.Cf; A1 = s.Sf; else, A0 = s.Ci; A1 = s.Si; end
w = v;
act = true(1,nc); if isfield(P,'active'), act = P.active; end
sz = size(w); sz(end+1:d+1) = 1;
for i = 1:d
    perm = [i, setdiff(1:d+1, i)];
    X = permute(w, perm);
    szX = size(X); szX(end+1:d+1) = 1;
    X = reshape(X, M+1, [], nc);
    for c = find(act)
        if bitget(c-1,i), A = A1; else, A = A0; end
        X(:,:,c) = A*X(:,:,c);       % fixed-shape GEMM per channel
    end
    w = ipermute(reshape(X, szX), perm);
end
w = reshape(w, sz);
end
