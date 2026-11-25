function Fc = iso_orbit(u, d, nc, F)
%ISO_ORBIT  Isotypic channels of F(u) for a general pointwise nonlinearity F
% (orbit evaluation).  Values of u on the 2^d sign-orbit points are
% synthesised from the channels by Walsh--Hadamard butterflies (p+q, p-q),
% F is applied pointwise, and the channels of F(u) are recovered by the
% butterflies ((p+q)/2, (p-q)/2).  If u lies in a symmetric subspace and F is
% odd-exact (fl(F(-v)) = -fl(F(v))), the forbidden channels of F(u) are exact
% zeros (each butterfly stage produces either 2a, 0 or a, exactly).
sz = size(u);
V = reshape(u, [], nc);
% synthesis: channel -> orbit values (orbit index g, bit i = reflect x_i)
for i = 1:d
    W = V;
    for c = 0:nc-1
        if ~bitget(c,i)
            p = V(:,c+1); q = V(:,c+1+2^(i-1));
            W(:,c+1) = p + q;            % value at x  (bit i of g = 0)
            W(:,c+1+2^(i-1)) = p - q;    % value at R_i x
        end
    end
    V = W;
end
V = F(V);
% analysis: orbit values -> channels
for i = 1:d
    W = V;
    for c = 0:nc-1
        if ~bitget(c,i)
            p = V(:,c+1); q = V(:,c+1+2^(i-1));
            W(:,c+1) = (p + q)/2;
            W(:,c+1+2^(i-1)) = (p - q)/2;
        end
    end
    V = W;
end
Fc = reshape(V, sz);   % odd channels vanish exactly on x_i in {0,L/2}
end

