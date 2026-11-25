function v = sp_r2r(v, dim, type, M, dir)
%SP_R2R  Real-to-real trigonometric transform along dimension DIM.
%   type 0 : DCT-I on points j=0..M      (even channel), coefficients k=0..M
%   type 1 : DST-I on points j=1..M-1    (odd channel),  coefficients k=1..M-1
%            (entries j=0,M resp. k=0,M are not degrees of freedom; they are
%             carried as structural zeros so that all channels share the
%             common (M+1)^d grid used for pointwise products.)
%   dir = +1 forward (values -> coefficients), -1 inverse.
% The transform of a channel depends on that channel only; hence rounding
% errors can never leak from one isotypic component into another.
if type == 0
    ext = cat(dim, v, flip(sub(v,dim,2:M),dim));
    if dir > 0
        V = fft(ext,[],dim);
    else
        V = ifft(ext,[],dim);
    end
    v = real(sub(V,dim,1:M+1));
else
    s   = sub(v,dim,2:M);
    z   = zeros(size(sub(v,dim,1)),'like',v);
    ext = cat(dim, z, s, z, -flip(s,dim));
    if dir > 0
        V = fft(ext,[],dim);
        v = cat(dim, z, -imag(sub(V,dim,2:M)), z);
    else
        V = ifft(ext,[],dim);
        v = cat(dim, z, imag(sub(V,dim,2:M)), z);
    end
end
end

function B = sub(A, dim, ix)
c = repmat({':'},1,max(ndims(A),dim));
c{dim} = ix;
B = A(c{:});
end
