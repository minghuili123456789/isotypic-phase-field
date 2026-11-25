function B = iso_spectral_backend(d, N, mode)
%ISO_SPECTRAL_BACKEND  Drop-in floating-point-equivariant spectral backend
% for MATLAB PDE solvers (ode45, ode15s, IMEX, ETDRK, CN, SAV, ESP, SBDF, SDC).
%
% Usage:
%   B = iso_spectral_backend(1, 128, 'AR');   % 1D reflection (even/odd channels)
%   B = iso_spectral_backend(2, 128, 'AR');   % 2D reflection (4 channels)
%   B = iso_spectral_backend(1, 128, 'T2');   % 1D half-period parity (arbitrary shifts)
%
% Returned struct B provides the exact same 3 primitives as a standard FFT solver:
%   B.K2          -- Laplacian eigenvalues |k|^2 matching channel array shape
%   B.toSpec(U)   -- channel-wise forward spectral transform
%   B.toPhys(Uh)  -- channel-wise inverse spectral transform
%   B.cube(U)     -- exact Z_2^d-graded cubic nonlinearity (u^3)_sigma
%   B.orbit(U, F) -- odd-exact Walsh-Hadamard orbit evaluation for general odd F(u)
%   B.split(u)    -- split full N^d periodic grid field into 2^d channels
%   B.assemble(U) -- reassemble full N^d periodic grid field from 2^d channels
%
% Algebraic identity: on ANY non-symmetric field u, B reproduces standard
% Fourier pseudo-spectral operators to ~1e-15 while keeping every symmetric
% subspace bitwise invariant (D === 0) without filters or projections.
if nargin < 3, mode = 'AR'; end
M = N/2; nc = 2^d;
B.d = d; B.N = N; B.M = M; B.nc = nc; B.mode = upper(mode);

switch B.mode
    case 'AR'
        % Fundamental grid [0, pi]^d of size (M+1)^d x 2^d (DCT-I / DST-I)
        k1 = (0:M)';
        if d == 1
            B.K2 = repmat(k1.^2, 1, 2);
        elseif d == 2
            B.K2 = repmat(k1.^2 + (k1.').^2, 1, 1, 4);
        else
            error('For d=3 see core/sp_setup.m and core/sp_ops.m.');
        end
        B.toSpec   = @(U)  ar_xform(U, d, M, +1);
        B.toPhys   = @(Uh) ar_xform(Uh, d, M, -1);
        B.split    = @(u)  ar_split(u, d, M);
        B.assemble = @(U)  ar_assem(U, d, M);

    case 'T2'
        % Half-box [0, pi)^d of size M^d x 2^d (even/odd-wavenumber twiddled FFT)
        j = (0:M-1)'; B.tw = exp(-1i*pi*j/M);
        ke = 2*[0:M/2-1, -M/2:-1]';
        ko = 2*[0:M/2-1, -M/2:-1]' + 1;
        K1 = {ke, ko};
        if d == 1
            B.K2 = [ke.^2, ko.^2];
        elseif d == 2
            B.K2 = zeros(M, M, 4);
            for c = 0:3
                kx = K1{bitget(c,1)+1}; ky = K1{bitget(c,2)+1};
                B.K2(:,:,c+1) = kx.^2 + (ky.').^2;
            end
        else
            error('For d=3 see core/sp_t2_setup.m and core/sp_t2_ops.m.');
        end
        B.toSpec   = @(U)  t2_xform(U, d, M, B.tw, +1);
        B.toPhys   = @(Uh) t2_xform(Uh, d, M, B.tw, -1);
        B.split    = @(u)  t2_split(u, d, M);
        B.assemble = @(U)  t2_assem(U, d, M);

    otherwise
        error('mode must be ''AR'' or ''T2''.');
end

B.cube  = @(U)    z2_cube(U, nc);
B.orbit = @(U, F) wh_orbit(U, d, nc, F);
end

% ======================== Internal helpers ===============================
function W = ar_xform(V, d, M, dir)
nc = 2^d; W = zeros(size(V),'like',V);
for c = 0:nc-1
    if d == 1
        W(:,c+1) = r2r_1d(V(:,c+1), bitget(c,1), M, dir);
    else
        A = r2r_1d(V(:,:,c+1), bitget(c,1), M, dir);
        W(:,:,c+1) = r2r_1d(A.', bitget(c,2), M, dir).';
    end
end
end

function Y = r2r_1d(X, isOdd, M, dir)
m = size(X,2); Z = zeros(1, m, 'like', X);
if ~isOdd
    ext = [X; X(M:-1:2,:)];
    if dir > 0, F = real(fft(ext,[],1)); else, F = real(ifft(ext,[],1)); end
    Y = F(1:M+1,:);
else
    ext = [Z; X(2:M,:); Z; -X(M:-1:2,:)];
    if dir > 0, F = -imag(fft(ext,[],1)); else, F = imag(ifft(ext,[],1)); end
    Y = F(1:M+1,:); Y(1,:) = 0; Y(M+1,:) = 0;
end
end

function U = ar_split(u, d, M)
N = 2*M; refl = [1, N:-1:2];
if d == 1
    ue = 0.5*(u + u(refl)); uo = 0.5*(u - u(refl));
    U = [ue(1:M+1), uo(1:M+1)]; U([1 end],2) = 0;
else
    ue = 0.5*(u + u(refl,:)); uo = 0.5*(u - u(refl,:));
    U = zeros(M+1, M+1, 4);
    U(:,:,1) = 0.5*(ue(1:M+1,1:M+1) + ue(1:M+1,refl(1:M+1)));
    U(:,:,2) = 0.5*(uo(1:M+1,1:M+1) + uo(1:M+1,refl(1:M+1)));
    U(:,:,3) = 0.5*(ue(1:M+1,1:M+1) - ue(1:M+1,refl(1:M+1)));
    U(:,:,4) = 0.5*(uo(1:M+1,1:M+1) - uo(1:M+1,refl(1:M+1)));
    U([1 end],:,[2 4]) = 0; U(:,[1 end],[3 4]) = 0;
end
end

function u = ar_assem(U, d, M)
map = [1:M+1, M:-1:2]; sgn = [ones(1,M+1), -ones(1,M-1)];
if d == 1
    u = U(map,1) + U(map,2).*sgn.';
else
    u = zeros(2*M, 2*M);
    for c = 0:3
        A = U(map, map, c+1);
        if bitget(c,1), A = A .* sgn.'; end
        if bitget(c,2), A = A .* sgn;   end
        u = u + A;
    end
end
end

function W = t2_xform(V, d, M, tw, dir)
nc = 2^d;
if dir > 0, W = complex(zeros(size(V))); else, W = zeros(size(V)); end
for c = 0:nc-1
    if d == 1
        vc = V(:,c+1);
        if bitget(c,1)
            if dir > 0, vc = fft(vc .* tw); else, vc = real(ifft(vc) .* conj(tw)); end
        else
            if dir > 0, vc = fft(vc);       else, vc = real(ifft(vc)); end
        end
        W(:,c+1) = vc;
    else
        vc = V(:,:,c+1);
        for ax = 1:2
            shp = [1 1]; shp(ax) = M;
            if bitget(c,ax)
                if dir > 0, vc = fft(vc .* reshape(tw,shp), [], ax);
                else,       vc = ifft(vc, [], ax) .* reshape(conj(tw),shp); end
            else
                if dir > 0, vc = fft(vc, [], ax); else, vc = ifft(vc, [], ax); end
            end
        end
        if dir < 0, vc = real(vc); end
        W(:,:,c+1) = vc;
    end
end
end

function U = t2_split(u, d, M)
if d == 1
    U = [0.5*(u(1:M)+u(M+1:2*M)), 0.5*(u(1:M)-u(M+1:2*M))];
else
    ue = 0.5*(u(1:M,:)+u(M+1:2*M,:)); uo = 0.5*(u(1:M,:)-u(M+1:2*M,:));
    U = zeros(M, M, 4);
    U(:,:,1) = 0.5*(ue(:,1:M)+ue(:,M+1:2*M));
    U(:,:,2) = 0.5*(uo(:,1:M)+uo(:,M+1:2*M));
    U(:,:,3) = 0.5*(ue(:,1:M)-ue(:,M+1:2*M));
    U(:,:,4) = 0.5*(uo(:,1:M)-uo(:,M+1:2*M));
end
end

function u = t2_assem(U, d, M)
map = [1:M, 1:M]; sgn = [ones(1,M), -ones(1,M)];
if d == 1
    u = U(map,1) + U(map,2).*sgn.';
else
    u = zeros(2*M, 2*M);
    for c = 0:3
        A = U(map, map, c+1);
        if bitget(c,1), A = A .* sgn.'; end
        if bitget(c,2), A = A .* sgn;   end
        u = u + A;
    end
end
end

function U3 = z2_cube(U, nc)
sz = size(U); A = reshape(U, [], nc);
A2 = zeros(size(A),'like',A); A3 = zeros(size(A),'like',A);
for r = 0:nc-1
    for a = 0:nc-1, A2(:,r+1) = A2(:,r+1) + A(:,a+1).*A(:,bitxor(a,r)+1); end
end
for s = 0:nc-1
    for r = 0:nc-1, A3(:,s+1) = A3(:,s+1) + A2(:,r+1).*A(:,bitxor(s,r)+1); end
end
U3 = reshape(A3, sz);
end

function Fc = wh_orbit(U, d, nc, F)
sz = size(U); V = reshape(U, [], nc);
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
            W(:,c+1) = 0.5*(p + q); W(:,c+1+2^(i-1)) = 0.5*(p - q);
        end
    end
    V = W;
end
Fc = reshape(V, sz);
end
