% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Full 16-channel (non-interventional!) A_R x T_2 isotypic ESP for Wang (2023)
% 2D Example 4.1: kappa = 0.1, u0 = sin(2x)*sin(2y), Np = 256, S = 1, Te = 1, Ns = 10, Niter = 20.
% 1) First verifies ~1e-15 algebraic equivalence with standard FFT on random non-symmetric 2D data.
% 2) Then runs the EXACT same Wang (2023) Example 4.1 ESP parameters to T = 120 and appends
%    to results/exp_wang41_2d.mat.
here = fileparts(mfilename('fullpath')); addpath(here);
resfile = fullfile(here,'..','..','results','exp_wang41_2d.mat');

% --- Step 1: Verify algebraic equivalence on random non-symmetric 2D data (N=32) ---
Nt = 32; kap = 0.1; S_stab = 1;
Pf_t = sp_setup(struct('d',2,'N',Nt,'eps2',kap^2,'gam',1,'S',S_stab,'rep','fourier'));
Of_t = sp_ops(Pf_t); If_t = sp_integrator(Pf_t, Of_t, 'etdrk2', 0.05);
[Pg_t, Og_t] = make_art2_2d(Nt, kap^2, 1, S_stab);
Ig_t = sp_integrator(Pg_t, Og_t, 'etdrk2', 0.05);
rng(42); u_ns = 0.4*randn(Nt, Nt);
U_ns = split_art2_2d(u_ns);
fprintf('Split/assemble round-trip err (N=%d): %.2e\n', Nt, max(abs(u_ns - assem_art2_2d(U_ns)),[],'all'));
uhf = Of_t.toSpec(u_ns); uhg = Og_t.toSpec(U_ns);
for k = 1:20, uhf = If_t.step(uhf); uhg = Ig_t.step(uhg); end
eq_err = max(abs(Of_t.toPhys(uhf) - assem_art2_2d(Og_t.toPhys(uhg))),[],'all');
fprintf('16-channel A_R x T_2 non-symmetric 2D equivalence (20 ETDRK2 steps): %.2e\n', eq_err);

% --- Step 2: Run Wang (2023) 2D Example 4.1 (N=256, EXACT same ESP params, T=120) ---
if exist('EQ_ONLY','var') && EQ_ONLY, return; end
S = load(resfile);
N = S.N; Te = S.Te; Ns = S.Ns; Nit = S.Nit; T = S.T; tsnaps = S.tsnaps;
[Pg, Og] = make_art2_2d(N, kap^2, 1, S_stab);
Ig = sp_integrator(Pg, Og, 'esp', Te, struct('Ns',Ns,'Niter',Nit));

% Exact initial datum u0(x,y) = sin(2x)*sin(2y) -> channel c = (rx=1,tx=0,ry=1,ty=0) = 1 + 4 = 5 (0-based) => index 6!
Q = N/4; xq = (0:Q)'*2*pi/N; [Xq, Yq] = ndgrid(xq, xq);
U0_g = zeros(Q+1, Q+1, 16);
U0_g(:,:,6) = sin(2*Xq).*sin(2*Yq); % exact odd-x, pi-per-x, odd-y, pi-per-y channel
% Enforce exact boundary zeros of DST-I
U0_g([1 end], :, 6) = 0; U0_g(:, [1 end], 6) = 0;

Pf = sp_setup(struct('d',2,'N',N,'eps2',kap^2,'gam',1,'S',S_stab,'rep','fourier'));
Of = sp_ops(Pf);
dx = 2*pi/N;
calcE = @(u, uh) 0.5*kap^2*(dx^2/N^2)*sum(Pf.K2(:).*abs(uh(:)).^2) + 0.25*dx^2*sum((u(:).^2 - 1).^2);

u0_check = assem_art2_2d(U0_g);
fprintf('Initial diff vs sin(2x)sin(2y) on N=256 grid: %.2e\n', ...
    max(abs(u0_check - sin(2*(0:N-1)'*dx).*sin(2*(0:N-1)*dx)),[],'all'));

uh_g = Og.toSpec(U0_g);
D_art2 = zeros(T+1,1); E_art2 = zeros(T+1,1);
E_art2(1) = calcE(u0_check, Of.toSpec(u0_check));
snap_art2 = cell(numel(tsnaps),1);
trans_idx = setdiff(1:16, 6); % all 15 transverse channels!

tic;
for t = 1:T
    uh_g = Ig.step(uh_g);
    Vg   = Og.toPhys(uh_g);
    D_art2(t+1) = max(abs(Vg(:,:,trans_idx)),[],'all');
    ug = assem_art2_2d(Vg);
    E_art2(t+1) = calcE(ug, Of.toSpec(ug));
    if mod(t,20)==0 || t==1 || t==5
        fprintf('  t=%3d | A_R x T_2 max transverse D = %g | E = %.4f (%.1fs)\n', ...
            t, D_art2(t+1), E_art2(t+1), toc);
    end
    k = find(tsnaps == t, 1);
    if ~isempty(k), snap_art2{k} = ug; end
end

% Also compute D_t2 (pi-periodicity defect) of the A_R-only run at snapshots!
D_ar_t2 = zeros(numel(tsnaps),1);
for k = 1:numel(tsnaps)
    u_ar = S.snap_iso{k};
    D_ar_t2(k) = 0.5*max(max(abs(u_ar - circshift(u_ar,[N/2 0])),[],'all'), ...
                         max(abs(u_ar - circshift(u_ar,[0 N/2])),[],'all'));
end
fprintf('A_R-only pi-periodicity defect at t=[%s]: %s\n', ...
    num2str(tsnaps), sprintf('%.2e ', D_ar_t2));

save(resfile, 'D_art2','E_art2','snap_art2','eq_err','D_ar_t2', '-append');
fprintf('Appended A_R x T_2 results to exp_wang41_2d.mat\n');

% =========================================================================
function [P, O] = make_art2_2d(N, eps2, gam, S)
% 16-channel A_R x T_2 (Z_2^4) representation on (Q+1)x(Q+1), Q = N/4
Q = N/4;
P.d = 2; P.N = N; P.Q = Q; P.nc = 16; P.eps2 = eps2; P.gam = gam; P.S = S;
P.rep = 'art2';
% Wavenumber table for each 1D channel c1 in {0,1,2,3}:
%  c1=0 (r=0,t=0): even k = 0,2,...,2Q
%  c1=1 (r=1,t=0): even k = 0,2,...,2Q (0 at endpoints)
%  c1=2 (r=0,t=1): odd  k = 1,3,...,2Q-1, (dummy 0 at Q+1)
%  c1=3 (r=1,t=1): odd  k = 1,3,...,2Q-1, (dummy 0 at 1)
ke = (0:Q)' * 2;
ko_c = [2*(0:Q-1)'+1; 0]; % for c1=2 (zero at last entry)
ko_s = [0; 2*(0:Q-1)'+1]; % for c1=3 (zero at first entry)
K1 = {ke, ke, ko_c, ko_s};
P.K2 = zeros(Q+1, Q+1, 16);
for c = 0:15
    cx = bitand(c, 3);       % rx + 2*tx
    cy = bitshift(c, -2);    % ry + 2*ty
    kx = K1{cx+1}; ky = K1{cy+1};
    P.K2(:,:,c+1) = kx.^2 + (ky.').^2;
end
P.Lk = eps2 * P.K2 + S;

O.toSpec = @(U) xform_art2_2d(U, +1);
O.toPhys = @(Uh) xform_art2_2d(Uh, -1);
O.cube   = @(U) cube_z2(U, 16);
O.N      = @(Uh) O.toSpec((gam+S)*O.toPhys(Uh) - gam*O.cube(O.toPhys(Uh)));
end

function W = xform_art2_2d(V, dir)
% Separable 2D transform over the 16 A_R x T_2 channels
Q = size(V,1) - 1;
W = zeros(size(V),'like',V);
for c = 0:15
    cx = bitand(c, 3);
    cy = bitshift(c, -2);
    A = V(:,:,c+1);
    if ~any(A(:)), continue; end % exact zero-channel preservation & huge speedup!
    A = r2r_g1(A, cx, dir);      % transform along dim 1
    A = r2r_g1(A.', cy, dir).';  % transform along dim 2
    W(:,:,c+1) = A;
end
end

function Y = r2r_g1(X, c1, dir)
% Exact 1D G1 = <R, T_pi> channel transform along dim 1 (columns of X)
% Normalized so forward followed by inverse is the exact identity!
Q = size(X,1) - 1; m = size(X,2);
switch c1
    case 0 % cos(2m x): DCT-I of length Q+1 via 2Q-point FFT
        ext = [X; X(Q:-1:2, :)];
        if dir > 0
            F = real(fft(ext, [], 1)); Y = F(1:Q+1, :);
        else
            F = real(ifft(ext, [], 1)); Y = F(1:Q+1, :);
        end
    case 1 % sin(2m x): DST-I of length Q+1 (zeros at 1 and Q+1)
        Z = zeros(1, m, 'like', X);
        ext = [Z; X(2:Q, :); Z; -X(Q:-1:2, :)];
        if dir > 0
            F = -imag(fft(ext, [], 1)); Y = F(1:Q+1, :);
        else
            F =  imag(ifft(ext, [], 1)); Y = F(1:Q+1, :);
        end
        Y(1,:) = 0; Y(Q+1,:) = 0;
    case 2 % cos((2m+1)x): zero at j=Q (row Q+1); spectrum at odd indices 2:2:2Q -> stored in rows 1..Q
        Z = zeros(1, m, 'like', X);
        if dir > 0
            ext = [X(1:Q,:); Z; -X(Q:-1:2,:); -X(1:Q,:); Z; X(Q:-1:2,:)];
            F = real(fft(ext, [], 1));
            Y = [F(2:2:2*Q, :); Z];
        else
            F = zeros(4*Q, m, 'like', X);
            F(2:2:2*Q, :)       = X(1:Q, :);
            F(4*Q:-2:2*Q+2, :)  = X(1:Q, :);
            v = real(ifft(F, [], 1));
            Y = [v(1:Q, :); Z];
        end
    case 3 % sin((2m+1)x): zero at j=0 (row 1); spectrum at odd indices 2:2:2Q -> stored in rows 2..Q+1
        Z = zeros(1, m, 'like', X);
        if dir > 0
            ext = [Z; X(2:Q+1,:); X(Q:-1:2,:); Z; -X(2:Q+1,:); -X(Q:-1:2,:)];
            F = -imag(fft(ext, [], 1));
            Y = [Z; F(2:2:2*Q, :)];
        else
            F = zeros(4*Q, m, 'like', X);
            F(2:2:2*Q, :)       = -1i * X(2:Q+1, :);
            F(4*Q:-2:2*Q+2, :)  =  1i * X(2:Q+1, :);
            v = real(ifft(F, [], 1));
            Y = [Z; v(2:Q+1, :)];
        end
end
end

function U3 = cube_z2(U, nc)
sz = size(U); A = reshape(U, [], nc);
nz = any(A ~= 0, 1);
A2 = zeros(size(A),'like',A); A3 = zeros(size(A),'like',A);
for a = find(nz)-1
    for b = find(nz)-1
        r = bitxor(a, b);
        A2(:,r+1) = A2(:,r+1) + A(:,a+1).*A(:,b+1);
    end
end
nz2 = any(A2 ~= 0, 1);
for r = find(nz2)-1
    for b = find(nz)-1
        s = bitxor(r, b);
        A3(:,s+1) = A3(:,s+1) + A2(:,r+1).*A(:,b+1);
    end
end
U3 = reshape(A3, sz);
end

function U = split_art2_2d(u)
% Split full N x N grid into 16 (Q+1)x(Q+1) channels of A_R x T_2
N = size(u,1); Q = N/4;
% First split along dim 1 into 4 channels -> (Q+1) x N x 4
V = zeros(Q+1, N, 4);
refl = [1, N:-1:2];
u_e = 0.5*(u + u(refl,:)); u_o = 0.5*(u - u(refl,:));
% Now split each by T_pi (shift by N/2 = 2Q)
V(:,:,1) = 0.5*(u_e(1:Q+1,:) + u_e(2*Q+1:3*Q+1,:)); % rx=0,tx=0
V(:,:,2) = 0.5*(u_o(1:Q+1,:) + u_o(2*Q+1:3*Q+1,:)); % rx=1,tx=0
V(:,:,3) = 0.5*(u_e(1:Q+1,:) - u_e(2*Q+1:3*Q+1,:)); % rx=0,tx=1
V(:,:,4) = 0.5*(u_o(1:Q+1,:) - u_o(2*Q+1:3*Q+1,:)); % rx=1,tx=1
% Now split along dim 2 into 4 channels -> (Q+1) x (Q+1) x 16
U = zeros(Q+1, Q+1, 16);
for cx = 0:3
    A = V(:,:,cx+1);
    Ae = 0.5*(A + A(:,refl)); Ao = 0.5*(A - A(:,refl));
    U(:,:,cx + 0*4 + 1) = 0.5*(Ae(:,1:Q+1) + Ae(:,2*Q+1:3*Q+1));
    U(:,:,cx + 1*4 + 1) = 0.5*(Ao(:,1:Q+1) + Ao(:,2*Q+1:3*Q+1));
    U(:,:,cx + 2*4 + 1) = 0.5*(Ae(:,1:Q+1) - Ae(:,2*Q+1:3*Q+1));
    U(:,:,cx + 3*4 + 1) = 0.5*(Ao(:,1:Q+1) - Ao(:,2*Q+1:3*Q+1));
end
end

function u = assem_art2_2d(U)
Q = size(U,1) - 1; N = 4*Q;
u = zeros(N, N, 'like', U);
for c = 0:15
    A = U(:,:,c+1);
    if ~any(A(:)), continue; end
    cx = bitand(c, 3); cy = bitshift(c, -2);
    A = ext_g1(A, cx);
    A = ext_g1(A.', cy).';
    u = u + A;
end
end

function E = ext_g1(X, c1)
Q = size(X,1) - 1;
switch c1
    case 0, E = [X(1:Q+1,:); X(Q:-1:2,:); X(1:Q+1,:); X(Q:-1:2,:)];
    case 1, E = [X(1:Q+1,:); -X(Q:-1:2,:); X(1:Q+1,:); -X(Q:-1:2,:)];
    case 2, E = [X(1:Q+1,:); -X(Q:-1:2,:); -X(1:Q+1,:); X(Q:-1:2,:)];
    case 3, E = [X(1:Q+1,:); X(Q:-1:2,:); -X(1:Q+1,:); -X(Q:-1:2,:)];
end
end
