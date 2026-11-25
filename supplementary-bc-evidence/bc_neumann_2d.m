% Supplementary BC Evidence 2: 2D Homogeneous Neumann X-Junction
%   du/dt = Delta u - eps^{-2}*(u^3 - u) on [0, pi]^2,  dn u = 0 on boundary
%   u0(x,y) = tanh(cos(x)*cos(y) / (sqrt(2)*eps)), eps = 0.05, N = 192 (mixed-radix!)
% Compares standard 2D DCT-I (and DCT-I + zero-mode suppression) vs.
% proposed 4-channel 2D midpoint-isotypic Neumann on [0, pi/2]^2 (D === 0 bitwise!).
clear; close all;
outfig = fullfile(fileparts(mfilename('fullpath')), 'fig_neumann_2d.png');
outpdf = fullfile(fileparts(mfilename('fullpath')), 'fig_neumann_2d.pdf');

N = 192; Q = N/2; h = pi/N; eps = 0.05; gam = 1/eps^2; S_stab = 2*gam;
tau = 1e-3; T = 2.0; nst = round(T/tau);
x = (0:N)'*h; xq = (0:Q)'*h;
[X, Y] = ndgrid(x, x); [Xq, Yq] = ndgrid(xq, xq);

% 1D DCT-I along columns (length N+1)
dct1_1d  = @(A) subsref(real(fft([A; A(end-1:-1:2,:)], [], 1)), struct('type','()','subs',{{1:N+1,':'}}));
idct1_1d = @(A) subsref(real(ifft([A; A(end-1:-1:2,:)], [], 1)), struct('type','()','subs',{{1:N+1,':'}}));
toSpec_f = @(U) dct1_1d(dct1_1d(U).').';
toPhys_f = @(Uh) idct1_1d(idct1_1d(Uh).').';
k1 = (0:N)'; Lk_f = (k1.^2 + (k1.').^2) + S_stab;

% 4-channel 2D Midpoint-Isotypic Neumann on (Q+1)x(Q+1):
%   c1=0 (even about pi/2): cos(2m x)   -> DCT-I on Q+1
%   c1=1 (odd  about pi/2): cos((2m+1)x)-> DCT-III on Q+1 (zero at row Q+1)
ke = (0:Q)'*2; ko = [2*(0:Q-1)'+1; 0]; K1 = {ke, ko};
Lk_i = zeros(Q+1, Q+1, 4);
for c = 0:3
    kx = K1{bitget(c,1)+1}; ky = K1{bitget(c,2)+1};
    Lk_i(:,:,c+1) = kx.^2 + (ky.').^2 + S_stab;
end
toSpec_i = @(U) xform2d_neum(U, +1);
toPhys_i = @(Uh) xform2d_neum(Uh, -1);

% 1. Verify non-symmetric 2D Neumann equivalence on N = 32 (20 steps)
Nt = 32; Qt = Nt/2; xt = (0:Nt)'*pi/Nt; [Xt, Yt] = ndgrid(xt, xt);
u_ns = 0.3*cos(Xt).*cos(2*Yt) + 0.25*cos(2*Xt).*cos(Yt) + 0.2*cos(3*Xt) + 0.1;
U_ns = split2d_neum(u_ns);
fprintf('2D Neumann split/assem round-trip err (N=%d): %.2e\n', Nt, max(abs(u_ns - assem2d_neum(U_ns)),[],'all'));

% 2. Symmetric 2D Neumann X-junction: odd about BOTH x=pi/2 and y=pi/2 (channel 4 = (1,1)!)
U0 = zeros(Q+1, Q+1, 4);
U0(:,:,4) = tanh(cos(Xq).*cos(Yq)/(sqrt(2)*eps));
U0(end,:,4) = 0; U0(:,end,4) = 0;
u0 = assem2d_neum(U0); % exact bitwise midpoint-odd-odd initial datum!

% ETDRK2 operators
[Ef, Q1f, Q2f] = etd2_ops(Lk_f, tau);
[Ei, Q1i, Q2i] = etd2_ops(Lk_i, tau);
Nop_f = @(uh) toSpec_f((gam+S_stab)*toPhys_f(uh) - gam*toPhys_f(uh).^3);
Nop_i = @(uh) toSpec_i((gam+S_stab)*toPhys_i(uh) - gam*cube4(toPhys_i(uh)));

uh_f = toSpec_f(u0); uh_z = uh_f; uh_i = toSpec_i(U0);
tsnaps = [0, 0.4, 1.0, 2.0]; s_steps = round(tsnaps/tau);
trec = (0:0.02:T)'; r_step = round(0.02/tau);
D_f = zeros(numel(trec),1); D_z = D_f; D_i = D_f;
Sf = cell(1,4); Si = cell(1,4);
Sf{1} = u0; Si{1} = u0;

for n = 1:nst
    % Standard DCT-I
    Nu = Nop_f(uh_f); a = Ef.*uh_f + Q1f.*Nu; uh_f = a + Q2f.*(Nop_f(a) - Nu);
    % DCT-I + zero-mode suppression
    Nu = Nop_f(uh_z); a = Ef.*uh_z + Q1f.*Nu; uh_z = a + Q2f.*(Nop_f(a) - Nu); uh_z(1,1) = 0;
    % Isotypic Neumann (4 channels)
    Nu = Nop_i(uh_i); a = Ei.*uh_i + Q1i.*Nu; uh_i = a + Q2i.*(Nop_i(a) - Nu);

    if mod(n, r_step) == 0
        idx = round(n/r_step) + 1;
        uf = toPhys_f(uh_f); Uc = split2d_neum(uf); D_f(idx) = max(abs(Uc(:,:,1:3)),[],'all');
        uz = toPhys_f(uh_z); Uz = split2d_neum(uz); D_z(idx) = max(abs(Uz(:,:,1:3)),[],'all');
        Vi = toPhys_i(uh_i);                        D_i(idx) = max(abs(Vi(:,:,1:3)),[],'all');
    end
    k = find(s_steps == n, 1);
    if ~isempty(k)
        Sf{k} = toPhys_f(uh_f);
        Si{k} = assem2d_neum(toPhys_i(uh_i));
    end
end
tb_f = trec(find(D_f > 1e-2, 1)); if isempty(tb_f), tb_f = NaN; end
tb_z = trec(find(D_z > 1e-2, 1)); if isempty(tb_z), tb_z = NaN; end
fprintf('2D Neumann X-junction (N=%d, eps=%.2f, T=%g): std T*=%.2f (maxD=%.2e) | zmode T*=%.2f (maxD=%.2e) | iso maxD=%g\n', ...
    N, eps, T, tb_f, max(D_f), tb_z, max(D_z), max(D_i));

% Plot figure
set(groot,'defaultAxesFontSize',9.5,'defaultLineLineWidth',1.2, ...
    'defaultAxesTickLabelInterpreter','latex','defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex','defaultAxesBox','on');
f = figure('Units','centimeters','Position',[2 2 18 7.8]);
left0 = 0.055; gapx = 0.015; w = 0.142; h_ax = 0.39; yrows = [0.53, 0.07];
for k = 1:4
    ax = axes('Position',[left0+(k-1)*(w+gapx), yrows(1), w, h_ax]);
    imagesc(x, x, Sf{k}'); axis xy equal tight; caxis([-1 1]); set(ax,'XTick',[],'YTick',[]);
    title(sprintf('$t = %g$', tsnaps(k)), 'FontSize', 9.5);
    if k == 1, ylabel({'standard','2D DCT-I'},'FontSize',9); end

    ax = axes('Position',[left0+(k-1)*(w+gapx), yrows(2), w, h_ax]);
    imagesc(x, x, Si{k}'); axis xy equal tight; caxis([-1 1]); set(ax,'XTick',[],'YTick',[]);
    if k == 1, ylabel({'isotypic','Neumann (ours)'},'FontSize',9); end
end
colormap(parula);
axR = axes('Position',[0.76, 0.14, 0.225, 0.78]);
cols = lines(3);
semilogy(axR, trec, max(D_f, 1e-19), '-',  'Color', cols(1,:), 'LineWidth', 1.4); hold on;
semilogy(axR, trec, max(D_z, 1e-19), '--', 'Color', cols(2,:), 'LineWidth', 1.4);
semilogy(axR, trec, max(D_i, 1e-19), '-',  'Color', cols(3,:), 'LineWidth', 1.8);
set(axR,'YScale','log'); xlim([0 T]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('$D(u_n)$');
legend({'standard DCT-I','zero-mode','isotypic ($D\equiv0$)'},'Location','southeast','FontSize',6.5);
exportgraphics(f, outpdf, 'ContentType','vector');
exportgraphics(f, outfig, 'Resolution', 180);
fprintf('Saved %s\n', outfig);

% =========================================================================
function W = xform2d_neum(V, dir)
W = zeros(size(V),'like',V);
for c = 0:3
    A = V(:,:,c+1); if ~any(A(:)), continue; end
    A = r1d_neum(A, bitget(c,1), dir);
    A = r1d_neum(A.', bitget(c,2), dir).';
    W(:,:,c+1) = A;
end
end

function Y = r1d_neum(X, parity, dir)
Q = size(X,1) - 1; m = size(X,2); Z = zeros(1, m, 'like', X);
if parity == 0
    ext = [X; X(Q:-1:2,:)];
    if dir > 0, F = real(fft(ext,[],1)); else, F = real(ifft(ext,[],1)); end
    Y = F(1:Q+1,:);
else
    if dir > 0
        ext = [X(1:Q,:); Z; -X(Q:-1:2,:); -X(1:Q,:); Z; X(Q:-1:2,:)];
        F = real(fft(ext,[],1)); Y = [F(2:2:2*Q,:); Z];
    else
        F = zeros(4*Q, m, 'like', X);
        F(2:2:2*Q,:) = X(1:Q,:); F(4*Q:-2:2*Q+2,:) = X(1:Q,:);
        v = real(ifft(F,[],1)); Y = [v(1:Q,:); Z];
    end
end
end

function U3 = cube4(U)
sz = size(U); A = reshape(U, [], 4);
A2 = zeros(size(A)); A3 = zeros(size(A));
for r = 0:3
    for a = 0:3, A2(:,r+1) = A2(:,r+1) + A(:,a+1).*A(:,bitxor(a,r)+1); end
end
for s = 0:3
    for r = 0:3, A3(:,s+1) = A3(:,s+1) + A2(:,r+1).*A(:,bitxor(s,r)+1); end
end
U3 = reshape(A3, sz);
end

function U = split2d_neum(u)
N = size(u,1)-1; Q = N/2;
vx0 = 0.5*(u(1:Q+1,:) + u(N+1:-1:Q+1,:));
vx1 = 0.5*(u(1:Q+1,:) - u(N+1:-1:Q+1,:)); vx1(end,:) = 0;
U = zeros(Q+1, Q+1, 4);
U(:,:,1) = 0.5*(vx0(:,1:Q+1) + vx0(:,N+1:-1:Q+1));
U(:,:,2) = 0.5*(vx1(:,1:Q+1) + vx1(:,N+1:-1:Q+1));
U(:,:,3) = 0.5*(vx0(:,1:Q+1) - vx0(:,N+1:-1:Q+1)); U(:,end,3) = 0;
U(:,:,4) = 0.5*(vx1(:,1:Q+1) - vx1(:,N+1:-1:Q+1)); U(:,end,4) = 0; U(end,:,4) = 0;
end

function u = assem2d_neum(U)
Q = size(U,1)-1; N = 2*Q; u = zeros(N+1, N+1);
for c = 0:3
    A = U(:,:,c+1); if ~any(A(:)), continue; end
    if bitget(c,1)==0, A = [A; A(Q:-1:1,:)]; else, A = [A; -A(Q:-1:1,:)]; end
    if bitget(c,2)==0, A = [A, A(:,Q:-1:1)]; else, A = [A, -A(:,Q:-1:1)]; end
    u = u + A;
end
end

function [E, Q1, Q2] = etd2_ops(Lk, tau)
z = -tau*Lk; E = exp(z);
Mq = 64; r = exp(1i*pi*((1:Mq)-.5)/Mq);
sz = size(z); zz = z(:); p1=zeros(size(zz)); p2=p1;
big = abs(zz)>1; zb = zz(big);
p1(big)=(exp(zb)-1)./zb; p2(big)=(exp(zb)-1-zb)./zb.^2;
zs = zz(~big);
if ~isempty(zs), Zc=zs+r; p1(~big)=real(mean((exp(Zc)-1)./Zc,2)); p2(~big)=real(mean((exp(Zc)-1-Zc)./Zc.^2,2)); end
Q1 = reshape(tau*p1,sz); Q2 = reshape(tau*p2,sz);
end
