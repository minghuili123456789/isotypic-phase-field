% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Comprehensive test of Direction B:
%  (1) Half-grid DCT-II/DST-II (fixed sign) in 1D/2D/3D
%  (2) Half-period anti-symmetric T2 in 1D/2D/3D with IRRATIONAL shifts
%  (3) Canonical-Frame Isotypic Spectral Method (CF-ISM) with zero-aliasing
%      Gauss-Newton gauge fixing on mixed-harmonic shifted data (Wang Ex 5.2 +
%      irrational shifts + 2D shifted X-junction + 2D shifted C2 + faithfulness
%      + non-symmetric equivalence!)
addpath(fileparts(mfilename('fullpath')));

fprintf('=== (1) Half-grid DCT-II / DST-II equivalence + invariance ===\n');
for d = 1:3
    N = [64 32 16]; N = N(d); M = N/2; h = 2*pi/N;
    x1 = ((0:N-1)' + 0.5)*h;
    switch d
        case 1, u_full = 0.5*sin(x1) + 0.3*cos(2*x1) + 0.2*sin(3*x1 + 0.4);
        case 2, [X,Y] = ndgrid(x1,x1); u_full = sin(X).*cos(Y) + 0.3*cos(2*X+Y) + 0.2*sin(X-3*Y);
        case 3, [X,Y,Z] = ndgrid(x1,x1,x1); u_full = sin(X).*cos(Y).*sin(Z) + 0.3*cos(2*X+Y-Z);
    end
    U = hg_split(u_full, d, M);
    Uh = hg_xform(U, d, M, +1);
    K2 = hg_k2(d, M);
    U2 = hg_xform(exp(-0.1*K2).*Uh, d, M, -1);
    u_hg = hg_assemble(U2, d, M);
    kf = [0:M, -M+1:-1]';
    switch d
        case 1, K2f = kf.^2; u_ref = real(ifft(exp(-0.1*K2f).*fft(u_full)));
        case 2, [a,b]=ndgrid(kf,kf); K2f=a.^2+b.^2; u_ref = real(ifftn(exp(-0.1*K2f).*fftn(u_full)));
        case 3, [a,b,c]=ndgrid(kf,kf,kf); K2f=a.^2+b.^2+c.^2; u_ref = real(ifftn(exp(-0.1*K2f).*fftn(u_full)));
    end
    fprintf('  d=%d N=%2d: |half-grid iso - Fourier| = %.2e\n', d, N, max(abs(u_hg(:)-u_ref(:))));
end

fprintf('\n=== (1b) T2 wavenumber-parity equivalence in 1D/2D/3D (non-symmetric data) ===\n');
for d = 1:3
    Nt = [128 64 32]; Nt = Nt(d); h = 2*pi/Nt; x1 = (0:Nt-1)'*h;
    switch d
        case 1, fun = @(x) 0.5*sin(x) + 0.3*cos(2*x) + 0.2*sin(3*x + 0.4); u_f = fun(x1);
        case 2, fun = @(x,y) sin(x).*cos(y) + 0.3*cos(2*x+y) + 0.2*sin(x-3*y); [X,Y]=ndgrid(x1,x1); u_f = fun(X,Y);
        case 3, fun = @(x,y,z) sin(x).*cos(y).*sin(z) + 0.3*cos(2*x+y-z); [X,Y,Z]=ndgrid(x1,x1,x1); u_f = fun(X,Y,Z);
    end
    Pt = sp_t2_setup(struct('d',d,'N',Nt,'eps2',1,'gam',1,'S',0)); Ot = sp_t2_ops(Pt);
    Ut = sp_t2_init(Pt, fun);
    u_t2 = sp_t2_assemble(Ot.toPhys(exp(-0.1*Pt.K2).*Ot.toSpec(Ut)), d, Nt/2);
    Pf = sp_setup(struct('d',d,'N',Nt,'eps2',1,'gam',1,'S',0,'rep','fourier')); Of = sp_ops(Pf);
    u_ref = Of.toPhys(exp(-0.1*Pf.K2).*Of.toSpec(u_f));
    fprintf('  d=%d N=%3d: |T2 iso - Fourier| = %.2e\n', d, Nt, max(abs(u_t2(:)-u_ref(:))));
end

fprintf('\n=== (2) CF-ISM with Gauss-Newton gauge on mixed-harmonic shifted data ===\n');
% Test all shifts: 0, pi/8, pi/4, pi/128 (half-grid), and 0.37123456789 (irrational)
N = 128; M = N/2; h = 2*pi/N;
shifts = [0, pi/8, pi/4, pi/N, 0.37123456789];
for q = 1:numel(shifts)
    sh = shifts(q);
    f_dbl = @(x) 0.5*(sin(x + sh) + sin(4*(x + sh)));
    f_vpa = @(x) 0.5*(sin(x + vpa(sh)) + sin(4*(x + vpa(sh))));
    [c_q, U0] = cf_init_1d(f_dbl, f_vpa, N, 34);
    % Run pure double-precision isotypic ETDRK4 (S=0, kappa=1, T=100)
    P = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));
    O = sp_ops(P); I = sp_integrator(P,O,'etdrk4',0.05);
    vh = O.toSpec(U0);
    maxD = max(abs(U0(:,1)));
    for n = 1:2000
        vh = I.step(vh);
        V = O.toPhys(vh);
        maxD = max(maxD, max(abs(V(:,1))));
    end
    fprintf('  sh=%11.8f: c_err=%.2e | init max|ev|=%g | T=100 maxD=%g | range=[%.4f,%.4f]\n', ...
        sh, double(min(abs(mod(c_q + vpa(sh) + vpa(pi)/2, vpa(pi)) - vpa(pi)/2))), ...
        max(abs(U0(:,1))), maxD, min(V(:,2)), max(V(:,2)));
end

fprintf('\n=== (3) CF-ISM faithfulness with irrational shift + genuine asymmetry delta ===\n');
sh = 0.37123456789;
for del = [1e-10, 1e-30]
    f_dbl = @(x) sin(x + sh) + del;
    f_vpa = @(x) sin(x + vpa(sh)) + vpa(del);
    [c_q, U0] = cf_init_1d(f_dbl, f_vpa, N, 60);
    P = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));
    O = sp_ops(P); I = sp_integrator(P,O,'etdrk4',0.05);
    vh = O.toSpec(U0); tb = NaN;
    for n = 1:2000
        vh = I.step(vh); V = O.toPhys(vh);
        if isnan(tb) && max(abs(V(:,1))) > 1e-2, tb = n*0.05; break; end
    end
    fprintf('  sh=%.5f delta=%.0e: init ev=%.3e (err=%.2e) | T*(D>1e-2)=%.2f\n', ...
        sh, del, max(abs(U0(:,1))), abs(max(abs(U0(:,1)))-del)/del, tb);
end

fprintf('\n=== (4) CF-ISM on NON-SYMMETRIC data (must match standard Fourier!) ===\n');
f_ns_d = @(x) 0.5*sin(x) + 0.3*cos(2*x) + 0.2*sin(3*x + 0.4);
f_ns_v = @(x) 0.5*sin(x) + 0.3*cos(2*x) + 0.2*sin(3*x + vpa(0.4));
[c_q, U0] = cf_init_1d(f_ns_d, f_ns_v, N, 34);
P_iso = sp_setup(struct('d',1,'N',N,'eps2',0.1,'gam',1,'S',1,'rep','iso'));
O_iso = sp_ops(P_iso); I_iso = sp_integrator(P_iso,O_iso,'etdrk4',0.05);
vh = O_iso.toSpec(U0);
for n = 1:40, vh = I_iso.step(vh); end
v_full = sp_assemble(O_iso.toPhys(vh), 1, M);
% Shift back by -c_q in Fourier space to compare on original grid x_j = j*h
c_d = double(c_q); kf = [0:M, -M+1:-1]';
u_cf_orig = real(ifft(fft(v_full) .* exp(-1i*kf*c_d)));
% Standard Fourier on original grid
P_f = sp_setup(struct('d',1,'N',N,'eps2',0.1,'gam',1,'S',1,'rep','fourier'));
O_f = sp_ops(P_f); I_f = sp_integrator(P_f,O_f,'etdrk4',0.05);
uh = O_f.toSpec(f_ns_d((0:N-1)'*h));
for n = 1:40, uh = I_f.step(uh); end
u_std = O_f.toPhys(uh);
fprintf('  non-symmetric data (40 steps): |CF-ISM - Fourier| = %.2e\n', max(abs(u_cf_orig - u_std)));

fprintf('\n=== (5) 2D shifted mixed-harmonic & 2D shifted C2 (irrational shifts!) ===\n');
sh1 = 0.37123456; sh2 = 0.71828182;
% 2D odd-odd mixed harmonic: (sin(x+sh1)+0.5*sin(4*(x+sh1))) * (sin(y+sh2)+0.5*sin(4*(y+sh2)))
N2 = 64; M2 = N2/2;
f2_d = @(x,y) 0.5*(sin(x+sh1)+0.5*sin(4*(x+sh1))).*(sin(y+sh2)+0.5*sin(4*(y+sh2)));
f2_v = @(x,y) 0.5*(sin(x+vpa(sh1))+0.5*sin(4*(x+vpa(sh1)))).*(sin(y+vpa(sh2))+0.5*sin(4*(y+vpa(sh2))));
[c2, U2_0] = cf_init_2d(f2_d, f2_v, N2, 34);
fprintf('  2D shifted odd-odd: c_err=(%.1e, %.1e), init forbidden max=%g\n', ...
    double(min(abs(mod(c2(1)+vpa(sh1)+vpa(pi)/2,vpa(pi))-vpa(pi)/2))), ...
    double(min(abs(mod(c2(2)+vpa(sh2)+vpa(pi)/2,vpa(pi))-vpa(pi)/2))), ...
    max(abs(reshape(U2_0(:,:,1:3),[],1))));
% 2D shifted C2 point-odd (companion benchmark shifted by off-grid sh1, sh2!)
% Note: sin((x+sh1) + 2.5*sin(3*(y+sh2))) has C2 centres at x+sh1 = m*pi, y+sh2 = n*pi/3!
ep = 0.05;
fc2_d = @(x,y) tanh(sin((x+sh1) + 2.5*sin(3*(y+sh2)))/(sqrt(2)*ep));
fc2_v = @(x,y) tanh(sin((x+vpa(sh1)) + vpa(2.5)*sin(3*(y+vpa(sh2))))/(sqrt(vpa(2))*vpa(ep)));
[cc2, Uc2_0] = cf_init_2d_pt(fc2_d, fc2_v, N2, 320);
Pc2 = sp_setup(struct('d',2,'N',N2,'eps2',1,'gam',1/ep^2,'S',2/ep^2,'rep','iso'));
Oc2 = sp_ops(Pc2); Ic2 = sp_integrator(Pc2, Oc2, 'etdrk2', 2e-3);
vhc = Oc2.toSpec(Uc2_0); maxDc2 = max([max(abs(Uc2_0(:,:,1)),[],'all'), max(abs(Uc2_0(:,:,4)),[],'all')]);
for n = 1:200
    vhc = Ic2.step(vhc); Vc = Oc2.toPhys(vhc);
    maxDc2 = max(maxDc2, max([max(abs(Vc(:,:,1)),[],'all'), max(abs(Vc(:,:,4)),[],'all')]));
end
fprintf('  2D shifted C2 (d0=320, 200 ETDRK2 steps): c_err=(%.1e, %.1e mod pi/3), init=(ee,oo)=%g, maxD=%g\n', ...
    double(min(abs(mod(cc2(1)+vpa(sh1)+vpa(pi)/2,vpa(pi))-vpa(pi)/2))), ...
    double(min(abs(mod(cc2(2)+vpa(sh2)+vpa(pi)/6,vpa(pi)/3)-vpa(pi)/6))), ...
    max([max(abs(Uc2_0(:,:,1)),[],'all'), max(abs(Uc2_0(:,:,4)),[],'all')]), maxDc2);

% -------- Helper functions --------
function [c_q, U0] = cf_init_1d(f_dbl, f_vpa, N, dig)
% Intrinsic translational gauge fixing (no knowledge of symmetry, no zeroing!)
M = N/2; h = 2*pi/N; xf = (0:N-1)'*h;
u = f_dbl(xf); uh = fft(u);
[~, idx] = min(real(ifft(uh.^2))/N);
c = (idx - 1)*h/2;
% Coarse double Newton on physical oddness J(c) = 0.5 sum (f(c+x)+f(c-x))^2
Mc = min(32, M); xc = (1:Mc-1)'*pi/Mc;
eps_fd = 1e-7;
for it = 1:6
    r = f_dbl(c + xc) + f_dbl(c - xc);
    dr = (f_dbl(c + eps_fd + xc) + f_dbl(c + eps_fd - xc) - (f_dbl(c - eps_fd + xc) + f_dbl(c - eps_fd - xc)))/(2*eps_fd);
    denom = dr'*dr; if denom < 1e-30, break; end
    dc = -(dr'*r)/denom; c = c + dc;
    if abs(dc) < 1e-14, break; end
end
% Map c to [-pi/2, pi/2) (canonical representative mod pi)
c = mod(c + pi/2, pi) - pi/2;
% Extended-precision Gauss-Newton polish on Mc points (at t=0 ONLY!)
digits(dig);
c_q = vpa(c); xc_q = (1:Mc-1)'*vpa(pi)/Mc;
eps_q = vpa(10)^(-floor(dig/2));
for it = 1:5
    r = f_vpa(c_q + xc_q) + f_vpa(c_q - xc_q);
    dr = (f_vpa(c_q + eps_q + xc_q) + f_vpa(c_q + eps_q - xc_q) - (f_vpa(c_q - eps_q + xc_q) + f_vpa(c_q - eps_q - xc_q)))/(2*eps_q);
    denom = sum(dr.^2); if double(denom) < 1e-30, break; end
    dc = -sum(dr.*r)/denom; c_q = c_q + dc;
    if double(abs(dc)) < 10^(-(dig-4)), break; end
end
% Evaluate channels on fundamental grid [0, pi] in extended precision -> round to double
xh_q = (0:M)' * (2*vpa(pi)/N);
vp = f_vpa(c_q + xh_q); vm = f_vpa(c_q - xh_q);
ev = double((vp + vm)/2); od = double((vp - vm)/2);
od([1 end]) = 0;   % structural boundary non-DOFs of DST-I
U0 = [ev, od];
end

function [c_q, U0] = cf_init_2d(f_dbl, f_vpa, N, dig)
% 2D reflection-odd gauge fixing along x and y separately
M = N/2; h = 2*pi/N; x1 = (0:N-1)'*h; [X,Y] = ndgrid(x1,x1);
U = f_dbl(X,Y);
% 1D auto-convolution summed over other axis
Ax = sum(real(ifft(fft(U,[],1).^2,[],1)), 2); [~,ix] = min(Ax); cx = (ix-1)*h/2;
Ay = sum(real(ifft(fft(U,[],2).^2,[],2)), 1); [~,iy] = min(Ay); cy = (iy-1)*h/2;
cx = mod(cx+pi/2,pi)-pi/2; cy = mod(cy+pi/2,pi)-pi/2;
digits(dig);
c_q = [vpa(cx); vpa(cy)];
Mc = 12; xc = (1:Mc-1)'*vpa(pi)/Mc; [Xc,Yc] = ndgrid(xc,xc);
Xc = Xc(:); Yc = Yc(:);
eps_q = vpa(10)^(-floor(dig/2));
for it = 1:6
    rx = f_vpa(c_q(1)+Xc, c_q(2)+Yc) + f_vpa(c_q(1)-Xc, c_q(2)+Yc);
    ry = f_vpa(c_q(1)+Xc, c_q(2)+Yc) + f_vpa(c_q(1)+Xc, c_q(2)-Yc);
    drx = (f_vpa(c_q(1)+eps_q+Xc, c_q(2)+Yc) + f_vpa(c_q(1)+eps_q-Xc, c_q(2)+Yc) ...
         - f_vpa(c_q(1)-eps_q+Xc, c_q(2)+Yc) - f_vpa(c_q(1)-eps_q-Xc, c_q(2)+Yc))/(2*eps_q);
    dry = (f_vpa(c_q(1)+Xc, c_q(2)+eps_q+Yc) + f_vpa(c_q(1)+Xc, c_q(2)+eps_q-Yc) ...
         - f_vpa(c_q(1)+Xc, c_q(2)-eps_q+Yc) - f_vpa(c_q(1)+Xc, c_q(2)-eps_q-Yc))/(2*eps_q);
    dc1 = -sum(drx.*rx)/sum(drx.^2); dc2 = -sum(dry.*ry)/sum(dry.^2);
    c_q = c_q + [dc1; dc2];
    if double(max(abs([dc1 dc2]))) < 10^(-(dig-4)), break; end
end
U0 = eval_iso_2d_vpa(f_vpa, c_q, N);
end

function [c_q, U0] = cf_init_2d_pt(f_dbl, f_vpa, N, dig)
% 2D point-reflection (C2) gauge fixing via 2D auto-convolution + Gauss-Newton
M = N/2; h = 2*pi/N; x1 = (0:N-1)'*h; [X,Y] = ndgrid(x1,x1);
U = f_dbl(X,Y);
A = real(ifftn(fftn(U).^2));
[~, idx] = min(A(:)); [ix, iy] = ind2sub([N N], idx);
cx = mod((ix-1)*h/2 + pi/2, pi) - pi/2;
cy = mod((iy-1)*h/2 + pi/2, pi) - pi/2;
digits(dig); c_q = [vpa(cx); vpa(cy)];
Mc = 10; xc = (1:Mc-1)'*vpa(pi)/Mc; [Xc,Yc] = ndgrid(xc, [xc; -xc]);
Xc = Xc(:); Yc = Yc(:);
eps_q = vpa(10)^(-floor(dig/2));
for it = 1:7
    r = f_vpa(c_q(1)+Xc, c_q(2)+Yc) + f_vpa(c_q(1)-Xc, c_q(2)-Yc);
    J1 = (f_vpa(c_q(1)+eps_q+Xc, c_q(2)+Yc) + f_vpa(c_q(1)+eps_q-Xc, c_q(2)-Yc) ...
        - f_vpa(c_q(1)-eps_q+Xc, c_q(2)+Yc) - f_vpa(c_q(1)-eps_q-Xc, c_q(2)-Yc))/(2*eps_q);
    J2 = (f_vpa(c_q(1)+Xc, c_q(2)+eps_q+Yc) + f_vpa(c_q(1)-Xc, c_q(2)+eps_q-Yc) ...
        - f_vpa(c_q(1)+Xc, c_q(2)-eps_q+Yc) - f_vpa(c_q(1)-Xc, c_q(2)-eps_q-Yc))/(2*eps_q);
    J = [J1, J2];
    dc = -(J.'*J) \ (J.'*r);
    c_q = c_q + dc;
    if double(max(abs(dc))) < 10^(-(dig-4)), break; end
end
U0 = eval_iso_2d_vpa(f_vpa, c_q, N);
end

function U0 = eval_iso_2d_vpa(f_vpa, c_q, N)
M = N/2; xh = (0:M)'*(2*vpa(pi)/N); [X,Y] = ndgrid(xh,xh);
vpp = f_vpa(c_q(1)+X, c_q(2)+Y);
vmp = f_vpa(c_q(1)-X, c_q(2)+Y);
vpm = f_vpa(c_q(1)+X, c_q(2)-Y);
vmm = f_vpa(c_q(1)-X, c_q(2)-Y);
% butterflies in VPA -> exact cancellation -> round to double
ep = (vpp + vmp)/2; op = (vpp - vmp)/2;
em = (vpm + vmm)/2; om = (vpm - vmm)/2;
U0 = zeros(M+1, M+1, 4);
U0(:,:,1) = double((ep + em)/2);
U0(:,:,2) = double((op + om)/2); U0([1 end],:,2) = 0;
U0(:,:,3) = double((ep - em)/2); U0(:,[1 end],3) = 0;
U0(:,:,4) = double((op - om)/2); U0([1 end],:,4) = 0; U0(:,[1 end],4) = 0;
end

function U = hg_split(u, d, M)
N = 2*M; nc = 2^d; vals = cell(1,nc); vals{1} = u; refl = N:-1:1;
for i = 1:d
    nv = cell(1,nc);
    for g = 0:2^(i-1)-1
        p = vals{g+1}; cc = repmat({':'},1,d); cc{i} = refl; q = p(cc{:});
        nv{g+1} = (p + q)/2; nv{g+1+2^(i-1)} = (p - q)/2;
    end
    vals(1:2^i) = nv(1:2^i);
end
cc = repmat({1:M},1,d); for c = 1:nc, vals{c} = vals{c}(cc{:}); end
U = cat(d+1, vals{:});
end

function uf = hg_assemble(U, d, M)
N = 2*M; nc = 2^d; map = [1:M, M:-1:1]; sgn = [ones(1,M), -ones(1,M)];
colons = repmat({':'},1,d); uf = [];
for c = 0:nc-1
    A = U(colons{:},c+1);
    for i = 1:d
        cc = repmat({':'},1,d); cc{i} = map; A = A(cc{:});
        if bitget(c,i), shp = ones(1,d); shp(i) = N; A = A .* reshape(sgn,[shp 1]); end
    end
    if isempty(uf), uf = A; else, uf = uf + A; end
end
end

function w = hg_xform(v, d, M, dir)
nc = 2^d; w = v; colons = repmat({':'},1,d);
for c = 0:nc-1
    vc = v(colons{:},c+1);
    for i = 1:d, vc = r2r_type2(vc, i, bitget(c,i), M, dir); end
    w(colons{:},c+1) = vc;
end
end

function v = r2r_type2(v, dim, isodd, M, dir)
nd = max(ndims(v), dim); shp = ones(1,nd); shp(dim) = 2*M;
kf = [0:M, -M+1:-1]';
ph_fwd = reshape(exp(-1i*kf*pi/(2*M)), shp);
ph_inv = reshape(exp(+1i*kf*pi/(2*M)), shp);
if ~isodd
    if dir > 0
        ext = cat(dim, v, flip(v,dim));
        V = fft(ext,[],dim) .* ph_fwd;
        v = real(sub_(V, dim, 1:M));
    else
        z = zeros(size(sub_(v,dim,1)),'like',v);
        spec = cat(dim, v, z, flip(sub_(v,dim,2:M),dim));
        V = ifft(spec .* ph_inv, [], dim);
        v = real(sub_(V, dim, 1:M));
    end
else
    if dir > 0
        ext = cat(dim, v, -flip(v,dim));
        V = fft(ext,[],dim) .* ph_fwd;
        v = -imag(sub_(V, dim, 2:M+1));
    else
        z = zeros(size(sub_(v,dim,1)),'like',v);
        % FIXED SIGN: sin(k x) = (-i/2) exp(i k x) + (i/2) exp(-i k x)
        spec = cat(dim, z, -1i*sub_(v,dim,1:M-1), -1i*sub_(v,dim,M), +1i*flip(sub_(v,dim,1:M-1),dim));
        V = ifft(spec .* ph_inv, [], dim);
        v = real(sub_(V, dim, 1:M));
    end
end
end

function K2 = hg_k2(d, M)
nc = 2^d; K2 = zeros([repmat(M,1,d) nc]);
for c = 0:nc-1
    ks = cell(1,d);
    for i = 1:d, if bitget(c,i), ks{i} = (1:M)'; else, ks{i} = (0:M-1)'; end, end
    if d == 1, K2(:,c+1) = ks{1}.^2;
    elseif d == 2, [a,b]=ndgrid(ks{:}); K2(:,:,c+1)=a.^2+b.^2;
    else, [a,b,e]=ndgrid(ks{:}); K2(:,:,:,c+1)=a.^2+b.^2+e.^2; end
end
end

function B = sub_(A, dim, ix)
c = repmat({':'},1,max(ndims(A),dim)); c{dim} = ix; B = A(c{:});
end
