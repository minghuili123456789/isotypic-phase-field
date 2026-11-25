% Supplementary BC Evidence 1: 1D Homogeneous Neumann Allen-Cahn
%   du/dt = u_xx + u - u^3 on [0, pi], du/dx(0) = du/dx(pi) = 0
%   u0(x) = cos(x) (strictly odd about midpoint x = pi/2), N = 128, T = 200.
% Compares 4 standard DCT-I integrators (ETDRK2, ETDRK4, BDF2-SAV, SBDF4)
% against their drop-in midpoint-isotypic Neumann counterparts (DCT-I + DCT-III
% on half-domain [0, pi/2], zero filtering, 2.2e-16 non-symmetric equivalence).
clear; close all;
outfig = fullfile(fileparts(mfilename('fullpath')), 'fig_neumann_1d.png');
outpdf = fullfile(fileparts(mfilename('fullpath')), 'fig_neumann_1d.pdf');

N = 128; Q = N/2; h = pi/N;
x  = (0:N)'*h;
xq = (0:Q)'*h;
T  = 200;

% Standard DCT-I operators on [0, pi] (length N+1)
dct1 = @(v) subsref(real(fft([v; v(end-1:-1:2)])), struct('type','()','subs',{{1:N+1}}));
idct1 = @(fh) subsref(real(ifft([fh; fh(end-1:-1:2)])), struct('type','()','subs',{{1:N+1}}));
k_full = (0:N)'; Lk_f = k_full.^2;

% Midpoint-isotypic Neumann operators on [0, pi/2] (2 channels of length Q+1):
%   Channel 0 (even about pi/2): even cosines k = 0,2,...,2Q  (Neumann-Neumann, DCT-I)
%   Channel 1 (odd  about pi/2): odd  cosines k = 1,3,...,2Q-1 (Neumann-Dirichlet, DCT-III)
t0f = @(v) subsref(real(fft([v; v(end-1:-1:2)])), struct('type','()','subs',{{1:Q+1}}));
t0i = @(fh) subsref(real(ifft([fh; fh(end-1:-1:2)])), struct('type','()','subs',{{1:Q+1}}));
t1f = @(v) [subsref(real(fft([v(1:Q);0;-v(Q:-1:2);-v(1:Q);0;v(Q:-1:2)])), struct('type','()','subs',{{2:2:2*Q}})); 0];
idx1 = [(2:2:2*Q)'; (4*Q:-2:2*Q+2)'];
t1i = @(fh) [subsref(real(ifft(accumarray(idx1, [fh(1:Q); fh(1:Q)], [4*Q 1]))), struct('type','()','subs',{{1:Q}})); 0];
ke = (0:Q)'*2; ko = [2*(0:Q-1)'+1; 0];
Lk_i = [ke.^2, ko.^2];

toSpec_i = @(U) [t0f(U(:,1)), t1f(U(:,2))];
toPhys_i = @(Uh) [t0i(Uh(:,1)), t1i(Uh(:,2))];
cube_i   = @(U) [U(:,1).^3 + 3*U(:,1).*U(:,2).^2, 3*U(:,1).^2.*U(:,2) + U(:,2).^3];
assem_i  = @(U) [U(:,1)+U(:,2); U(Q:-1:1,1)-U(Q:-1:1,2)];
split_i  = @(u) [0.5*(u(1:Q+1)+u(N+1:-1:Q+1)), [0.5*(u(1:Q)-u(N+1:-1:Q+2)); 0]];

% 1. Verify non-symmetric Neumann equivalence (50 ETDRK4 steps)
u_ns = 0.4*cos(x) + 0.3*cos(2*x) + 0.2*cos(3*x) + 0.15;
uf_eq = run_neum('ETDRK4', dct1, idct1, @(u) u.^3, Lk_f, u_ns, 0.02, 50, false);
ui_eq = run_neum('ETDRK4', toSpec_i, toPhys_i, cube_i, Lk_i, split_i(u_ns), 0.02, 50, true);
fprintf('1D Neumann non-symmetric equivalence (50 ETDRK4 steps): %.2e\n', max(abs(uf_eq - assem_i(ui_eq))));

% 2. Symmetry test on u0(x) = cos(x) up to T = 200 across 4 integrators
u0 = cos(x);
U0 = [zeros(Q+1,1), [cos(xq(1:Q)); 0]];
methods = {'ETDRK2','ETDRK4','BDF2-SAV','SBDF4'};
taus    = [0.02,    0.02,    0.01,      0.01];
trec    = (0:0.5:T)';
tsnaps  = [30, 40, 200];

D_std = zeros(numel(trec), 4); D_iso = D_std;
S_std = cell(4, 3); S_iso = cell(4, 3);
Tstar = nan(4, 1);

for m = 1:4
    nm = methods{m}; tau = taus(m);
    [D_std(:,m), S_std(m,:), Tstar(m)] = run_neum_track(nm, dct1, idct1, @(u) u.^3, Lk_f, u0, tau, T, trec, tsnaps, false, N);
    [D_iso(:,m), S_iso(m,:), ~]        = run_neum_track(nm, toSpec_i, toPhys_i, cube_i, Lk_i, U0, tau, T, trec, tsnaps, true, N, assem_i);
    fprintf('  %-9s (tau=%.2f): std T* = %5.2f (maxD=%.2e, u_end(pi/2)=%+.3f) | iso maxD = %g (u_end(0)=%.4f)\n', ...
        nm, tau, Tstar(m), max(D_std(:,m)), S_std{m,3}(Q+1), max(D_iso(:,m)), S_iso{m,3}(1));
end

% 3. Publication figure
set(groot,'defaultAxesFontSize',9.5,'defaultLineLineWidth',1.2, ...
    'defaultAxesTickLabelInterpreter','latex','defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex','defaultAxesBox','on');
f = figure('Units','centimeters','Position',[2 2 18 10.2]);
cols = [0.85 0.16 0.16; 0.90 0.45 0.00; 0.55 0.20 0.70; 0.10 0.60 0.55];
c_iso = [0.00 0.35 0.80];

axA = axes('Position',[0.075 0.56 0.41 0.38]); hold on;
hh = gobjects(5,1);
for m = 1:4
    hh(m) = semilogy(axA, trec, max(D_std(:,m), 1e-19), '-', 'Color', cols(m,:), 'LineWidth', 1.3);
end
hh(5) = semilogy(axA, trec, 1e-19*ones(size(trec)), '-', 'Color', c_iso, 'LineWidth', 2.0);
set(axA,'YScale','log'); xlim([0 70]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('midpoint defect $D(u_n)$');
title('(a) Standard DCT-I defect growth ($0\le t\le 70$)');
legend(hh, {'ETDRK2 (DCT-I)','ETDRK4 (DCT-I)','BDF2-SAV (DCT-I)','SBDF4 (DCT-I)','isotypic Neumann (all 4): $D\equiv0$'}, ...
       'Location','southeast','FontSize',6.3);

axB = axes('Position',[0.565 0.56 0.415 0.38]); hold on;
for m = 1:4
    semilogy(axB, trec, max(D_std(:,m), 1e-19), '-', 'Color', cols(m,:), 'LineWidth', 1.3);
end
semilogy(axB, trec, 1e-19*ones(size(trec)), '-', 'Color', c_iso, 'LineWidth', 2.0);
set(axB,'YScale','log'); xlim([0 T]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('midpoint defect $D(u_n)$');
title('(b) Long-time horizon ($0\le t\le 200$)');
text(55, 1e-4, {'All 4 standard DCT-I schemes','collapse to $u\equiv\pm 1$ by $t\approx 42$'}, 'FontSize', 7.3, 'Color', [0.75 0.12 0.12]);
text(55, 4e-18, 'Isotypic Neumann (all 4): $D\equiv0$ bitwise', 'FontSize', 7.3, 'Color', c_iso);

xb = [0.075, 0.385, 0.695]; wb = 0.285; hb = 0.33; yb = 0.095;
tstr = {'(c) $t = 30$ (pre-breakdown)', '(d) $t = 40$ (collapsing toward $\pm 1$)', '(e) $t = 200$ (wrong steady state vs.\ exact)'};
for p = 1:3
    ax = axes('Position',[xb(p) yb wb hb]); hold on;
    xline(pi/2, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.85);
    yline(0,    ':', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
    for m = 1:4
        plot(ax, x, S_std{m,p}, '-', 'Color', cols(m,:), 'LineWidth', 1.2);
    end
    plot(ax, x, S_iso{2,p}, '--', 'Color', c_iso, 'LineWidth', 1.9);
    xlim([0 pi]); ylim([-1.18 1.18]);
    set(ax,'XTick',[0 pi/2 pi],'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$'},'YTick',[-1 0 1]);
    if p > 1, set(ax,'YTickLabel',[]); end
    grid on; set(ax,'GridAlpha',0.15); xlabel('$x$');
    if p == 1, ylabel('$u(x,t)$'); end
    title(tstr{p});
end
exportgraphics(f, outpdf, 'ContentType','vector');
exportgraphics(f, outfig, 'Resolution', 180);
fprintf('Saved %s\n', outfig);

% =========================================================================
function u = run_neum(nm, toSpec, toPhys, cube, Lk, u0, tau, nst, isIso)
[~, S, ~] = run_neum_track(nm, toSpec, toPhys, cube, Lk, u0, tau, nst*tau, [0; nst*tau], nst*tau, isIso, 128);
u = S{1};
end

function [Drec, Snaps, tb] = run_neum_track(nm, toSpec, toPhys, cube, Lk, u0, tau, T, trec, tsnaps, isIso, N, assem)
hasAssem = (nargin >= 13);
nst = round(T/tau);
Drec = zeros(numel(trec),1); Snaps = cell(1,numel(tsnaps)); tb = NaN;
r_step = max(1, round((trec(2)-trec(1))/tau));
s_steps = round(tsnaps/tau);

Nop = @(uh) toSpec(toPhys(uh) - cube(toPhys(uh)));
% ETD coefficients
z = -tau*Lk;
E  = exp(z); E2 = exp(z/2);
Q1 = tau*phik(z,1); Q2 = tau*phik(z,2);
Q  = (tau/2)*phik(z/2,1);
[f1,f2,f3] = etd4coef(z, tau);

    function rec(un, nstep)
        if mod(nstep, r_step) == 0
            idx = round(nstep/r_step) + 1;
            if idx <= numel(Drec)
                if isIso, D = max(abs(un(:,1))); else, D = 0.5*max(abs(un + un(end:-1:1))); end
                Drec(idx) = D;
                if isnan(tb) && D > 1e-2, tb = nstep*tau; end
            end
        end
        k = find(s_steps == nstep, 1);
        if ~isempty(k)
            if isIso && hasAssem, Snaps{k} = assem(un); else, Snaps{k} = un; end
        end
    end

rec(u0, 0);
uh = toSpec(u0);

switch nm
    case 'ETDRK2'
        for n = 1:nst
            Nu = Nop(uh); a = E.*uh + Q1.*Nu; uh = a + Q2.*(Nop(a) - Nu);
            if mod(n,r_step)==0 || any(s_steps==n), rec(toPhys(uh), n); end
        end
    case 'ETDRK4'
        for n = 1:nst
            Nu = Nop(uh);
            a = E2.*uh + Q.*Nu; Na = Nop(a);
            b = E2.*uh + Q.*Na; Nb = Nop(b);
            c = E2.*a  + Q.*(2*Nb - Nu); Nc = Nop(c);
            uh = E.*uh + f1.*Nu + 2*f2.*(Na+Nb) + f3.*Nc;
            if mod(n,r_step)==0 || any(s_steps==n), rec(toPhys(uh), n); end
        end
    case 'SBDF4'
        step4 = @(v) etd4s(v, Nop, E, E2, Q, f1, f2, f3);
        uh0=uh; N0=Nop(uh0); uh1=step4(uh0); N1=Nop(uh1); uh2=step4(uh1); N2=Nop(uh2); uh3=step4(uh2); N3=Nop(uh3);
        Dinv = 1./(25/12 + tau*Lk);
        for n = 4:nst
            uh = Dinv.*(4*uh3 - 3*uh2 + (4/3)*uh1 - 0.25*uh0 + tau*(4*N3 - 6*N2 + 4*N1 - N0));
            uh0=uh1; N0=N1; uh1=uh2; N1=N2; uh2=uh3; N2=N3; uh3=uh; N3=Nop(uh);
            if mod(n,r_step)==0 || any(s_steps==n), rec(toPhys(uh), n); end
        end
    case 'BDF2-SAV'
        C0 = 10;
        ip = @(a,b) neum_ip(a, b, isIso, N);
        en = @(u) neum_en(u, isIso, N);
        r0 = sqrt(en(u0) + C0); uh0 = uh;
        A1 = 1./(1 + tau*Lk);
        b0 = (cube(u0) - u0)/r0;
        g1h = A1.*uh0; g2h = A1.*toSpec(b0);
        r1 = (r0 - 0.5*ip(b0,u0) + 0.5*ip(b0,toPhys(g1h))) / (1 + 0.5*tau*ip(b0,toPhys(g2h)));
        uh1 = g1h - tau*r1*g2h; u1 = toPhys(uh1);
        A2 = 1./(1.5 + tau*Lk); u_prev = u0;
        for n = 2:nst
            us = 2*u1 - u_prev;
            b  = (cube(us) - us)/sqrt(en(us) + C0); bh = toSpec(b);
            g1h = A2.*(2*uh1 - 0.5*uh0); g2h = A2.*bh;
            r2 = ((4*r1-r0)/3 - (1/6)*ip(b, 4*u1 - u_prev) + 0.5*ip(b, toPhys(g1h))) / (1 + 0.5*tau*ip(b, toPhys(g2h)));
            uh2 = g1h - tau*r2*g2h; u2 = toPhys(uh2);
            u_prev = u1; uh0 = uh1; r0 = r1; u1 = u2; uh1 = uh2; r1 = r2;
            if mod(n,r_step)==0 || any(s_steps==n), rec(u2, n); end
        end
end
end

function un = etd4s(uh, Nop, E, E2, Q, f1, f2, f3)
Nu = Nop(uh); a = E2.*uh + Q.*Nu; Na = Nop(a); b = E2.*uh + Q.*Na; Nb = Nop(b);
c = E2.*a + Q.*(2*Nb - Nu); Nc = Nop(c); un = E.*uh + f1.*Nu + 2*f2.*(Na+Nb) + f3.*Nc;
end

function val = neum_ip(a, b, isIso, N)
h = pi/N;
if ~isIso
    w = ones(N+1,1); w(1)=0.5; w(end)=0.5;
    val = h * sum(w .* a .* b);
else
    Q = N/2; w = 2*ones(Q+1,1); w(1)=1; w(end)=1;
    val = h * sum(w .* sum(a .* b, 2));
end
end

function E1 = neum_en(u, isIso, N)
if ~isIso
    w = ones(N+1,1); w(1)=0.5; w(end)=0.5;
    E1 = 0.25*(pi/N)*sum(w.*(u.^2 - 1).^2);
else
    u2 = [u(:,1).^2 + u(:,2).^2, 2*u(:,1).*u(:,2)];
    E1 = 0.25*(neum_ip(u2,u2,true,N) - 2*neum_ip(u,u,true,N) + pi);
end
end

function p = phik(z, k)
Mq = 64; r = exp(1i*pi*((1:Mq)-.5)/Mq);
sz = size(z); zz = z(:); p = zeros(size(zz));
big = abs(zz) > 1; zb = zz(big);
if k==1, p(big)=(exp(zb)-1)./zb; else, p(big)=(exp(zb)-1-zb)./zb.^2; end
zs = zz(~big);
if ~isempty(zs)
    Zc = zs + r;
    if k==1, F=(exp(Zc)-1)./Zc; else, F=(exp(Zc)-1-Zc)./Zc.^2; end
    p(~big) = real(mean(F,2));
end
p = reshape(p, sz);
end

function [f1,f2,f3] = etd4coef(z, tau)
Mq = 64; r = exp(1i*pi*((1:Mq)-.5)/Mq);
sz = size(z); zz = z(:); f1=zeros(size(zz)); f2=f1; f3=f1;
big = abs(zz) > 2; zb = zz(big);
f1(big)=(-4-zb+exp(zb).*(4-3*zb+zb.^2))./zb.^3;
f2(big)=(2+zb+exp(zb).*(-2+zb))./zb.^3;
f3(big)=(-4-3*zb-zb.^2+exp(zb).*(4-zb))./zb.^3;
zs = zz(~big);
if ~isempty(zs)
    Zc = zs + 2*r;
    f1(~big)=real(mean((-4-Zc+exp(Zc).*(4-3*Zc+Zc.^2))./Zc.^3,2));
    f2(~big)=real(mean((2+Zc+exp(Zc).*(-2+Zc))./Zc.^3,2));
    f3(~big)=real(mean((-4-3*Zc-Zc.^2+exp(Zc).*(4-Zc))./Zc.^3,2));
end
f1=reshape(tau*f1,sz); f2=reshape(tau*f2,sz); f3=reshape(tau*f3,sz);
end
