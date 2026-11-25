% Supplementary BC Evidence 3: 1D Homogeneous Dirichlet Allen-Cahn
%   du/dt = kappa^2 * u_xx + u - u^3 on [0, pi], u(0) = u(pi) = 0
%   kappa = 0.5, u0(x) = sin(2x) (strictly odd about midpoint x = pi/2), N = 128, T = 200.
% Note: there is NO k=0 constant mode here! The unstable transverse mode is sin(x)
% (even about pi/2, eigenvalue 1 - kappa^2 = +0.75 > 0).
% Standard DST-I (ETDRK2, ETDRK4, BDF2-SAV, SBDF4) breaks midpoint oddness and collapses
% to a single-sign positive/negative arch, whereas our 2-channel midpoint-isotypic
% Dirichlet backend (DST-I for even k + DST-III for odd k on [0, pi/2]) keeps D === 0!
clear; close all;
outfig = fullfile(fileparts(mfilename('fullpath')), 'fig_dirichlet_1d.png');
outpdf = fullfile(fileparts(mfilename('fullpath')), 'fig_dirichlet_1d.pdf');

N = 128; Q = N/2; h = pi/N; kap = 0.5;
x  = (0:N)'*h; xq = (0:Q)'*h; T = 200;

% Standard DST-I on [0, pi] (length N+1, zeros at 1 and N+1)
dst1  = @(v) [0; subsref(-imag(fft([0; v(2:N); 0; -v(N:-1:2)])), struct('type','()','subs',{{2:N}})); 0];
idst1 = @(fh) [0; subsref(imag(ifft([0; fh(2:N); 0; -fh(N:-1:2)])), struct('type','()','subs',{{2:N}})); 0];
k_f = (0:N)'; Lk_f = kap^2 * k_f.^2;

% 2-channel Midpoint-Isotypic Dirichlet on [0, pi/2] (length Q+1):
%   Channel 0 (EVEN about pi/2): odd  sines sin((2m+1)x), 0 at x=0, symmetric at x=pi/2 (DST-III)
%   Channel 1 (ODD  about pi/2): even sines sin(2m x),    0 at x=0 and x=pi/2           (DST-I)
t1f = @(v) [0; subsref(-imag(fft([0; v(2:Q); 0; -v(Q:-1:2)])), struct('type','()','subs',{{2:Q}})); 0];
t1i = @(fh) [0; subsref(imag(ifft([0; fh(2:Q); 0; -fh(Q:-1:2)])), struct('type','()','subs',{{2:Q}})); 0];
t0f = @(v) [0; subsref(-imag(fft([0; v(2:Q+1); v(Q:-1:2); 0; -v(2:Q+1); -v(Q:-1:2)])), struct('type','()','subs',{{2:2:2*Q}}))];
idx0 = [(2:2:2*Q)'; (4*Q:-2:2*Q+2)'];
t0i = @(fh) [0; subsref(real(ifft(accumarray(idx0, [-1i*fh(2:Q+1); 1i*fh(2:Q+1)], [4*Q 1]))), struct('type','()','subs',{{2:Q+1}}))];

ko = [0; 2*(0:Q-1)'+1]; ke = (0:Q)'*2;
Lk_i = kap^2 * [ko.^2, ke.^2];
toSpec_i = @(U) [t0f(U(:,1)), t1f(U(:,2))];
toPhys_i = @(Uh) [t0i(Uh(:,1)), t1i(Uh(:,2))];
cube_i   = @(U) [U(:,1).^3 + 3*U(:,1).*U(:,2).^2, 3*U(:,1).^2.*U(:,2) + U(:,2).^3];
assem_i  = @(U) [U(:,1)+U(:,2); U(Q:-1:1,1)-U(Q:-1:1,2)];
split_i  = @(u) [0.5*(u(1:Q+1)+u(N+1:-1:Q+1)), [0.5*(u(1:Q)-u(N+1:-1:Q+2)); 0]];

% 1. Verify non-symmetric Dirichlet equivalence (50 ETDRK4 steps)
u_ns = 0.4*sin(x) + 0.3*sin(2*x) + 0.2*sin(3*x);
uf_eq = run_dir('ETDRK4', dst1, idst1, @(u) u.^3, Lk_f, u_ns, 0.02, 50, false, N);
ui_eq = run_dir('ETDRK4', toSpec_i, toPhys_i, cube_i, Lk_i, split_i(u_ns), 0.02, 50, true, N);
fprintf('1D Dirichlet non-symmetric equivalence (50 ETDRK4 steps): %.2e\n', max(abs(uf_eq - assem_i(ui_eq))));

% 2. Symmetry test on u0(x) = sin(2x), kappa = 0.5, T = 200
u0 = sin(2*x); u0([1 end]) = 0;
U0 = [zeros(Q+1,1), [0; sin(2*xq(2:Q)); 0]];
methods = {'ETDRK2','ETDRK4','SBDF4'};
taus    = [0.02,    0.02,    0.01];
trec    = (0:0.5:T)'; tsnaps = [40, 58, 200];

D_std = zeros(numel(trec), 3); D_iso = D_std;
S_std = cell(3, 3); S_iso = cell(3, 3); Tstar = nan(3, 1);
for m = 1:3
    nm = methods{m}; tau = taus(m);
    [D_std(:,m), S_std(m,:), Tstar(m)] = run_dir_track(nm, dst1, idst1, @(u) u.^3, Lk_f, u0, tau, T, trec, tsnaps, false, N);
    [D_iso(:,m), S_iso(m,:), ~]        = run_dir_track(nm, toSpec_i, toPhys_i, cube_i, Lk_i, U0, tau, T, trec, tsnaps, true, N, assem_i);
    fprintf('  %-7s (tau=%.2f): std T* = %5.2f (maxD=%.2e, u_end(pi/2)=%+.3f) | iso maxD = %g\n', ...
        nm, tau, Tstar(m), max(D_std(:,m)), S_std{m,3}(Q+1), max(D_iso(:,m)));
end

% Plot
set(groot,'defaultAxesFontSize',9.5,'defaultLineLineWidth',1.2, ...
    'defaultAxesTickLabelInterpreter','latex','defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex','defaultAxesBox','on');
f = figure('Units','centimeters','Position',[2 2 18 9.8]);
cols = [0.85 0.16 0.16; 0.90 0.45 0.00; 0.55 0.20 0.70];
c_iso = [0.00 0.35 0.80];

axA = axes('Position',[0.075 0.56 0.89 0.37]); hold on;
hh = gobjects(4,1);
for m = 1:3
    hh(m) = semilogy(axA, trec, max(D_std(:,m), 1e-19), '-', 'Color', cols(m,:), 'LineWidth', 1.4);
end
hh(4) = semilogy(axA, trec, 1e-19*ones(size(trec)), '-', 'Color', c_iso, 'LineWidth', 2.0);
set(axA,'YScale','log'); xlim([0 T]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('midpoint defect $D(u_n)$');
title('(a) Homogeneous Dirichlet ($\kappa=0.5$, $u_0=\sin 2x$): standard DST-I vs.\ isotypic Dirichlet');
legend(hh, {'ETDRK2 (DST-I)','ETDRK4 (DST-I)','SBDF4 (DST-I)','isotypic Dirichlet (all 3): $D\equiv0$'}, ...
       'Location','east','FontSize',7.2);

xb = [0.075, 0.385, 0.695]; wb = 0.285; hb = 0.33; yb = 0.095;
tstr = {'(b) $t = 40$ (pre-breakdown)', '(c) $t = 58$ (breaking into single arch)', '(d) $t = 200$ (wrong arch vs.\ exact odd)'};
for p = 1:3
    ax = axes('Position',[xb(p) yb wb hb]); hold on;
    xline(pi/2, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.85);
    yline(0,    ':', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
    for m = 1:3, plot(ax, x, S_std{m,p}, '-', 'Color', cols(m,:), 'LineWidth', 1.25); end
    plot(ax, x, S_iso{2,p}, '--', 'Color', c_iso, 'LineWidth', 1.9);
    xlim([0 pi]); ylim([-1.15 1.15]);
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
function u = run_dir(nm, toSpec, toPhys, cube, Lk, u0, tau, nst, isIso, N)
[~, S, ~] = run_dir_track(nm, toSpec, toPhys, cube, Lk, u0, tau, nst*tau, [0; nst*tau], nst*tau, isIso, N);
u = S{1};
end

function [Drec, Snaps, tb] = run_dir_track(nm, toSpec, toPhys, cube, Lk, u0, tau, T, trec, tsnaps, isIso, N, assem)
hasAssem = (nargin >= 13);
nst = round(T/tau);
Drec = zeros(numel(trec),1); Snaps = cell(1,numel(tsnaps)); tb = NaN;
r_step = max(1, round((trec(2)-trec(1))/tau)); s_steps = round(tsnaps/tau);
Nop = @(uh) toSpec(toPhys(uh) - cube(toPhys(uh)));
z = -tau*Lk; E = exp(z); E2 = exp(z/2);
Q1 = tau*phik(z,1); Q2 = tau*phik(z,2); Q = (tau/2)*phik(z/2,1);
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
rec(u0, 0); uh = toSpec(u0);
switch nm
    case 'ETDRK2'
        for n = 1:nst
            Nu = Nop(uh); a = E.*uh + Q1.*Nu; uh = a + Q2.*(Nop(a) - Nu);
            if mod(n,r_step)==0 || any(s_steps==n), rec(toPhys(uh), n); end
        end
    case 'ETDRK4'
        for n = 1:nst
            Nu = Nop(uh); a=E2.*uh+Q.*Nu; Na=Nop(a); b=E2.*uh+Q.*Na; Nb=Nop(b);
            c=E2.*a+Q.*(2*Nb-Nu); Nc=Nop(c); uh=E.*uh+f1.*Nu+2*f2.*(Na+Nb)+f3.*Nc;
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
end
end

function un = etd4s(uh, Nop, E, E2, Q, f1, f2, f3)
Nu=Nop(uh); a=E2.*uh+Q.*Nu; Na=Nop(a); b=E2.*uh+Q.*Na; Nb=Nop(b);
c=E2.*a+Q.*(2*Nb-Nu); Nc=Nop(c); un=E.*uh+f1.*Nu+2*f2.*(Na+Nb)+f3.*Nc;
end

function p = phik(z, k)
Mq=64; r=exp(1i*pi*((1:Mq)-.5)/Mq); sz=size(z); zz=z(:); p=zeros(size(zz));
big=abs(zz)>1; zb=zz(big);
if k==1, p(big)=(exp(zb)-1)./zb; else, p(big)=(exp(zb)-1-zb)./zb.^2; end
zs=zz(~big);
if ~isempty(zs), Zc=zs+r; if k==1, F=(exp(Zc)-1)./Zc; else, F=(exp(Zc)-1-Zc)./Zc.^2; end; p(~big)=real(mean(F,2)); end
p=reshape(p,sz);
end

function [f1,f2,f3] = etd4coef(z, tau)
Mq=64; r=exp(1i*pi*((1:Mq)-.5)/Mq); sz=size(z); zz=z(:); f1=zeros(size(zz)); f2=f1; f3=f1;
big=abs(zz)>2; zb=zz(big);
f1(big)=(-4-zb+exp(zb).*(4-3*zb+zb.^2))./zb.^3;
f2(big)=(2+zb+exp(zb).*(-2+zb))./zb.^3;
f3(big)=(-4-3*zb-zb.^2+exp(zb).*(4-zb))./zb.^3;
zs=zz(~big);
if ~isempty(zs), Zc=zs+2*r; f1(~big)=real(mean((-4-Zc+exp(Zc).*(4-3*Zc+Zc.^2))./Zc.^3,2)); f2(~big)=real(mean((2+Zc+exp(Zc).*(-2+Zc))./Zc.^3,2)); f3(~big)=real(mean((-4-3*Zc-Zc.^2+exp(Zc).*(4-Zc))./Zc.^3,2)); end
f1=reshape(tau*f1,sz); f2=reshape(tau*f2,sz); f3=reshape(tau*f3,sz);
end
