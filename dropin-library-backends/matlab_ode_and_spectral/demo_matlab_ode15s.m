% Worked Example: Plugging iso_spectral_backend into MATLAB's built-in ode15s!
%   PDE: du/dt = u_xx + u - u^3 on [-pi, pi), periodic BC, u0 = sin(x), N = 128, T = 150.
% Demonstrates that MATLAB's own stiff ODE solver ode15s (RelTol=1e-10, AbsTol=1e-12)
% breaks odd symmetry at T* ~ 42 and collapses to u = +/-1 when using standard fft,
% whereas changing ONLY the 3-line spatial RHS to iso_spectral_backend preserves D === 0!
clear; close all;
addpath(fileparts(mfilename('fullpath')));
outfig = fullfile(fileparts(mfilename('fullpath')), 'fig_demo_ode15s.png');

N = 128; M = N/2; T = 150;
x = (0:N-1)'*2*pi/N; xs = x; xs(xs>=pi) = xs(xs>=pi)-2*pi; [xs, ix] = sort(xs);
k = [0:N/2-1, -N/2:-1]'; K2_f = k.^2;

B = iso_spectral_backend(1, N, 'AR');

% Standard FFT RHS for MATLAB ode15s
rhs_std = @(t, u) real(ifft(-K2_f .* fft(u))) + u - u.^3;

% Drop-in Isotypic RHS for MATLAB ode15s (ONLY 3 lines changed!)
rhs_iso = @(t, y) reshape( ...
    B.toPhys(-B.K2 .* B.toSpec(reshape(y, M+1, 2))) + reshape(y, M+1, 2) - B.cube(reshape(y, M+1, 2)), ...
    [], 1);

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
tgrid = 0:1:T;

% 1. Verify non-symmetric equivalence with MATLAB ode15s
u_ns = 0.4*sin(x) + 0.3*cos(2*x) + 0.15*sin(3*x + 0.7);
[~, Yf_ns] = ode15s(rhs_std, [0 2], u_ns, opts);
[~, Yi_ns] = ode15s(rhs_iso, [0 2], reshape(B.split(u_ns), [], 1), opts);
u_iso_ns = B.assemble(reshape(Yi_ns(end,:).', M+1, 2));
fprintf('MATLAB ode15s non-symmetric equivalence (t=2): %.2e\n', max(abs(Yf_ns(end,:).' - u_iso_ns)));

% 2. Symmetric run u0 = sin(x) to T = 150
U0 = [zeros(M+1,1), sin((0:M)'*2*pi/N)]; U0([1 end],2) = 0;
u0 = B.assemble(U0);

[~, Yf] = ode15s(rhs_std, tgrid, u0, opts);
[~, Yi] = ode15s(rhs_iso, tgrid, reshape(U0, [], 1), opts);

D_std = zeros(numel(tgrid), 1); D_iso = D_std;
for j = 1:numel(tgrid)
    Uc = B.split(Yf(j,:).');  D_std(j) = max(abs(Uc(:,1)));
    Vi = reshape(Yi(j,:).', M+1, 2); D_iso(j) = max(abs(Vi(:,1)));
end
tb = tgrid(find(D_std > 1e-2, 1));
fprintf('MATLAB ode15s (u0=sin x, N=%d, T=%d): std fft T* = %g (maxD=%.2e) | drop-in iso maxD = %g\n', ...
    N, T, tb, max(D_std), max(D_iso));

% Plot result
set(groot,'defaultAxesFontSize',9.5,'defaultLineLineWidth',1.3, ...
    'defaultAxesTickLabelInterpreter','latex','defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex','defaultAxesBox','on');
f = figure('Units','centimeters','Position',[2 2 16 6.0]);
subplot(1,2,1);
semilogy(tgrid, max(D_std, 1e-19), 'r-', 'LineWidth', 1.5); hold on;
semilogy(tgrid, max(D_iso, 1e-19), 'b-', 'LineWidth', 2.0);
xlim([0 T]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('odd-symmetry defect $D(u)$');
title('(a) MATLAB \texttt{ode15s}: defect history');
legend({'standard \texttt{fft}\ RHS','drop-in \texttt{iso\_spectral\_backend}\ ($D\equiv0$)'}, ...
       'Location','southeast','FontSize',7.0);

subplot(1,2,2);
u_end_f = Yf(end,:).';
u_end_i = B.assemble(reshape(Yi(end,:).', M+1, 2));
plot(xs, u_end_f(ix), 'r-', 'LineWidth', 1.6); hold on;
plot(xs, u_end_i(ix), 'b--', 'LineWidth', 1.9);
xlim([-pi pi]); ylim([-1.15 1.15]); grid on;
set(gca,'XTick',[-pi 0 pi],'XTickLabel',{'$-\pi$','$0$','$\pi$'});
xlabel('$x$'); ylabel('$u(x,150)$');
title('(b) Final profile at $T = 150$');
legend({'standard \texttt{ode15s}\ (collapsed)','isotypic \texttt{ode15s}\ (exact odd)'}, ...
       'Location','best','FontSize',7.2);

exportgraphics(f, outfig, 'Resolution', 180);
fprintf('Saved %s\n', outfig);
