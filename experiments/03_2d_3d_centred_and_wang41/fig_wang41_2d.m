% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: Exact 2D comparison on Wang (2023) Example 4.1 (u0 = sin(2x)sin(2y),
% kappa = 0.1, Np = 256, S = 1, Te = 1, Ns = 10, Niter = 20 -- ALL conditions identical!).
% 3 rows x 4 snapshots (t = 20, 80, 100, 120) + 2 right-hand panels (defect & energy).
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp_wang41_2d.mat'));

N = S.N; x = (0:N-1)*2*pi/N; x(x>=pi) = x(x>=pi) - 2*pi; [xs, ix] = sort(x);
tsel = [20, 80, 100, 120];
isel = arrayfun(@(t) find(S.tsnaps == t, 1), tsel);

f = figure('Units','centimeters','Position',[2 2 18 10.8]);
left0 = 0.068; gapx = 0.014; w = 0.140; h = 0.275;
yrows = [0.665, 0.355, 0.045];
Snaps = {S.snap_std, S.snap_iso, S.snap_art2};
Rlabs = {{'Wang ESP','(standard FFT)'}, ...
         {'$A_{\mathrm R}$ isotypic ESP','(odd--odd only)'}, ...
         {'$A_{\mathrm R}\!\times\!\mathcal{T}_2$ ESP','(ours, same params)'}};

for r = 1:3
    for k = 1:numel(tsel)
        s = isel(k);
        ax = axes('Position',[left0 + (k-1)*(w+gapx), yrows(r), w, h]);
        imagesc(xs, xs, Snaps{r}{s}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
        set(ax,'XTick',[],'YTick',[]);
        if r == 1, title(sprintf('$t = %g$', tsel(k)), 'FontSize', 9.5); end
        if k == 1, ylabel(Rlabs{r}, 'FontSize', 8.3); end
    end
end
colormap(parula);

c_std  = [0.85 0.16 0.16];
c_zm   = [0.90 0.45 0.00];
c_ar   = [0.00 0.60 0.35];
c_art2 = [0.00 0.35 0.80];

% Top-right: Defect evolution 0 <= t <= 120
axRT = axes('Position',[0.745, 0.565, 0.240, 0.375]); hold on;
h1 = semilogy(axRT, S.tvec, max(S.D_std,  1e-19), '-',  'Color', c_std,  'LineWidth', 1.4);
h2 = semilogy(axRT, S.tvec, max(S.D_zm,   1e-19), '--', 'Color', c_zm,   'LineWidth', 1.35);
h3 = semilogy(axRT, S.tvec, max(S.D_iso,  1e-19), '-',  'Color', c_ar,   'LineWidth', 1.6);
h4 = semilogy(axRT, S.tvec, max(S.D_art2, 1e-19), '--', 'Color', c_art2, 'LineWidth', 1.6);
xline(20, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility','off');
text(23, 3e-7, {'Wang Fig.\,5','stopped','at $t=20$'}, 'FontSize', 6.5, 'Color', [0.35 0.35 0.35]);
set(axRT,'YScale','log');
xlim([0 120]); ylim([1e-20 10]); grid on;
xlabel('$t$','FontSize',9); ylabel('defect $D(u_n)$','FontSize',9);
title('Symmetry defect','FontSize',9.5);
legend([h1 h2 h3 h4], {'Wang ESP','ESP+zero-mode','$A_{\mathrm R}$ ($D_{\mathrm R}\!\equiv\!0$)','$A_{\mathrm R}\!\times\!\mathcal{T}_2$ ($D\!\equiv\!0$)'}, ...
       'Location','southeast','FontSize',6.0);

% Bottom-right: Energy evolution 0 <= t <= 120 (extending Wang Fig. 7 left!)
axRB = axes('Position',[0.745, 0.095, 0.240, 0.355]); hold on;
plot(axRB, S.tvec, S.E_std,  '-',  'Color', c_std,  'LineWidth', 1.4);
plot(axRB, S.tvec, S.E_zm,   '--', 'Color', c_zm,   'LineWidth', 1.3);
plot(axRB, S.tvec, S.E_iso,  '-.', 'Color', c_ar,   'LineWidth', 1.4);
plot(axRB, S.tvec, S.E_art2, '-',  'Color', c_art2, 'LineWidth', 1.7);
xline(20, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility','off');
xlim([0 120]); ylim([3.2 5.2]); grid on;
xlabel('$t$','FontSize',9); ylabel('energy $E(u_n)$','FontSize',9);
title('Energy (cf.\ Wang Fig.\,7)','FontSize',9.5);
text(23, 3.65, {'false drop','after $t\approx 85$'}, 'FontSize', 6.5, 'Color', c_std);
text(52, 4.45, {'exact plateau','$E=4.2351$'}, 'FontSize', 6.5, 'Color', c_art2);

exportgraphics(f, fullfile(here,'..','..','figs','fig_wang41_2d.pdf'),'ContentType','vector');
fprintf('Wrote 3-row fig_wang41_2d.pdf\n');
