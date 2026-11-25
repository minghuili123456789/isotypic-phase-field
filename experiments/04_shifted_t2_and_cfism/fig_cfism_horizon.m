% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: CF-ISM horizon scaling T*(d0) (polished layout)
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
L1 = load(fullfile(here,'..','..','results','exp_shifted_1d.mat'),'horizon');
H = L1.horizon;

f = figure('Units','centimeters','Position',[2 2 17.5 6.8]);
d0 = [H.d0]; Ts = [H.Tstar]; fin = ~isnan(Ts);
pp = polyfit(d0(fin), Ts(fin), 1);

% (a) Breakdown time T*(d0)
subplot(1,2,1);
plot(d0(fin), Ts(fin), 'bo', 'MarkerFaceColor',[.2 .4 .9], 'MarkerSize', 5.5); hold on;
xx = linspace(0, 310, 100);
plot(xx, polyval(pp, xx), 'b-', 'LineWidth', 1.2);
plot(d0(~fin), 740, 'r^', 'MarkerFaceColor','r', 'MarkerSize', 7.5);
plot([0 340], [700 700], 'k:');
text(15, 735, '$T_\ast=\infty$ ($D\equiv0$)', 'FontSize', 7.5, 'Color', [0.75 0 0]);
plot([0 340], [37.4 37.4], 'k--', 'LineWidth', 1.0);
text(125, 70, 'standard / ESP ($T_\ast\approx 37$--$41$)', 'FontSize', 7);
xlim([0 340]); ylim([0 790]); grid on;
xlabel('$t=0$ precision $d_0$ (double for $t>0$)');
ylabel('$T_\ast$ ($D>10^{-2}$)');
legend({'CF-ISM ($T_\ast < \infty$)', sprintf('fit: slope $%.3f \\approx \\ln 10$', pp(1)), ...
    'CF-ISM ($d_0=320$, $D\equiv0$)'}, 'Position',[0.23 0.26 0.21 0.17], 'FontSize', 6.5);
title('(a) CF-ISM horizon law (irrational shift)', 'FontSize', 8.5);

% (b) Initial transverse defect vs d0 (log10 scale so underflow at 0 is cleanly shown)
subplot(1,2,2);
ev_log = log10([H(fin).init_ev]);
plot(d0(fin), ev_log, 'ks-', 'MarkerFaceColor',[.3 .7 .4], 'MarkerSize', 5); hold on;
plot([0 340], [-307.65 -307.65], 'r--', 'LineWidth', 1.2);
text(20, -292, 'IEEE double underflow ($\approx 2.2\times10^{-308}$)', 'FontSize', 7, 'Color', 'r');
plot(320, -328, 'r^', 'MarkerFaceColor','r', 'MarkerSize', 7.5);
text(225, -328, 'exact $0$', 'FontSize', 7.5, 'Color', [0.75 0 0]);
xlim([0 340]); ylim([-345 0]); grid on;
xlabel('$t=0$ precision $d_0$');
ylabel('$\log_{10}\|v^0_\perp\|_\infty$');
title('(b) initial transverse defect at $t=0$', 'FontSize', 8.5);
exportgraphics(f, fullfile(here,'..','..','figs','fig_cfism_horizon.pdf'),'ContentType','vector');
