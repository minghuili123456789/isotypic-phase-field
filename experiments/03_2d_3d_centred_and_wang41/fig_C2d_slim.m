% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Compact 1x3 figure for 2D C2 point-symmetric pattern (N=192, T=200):
% Replaces visually-identical 2-row grid with sleek 1-row (t=0, t=200, defect history)
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp2d3d_b.mat'));
R = struct;
for j = 1:numel(S.out)
    if strcmp(S.out{j}.name,'C2d'), R.(S.out{j}.variant) = S.out{j}; end
end
N = size(R.iso.snap{1},1); x = (0:N-1)*2*pi/N; x(x>=pi) = x(x>=pi)-2*pi; [xs,ix] = sort(x);

f = figure('Units','centimeters','Position',[2 2 18 4.9]);
% (a) t = 0
ax1 = axes('Position',[0.04, 0.14, 0.21, 0.74]);
imagesc(xs, xs, R.iso.snap{1}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
set(ax1,'XTick',[],'YTick',[]); title('(a) $t=0$ ($C_2$ datum)','FontSize',9);

% (b) t = 200
s200 = find(abs(R.iso.snapT - 200) < 1e-9, 1);
ax2 = axes('Position',[0.28, 0.14, 0.21, 0.74]);
imagesc(xs, xs, R.iso.snap{s200}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
set(ax2,'XTick',[],'YTick',[]); title('(b) $t=200$ ($C_2$ stripe)','FontSize',9);
colormap(parula);

% (c) Defect history 0 <= t <= 200
ax3 = axes('Position',[0.575, 0.18, 0.405, 0.70]); hold on;
cols = lines(3);
semilogy(ax3, R.fourier.t, max(R.fourier.D, 1e-19), '-',  'Color', cols(1,:), 'LineWidth', 1.4);
semilogy(ax3, R.zmode.t,   max(R.zmode.D,   1e-19), '--', 'Color', cols(2,:), 'LineWidth', 1.4);
semilogy(ax3, R.iso.t,     max(R.iso.D,     1e-19), '-',  'Color', cols(3,:), 'LineWidth', 1.8);
set(ax3,'YScale','log');
xlim([0 200]); ylim([1e-20 1e-10]); grid on;
xlabel('$t$','FontSize',9); ylabel('$D(u_n)$','FontSize',9);
title('(c) Neutral linear drift vs.\ exact bitwise invariance','FontSize',9);
legend({'standard FFT','zero-mode suppression','isotypic ($D\equiv0$)'}, ...
       'Location','east','FontSize',6.8);

exportgraphics(f, fullfile(here,'..','..','figs','fig_C2d.pdf'),'ContentType','vector');
fprintf('Wrote slim 1x3 fig_C2d.pdf\n');
