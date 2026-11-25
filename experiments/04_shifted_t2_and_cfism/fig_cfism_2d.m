% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: 2D irrational-shift mixed-harmonic X-junction (Standard FFT vs CF-ISM)
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp_cfism_2d.mat'));
tsel = S.snaps; nt = numel(tsel);
N = size(S.snap_cf{1},1); x = (0:N-1)*2*pi/N; x(x>=pi) = x(x>=pi)-2*pi; [xs,ix] = sort(x);

f = figure('Units','centimeters','Position',[2 2 18 7.8]);
left0 = 0.045; gapx = 0.015; w = 0.145; h = 0.39;
yrows = [0.53, 0.07];
for k = 1:nt
    ax = axes('Position',[left0 + (k-1)*(w+gapx), yrows(1), w, h]);
    imagesc(xs, xs, S.snap_std{k}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
    set(ax,'XTick',[],'YTick',[]); title(sprintf('$t=%g$',tsel(k)),'FontSize',9.5);
    if k == 1, ylabel('standard','FontSize',9.5); end

    ax = axes('Position',[left0 + (k-1)*(w+gapx), yrows(2), w, h]);
    imagesc(xs, xs, S.snap_cf{k}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
    set(ax,'XTick',[],'YTick',[]);
    if k == 1, ylabel('CF-ISM','FontSize',9.5); end
end
colormap(parula);

axR = axes('Position',[0.76, 0.14, 0.225, 0.78]);
cols = lines(3);
semilogy(axR, S.t, max(S.D_std, 1e-19),    '-',  'Color', cols(1,:), 'LineWidth', 1.3); hold on;
semilogy(axR, S.t, max(S.mean_std, 1e-19), '--', 'Color', cols(2,:), 'LineWidth', 1.3);
semilogy(axR, S.t, max(S.D_cf, 1e-19),     '-',  'Color', cols(3,:), 'LineWidth', 1.6);
xlim([0 4]); ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('defect / $|\mathrm{mean}|$');
legend({'standard $D$','standard $|\mathrm{mean}|$','CF-ISM ($D\equiv0$)'},'Location','southeast','FontSize',6.8);
exportgraphics(f, fullfile(here,'..','..','figs','fig_cfism_2d.pdf'),'ContentType','vector');
