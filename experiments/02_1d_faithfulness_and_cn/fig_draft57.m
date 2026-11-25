% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: Exact reproduction of user's Draft Sec 5.7 (eps=1, N=4096, CN, T=100)
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
load(fullfile(here,'..','..','results','exp_draft57.mat'),'R_std','R_zm','R_iso','snaps');
N = numel(R_std.ufinal);
x = (0:N-1)'*2*pi/N; x(x>=pi) = x(x>=pi) - 2*pi; [xs, ix] = sort(x);

f = figure('Units','centimeters','Position',[2 2 18 6.8]);
tsel = [30 45 100];
for k = 1:3
    s = find(abs(snaps - tsel(k)) < 1e-9);
    subplot(1,4,k);
    plot(xs, R_std.snap{s}(ix), 'r-', 'LineWidth', 1.3); hold on;
    plot(xs, R_zm.snap{s}(ix), 'g-.', 'LineWidth', 1.2);
    plot(xs, R_iso.snap{s}(ix), 'b--', 'LineWidth', 1.3);
    xlim([-pi pi]); ylim([-1.08 1.08]); grid on;
    xlabel('$x$'); if k == 1, ylabel('$u(x,t)$'); end
    title(sprintf('$t=%g$', tsel(k)), 'FontSize', 9);
    if k == 3
        legend({'standard CN','zero-mode CN','isotypic CN'},'Location','southeast','FontSize',7);
    end
end
subplot(1,4,4);
semilogy(R_std.t, max(R_std.D,1e-20), 'r-', 'LineWidth', 1.3); hold on;
semilogy(R_zm.t,  max(R_zm.D,1e-20),  'g-.', 'LineWidth', 1.2);
semilogy(R_iso.t, max(R_iso.D, 10^-19.5), 'b--', 'LineWidth', 1.5);
xlim([0 100]); ylim([1e-20 10]); grid on;
xlabel('$t$'); ylabel('$D(u_n)$'); title('defect ($D_{\rm iso}\equiv0$)','FontSize',9);
exportgraphics(f, fullfile(here,'..','..','figs','fig_draft57.pdf'),'ContentType','vector');
