% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: Time-evolving departure of 7 high-order integrators vs exact
% long-time invariance of their drop-in isotypic counterparts (0 <= t <= 1000).
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp_highorder_evol.mat'));

x = S.xs; H = S.H;
trec = S.trec; tlong = S.tlong;
nm = numel(H);
cols = [0.85 0.16 0.16;   % BDF2-SAV (red)
        0.90 0.45 0.00;   % ARS443   (orange)
        0.72 0.52 0.04;   % SBDF3    (ochre)
        0.55 0.20 0.70;   % SBDF4    (purple)
        0.10 0.60 0.55;   % IFRK4    (teal)
        0.80 0.20 0.55;   % SDC4     (magenta)
        0.45 0.45 0.45];  % Strang2  (gray)
c_iso = [0.00 0.35 0.80]; % deep blue for isotypic
labs = {'BDF2-SAV','ARS(4,4,3)','SBDF3','SBDF4','IFRK4','SDC4','Strang2'};

f = figure('Units','centimeters','Position',[2 2 18 10.6]);

% Top-left (a): Defect growth D(u_n) on [0, 75]
axA = axes('Position',[0.075, 0.56, 0.41, 0.38]); hold on;
hh = gobjects(nm+2, 1);
for m = 1:nm
    ls = '-'; if m == 7, ls = ':'; end
    hh(m) = semilogy(axA, trec, max(H{m}.D_f, 1e-19), ls, 'Color', cols(m,:), 'LineWidth', 1.25);
end
hh(nm+1) = semilogy(axA, trec, 2e-19*exp(trec), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
hh(nm+2) = semilogy(axA, trec, 1e-19*ones(size(trec)), '-', 'Color', c_iso, 'LineWidth', 2.0);
set(axA,'YScale','log');
xlim([0 75]); ylim([1e-20 10]); grid on;
xlabel('$t$','FontSize',9.5); ylabel('defect $D(u_n)$','FontSize',9.5);
title('(a) Symmetry defect growth ($0 \le t \le 75$)','FontSize',10);
legend(hh, [labs, {'$\propto \mathrm{e}^t$','isotypic ($D\equiv0$)'}], ...
       'Location','southeast','NumColumns',2,'FontSize',6.0);

% Top-right (b): Long-time horizon [0, 1000]
axB = axes('Position',[0.565, 0.56, 0.415, 0.38]); hold on;
for m = 1:nm
    ls = '-'; if m == 7, ls = ':'; end
    semilogy(axB, tlong, max(H{m}.D_fl, 1e-19), ls, 'Color', cols(m,:), 'LineWidth', 1.25);
end
semilogy(axB, tlong, 1e-19*ones(size(tlong)), '-', 'Color', c_iso, 'LineWidth', 2.0);
set(axB,'YScale','log');
xlim([0 1000]); ylim([1e-20 10]); grid on;
xlabel('$t$','FontSize',9.5); ylabel('defect $D(u_n)$','FontSize',9.5);
title('(b) Long-time horizon ($0 \le t \le 1000$)','FontSize',10);
text(140, 2e-4, {'All 6 standard high-order schemes','collapse to $u\equiv -1$ by $t\approx 48$'}, ...
     'FontSize',7.5,'Color',[0.75 0.12 0.12]);
text(140, 4e-18, 'All 7 isotypic schemes: $D \equiv 0$ through $T=1000$', ...
     'FontSize',7.4,'Color',c_iso);

% Bottom row: 3 solution snapshots at t = 35, t = 47 (mid-departure), t = 1000
% Note: S.S47_f / S.S47_i holds the mid-collapse snapshot at t = 47
xb = [0.075, 0.385, 0.695]; wb = 0.285; hb = 0.33; yb = 0.095;
tstr = {'(c) $t = 35$ (pre-breakdown)', '(d) $t = 47$ (peeling off toward $-1$)', '(e) $t = 1000$ ($20\times$ farther in time)'};
for p = 1:3
    ax = axes('Position',[xb(p), yb, wb, hb]); hold on;
    xline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    yline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    for m = 1:6
        if p == 1, uf = H{m}.S_f{1};
        elseif p == 2, uf = S.S47_f{m};
        else, uf = H{m}.S_f{4};
        end
        plot(ax, x, uf, '-', 'Color', cols(m,:), 'LineWidth', 1.2);
    end
    if p == 1, ui = H{4}.S_i{1};
    elseif p == 2, ui = S.S47_i;
    else, ui = H{4}.S_i{4};
    end
    plot(ax, x, ui, '--', 'Color', c_iso, 'LineWidth', 1.9);
    xlim([-pi pi]); ylim([-1.18 1.18]);
    set(ax,'XTick',[-pi 0 pi],'XTickLabel',{'$-\pi$','$0$','$\pi$'},'YTick',[-1 0 1]);
    if p > 1, set(ax,'YTickLabel',[]); end
    grid on; set(ax,'GridAlpha',0.15);
    xlabel('$x$','FontSize',9.5);
    if p == 1, ylabel('$u(x,t)$','FontSize',9.5); end
    title(tstr{p},'FontSize',9.5);
    if p == 3
        % Clean zoomed inset in the empty upper half of panel (e)
        axIn = axes('Position',[xb(3)+0.055, yb+0.165, 0.195, 0.125]); hold on;
        xline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7);
        yline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7);
        plot(axIn, x, ui, '-', 'Color', c_iso, 'LineWidth', 1.35);
        xlim([-pi pi]); ylim([-0.035 0.035]);
        set(axIn,'XTick',[-pi 0 pi],'XTickLabel',{'$-\pi$','$0$','$\pi$'}, ...
                 'YTick',[-0.03 0 0.03],'FontSize',6.3,'Box','on','Color',[0.98 0.99 1]);
        title('zoom: isotypic ($D\equiv0$)','FontSize',6.6,'Color',c_iso);
    end
end

exportgraphics(f, fullfile(here,'..','..','figs','fig_highorder_1d.pdf'),'ContentType','vector');
fprintf('Wrote polished fig_highorder_1d.pdf\n');
