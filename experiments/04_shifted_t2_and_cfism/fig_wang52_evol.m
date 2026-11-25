% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: 1D long-time shifted evolution comparison on Wang (2023) Ex. 5.2(a)
% u0(x) = sin(x + pi/4), kappa = 1, N = 512, T = 1000.
% Visually demonstrates that ALL existing methods fail while our non-interventional
% T2 isotypic ESP & CF-ISM preserve the exact shifted odd symmetry bitwise (D === 0).
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp_wang52_evol.mat'));

x = S.x; tsnaps = S.tsnaps;
tshow = [30, 43, 100];
ishow = arrayfun(@(t) find(tsnaps == t, 1), tshow);
i1000 = find(tsnaps == 1000, 1);

f = figure('Units','centimeters','Position',[2 2 18 11.6]);
c_std  = [0.85 0.16 0.16];   % red
c_esp  = [0.90 0.45 0.00];   % orange
c_lqty = [0.55 0.20 0.70];   % purple
c_t2   = [0.00 0.35 0.80];   % deep blue
c_cf   = [0.00 0.62 0.35];   % emerald green
c_ref  = [0.60 0.60 0.60];   % reference gray

Snaps = {S.S1, S.S2, S.S3, S.S4};
Cols  = {c_std, c_esp, c_lqty, c_t2};
Rlabs = {{'standard','ETDRK4'}, {'Wang ESP','($S\!=\!1,T_e\!=\!1$)'}, ...
         {'LQTY','parity filter'}, {'$\mathcal{T}_2$ isotypic','ESP (ours)'}};

% Left 4x3 grid of snapshots
left0 = 0.095; w = 0.162; gapx = 0.017;
top0  = 0.935; h = 0.180; gapy = 0.036;
for r = 1:4
    yb = top0 - r*h - (r-1)*gapy;
    for c = 1:3
        xb = left0 + (c-1)*(w+gapx);
        ax = axes('Position',[xb, yb, w, h]); hold on;
        u_ref = S.S4{ishow(c)};
        u_cur = Snaps{r}{ishow(c)};
        xline(-pi/4, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.85, 'HandleVisibility','off');
        yline(0,     ':', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.75, 'HandleVisibility','off');
        if r < 4
            plot(x, u_ref, '--', 'Color', c_ref, 'LineWidth', 0.95);
        end
        plot(x, u_cur, '-', 'Color', Cols{r}, 'LineWidth', 1.45);
        xlim([-pi pi]); ylim([-1.18 1.18]);
        set(ax,'YTick',[-1 0 1],'XTick',[-pi 0 pi]);
        if c > 1, set(ax,'YTickLabel',[]); end
        if r < 4
            set(ax,'XTickLabel',[]);
        else
            set(ax,'XTickLabel',{'$-\pi$','$0$','$\pi$'});
            xlabel('$x$','FontSize',9);
        end
        if r == 1, title(sprintf('$t = %g$', tshow(c)), 'FontSize', 9.5); end
        if c == 1, ylabel(Rlabs{r}, 'FontSize', 8.3); end
        grid on; set(ax,'GridAlpha',0.15);
    end
end

% Top-right: Shifted-symmetry defect history 0 <= t <= 1000
axRT = axes('Position',[0.705, 0.580, 0.280, 0.355]); hold on;
t = S.tgrid;
hD1 = semilogy(axRT, t, max(S.D1, 1e-19),    '-',  'Color', c_std,  'LineWidth', 1.3);
hD2 = semilogy(axRT, t, max(S.D2, 1e-19),    '-',  'Color', c_esp,  'LineWidth', 1.4);
hD3 = semilogy(axRT, t, max(S.D3, 1e-19),    '-.', 'Color', c_lqty, 'LineWidth', 1.3);
hD4 = semilogy(axRT, t, max(S.D4_t2, 1e-19), '-',  'Color', c_t2,   'LineWidth', 1.7);
hD5 = semilogy(axRT, t, max(S.D5, 1e-19),    '--', 'Color', c_cf,   'LineWidth', 1.4);
set(axRT,'YScale','log');
xlim([0 1000]); ylim([1e-20 10]); grid on;
xlabel('$t$','FontSize',9); ylabel('defect $D(u_n)$','FontSize',9);
title('Long-time defect ($0 \le t \le 1000$)','FontSize',9.5);
legend([hD1 hD2 hD3 hD4 hD5], ...
       {'std ETDRK4','Wang ESP ($S=1$)','LQTY filter', ...
        '$\mathcal{T}_2$ ESP ($D\equiv0$)','CF-ISM ($D\equiv0$)'}, ...
        'Location','east','FontSize',6.4);

% Bottom-right: Zoomed long-time snapshot at T = 1000
axRB = axes('Position',[0.705, 0.095, 0.280, 0.345]); hold on;
xline(-pi/4, ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9, 'HandleVisibility','off');
yline(0,     ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility','off');
hZ3 = plot(x, S.S3{i1000}, '-.', 'Color', c_lqty, 'LineWidth', 1.3);
hZ4 = plot(x, S.S4{i1000}, '-',  'Color', c_t2,   'LineWidth', 1.7);
hZ5 = plot(x, S.S5{i1000}, '--', 'Color', c_cf,   'LineWidth', 1.35);
xlim([-pi pi]); ylim([-0.044 0.044]);
set(axRB,'XTick',[-pi -pi/4 0 pi],'XTickLabel',{'$-\pi$','$-\frac{\pi}{4}$','$0$','$\pi$'});
grid on; xlabel('$x$','FontSize',9); ylabel('$u(x,1000)$','FontSize',9);
title('Zoom at $T = 1000$','FontSize',9.5);
text(-2.98, 0.035, 'Std \& ESP: $u\equiv +1$ (off scale)', 'FontSize', 7.0, ...
     'Color', c_std, 'BackgroundColor', [1 1 1 0.88]);
legend([hZ3 hZ4 hZ5], {'LQTY (wrong phase)','$\mathcal{T}_2$ ESP (ours)','CF-ISM (ours)'}, ...
       'Location','southeast','FontSize',6.3);

exportgraphics(f, fullfile(here,'..','..','figs','fig_wang52_evol.pdf'),'ContentType','vector');
fprintf('Wrote polished fig_wang52_evol.pdf\n');
