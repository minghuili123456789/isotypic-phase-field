% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Unified 1x3 Faithfulness figure merging fig_delta + fig_vpa (zero duplication!):
%  (a) Computed breakdown time T*(delta) over delta = 1e-2..1e-200 (iso vs std double/single vs LQTY)
%  (b) Defect histories vs 120-digit VPA reference for delta = 1e-10, 1e-30, 1e-60
%  (c) Relative deviation from 120-digit VPA reference in the linear regime
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
L = load(fullfile(here,'..','..','results','exp_delta.mat'));
extra_file = fullfile(here,'..','..','results','exp_delta_extra.mat');
if exist(extra_file,'file')
    load(extra_file,'extra');
    L.out = [L.out(~cellfun(@(o) strcmp(o.variant,'iso') && o.delta==1e-200, L.out)), {extra}];
end
load(fullfile(here,'..','..','results','ref_vpa.mat'),'ref');

f = figure('Units','centimeters','Position',[2 2 18 5.6]);

% (a) T*(delta) vs -log10(delta)
ax1 = axes('Position',[0.062, 0.175, 0.255, 0.72]); hold on;
vars  = {'iso','fourier','single','filter'};
names = {'isotypic (double)','standard (double)','standard (single)','LQTY filter'};
mk    = {'o','s','^','x'}; cols4 = lines(4);
for v = 1:4
    d = []; T = [];
    for j = 1:numel(L.out)
        o = L.out{j};
        if strcmp(o.case,'k1') && strcmp(o.variant,vars{v}) && o.delta > 0
            tb = o.t(find(o.D > 0.1,1)); if isempty(tb), tb = inf; end
            d(end+1) = o.delta; T(end+1) = tb; %#ok<SAGROW>
        end
    end
    [d,ix] = sort(d,'descend'); T = T(ix);
    Tp = T; Tp(isinf(Tp)) = 520 + 12*v;
    plot(ax1, -log10(d), Tp, [mk{v} '-'], 'Color', cols4(v,:), 'MarkerSize', 4.2, 'LineWidth', 1.1);
end
plot(ax1, [0 205],[500 500],'k:','HandleVisibility','off');
text(ax1, 30, 512, 'not reached by $T=400$', 'FontSize', 6.3);
xlim([0 205]); ylim([0 585]); grid on;
xlabel('$-\log_{10}\delta$','FontSize',9); ylabel('$T_\ast$','FontSize',9);
legend(names,'Location','southeast','FontSize',6.0);
title('(a) Breakdown time vs.\ $\delta$','FontSize',9.2);

% (b) Defect histories vs 120-digit VPA
ax2 = axes('Position',[0.395, 0.175, 0.255, 0.72]); hold on;
cols3 = lines(3);
hV = gobjects(3,1);
for q = 1:numel(ref)
    hV(q) = semilogy(ax2, ref(q).t, ref(q).D, '-', 'Color', cols3(q,:), 'LineWidth', 1.3);
    semilogy(ax2, ref(q).t, ref(q).Diso, 'o', 'Color', cols3(q,:), 'MarkerSize', 3.0);
    semilogy(ax2, ref(q).t, max(ref(q).Dfou,1e-300), ':', 'Color', cols3(q,:), 'LineWidth', 1.5);
end
set(ax2,'YScale','log');
ylim([1e-68 10]); grid on;
xlabel('$t$','FontSize',9); ylabel('$D(u_n)$','FontSize',9);
title('(b) VPA (---), iso ($\circ$), std ($\cdot\cdot$)','FontSize',9.2);
lg = arrayfun(@(r) sprintf('$\\delta=10^{%d}$',round(log10(r.delta))), ref, 'UniformOutput', false);
legend(hV, lg, 'Location','southeast','FontSize',6.4);

% (c) Relative deviation from 120-digit VPA
ax3 = axes('Position',[0.725, 0.175, 0.255, 0.72]); hold on;
hR = gobjects(3,1);
for q = 1:numel(ref)
    lin = ref(q).D < 1e-2;
    hR(q) = semilogy(ax3, ref(q).t(lin), abs(ref(q).Diso(lin)-ref(q).D(lin))./ref(q).D(lin), ...
        'o-', 'Color', cols3(q,:), 'MarkerSize', 3.0, 'LineWidth', 1.1);
    semilogy(ax3, ref(q).t(lin), abs(ref(q).Dfou(lin)-ref(q).D(lin))./ref(q).D(lin), ...
        ':', 'Color', cols3(q,:), 'LineWidth', 1.5);
end
set(ax3,'YScale','log');
ylim([1e-17 1e50]); grid on;
xlabel('$t$','FontSize',9); ylabel('rel.\ dev.\ from VPA','FontSize',9);
title('(c) Relative error vs.\ VPA','FontSize',9.2);
text(ax3, 35, 1e33, 'std ($\cdot\cdot$)', 'FontSize', 6.8, 'Color', [0.3 0.3 0.3]);
text(ax3, 35, 1e-9, 'iso ($\circ$)', 'FontSize', 6.8, 'Color', [0.1 0.3 0.7]);

exportgraphics(f, fullfile(here,'..','..','figs','fig_faithful.pdf'),'ContentType','vector');
fprintf('Wrote unified fig_faithful.pdf\n');
