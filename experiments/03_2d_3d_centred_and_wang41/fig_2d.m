function fig_2d(matfile, name, outname, tsel, tmaxD)
% snapshots (standard vs isotypic) + defect history for a 2D experiment
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'..','..','core')); addpath(here); fig_common;
if ~exist(fullfile(here,'..','..','figs'),'dir'), mkdir(fullfile(here,'..','..','figs')); end
S = load(fullfile(here,'..','..','results',matfile));
R = struct;
for j = 1:numel(S.out)
    if strcmp(S.out{j}.name,name), R.(S.out{j}.variant) = S.out{j}; end
end
nt = numel(tsel);
f = figure('Units','centimeters','Position',[2 2 18 7.8]);
vars = {'fourier','iso'}; lab = {'standard','isotypic'};
N = size(R.iso.snap{1},1); x = (0:N-1)*2*pi/N; x(x>=pi) = x(x>=pi)-2*pi; [xs,ix] = sort(x);

% Explicit layout so the right-hand defect plot has plenty of room for ylabel & legend
left0 = 0.045; gapx = 0.015; w = 0.145; h = 0.39;
yrows = [0.53, 0.07];
for v = 1:2
    for k = 1:nt
        s = find(abs(R.(vars{v}).snapT - tsel(k)) < 1e-9);
        ax = axes('Position',[left0 + (k-1)*(w+gapx), yrows(v), w, h]);
        U = R.(vars{v}).snap{s};
        imagesc(xs, xs, U(ix,ix)'); axis xy equal tight; caxis([-1 1]);
        set(ax,'XTick',[],'YTick',[]);
        if v == 1, title(sprintf('$t=%g$',tsel(k)),'FontSize',9.5); end
        if k == 1, ylabel(lab{v},'FontSize',9.5); end
    end
end
colormap(parula);

axR = axes('Position',[0.76, 0.14, 0.225, 0.78]);
cols = lines(3); vv = {'fourier','zmode','iso'};
ll = {'standard','zero-mode','isotypic ($D\equiv0$)'};
for v = 1:3
    if ~isfield(R,vv{v}), continue; end
    D = R.(vv{v}).D; if strcmp(vv{v},'iso'), D = max(D,1e-19); end
    lw = 1.3; ls = '-';
    if v == 2, ls = '--'; lw = 1.4; end
    if v == 3, lw = 1.6; end
    semilogy(axR, R.(vv{v}).t, max(D,1e-19), ls, 'Color', cols(v,:), 'LineWidth', lw); hold on;
end
ylim([1e-20 10]); if nargin > 4, xlim([0 tmaxD]); end
grid on; xlabel('$t$'); ylabel('$D(u_n)$');
legend(ll,'Location','southeast','FontSize',7.2);
exportgraphics(f, fullfile(here,'..','..','figs',outname),'ContentType','vector');
end
