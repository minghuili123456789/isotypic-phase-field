% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure: Irrational-shift 2D X-junction & 3D octant saddle (standard vs T2 isotypic)
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results','exp_shifted_2d3d.mat'));
X = S.X2d_sh; O3 = S.O3d_sh;

% --- Fig A: 2D irrational-shift X-junction ---
f = figure('Units','centimeters','Position',[2 2 18 7.8]);
tsel = X.snapT; nt = numel(tsel);
N = size(X.snap_t2{1},1); x = (0:N-1)*2*pi/N; x(x>=pi) = x(x>=pi)-2*pi; [xs,ix] = sort(x);
for k = 1:nt
    ax = subplot(2, nt+1, k);
    imagesc(xs, xs, X.snap_std{k}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
    set(ax,'XTick',[],'YTick',[]); title(sprintf('$t=%g$',tsel(k)),'FontSize',9);
    if k == 1, ylabel('standard'); end
    ax = subplot(2, nt+1, nt+1+k);
    imagesc(xs, xs, X.snap_t2{k}(ix,ix)'); axis xy equal tight; caxis([-1 1]);
    set(ax,'XTick',[],'YTick',[]);
    if k == 1, ylabel('$\mathcal T_2$ isotypic'); end
end
colormap(parula);
subplot(2, nt+1, [nt+1, 2*nt+2]);
cols = lines(3);
semilogy(X.t, max(X.D_std,1e-19),  '-',  'Color', cols(1,:), 'LineWidth', 1.3); hold on;
semilogy(X.t, max(X.D_iso0,1e-19), '--', 'Color', cols(2,:), 'LineWidth', 1.3);
semilogy(X.t, max(X.D_t2,1e-19),   '-',  'Color', cols(3,:), 'LineWidth', 1.6);
xlim([0 2]); ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('$D(u_n)$');
legend({'standard FFT','fixed-origin iso','$\mathcal T_2$ iso ($D\equiv0$)'},'Location','southeast','FontSize',6.5);
exportgraphics(f, fullfile(here,'..','..','figs','fig_X2d_shifted.pdf'),'ContentType','vector');

% --- Fig B: 3D irrational-shift octant saddle ---
f3 = figure('Units','centimeters','Position',[2 2 18 8.2],'Color','w');
tsel3 = O3.snapT; nt3 = numel(tsel3);
for v = 1:2
    for k = 1:nt3
        if v == 1, U = O3.snap_std{k}; else, U = O3.snap_t2{k}; end
        N3 = size(U,1); U = circshift(U, [N3/2 N3/2 N3/2]);
        x3 = linspace(-pi, pi - 2*pi/N3, N3); [X3,Y3,Z3] = meshgrid(x3,x3,x3);
        ax = subplot(2, nt3+1, (v-1)*(nt3+1)+k);
        p = patch(isosurface(X3,Y3,Z3,permute(U,[2 1 3]),0.5));
        isonormals(X3,Y3,Z3,permute(U,[2 1 3]),p);
        p.FaceColor = [0.85 0.33 0.1]; p.EdgeColor = 'none';
        daspect([1 1 1]); view(35,25); axis([-pi pi -pi pi -pi pi]); camlight; lighting gouraud;
        set(ax,'XTick',[],'YTick',[],'ZTick',[]); box on;
        if v == 1, title(sprintf('$t=%g$',tsel3(k)),'FontSize',9); end
        if k == 1
            if v == 1, zlabel('standard'); else, zlabel('$\mathcal T_2$ isotypic'); end
        end
    end
end
subplot(2, nt3+1, [nt3+1, 2*nt3+2]);
semilogy(O3.t, max(O3.D_std,1e-19), '-', 'Color', cols(1,:), 'LineWidth', 1.3); hold on;
semilogy(O3.t, max(O3.D_t2,1e-19),  '-', 'Color', cols(3,:), 'LineWidth', 1.6);
ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('$D(u_n)$');
legend({'standard FFT','$\mathcal T_2$ iso ($D\equiv0$)'},'Location','southeast','FontSize',7);
exportgraphics(f3, fullfile(here,'..','..','figs','fig_O3d_shifted.png'),'ContentType','image','Resolution',200);
