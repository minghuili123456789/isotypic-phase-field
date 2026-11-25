function fig_3d(matfile, name, outname, tsel)
% zero-level isosurfaces (standard vs isotypic) and defect histories, 3D
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'..','..','core')); addpath(here); fig_common;
S = load(fullfile(here,'..','..','results',matfile));
R = struct;
for j = 1:numel(S.out)
    if strcmp(S.out{j}.name,name), R.(S.out{j}.variant) = S.out{j}; end
end
nt = numel(tsel);
f = figure('Units','centimeters','Position',[2 2 18 8.5],'Color','w');
vars = {'fourier','iso'}; lab = {'standard','isotypic'};
for v = 1:2
    for k = 1:nt
        s = find(abs(R.(vars{v}).snapT - tsel(k)) < 1e-9);
        U = R.(vars{v}).snap{s}; N = size(U,1);
        U = circshift(U, [N/2 N/2 N/2]);          % centre the box at the origin
        x = linspace(-pi, pi - 2*pi/N, N);
        [X,Y,Z] = meshgrid(x,x,x);
        ax = subplot(2,nt+1,(v-1)*(nt+1)+k);
        p = patch(isosurface(X,Y,Z,permute(U,[2 1 3]),0.5));
        isonormals(X,Y,Z,permute(U,[2 1 3]),p);
        p.FaceColor = [0.85 0.33 0.1]; p.EdgeColor = 'none';
        daspect([1 1 1]); view(35,25); axis([-pi pi -pi pi -pi pi]); camlight; lighting gouraud;
        set(ax,'XTick',[],'YTick',[],'ZTick',[]); box on;
        if v == 1, title(sprintf('$t=%g$',tsel(k)),'FontSize',9); end
        if k == 1, zlabel(lab{v}); end
    end
end
subplot(2,nt+1,[nt+1, 2*nt+2]);
cols = lines(3); vv = {'fourier','zmode','iso'}; ll = {'standard','zero-mode suppression','isotypic ($D\equiv0$)'};
for v = 1:3
    if ~isfield(R,vv{v}), continue; end
    semilogy(R.(vv{v}).t, max(R.(vv{v}).D,1e-19), 'Color', cols(v,:)); hold on;
end
ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('$D(u_n)$');
legend(ll,'Location','southeast','FontSize',7);
exportgraphics(f, fullfile(here,'..','..','figs',outname),'ContentType','image','Resolution',200);
end

