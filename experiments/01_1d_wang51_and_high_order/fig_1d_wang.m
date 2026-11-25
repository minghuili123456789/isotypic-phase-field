% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Figure 1D: Wang's benchmark kappa=1, u0 = sin x
here = fileparts(mfilename('fullpath')); addpath(here); fig_common;
L = load(fullfile(here,'..','..','results','exp_law.mat')); E = load(fullfile(here,'..','..','results','explore1d.mat'));
f = figure('Units','centimeters','Position',[2 2 17 11.5]);
% (a) defect histories, S=0, double, various schemes
subplot(2,2,1); cols = lines(7); k = 0; lg = {};
for j = 1:numel(L.out)
    J = L.out{j};
    if strcmp(J.type,'fourier') && strcmp(J.prec,'double') && (J.tau==0.01 || strcmp(J.sch,'esp'))
        if strcmp(J.sch,'esp') && J.tau==5, continue; end
        k = k+1; semilogy(J.t, max(J.D,1e-20), 'Color', cols(k,:)); hold on;
        if strcmp(J.sch,'esp'), lg{end+1} = sprintf('ESP, $T_e=%g$',J.tau); %#ok<SAGROW>
        else, lg{end+1} = upper(J.sch); end %#ok<SAGROW>
    end
end
tt = 0:60; semilogy(tt, 1e-19*exp(tt), 'k--'); lg{end+1} = '$\propto e^{t}$';
semilogy([0 60],[10^-19.5 10^-19.5],'k-','LineWidth',2); lg{end+1} = 'isotypic: $D\equiv0$';
xlim([0 60]); ylim([1e-20 10]); xlabel('$t$'); ylabel('$D(u_n)$'); grid on;
legend(lg,'Location','southeast','NumColumns',2,'FontSize',6.0); title('(a) $S=0$: standard schemes, double');
% (b) breakdown-time law
subplot(2,2,2);
dl = []; ts = [];
for j = 1:numel(L.out)
    J = L.out{j}; if strcmp(J.type,'iso'), dl(end+1) = J.del; ts(end+1) = J.Tstar; end %#ok<SAGROW>
end
ok = dl < 1e-3;
pp = polyfit(-log(dl(ok)), ts(ok), 1);
semilogx(dl, ts, 'ko', 'MarkerFaceColor', [.3 .3 .9]); hold on;
xx = logspace(-42,-2,50); semilogx(xx, polyval(pp,-log(xx)), 'b-');
Tf = []; for j = 1:numel(L.out), J = L.out{j};
    if strcmp(J.type,'fourier') && strcmp(J.prec,'double'), Tf(end+1) = J.Tstar; end, end %#ok<SAGROW>
fill([1e-42 1e-2 1e-2 1e-42],[min(Tf) min(Tf) max(Tf) max(Tf)],[1 .6 .6],'EdgeColor','none','FaceAlpha',.5);
set(gca,'XDir','reverse'); xlim([1e-42 1e-2]); grid on;
xlabel('asymmetry $\delta$ of the data'); ylabel('$T_\ast$');
legend({'isotypic, $T_\ast(\delta)$', sprintf('fit: slope $1/\\lambda_\\perp$, $\\lambda_\\perp=%.3f$',1/pp(1)), ...
    'standard (all schemes, $\tau$), $\delta=0$'},'Location','northwest','FontSize',7);
title('(b) breakdown time law');
% (c) S=1: erratic behaviour of standard schemes (to t=3000/6000)
subplot(2,2,3); k = 0; lg = {};
for q = 1:numel(E.out)
    nm = E.out{q}.cfg.name;
    if startsWith(nm,'k1_') && (contains(nm,'_S1') || contains(nm,'esp'))
        F = E.out{q}.fourier; k = k+1;
        semilogy(F.t, max(F.D,1e-20), 'Color', cols(k,:)); hold on;
        lbl = strrep(strrep(nm,'k1_',''),'_S1','');
        if strcmp(lbl,'esp_Te1'), lbl = 'ESP ($T_e=1$)';
        elseif strcmp(lbl,'esp_Te20'), lbl = 'ESP ($T_e=20$)';
        else, lbl = upper(lbl); end
        lg{end+1} = lbl; %#ok<SAGROW>
    end
end
xlim([0 3000]); ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('$D(u_n)$');
legend(lg,'Location','southeast','FontSize',7); title('(c) $S=1$: standard schemes (erratic)');
% (d) final states, ESP Te=20
subplot(2,2,4);
q = find(cellfun(@(o) strcmp(o.cfg.name,'k1_esp_Te20'), E.out));
xg = (0:127)'*2*pi/128; xg(xg>=pi) = xg(xg>=pi)-2*pi; [xs,ix] = sort(xg);
plot(xs, E.out{q}.fourier.ufinal(ix), 'r-', xs, E.out{q}.iso.ufinal(ix), 'b--');
xlim([-pi pi]); ylim([-1.1 1.1]); grid on; xlabel('$x$'); ylabel('$u(x,T)$');
legend({'ESP standard','ESP isotypic'},'Location','southeast');
title(sprintf('(d) ESP ($T_e=20$, $S=1$), $T=%g$', E.out{q}.cfg.T));
exportgraphics(f, fullfile(here,'..','..','figs','fig_1d_wang.pdf'),'ContentType','vector');
fprintf('fit lambda_perp = %.4f, intercept %.3f\n', 1/pp(1), pp(2));
