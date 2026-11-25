% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Exact reproduction & long-time continuation of Wang (2023) Example 4.1 (2D Allen-Cahn):
%   du/dt = kappa^2 * Delta u + u - u^3,  kappa = 0.1,  Omega = [-pi, pi]^2
%   u0(x,y) = sin(2x)*sin(2y),  Np = 256 x 256
%   Wang's EXACT ESP parameters: S = 1, Te = 1, Ns = 10, Niter = 20
% Wang (2023) Fig. 5 stopped at T = 20! Here we run under 100% IDENTICAL conditions
% up to T = 120 (or until breakdown + far beyond) comparing:
%   1. Wang's standard FFT ESP (S=1, Te=1, Ns=10, Niter=20)
%   2. Standard FFT ESP + a posteriori zero-mode suppression
%   3. Proposed drop-in ISOTYPIC ESP (EXACT SAME S=1, Te=1, Ns=10, Niter=20!)
here = fileparts(mfilename('fullpath')); addpath(here);
res  = fullfile(here,'..','..','results');

N = 256; kap = 0.1; S_stab = 1; Te = 1; Ns = 10; Nit = 20; T = 120;
Pf = sp_setup(struct('d',2,'N',N,'eps2',kap^2,'gam',1,'S',S_stab,'rep','fourier'));
Pi = sp_setup(struct('d',2,'N',N,'eps2',kap^2,'gam',1,'S',S_stab,'rep','iso'));
Of = sp_ops(Pf); Oi = sp_ops(Pi);

If = sp_integrator(Pf, Of, 'esp', Te, struct('Ns',Ns,'Niter',Nit));
Ii = sp_integrator(Pi, Oi, 'esp', Te, struct('Ns',Ns,'Niter',Nit));

% Exact bitwise-symmetric initial datum u0(x,y) = sin(2x)*sin(2y)
U0 = sp_init(Pi, @(x,y) sin(2*x).*sin(2*y));
u0 = sp_assemble(U0, 2, N/2);
fprintf('Initial max transverse channel of U0 = %g\n', max(abs(U0(:,:,[1 2 3])),[],'all'));

% Also check energy E(u) = (kap^2/2)*int |grad u|^2 + (1/4)*int (u^2-1)^2
dx = 2*pi/N;
calcE = @(u, uh) 0.5*kap^2*(dx^2/N^2)*sum(Pf.K2(:).*abs(uh(:)).^2) + 0.25*dx^2*sum((u(:).^2 - 1).^2);

tsnaps = [1, 5, 20, 40, 60, 80, 100, 120];
tvec   = (0:T)';
D_std = zeros(T+1,1); D_zm = D_std; D_iso = D_std;
E_std = zeros(T+1,1); E_zm = E_std; E_iso = E_std;
snap_std = cell(numel(tsnaps),1);
snap_zm  = cell(numel(tsnaps),1);
snap_iso = cell(numel(tsnaps),1);

uh_f = Of.toSpec(u0);
uh_z = uh_f;
uh_i = Oi.toSpec(U0);

D_std(1) = 0; D_zm(1) = 0; D_iso(1) = 0;
E0 = calcE(u0, uh_f); E_std(1) = E0; E_zm(1) = E0; E_iso(1) = E0;

fprintf('t    | D_std      E_std   | D_zmode    E_zm    | D_iso      E_iso\n');
fprintf('--------------------------------------------------------------------\n');
for t = 1:T
    % 1. Wang's standard FFT ESP
    uh_f = If.step(uh_f);
    uf   = Of.toPhys(uh_f);
    Uc_f = sp_channels(uf, Pf);
    D_std(t+1) = max(abs(Uc_f(:,:,[1 2 3])),[],'all');
    E_std(t+1) = calcE(uf, uh_f);

    % 2. Standard FFT ESP + zero-mode suppression
    uh_z = If.step(uh_z); uh_z(1,1) = 0;
    uz   = Of.toPhys(uh_z);
    Uc_z = sp_channels(uz, Pf);
    D_zm(t+1) = max(abs(Uc_z(:,:,[1 2 3])),[],'all');
    E_zm(t+1) = calcE(uz, uh_z);

    % 3. Drop-in Isotypic ESP (EXACT same Te=1, Ns=10, Niter=20, S=1!)
    uh_i = Ii.step(uh_i);
    Vi   = Oi.toPhys(uh_i);
    D_iso(t+1) = max(abs(Vi(:,:,[1 2 3])),[],'all');
    ui   = sp_assemble(Vi, 2, N/2);
    E_iso(t+1) = calcE(ui, Of.toSpec(ui));

    if mod(t, 10) == 0 || t == 1 || t == 5 || t == 20 || (D_std(t+1) > 1e-3 && D_std(t) <= 1e-3)
        fprintf('%3d  | %9.2e  %7.4f | %9.2e  %7.4f | %9.2e  %7.4f\n', ...
            t, D_std(t+1), E_std(t+1), D_zm(t+1), E_zm(t+1), D_iso(t+1), E_iso(t+1));
    end
    k = find(tsnaps == t, 1);
    if ~isempty(k)
        snap_std{k} = uf; snap_zm{k} = uz; snap_iso{k} = ui;
    end
end

save(fullfile(res,'exp_wang41_2d.mat'), 'N','kap','S_stab','Te','Ns','Nit','T', ...
    'tvec','tsnaps','D_std','D_zm','D_iso','E_std','E_zm','E_iso', ...
    'snap_std','snap_zm','snap_iso');
fprintf('Saved exp_wang41_2d.mat\n');
