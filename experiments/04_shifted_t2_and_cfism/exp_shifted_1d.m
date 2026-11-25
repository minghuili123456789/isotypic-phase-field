% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Direction B 1D experiments:
% (A) CF-ISM breakdown horizon T*(d0) vs t=0 initialization digits d0
%     (pure IEEE double time stepping for ALL t > 0!) on mixed-harmonic
%     irrational-shift u0(x) = 0.5*(sin(x+th) + sin(4*(x+th))), th = 0.37123456789
% (B) Full Wang (2023) Example 5.2 (a)-(d) with shifts 0, pi/8, pi/4, 0.37123456789
%     comparing standard ESP (S=1), standard ETDRK4 (S=0), T2 isotypic, and CF-ISM.
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');

fprintf('=== (A) CF-ISM horizon law T*(d0) for irrational shift + mixed harmonics ===\n');
N = 128; tau = 0.05; th = 0.37123456789;
f_d = @(x) 0.5*(sin(x + th) + sin(4*(x + th)));
f_v = @(x) 0.5*(sin(x + vpa(th)) + sin(4*(x + vpa(th))));
d0_list = [16, 25, 40, 60, 100, 150, 220, 320];
P_iso = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));
O_iso = sp_ops(P_iso); I_iso = sp_integrator(P_iso,O_iso,'etdrk4',tau);
horizon = struct('d0',{},'init_ev',{},'Tstar',{},'maxD',{});
for q = 1:numel(d0_list)
    d0 = d0_list(q);
    [~, U0] = cf_gauge_init(1, N, f_d, f_v, d0, 'refl');
    Tmax = 800; nst = round(Tmax/tau);
    vh = O_iso.toSpec(U0);
    tb = NaN; maxD = max(abs(U0(:,1)));
    for n = 1:nst
        vh = I_iso.step(vh);
        if mod(n,2) == 0
            V = O_iso.toPhys(vh); D = max(abs(V(:,1)));
            maxD = max(maxD, D);
            if isnan(tb) && D > 1e-2, tb = n*tau; break; end
        end
    end
    horizon(q) = struct('d0',d0,'init_ev',max(abs(U0(:,1))),'Tstar',tb,'maxD',maxD);
    fprintf('  d0=%3d: init_ev=%9.2e | T*(D>1e-2)=%7.2f | maxD=%9.2e\n', d0, horizon(q).init_ev, tb, maxD);
end

fprintf('\n=== (B) Wang Ex 5.2 (a)-(d) with shifts 0, pi/8, pi/4, 0.37123456789 ===\n');
% (a) kappa=1, sin(x~), T=200; (b) kappa=0.9, sin(x~), T=200;
% (c) kappa=0.1, 0.5(sin(x~)+sin(4x~)), T=200; (d) kappa=sqrt(0.001), 0.5(sin(x~)+sin(8x~)), T=600
Nw = 512;
shs = [0, pi/8, pi/4, 0.37123456789];
cases = {
  struct('id','52a','eps2',1.0,  'T',200,'Te',1, 'tau',0.02,'odd_harm',true, ...
    'fd',@(x,s) sin(x+s), 'fv',@(x,s) sin(x+vpa(s))), ...
  struct('id','52b','eps2',0.81, 'T',200,'Te',2, 'tau',0.02,'odd_harm',true, ...
    'fd',@(x,s) sin(x+s), 'fv',@(x,s) sin(x+vpa(s))), ...
  struct('id','52c','eps2',0.01, 'T',200,'Te',2, 'tau',0.02,'odd_harm',false, ...
    'fd',@(x,s) 0.5*(sin(x+s)+sin(4*(x+s))), 'fv',@(x,s) 0.5*(sin(x+vpa(s))+sin(4*(x+vpa(s))))), ...
  struct('id','52d','eps2',0.001,'T',600,'Te',10,'tau',0.05,'odd_harm',false, ...
    'fd',@(x,s) 0.5*(sin(x+s)+sin(8*(x+s))), 'fv',@(x,s) 0.5*(sin(x+vpa(s))+sin(8*(x+vpa(s)))))
};
wang52 = cell(numel(cases), numel(shs));
xf = -pi + (0:Nw-1)'*2*pi/Nw;
for ic = 1:numel(cases)
    C = cases{ic};
    for is = 1:numel(shs)
        sh = shs(is);
        res_entry = struct('id',C.id,'sh',sh);
        % 1. Standard ESP (Wang's exact S=1, Ns=24, Niter=20)
        Pf = sp_setup(struct('d',1,'N',Nw,'eps2',C.eps2,'gam',1,'S',1,'rep','fourier'));
        Of = sp_ops(Pf); Iesp = sp_integrator(Pf,Of,'esp',C.Te,struct('Ns',24,'Niter',20));
        uh = Of.toSpec(C.fd(xf, sh));
        tb_esp = NaN; maxM_esp = 0;
        for n = 1:round(C.T/C.Te)
            uh = Iesp.step(uh);
            mu = abs(real(uh(1)))/Nw; maxM_esp = max(maxM_esp, mu);
            if isnan(tb_esp) && mu > 1e-2, tb_esp = n*C.Te; end
        end
        res_entry.esp_S1_Tstar = tb_esp; res_entry.esp_S1_maxMean = maxM_esp;
        % 2. Standard ETDRK4 (S=0)
        Pf0 = sp_setup(struct('d',1,'N',Nw,'eps2',C.eps2,'gam',1,'S',0,'rep','fourier'));
        Of0 = sp_ops(Pf0); I4 = sp_integrator(Pf0,Of0,'etdrk4',C.tau);
        uh = Of0.toSpec(C.fd(xf, sh));
        tb_f = NaN; maxM_f = 0;
        for n = 1:round(C.T/C.tau)
            uh = I4.step(uh);
            if mod(n,10)==0
                mu = abs(real(uh(1)))/Nw; maxM_f = max(maxM_f, mu);
                if isnan(tb_f) && mu > 1e-2, tb_f = n*C.tau; end
            end
        end
        res_entry.std_S0_Tstar = tb_f; res_entry.std_S0_maxMean = maxM_f;
        % 3. T2 isotypic (pure double, zero VPA) for odd-harmonic cases (a),(b)
        if C.odd_harm
            Pt2 = sp_t2_setup(struct('d',1,'N',Nw,'eps2',C.eps2,'gam',1,'S',0));
            Ot2 = sp_t2_ops(Pt2); It2 = sp_integrator(Pt2,Ot2,'etdrk4',C.tau);
            U0_t2 = sp_t2_init(Pt2, @(x) C.fd(x, sh), 1);
            vh = Ot2.toSpec(U0_t2); maxD_t2 = 0;
            for n = 1:round(C.T/C.tau)
                vh = It2.step(vh);
                if mod(n,50)==0, V = Ot2.toPhys(vh); maxD_t2 = max(maxD_t2, max(abs(V(:,1)))); end
            end
            res_entry.t2_maxD = maxD_t2;
        else
            res_entry.t2_maxD = NaN;
        end
        % 4. CF-ISM (d0=320 at t=0, pure double ETDRK4 S=0 for all t>0)
        [~, U0_cf] = cf_gauge_init(1, Nw, @(x) C.fd(x,sh), @(x) C.fv(x,sh), 320, 'refl');
        Pi = sp_setup(struct('d',1,'N',Nw,'eps2',C.eps2,'gam',1,'S',0,'rep','iso'));
        Oi = sp_ops(Pi); Ii = sp_integrator(Pi,Oi,'etdrk4',C.tau);
        vh = Oi.toSpec(U0_cf); maxD_cf = max(abs(U0_cf(:,1)));
        for n = 1:round(C.T/C.tau)
            vh = Ii.step(vh);
            if mod(n,50)==0, V = Oi.toPhys(vh); maxD_cf = max(maxD_cf, max(abs(V(:,1)))); end
        end
        res_entry.cf_maxD = maxD_cf;
        wang52{ic,is} = res_entry;
        fprintf('  %s sh=%7.4f | ESP(S=1) T*=%5g (maxMean=%.1e) | Std(S=0) T*=%6.1f | T2 maxD=%g | CF-ISM maxD=%g\n', ...
            C.id, sh, tb_esp, maxM_esp, tb_f, res_entry.t2_maxD, maxD_cf);
    end
end
save(fullfile(res,'exp_shifted_1d.mat'),'horizon','wang52','shs');
