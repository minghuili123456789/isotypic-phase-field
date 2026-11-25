% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Quickly compute the mid-collapse snapshot at t = 47 (T = 47 only, ~1 sec)
% and append S47_f, S47_i into results/exp_highorder_evol.mat, then regenerate figure.
here = fileparts(mfilename('fullpath')); addpath(here);
resfile = fullfile(here,'..','..','results','exp_highorder_evol.mat');
S = load(resfile);

N = 128; T47 = 47;
Pf = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','fourier')); Of = sp_ops(Pf);
Pi = sp_setup(struct('d',1,'N',N,'eps2',1,'gam',1,'S',0,'rep','iso'));     Oi = sp_ops(Pi);
x_raw = (0:N-1)'*2*pi/N; x = x_raw; x(x>=pi) = x(x>=pi) - 2*pi; [~, ix] = sort(x);
U0 = sp_init(Pi, @(xx) sin(xx)); u0 = sp_assemble(U0, 1, N/2);

methods = {'BDF2-SAV','ARS443','SBDF3','SBDF4','IFRK4','SDC4'};
taus    = [0.01,      0.02,    0.01,   0.01,   0.02,   0.05];
S47_f = cell(6,1);
for m = 1:6
    u = run47(methods{m}, Pf, Of, u0, taus(m), round(T47/taus(m)));
    S47_f{m} = u(ix);
    fprintf('t=47 %-10s: mean=%+.3f min=%.3f max=%.3f\n', methods{m}, mean(u), min(u), max(u));
end
ui = run47('SBDF4', Pi, Oi, U0, 0.01, round(T47/0.01));
uif = sp_assemble(ui, 1, N/2); S47_i = uif(ix);

save(resfile, 'S47_f', 'S47_i', '-append');
fig_highorder_1d;

function u = run47(nm, P, O, u0, tau, nst)
Lk = P.eps2 * P.K2;
switch nm
    case 'IFRK4'
        E2 = exp(-0.5*tau*Lk); E = exp(-tau*Lk); uh = O.toSpec(u0);
        for n = 1:nst
            k1=O.N(uh); u2=E2.*(uh+0.5*tau*k1); k2=O.N(u2);
            u3=E2.*uh+0.5*tau*k2; k3=O.N(u3); u4=E.*uh+tau*(E2.*k3); k4=O.N(u4);
            uh = E.*uh + (tau/6)*(E.*k1 + 2*E2.*(k2+k3) + k4);
        end
        u = O.toPhys(uh);
    case 'SBDF3'
        I4 = sp_integrator(P,O,'etdrk4',tau);
        uh0=O.toSpec(u0); N0=O.N(uh0); uh1=I4.step(uh0); N1=O.N(uh1); uh2=I4.step(uh1); N2=O.N(uh2);
        Dinv = 1./(11/6+tau*Lk);
        for n=3:nst, uh=Dinv.*(3*uh2-1.5*uh1+(1/3)*uh0+tau*(3*N2-3*N1+N0)); uh0=uh1;N0=N1;uh1=uh2;N1=N2;uh2=uh;N2=O.N(uh); end
        u = O.toPhys(uh2);
    case 'SBDF4'
        I4 = sp_integrator(P,O,'etdrk4',tau);
        uh0=O.toSpec(u0); N0=O.N(uh0); uh1=I4.step(uh0); N1=O.N(uh1); uh2=I4.step(uh1); N2=O.N(uh2); uh3=I4.step(uh2); N3=O.N(uh3);
        Dinv = 1./(25/12+tau*Lk);
        for n=4:nst, uh=Dinv.*(4*uh3-3*uh2+(4/3)*uh1-0.25*uh0+tau*(4*N3-6*N2+4*N1-N0)); uh0=uh1;N0=N1;uh1=uh2;N1=N2;uh2=uh3;N2=N3;uh3=uh;N3=O.N(uh); end
        u = O.toPhys(uh3);
    case 'ARS443'
        a=zeros(5,5); ah=zeros(5,5);
        a(2,2)=1/2; a(3,2)=1/6; a(3,3)=1/2; a(4,2)=-1/2; a(4,3)=1/2; a(4,4)=1/2; a(5,2)=3/2; a(5,3)=-3/2; a(5,4)=1/2; a(5,5)=1/2;
        ah(2,1)=1/2; ah(3,1)=11/18; ah(3,2)=1/18; ah(4,1)=5/6; ah(4,2)=-5/6; ah(4,3)=1/2; ah(5,1)=1/4; ah(5,2)=7/4; ah(5,3)=3/4; ah(5,4)=-7/4;
        Dinv=1./(1+0.5*tau*Lk); uh=O.toSpec(u0);
        for n=1:nst
            L_st=cell(5,1); N_st=cell(5,1); N_st{1}=O.N(uh); L_st{1}=-Lk.*uh;
            for i=2:5
                rhs=uh; for j=1:i-1, if a(i,j)~=0, rhs=rhs+(tau*a(i,j))*L_st{j}; end; if ah(i,j)~=0, rhs=rhs+(tau*ah(i,j))*N_st{j}; end; end
                Ui=Dinv.*rhs; L_st{i}=-Lk.*Ui; if i<5, N_st{i}=O.N(Ui); end
            end
            uh=Ui;
        end
        u = O.toPhys(uh);
    case 'BDF2-SAV'
        C0=10; h=2*pi/P.N; E1=(P.gam/4)*h*sum((u0.^2-1).^2); r0=sqrt(E1+C0);
        uh0=O.toSpec(u0); A1=1./(1+tau*Lk); b0=P.gam*(u0.^3-u0)/sqrt(E1+C0);
        g1h=A1.*uh0; g2h=A1.*O.toSpec(b0); ip=@(p,q) h*sum(p.*q);
        r1=(r0-0.5*ip(b0,u0)+0.5*ip(b0,O.toPhys(g1h)))/(1+0.5*tau*ip(b0,O.toPhys(g2h)));
        uh1=g1h-tau*r1*g2h; u1=O.toPhys(uh1); A2=1./(1.5+tau*Lk);
        for n=2:nst
            us=2*u1-u0; E1=(P.gam/4)*h*sum((us.^2-1).^2); b=P.gam*(us.^3-us)/sqrt(E1+C0); bh=O.toSpec(b);
            g1h=A2.*(2*uh1-0.5*uh0); g2h=A2.*bh;
            r2=((4*r1-r0)/3-(1/6)*ip(b,4*u1-u0)+0.5*ip(b,O.toPhys(g1h)))/(1+0.5*tau*ip(b,O.toPhys(g2h)));
            uh2=g1h-tau*r2*g2h; u2=O.toPhys(uh2); u0=u1;uh0=uh1;r0=r1; u1=u2;uh1=uh2;r1=r2;
        end
        u = u1;
    case 'SDC4'
        Sq=[5/24,8/24,-1/24; -1/24,8/24,5/24]; dt=tau/2; Dinv=1./(1+dt*Lk); uh=O.toSpec(u0);
        for n=1:nst
            U=cell(3,1); L=cell(3,1); Nv=cell(3,1); U{1}=uh; L{1}=-Lk.*uh; Nv{1}=O.N(uh);
            for m=1:2, U{m+1}=Dinv.*(U{m}+dt*Nv{m}); L{m+1}=-Lk.*U{m+1}; Nv{m+1}=O.N(U{m+1}); end
            for sw=1:5
                Un=cell(3,1); Ln=cell(3,1); Nn=cell(3,1); Un{1}=uh; Ln{1}=L{1}; Nn{1}=Nv{1};
                for m=1:2
                    Iint=tau*(Sq(m,1)*(L{1}+Nv{1})+Sq(m,2)*(L{2}+Nv{2})+Sq(m,3)*(L{3}+Nv{3}));
                    Un{m+1}=Dinv.*(Un{m}+dt*(Nn{m}-Nv{m})-dt*L{m+1}+Iint);
                    Ln{m+1}=-Lk.*Un{m+1}; Nn{m+1}=O.N(Un{m+1});
                end
                U=Un; L=Ln; Nv=Nn;
            end
            uh=U{3};
        end
        u = O.toPhys(uh);
end
end
