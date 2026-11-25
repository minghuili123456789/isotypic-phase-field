% Standalone verification of the Chebfun spin / spin2 (Cox-Matthews ETDRK4)
% drop-in patch: compares standard fft/fft2 ETDRK4 vs. iso_spectral_backend
% on 1D Allen-Cahn (u0 = sin x, T = 200) and 2D X-junction (N = 256, T = 2).
clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'matlab_ode_and_spectral'));
outfig = fullfile(here, 'fig_demo_chebfun_spin.png');

% --- Part 1: 1D Chebfun spin-style ETDRK4 (u0 = sin x, N = 128, tau = 0.05, T = 200) ---
N1 = 128; tau1 = 0.05; T1 = 200; nst1 = round(T1/tau1);
B1 = iso_spectral_backend(1, N1, 'AR');
k1 = [0:N1/2-1, -N1/2:-1]'; K2_1f = k1.^2;
[Ef1,E2f1,Qf1,f1f1,f2f1,f3f1] = etd4_weights(-tau1*K2_1f, tau1);
[Ei1,E2i1,Qi1,f1i1,f2i1,f3i1] = etd4_weights(-tau1*B1.K2, tau1);

U0_1 = [zeros(B1.M+1,1), sin((0:B1.M)'*2*pi/N1)]; U0_1([1 end],2) = 0;
u0_1 = B1.assemble(U0_1);
uhf = fft(u0_1); uhi = B1.toSpec(U0_1);
t1 = (0:1:T1)'; D1f = zeros(numel(t1),1); D1i = D1f; r1 = round(1/tau1);
for n = 1:nst1
    % Standard Chebfun spin step
    Nf = @(v) fft(real(ifft(v)) - real(ifft(v)).^3);
    Nu=Nf(uhf); a=E2f1.*uhf+Qf1.*Nu; Na=Nf(a); b=E2f1.*uhf+Qf1.*Na; Nb=Nf(b);
    c=E2f1.*a+Qf1.*(2*Nb-Nu); Nc=Nf(c); uhf=Ef1.*uhf+f1f1.*Nu+2*f2f1.*(Na+Nb)+f3f1.*Nc;
    % Drop-in isotypic spin step (ONLY Nf changed!)
    Ni = @(V) B1.toSpec(B1.toPhys(V) - B1.cube(B1.toPhys(V)));
    Nu=Ni(uhi); a=E2i1.*uhi+Qi1.*Nu; Na=Ni(a); b=E2i1.*uhi+Qi1.*Na; Nb=Ni(b);
    c=E2i1.*a+Qi1.*(2*Nb-Nu); Nc=Ni(c); uhi=Ei1.*uhi+f1i1.*Nu+2*f2i1.*(Na+Nb)+f3i1.*Nc;
    if mod(n,r1)==0
        idx = round(n/r1)+1;
        Uc = B1.split(real(ifft(uhf))); D1f(idx) = max(abs(Uc(:,1)));
        Vi = B1.toPhys(uhi);            D1i(idx) = max(abs(Vi(:,1)));
    end
end
fprintf('Chebfun spin 1D (T=%d): std T*=%g (maxD=%.2e) | iso maxD=%g\n', ...
    T1, t1(find(D1f>1e-2,1)), max(D1f), max(D1i));

% --- Part 2: 2D Chebfun spin2-style ETDRK4 (X-junction, N = 256, eps = 0.05, T = 1.0) ---
N2 = 256; eps2 = 0.05; gam2 = 1/eps2^2; S2 = 2*gam2; tau2 = 2e-3; T2 = 1.0; nst2 = round(T2/tau2);
B2 = iso_spectral_backend(2, N2, 'AR');
k2 = [0:N2/2-1, -N2/2:-1]'; K2_2f = k2.^2 + (k2.').^2;
[Ef2,E2f2,Qf2,f1f2,f2f2,f3f2] = etd4_weights(-tau2*(K2_2f + S2), tau2);
[Ei2,E2i2,Qi2,f1i2,f2i2,f3i2] = etd4_weights(-tau2*(B2.K2 + S2), tau2);

xh = (0:B2.M)'*2*pi/N2; [Xh, Yh] = ndgrid(xh, xh);
U0_2 = zeros(B2.M+1, B2.M+1, 4);
U0_2(:,:,4) = tanh(sin(Xh).*sin(Yh)/(sqrt(2)*eps2));
U0_2([1 end],:,4) = 0; U0_2(:,[1 end],4) = 0;
u0_2 = B2.assemble(U0_2);
uhf2 = fft2(u0_2); uhi2 = B2.toSpec(U0_2);
t2 = (0:0.04:T2)'; D2f = zeros(numel(t2),1); D2i = D2f; r2 = round(0.04/tau2);
Nf2 = @(v) fft2((gam2+S2)*real(ifft2(v)) - gam2*real(ifft2(v)).^3);
Ni2 = @(V) B2.toSpec((gam2+S2)*B2.toPhys(V) - gam2*B2.cube(B2.toPhys(V)));
for n = 1:nst2
    Nu=Nf2(uhf2); a=E2f2.*uhf2+Qf2.*Nu; Na=Nf2(a); b=E2f2.*uhf2+Qf2.*Na; Nb=Nf2(b);
    c=E2f2.*a+Qf2.*(2*Nb-Nu); Nc=Nf2(c); uhf2=Ef2.*uhf2+f1f2.*Nu+2*f2f2.*(Na+Nb)+f3f2.*Nc;
    Nu=Ni2(uhi2); a=E2i2.*uhi2+Qi2.*Nu; Na=Ni2(a); b=E2i2.*uhi2+Qi2.*Na; Nb=Ni2(b);
    c=E2i2.*a+Qi2.*(2*Nb-Nu); Nc=Ni2(c); uhi2=Ei2.*uhi2+f1i2.*Nu+2*f2i2.*(Na+Nb)+f3i2.*Nc;
    if mod(n,r2)==0
        idx = round(n/r2)+1;
        Uc = B2.split(real(ifft2(uhf2))); D2f(idx) = max(abs(Uc(:,:,1:3)),[],'all');
        Vi = B2.toPhys(uhi2);             D2i(idx) = max(abs(Vi(:,:,1:3)),[],'all');
    end
end
fprintf('Chebfun spin2 2D (N=%d, T=%g): std T*=%.2f (maxD=%.2e) | iso maxD=%g\n', ...
    N2, T2, t2(find(D2f>1e-2,1)), max(D2f), max(D2i));

% Plot
set(groot,'defaultAxesFontSize',9.5,'defaultLineLineWidth',1.3, ...
    'defaultAxesTickLabelInterpreter','latex','defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex','defaultAxesBox','on');
f = figure('Units','centimeters','Position',[2 2 18 5.6]);
subplot(1,3,1);
semilogy(t1, max(D1f,1e-19),'r-','LineWidth',1.5); hold on;
semilogy(t1, max(D1i,1e-19),'b-','LineWidth',2.0);
xlim([0 T1]); ylim([1e-20 10]); grid on; xlabel('$t$'); ylabel('$D(u_n)$');
title('(a) 1D \texttt{spin}\ ETDRK4');
legend({'standard FFT','isotypic ($D\equiv0$)'},'Location','southeast','FontSize',6.8);

x2 = (0:N2-1)*2*pi/N2; x2(x2>=pi)=x2(x2>=pi)-2*pi; [xs2,ix2]=sort(x2);
subplot(1,3,2);
uf_end = real(ifft2(uhf2));
imagesc(xs2, xs2, uf_end(ix2,ix2)'); axis xy equal tight; caxis([-1 1]); set(gca,'XTick',[],'YTick',[]);
title('(b) 2D \texttt{spin2}\ standard ($t=1$, broken)');

subplot(1,3,3);
ui_end = B2.assemble(B2.toPhys(uhi2));
imagesc(xs2, xs2, ui_end(ix2,ix2)'); axis xy equal tight; caxis([-1 1]); set(gca,'XTick',[],'YTick',[]);
title('(c) 2D \texttt{spin2}\ isotypic ($t=1$, $D\equiv0$)');
colormap(parula);
exportgraphics(f, outfig, 'Resolution', 180);
fprintf('Saved %s\n', outfig);

function [E,E2,Q,f1,f2,f3] = etd4_weights(z, tau)
E = exp(z); E2 = exp(z/2);
Mq=64; r=exp(1i*pi*((1:Mq)-.5)/Mq); sz=size(z); zz=z(:);
% Q = (tau/2)*phi1(z/2)
z2 = zz/2; q = zeros(size(z2)); bq = abs(z2)>1;
q(bq) = (exp(z2(bq))-1)./z2(bq);
if any(~bq), Zc=z2(~bq)+r; q(~bq)=real(mean((exp(Zc)-1)./Zc,2)); end
Q = reshape((tau/2)*q, sz);
% f1, f2, f3
f1=zeros(size(zz)); f2=f1; f3=f1; big=abs(zz)>2; zb=zz(big);
f1(big)=(-4-zb+exp(zb).*(4-3*zb+zb.^2))./zb.^3;
f2(big)=(2+zb+exp(zb).*(-2+zb))./zb.^3;
f3(big)=(-4-3*zb-zb.^2+exp(zb).*(4-zb))./zb.^3;
if any(~big)
    Zc=zz(~big)+2*r;
    f1(~big)=real(mean((-4-Zc+exp(Zc).*(4-3*Zc+Zc.^2))./Zc.^3,2));
    f2(~big)=real(mean((2+Zc+exp(Zc).*(-2+Zc))./Zc.^3,2));
    f3(~big)=real(mean((-4-3*Zc-Zc.^2+exp(Zc).*(4-Zc))./Zc.^3,2));
end
f1=reshape(tau*f1,sz); f2=reshape(tau*f2,sz); f3=reshape(tau*f3,sz);
end
