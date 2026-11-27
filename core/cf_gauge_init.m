function [c_q, U0] = cf_gauge_init(d, N, f_dbl, f_vpa, dig, mode)
%CF_GAUGE_INIT  Canonical-Frame Isotypic Spectral Method (CF-ISM) t=0 setup.
% Finds the canonical torus origin c in [-pi/2, pi/2)^d minimizing the
% alias-free oddness functional via 1 FFT auto-convolution + Gauss-Newton
% polish in `dig`-digit precision (at t=0 ONLY), then evaluates the initial
% isotypic channels U0 in `dig`-digit precision and rounds to IEEE double.
% ZERO filtering, ZERO mode zeroing, ZERO knowledge of whether u0 is symmetric!
% mode: 'refl' (coordinate-wise odd reflection, default) | 'point' (C2 point odd)
if nargin < 5 || isempty(dig), dig = 320; end
if nargin < 6 || isempty(mode), mode = 'refl'; end
M = N/2; h = 2*pi/N; x1 = (0:N-1)'*h;
if d == 1
    u = f_dbl(x1); uh = fft(u);
    [~, idx] = min(real(ifft(uh.^2))/N);
    c = mod((idx - 1)*h/2 + pi/2, pi) - pi/2;
    Mc = min(32, M); xc = (1:Mc-1)'*pi/Mc;
    % Double Gauss-Newton on mean-centred alias-free oddness residual
    eps_fd = 1e-7;
    for it = 1:8
        r = f_dbl(c + xc) + f_dbl(c - xc); r = r - mean(r);
        dr = (f_dbl(c+eps_fd+xc)+f_dbl(c+eps_fd-xc)-f_dbl(c-eps_fd+xc)-f_dbl(c-eps_fd-xc))/(2*eps_fd);
        dr = dr - mean(dr);
        den = dr'*dr; if den < 1e-30, break; end
        dc = -(dr'*r)/den; c = c + dc;
        if abs(dc) < 1e-14, break; end
    end
    c = mod(c + pi/2, pi) - pi/2;
    % High-precision Gauss-Newton polish on Mc points (t=0 only)
    digits(dig);
    c_q = vpa(c); xc_q = (1:Mc-1)'*vpa(pi)/Mc;
    eps_q = vpa(10)^(-floor(dig/2));
    for it = 1:8
        r = f_vpa(c_q + xc_q) + f_vpa(c_q - xc_q); r = r - sum(r)/numel(r);
        dr = (f_vpa(c_q+eps_q+xc_q)+f_vpa(c_q+eps_q-xc_q)-f_vpa(c_q-eps_q+xc_q)-f_vpa(c_q-eps_q-xc_q))/(2*eps_q);
        dr = dr - sum(dr)/numel(dr);
        den = sum(dr.^2); if double(den) < 1e-30, break; end
        dc = -sum(dr.*r)/den; c_q = c_q + dc;
        if double(abs(dc)) < 10^(-(dig-6)), break; end
    end
    xh_q = (0:M)'*(2*vpa(pi)/N);
    vp = f_vpa(c_q + xh_q); vm = f_vpa(c_q - xh_q);
    ev = double((vp + vm)/2); od = double((vp - vm)/2);
    od([1 end]) = 0;
    U0 = [ev, od];
elseif d == 2
    [X,Y] = ndgrid(x1,x1); U = f_dbl(X,Y);
    if strcmp(mode,'refl')
        Ax = sum(real(ifft(fft(U,[],1).^2,[],1)), 2); [~,ix] = min(Ax);
        Ay = sum(real(ifft(fft(U,[],2).^2,[],2)), 1); [~,iy] = min(Ay);
        cx = mod((ix-1)*h/2 + pi/2, pi) - pi/2;
        cy = mod((iy-1)*h/2 + pi/2, pi) - pi/2;
        Mc = 10; xc = (1:Mc-1)'*pi/Mc; [Xc,Yc] = ndgrid(xc,xc); Xc=Xc(:); Yc=Yc(:);
        eps_fd = 1e-7;
        for it = 1:6
            rx = f_dbl(cx+Xc,cy+Yc)+f_dbl(cx-Xc,cy+Yc); rx = rx - mean(rx);
            ry = f_dbl(cx+Xc,cy+Yc)+f_dbl(cx+Xc,cy-Yc); ry = ry - mean(ry);
            drx = (f_dbl(cx+eps_fd+Xc,cy+Yc)+f_dbl(cx+eps_fd-Xc,cy+Yc)-f_dbl(cx-eps_fd+Xc,cy+Yc)-f_dbl(cx-eps_fd-Xc,cy+Yc))/(2*eps_fd);
            dry = (f_dbl(cx+Xc,cy+eps_fd+Yc)+f_dbl(cx+Xc,cy+eps_fd-Yc)-f_dbl(cx+Xc,cy-eps_fd+Yc)-f_dbl(cx+Xc,cy-eps_fd-Yc))/(2*eps_fd);
            drx = drx - mean(drx); dry = dry - mean(dry);
            cx = cx - (drx'*rx)/(drx'*drx + 1e-30);
            cy = cy - (dry'*ry)/(dry'*dry + 1e-30);
        end
        digits(dig); c_q = [vpa(cx); vpa(cy)];
        xc_q = (1:Mc-1)'*vpa(pi)/Mc; [Xq,Yq] = ndgrid(xc_q,xc_q); Xq=Xq(:); Yq=Yq(:);
        eps_q = vpa(10)^(-floor(dig/2));
        for it = 1:8
            rx = f_vpa(c_q(1)+Xq,c_q(2)+Yq)+f_vpa(c_q(1)-Xq,c_q(2)+Yq); rx = rx - sum(rx)/numel(rx);
            ry = f_vpa(c_q(1)+Xq,c_q(2)+Yq)+f_vpa(c_q(1)+Xq,c_q(2)-Yq); ry = ry - sum(ry)/numel(ry);
            drx = (f_vpa(c_q(1)+eps_q+Xq,c_q(2)+Yq)+f_vpa(c_q(1)+eps_q-Xq,c_q(2)+Yq)-f_vpa(c_q(1)-eps_q+Xq,c_q(2)+Yq)-f_vpa(c_q(1)-eps_q-Xq,c_q(2)+Yq))/(2*eps_q);
            dry = (f_vpa(c_q(1)+Xq,c_q(2)+eps_q+Yq)+f_vpa(c_q(1)+Xq,c_q(2)+eps_q-Yq)-f_vpa(c_q(1)+Xq,c_q(2)-eps_q+Yq)-f_vpa(c_q(1)+Xq,c_q(2)-eps_q-Yq))/(2*eps_q);
            drx = drx - sum(drx)/numel(drx); dry = dry - sum(dry)/numel(dry);
            dc1 = -sum(drx.*rx)/sum(drx.^2); dc2 = -sum(dry.*ry)/sum(dry.^2);
            c_q = c_q + [dc1; dc2];
            if double(max(abs([dc1 dc2]))) < 10^(-(dig-6)), break; end
        end
    else
        A = real(ifftn(fftn(U).^2)); [~,idx] = min(A(:)); [ix,iy] = ind2sub([N N], idx);
        cx = mod((ix-1)*h/2 + pi/2, pi) - pi/2; cy = mod((iy-1)*h/2 + pi/2, pi) - pi/2;
        digits(dig); c_q = [vpa(cx); vpa(cy)];
        Mc = 8; xc_q = (1:Mc-1)'*vpa(pi)/Mc; [Xq,Yq] = ndgrid(xc_q, [xc_q; -xc_q]); Xq=Xq(:); Yq=Yq(:);
        eps_q = vpa(10)^(-floor(dig/2));
        for it = 1:9
            r = f_vpa(c_q(1)+Xq,c_q(2)+Yq)+f_vpa(c_q(1)-Xq,c_q(2)-Yq); r = r - sum(r)/numel(r);
            J1 = (f_vpa(c_q(1)+eps_q+Xq,c_q(2)+Yq)+f_vpa(c_q(1)+eps_q-Xq,c_q(2)-Yq)-f_vpa(c_q(1)-eps_q+Xq,c_q(2)+Yq)-f_vpa(c_q(1)-eps_q-Xq,c_q(2)-Yq))/(2*eps_q);
            J2 = (f_vpa(c_q(1)+Xq,c_q(2)+eps_q+Yq)+f_vpa(c_q(1)-Xq,c_q(2)+eps_q-Yq)-f_vpa(c_q(1)+Xq,c_q(2)-eps_q+Yq)-f_vpa(c_q(1)-Xq,c_q(2)-eps_q-Yq))/(2*eps_q);
            J1 = J1 - sum(J1)/numel(J1); J2 = J2 - sum(J2)/numel(J2);
            J = [J1, J2]; dc = -(J.'*J)\(J.'*r); c_q = c_q + dc;
            if double(max(abs(dc))) < 10^(-(dig-6)), break; end
        end
    end
    xh = (0:M)'*(2*vpa(pi)/N); [X,Y] = ndgrid(xh,xh);
    vpp = f_vpa(c_q(1)+X, c_q(2)+Y); vmp = f_vpa(c_q(1)-X, c_q(2)+Y);
    vpm = f_vpa(c_q(1)+X, c_q(2)-Y); vmm = f_vpa(c_q(1)-X, c_q(2)-Y);
    ep = (vpp+vmp)/2; op = (vpp-vmp)/2; em = (vpm+vmm)/2; om = (vpm-vmm)/2;
    U0 = zeros(M+1,M+1,4);
    U0(:,:,1) = double((ep+em)/2);
    U0(:,:,2) = double((op+om)/2); U0([1 end],:,2) = 0;
    U0(:,:,3) = double((ep-em)/2); U0(:,[1 end],3) = 0;
    U0(:,:,4) = double((op-om)/2); U0([1 end],:,4) = 0; U0(:,[1 end],4) = 0;
end
end
