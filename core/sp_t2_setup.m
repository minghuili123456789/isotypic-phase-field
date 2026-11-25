function P = sp_t2_setup(opt)
%SP_T2_SETUP  Half-period translation (T_2 = <T_{pi e_i}> = Z_2^d) isotypic
% spectral discretisation of u_t = eps2*Lap(u) + gam*(u - u^3) on [0,L)^d.
% Channels are indexed by wavenumber parity pi in {0,1}^d (0: even k_i,
% 1: odd k_i).  Each channel lives on the half-domain [0,L/2)^d of size M^d
% (M = N/2), so the total number of grid points is 2^d * M^d = N^d.
% Because T_{pi e_i} commutes with ALL continuous translations T_a, any
% anti-periodic / odd-harmonic symmetry u(x + pi e_i) = -u(x) is preserved
% BITWISE in pure IEEE double precision for ANY arbitrary (even irrational)
% spatial shift x_0 in R^d, without shift detection or filtering.
if ~isfield(opt,'L'),     opt.L = 2*pi;        end
if ~isfield(opt,'prec'),  opt.prec = 'double'; end
if ~isfield(opt,'S'),     opt.S = 0;           end
if ~isfield(opt,'model'), opt.model = 'AC';    end
P = opt; d = opt.d; N = opt.N; M = N/2; L = opt.L;
assert(mod(N,4)==0, 'N must be divisible by 4 for T2 split');
P.M = M; P.h = L/N; P.nc = 2^d; P.rep = 't2';
cls = opt.prec;
% Wavenumbers in channel parity 0 (even) and 1 (odd):
% On [0,L/2) with M points, exp(i*2*pi*m*j/M) has physical wavenumber:
%   even channel: 2*(0:M/2, -M/2+1:-1) * (2*pi/L)
%   odd  channel: (2*(0:M/2-1, -M/2:-1) + 1) * (2*pi/L)
ke = (2*pi/L) * (2*[0:M/2, -M/2+1:-1]');
ko = (2*pi/L) * (2*[0:M/2-1, -M/2:-1]' + 1);
K2 = zeros([repmat(M,1,d) P.nc], cls);
P.par = zeros(P.nc, d);
for c = 0:P.nc-1
    ks = cell(1,d);
    for i = 1:d
        P.par(c+1,i) = bitget(c,i);
        if bitget(c,i), ks{i} = ko; else, ks{i} = ke; end
    end
    switch d
        case 1, K2(:,c+1) = cast(ks{1}.^2, cls);
        case 2, [a,b] = ndgrid(ks{:}); K2(:,:,c+1) = cast(a.^2 + b.^2, cls);
        case 3, [a,b,e] = ndgrid(ks{:}); K2(:,:,:,c+1) = cast(a.^2 + b.^2 + e.^2, cls);
    end
end
P.K2 = K2;
P.Lk = opt.eps2*K2 + opt.S;
P.tw = cast(exp(-1i*pi*(0:M-1)'/M), cls);   % anti-periodic twiddle factor
end
