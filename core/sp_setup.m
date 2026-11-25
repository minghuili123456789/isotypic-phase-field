function P = sp_setup(opt)
%SP_SETUP  Build a spectral discretisation of the (stabilised) Allen--Cahn flow
%      u_t = eps2*Lap(u) + gam*(u - u^3)
% on the periodic box [-L/2,L/2)^d (equivalently [0,L)^d), in one of two
% algebraically equivalent representations:
%   rep = 'fourier' : standard full-grid Fourier pseudo-spectral (FFT)
%   rep = 'iso'     : isotypic (symmetry-adapted) representation w.r.t. the
%                     reflection group Z_2^d, i.e. 2^d channels, each a
%                     tensor product of DCT-I / DST-I series on [0,L/2]^d.
% Both representations define exactly the same discrete operator in exact
% arithmetic; they differ only in how rounding errors are distributed.
%
% opt fields: d, N (even), L (default 2*pi), eps2, gam, S (stabiliser),
%             rep, prec ('double'|'single'), gpu (false), model ('AC'|'CAC')

if ~isfield(opt,'L'),     opt.L = 2*pi;        end
if ~isfield(opt,'prec'),  opt.prec = 'double'; end
if ~isfield(opt,'gpu'),   opt.gpu = false;     end
if ~isfield(opt,'S'),     opt.S = 0;           end
if ~isfield(opt,'model'), opt.model = 'AC';    end
P = opt;
d = opt.d; N = opt.N; M = N/2; L = opt.L;
P.M = M; P.h = L/N;
cls = opt.prec;

if strcmp(opt.rep,'fourier')
    kf = (2*pi/L)*[0:M, -M+1:-1]';
    switch d
        case 1, K2 = kf.^2;
        case 2, [a,b] = ndgrid(kf,kf); K2 = a.^2 + b.^2;
        case 3, [a,b,c] = ndgrid(kf,kf,kf); K2 = a.^2 + b.^2 + c.^2;
    end
    P.nc = 1;
    P.sig = zeros(1,d);
else
    kk = (2*pi/L)*(0:M)';
    switch d
        case 1, K2 = kk.^2;
        case 2, [a,b] = ndgrid(kk,kk); K2 = a.^2 + b.^2;
        case 3, [a,b,c] = ndgrid(kk,kk,kk); K2 = a.^2 + b.^2 + c.^2;
    end
    nc = 2^d; P.nc = nc;
    P.sig = zeros(nc,d);           % sig(c,i)=1 : channel c is odd in x_i
    for c = 0:nc-1
        for i = 1:d, P.sig(c+1,i) = bitget(c,i); end
    end
    K2 = repmat(K2,[ones(1,d) nc]);
end
P.K2 = cast(K2,cls);
if opt.gpu, P.K2 = gpuArray(P.K2); end
P.Lk = opt.eps2*P.K2 + opt.S;          % symbol of the stabilised linear part
end
