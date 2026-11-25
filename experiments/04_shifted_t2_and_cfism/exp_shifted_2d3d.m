% Auto-add core library path (repository root/core)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','core'));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','results'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','results')); end
if ~exist(fullfile(fileparts(mfilename('fullpath')),'..','..','figs'),'dir'), mkdir(fullfile(fileparts(mfilename('fullpath')),'..','..','figs')); end
% Direction B 2D & 3D irrational-shift experiments + T2 cost benchmarks
addpath(fileparts(mfilename('fullpath')));
res = fullfile(fileparts(mfilename('fullpath')),'..','..','results');

th1 = 0.37123456; th2 = 0.71828182; th3 = 1.13141592;

fprintf('=== (1) 2D irrational-shift X-junction (N=256, eps=0.05, T=10) ===\n');
ep = 0.05; g = 1/ep^2; S = 2*g; tau = 1e-3; T2d = 10; N2 = 256;
f_x2d = @(x,y) tanh(sin(x + th1).*sin(y + th2)/(sqrt(2)*ep));
Pt2 = sp_t2_setup(struct('d',2,'N',N2,'eps2',1,'gam',g,'S',S));
Ot2 = sp_t2_ops(Pt2); It2 = sp_integrator(Pt2,Ot2,'etdrk2',tau);
[U0_t2, uf0] = sp_t2_init(Pt2, f_x2d, 3);   % channel 3 = odd in x AND odd in y
snaps2 = [0 0.5 2 10];
% Run T2 isotypic (pure double, zero VPA, zero shift detection!)
vh = Ot2.toSpec(U0_t2);
nst = round(T2d/tau); nout = 100;
t_t2 = (0:nout:nst)'*tau; D_t2 = zeros(size(t_t2)); snap_t2 = cell(size(snaps2));
km = 0; tic;
for n = 0:nst
    if mod(n,nout)==0
        km = km+1; V = Ot2.toPhys(vh);
        D_t2(km) = max(abs(reshape(V(:,:,1:3),[],1)));
    end
    s = find(abs(snaps2 - n*tau) < tau/2);
    if ~isempty(s), snap_t2{s} = sp_t2_assemble(Ot2.toPhys(vh), 2, Pt2.M); end
    if n < nst, vh = It2.step(vh); end
end
cpu_t2 = toc;
% Compare against standard Fourier ETDRK2 on EXACT SAME bitwise initial vector uf0!
Pf = sp_setup(struct('d',2,'N',N2,'eps2',1,'gam',g,'S',S,'rep','fourier'));
Of = sp_ops(Pf); If = sp_integrator(Pf,Of,'etdrk2',tau);
uh = Of.toSpec(uf0);
D_f = zeros(size(t_t2)); snap_f = cell(size(snaps2)); km = 0; tic;
for n = 0:nst
    if mod(n,nout)==0
        km = km+1; u = Of.toPhys(uh);
        Ut = t2_channels_local(u, 2, Pt2.M);
        D_f(km) = max(abs(reshape(Ut(:,:,1:3),[],1)));
    end
    s = find(abs(snaps2 - n*tau) < tau/2);
    if ~isempty(s), snap_f{s} = Of.toPhys(uh); end
    if n < nst, uh = If.step(uh); end
end
cpu_f = toc;
% Compare against fixed-origin reflection isotypic on EXACT SAME uf0 (must break!)
Pi = sp_setup(struct('d',2,'N',N2,'eps2',1,'gam',g,'S',S,'rep','iso'));
Oi = sp_ops(Pi); Ii = sp_integrator(Pi,Oi,'etdrk2',tau);
Ui0 = sp_channels(uf0, Pf);
vh = Oi.toSpec(Ui0); D_iso0 = zeros(size(t_t2)); km = 0;
for n = 0:nst
    if mod(n,nout)==0
        km = km+1; u = sp_assemble(Oi.toPhys(vh), 2, Pi.M);
        Ut = t2_channels_local(u, 2, Pt2.M);
        D_iso0(km) = max(abs(reshape(Ut(:,:,1:3),[],1)));
    end
    if n < nst, vh = Ii.step(vh); end
end
tb_f = t_t2(find(D_f>1e-2,1)); tb_i0 = t_t2(find(D_iso0>1e-2,1));
fprintf('  2D shifted X-junction: Std T*=%.2f (cpu=%.1fs) | Fixed-origin iso T*=%.2f | T2 iso maxD=%g (cpu=%.1fs)\n', ...
    tb_f, cpu_f, tb_i0, max(D_t2), cpu_t2);
X2d_sh = struct('t',t_t2,'D_std',D_f,'D_iso0',D_iso0,'D_t2',D_t2,'snapT',snaps2, ...
    'snap_std',{snap_f},'snap_t2',{snap_t2},'cpu_std',cpu_f,'cpu_t2',cpu_t2);

fprintf('\n=== (2) 2D shifted MIXED-harmonic X-junction via CF-ISM (N=128, eps=0.1, T=4) ===\n');
epm = 0.1; gm = 1/epm^2; Sm = 2*gm; taum = 1e-3; Tm = 4; Nm = 128;
fm_d = @(x,y) tanh((sin(x+th1)+0.4*sin(4*(x+th1))).*(sin(y+th2)+0.4*sin(4*(y+th2)))/(sqrt(2)*epm));
fm_v = @(x,y) tanh((sin(x+vpa(th1))+vpa(0.4)*sin(4*(x+vpa(th1)))).*(sin(y+vpa(th2))+vpa(0.4)*sin(4*(y+vpa(th2))))/(sqrt(vpa(2))*vpa(epm)));
tic; [c2_q, U0_cf2] = cf_gauge_init(2, Nm, fm_d, fm_v, 320, 'refl'); t_init2 = toc;
Pm = sp_setup(struct('d',2,'N',Nm,'eps2',1,'gam',gm,'S',Sm,'rep','iso')); Pm.xform = 'mat';
Om = sp_ops(Pm); Im = sp_integrator(Pm,Om,'etdrk2',taum);
vh = Om.toSpec(U0_cf2); maxD_cf2 = max(abs(reshape(U0_cf2(:,:,1:3),[],1)));
for n = 1:round(Tm/taum)
    vh = Im.step(vh);
    if mod(n,100)==0, V = Om.toPhys(vh); maxD_cf2 = max(maxD_cf2, max(abs(reshape(V(:,:,1:3),[],1)))); end
end
V_cf2_end = sp_assemble(Om.toPhys(vh), 2, Pm.M);
% Standard Fourier (N=128 is equivariant in 2D -> use N=192 or unaligned grid! Even on N=128 shifted data has non-zero defect!)
Pmf = sp_setup(struct('d',2,'N',Nm,'eps2',1,'gam',gm,'S',Sm,'rep','fourier'));
Omf = sp_ops(Pmf); Imf = sp_integrator(Pmf,Omf,'etdrk2',taum);
xm = (0:Nm-1)'*2*pi/Nm; [Xm,Ym] = ndgrid(xm,xm);
uh = Omf.toSpec(fm_d(Xm,Ym));
for n = 1:round(Tm/taum), uh = Imf.step(uh); end
u_mf_end = Omf.toPhys(uh);
fprintf('  2D shifted mixed-harm: t_init=%.2fs, init_forbid=%g, CF-ISM maxD=%g | Std mean=%.3e\n', ...
    t_init2, max(abs(reshape(U0_cf2(:,:,1:3),[],1))), maxD_cf2, mean(u_mf_end(:)));
M2d_sh = struct('t_init',t_init2,'maxD_cf',maxD_cf2,'u_cf_end',V_cf2_end,'u_std_end',u_mf_end,'c',double(c2_q));

fprintf('\n=== (3) 3D irrational-shift octant saddle (N=96, eps=0.1, T=5) ===\n');
ep3 = 0.1; g3 = 1/ep3^2; S3 = 2*g3; tau3 = 2e-3; T3 = 5; N3 = 96;
f_o3d = @(x,y,z) tanh(sin(x+th1).*sin(y+th2).*sin(z+th3)/(sqrt(2)*ep3));
Pt3 = sp_t2_setup(struct('d',3,'N',N3,'eps2',1,'gam',g3,'S',S3));
Ot3 = sp_t2_ops(Pt3); It3 = sp_integrator(Pt3,Ot3,'etdrk2',tau3);
[U0_3, uf3_0] = sp_t2_init(Pt3, f_o3d, 7);   % channel 7 = odd in x,y,z
snaps3 = [0 1 3 5];
nst3 = round(T3/tau3); nout3 = 50;
t3 = (0:nout3:nst3)'*tau3; D3_t2 = zeros(size(t3)); snap3_t2 = cell(size(snaps3));
vh = Ot3.toSpec(U0_3); km = 0; tic;
for n = 0:nst3
    if mod(n,nout3)==0
        km = km+1; V = Ot3.toPhys(vh);
        D3_t2(km) = max(abs(reshape(V(:,:,:,1:7),[],1)));
    end
    s = find(abs(snaps3 - n*tau3) < tau3/2);
    if ~isempty(s), snap3_t2{s} = sp_t2_assemble(Ot3.toPhys(vh), 3, Pt3.M); end
    if n < nst3, vh = It3.step(vh); end
end
cpu3_t2 = toc;
Pf3 = sp_setup(struct('d',3,'N',N3,'eps2',1,'gam',g3,'S',S3,'rep','fourier'));
Of3 = sp_ops(Pf3); If3 = sp_integrator(Pf3,Of3,'etdrk2',tau3);
uh = Of3.toSpec(uf3_0); D3_f = zeros(size(t3)); snap3_f = cell(size(snaps3)); km = 0; tic;
for n = 0:nst3
    if mod(n,nout3)==0
        km = km+1; u = Of3.toPhys(uh);
        Ut = t2_channels_local(u, 3, Pt3.M);
        D3_f(km) = max(abs(reshape(Ut(:,:,:,1:7),[],1)));
    end
    s = find(abs(snaps3 - n*tau3) < tau3/2);
    if ~isempty(s), snap3_f{s} = Of3.toPhys(uh); end
    if n < nst3, uh = If3.step(uh); end
end
cpu3_f = toc;
tb3_f = t3(find(D3_f>1e-2,1));
fprintf('  3D shifted octant (N=96): Std T*=%.2f (cpu=%.1fs) | T2 iso maxD=%g (cpu=%.1fs)\n', ...
    tb3_f, cpu3_f, max(D3_t2), cpu3_t2);
O3d_sh = struct('t',t3,'D_std',D3_f,'D_t2',D3_t2,'snapT',snaps3, ...
    'snap_std',{snap3_f},'snap_t2',{snap3_t2},'cpu_std',cpu3_f,'cpu_t2',cpu3_t2);

fprintf('\n=== (4) Wall-clock cost of T2 isotypic vs standard Fourier (ETDRK4, 4 threads) ===\n');
maxNumCompThreads(4);
cost_cases = {struct('d',1,'N',4096,'odd',1,'fun',@(x) sin(x+th1)), ...
              struct('d',2,'N',256, 'odd',3,'fun',@(x,y) sin(x+th1).*sin(y+th2)), ...
              struct('d',3,'N',64,  'odd',7,'fun',@(x,y,z) sin(x+th1).*sin(y+th2).*sin(z+th3)), ...
              struct('d',3,'N',128, 'odd',7,'fun',@(x,y,z) sin(x+th1).*sin(y+th2).*sin(z+th3))};
Tcost = table();
for q = 1:numel(cost_cases)
    c = cost_cases{q}; nst_c = 40; if c.d==3 && c.N==128, nst_c = 10; end
    Pt = sp_t2_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200));
    [Ut0, uf_c] = sp_t2_init(Pt, c.fun, c.odd);
    % Standard Fourier
    Pfc = sp_setup(struct('d',c.d,'N',c.N,'eps2',1,'gam',100,'S',200,'rep','fourier'));
    Ofc = sp_ops(Pfc); Ifc = sp_integrator(Pfc,Ofc,'etdrk4',1e-3);
    uh = Ofc.toSpec(uf_c); uh = Ifc.step(uh);
    tic; for n=1:nst_c, uh = Ifc.step(uh); end; t_std = toc/nst_c;
    % T2 full
    Otc = sp_t2_ops(Pt); Itc = sp_integrator(Pt,Otc,'etdrk4',1e-3);
    vh = Otc.toSpec(Ut0); vh = Itc.step(vh);
    tic; for n=1:nst_c, vh = Itc.step(vh); end; t_t2f = toc/nst_c;
    U_full = Otc.toPhys(vh);
    % T2 reduced
    Ptr = Pt; Ptr.active = false(1,Pt.nc); Ptr.active(c.odd+1) = true;
    Otr = sp_t2_ops(Ptr); Itr = sp_integrator(Ptr,Otr,'etdrk4',1e-3);
    vh = Otr.toSpec(Ut0); vh = Itr.step(vh);
    tic; for n=1:nst_c, vh = Itr.step(vh); end; t_t2r = toc/nst_c;
    U_red = Otr.toPhys(vh);
    same = isequal(abs(U_full), abs(U_red));
    Tcost = [Tcost; table(c.d, c.N, t_std, t_t2f, t_t2r, same, ...
        'VariableNames',{'d','N','fourier','t2_full','t2_reduced','bitwise_identical'})]; %#ok<AGROW>
    disp(Tcost(end,:));
end

save(fullfile(res,'exp_shifted_2d3d.mat'),'X2d_sh','M2d_sh','O3d_sh','Tcost','-v7.3');

function U = t2_channels_local(u, d, M)
N = 2*M; nc = 2^d; vals = cell(1,nc); vals{1} = u; sh = [M+1:N, 1:M];
for i = 1:d
    nv = cell(1,nc);
    for g = 0:2^(i-1)-1
        p = vals{g+1}; cc = repmat({':'},1,d); cc{i} = sh; q = p(cc{:});
        nv{g+1} = (p + q)/2; nv{g+1+2^(i-1)} = (p - q)/2;
    end
    vals(1:2^i) = nv(1:2^i);
end
cc = repmat({1:M},1,d); for c = 1:nc, vals{c} = vals{c}(cc{:}); end
U = cat(d+1, vals{:});
end
