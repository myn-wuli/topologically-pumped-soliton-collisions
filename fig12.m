%% pengzhuang_workspace_peaks.m
% 默认只搜索并演化严格固定峰值的无节点基态；保留可选的激发态碰撞功能。
% 驱动势 p1=p2=25、扫描速率 nu=0.05、相互作用 g=-1 保持不变。
% 定态搜索收敛不等价于有限速率泵浦保持局域；不人为截去晶格尾部。
%
% 默认示例：第二能带 bandIndex = 2。
% 改 bandIndex 可更换 Wannier 初猜；sol.ampFactor 指定激发态的峰值振幅。
%
% 方程：

%   i psi_t = [-1/2 d_xx + V(x,t) + g |psi|^2] psi,
%   V(x,t) = -p1 cos^2(pi*x/d1) - p2 cos^2(pi*x/d2 - nu*t).
%
% 流程：
%   1. Bloch 平面波基底计算 E_n(kx,ky)；
%   2. 在起始晶格相位 ky=nu*t0 构造第 bandIndex 条能带 Wannier 态；
%   3. 固定峰值振幅 max(abs(psi))=sol.ampFactor，增广 Newton 同时求解 psi 和 mu；
%      化学势 mu 与粒子数 NB 均为输出；搜索后不再缩放振幅。
%   4. 以组合初态在滑动晶格中做 split-step 泵浦演化。
%   5. 计算完成后保存整个脚本工作区，文件名标明基态和激发态的实际初始峰值振幅。
%   6. 保存成功后，在主脚本中统一绘图。

clear; clc; close all;

%% ========================== 总开关 ==========================
job.bandIndex = 2;              % 第二能带为例
job.fastMode = true;            % true: 先用较小网格跑通；正式计算可改 false
job.calcChern = false;          % Chern 对泵浦演化不是必须；需要时改 true
job.doPump = true;             % false: 只搜索两个定态并作初态图
job.excitedCenter = 3;          % Wannier 晶胞平移 R；实际中心还含带内位置偏移
job.groundCenter  = -2;         % 基态初猜中心
job.groundWeight   = 1;         % 基态孤子权重；设为0时，只演化激发态，结果应与单独跑激发态一致
job.excitWeight   = 1;          % 默认只搜索并演化基态；设为1恢复碰撞
% 激发态给定振幅对应未缩放的 psi；实际输入的激发态振幅为 abs(excitWeight)*ampFactor。
% excitWeight=1 保持求得的定态，excitWeight=0 关闭激发态。
job.solitonCenter = job.excitedCenter;  % 兼容原有变量名：Wannier 晶胞平移

%% ========================== 晶格与非线性参数 ==========================
lat.p1 = 25;
lat.p2 = 25;
lat.d1 = 1/2;
lat.d2 = 1;
lat.a  = 1;
lat.g  = -1;                    % 吸引非线性
lat.nu = 0.05;                   % 泵浦速度；周期 T = pi/nu

%% ========================== 能带/Wannier 参数 ==========================
bpar.p1 = lat.p1;
bpar.p2 = lat.p2;
bpar.d1 = lat.d1;
bpar.d2 = lat.d2;
bpar.a  = lat.a;

if job.fastMode
    bpar.NG = 61;               % 调试/快速运行
    bpar.Nkx = 128;
    bpar.Nky = 81;
else
    bpar.NG = 81;               % 正式计算可提高到 81/101
    bpar.Nkx = 256;
    bpar.Nky = 121;
end
bpar.Nband = max(6, job.bandIndex);                 % 计算前几条能带
bpar.kx = (-bpar.Nkx/2 : bpar.Nkx/2-1).' * (2*pi/(bpar.a*bpar.Nkx));
bpar.ky = (0:bpar.Nky-1).' * (pi/bpar.Nky);

bpar.iyWannier = 1;             % 起始晶格相位处构造 Wannier 态
bpar.dx = 0.05;                 % 和非线性泵浦网格保持一致
bpar.NxCell = round(bpar.a/bpar.dx);
bpar.NcellWannier = bpar.Nkx;
bpar.R0 = job.solitonCenter;     % 晶胞平移，实际中心另行计算

%% ========================== 孤子求解参数 ==========================
% FFT 周期边界要求区间长度为晶格共同周期的整数倍。
sol.L = 48;                    % 增大计算区间以检查弱态尾部；不改变晶格与驱动
sol.N = round(sol.L/bpar.dx);    % N=960；FFT 不要求 N 是2的整数次幂
sol.t0 = 0;
vphi = 0;                      % 基态相对于激发态的初始相位

% ================= 用户指定激发态初始峰值振幅：只改这一行 =================
sol.ampFactor = sqrt(10);        % A=max(abs(psi))；峰值密度为 A^2
% 增广定态方程同时确定化学势 mu 和波函数 psi，不再指定 mu 或 NB。
% 下面的宽度和化学势偏移仅用于生成不同初猜，不约束最终 mu。
sol.seedWidthFactors = [1 0.7 1.4];
sol.muGuessOffsets = [0 -0.2 0.2]; % 乘 max(1,abs(g)*A^2)，加到初猜 Rayleigh 商上
sol.maxPeakError = 1e-9;        % 峰值绝对容差=此值*max(1,A)
sol.localWindow = 4.0;
sol.minLocalMass = 0.86;
sol.minNorm = 1e-6;             % 排除 Newton 收敛到零解
sol.minWannierOverlap = 0.50;   % 防止明显转到不相关的态；不等价于严格带投影
sol.maxCenterShift = 0.5*lat.a;
sol.boundaryWidth = 2*lat.a;
sol.maxBoundaryMass = 1e-6;
sol.maxResidualAccept = 1e-8;

sol.useEnvelope = true;
sol.wEnvelope = 3.0;                    % 第二能带不要太窄；太窄会破坏 Wannier 节点结构
sol.R0 = job.solitonCenter;      % 构造 Wannier 后更新为实际 Wannier 中心
sol.keepSymmetry = true;
sol.symmetry = 'auto';                  % 'auto','none','even','odd'
sol.lockSignPattern = false;            % 不建议默认锁全部尾部符号；会制造扩展振荡
sol.nodeTol = 1e-3;

sol.NewtonTol = 1e-10;
sol.NewtonMaxIt = 80;
sol.LineSearchMax = 18;
sol.verbose = false;                 % true 会输出大量 Newton 迭代信息

%% ========================== 泵浦演化参数 ==========================
pump.dt = 1e-3;
pump.tmax = 5*pi/lat.nu;
pump.nSave = 1000;
pump.useFourthOrder = true;
pump.plotXLim = [-8, 8];

%% ========================== 基态孤子虚时间搜索参数 ==========================
ground.center = job.groundCenter;
ground.dt = 1e-3;
ground.tmax = 300;
ground.tol = 1e-11;
ground.printEvery = 1000;

% 严格固定初始波函数峰值振幅，不固定粒子数；不截尾、不加人为局域势。
ground.normalizeMode = 'fixedPeak';
ground.ampFactor = sqrt(5);      % A=max(abs(psiGround))；峰值密度为 A^2
ground.seedWidth = 0.35;        % 仅初猜宽度；最终宽度由定态方程决定
ground.targetNB = 3.4;          % 仅在手动选择 fixedNorm 时使用
ground.refineNewton = true;     % 精修时仍固定峰值，mu 和 NB 均可变化


%% ========================== 非线性实空间网格 ==========================
validateattributes(sol.ampFactor, {'numeric'}, {'real','finite','positive','scalar'});
if abs(sol.L/lat.a-round(sol.L/lat.a)) > 1e-10
    error('sol.L 必须是晶格周期 lat.a 的整数倍。');
end
if isempty(sol.seedWidthFactors) || any(sol.seedWidthFactors <= 0) || any(~isfinite(sol.seedWidthFactors))
    error('sol.seedWidthFactors 必须为非空有限正数列表。');
end
if isempty(sol.muGuessOffsets) || any(~isfinite(sol.muGuessOffsets))
    error('sol.muGuessOffsets 必须为非空有限数列表。');
end
pump.t0 = sol.t0;
ground.t0 = sol.t0;
% 能带、定态势、虚时间势和泵浦起始时刻保持一致。
bpar.ky = mod(lat.nu*sol.t0 + (0:bpar.Nky-1).'*pi/bpar.Nky, pi);
[x, k, K2, dx] = make_grid(sol.L, sol.N);
wt = dx * ones(sol.N,1);
if abs(dx - bpar.dx) > 1e-12
    error('非线性网格 dx=%.16g 与 Wannier 网格 dx=%.16g 不一致。', dx, bpar.dx);
end

% 关闭激发态时跳过对应搜索，避免无关的初猜影响运行流程。
V0 = external_potential(x,sol.t0,lat);
Tmat = build_kinetic_matrix(K2);
Eband = []; Coeff = []; Chern = []; Berry = [];
xWannier = []; WannierAll = []; Wcenter = []; Wspread = [];
Wapsi = zeros(size(x)); psi0 = zeros(size(x)); psi = zeros(size(x));
mu = NaN; NB = 0; Energy = 0; R = zeros(size(x));
nodeNumber = 0; overlapW = NaN; peakAmplitude = 0; peakDensity = 0;
outAll = {}; searchInfo = struct('skipped',true);
if job.excitWeight ~= 0 || job.calcChern
fprintf('\n========== Stage 1: Bloch bands and Wannier states ==========%s', newline);
[Eband, Coeff, bpar] = calc_bloch_bands(bpar);

if job.calcChern
    Chern = zeros(bpar.Nband,1);
    Berry = zeros(bpar.Nkx,bpar.Nky,bpar.Nband);
    for ib = 1:bpar.Nband
        [Chern(ib), Berry(:,:,ib)] = chern_fukui_sewing(squeeze(Coeff(:,:,:,ib)), bpar);
        fprintf('Band %d: Chern = %.12f\n', ib, Chern(ib));
    end
else
    Chern = [];
    Berry = [];
end

[xWannier, WannierAll, Wcenter, Wspread] = calc_wannier_all_bands(Coeff, bpar, bpar.iyWannier);
Wfull = WannierAll(:, job.bandIndex);
Wapsi = crop_wannier_to_grid(Wfull, xWannier, x);
Wapsi = preprocess_wannier(Wapsi);
sol.R0 = Wcenter(job.bandIndex);
if abs(sol.R0)+sol.localWindow+sol.boundaryWidth >= sol.L/2
    error('孤子中心距边界过近，请增大 sol.L 或调整 job.excitedCenter。');
end
fprintf('Wannier cell translation R = %.8f; actual envelope center = %.8f\n', ...
        job.excitedCenter, sol.R0);

fprintf('Selected band = %d\n', job.bandIndex);
fprintf('Wannier center(full grid) = %.12g, spread = %.12g\n', Wcenter(job.bandIndex), Wspread(job.bandIndex));

%% ========================== 初猜与化学势初值 ==========================

if sol.keepSymmetry
    sol.symmetry = detect_or_set_symmetry(Wapsi, x, sol.symmetry, sol.R0, V0);
else
    sol.symmetry = 'none';
end
fprintf('Using symmetry = %s about x = %.12g\n', sol.symmetry, sol.R0);
fprintf('Specified peak amplitude A = %.12g; peak density A^2 = %.12g\n', ...
    sol.ampFactor,sol.ampFactor^2);

if sol.useEnvelope
    envelope = 1./cosh((x-sol.R0)/sol.wEnvelope);
else
    envelope = ones(size(x));
end
psi0 = project_solution(Wapsi .* envelope, Wapsi, x, sol);
psi0 = normalize_to_peak(psi0,sol.ampFactor);
NB0 = calc_NB(psi0,wt);
mu0 = real(sum(wt.*psi0.*(Tmat*psi0+(V0+lat.g*psi0.^2).*psi0)))/NB0;

fprintf('\n========== Stage 2: fixed-amplitude excited soliton, band %d seed ==========\n', job.bandIndex);
[psi, mu, outAll, searchInfo] = search_fixed_amplitude_soliton( ...
    Wapsi, x, Tmat, V0, wt, Eband, bpar, lat, sol);
% 峰值约束已经在 Newton 方程中求解；这里不再对解归一化。
peakAmplitude = max(abs(psi));
peakDensity = peakAmplitude^2;
assert(abs(peakAmplitude-sol.ampFactor) <= sol.maxPeakError*max(1,sol.ampFactor), ...
    '内部错误：最终激发态没有满足指定峰值振幅。');

R = residual_stationary(psi, mu, Tmat, V0, lat.g);
NB = calc_NB(psi, wt);
Energy = calc_energy(psi, x, k, wt, V0, lat.g);
nodeNumber = count_nodes(psi);
overlapW = calc_overlap(psi, Wapsi, wt);

fprintf('\n========== Soliton result ==========%s', newline);
fprintf('mu             = %.15g\n', mu);
fprintf('NB             = %.15g\n', NB);
fprintf('peakAmplitude  = %.15g (target %.15g)\n', peakAmplitude,sol.ampFactor);
fprintf('peakDensity    = %.15g\n', peakDensity);
fprintf('Energy         = %.15g\n', Energy);
fprintf('maxResidual    = %.3e\n', max(abs(R)));
fprintf('nodeNumber     = %d\n', nodeNumber);
fprintf('overlap(W,psi) = %.8f\n', overlapW);
locFinal = calc_localization(psi, x, wt, sol.localWindow, sol.R0);
fprintf('localMass      = %.8f within |x-R0| < %.3g\n', locFinal.massFraction, sol.localWindow);
if is_mu_inside_any_band(mu, Eband, bpar.iyWannier)
    warning('最终 mu 落在线性能带内，尾部通常不会局域。请使用 searchInfo 查看候选结果。');
end

else
    fprintf('Excited-state search skipped: job.excitWeight=0.\n');
end

%% ========================== Stage 2b: ground-state soliton ==========================
fprintf('\n========== Stage 2b: ground-state soliton at x = %.3f ==========%s', ground.center, newline);
if job.groundWeight ~= 0
    [psiGround,groundInfo] = search_ground_soliton_imag_time(x,k,wt,lat,ground,Tmat,sol);
else
    psiGround = zeros(size(x)); groundInfo = struct('skipped',true);
end
if job.groundWeight==0 && job.excitWeight==0
    error('两个孤子权重不能同时为0。');
end

psiPump0 = job.excitWeight*psi + job.groundWeight * psiGround*exp(1i*vphi);

fprintf('\n========== Initial combined state ==========%s', newline);
fprintf('Excited soliton center = %.8f, NB = %.12g\n', calc_center(x, psi, wt), calc_NB(psi, wt));
fprintf('Ground  soliton center = %.8f, NB = %.12g\n', calc_center(x, psiGround, wt), calc_NB(psiGround, wt));
fprintf('Ground weight          = %.8f\n', job.groundWeight);
fprintf('Total initial center   = %.8f, NB = %.12g\n', calc_center(x, psiPump0, wt), calc_NB(psiPump0, wt));

%% ========================== Stage 3: nonlinear pumping ==========================
if job.doPump
    fprintf('\n========== Stage 3: nonlinear pumping ==========\n');
    [tdata, udata, densityData, centerX, pumpInfo] = ...
        pump_split_step(psiPump0, x, k, lat, pump, wt);
else
    tdata = 0;
    udata = psiPump0;
    densityData = abs(psiPump0).'.^2;
    centerX = calc_center(x, psiPump0, wt);
    pumpInfo = struct('ran',false);
end

%% ========================== 保存全部工作区 ==========================
% 峰值按各分量实际进入泵浦的初始振幅 max(abs(psi)) 计算（含权重）。
% 关闭的分量峰值记为 0；峰值密度则为对应振幅的平方。
groundPeakAmplitude = max(abs(job.groundWeight*psiGround));
excitedPeakAmplitude = max(abs(job.excitWeight*psi));
groundPeakDensity = groundPeakAmplitude^2;
excitedPeakDensity = excitedPeakAmplitude^2;
workspaceFile = fullfile(pwd, sprintf( ...
    'GroundPeakA_%.12g_ExcitedPeakA_%.12g.mat', ...
    groundPeakAmplitude, excitedPeakAmplitude));
% 不指定变量列表：保存此刻脚本工作区中的全部变量；同名文件会覆盖。
save(workspaceFile, '-v7.3');
fprintf('\nSaved complete workspace: %s\n', workspaceFile);

%% ========================== 统一绘图（数据保存后） ==========================
% 可 load 对应 MAT 文件后单独运行本节；绘图直接位于主脚本中。
xcPlot = real(sum(x(:).*abs(psi(:)).^2)/(sum(abs(psi(:)).^2) + 1e-14));

if ~isempty(Eband)
    figure('Color','w','Name','Bands at initial lattice phase'); hold on; box on;
    for ib = 1:bpar.Nband
        plot(bpar.kx, squeeze(Eband(:,bpar.iyWannier,ib)), 'LineWidth', 1.1);
    end
    xlabel('k_x'); ylabel('E_n(k_x,k_y(t_0))');
    title(sprintf('Bloch bands, selected band = %d', job.bandIndex));
    set(gca,'FontSize',14);
end

if job.excitWeight ~= 0 && max(abs(psi))>0
    figure('Color','w','Name','Wannier and excited soliton');
    subplot(1,2,1); hold on; box on;
    yyaxis left; plot(x, V0, 'LineWidth',1.4); ylabel('V(x,0)');
    yyaxis right; plot(x, abs(psi).^2, 'LineWidth',1.6); ylabel('|\psi|^2');
    xlabel('x'); xlim([xcPlot-5, xcPlot+5]); title('Localized excited soliton'); set(gca,'FontSize',14);
    subplot(1,2,2); hold on; box on;
    plot(x, Wapsi/max(abs(Wapsi)), '--', 'LineWidth',1.4);
    plot(x, real(psi)/max(abs(psi)), '-', 'LineWidth',1.6);
    legend('band Wannier','soliton'); xlabel('x'); ylabel('normalized amplitude');
    xlim([-20,20]); set(gca,'FontSize',14);
end

if numel(tdata)>1
    figure('Position',[100 100 900 900]);
    subplot(1,2,1);
    plot(centerX, tdata, 'LineWidth',1.5);
    xlabel('center x'); ylabel('t'); title('Total-field center'); set(gca,'FontSize',14);
    subplot(1,2,2);
    surf(x, tdata, densityData, 'EdgeColor','none');
    shading interp; view(0,90); colormap("jet"); colorbar;
    xlabel('x'); ylabel('t'); title('Pumping density |\psi(x,t)|^2');
    axis([-20,20, 0, tmax]);
    set(gca,'FontSize',14);
end % pumping plots

if ~isempty(Chern)
    figure('Color','w','Name','Berry curvature');
    ibMax = min(bpar.Nband,4);
    for ib = 1:ibMax
        subplot(2,2,ib);
        surf(bpar.kx, bpar.ky, Berry(:,:,ib).', 'EdgeColor','none');
        shading interp; view(45,30);
        xlabel('k_x'); ylabel('k_y'); zlabel('F');
        title(sprintf('Band %d, C=%.4f', ib, Chern(ib)));
    end
end

figure('Color','w','Name','Initial two-soliton state'); hold on; box on;
plot(x,abs(job.excitWeight*psi).^2,'LineWidth',1.5);
plot(x,abs(job.groundWeight*psiGround).^2,'LineWidth',1.5);
plot(x,abs(psiPump0).^2,'--','LineWidth',1.5);
xlabel('x'); ylabel('|\psi|^2');
legend('excited','ground','combined'); xlim(pump.plotXLim); set(gca,'FontSize',14);

%% ========================================================================
%% 函数区：基础工具
%% ========================================================================
function [x, k, K2, dx] = make_grid(L, N)
    if mod(N,2) ~= 0
        error('FFT 网格点数 sol.N 必须为偶数。');
    end
    dx = L/N;
    x = linspace(-L/2, L/2-dx, N).';
    k = [0:N/2-1, -N/2:-1].' * 2*pi/L;
    K2 = k.^2;
end

function V = external_potential(x, t, lat)
    V = -lat.p1*cos(pi*x/lat.d1).^2 - lat.p2*cos(pi*x/lat.d2 - lat.nu*t).^2;
end

function Tmat = build_kinetic_matrix(K2)
    N = length(K2);
    I = eye(N);
    FI = fft(I, [], 1);
    Tmat = real(0.5 * ifft(bsxfun(@times, K2, FI), [], 1));
    Tmat = 0.5*(Tmat + Tmat');
end

function W = crop_wannier_to_grid(Wfull, xFull, xTarget)
    W = interp1(xFull(:), Wfull(:), xTarget(:), 'pchip', 0);
    if max(abs(W)) < 1e-12
        error('裁剪后的 Wannier 态几乎为零，请检查 xWannier 与非线性网格是否匹配。');
    end
end

function W = preprocess_wannier(W)
    W = W(:);
    [~, imax] = max(abs(W));
    W = W * exp(-1i*angle(W(imax)));
    W = real(W);
    [~, imax] = max(abs(W));
    if W(imax) < 0
        W = -W;
    end
    W = W / max(abs(W));
end

function NB = calc_NB(u, wt)
    NB = real(sum(wt.*abs(u).^2));
end

function u = normalize_to_NB(u, NB, wt)
    nrm = sum(wt.*abs(u).^2);
    if nrm < 1e-14
        error('初态范数太小，无法归一化。');
    end
    u = u * sqrt(NB/nrm);
end

function u = normalize_to_peak(u,A)
    peak = max(abs(u));
    if ~isfinite(peak) || peak <= 0
        error('初猜峰值为零或非有限值，无法设置峰值振幅。');
    end
    u = u*(A/peak);
end

function R = residual_stationary(u, mu, Tmat, Vext, g)
    u = real(u(:));
    R = Tmat*u + (Vext - mu + g*u.^2).*u;
end

function E = calc_energy(psi, x, k, wt, Vext, g)
    psix = ifft(1i*k.*fft(psi));
    Ekin = sum(wt.*0.5.*abs(psix).^2);
    Epot = sum(wt.*Vext.*abs(psi).^2);
    Enl  = sum(wt.*0.5*g.*abs(psi).^4);
    E = real(Ekin + Epot + Enl);
    %#ok<INUSD>
end

function ov = calc_overlap(u, v, wt)
    nu = sqrt(sum(wt.*abs(u).^2));
    nv = sqrt(sum(wt.*abs(v).^2));
    ov = abs(sum(wt.*conj(v).*u))/(nu*nv + 1e-14);
end

function nnode = count_nodes(u)
    u = real(u(:));
    tol = 1e-5 * max(abs(u));
    s = sign(u(abs(u) > tol));
    if isempty(s)
        nnode = 0;
    else
        nnode = sum(s(1:end-1).*s(2:end) < 0);
    end
end

%% ========================================================================
%% 函数区：Bloch 能带、Wannier、Chern
%% ========================================================================
function [Eband, Coeff, par] = calc_bloch_bands(par)
    if mod(par.NG,2) == 0
        error('par.NG 必须取奇数。');
    end
    if mod(par.Nkx,2) ~= 0
        error('par.Nkx 建议取偶数。');
    end

    M = (par.NG-1)/2;
    par.mList = (-M:M).';
    par.G = 2*pi/par.a * par.mList;
    par.qMat = par.mList - par.mList.';

    Eband = zeros(par.Nkx, par.Nky, par.Nband);
    Coeff = zeros(par.NG, par.Nkx, par.Nky, par.Nband);

    for iy = 1:par.Nky
        ky = par.ky(iy);
        Vmat = potential_matrix_analytic(ky, par);
        for ik = 1:par.Nkx
            kx = par.kx(ik);
            H = diag(0.5 * (kx + par.G).^2) + Vmat;
            H = 0.5*(H + H');
            [U,D] = eig(H, 'vector');
            [Eval, idx] = sort(real(D), 'ascend');
            U = U(:,idx);
            Eband(ik,iy,:) = Eval(1:par.Nband);
            Coeff(:,ik,iy,:) = U(:,1:par.Nband);
        end
        fprintf('Bloch spectrum: ky %d / %d\n', iy, par.Nky);
    end
end

function Vmat = potential_matrix_analytic(ky, par)
    qMax = max(abs(par.qMat(:)));
    c = zeros(2*qMax+1,1);

    q1raw = par.a/par.d1;
    q2raw = par.a/par.d2;
    q1 = round(q1raw);
    q2 = round(q2raw);
    if abs(q1-q1raw) > 1e-12 || abs(q2-q2raw) > 1e-12
        error('解析 Fourier 势要求 a/d1 和 a/d2 为整数。');
    end

    add_c(0, -(par.p1+par.p2)/2);
    add_c(+q1, -par.p1/4);
    add_c(-q1, -par.p1/4);
    add_c(+q2, -par.p2/4 * exp(-2i*ky));
    add_c(-q2, -par.p2/4 * exp(+2i*ky));

    Vmat = c(par.qMat + qMax + 1);
    Vmat = 0.5*(Vmat + Vmat');

    function add_c(q,val)
        if abs(q) <= qMax
            c(q + qMax + 1) = c(q + qMax + 1) + val;
        end
    end
end

function [Chern, Berry] = chern_fukui_sewing(Cband, par)
    Berry = zeros(par.Nkx, par.Nky);
    for ik = 1:par.Nkx
        ikp = mod(ik, par.Nkx) + 1;
        wrapK = (ik == par.Nkx);
        for iy = 1:par.Nky
            iyp = mod(iy, par.Nky) + 1;
            u00 = Cband(:,ik,iy);
            u01 = Cband(:,ik,iyp);
            u10 = Cband(:,ikp,iy);
            u11 = Cband(:,ikp,iyp);
            if wrapK
                u10 = shift_plus_G(u10);
                u11 = shift_plus_G(u11);
            end
            Ux00 = link_var(u00,u10);
            Uy10 = link_var(u10,u11);
            Ux01 = link_var(u01,u11);
            Uy00 = link_var(u00,u01);
            Berry(ik,iy) = angle(Ux00 * Uy10 * conj(Ux01) * conj(Uy00));
        end
    end
    Chern = sum(Berry(:))/(2*pi);
end

function z = link_var(v,w)
    ov = v' * w;
    if abs(ov) < 1e-14
        z = 1;
    else
        z = ov / abs(ov);
    end
end

function vS = shift_plus_G(v)
    vS = zeros(size(v));
    vS(1:end-1) = v(2:end);
    vS(end) = 0;
    nrm = norm(vS);
    if nrm > 1e-14
        vS = vS / nrm;
    end
end

function [xWannier, Wannier, Wcenter, Wspread] = calc_wannier_all_bands(Coeff, par, iy)
    NxCell = round(par.a/par.dx);
    if abs(NxCell*par.dx - par.a) > 1e-12
        error('要求 par.a/par.dx 为整数。');
    end
    Ncell = par.NcellWannier;
    Ntot = Ncell * NxCell;
    xCell = -par.a/2 + (0:NxCell-1).' * par.dx;
    cellList = (-Ncell/2 : Ncell/2-1).';

    xWannier = zeros(Ntot,1);
    for ic = 1:Ncell
        idx = (ic-1)*NxCell + (1:NxCell);
        xWannier(idx) = cellList(ic)*par.a + xCell;
    end

    Wannier = zeros(Ntot, par.Nband);
    Wcenter = zeros(par.Nband,1);
    Wspread = zeros(par.Nband,1);

    for ib = 1:par.Nband
        C = squeeze(Coeff(:,:,iy,ib));
        C = smooth_gauge_1d_sewing(C, par);
        W = zeros(Ntot,1);
        for ik = 1:par.Nkx
            kx = par.kx(ik);
            ux = periodic_part_realspace(C(:,ik), xCell, par);
            for ic = 1:Ncell
                idx = (ic-1)*NxCell + (1:NxCell);
                xg = xWannier(idx);
                psiK = exp(1i*kx*xg).*ux;
                W(idx) = W(idx) + exp(-1i*kx*par.R0)*psiK;
            end
        end
        W = W / sqrt(par.Nkx);
        W = W / sqrt(trapz(xWannier, abs(W).^2) + eps);
        rho = abs(W).^2;
        xc = trapz(xWannier, xWannier.*rho);
        sp = trapz(xWannier, (xWannier-xc).^2.*rho);
        Wannier(:,ib) = W;
        Wcenter(ib) = real(xc);
        Wspread(ib) = real(sp);
        fprintf('Wannier band %d: center = %.8g, spread = %.8g\n', ib, Wcenter(ib), Wspread(ib));
    end
end

function Csm = smooth_gauge_1d_sewing(C, par)
    Csm = zeros(size(C));
    Csm(:,1) = C(:,1) / norm(C(:,1));
    for ik = 2:par.Nkx
        v = C(:,ik) / norm(C(:,ik));
        ov = Csm(:,ik-1)' * v;
        Csm(:,ik) = v * exp(-1i*angle(ov));
    end
    firstShift = shift_plus_G(Csm(:,1));
    ovEnd = Csm(:,end)' * firstShift;
    phi = angle(ovEnd);
    if abs(abs(phi)-pi) < 1e-10
        phi = pi; % 统一 Zak 相位为 +/-pi 时的晶胞代表，避免中心随机跳一个晶胞
    end
    for ik = 1:par.Nkx
        alpha = (ik-1)/par.Nkx * phi; % Nk 条链路含首尾拼接链路
        Csm(:,ik) = Csm(:,ik) * exp(1i*alpha);
    end
end

function ux = periodic_part_realspace(Ck, xCell, par)
    ux = exp(1i*xCell*par.G.') * Ck / sqrt(par.a);
end

%% ========================================================================
%% 函数区：激发态孤子 Newton
%% ========================================================================
function symOut = detect_or_set_symmetry(W, x, symIn, R0, Vext)
    if ~any(strcmpi(symIn, {'auto','none','even','odd'}))
        error('sol.symmetry 必须为 auto/none/even/odd。');
    end
    if strcmpi(symIn,'none')
        symOut = 'none'; return;
    end
    dx = x(2)-x(1);
    mapped = (2*R0-x-x(1))/dx;
    if max(abs(mapped-round(mapped))) > 1e-7
        if ~strcmpi(symIn,'auto')
            error('指定对称中心的反射不能精确映射到当前网格。');
        end
        symOut = 'none'; return;
    end
    Vf = reflect_on_grid(Vext,x,R0);
    potentialEven = norm(Vext-Vf,inf) < 1e-9*max(1,norm(Vext,inf));
    if ~potentialEven
        if ~strcmpi(symIn,'auto')
            error('外势关于 sol.R0 不对称，不能强制奇偶对称性。');
        end
        symOut = 'none'; return;
    end
    if ~strcmpi(symIn,'auto')
        symOut = lower(symIn); return;
    end
    Wf = reflect_on_grid(W,x,R0);
    errEven = norm(W-Wf)/max(norm(W),eps);
    errOdd = norm(W+Wf)/max(norm(W),eps);
    if errEven < 0.20 && errEven < errOdd
        symOut = 'even';
    elseif errOdd < 0.20 && errOdd < errEven
        symOut = 'odd';
    else
        symOut = 'none';
    end
end

function uf = reflect_on_grid(u,x,R0)
    % 关于实际中心反射；周期索引，避免 interp1(...,-x,0) 的错位和截断。
    idx = mod(round((2*R0-x-x(1))/(x(2)-x(1))),numel(x))+1;
    uf = u(idx);
end

function u = project_solution(u,Wref,x,sol)
    u = real(u(:));
    if sol.keepSymmetry && ~strcmpi(sol.symmetry,'none')
        uf = reflect_on_grid(u,x,sol.R0);
        if strcmpi(sol.symmetry,'even')
            u = 0.5*(u+uf);
        elseif strcmpi(sol.symmetry,'odd')
            u = 0.5*(u-uf);
        end
    end
    if sol.lockSignPattern
        mask = abs(Wref)>sol.nodeTol*max(abs(Wref));
        sgn = sign(Wref); sgn(sgn==0) = 1;
        u(mask) = sgn(mask).*abs(u(mask));
        u(~mask) = 0;
    end
end

function [psiBest,muBest,outAll,searchInfo] = search_fixed_amplitude_soliton( ...
    Wapsi,x,Tmat,Vext,wt,Eband,bpar,lat,sol)
    A = sol.ampFactor;
    for ib = 1:size(Eband,3)
        fprintf('Band %d at t0: [%.12g, %.12g]\n',ib, ...
            min(Eband(:,bpar.iyWannier,ib)),max(Eband(:,bpar.iyWannier,ib)));
    end
    widthFactors = sol.seedWidthFactors;
    if ~sol.useEnvelope, widthFactors = 1; end
    nTrial = numel(widthFactors)*numel(sol.muGuessOffsets);
    fprintf('Fixed peak amplitude A=%.12g; %d initial-guess trials.\n',A,nTrial);
    psiBest = []; muBest = NaN; bestScore = -inf; bestIndex = NaN;
    outAll = cell(nTrial,1);
    records = struct('seedWidth',{},'muInitial',{},'mu',{},'residual',{}, ...
        'peakAmplitude',{},'peakError',{},'NB',{},'localMass',{},'boundaryMass',{}, ...
        'overlap',{},'nodes',{},'center',{},'inBand',{},'score',{},'accepted',{},'flag',{});
    edge = (x<x(1)+sol.boundaryWidth) | ...
           (x>=x(1)+numel(x)*(x(2)-x(1))-sol.boundaryWidth);
    iTrial = 0;
    for iw = 1:numel(widthFactors)
        width = sol.wEnvelope*widthFactors(iw);
        if sol.useEnvelope
            seedShape = Wapsi./cosh((x-sol.R0)/width);
        else
            seedShape = Wapsi;
        end
        seedShape = project_solution(seedShape,Wapsi,x,sol);
        psi0 = normalize_to_peak(seedShape,A);
        seedNB = calc_NB(psi0,wt);
        muRayleigh = real(sum(wt.*psi0.*(Tmat*psi0+(Vext+lat.g*psi0.^2).*psi0)))/seedNB;
        for im = 1:numel(sol.muGuessOffsets)
            iTrial = iTrial+1;
            mu0 = muRayleigh+sol.muGuessOffsets(im)*max(1,abs(lat.g)*A^2);
            [u,mu,out] = newton_fixed_peak(psi0,mu0,A,Tmat,Vext,wt,Wapsi,x,lat.g,sol);
            res = norm(residual_stationary(u,mu,Tmat,Vext,lat.g),inf);
            NBnow = calc_NB(u,wt);
            peak = max(abs(u)); peakError = abs(peak-A);
            loc = calc_localization(u,x,wt,sol.localWindow,sol.R0);
            ov = calc_overlap(u,Wapsi,wt);
            edgeMass = sum(wt(edge).*abs(u(edge)).^2)/max(NBnow,realmin);
            xc = calc_center(x,u,wt);
            nd = count_nodes_in_window(u,x,sol.localWindow,sol.R0);
            inBand = is_mu_inside_any_band(mu,Eband,bpar.iyWannier);
            covered = mu<max(Eband(:,bpar.iyWannier,end));
            accepted = all(isfinite(u)) && isfinite(mu) && isfinite(res) && ...
                res<sol.maxResidualAccept && peakError<=sol.maxPeakError*max(1,A) && ...
                NBnow>sol.minNorm && loc.massFraction>sol.minLocalMass && ...
                edgeMass<sol.maxBoundaryMass && ov>=sol.minWannierOverlap && ...
                abs(xc-sol.R0)<=sol.maxCenterShift && ~inBand && covered;
            score = ov-0.25*(1-loc.massFraction)-0.05*abs(xc-sol.R0)/lat.a;
            records(iTrial) = struct('seedWidth',width,'muInitial',mu0,'mu',mu, ...
                'residual',res,'peakAmplitude',peak,'peakError',peakError,'NB',NBnow, ...
                'localMass',loc.massFraction,'boundaryMass',edgeMass,'overlap',ov, ...
                'nodes',nd,'center',xc,'inBand',inBand,'score',score, ...
                'accepted',accepted,'flag',out.flag);
            outAll{iTrial} = out;
            fprintf(['trial %d/%d: mu=%.10g, NB=%.9g, peak=%.9g, peakError=%.2e, ' ...
                'residual=%.2e, localMass=%.6f, overlap=%.5f, accepted=%d\n'], ...
                iTrial,nTrial,mu,NBnow,peak,peakError,res,loc.massFraction,ov,accepted);
            if accepted && score>bestScore
                psiBest = u; muBest = mu; bestScore = score; bestIndex = iTrial;
            end
        end
    end
    searchInfo = struct('amplitudeTarget',A,'records',records, ...
        'bestIndex',bestIndex,'bestScore',bestScore);
    if isempty(psiBest)
        disp(struct2table(records));
        error(['给定峰值振幅 A=%.12g 下未找到满足检查条件的局域激发态。' ...
               '可调整 seedWidthFactors、wEnvelope、muGuessOffsets 或扩大网格；' ...
               '本次搜索失败不证明解不存在，程序不会缩放失败态继续演化。'],A);
    end
    searchInfo.outBest = outAll{bestIndex};
    searchInfo.localizationBest = calc_localization(psiBest,x,wt,sol.localWindow,sol.R0);
    fprintf('Accepted fixed-amplitude state: A=%.12g, mu=%.12g, NB=%.12g\n', ...
        max(abs(psiBest)),muBest,calc_NB(psiBest,wt));
end

function tf = is_mu_inside_any_band(mu, Eband, iy)
    tf = false;
    for ib = 1:size(Eband,3)
        Emin = min(squeeze(Eband(:,iy,ib)));
        Emax = max(squeeze(Eband(:,iy,ib)));
        if mu >= Emin && mu <= Emax
            tf = true;
            return;
        end
    end
end

function loc = calc_localization(u, x, wt, window, R0)
    den = sum(wt.*abs(u).^2) + 1e-300;
    mask = abs(x-R0) <= window;
    massIn = sum(wt(mask).*abs(u(mask)).^2);
    loc = struct();
    loc.massFraction = real(massIn/den);
    loc.tailFraction = real(1 - massIn/den);
    loc.window = window;
    loc.R0 = R0;
end

function nnode = count_nodes_in_window(u, x, window, R0)
    mask = abs(x-R0) <= window;
    uu = real(u(mask));
    if isempty(uu)
        nnode = 0;
        return;
    end
    tol = 1e-5 * max(abs(uu));
    s = sign(uu(abs(uu)>tol));
    if isempty(s)
        nnode = 0;
    else
        nnode = sum(s(1:end-1).*s(2:end) < 0);
    end
end

function [psi,mu,out] = newton_fixed_peak(psi0,mu0,A,Tmat,Vext,wt,Wref,x,g,sol)
    % 增广方程：R(u,mu)=0，C(u)=max(abs(u))-A=0。
    % 每步以当前峰值格点构造 C 的一行导数；若峰值位置变化则自动更新。
    u = normalize_to_peak(project_solution(psi0,Wref,x,sol),A);
    mu = mu0; N = numel(u);
    out = struct('errHistory',[],'NBHistory',[],'muHistory',[],'peakHistory',[],'flag',1);
    for it = 1:sol.NewtonMaxIt
        R = residual_stationary(u,mu,Tmat,Vext,g);
        [peak,ip] = max(abs(u)); C = peak-A;
        err = max(norm(R,inf),abs(C));
        out.errHistory(end+1,1) = err; %#ok<AGROW>
        out.NBHistory(end+1,1) = calc_NB(u,wt); %#ok<AGROW>
        out.muHistory(end+1,1) = mu; %#ok<AGROW>
        out.peakHistory(end+1,1) = peak; %#ok<AGROW>
        if err<sol.NewtonTol
            out.flag = 0; break;
        end
        peakRow = zeros(1,N); peakRow(ip) = sign(u(ip));
        J = [Tmat+diag(Vext-mu+3*g*u.^2), -u; peakRow,0];
        step = J\[-R;-C];
        if any(~isfinite(step))
            out.flag = 3; break;
        end
        alpha = 1; accepted = false; merit = norm([R;C],2);
        for ls = 1:sol.LineSearchMax
            unew = project_solution(u+alpha*step(1:N),Wref,x,sol);
            munew = mu+alpha*step(N+1);
            Rnew = residual_stationary(unew,munew,Tmat,Vext,g);
            Cnew = max(abs(unew))-A;
            if all(isfinite(Rnew)) && isfinite(munew) && ...
                    norm([Rnew;Cnew],2)<(1-1e-4*alpha)*merit
                accepted = true; break;
            end
            alpha = alpha/2;
        end
        if ~accepted
            out.flag = 2; break;
        end
        u = unew; mu = munew;
        if sol.verbose
            fprintf('fixed-peak Newton: it=%d, residual=%.3e, mu=%.10g, peak=%.10g\n', ...
                it,norm(Rnew,inf),mu,max(abs(u)));
        end
    end
    % 不在这里重新归一化：峰值必须由增广定态方程自身满足。
    psi = u;
    out.iter = it;
    out.finalResidual = norm(residual_stationary(u,mu,Tmat,Vext,g),inf);
    out.finalPeakError = abs(max(abs(u))-A);
    out.mu = mu;
    if max(out.finalResidual,out.finalPeakError)<sol.NewtonTol
        out.flag = 0;
    end
end

%% ========================================================================
%% 函数区：基态孤子虚时间搜索
%% ========================================================================
function [hw1, info] = search_ground_soliton_imag_time(x, k, wt, lat, ground, Tmat, sol)
    % 虚时间演化搜索静态晶格中的基态孤子，中心位于 ground.center。
    %
    % fixedPeak 在虚时间和 Newton 两个阶段均严格固定 max(abs(hw1))。
    % originalPeak 仅为兼容旧程序保留，不保证峰值等于 ampFactor。
    % 如果需要固定粒子数，可设置：
    %   ground.normalizeMode = 'fixedNorm';
    %   ground.targetNB = 3.4;

    hdt = ground.dt;
    Nt = round(ground.tmax / hdt);

    if isfield(ground,'seedWidth')
        w0 = ground.seedWidth;
    else
        w0 = 1;
    end
    hw1 = exp(-0.5*((x-ground.center)/w0).^2);
    if strcmpi(ground.normalizeMode,'fixedPeak')
        validateattributes(ground.ampFactor,{'numeric'},{'positive','finite','real','scalar'});
        hw1 = normalize_to_peak(hw1,ground.ampFactor);
    end
    hw1 = real(hw1(:));
    ut1 = hw1;

    Vhw = external_potential(x,ground.t0,lat);
    if strcmpi(ground.normalizeMode,'originalPeak')
        fprintf(['Ground uses originalPeak compatibility scaling; ampFactor is not ' ...
                 'an exactly enforced peak.\n']);
    end

    MaxErr = inf;
    for it = 1:Nt
        hw1 = ifft(exp(-hdt/2*k.^2/2).*fft(hw1));
        hw1 = exp(-hdt*(lat.g*abs(hw1).^2 + Vhw)).*hw1;
        hw1 = ifft(exp(-hdt/2*k.^2/2).*fft(hw1));
        hw1 = real(hw1);

        switch lower(ground.normalizeMode)
            case 'originalpeak'
                peak = max(abs(hw1));
                hw1 = hw1 * sqrt(ground.ampFactor/(peak + eps));
            case 'fixedpeak'
                hw1 = normalize_to_peak(hw1,ground.ampFactor);
            case 'fixednorm'
                hw1 = normalize_to_NB(hw1, ground.targetNB, wt);
            otherwise
                error('未知 ground.normalizeMode = %s。可选 fixedPeak、originalPeak 或 fixedNorm。', ground.normalizeMode);
        end

        MaxErr = max(abs(hw1 - ut1));
        ut1 = hw1;

        if isfield(ground,'printEvery') && ground.printEvery > 0 && mod(it, ground.printEvery) == 0
            fprintf('ground imaginary time: it = %d / %d, MaxErr = %.3e, center = %.6f, NB = %.8g\n', ...
                    it, Nt, MaxErr, calc_center(x, hw1, wt), calc_NB(hw1, wt));
        end

        if MaxErr <= ground.tol
            break;
        end
    end

    if MaxErr > ground.tol
        warning('基态虚时间迭代达到上限，MaxErr=%.3e，检查后续定态残差。',MaxErr);
    end
    NBseed = calc_NB(hw1,wt);
    if ~isfinite(NBseed) || NBseed <= 0 || any(~isfinite(hw1))
        error('基态虚时间搜索产生零解或非有限值，请检查参数。');
    end
    muGround = real(sum(wt.*hw1.*(Tmat*hw1+(Vhw+lat.g*hw1.^2).*hw1)))/NBseed;
    initialResidual = norm(residual_stationary(hw1,muGround,Tmat,Vhw,lat.g),inf);
    refineInfo = struct('flag',-1);
    if ground.refineNewton
        if strcmpi(ground.normalizeMode,'fixedPeak')
            groundSol = sol;
            groundSol.keepSymmetry = false;
            groundSol.lockSignPattern = false;
            [hw1,muGround,refineInfo] = newton_fixed_peak(hw1,muGround,ground.ampFactor, ...
                Tmat,Vhw,wt,hw1,x,lat.g,groundSol);
        else
            [hw1,muGround,refineInfo] = refine_ground_fixed_NB(hw1,muGround, ...
                NBseed,Tmat,Vhw,wt,lat.g,sol);
        end
    end
    if strcmpi(ground.normalizeMode,'fixedPeak') && ...
            abs(max(abs(hw1))-ground.ampFactor)>sol.maxPeakError*max(1,ground.ampFactor)
        error('基态没有达到指定峰值振幅。');
    end
    finalResidual = norm(residual_stationary(hw1,muGround,Tmat,Vhw,lat.g),inf);
    if ~isfinite(finalResidual) || finalResidual > sol.maxResidualAccept
        error('基态未通过定态残差检查：%.3e。请减小 ground.dt 或增加迭代次数。',finalResidual);
    end
    info = struct();
    info.mu = muGround;
    info.initialResidual = initialResidual;
    info.residual = finalResidual;
    info.refineInfo = refineInfo;
    info.iter = it;
    info.maxErr = MaxErr;
    info.center = calc_center(x, hw1, wt);
    info.NB = calc_NB(hw1, wt);
    info.peak = max(abs(hw1));
    info.normalizeMode = ground.normalizeMode;
    centerMask = abs(x-info.center)<=lat.a/2;
    tailMask = ~centerMask;
    info.massOutsideCentralCell = sum(wt(tailMask).*abs(hw1(tailMask)).^2)/info.NB;
    info.rmsWidth = sqrt(sum(wt.*(x-info.center).^2.*abs(hw1).^2)/info.NB);
    sig = hw1(abs(hw1)>1e-6*max(abs(hw1)));
    info.nodeless = all(sig>=0) || all(sig<=0);
    if ~info.nodeless
        error('搜索结果具有显著节点，不能作为目标无节点基态使用。');
    end

    fprintf('Ground soliton done: iter = %d, MaxErr = %.3e, center = %.8f, NB = %.12g, peak = %.8g\n', ...
            info.iter, info.maxErr, info.center, info.NB, info.peak);
    fprintf('Ground mu=%.12g; residual before/after refinement=%.3e / %.3e\n', ...
        info.mu,info.initialResidual,info.residual);
    fprintf('Ground RMS width=%.8g; mass outside central cell=%.6f; nodeless=%d\n', ...
        info.rmsWidth,info.massOutsideCentralCell,info.nodeless);
    fprintf('The decaying lattice tails are retained; stationary residual is checked without a cutoff.\n');
end

%% ========================================================================
%% 函数区：Split-step 泵浦演化
%% ========================================================================
function [u,mu,out] = refine_ground_fixed_NB(u,mu,NBtarget,Tmat,Vext,wt,g,sol)
    % 只精修基态，固定虚时间所得粒子数；不用于激发态搜索。
    out = struct('flag',1,'errHistory',[]);
    for it = 1:sol.NewtonMaxIt
        R = residual_stationary(u,mu,Tmat,Vext,g);
        C = calc_NB(u,wt)-NBtarget;
        out.errHistory(end+1,1) = max(norm(R,inf),abs(C)); %#ok<AGROW>
        if out.errHistory(end)<sol.NewtonTol
            out.flag = 0; break;
        end
        J = [Tmat+diag(Vext-mu+3*g*u.^2), -u; (2*wt.*u).',0];
        step = J\[-R;-C];
        if any(~isfinite(step)), break; end
        alpha = 1; accepted = false;
        merit = norm([R;C],2);
        for ls = 1:sol.LineSearchMax
            unew = u+alpha*step(1:end-1); munew = mu+alpha*step(end);
            Rnew = residual_stationary(unew,munew,Tmat,Vext,g);
            Cnew = calc_NB(unew,wt)-NBtarget;
            if all(isfinite(Rnew)) && norm([Rnew;Cnew],2)<(1-1e-4*alpha)*merit
                accepted = true; break;
            end
            alpha = alpha/2;
        end
        if ~accepted, break; end
        u = unew; mu = munew;
    end
    out.iter = it;
    out.normError = abs(calc_NB(u,wt)-NBtarget);
    if ~isfinite(out.normError) || out.normError > sol.maxResidualAccept
        error('基态 Newton 精修未保持粒子数：误差 %.3e。',out.normError);
    end
end

function [tdata,udata,densityData,centerX,info] = pump_split_step(psi0,x,k,lat,pump,wt)
    validateattributes(pump.dt,{'numeric'},{'positive','finite','scalar'});
    validateattributes(pump.tmax,{'numeric'},{'positive','finite','scalar'});
    nmax = ceil(pump.tmax/pump.dt);
    dt = pump.tmax/nmax; % 恰好到达 tmax
    saveStep = max(1,round(nmax/pump.nSave));
    saveIndices = unique([0:saveStep:nmax,nmax]);
    nFrame = numel(saveIndices); N = numel(x);
    u = psi0(:);
    udata = complex(zeros(N,nFrame));
    tdata = zeros(1,nFrame); densityData = zeros(nFrame,N);
    centerX = zeros(1,nFrame); normData = zeros(1,nFrame);
    frame = 1; udata(:,frame) = u;
    densityData(frame,:) = abs(u).'.^2;
    centerX(frame) = calc_center(x,u,wt); normData(frame) = calc_NB(u,wt);
    if pump.useFourthOrder
        c = 1/(2-2^(1/3));
        a = [c/2,(1-c)/2,(1-c)/2,c/2]; b = [c,1-2*c,c];
    else
        a = [1/2,1/2]; b = 1;
    end
    E = cell(numel(a),1);
    for j = 1:numel(a)
        E{j} = exp(-1i*a(j)*dt*k.^2/2);
    end
    stageTimes = cumsum(a);
    for nn = 1:nmax
        tStart = pump.t0+(nn-1)*dt;
        u = ifft(fft(u).*E{1});
        for j = 1:numel(b)
            % 在扩展时间空间作对称分裂：每个势子步使用自己的阶段时刻。
            V = external_potential(x,tStart+stageTimes(j)*dt,lat);
            u = u.*exp(-1i*b(j)*dt*(V+lat.g*abs(u).^2));
            u = ifft(fft(u).*E{j+1});
        end
        if frame<nFrame && nn==saveIndices(frame+1)
            frame = frame+1; udata(:,frame) = u; tdata(frame) = nn*dt;
            densityData(frame,:) = abs(u).'.^2;
            centerX(frame) = calc_center(x,u,wt); normData(frame) = calc_NB(u,wt);
            if mod(frame,50)==0 || frame==nFrame
                fprintf('pump: %.6f / %.6f, relative norm drift=%.2e\n', ...
                    nn*dt,pump.tmax,abs(normData(frame)-normData(1))/max(normData(1),eps));
            end
        end
    end
    edgeMask = abs(x)>max(abs(x))-2*lat.a;
    edgeMass = densityData(:,edgeMask)*wt(edgeMask)./max(normData(:),eps);
    firstEdge = find(edgeMass>0.01,1);
    info = struct('ran',true,'dtActual',dt,'normData',normData, ...
        'edgeMass',edgeMass,'maxRelativeNormDrift',max(abs(normData-normData(1)))/max(normData(1),eps));
    if ~isempty(firstEdge)
        warning('t=%.6g 时边缘区域已含超过1%%的粒子数，之后的总质心可能受周期边界回绕影响。',tdata(firstEdge));
    end
    % centerX 为整个场的质心，碰撞期间不是每个孤子的独立轨迹。
end

function xc = calc_center(x, u, wt)
    den = sum(wt.*abs(u).^2);
    xc = real(sum(wt.*x.*abs(u).^2)/(den + 1e-14));
end
