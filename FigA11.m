%% Integrated_SecondBand_ExcitedSoliton_Pump_Localized_Fixed.m
% 一体化程序：线性能带/Wannier 态 -> 第二能带激发态孤子 -> 非线性泵浦动力学
%
% 默认示例：第二能带 bandIndex = 2。
% 你只需要改 bandIndex，就可以换成第一、第三、第四能带。
%
% 方程：

%   i psi_t = [-1/2 d_xx + V(x,t) + g |psi|^2] psi,
%   V(x,t) = -p1 cos^2(pi*x/d1) - p2 cos^2(pi*x/d2 - nu*t).
%
% 流程：
%   1. Bloch 平面波基底计算 E_n(kx,ky)；
%   2. 固定 ky=0 构造第 bandIndex 条能带 Wannier 态；
%   3. 用该 Wannier 态作为初猜，固定粒子数增广 Newton 搜索激发态孤子；
%   4. 以该激发态孤子为初态，在滑动晶格中做 split-step 泵浦演化。

clear; clc; close all;

%% ========================== 总开关 ==========================
job.bandIndex = 2;              % 第二能带为例
job.fastMode = true;            % true: 先用较小网格跑通；正式计算可改 false
job.calcChern = false;          % Chern 对泵浦演化不是必须；需要时改 true
job.savePrefix = 'SecondBand_ExcitedSoliton_x4_GroundSoliton_xm3p5_Pump';
job.excitedCenter = 3;         % 第二能带激发态孤子中心位置：x = 4
job.groundCenter  = -2;      % 虚时间搜索得到的基态孤子中心位置：x = -3.5
job.groundWeight   = 1;         % 基态孤子权重；设为0时，只演化激发态，结果应与单独跑激发态一致
job.excitWeight   = 1;    % 激发态孤子权重；设为0时，只演化激发态，结果应与单独跑激发态一致
job.solitonCenter = job.excitedCenter;  % 兼容原有变量名：激发态孤子中心

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

bpar.iyWannier = 1;             % ky=0 处构造 Wannier 态
bpar.dx = 0.05;                 % 和非线性泵浦网格保持一致
bpar.NxCell = round(bpar.a/bpar.dx);
bpar.NcellWannier = bpar.Nkx;
bpar.R0 = job.solitonCenter;     % Wannier 态中心位置

%% ========================== 孤子求解参数 ==========================
sol.L = 25.6;
sol.N = round(sol.L/bpar.dx);    % dx=0.05 时 N=512
sol.t0 = 0; vphi = pi*0;

% 关键改动：激发态孤子必须把化学势放在相应线性能隙中。
% 原先 fixedNB 会让 mu 自由跑，容易落入线性能带，得到扩展态/多峰态。
sol.solitonSearchMode = 'gapFixedMu';   % 推荐：'gapFixedMu'；不再默认 fixedNB
sol.gapChoice = 'both';                 % 'below','above','both'，第二能带建议 both 自动尝试
sol.muSearchN = 5;                      % 每个能隙中试探几个 mu；正式计算可改 9
sol.gapMargin = 0.08;                   % 避开能带边缘的比例，防止病态
sol.seedNBList = [0.1 0.3 0.7 1.5 3 5];  % 初猜范数扫描；不宜太多，否则很慢
sol.localWindow = 4.0;                  % 判断局域性的窗口 |x-R0|<localWindow
sol.minLocalMass = 0.86;                % 至少多少粒子数留在局域窗口内
sol.maxResidualAccept = 1e-7;           % 固定 mu Newton 可接受残差
sol.targetNB = [];                      % 默认不固定 NB；要指定粒子数，先找到局域分支后再延拓
sol.doFixedNBAfterSearch = false;       % 若一定要指定 NB，可改 true
sol.NBPath = [];                        % 如 doFixedNBAfterSearch=true，可设 linspace(NB_found,targetNB,n)

sol.muGuessMode = 'bandedge';           % 只作为 fallback 使用
sol.muOffset = 0.3;
sol.muManual = -50;

sol.useEnvelope = true;
sol.wEnvelope = 3.0;                    % 第二能带不要太窄；太窄会破坏 Wannier 节点结构
sol.R0 = job.solitonCenter;       % 孤子包络和局域窗口中心位置
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
T = pi/lat.nu;
pump.nSave = 1000;
pump.useFourthOrder = true;
pump.plotXLim = [-8, 8];

%% ========================== 基态孤子虚时间搜索参数 ==========================
ground.center = job.groundCenter;      % 基态孤子中心：x = -3.5
ground.dt = 1e-3;
ground.tmax = 300;
ground.tol = 1e-11;
ground.printEvery = 1000;

% 保留你原程序的归一化方式：
% hw1 = hw1 * sqrt(af/max(abs(hw1))), af = sqrt(3.4)。
% 若后续想固定粒子数，可把 normalizeMode 改为 'fixedNorm'，并设置 targetNB。
ground.normalizeMode = 'originalPeak';
ground.ampFactor = sqrt(10);
ground.targetNB = 3.4;


%% ========================== 非线性实空间网格 ==========================
[x, k, K2, dx] = make_grid(sol.L, sol.N);
wt = dx * ones(sol.N,1);
if abs(dx - bpar.dx) > 1e-12
    error('非线性网格 dx=%.16g 与 Wannier 网格 dx=%.16g 不一致。', dx, bpar.dx);
end

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

fprintf('Selected band = %d\n', job.bandIndex);
fprintf('Wannier center(full grid) = %.12g, spread = %.12g\n', Wcenter(job.bandIndex), Wspread(job.bandIndex));

%% ========================== 初猜与化学势初值 ==========================
V0 = external_potential(x, sol.t0, lat);
Tmat = build_kinetic_matrix(K2);

if strcmpi(sol.muGuessMode, 'bandedge')
    Eb = squeeze(Eband(:, bpar.iyWannier, job.bandIndex));
    mu0 = min(Eb) - sol.muOffset;
else
    mu0 = sol.muManual;
end

if sol.keepSymmetry
    sol.symmetry = detect_or_set_symmetry(Wapsi, x, sol.symmetry);
end
fprintf('Using symmetry = %s\n', sol.symmetry);
fprintf('Initial mu guess = %.12g\n', mu0);

if sol.useEnvelope
    envelope = 1./cosh((x-sol.R0)/sol.wEnvelope);
else
    envelope = ones(size(x));
end
psi0 = Wapsi .* envelope;
% 仅用于显示和保存的初猜。正式搜索在 search_localized_gap_soliton 中会按 seedNBList 逐个归一化。
if isempty(sol.NBPath)
    psi0NB = sol.seedNBList(1);
else
    psi0NB = sol.NBPath(1);
end
psi0 = normalize_to_NB(psi0, psi0NB, wt);
psi0 = project_solution(psi0, Wapsi, x, sol);

fprintf('\n========== Stage 2: localized excited soliton from band-%d Wannier ==========%s', job.bandIndex, newline);

% 重要：先在对应线性能隙中固定 mu 搜索局域态。
% 若 mu 落在线性能带内，尾部不会指数衰减，Newton 很容易得到扩展多峰态。
[psi, mu, outAll, searchInfo] = search_localized_gap_soliton(Wapsi, x, Tmat, V0, wt, Eband, bpar, job, lat, sol);

% 可选：已经找到局域分支后，再用固定 NB 做短延拓。
% 注意：如果固定 NB 延拓把 mu 推回线性能带，态会重新变成扩展态。
if isfield(sol,'doFixedNBAfterSearch') && sol.doFixedNBAfterSearch
    if isempty(sol.targetNB)
        error('doFixedNBAfterSearch=true 时必须设置 sol.targetNB。');
    end
    NBstart = calc_NB(psi, wt);
    if isempty(sol.NBPath)
        NBPath = linspace(NBstart, sol.targetNB, 6);
    else
        NBPath = sol.NBPath(:).';
        if abs(NBPath(1)-NBstart) > 1e-12
            NBPath = [NBstart NBPath];
        end
        if abs(NBPath(end)-sol.targetNB) > 1e-12
            NBPath(end+1) = sol.targetNB;
        end
    end
    for j = 2:numel(NBPath)
        NBnow = NBPath(j);
        psi = normalize_to_NB(psi, NBnow, wt);
        psi = project_solution(psi, Wapsi, x, sol);
        fprintf('\n---- Local branch fixed-NB continuation %d/%d: target NB = %.12g ----\n', j-1, numel(NBPath)-1, NBnow);
        [psi, mu, outNB] = newton_fixed_NB(psi, mu, NBnow, Tmat, V0, wt, Wapsi, x, lat.g, sol);
        outAll{end+1} = outNB; %#ok<SAGROW>
        locTmp = calc_localization(psi, x, wt, sol.localWindow, sol.R0);
        fprintf('After fixed-NB continuation: mu = %.12g, localMass = %.6f\n', mu, locTmp.massFraction);
        if is_mu_inside_any_band(mu, Eband, bpar.iyWannier)
            warning('当前 mu=%.12g 已落入线性能带，继续延拓可能得到非局域态。', mu);
        end
    end
end

R = residual_stationary(psi, mu, Tmat, V0, lat.g);
NB = calc_NB(psi, wt);
Energy = calc_energy(psi, x, k, wt, V0, lat.g);
nodeNumber = count_nodes(psi);
overlapW = calc_overlap(psi, Wapsi, wt);

fprintf('\n========== Soliton result ==========%s', newline);
fprintf('mu             = %.15g\n', mu);
fprintf('NB             = %.15g\n', NB);
fprintf('Energy         = %.15g\n', Energy);
fprintf('maxResidual    = %.3e\n', max(abs(R)));
fprintf('nodeNumber     = %d\n', nodeNumber);
fprintf('overlap(W,psi) = %.8f\n', overlapW);
locFinal = calc_localization(psi, x, wt, sol.localWindow, sol.R0);
fprintf('localMass      = %.8f within |x-R0| < %.3g\n', locFinal.massFraction, sol.localWindow);
if is_mu_inside_any_band(mu, Eband, bpar.iyWannier)
    warning('最终 mu 落在线性能带内，尾部通常不会局域。请使用 searchInfo 查看候选结果。');
end

%% ========================== Stage 2b: ground-state soliton at x = -3.5 ==========================
fprintf('\n========== Stage 2b: ground-state soliton at x = %.3f ==========%s', ground.center, newline);
[psiGround, groundInfo] = search_ground_soliton_imag_time(x, k, wt, lat, ground);

psiPump0 = job.excitWeight*psi + job.groundWeight * psiGround*exp(1i*vphi);

fprintf('\n========== Initial combined state ==========%s', newline);
fprintf('Excited soliton center = %.8f, NB = %.12g\n', calc_center(x, psi, wt), calc_NB(psi, wt));
fprintf('Ground  soliton center = %.8f, NB = %.12g\n', calc_center(x, psiGround, wt), calc_NB(psiGround, wt));
fprintf('Ground weight          = %.8f\n', job.groundWeight);
fprintf('Total initial center   = %.8f, NB = %.12g\n', calc_center(x, psiPump0, wt), calc_NB(psiPump0, wt));
% 
% figure('Color','w','Name','Initial two-soliton state'); hold on; box on;
% plot(x, (psi), 'LineWidth', 1.5);
% plot(x, abs(job.groundWeight*psiGround).^2, 'LineWidth', 1.5);
% plot(x, (psiPump0), '--', 'LineWidth', 1.5);
% xlabel('x'); ylabel('|\psi|^2');
% legend('excited soliton at x=4', 'ground soliton at x=-3.5', 'combined initial state');
% xlim([-8, 8]); set(gca,'FontSize',14);

%% ========================== Stage 3: nonlinear pumping ==========================
fprintf('\n========== Stage 3: nonlinear pumping ==========%s', newline);
[tdata, udata, densityData, centerX] = pump_split_step(psiPump0, x, k, lat, pump, wt);

%% ========================== 作图 ==========================
plot_all_results(x, k, bpar, Eband, xWannier, Wapsi, psi, V0, tdata, densityData, centerX, Chern, Berry, job, pump);

%% ========================== 保存 ==========================
% save([job.savePrefix '_Result.mat'], ...
%     'job','lat','bpar','sol','pump','x','k','dx','Eband','Coeff','Chern','Berry', ...
%     'xWannier','WannierAll','Wcenter','Wspread','Wapsi','psi0','psi','mu','NB','Energy', ...
%     'nodeNumber','overlapW','R','outAll','searchInfo','ground','psiGround','psiPump0','groundInfo', ...
%     'tdata','udata','densityData','centerX','-v7.3');
% 
% save('Wapsi.mat','Wapsi');
% save('psi.mat','psi','mu','NB','Energy','nodeNumber');
% save('psiGround.mat','psiGround','groundInfo');
% 
% fprintf('\nAll done. Result saved to %s_Result.mat\n', job.savePrefix);

%% ========================================================================
%% 函数区：基础工具
%% ========================================================================
function [x, k, K2, dx] = make_grid(L, N)
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
    for ik = 1:par.Nkx
        alpha = (ik-1)/(par.Nkx-1) * phi;
        Csm(:,ik) = Csm(:,ik) * exp(1i*alpha);
    end
end

function ux = periodic_part_realspace(Ck, xCell, par)
    ux = exp(1i*xCell*par.G.') * Ck / sqrt(par.a);
end

%% ========================================================================
%% 函数区：激发态孤子 Newton
%% ========================================================================
function symOut = detect_or_set_symmetry(W, x, symIn)
    if ~strcmpi(symIn, 'auto')
        symOut = lower(symIn);
        return;
    end
    Wf = reflect_on_grid(W, x);
    nrm = norm(W) + 1e-14;
    errEven = norm(W - Wf)/nrm;
    errOdd  = norm(W + Wf)/nrm;
    if errEven < 0.20 && errEven < errOdd
        symOut = 'even';
    elseif errOdd < 0.20 && errOdd < errEven
        symOut = 'odd';
    else
        symOut = 'none';
    end
end

function uf = reflect_on_grid(u, x)
    uf = interp1(x(:), u(:), -x(:), 'linear', 0);
end

function u = project_solution(u, Wref, x, sol)
    u = real(u(:));
    if sol.keepSymmetry
        uf = reflect_on_grid(u, x);
        switch lower(sol.symmetry)
            case 'even'
                u = 0.5*(u + uf);
            case 'odd'
                u = 0.5*(u - uf);
            otherwise
        end
    end
    if sol.lockSignPattern
        mask = abs(Wref) > sol.nodeTol*max(abs(Wref));
        s = sign(Wref);
        s(s==0) = 1;
        u(mask) = s(mask).*abs(u(mask));
        u(~mask) = 0;
    end
end

function [psi, mu, out] = newton_fixed_NB(psi0, mu0, targetNB, Tmat, Vext, wt, Wref, x, g, sol)
    u = project_solution(psi0, Wref, x, sol);
    mu = mu0;
    N = length(u);
    out.errHistory = [];
    out.muHistory = [];
    out.NBHistory = [];
    out.flag = 1;

    for it = 1:sol.NewtonMaxIt
        R = residual_stationary(u, mu, Tmat, Vext, g);
        C = calc_NB(u, wt) - targetNB;
        err = max([max(abs(R)), abs(C)]);
        out.errHistory(end+1,1) = err; %#ok<AGROW>
        out.muHistory(end+1,1) = mu; %#ok<AGROW>
        out.NBHistory(end+1,1) = calc_NB(u, wt); %#ok<AGROW>

        if sol.verbose
            fprintf('it = %3d, err = %.3e, mu = %.12g, NB = %.12g, nodes = %d\n', ...
                    it, err, mu, calc_NB(u, wt), count_nodes(u));
        end
        if err < sol.NewtonTol
            out.flag = 0;
            break;
        end

        Lmat = Tmat + diag(Vext - mu + 3*g*u.^2);
        Aug = [Lmat, -u; (2*wt.*u).', 0];
        rhs = [-R; -C];
        solvec = Aug \ rhs;
        du = solvec(1:N);
        dmu = solvec(N+1);

        alpha = 1.0;
        accepted = false;
        for ls = 1:sol.LineSearchMax
            unew = project_solution(u + alpha*du, Wref, x, sol);
            munew = mu + alpha*dmu;
            Rnew = residual_stationary(unew, munew, Tmat, Vext, g);
            Cnew = calc_NB(unew, wt) - targetNB;
            errNew = max([max(abs(Rnew)), abs(Cnew)]);
            if errNew < err || alpha < 1e-4
                accepted = true;
                break;
            end
            alpha = alpha/2;
        end
        if ~accepted
            warning('Line search 未找到明显下降步长，使用 alpha = %.3e。', alpha);
        end
        u = unew;
        mu = munew;
    end
    psi = u;
end


function [psiBest, muBest, outAll, searchInfo] = search_localized_gap_soliton(Wapsi, x, Tmat, Vext, wt, Eband, bpar, job, lat, sol)
    % 在目标能带上下能隙中固定 mu 搜索局域激发态。
    % 选择标准：残差小、窗口内粒子数比例高、与 Wapsi 重叠高。

    if sol.keepSymmetry
        sol.symmetry = detect_or_set_symmetry(Wapsi, x, sol.symmetry);
    end

    if sol.useEnvelope
        envelope = 1./cosh((x-sol.R0)/sol.wEnvelope);
    else
        envelope = ones(size(x));
    end
    seedShape = project_solution(Wapsi(:).*envelope, Wapsi, x, sol);

    muList = make_gap_mu_list(Eband, bpar, job.bandIndex, sol);
    if isempty(muList)
        warning('没有找到可用线性能隙，退回 bandedge/manual 初猜。');
        if strcmpi(sol.muGuessMode,'bandedge')
            Eb = squeeze(Eband(:, bpar.iyWannier, job.bandIndex));
            muList = min(Eb) - sol.muOffset;
        else
            muList = sol.muManual;
        end
    end

    fprintf('Trying %d mu candidates in linear gaps.\n', numel(muList));
    fprintf('mu candidates: '); fprintf('%.6g ', muList); fprintf('\n');

    bestScore = inf;
    psiBest = [];
    muBest = NaN;
    outBest = [];
    outAll = {};
    records = struct('mu',{},'seedNB',{},'residual',{},'NB',{},'localMass',{},'overlap',{},'nodes',{},'score',{},'accepted',{});

    recId = 0;
    for im = 1:numel(muList)
        muTry = muList(im);
        if is_mu_inside_any_band(muTry, Eband, bpar.iyWannier)
            fprintf('Skip mu %.12g because it is inside a linear band.\n', muTry);
            continue;
        end
        for iseed = 1:numel(sol.seedNBList)
            seedNB = sol.seedNBList(iseed);
            psi0 = normalize_to_NB(seedShape, seedNB, wt);
            psi0 = project_solution(psi0, Wapsi, x, sol);

            fprintf('\n---- fixed-mu trial: mu = %.12g, seedNB = %.12g ----\n', muTry, seedNB);
            [psiTry, outTry] = newton_fixed_mu(psi0, muTry, Tmat, Vext, wt, Wapsi, x, lat.g, sol);
            Rtry = residual_stationary(psiTry, muTry, Tmat, Vext, lat.g);
            res = max(abs(Rtry));
            NBtry = calc_NB(psiTry, wt);
            loc = calc_localization(psiTry, x, wt, sol.localWindow, sol.R0);
            ov = calc_overlap(psiTry, Wapsi, wt);
            nd = count_nodes_in_window(psiTry, x, sol.localWindow, sol.R0);

            % 分数越小越好。首先要求局域，然后再看重叠和残差。
            score = log10(res + 1e-300) ...
                    + 8*(1-loc.massFraction) ...
                    + 0.5*(1-ov) ...
                    + 0.02*abs(NBtry-seedNB);

            accepted = (res < sol.maxResidualAccept) && (loc.massFraction > sol.minLocalMass) && (NBtry > 1e-8);
            recId = recId + 1;
            records(recId).mu = muTry;
            records(recId).seedNB = seedNB;
            records(recId).residual = res;
            records(recId).NB = NBtry;
            records(recId).localMass = loc.massFraction;
            records(recId).overlap = ov;
            records(recId).nodes = nd;
            records(recId).score = score;
            records(recId).accepted = accepted;
            outAll{end+1} = outTry; %#ok<AGROW>

            fprintf('trial result: residual=%.3e, NB=%.6g, localMass=%.4f, overlap=%.4f, nodes(win)=%d\n', ...
                    res, NBtry, loc.massFraction, ov, nd);

            if accepted && score < bestScore
                bestScore = score;
                psiBest = psiTry;
                muBest = muTry;
                outBest = outTry;
            end
        end
    end

    if isempty(psiBest)
        % 如果严格标准下没有通过，选 localMass 最高且残差相对较小的作为诊断结果。
        warning('没有找到满足局域性阈值的固定-mu 解。将返回最佳诊断态；请查看 searchInfo.records。');
        if isempty(records)
            error('所有 mu 搜索均失败，没有可返回的态。');
        end
        vals = [[records.localMass].' , -log10([records.residual].'+1e-300), [records.overlap].'];
        metric = 6*vals(:,1) + vals(:,2) + vals(:,3);
        [~,ii] = max(metric);
        muBest = records(ii).mu;
        seedNB = records(ii).seedNB;
        psi0 = normalize_to_NB(seedShape, seedNB, wt);
        [psiBest, outBest] = newton_fixed_mu(psi0, muBest, Tmat, Vext, wt, Wapsi, x, lat.g, sol);
    end

    locBest = calc_localization(psiBest, x, wt, sol.localWindow, sol.R0);
    fprintf('\nSelected localized branch: mu=%.12g, NB=%.12g, localMass=%.6f, residual=%.3e\n', ...
            muBest, calc_NB(psiBest, wt), locBest.massFraction, max(abs(residual_stationary(psiBest, muBest, Tmat, Vext, lat.g))));

    searchInfo = struct();
    searchInfo.records = records;
    searchInfo.bestScore = bestScore;
    searchInfo.outBest = outBest;
    searchInfo.localizationBest = locBest;
    searchInfo.muList = muList;
    searchInfo.seedNBList = sol.seedNBList;
end

function muList = make_gap_mu_list(Eband, bpar, bandIndex, sol)
    iy = bpar.iyWannier;
    muList = [];
    nband = size(Eband,3);

    if bandIndex > 1 && (strcmpi(sol.gapChoice,'below') || strcmpi(sol.gapChoice,'both'))
        lo = max(squeeze(Eband(:,iy,bandIndex-1)));
        hi = min(squeeze(Eband(:,iy,bandIndex)));
        muList = [muList; points_in_gap(lo, hi, sol.muSearchN, sol.gapMargin)]; %#ok<AGROW>
        fprintf('Gap below band %d: (%.12g, %.12g), width %.6g\n', bandIndex, lo, hi, hi-lo);
    end

    if bandIndex < nband && (strcmpi(sol.gapChoice,'above') || strcmpi(sol.gapChoice,'both'))
        lo = max(squeeze(Eband(:,iy,bandIndex)));
        hi = min(squeeze(Eband(:,iy,bandIndex+1)));
        muList = [muList; points_in_gap(lo, hi, sol.muSearchN, sol.gapMargin)]; %#ok<AGROW>
        fprintf('Gap above band %d: (%.12g, %.12g), width %.6g\n', bandIndex, lo, hi, hi-lo);
    end

    muList = unique(muList(:));
end

function pts = points_in_gap(lo, hi, n, margin)
    if hi <= lo
        pts = [];
        return;
    end
    width = hi - lo;
    a = lo + margin*width;
    b = hi - margin*width;
    if b <= a
        pts = 0.5*(lo+hi);
    else
        pts = linspace(a, b, n).';
    end
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

function [psi, out] = newton_fixed_mu(psi0, mu, Tmat, Vext, wt, Wref, x, g, sol)
    u = project_solution(psi0, Wref, x, sol);
    out.errHistory = [];
    out.NBHistory = [];
    out.flag = 1;

    for it = 1:sol.NewtonMaxIt
        R = residual_stationary(u, mu, Tmat, Vext, g);
        err = max(abs(R));
        out.errHistory(end+1,1) = err; %#ok<AGROW>
        out.NBHistory(end+1,1) = calc_NB(u, wt); %#ok<AGROW>

        if sol.verbose
            loc = calc_localization(u, x, wt, sol.localWindow, sol.R0);
            fprintf('it = %3d, err = %.3e, NB = %.12g, localMass = %.4f, nodes(win) = %d\n', ...
                    it, err, calc_NB(u, wt), loc.massFraction, count_nodes_in_window(u,x,sol.localWindow,sol.R0));
        end

        if err < sol.NewtonTol
            out.flag = 0;
            break;
        end

        Lmat = Tmat + diag(Vext - mu + 3*g*u.^2);
        du = Lmat \ (-R);

        alpha = 1.0;
        errOld = err;
        accepted = false;
        for ls = 1:sol.LineSearchMax
            unew = project_solution(u + alpha*du, Wref, x, sol);
            Rnew = residual_stationary(unew, mu, Tmat, Vext, g);
            errNew = max(abs(Rnew));
            if errNew < errOld || alpha < 1e-4
                accepted = true;
                break;
            end
            alpha = alpha/2;
        end
        if ~accepted
            warning('fixed-mu line search 未找到明显下降步长，使用 alpha = %.3e。', alpha);
        end
        u = unew;
    end
    psi = u;
end


%% ========================================================================
%% 函数区：基态孤子虚时间搜索
%% ========================================================================
function [hw1, info] = search_ground_soliton_imag_time(x, k, wt, lat, ground)
    % 虚时间演化搜索静态晶格中的基态孤子，中心位于 ground.center。
    %
    % 默认保留你原来的峰值归一化方式：
    %   hw1 = hw1 * sqrt(af/max(abs(hw1))), af = sqrt(3.4)。
    % 如果需要固定粒子数，可设置：
    %   ground.normalizeMode = 'fixedNorm';
    %   ground.targetNB = 3.4;

    hdt = ground.dt;
    Nt = round(ground.tmax / hdt);

    hw1 = sech(x - ground.center);
    hw1 = real(hw1(:));
    ut1 = hw1;

    Vhw = -lat.p1*cos(pi*x/lat.d1).^2 - lat.p2*cos(pi*x/lat.d2).^2;

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
            case 'fixednorm'
                hw1 = normalize_to_NB(hw1, ground.targetNB, wt);
            otherwise
                error('未知 ground.normalizeMode = %s。可选 originalPeak 或 fixedNorm。', ground.normalizeMode);
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

    info = struct();
    info.iter = it;
    info.maxErr = MaxErr;
    info.center = calc_center(x, hw1, wt);
    info.NB = calc_NB(hw1, wt);
    info.peak = max(abs(hw1));
    info.normalizeMode = ground.normalizeMode;

    fprintf('Ground soliton done: iter = %d, MaxErr = %.3e, center = %.8f, NB = %.12g, peak = %.8g\n', ...
            info.iter, info.maxErr, info.center, info.NB, info.peak);
end

%% ========================================================================
%% 函数区：Split-step 泵浦演化
%% ========================================================================
function [tdata, udata, densityData, centerX] = pump_split_step(psi0, x, k, lat, pump, wt)
    dt = pump.dt;
    nmax = round(pump.tmax/dt);
    saveStep = max(1, round(nmax/pump.nSave));
    nFrame = floor(nmax/saveStep) + 1;
    N = length(x);

    u = psi0(:);

% figure
% plot(x,udata(:,1))
    udata = zeros(N, nFrame);
    tdata = zeros(1, nFrame);
    densityData = zeros(nFrame, N);
    centerX = zeros(1, nFrame);

    frame = 1;
    udata(:,frame) = u;
    densityData(frame,:) = abs(u).'.^2;
    centerX(frame) = calc_center(x, u, wt);

    if pump.useFourthOrder
        c = 1/(2-2^(1/3));
        a = [c/2, (1-c)/2, (1-c)/2, c/2];
        b = [c, 1-2*c, c];
    else
        a = [1/2, 1/2];
        b = 1;
    end

    E = cell(numel(a),1);
    for j = 1:numel(a)
        E{j} = exp(-1i*a(j)*dt*(k.^2)/2);
    end

    for nn = 1:nmax
        tMid = (nn-0.5)*dt;
        V = external_potential(x, tMid, lat);

        if pump.useFourthOrder
            u = ifft(fft(u).*E{1});
            u = u .* exp(-1i*b(1)*dt*(V + lat.g*abs(u).^2));
            u = ifft(fft(u).*E{2});
            u = u .* exp(-1i*b(2)*dt*(V + lat.g*abs(u).^2));
            u = ifft(fft(u).*E{3});
            u = u .* exp(-1i*b(3)*dt*(V + lat.g*abs(u).^2));
            u = ifft(fft(u).*E{4});
        else
            u = ifft(fft(u).*E{1});
            u = u .* exp(-1i*dt*(V + lat.g*abs(u).^2));
            u = ifft(fft(u).*E{2});
        end

        if mod(nn, saveStep) == 0
            frame = frame + 1;
            udata(:,frame) = u;
            tdata(frame) = nn*dt;
            densityData(frame,:) = abs(u).'.^2;
            centerX(frame) = calc_center(x, u, wt);
            fprintf('pump: %8.4f / %8.4f\n', nn*dt, pump.tmax);
        end
    end
end

function xc = calc_center(x, u, wt)
    den = sum(wt.*abs(u).^2);
    xc = real(sum(wt.*x.*abs(u).^2)/(den + 1e-14));
end

%% ========================================================================
%% 函数区：作图
%% ========================================================================
function plot_all_results(x, k, bpar, Eband, xWannier, Wapsi, psi, V0, tdata, densityData, centerX, Chern, Berry, job, pump)
    %#ok<INUSD>
    cmap = soft_colormap();
    xcPlot = real(sum(x(:).*abs(psi(:)).^2)/(sum(abs(psi(:)).^2) + 1e-14));

    figure('Color','w','Name','Bands at ky=0'); hold on; box on;
    for ib = 1:bpar.Nband
        plot(bpar.kx, squeeze(Eband(:,bpar.iyWannier,ib)), 'LineWidth', 1.1);
    end
    xlabel('k_x'); ylabel('E_n(k_x,k_y=0)');
    title(sprintf('Bloch bands, selected band = %d', job.bandIndex));
    set(gca,'FontSize',14);

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

     figure('Position',[100 100 900 900])
    subplot(1,2,1);
    plot(centerX, tdata, 'LineWidth',1.5);
    xlabel('center x'); ylabel('t'); title('Soliton center'); set(gca,'FontSize',14);
    subplot(1,2,2);
    surf(x, tdata, densityData, 'EdgeColor','none');
    shading interp; view(0,90); colormap(cmap); colorbar;
    xlabel('x'); ylabel('t'); title('Pumping density |\psi(x,t)|^2');
    axis([pump.plotXLim(1), pump.plotXLim(2), 0, tdata(end)]);
    set(gca,'FontSize',14);

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
end

function bR = soft_colormap()
    Cx = linspace(0,1,1001)';
    Cy = fliplr(linspace(0,1,1001))';
    Cz = zeros(1001,1);
    Ch = ones(1001,1);
    CR = [Cz;Cz;Cx;Ch];
    CG = [Cx;Ch;Ch;Cy];
    CB = [Ch;Cy;Cz;Cz];
    bR = [CR CG CB];
    brightness_boost = 0.2;
    soften_factor = 0.7;
    bR = min(bR*soften_factor + brightness_boost, 1);
end
