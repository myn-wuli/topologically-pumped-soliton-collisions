%% vpump.m: fixed-peak vector ground/excited solitons
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
%   3. 用 Wannier 初猜和固定峰值增广 Newton 搜索激发态；基态也按固定峰值精修。
%      A=max(abs(psi))；基态、激发态默认均为 sqrt(5)，mu 和 NB 由求解得到。
%   4. 仅在泵浦阶段构造双组分矢量初态并做 split-step 演化。
%      左侧为矢量基态孤子，右侧为矢量激发态孤子。
%   5. 第一组分中两孤子的相对相位分别取 0、pi，第二组分保持 0。
%      两次演化分别保存全部工作区，最后统一绘图。

clear; clc; close all;

%% ========================== 总开关 ==========================
job.bandIndex = 2;              % 第二能带为例
job.fastMode = true;            % true: 先用较小网格跑通；正式计算可改 false
job.calcChern = false;          % Chern 对泵浦演化不是必须；需要时改 true
job.excitedCenter = 3;         % 第二能带激发态孤子中心位置；需要改位置只改这里
job.groundCenter  = -2;      % 虚时间搜索得到的基态孤子中心位置；需要改位置只改这里
job.groundWeight = 1;           % 左侧矢量基态孤子的总权重
job.excitWeight  = 1;           % 右侧矢量激发态孤子的总权重

% 矢量孤子的组分配比。
% 默认 [1;1]/sqrt(2) 表示两个组分等幅，且保持总密度与对应单组分孤子一致：
%   |chi_1 psi|^2 + |chi_2 psi|^2 = |psi|^2。
job.groundSpinor = [1; 1] / sqrt(2);     % 左侧矢量基态孤子的二组分 spinor
job.excitedSpinor = [1; 1] / sqrt(2);    % 右侧矢量激发态孤子的二组分 spinor
job.groundPhase = [0; 0];                % 左侧矢量基态孤子两个组分的相位
job.excitedPhase = [0; 0];               % 右侧矢量激发态孤子两个组分的相位

job.solitonCenter = job.excitedCenter;  % 兼容原有变量名：激发态孤子中心

%% ========================== 晶格与非线性参数 ==========================
lat.p1 = 25;
lat.p2 = 25;
lat.d1 = 1/2;
lat.d2 = 1;
lat.a  = 1;

% 单组分孤子搜索使用的非线性系数，保持原程序不变。
lat.g  = -1;

% 只在泵浦演化阶段使用的双组分非线性系数。
% 第一组分自相互作用：g11；第二组分自相互作用：g22；组分间耦合：g12。
lat.g11 = -1;
lat.g22 = -1;
lat.g12 = -1;

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
sol.t0 = 0;
phaseList = [0, pi];
phaseTags = {'0','pi'};  % 文件名相位标记

% 指定孤子包络峰值；归一化 spinor 下，矢量孤子的总峰值密度为 A^2。
sol.ampFactor = sqrt(12);
sol.solitonSearchMode = 'fixedPeak';
sol.seedWidthFactors = [1 0.7 1.4];
sol.muGuessOffsets = [0 -0.2 0.2];   % 仅改变 Newton 初猜，最终 mu 不固定
sol.maxPeakError = 1e-9;
sol.localWindow = 4.0;
sol.minLocalMass = 0.86;
sol.minNorm = 1e-6;
sol.minWannierOverlap = 0.50;
sol.maxCenterShift = 0.5*lat.a;
sol.boundaryWidth = 2*lat.a;
sol.maxBoundaryMass = 1e-6;
sol.maxResidualAccept = 1e-8;

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
pump.tmax = 7*pi/lat.nu;
pump.nSave = 1000;
pump.useFourthOrder = true;
pump.plotXLim = [-8, 8];

%% ========================== 基态孤子虚时间搜索参数 ==========================
ground.center = job.groundCenter;      % 基态孤子中心
ground.dt = 1e-3;
ground.tmax = 300;
ground.tol = 1e-11;
ground.printEvery = 1000;

% 固定峰值虚时间初猜，再用增广 Newton 同时求波函数和化学势。
ground.normalizeMode = 'fixedPeak';
ground.ampFactor = sqrt(6);
ground.seedWidth = 0.35;
ground.t0 = sol.t0;
validateattributes(sol.ampFactor,{'numeric'},{'real','finite','positive','scalar'});
validateattributes(ground.ampFactor,{'numeric'},{'real','finite','positive','scalar'});

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

%% ========================== 固定峰值激发态搜索 ==========================
V0 = external_potential(x, sol.t0, lat);
Tmat = build_kinetic_matrix(K2);
sol.R0 = Wcenter(job.bandIndex);  % 实际 Wannier 中心，包含带内位置偏移
if sol.keepSymmetry
    sol.symmetry = detect_or_set_symmetry(Wapsi,x,sol.symmetry,sol.R0,V0);
else
    sol.symmetry = 'none';
end
if sol.useEnvelope
    envelope = 1./cosh((x-sol.R0)/sol.wEnvelope);
else
    envelope = ones(size(x));
end
psi0 = normalize_to_peak(project_solution(Wapsi.*envelope,Wapsi,x,sol),sol.ampFactor);
fprintf('\n========== Stage 2: fixed-peak excited soliton ==========\n');
[psi,mu,outAll,searchInfo] = search_fixed_amplitude_soliton( ...
    Wapsi,x,Tmat,V0,wt,Eband,bpar,lat,sol);
peakAmplitude = max(abs(psi));
peakDensity = peakAmplitude^2;
assert(abs(peakAmplitude-sol.ampFactor)<=sol.maxPeakError*max(1,sol.ampFactor), ...
    '激发态峰值没有满足指定值。');
fprintf('Excited peak amplitude = %.12g; peak density = %.12g\n',peakAmplitude,peakDensity);

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

%% ========================== Stage 2b: ground-state soliton ==========================
fprintf('\n========== Stage 2b: ground-state soliton at x = %.3f ==========%s', ground.center, newline);
[psiGround, groundInfo] = search_ground_soliton_imag_time(x, k, wt, lat, ground, Tmat, sol);

% 只在泵浦演化阶段改成双组分矢量孤子初态；前面的两个孤子搜索仍保持单组分。
% 左侧：矢量基态孤子，两个组分都含有 psiGround；
% 右侧：矢量激发态孤子，两个组分都含有 psi。
%
% 形式：
%   psi_1(x,0) = cG_1 psiGround(x) + cE_1 psiExcited(x),
%   psi_2(x,0) = cG_2 psiGround(x) + cE_2 psiExcited(x).
% 默认 cG=cE=[1;1]/sqrt(2)，这样每个矢量孤子的总密度保持等于原来的单组分密度。
chiG = job.groundSpinor(:) .* exp(1i*job.groundPhase(:));
chiE = job.excitedSpinor(:) .* exp(1i*job.excitedPhase(:));
if numel(chiG) ~= 2 || numel(chiE) ~= 2
    error('job.groundSpinor 和 job.excitedSpinor 必须是二分量列向量，例如 [1;1]/sqrt(2)。');
end
if any(~isfinite(chiG)) || any(~isfinite(chiE)) || norm(chiG)==0 || norm(chiE)==0
    error('两个 spinor 均须为有限非零向量。');
end
chiG = chiG / norm(chiG);
chiE = chiE / norm(chiE);

psiExcited = psi;
psiGroundVecTotalDensity = abs(job.groundWeight*chiG(1)*psiGround).^2 + abs(job.groundWeight*chiG(2)*psiGround).^2;
psiExcitedVecTotalDensity = abs(job.excitWeight*chiE(1)*psiExcited).^2 + abs(job.excitWeight*chiE(2)*psiExcited).^2;
groundPeakAmplitude = max(abs(psiGround));
excitedPeakAmplitude = max(abs(psiExcited));
groundPeakDensity = max(psiGroundVecTotalDensity);
excitedPeakDensity = max(psiExcitedVecTotalDensity);
fprintf('Ground/excited envelope peaks = %.12g / %.12g\n', ...
    groundPeakAmplitude,excitedPeakAmplitude);
fprintf('Ground/excited vector total peak densities = %.12g / %.12g\n', ...
    groundPeakDensity,excitedPeakDensity);

%% ========================== 相位 0、pi 分别演化并保存 ==========================
nPhase = numel(phaseList);
dataFiles = cell(nPhase,1);
for iPhase = 1:nPhase
    dataFiles{iPhase} = sprintf('g12=%.6g,G=%.6g,E=%.6g,p=%s.mat', ...
        lat.g12,groundPeakDensity,excitedPeakDensity,phaseTags{iPhase});
end
for iPhase = 1:nPhase
    vphi = phaseList(iPhase);
    phaseTag = phaseTags{iPhase};
    % 每次重新组装初态；相位仅乘第一组分中的基态贡献。
    psi1_0 = job.groundWeight*chiG(1)*psiGround*exp(1i*vphi) ...
           + job.excitWeight*chiE(1)*psiExcited;
    psi2_0 = job.groundWeight*chiG(2)*psiGround ...
           + job.excitWeight*chiE(2)*psiExcited;
    fprintf('\n========== Phase %s (%d/%d) ==========\n',phaseTag,iPhase,nPhase);
    fprintf('Initial component norms = %.12g / %.12g\n', ...
        calc_NB(psi1_0,wt),calc_NB(psi2_0,wt));
    [tdata,udata1,udata2,densityData1,densityData2,centerX] = ...
        pump_split_step(psi1_0,psi2_0,x,k,lat,pump,wt);
    workspaceFile = fullfile(pwd,dataFiles{iPhase});
    save(workspaceFile,'-v7.3');  % 保存本例全部工作区，包括 vphi 和初态
    fprintf('Saved: %s\n',dataFiles{iPhase});
end

%% ========================== 全部保存后统一绘图 ==========================
% 两列分别为相位 0、pi；上排第一组分，下排第二组分。
plotData = cell(nPhase,1);
densityMax = zeros(1,2);
for iPhase = 1:nPhase
    plotData{iPhase} = load(dataFiles{iPhase}, ...
        'x','tdata','densityData1','densityData2','vphi','pump');
    plotMask = plotData{iPhase}.x>=plotData{iPhase}.pump.plotXLim(1) & ...
               plotData{iPhase}.x<=plotData{iPhase}.pump.plotXLim(2);
    if ~any(plotMask)
        error('pump.plotXLim 与计算网格没有交集。');
    end
    densityMax(1) = max(densityMax(1),max(max(plotData{iPhase}.densityData1(:,plotMask))));
    densityMax(2) = max(densityMax(2),max(max(plotData{iPhase}.densityData2(:,plotMask))));
end
densityMax = max(densityMax,eps);
figure('Color','w','Name','0, pi','Position',[80,60,1100,850]);
tiledlayout(2,nPhase,'TileSpacing','compact','Padding','compact');
for iPhase = 1:nPhase
    D = plotData{iPhase};
    for iComp = 1:2
        nexttile((iComp-1)*nPhase+iPhase);
        if iComp==1
            rhoPlot = D.densityData1;
        else
            rhoPlot = D.densityData2;
        end
        surf(D.x,D.tdata,rhoPlot,'EdgeColor','none');
        shading interp; view(2);
        set(gca,'YDir','normal','FontName','Times New Roman', ...
            'FontSize',18,'LineWidth',1.5);
        axis([D.pump.plotXLim(1),D.pump.plotXLim(2),0,D.tdata(end)]);
        caxis([0,densityMax(iComp)]); colormap(jet(256));
        cb = colorbar;
        cb.Label.String = sprintf('|\\psi_%d|^2',iComp);
        xlabel('x'); ylabel('t');
        if abs(D.vphi)<1e-12
            phaseText = '0';
        elseif abs(D.vphi-pi)<1e-12
            phaseText = '\pi';
        else
            phaseText = sprintf('%g\\pi',D.vphi/pi);
        end
        title(['\Delta\phi = ',phaseText],'FontSize',22);
    end
end

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
    % 固定峰值虚时间提供初猜，增广 Newton 同时确定波函数和化学势。
    hdt = ground.dt;
    Nt = round(ground.tmax / hdt);

    if isfield(ground,'seedWidth')
        w0 = ground.seedWidth;
    else
        w0 = 1;
    end
    hw1 = exp(-0.5*((x-ground.center)/w0).^2);
    hw1 = normalize_to_peak(real(hw1(:)),ground.ampFactor);
    ut1 = hw1;
    Vhw = external_potential(x,ground.t0,lat);

    MaxErr = inf;
    for it = 1:Nt
        hw1 = ifft(exp(-hdt/2*k.^2/2).*fft(hw1));
        hw1 = exp(-hdt*(lat.g*abs(hw1).^2 + Vhw)).*hw1;
        hw1 = ifft(exp(-hdt/2*k.^2/2).*fft(hw1));
        hw1 = real(hw1);

        hw1 = normalize_to_peak(hw1,ground.ampFactor);

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
    groundSol = sol;
    groundSol.keepSymmetry = false;
    groundSol.lockSignPattern = false;
    [hw1,muGround,refineInfo] = newton_fixed_peak(hw1,muGround,ground.ampFactor, ...
        Tmat,Vhw,wt,hw1,x,lat.g,groundSol);
    if abs(max(abs(hw1))-ground.ampFactor)>sol.maxPeakError*max(1,ground.ampFactor)
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
    info.peakDensity = info.peak^2;
    info.peakError = abs(info.peak-ground.ampFactor);
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
function [tdata, udata1, udata2, densityData1, densityData2, centerX] = pump_split_step(psi1_0, psi2_0, x, k, lat, pump, wt)
    % 双组分泵浦演化。
    % u1,u2 初态均包含左侧矢量基态孤子和右侧矢量激发态孤子。
    % 演化方程采用独立的双组分非线性系数：
    %   i d_t u1 = [-1/2 d_xx + V(x,t) + g11 |u1|^2 + g12 |u2|^2] u1,
    %   i d_t u2 = [-1/2 d_xx + V(x,t) + g22 |u2|^2 + g12 |u1|^2] u2.
    %
    % 前面的基态/激发态标量孤子搜索仍保持单组分流程；仅泵浦初态和演化是双组分。

    dt = pump.dt;
    nmax = round(pump.tmax/dt);
    saveStep = max(1, round(nmax/pump.nSave));
    nFrame = floor(nmax/saveStep) + 1;
    N = length(x);

    u1 = psi1_0(:);
    u2 = psi2_0(:);

    udata1 = zeros(N, nFrame);
    udata2 = zeros(N, nFrame);
    tdata = zeros(1, nFrame);
    densityData1 = zeros(nFrame, N);
    densityData2 = zeros(nFrame, N);

    % centerX(1,:)：第一组分中心；centerX(2,:)：第二组分中心；centerX(3,:)：总密度中心。
    centerX = zeros(3, nFrame);

    frame = 1;
    udata1(:,frame) = u1;
    udata2(:,frame) = u2;
    densityData1(frame,:) = abs(u1).'.^2;
    densityData2(frame,:) = abs(u2).'.^2;
    centerX(1,frame) = calc_center(x, u1, wt);
    centerX(2,frame) = calc_center(x, u2, wt);
    centerX(3,frame) = calc_center_total(x, u1, u2, wt);

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
            % kinetic a1
            u1 = ifft(fft(u1).*E{1});
            u2 = ifft(fft(u2).*E{1});

            % nonlinear b1，两个组分分别使用 g11/g22/g12 相位，必须同步更新
            rho1 = abs(u1).^2;
            rho2 = abs(u2).^2;
            phase1 = exp(-1i*b(1)*dt*(V + lat.g11*rho1 + lat.g12*rho2));
            phase2 = exp(-1i*b(1)*dt*(V + lat.g22*rho2 + lat.g12*rho1));
            u1 = u1 .* phase1;
            u2 = u2 .* phase2;

            % kinetic a2
            u1 = ifft(fft(u1).*E{2});
            u2 = ifft(fft(u2).*E{2});

            % nonlinear b2
            rho1 = abs(u1).^2;
            rho2 = abs(u2).^2;
            phase1 = exp(-1i*b(2)*dt*(V + lat.g11*rho1 + lat.g12*rho2));
            phase2 = exp(-1i*b(2)*dt*(V + lat.g22*rho2 + lat.g12*rho1));
            u1 = u1 .* phase1;
            u2 = u2 .* phase2;

            % kinetic a3
            u1 = ifft(fft(u1).*E{3});
            u2 = ifft(fft(u2).*E{3});

            % nonlinear b3
            rho1 = abs(u1).^2;
            rho2 = abs(u2).^2;
            phase1 = exp(-1i*b(3)*dt*(V + lat.g11*rho1 + lat.g12*rho2));
            phase2 = exp(-1i*b(3)*dt*(V + lat.g22*rho2 + lat.g12*rho1));
            u1 = u1 .* phase1;
            u2 = u2 .* phase2;

            % kinetic a4
            u1 = ifft(fft(u1).*E{4});
            u2 = ifft(fft(u2).*E{4});
        else
            % 二阶 split-step
            u1 = ifft(fft(u1).*E{1});
            u2 = ifft(fft(u2).*E{1});

            rho1 = abs(u1).^2;
            rho2 = abs(u2).^2;
            phase1 = exp(-1i*dt*(V + lat.g11*rho1 + lat.g12*rho2));
            phase2 = exp(-1i*dt*(V + lat.g22*rho2 + lat.g12*rho1));
            u1 = u1 .* phase1;
            u2 = u2 .* phase2;

            u1 = ifft(fft(u1).*E{2});
            u2 = ifft(fft(u2).*E{2});
        end

        if mod(nn, saveStep) == 0
            frame = frame + 1;
            udata1(:,frame) = u1;
            udata2(:,frame) = u2;
            tdata(frame) = nn*dt;
            densityData1(frame,:) = abs(u1).'.^2;
            densityData2(frame,:) = abs(u2).'.^2;
            centerX(1,frame) = calc_center(x, u1, wt);
            centerX(2,frame) = calc_center(x, u2, wt);
            centerX(3,frame) = calc_center_total(x, u1, u2, wt);
            fprintf('pump: %8.4f / %8.4f\n', nn*dt, pump.tmax);
        end
    end
end

function xc = calc_center(x, u, wt)
    den = sum(wt.*abs(u).^2);
    xc = real(sum(wt.*x.*abs(u).^2)/(den + 1e-14));
end

function xc = calc_center_total(x, u1, u2, wt)
    den = sum(wt.*(abs(u1).^2 + abs(u2).^2));
    xc = real(sum(wt.*x.*(abs(u1).^2 + abs(u2).^2))/(den + 1e-14));
end

function xc = calc_center_density(x, rho, wt)
    den = sum(wt.*rho);
    xc = real(sum(wt.*x.*rho)/(den + 1e-14));
end
