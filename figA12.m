function figA12(mode, caseFile)
% 粗筛+临界区细扫，轻量运行，仅保存最终相图pump_displacement_scan.png。
% 在adaptive_scan_settings函数中调整粗扫点、细扫步长和位移变化阈值。
% 运行：figA12
% 读取旧版已有结果并完成剩余扫描：figA12('resume')
% 指定旧版目录：figA12('resume','D:/your_folder/GroundAmplitude_PumpScan')
% 新扫描不保存断点、逐组MAT、CSV、密度图或完整时空密度矩阵。
% 演化时仅暂存当前组的波包特征，分析后立即释放；各组只累计少量最终指标。
% 合成轨迹自检：figA12('selftest')
% q=[X(t0+T)-X(t0)]/a；右向为正，左向为负；无法辨识记为NaN。
% packetID按测量起点从左到右编号，不强行对应入射孤子。
if nargin>0
    if strcmpi(mode,'selftest')
        scan_selftest(); return;
    elseif strcmpi(mode,'resume')
        if nargin<2, caseFile=fullfile(pwd,'GroundAmplitude_PumpScan'); end
        resume_saved_scan(caseFile); return;
    elseif strcmpi(mode,'reanalyze')
        if nargin<2, error('请提供case MAT路径。'); end
        if nargin<3, collisionEndT=[]; end
        reanalyze_saved_case(caseFile,collisionEndT); return;
    else
        error('mode可省略，或使用selftest/resume。');
    end
end
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

clc; close all;

%% ========================== 总开关 ==========================
job.bandIndex = 2;              % 第二能带为例
job.fastMode = true;            % true: 先用较小网格跑通；正式计算可改 false
job.calcChern = false;          % Chern 对泵浦演化不是必须；需要时改 true
job.savePrefix = 'SecondBand_ExcitedSoliton_x4_GroundSoliton_xm3p5_Pump';
job.excitedCenter = 3;         % 第二能带激发态孤子中心位置
job.groundCenter  = -2;      % 虚时间搜索得到的基态孤子中心位置
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
pump.nSave = 1000;
pump.useFourthOrder = true;
pump.plotXLim = [-8, 8];

%% ========================== 基态孤子虚时间搜索参数 ==========================
ground.center = job.groundCenter;      % 基态孤子中心：由job.groundCenter指定
ground.dt = 1e-3;
ground.tmax = 300;
ground.tol = 1e-11;
ground.printEvery = 1000;

% 保留你原程序的归一化方式：
% hw1 = hw1 * sqrt(af/max(abs(hw1))), af = sqrt(3.4)。
% 若后续想固定粒子数，可把 normalizeMode 改为 'fixedNorm'，并设置 targetNB。
ground.normalizeMode = 'originalPeak';
ground.ampFactor = sqrt(2.5);
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

%% ========================== 扫描设置（主要修改这里） ==========================
scan = adaptive_scan_settings(); % 粗扫21点，自动选择变化区间，局部细扫
pump.nPeriods = 5;            % 若碰撞后不足一周期，可增加这里
pump.framesPerPeriod = 200;   % 每周期等间隔保存，确保比较时间严格相差T
pump.tmax = pump.nPeriods*pi/lat.nu;

% 下列均为可调的波包识别阈值，不是对“拓扑孤子”的严格数学认证。
trackOpt.quantTolerance = 0.05;   % |Delta x/a - 1| <= 0.05
trackOpt.smoothSigma = 0.35*lat.a; % 仅检测包络峰；质心使用未平滑的密度
trackOpt.minPeakDistance = 0.9*lat.a;
trackOpt.coreRadius = 0.9*lat.a;
trackOpt.maxWidth = 0.65*lat.a;    % 局域密度的均方根宽度上限
trackOpt.minConcentration = 0.65; % N(|x-xc|<R)/N(|x-xc|<2R)
trackOpt.minMassRelativeSmall = 0.10; % 相对于较弱初始孤子的粒子数
trackOpt.minMassRelativeTotal = 1e-4;
trackOpt.minPeakRelative = 1e-5;
trackOpt.isolationDistance = 1.8*lat.a;
trackOpt.contactDistance = 1.8*lat.a;
trackOpt.minApproach = 0.3*lat.a;
trackOpt.maxJump = 0.35*lat.a;     % 相邻保存帧的最大允许位移
trackOpt.maxMassVariation = 0.30; % 整个跟踪区间局域粒子数相对首帧变化上限
trackOpt.maxWidthRatio = 1.8;
trackOpt.edgeMargin = 2.5*lat.a;
trackOpt.edgeMassTolerance = 1e-3;
trackOpt.normTolerance = 1e-7;
trackOpt.relaxPeriods = 0.25;     % 首次接触后至少等待0.25T
trackOpt.holdPeriods = 0.15;      % 连续稳定0.15T后开始测量
% 手动指定碰撞结束时刻（单位T）：[]自动；标量所有组共用。
% 也可设函数句柄，如 @(ampSquared) 3.2，供新增的细扫点使用。
% 可根据已知碰撞结束时刻手动设置；轻量版不保存中间密度图。
trackOpt.collisionEndT = [];

if ~strcmpi(ground.normalizeMode, 'originalPeak')
    error('振幅扫描需要 ground.normalizeMode = originalPeak；fixedNorm不会使用ampFactor。');
end
if job.excitWeight == 0 || job.groundWeight == 0
    error('碰撞扫描要求两种孤子的权重均非零。');
end
if max(abs(R)) > sol.maxResidualAccept || locFinal.massFraction < sol.minLocalMass
    error('激发态仅为未通过阈值的诊断态，请先得到合格的局域定态再扫描。');
end
if ~isempty(trackOpt.collisionEndT) && ~isa(trackOpt.collisionEndT,'function_handle') && ...
        ~isscalar(trackOpt.collisionEndT)
    error('粗细扫描的collisionEndT应为空、标量或关于ampSquared的函数句柄。');
end
if ~exist(scan.outputDir, 'dir'), mkdir(scan.outputDir); end
psiExcitedFixed = psi;        % 在循环外只计算一次，不随基态振幅改变
% 定态求解结束后释放不再使用的能带/Wannier/Newton中间量。
clear Eband Coeff Berry Chern WannierAll xWannier Wapsi psi0 outAll searchInfo Tmat V0 R psi
run_amplitude_scan(psiExcitedFixed,x,k,wt,lat,ground,pump,scan,trackOpt,job,vphi,false);
end % 主函数结束；以下局部函数无需另存

function resume_saved_scan(folder)
    commonFile=fullfile(folder,'common_initial_state.mat');
    if ~exist(commonFile,'file')
        error('找不到%s，请指定原扫描结果目录。',commonFile);
    end
    required={'psiExcitedFixed','x','k','wt','lat','ground','pump','scan','trackOpt','job','vphi'};
    C=load(commonFile,required{:});
    for j=1:numel(required)
        if ~isfield(C,required{j}), error('保存文件缺少变量%s。',required{j}); end
    end
    oldScan=C.scan;
    C.scan=adaptive_scan_settings();
    C.scan.outputDir=folder;
    C.scan.legacyAmpSquared=oldScan.ampSquared;
    % 旧版按固定网格逐点指定的碰撞时间，只在原网格范围内插值。
    if ~isempty(C.trackOpt.collisionEndT) && ~isa(C.trackOpt.collisionEndT,'function_handle') && ...
            ~isscalar(C.trackOpt.collisionEndT)
        oldA2=oldScan.ampSquared(:); oldTimes=C.trackOpt.collisionEndT(:);
        assert(numel(oldA2)==numel(oldTimes),'旧版碰撞时间与扫描网格长度不一致。');
        C.trackOpt.collisionEndT=@(a2) interp1(oldA2,oldTimes,a2,'linear',NaN);
    end
    fprintf('读取旧初态和已有结果，按新的粗筛/细扫网格完成扫描。\n');
    run_amplitude_scan(C.psiExcitedFixed,C.x,C.k,C.wt,C.lat,C.ground, ...
        C.pump,C.scan,C.trackOpt,C.job,C.vphi,true);
end

function scan = adaptive_scan_settings()
    % 所有步长均针对A2=ground.ampFactor^2，不是ampFactor本身。
    scan.ampSquared = [0.1,1:1:20]; % 粗筛：21点，覆盖sqrt(0.1)至sqrt(20)
    scan.ampFactor = sqrt(scan.ampSquared);
    scan.fineStep = 0.1;           % 变化区间内最大A2步长；需要时改为0.05/0.02
    scan.displacementJump = 0.25;  % 相邻点Q变化超过此值时也加密
    scan.extraFineRanges = [];    % 可手动补充区间，如[4 5; 8 8.5]
    scan.outputDir = fullfile(pwd,'GroundAmplitude_PumpScan');
    scan.legacyAmpSquared = [];   % 仅resume时用来匹配旧版文件
end

function run_amplitude_scan(psiExcitedFixed,x,k,wt,lat,ground,pump,scan,trackOpt,job,vphi,resumeMode)
coarseA2=unique(scan.ampSquared(:).');
assert(~isempty(coarseA2) && all(isfinite(coarseA2)) && all(coarseA2>0),'粗扫点必须为正数。');
assert(isscalar(scan.fineStep) && isfinite(scan.fineStep) && scan.fineStep>0,'fineStep必须为正数。');
assert(isscalar(scan.displacementJump) && isfinite(scan.displacementJump) && ...
    scan.displacementJump>0,'displacementJump必须为正数。');
allRows=repmat(empty_scan_row(),0,1);
coarseRows=cell(numel(coarseA2),1);
fprintf('\n===== 第一阶段：粗筛，共%d点 =====\n',numel(coarseA2));
for j=1:numel(coarseA2)
    rows=evaluate_amplitude(coarseA2(j),j,psiExcitedFixed,x,k,wt,lat,ground, ...
        pump,scan,trackOpt,job,vphi,resumeMode);
    coarseRows{j}=rows;
    allRows=[allRows;rows(:)]; %#ok<AGROW>
end
ranges=select_fine_ranges(coarseA2,coarseRows,scan);
fineA2=make_fine_points(ranges,coarseA2,scan.fineStep);
fprintf('\n===== 第二阶段：%d个待加密区间，新增%d点 =====\n',size(ranges,1),numel(fineA2));
for j=1:size(ranges,1)
    fprintf('  A2区间 [%.6g, %.6g]\n',ranges(j,1),ranges(j,2));
end
for j=1:numel(fineA2)
    rows=evaluate_amplitude(fineA2(j),numel(coarseA2)+j,psiExcitedFixed,x,k,wt,lat, ...
        ground,pump,scan,trackOpt,job,vphi,resumeMode);
    allRows=[allRows;rows(:)]; %#ok<AGROW>
end
resultTable=struct2table(allRows,'AsArray',true);
resultTable=sortrows(resultTable,{'ampSquared','packetID'});
plot_scan_summary(resultTable,scan.outputDir,trackOpt.quantTolerance);
fprintf('粗筛%d点，细扫新增%d点，总计%d点。\n',numel(coarseA2),numel(fineA2),numel(coarseA2)+numel(fineA2));
fprintf('仅保存最终相图：%s\n',fullfile(scan.outputDir,'pump_displacement_scan.png'));
end

function ranges=select_fine_ranges(a2,rowsByPoint,scan)
    % 候选变化区间，不将识别失败边界直接解释为物理临界点。
    % 不同的出射孤子不作跨参数身份匹配；比较各类波包数量及排序后的Q。
    ranges=zeros(0,2);
    for j=1:numel(a2)-1
        r1=rowsByPoint{j}; r2=rowsByPoint{j+1};
        sig1=packet_phase_signature(r1); sig2=packet_phase_signature(r2);
        changed=any(sig1~=sig2);
        q1=sort([r1([r1.valid]).deltaXOverA]);
        q2=sort([r2([r2.valid]).deltaXOverA]);
        jump=~isempty(q1) && numel(q1)==numel(q2) && ...
            any(abs(q1-q2)>scan.displacementJump);
        if changed || jump
            ranges(end+1,:)=[a2(j),a2(j+1)]; %#ok<AGROW>
        end
    end
    if ~isempty(scan.extraFineRanges)
        extra=scan.extraFineRanges;
        assert(size(extra,2)==2 && all(isfinite(extra(:))) && ...
            all(extra(:,2)>extra(:,1)),'extraFineRanges每行必须是递增的两个有限数。');
        assert(all(extra(:,1)>=min(a2)) && all(extra(:,2)<=max(a2)), ...
            '额外细扫区间必须位于粗扫范围内。');
        ranges=[ranges;extra];
    end
    ranges=unique(ranges,'rows');
end

function sig=packet_phase_signature(rows)
    valid=[rows.valid];
    sig=[nnz(valid & [rows.isPlusOne]==1), ...
         nnz(valid & [rows.isMinusOne]==1), ...
         nnz(valid & [rows.isUnitQuantized]==0),nnz(~valid)];
end

function fineA2=make_fine_points(ranges,coarseA2,step)
    fineA2=[];
    for j=1:size(ranges,1)
        % 用linspace覆盖端点，避免最后一段超过指定步长。
        segments=max(1,ceil((ranges(j,2)-ranges(j,1))/step-1e-12));
        fineA2=[fineA2,linspace(ranges(j,1),ranges(j,2),segments+1)]; %#ok<AGROW>
    end
    fineA2=sort(fineA2);
    if isempty(fineA2), return; end
    fineA2=fineA2([true,diff(fineA2)>1e-10]);
    keep=true(size(fineA2));
    for j=1:numel(fineA2)
        keep(j)=all(abs(coarseA2-fineA2(j))>1e-10);
    end
    fineA2=fineA2(keep); % 粗扫点/重叠区间不重复求解
end

function rows=evaluate_amplitude(a2,iEval,psiExcitedFixed,x,k,wt,lat,ground,pump,scan,trackOpt,job,vphi,resumeMode)
    % 每次调用结束即释放当前点的波函数、波包特征和诊断数组。
    if resumeMode && ~isempty(scan.legacyAmpSquared)
        oldIndex=find(abs(scan.legacyAmpSquared-a2)<1e-10,1);
        if ~isempty(oldIndex)
            stem=sprintf('case_%03d_A2_%07.3f',oldIndex,scan.legacyAmpSquared(oldIndex));
            casePath=fullfile(scan.outputDir,[stem '.mat']);
            if exist(casePath,'file')
                C=load(casePath,'rows','ground');
                if isfield(C,'ground') && abs(C.ground.ampFactor-sqrt(a2))<1e-10 && ...
                        isfield(C,'rows') && ~isempty(C.rows) && ...
                        ~any(strcmp({C.rows.status},'computation_failed'))
                    rows=C.rows;
                    for j=1:numel(rows), rows(j).scanIndex=iEval; end
                    fprintf('复用已有结果：A2=%.6g\n',a2); return;
                end
            end
        end
    end
    ground.ampFactor=sqrt(a2);
    fprintf('\n===== 计算点%d：ampFactor=sqrt(%.6g) =====\n',iEval,a2);
    rowBase=empty_scan_row();
    rowBase.scanIndex=iEval; rowBase.ampSquared=a2; rowBase.ampFactor=sqrt(a2);
    rowBase.excitedNB=calc_NB(job.excitWeight*psiExcitedFixed,wt);
    optNow=trackOpt;
    try
        if isa(trackOpt.collisionEndT,'function_handle')
            optNow.collisionEndT=trackOpt.collisionEndT(a2);
        end
        if ~isempty(optNow.collisionEndT)
            assert(isscalar(optNow.collisionEndT) && isfinite(optNow.collisionEndT) && ...
                optNow.collisionEndT>=0,'指定的碰撞结束时间必须为非负有限标量。');
        end
        [psiGround,groundInfo]=search_ground_soliton_imag_time(x,k,wt,lat,ground);
        rowBase.groundNB=calc_NB(job.groundWeight*psiGround,wt);
        rowBase.groundPeak=max(abs(psiGround));
        rowBase.groundIterationError=groundInfo.maxErr;
        if groundInfo.maxErr>ground.tol
            error('基态虚时间迭代尚未收敛：maxErr=%.3e。',groundInfo.maxErr);
        end
        psiPump0=job.excitWeight*psiExcitedFixed+job.groundWeight*psiGround*exp(1i*vphi);
        optNow.minMass=max(trackOpt.minMassRelativeSmall*min(rowBase.excitedNB,rowBase.groundNB), ...
            trackOpt.minMassRelativeTotal*calc_NB(psiPump0,wt));
        optNow.initialCenters=[calc_center(x,psiGround,wt),calc_center(x,psiExcitedFixed,wt)];
        [tdata,packetData,diagPump]=pump_scan_exact_period(psiPump0,x,k,lat,pump,wt,optNow);
        rows=analyze_scan_case(tdata,packetData,x,lat,pump.framesPerPeriod,optNow,rowBase,diagPump);
    catch ME
        rows=rowBase; rows.status='computation_failed'; rows.note=ME.message;
        warning('A2=%.6g计算失败：%s',a2,ME.message);
    end
    fprintf('A2=%.6g完成；有效波包%d个。\n',a2,nnz([rows.valid]));
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
    groundKinetic = exp(-hdt*k.^2/4);

    hw1 = sech(x - ground.center);
    hw1 = real(hw1(:));
    ut1 = hw1;

    Vhw = -lat.p1*cos(pi*x/lat.d1).^2 - lat.p2*cos(pi*x/lat.d2).^2;

    MaxErr = inf;
    for it = 1:Nt
        hw1 = ifft(groundKinetic.*fft(hw1));
        hw1 = exp(-hdt*(lat.g*abs(hw1).^2 + Vhw)).*hw1;
        hw1 = ifft(groundKinetic.*fft(hw1));
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
function xc = calc_center(x, u, wt)
    den = sum(wt.*abs(u).^2);
    xc = real(sum(wt.*x.*abs(u).^2)/(den + 1e-14));
end

%% ========================================================================
%% 函数区：作图
%% ========================================================================
%% ========================== 振幅扫描专用函数 ==========================
function r = empty_scan_row()
    r = struct('scanIndex',NaN,'ampSquared',NaN,'ampFactor',NaN, ...
        'groundPeak',NaN,'groundNB',NaN,'excitedNB',NaN, ...
        'groundIterationError',NaN,'packetID',NaN,'contactTimeOverT',NaN, ...
        'tStartOverT',NaN,'tEndOverT',NaN,'xStart',NaN,'xEnd',NaN, ...
        'deltaX',NaN,'deltaXOverA',NaN,'errorToPlusOne',NaN, ...
        'errorToUnitMagnitude',NaN,'isPlusOne',NaN,'isMinusOne',NaN, ...
        'isUnitQuantized',NaN,'massStart',NaN,'massEnd',NaN, ...
        'maxMassVariation',NaN,'maxNormRelativeError',NaN, ...
        'maxEdgeMassFraction',NaN,'valid',false,'status','','note','');
end

function [tdata,packetData,dg] = pump_scan_exact_period(u,x,k,lat,pump,wt,opt)
    % 保留原实时间推进算法，仅使分析采样帧严格对齐周期。
    % 原useFourthOrder把V冻结在整步中点；对含时势不宣称整体四阶精度。
    assert(lat.nu>0 && pump.dt>0,'nu和dt必须为正。');
    assert(pump.nPeriods>=1 && mod(pump.nPeriods,1)==0,'nPeriods必须为正整数。');
    assert(pump.framesPerPeriod>=2 && mod(pump.framesPerPeriod,1)==0,'framesPerPeriod必须为至少2的整数。');
    T = pi/lat.nu;
    nStepsPerFrame = max(1,ceil(T/(pump.dt*pump.framesPerPeriod)));
    nStepsPerPeriod = nStepsPerFrame*pump.framesPerPeriod;
    dt = T/nStepsPerPeriod; % 与原dt非常接近，且不超过原dt
    nFrames = pump.nPeriods*pump.framesPerPeriod+1;
    nSteps = (nFrames-1)*nStepsPerFrame;
    tdata = (0:nFrames-1)*(T/pump.framesPerPeriod);
    packetData = cell(nFrames,1);
    packetData{1} = detect_packets(abs(u(:)).^2,x,opt);
    norm0 = calc_NB(u,wt);
    normTrace = zeros(nFrames,1); normTrace(1) = norm0;
    edge = x < x(1)+opt.edgeMargin | x > x(end)-opt.edgeMargin;
    edgeTrace = zeros(nFrames,1);
    edgeTrace(1) = sum(wt(edge).*abs(u(edge)).^2)/norm0;
    if pump.useFourthOrder
        c = 1/(2-2^(1/3));
        a = [c/2,(1-c)/2,(1-c)/2,c/2]; b = [c,1-2*c,c];
    else
        a = [1/2,1/2]; b = 1;
    end
    E = exp(-1i*(k(:).^2)*a*dt/2);
    frame = 1;
    for nn = 1:nSteps
        V = external_potential(x,(nn-0.5)*dt,lat);
        u = ifft(fft(u).*E(:,1));
        for j = 1:numel(b)
            u = u.*exp(-1i*b(j)*dt*(V+lat.g*abs(u).^2));
            u = ifft(fft(u).*E(:,j+1));
        end
        if mod(nn,nStepsPerFrame) == 0
            frame = frame+1;
            rhoNow = abs(u(:)).^2;
            normTrace(frame) = calc_NB(u,wt);
            edgeTrace(frame) = sum(wt(edge).*abs(u(edge)).^2)/normTrace(frame);
            if any(~isfinite(rhoNow))
                error('实时间演化出现NaN/Inf。');
            end
            packetData{frame} = detect_packets(rhoNow,x,opt);
            if mod(frame-1,pump.framesPerPeriod) == 0
                fprintf('  pump: %d/%d periods\n',(frame-1)/pump.framesPerPeriod,pump.nPeriods);
            end
        end
    end
    dg = struct('dtActual',dt,'T',T,'normTrace',normTrace, ...
        'normRelativeError',abs(normTrace/norm0-1), ...
        'edgeMassFraction',edgeTrace);
end

function p = detect_packets(rho,x,opt)
    % p每行：[局域质心,局域粒子数,RMS宽度,局域质量比,是否可辨识]。
    % 高斯平滑仅消除晶胞内多峰；质量/质心/宽度均由原始密度计算。
    rho = max(real(rho(:)),0); x = x(:); dx = x(2)-x(1);
    n = ceil(4*opt.smoothSigma/dx);
    ker = exp(-0.5*((-n:n)'*dx/opt.smoothSigma).^2);
    sm = conv(rho,ker/sum(ker),'same');
    idx = find(sm(2:end-1)>=sm(1:end-2) & sm(2:end-1)>sm(3:end))+1;
    idx = idx(sm(idx)>=opt.minPeakRelative*max(sm));
    [~,order] = sort(sm(idx),'descend'); idx = idx(order);
    kept = [];
    for j = 1:numel(idx)
        if isempty(kept) || all(abs(x(idx(j))-x(kept))>=opt.minPeakDistance)
            kept(end+1) = idx(j); %#ok<AGROW>
        end
    end
    p = zeros(0,5);
    for j = 1:numel(kept)
        xc = x(kept(j));
        for it = 1:3
            core = abs(x-xc)<=opt.coreRadius;
            mass = dx*sum(rho(core));
            if mass<=eps, break; end
            xc = dx*sum(x(core).*rho(core))/mass;
        end
        core = abs(x-xc)<=opt.coreRadius;
        mass = dx*sum(rho(core));
        if mass<opt.minMass, continue; end
        xc = dx*sum(x(core).*rho(core))/mass;
        width = sqrt(dx*sum((x(core)-xc).^2.*rho(core))/mass);
        outer = abs(x-xc)<=2*opt.coreRadius;
        concentration = mass/(dx*sum(rho(outer))+eps);
        good = width<=opt.maxWidth && concentration>=opt.minConcentration && ...
            xc>x(1)+opt.edgeMargin && xc<x(end)-opt.edgeMargin;
        % 两个原始峰收敛到同一个包络质心时，仅保留一个候选。
        if ~isempty(p) && any(abs(p(:,1)-xc)<opt.minPeakDistance/2), continue; end
        p(end+1,:) = [xc,mass,width,concentration,good]; %#ok<AGROW>
    end
    if ~isempty(p), p = sortrows(p,1); end
    for j = 1:size(p,1)
        others = [1:j-1,j+1:size(p,1)];
        if any(abs(p(others,1)-p(j,1))<opt.isolationDistance), p(j,5)=0; end
    end
end

function [ic,reason] = find_contact(packets,opt)
    % 只在入射两波包已明显靠近的情况下标记接触；远处消失不算碰撞。
    ic = NaN; reason = 'initial_packets_unresolved';
    prev = sort(opt.initialCenters(:));
    p = packets{1};
    if size(p,1)<2, return; end
    pick = zeros(2,1);
    for j = 1:2
        [dist,pick(j)] = min(abs(p(:,1)-prev(j)));
        if dist>opt.coreRadius, return; end
    end
    if pick(1)==pick(2), return; end
    prev = p(pick,1); d0 = abs(diff(prev));
    if d0<=opt.contactDistance, return; end
    reason = 'no_collision_detected';
    for it = 2:numel(packets)
        p = packets{it}; now = nan(2,1); pick = zeros(2,1);
        for j = 1:2
            ids = find(abs(p(:,1)-prev(j))<=opt.maxJump);
            if numel(ids)==1, pick(j)=ids; now(j)=p(ids,1); end
        end
        if any(~isfinite(now)) || pick(1)==pick(2)
            if abs(diff(prev))<=opt.contactDistance+2*opt.maxJump && ...
                    d0-abs(diff(prev))>=opt.minApproach
                ic=it; reason='approach_then_unresolved_contact';
            else
                reason='incoming_track_lost_before_contact';
            end
            return;
        end
        sep = abs(diff(now));
        if sep<=opt.contactDistance && d0-sep>=opt.minApproach
            ic=it; reason='approaching_packets_contact'; return;
        end
        prev=now;
    end
end

function [path,ok,reason] = follow_packet(packets,startFrame,endFrame,id,opt)
    path = nan(endFrame-startFrame+1,4);
    first = packets{startFrame}; p = first(id,:);
    ok=false; reason='packet_not_localized_or_isolated';
    if ~p(5), return; end
    path(1,:) = p(1:4); mass0=p(2); width0=max(p(3),eps);
    for it = startFrame+1:endFrame
        q = packets{it}; ids=find(abs(q(:,1)-p(1))<=opt.maxJump);
        if isempty(ids), reason='packet_track_lost'; return; end
        if numel(ids)>1, reason='ambiguous_packet_assignment'; return; end
        p=q(ids,:);
        if ~p(5), reason='packet_lost_localization_or_isolation'; return; end
        path(it-startFrame+1,:)=p(1:4);
        if abs(p(2)/mass0-1)>opt.maxMassVariation
            reason='local_mass_not_stable'; return;
        end
        if p(3)/width0>opt.maxWidthRatio || width0/max(p(3),eps)>opt.maxWidthRatio
            reason='packet_width_not_stable'; return;
        end
    end
    ok=true; reason='ok';
end

function [rows,out] = analyze_scan_case(t,rho,x,lat,fpp,opt,base,dg)
    T=pi/lat.nu; nt=numel(t);
    if iscell(rho)
        % 轻量演化入口：已在线提取波包特征，没有完整密度矩阵。
        packets=rho;
        assert(numel(packets)==nt,'波包特征帧数与时间采样不一致。');
    else
        % 仅供小规模合成轨迹自检。
        packets=cell(nt,1);
        for it=1:nt, packets{it}=detect_packets(rho(it,:),x,opt); end
    end
    out=struct('packets',{packets},'contactFrame',NaN,'startFrame',NaN, ...
        'paths',{{}},'pathFrames',{{}},'contactReason','','startReason','');
    rows=base;
    if isempty(opt.collisionEndT)
        [ic,why]=find_contact(packets,opt);
        out.contactReason=why;
        if isnan(ic), rows.status=why; return; end
        out.contactFrame=ic; base.contactTimeOverT=t(ic)/T;
        earliest=ic+ceil(opt.relaxPeriods*fpp);
        out.startReason='automatic_localization_and_isolation_heuristic';
    else
        % 手动值代表已确认碰撞结束，而非首次接触时刻。
        earliest=find(t/T>=opt.collisionEndT-1e-12,1);
        if isempty(earliest), rows.status='manual_time_outside_simulation'; return; end
        out.contactReason='manual_collision_end';
        out.startReason='manual_end_then_stability_hold';
    end
    holdFrames=max(1,ceil(opt.holdPeriods*fpp));
    start=NaN;
    % 窗口选择只依据局域/分离/质量稳定，不使用最终位移是否接近1。
    for trial=earliest:nt-fpp-holdFrames
        p=packets{trial};
        for id=1:size(p,1)
            % 自动模式要求整个待测周期也连续可辨识，避免在碰撞的暂时合峰阶段启动。
            checkEnd=trial+holdFrames;
            if isempty(opt.collisionEndT), checkEnd=checkEnd+fpp; end
            [~,ok]=follow_packet(packets,trial,checkEnd,id,opt);
            if ok, start=trial+holdFrames; break; end
        end
        if isfinite(start), break; end
    end
    rows=base;
    if isnan(start)
        rows.status='no_stable_postcollision_window';
        rows.note='未找到可辨识的碰撞后波包或剩余时间不足T；查看密度图，必要时延长演化/手动指定结束时间。';
        return;
    end
    out.startFrame=start; finish=start+fpp;
    if abs((t(finish)-t(start))/T-1)>1e-10
        error('分析窗口长度必须严格等于驱动周期T。');
    end
    p=packets{start}; ids=find(p(:,5)>0);
    rows=repmat(base,numel(ids),1);
    for j=1:numel(ids)
        [path,ok,why]=follow_packet(packets,start,finish,ids(j),opt);
        out.paths{j}=path; out.pathFrames{j}=start:finish;
        r=base; r.packetID=j; r.tStartOverT=t(start)/T; r.tEndOverT=t(finish)/T;
        r.maxNormRelativeError=max(dg.normRelativeError(1:finish));
        r.maxEdgeMassFraction=max(dg.edgeMassFraction(1:finish));
        if r.maxNormRelativeError>opt.normTolerance, ok=false; why='norm_drift'; end
        if r.maxEdgeMassFraction>opt.edgeMassTolerance, ok=false; why='boundary_contamination'; end
        r.valid=ok; r.status=why;
        if ok
            r.xStart=path(1,1); r.xEnd=path(end,1);
            r.deltaX=r.xEnd-r.xStart; r.deltaXOverA=r.deltaX/lat.a;
            r.massStart=path(1,2); r.massEnd=path(end,2);
            r.maxMassVariation=max(abs(path(:,2)/path(1,2)-1));
            r.errorToPlusOne=abs(r.deltaXOverA-1);
            r.errorToUnitMagnitude=abs(abs(r.deltaXOverA)-1);
            r.isPlusOne=double(r.errorToPlusOne<=opt.quantTolerance);
            r.isMinusOne=double(abs(r.deltaXOverA+1)<=opt.quantTolerance);
            r.isUnitQuantized=double(r.errorToUnitMagnitude<=opt.quantTolerance);
            if r.isPlusOne, r.status='right_unit_quantized';
            elseif r.isMinusOne, r.status='left_unit_quantized';
            else, r.status='resolved_nonunit_displacement'; end
        else
            r.note='位移及量子化标记为NaN；此行不能解释为零位移或不量子化。';
        end
        rows(j)=r;
    end
end

function plot_scan_summary(tbl,outDir,tol)
    fig=figure('Color','w'); hold on; box on;
    valid=tbl.valid; unit=valid & tbl.isUnitQuantized==1;
    right=valid & tbl.isPlusOne==1; left=valid & tbl.isMinusOne==1;
    h1=scatter(tbl.ampSquared(right),tbl.deltaXOverA(right),30,[0 .4 .85],'filled');
    h2=scatter(tbl.ampSquared(left),tbl.deltaXOverA(left),30,[.85 .25 .1],'filled');
    h3=scatter(tbl.ampSquared(valid & ~unit),tbl.deltaXOverA(valid & ~unit),22,'k','filled');
    legend([h1,h2,h3],{'Right: +1','Left: -1','Other displacement'},'Location','best');
    xl=[min(tbl.ampSquared),max(tbl.ampSquared)];
    if diff(xl)==0, xl=xl+[-.05,.05]; end
    for q=[-1,1]
        plot(xl,[q,q],'--','Color',[.4 .4 .4]);
        plot(xl,[q-tol,q-tol],':','Color',[.7 .7 .7]);
        plot(xl,[q+tol,q+tol],':','Color',[.7 .7 .7]);
    end
    xlabel('ground.ampFactor^2'); ylabel('\Delta x/a'); xlim(xl);
    title('Post-collision displacement over one drive period');
    set(gca,'FontName','Times New Roman','FontSize',13);
    print(fig,fullfile(outDir,'pump_displacement_scan.png'),'-dpng','-r300');
    fprintf('有效波包行数：%d；接近+1：%d；接近-1：%d；无效行数：%d。\n', ...
        nnz(valid),nnz(tbl.isPlusOne==1),nnz(tbl.isMinusOne==1),nnz(~valid));
end

function scan_selftest()
    % 调度回归：稳定区不加密，状态改变区加密，避免重复端点。
    cfg=adaptive_scan_settings();
    right=empty_scan_row(); right.valid=true; right.deltaXOverA=1;
    right.isPlusOne=1; right.isMinusOne=0; right.isUnitQuantized=1;
    left=right; left.deltaXOverA=-1; left.isPlusOne=0; left.isMinusOne=1;
    assert(isempty(select_fine_ranges([0.1,1,2],{right,right,right},cfg)));
    ranges=select_fine_ranges([0.1,1,2],{right,right,left},cfg);
    assert(isequal(ranges,[1,2]));
    pts=make_fine_points([1,2;1,2],[0.1,1,2],0.1);
    assert(numel(pts)==9 && all(pts>1 & pts<2));
    assert(max(diff([1,pts,2]))<=0.1+1e-10);
    unknown=empty_scan_row();
    assert(isequal(select_fine_ranges([1,2],{right,unknown},cfg),[1,2]));
    % 回归检查：单条结构体的空文本/非空文本应均保留在同一行。
    testRow=empty_scan_row(); testRow.status='right_unit_quantized';
    oneTable=struct2table(testRow,'AsArray',true);
    assert(height(oneTable)==1 && isnan(oneTable.deltaXOverA));
    secondRow=testRow; secondRow.status='ok'; secondRow.note='nonempty note';
    twoTable=struct2table([testRow;secondRow],'AsArray',true);
    assert(height(twoTable)==2 && all(isnan(twoTable.deltaXOverA)));
    % 不求NLSE，使用已知中心轨迹验证检测、方向、周期差和无效值处理。
    x=(-12:0.025:12)'; fpp=100; t=(0:4*fpp)/fpp;
    lat.a=1; lat.nu=pi; % T=1
    o=struct('quantTolerance',0.05,'smoothSigma',0.35,'minPeakDistance',0.9, ...
        'coreRadius',0.9,'maxWidth',0.65,'minConcentration',0.65, ...
        'minMass',0.03,'minPeakRelative',1e-5,'isolationDistance',1.8, ...
        'contactDistance',1.8,'minApproach',0.3,'maxJump',0.35, ...
        'maxMassVariation',0.30,'maxWidthRatio',1.8,'edgeMargin',2.5, ...
        'edgeMassTolerance',1e-3,'normTolerance',1e-7, ...
        'relaxPeriods',0.25,'holdPeriods',0.15,'collisionEndT',0, ...
        'initialCenters',[-2,2]);
    dg=struct('normRelativeError',zeros(numel(t),1),'edgeMassFraction',zeros(numel(t),1));
    base=empty_scan_row();
    for v=[1,-1,0,0.6]
        rho=exp(-0.5*((x'-v*t')/0.22).^2);
        [r,~]=analyze_scan_case(t,rho,x,lat,fpp,o,base,dg);
        compactPackets=cell(numel(t),1);
        for it=1:numel(t), compactPackets{it}=detect_packets(rho(it,:),x,o); end
        rCompact=analyze_scan_case(t,compactPackets,x,lat,fpp,o,base,dg);
        assert(isequaln(r,rCompact),'在线特征与完整密度分析的结果不一致');
        assert(numel(r)==1 && r.valid,'单波包未被识别');
        assert(abs(r.deltaXOverA-v)<0.01,'有符号位移计算错误');
        assert(r.isPlusOne==double(v==1),'右向标记错误');
        assert(r.isMinusOne==double(v==-1),'左向标记错误');
        assert(r.isUnitQuantized==double(abs(v)==1),'单位位移标记错误');
        assert(abs(r.tEndOverT-r.tStartOverT-1)<1e-12,'窗口不是一个周期');
    end
    % 双峰内部结构不能被误识别成两只孤子。
    rho1=exp(-0.5*((x-0.22)/0.10).^2)+exp(-0.5*((x+0.22)/0.10).^2);
    p=detect_packets(rho1,x,o);
    assert(size(p,1)==1 && p(1,5),'内部双峰被误分裂');
    % 两个质量差别很大的波包仍应分别输出。
    rho=exp(-0.5*((x'-(-3+t'))/0.22).^2)+ ...
        0.1*exp(-0.5*((x'-(3-t'))/0.22).^2);
    [r,~]=analyze_scan_case(t,rho,x,lat,fpp,o,base,dg);
    assert(numel(r)==2 && all([r.valid]),'不等质量两波包跟踪失败');
    assert(abs(r(1).deltaXOverA-1)<0.01 && abs(r(2).deltaXOverA+1)<0.01);
    % 窗口中途消失不能报告零位移或量子化。
    rho=exp(-0.5*((x'-t')/0.22).^2); rho(t>0.7,:)=0;
    [r,~]=analyze_scan_case(t,rho,x,lat,fpp,o,base,dg);
    assert(~r.valid && isnan(r.deltaXOverA) && isnan(r.isPlusOne));
    % 边界污染与范数漂移均使判定失效。
    rho=exp(-0.5*((x'-t')/0.22).^2); bad=dg;
    bad.edgeMassFraction(:)=0.01;
    r=analyze_scan_case(t,rho,x,lat,fpp,o,base,bad);
    assert(~r.valid && strcmp(r.status,'boundary_contamination'));
    bad=dg; bad.normRelativeError(:)=1e-4;
    r=analyze_scan_case(t,rho,x,lat,fpp,o,base,bad);
    assert(~r.valid && strcmp(r.status,'norm_drift'));
    % 自动接触检测：先靠近、碰撞，随后左右各移一晶胞/周期。
    sep=abs(3-2*t');
    rho=exp(-0.5*((x'-sep/2)/0.22).^2)+exp(-0.5*((x'+sep/2)/0.22).^2);
    o.collisionEndT=[]; o.initialCenters=[-1.5,1.5];
    [r,a]=analyze_scan_case(t,rho,x,lat,fpp,o,base,dg);
    assert(isfinite(a.contactFrame) && all([r.valid]),'自动碰撞窗口失败');
    assert(all(abs(abs([r.deltaXOverA])-1)<0.01),'自动窗口位移错误');
    % 未碰撞的两孤子不能当作“碰撞后”结果。
    rho=repmat((exp(-0.5*((x-3)/0.22).^2)+exp(-0.5*((x+3)/0.22).^2))',numel(t),1);
    o.initialCenters=[-3,3];
    r=analyze_scan_case(t,rho,x,lat,fpp,o,base,dg);
    assert(~r.valid && strcmp(r.status,'no_collision_detected'));
    fprintf('All synthetic packet-tracking tests passed.\n');
end
