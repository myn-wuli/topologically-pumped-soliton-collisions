%% analyze_two_vector_collision_phase.m
% 两个矢量孤子（基态 + 激发态）碰撞的后处理：
%   1) 分别计算两个组分中基态/激发态波包的空间窗口粒子数；
%   2) 跟踪两个内部模的中心和密度交叠；
%   3) 根据出射轨迹区分俘获、穿越和反射；
%   4) 使用原程序求得的 psiGround 与 psi，扫描相对相位和间距，
%      计算相干密度干涉与双组分相互作用能。
% 中心跟踪默认先把两个矢量孤子分别在同一泵浦势中单独演化，得到无碰撞
% 动态参考态和参考轨迹，再在真实碰撞数据附近做锚定匹配。
%
% 推荐用法：
%   先运行 pengzhuang_two_vector_soliton_pump_g11_g22_g12(1).m，
%   保留工作区变量，再运行本程序。
% 或者取消原程序 save 区域的注释，并在 ana.dataFile 中填写 MAT 文件。
%
% 主粒子数曲线采用两个随孤子移动的软窗口；窗口重叠时用归一化分区，
% 保证同一部分密度不重复计数。程序还用 Lowdin 正交投影给出内部模权重，
% 但固定参考态在泵浦中会失配，所以投影权重只作为辅助诊断，不称为严格粒子数。

%% ============================= 用户设置 ==============================
ana.dataFile = '';            % 空：读取当前工作区；也可填写结果 MAT 文件
ana.caseName = 'current_phase';

% 中心跟踪：用基态/激发态参考密度做局域匹配。
% 原来的 searchHalfWidth=2 和 maxStep=1.2 相对于相邻保存帧过大，
% 会使激发态跟踪器跳到基态峰。下面的数值同时作为单帧硬位移限制。
ana.trackSearchHalfWidth = 0.45;
ana.trackWindowHalfWidth = 1.50;
ana.trackTaperWidth = 0.25;
ana.trackMaxStep = 0.25;
ana.trackMotionPenalty = 0.50;
ana.centerSmoothFrames = 3;

% 先分别演化孤立的矢量基态和矢量激发态，生成无碰撞动态参考轨迹。
% 这是多波瓣激发态的默认、也是最可靠的中心锚定方法。
ana.useIsolatedReferenceEvolution = true;
ana.referenceDt = [];                    % 空：优先读取原程序 pump.dt
ana.referenceUseFourthOrder = [];        % 空：优先读取 pump.useFourthOrder
ana.actualGroundExtraPhase = [];         % 空：读取工作区 vphi；仅影响参考态整体相位
ana.reuseReferenceCache = true;
ana.referenceCacheFile = 'isolated_vector_dynamic_references.mat';
ana.referenceAnchorHalfWidth = 0.65;
ana.trackMinMassFraction = 0.08;
ana.trackMinQuality = 0.05;
ana.survivalThreshold = 0.35;
ana.captureReceiverThreshold = 0.70;

% 粒子数窗口。半宽必须小于碰撞前两中心间距的一半，否则尚未碰撞时
% 两个窗口就会重叠并人为交换粒子数。
ana.modeWindowHalfWidth = 1.40;
ana.modeWindowTaperWidth = 0.25;
ana.removeTranslationalPhase = true;
ana.gramEigenvalueFloor = 1e-9;
ana.populationSmoothFrames = 5;

% 动力学俘获/穿越判据。
ana.ignoreInitialFraction = 0.05;
ana.postFraction = 0.20;
ana.captureDistance = 1.20;
ana.captureOverlap = 0.65;
ana.separatedDistance = 2.00;

% 静态相位--间距扫描。默认只改变第一组分中基态相对于激发态的相位；
% 第二组分的相对相位固定为 ana.fixedRelativePhase(2)。
ana.doPhaseScan = true;
ana.phaseScanComponent = 1;                 % 1 或 2
ana.fixedRelativePhase = [0; 0];
ana.phaseCases = [0, pi];
ana.phaseGrid = linspace(0,2*pi,181);
ana.separationGrid = linspace(0,6,81);
ana.scanCenter = [];                        % 空：取实际碰撞中心
ana.scanTime = [];                          % 空：取实际碰撞时刻
ana.scanSelectedDistance = [];              % 空：取实际最小间距
ana.scanWindowHalfWidth = 3.0;
ana.scanWindowTaperWidth = 0.4;

ana.saveOutput = true;
ana.outputPrefix = 'two_vector_collision_mode_population';

%% =============================== 读数据 ==============================
if isempty(ana.dataFile)
    needed = {'x','tdata','udata1','udata2','lat','psiGround'};
    for jj = 1:numel(needed)
        if exist(needed{jj},'var')~=1
            error(['当前工作区缺少变量 %s。请先运行碰撞演化程序，', ...
                   '或设置 ana.dataFile。'],needed{jj});
        end
    end
    S = struct();
    S.x = x;
    S.tdata = tdata;
    S.udata1 = udata1;
    S.udata2 = udata2;
    S.lat = lat;
    S.psiGround = psiGround;
    if exist('psi','var')==1, S.psi = psi; end
    if exist('psiExcited','var')==1, S.psiExcited = psiExcited; end
    if exist('job','var')==1, S.job = job; end
    if exist('pump','var')==1, S.pump = pump; end
    if exist('k','var')==1, S.k = k; end
    if exist('vphi','var')==1, S.vphi = vphi; end
else
    S = load(ana.dataFile);
end

if ~isfield(S,'psi')
    if isfield(S,'psiExcited')
        S.psi = S.psiExcited;
    else
        error('数据中缺少激发态参考解 psi（或 psiExcited）。');
    end
end

requiredFields = {'x','tdata','udata1','udata2','lat','psiGround','psi'};
for jj = 1:numel(requiredFields)
    if ~isfield(S,requiredFields{jj})
        error('输入数据中缺少变量 %s。',requiredFields{jj});
    end
end

x = S.x(:);
tdata = S.tdata(:).';
u1Data = S.udata1;
u2Data = S.udata2;
phiGround = S.psiGround(:);
phiExcited = S.psi(:);
lat = S.lat;
Nx = numel(x);
Nt = numel(tdata);
dx = x(2)-x(1);
Lx = Nx*dx;

if size(u1Data,1)~=Nx && size(u1Data,2)==Nx, u1Data=u1Data.'; end
if size(u2Data,1)~=Nx && size(u2Data,2)==Nx, u2Data=u2Data.'; end
if size(u1Data,1)~=Nx || size(u1Data,2)~=Nt || ...
   size(u2Data,1)~=Nx || size(u2Data,2)~=Nt
    error('udata1 和 udata2 的尺寸必须为 length(x) x length(tdata)。');
end
if numel(phiGround)~=Nx || numel(phiExcited)~=Nx
    error('psiGround 和 psi 的长度必须等于 length(x)。');
end

if isfield(S,'k') && numel(S.k)==Nx
    k = S.k(:);
elseif mod(Nx,2)==0
    k = [0:Nx/2-1,-Nx/2:-1].'*(2*pi/Lx);
else
    k = [0:(Nx-1)/2,-(Nx-1)/2:-1].'*(2*pi/Lx);
end

ana.centerSmoothFrames = odd_positive_collision(ana.centerSmoothFrames);
ana.populationSmoothFrames = odd_positive_collision(ana.populationSmoothFrames);

%% ======================== 跟踪基态和激发态中心 ======================
rhoTotal = abs(u1Data).^2+abs(u2Data).^2;
xGroundRef = density_center_collision(x,phiGround,dx);
xExcitedRef = density_center_collision(x,phiExcited,dx);

wGroundInitial = soft_window_collision(x,xGroundRef, ...
    ana.modeWindowHalfWidth,ana.modeWindowTaperWidth);
wExcitedInitial = soft_window_collision(x,xExcitedRef, ...
    ana.modeWindowHalfWidth,ana.modeWindowTaperWidth);
groundMassInWindow = real(dx*sum(wGroundInitial.*abs(phiGround).^2))/ ...
    (real(dx*sum(abs(phiGround).^2))+eps);
excitedMassInWindow = real(dx*sum(wExcitedInitial.*abs(phiExcited).^2))/ ...
    (real(dx*sum(abs(phiExcited).^2))+eps);
fprintf('reference mass inside particle window: ground=%.6f, excited=%.6f\n', ...
        groundMassInWindow,excitedMassInWindow);
if min(groundMassInWindow,excitedMassInWindow)<0.98
    warning(['当前粒子数窗口包含的参考态质量低于 98%%。可适当增大 ', ...
             'ana.modeWindowHalfWidth，但必须保持碰撞前两个窗口不重叠。']);
end

fprintf('\n========== Tracking ground/excited packets ==========\n');
if ana.useIsolatedReferenceEvolution
    referenceCacheSignature=reference_signature_collision( ...
        x,tdata,phiGround,phiExcited,lat,S,ana);
    loadedReferenceCache=false;
    if ana.reuseReferenceCache && isfile(ana.referenceCacheFile)
        cacheData=load(ana.referenceCacheFile,'isolatedReference', ...
                       'referenceCacheSignature');
        if isfield(cacheData,'isolatedReference') && ...
           isfield(cacheData,'referenceCacheSignature') && ...
           isequal(size(cacheData.referenceCacheSignature), ...
                   size(referenceCacheSignature)) && ...
           max(abs(cacheData.referenceCacheSignature-referenceCacheSignature))<1e-10
            isolatedReference=cacheData.isolatedReference;
            loadedReferenceCache=true;
            fprintf('Loaded isolated reference cache: %s\n',ana.referenceCacheFile);
        end
    end
    if ~loadedReferenceCache
        isolatedReference = evolve_isolated_references_collision( ...
            x,k,tdata,phiGround,phiExcited,lat,S,ana);
        if ana.reuseReferenceCache
            save(ana.referenceCacheFile,'isolatedReference', ...
                 'referenceCacheSignature','-v7.3');
            fprintf('Saved isolated reference cache: %s\n',ana.referenceCacheFile);
        end
    end
    rhoGroundReferenceData = abs(isolatedReference.g1).^2 ...
                           + abs(isolatedReference.g2).^2;
    rhoExcitedReferenceData = abs(isolatedReference.e1).^2 ...
                            + abs(isolatedReference.e2).^2;
    XgroundIsolated = isolatedReference.Xground;
    XexcitedIsolated = isolatedReference.Xexcited;

    [XgroundRaw,scoreGround,trackMassGround,validGround] = ...
        track_anchored_reference_collision(x,k,rhoTotal, ...
        rhoGroundReferenceData,XgroundIsolated,isolatedReference.normGround,ana);
    [XexcitedRaw,scoreExcited,trackMassExcited,validExcited] = ...
        track_anchored_reference_collision(x,k,rhoTotal, ...
        rhoExcitedReferenceData,XexcitedIsolated,isolatedReference.normExcited,ana);
    referenceApproachOverlap = reference_approach_overlap_collision( ...
        rhoGroundReferenceData,rhoExcitedReferenceData,dx);
else
    isolatedReference = struct();
    XgroundIsolated = nan(1,Nt);
    XexcitedIsolated = nan(1,Nt);
    [XgroundRaw,scoreGround] = track_template_center_collision( ...
        x,k,rhoTotal,phiGround,xGroundRef,ana);
    [XexcitedRaw,scoreExcited] = track_template_center_collision( ...
        x,k,rhoTotal,phiExcited,xExcitedRef,ana);
    trackMassGround = zeros(1,Nt);
    trackMassExcited = zeros(1,Nt);
    validGround = true(1,Nt);
    validExcited = true(1,Nt);
    referenceApproachOverlap = nan(1,Nt);
end

Xground = movmean(XgroundRaw,ana.centerSmoothFrames,'Endpoints','shrink');
Xexcited = movmean(XexcitedRaw,ana.centerSmoothFrames,'Endpoints','shrink');
XgroundDisplay = Xground;
XexcitedDisplay = Xexcited;
XgroundDisplay(~validGround) = nan;
XexcitedDisplay(~validExcited) = nan;
signedSeparation = Xexcited-Xground;
distanceGE = abs(signedSeparation);

%% =================== 四个模态粒子数与空间粒子数 ====================
N1ground = zeros(1,Nt); N1excited = zeros(1,Nt);
N2ground = zeros(1,Nt); N2excited = zeros(1,Nt);
N1groundSpatial = zeros(1,Nt); N1excitedSpatial = zeros(1,Nt);
N2groundSpatial = zeros(1,Nt); N2excitedSpatial = zeros(1,Nt);
N1local = zeros(1,Nt); N2local = zeros(1,Nt);
N1total = real(dx*sum(abs(u1Data).^2,1));
N2total = real(dx*sum(abs(u2Data).^2,1));
N1projectionResidual = zeros(1,Nt);
N2projectionResidual = zeros(1,Nt);
gramCondition1 = zeros(1,Nt);
gramCondition2 = zeros(1,Nt);
densityOverlap = zeros(1,Nt);

for it = 1:Nt
    xg = Xground(it);
    xe = Xexcited(it);
    wg = soft_window_collision(x,xg,ana.modeWindowHalfWidth, ...
                               ana.modeWindowTaperWidth);
    we = soft_window_collision(x,xe,ana.modeWindowHalfWidth, ...
                               ana.modeWindowTaperWidth);
    wUnion = min(1,wg+we);

    % 空间窗口粒子数采用归一化分区，保证两个窗口重叠时不重复计数。
    cover = min(1,wg+we);
    shareDen = wg+we+eps;
    wGroundPartition = cover.*wg./shareDen;
    wExcitedPartition = cover.*we./shareDen;
    rho1 = abs(u1Data(:,it)).^2;
    rho2 = abs(u2Data(:,it)).^2;
    N1groundSpatial(it) = real(dx*sum(wGroundPartition.*rho1));
    N1excitedSpatial(it) = real(dx*sum(wExcitedPartition.*rho1));
    N2groundSpatial(it) = real(dx*sum(wGroundPartition.*rho2));
    N2excitedSpatial(it) = real(dx*sum(wExcitedPartition.*rho2));
    N1local(it) = real(dx*sum(wUnion.*rho1));
    N2local(it) = real(dx*sum(wUnion.*rho2));

    if ana.useIsolatedReferenceEvolution
        phiG1Now = isolatedReference.g1(:,it);
        phiG2Now = isolatedReference.g2(:,it);
        phiE1Now = isolatedReference.e1(:,it);
        phiE2Now = isolatedReference.e2(:,it);
        xgReferenceNow = XgroundIsolated(it);
        xeReferenceNow = XexcitedIsolated(it);
        rhoGReferenceNow = abs(phiG1Now).^2+abs(phiG2Now).^2;
        rhoEReferenceNow = abs(phiE1Now).^2+abs(phiE2Now).^2;
    else
        phiG1Now = phiGround; phiG2Now = phiGround;
        phiE1Now = phiExcited; phiE2Now = phiExcited;
        xgReferenceNow = xGroundRef;
        xeReferenceNow = xExcitedRef;
        rhoGReferenceNow = abs(phiGround).^2;
        rhoEReferenceNow = abs(phiExcited).^2;
    end

    % 将孤立动态参考态仅作小幅平移到真实数据的局域中心，再计算几何交叠。
    rhoGOverlap = real(spectral_shift_collision( ...
        rhoGReferenceNow,k,xg-xgReferenceNow));
    rhoEOverlap = real(spectral_shift_collision( ...
        rhoEReferenceNow,k,xe-xeReferenceNow));
    rhoGOverlap=max(0,rhoGOverlap);
    rhoEOverlap=max(0,rhoEOverlap);
    rhoECoCentered = real(spectral_shift_collision( ...
        rhoEReferenceNow,k,xgReferenceNow-xeReferenceNow));
    rhoECoCentered=max(0,rhoECoCentered);
    overlapRaw = density_bhattacharyya_collision(rhoGOverlap,rhoEOverlap,dx);
    overlapMaximumNow = density_bhattacharyya_collision( ...
        rhoGReferenceNow,rhoECoCentered,dx);
    densityOverlap(it)=min(1,max(0,overlapRaw/(overlapMaximumNow+eps)));
    if ana.useIsolatedReferenceEvolution && ...
       ~(validGround(it) && validExcited(it))
        densityOverlap(it)=nan;
    end

    [N1ground(it),N1excited(it),gramCondition1(it)] = ...
        mode_population_collision(u1Data(:,it),x,k,dx,wUnion,wg,we, ...
        phiG1Now,phiE1Now,xgReferenceNow,xeReferenceNow,xg,xe,ana);
    [N2ground(it),N2excited(it),gramCondition2(it)] = ...
        mode_population_collision(u2Data(:,it),x,k,dx,wUnion,wg,we, ...
        phiG2Now,phiE2Now,xgReferenceNow,xeReferenceNow,xg,xe,ana);

    % 保留符号：若该量明显为负，说明窗口、中心或 Gram 矩阵条件数需检查。
    N1projectionResidual(it) = N1local(it)-N1ground(it)-N1excited(it);
    N2projectionResidual(it) = N2local(it)-N2ground(it)-N2excited(it);
end

N1groundPlot = movmean(N1ground,ana.populationSmoothFrames,'Endpoints','shrink');
N1excitedPlot = movmean(N1excited,ana.populationSmoothFrames,'Endpoints','shrink');
N2groundPlot = movmean(N2ground,ana.populationSmoothFrames,'Endpoints','shrink');
N2excitedPlot = movmean(N2excited,ana.populationSmoothFrames,'Endpoints','shrink');
N1groundSpatialPlot = movmean(N1groundSpatial,ana.populationSmoothFrames, ...
                              'Endpoints','shrink');
N1excitedSpatialPlot = movmean(N1excitedSpatial,ana.populationSmoothFrames, ...
                               'Endpoints','shrink');
N2groundSpatialPlot = movmean(N2groundSpatial,ana.populationSmoothFrames, ...
                              'Endpoints','shrink');
N2excitedSpatialPlot = movmean(N2excitedSpatial,ana.populationSmoothFrames, ...
                               'Endpoints','shrink');

%% ========================== 碰撞结果分类 ============================
iStart = min(Nt,max(1,round(ana.ignoreInitialFraction*Nt)));
if ana.useIsolatedReferenceEvolution && any(isfinite(referenceApproachOverlap))
    collisionIndicator=referenceApproachOverlap;
else
    collisionIndicator=densityOverlap;
end
[overlapMaximum,iOverlapRelative] = max(collisionIndicator(iStart:end));
if ana.useIsolatedReferenceEvolution
    % 孤立参考轨迹不会因俘获而形成长期平台，因此最大交叠点就是无歧义的
    % 预期碰撞中心；不再取 85% 阈值的碰撞前沿。
    iCollision=iStart+iOverlapRelative-1;
else
    overlapThreshold=max(0.40,0.85*overlapMaximum);
    iFirstOverlapRelative=find( ...
        collisionIndicator(iStart:end)>=overlapThreshold,1,'first');
    if isempty(iFirstOverlapRelative)
        [~,iClosestRelative] = min(distanceGE(iStart:end));
        iCollision=iStart+iClosestRelative-1;
    else
        % 对俘获形成的长期交叠平台，取第一次进入强交叠区。
        iCollision=iStart+iFirstOverlapRelative-1;
    end
end

if ana.useIsolatedReferenceEvolution
    [collisionClass,collisionInfo] = classify_collision_survival_collision( ...
        tdata,Xground,Xexcited,trackMassGround,trackMassExcited, ...
        validGround,validExcited,iCollision,ana);
else
    [collisionClass,collisionInfo] = classify_collision_collision( ...
        tdata,Xground,Xexcited,densityOverlap,iCollision,ana);
end

fprintf('\n========== Dynamic collision result ==========\n');
fprintf('collision frame/time     = %d / %.12g\n',iCollision,tdata(iCollision));
fprintf('minimum separation       = %.8g\n',min(distanceGE));
fprintf('post-collision separation= %.8g\n',collisionInfo.postDistance);
if isfield(collisionInfo,'postOverlap')
    fprintf('post-collision overlap   = %.8g\n',collisionInfo.postOverlap);
end
if isfield(collisionInfo,'groundSurvival')
    fprintf('ground/excited survival  = %.8g / %.8g\n', ...
        collisionInfo.groundSurvival,collisionInfo.excitedSurvival);
end
fprintf('trajectory sign changed  = %d\n',collisionInfo.signChanged);
fprintf('classification           = %s\n',collisionClass);

%% ============================== 作图 ================================
if isfield(lat,'nu') && abs(lat.nu)>0
    pumpPeriod = pi/abs(lat.nu);
    tPlot = tdata/pumpPeriod;
    timeLabel = '$t/T$';
else
    pumpPeriod = nan;
    tPlot = tdata;
    timeLabel = '$t$';
end

normDrift1=max(abs(N1total-N1total(1)))/(abs(N1total(1))+eps);
normDrift2=max(abs(N2total-N2total(1)))/(abs(N2total(1))+eps);
if isfinite(pumpPeriod)
    earlyIndex=find(tPlot<=0.5);
else
    earlyIndex=find(tdata<=0.5);
end
if numel(earlyIndex)<2, earlyIndex=1:min(Nt,max(2,round(0.1*Nt))); end
earlyVariation=[relative_span_collision(N1groundSpatial(earlyIndex)), ...
                relative_span_collision(N1excitedSpatial(earlyIndex)); ...
                relative_span_collision(N2groundSpatial(earlyIndex)), ...
                relative_span_collision(N2excitedSpatial(earlyIndex))];
fprintf('\n========== Conservation/window diagnostics ==========\n');
fprintf('total norm relative drift: N1 = %.3e, N2 = %.3e\n', ...
        normDrift1,normDrift2);
fprintf('relative window variation before t/T=0.5:\n');
fprintf('  component 1: ground = %.3e, excited = %.3e\n', ...
        earlyVariation(1,1),earlyVariation(1,2));
fprintf('  component 2: ground = %.3e, excited = %.3e\n', ...
        earlyVariation(2,1),earlyVariation(2,2));
minimumEarlyDistance=min(distanceGE(earlyIndex));
fprintf('minimum center distance before t/T=0.5 = %.6g; 2L_w = %.6g\n', ...
        minimumEarlyDistance,2*ana.modeWindowHalfWidth);
if minimumEarlyDistance<=2*ana.modeWindowHalfWidth
    warning(['t/T<0.5 时两个粒子数窗口已经发生几何重叠。请继续减小 ', ...
             'ana.modeWindowHalfWidth，或把无碰撞检验区间缩短。']);
end

figure('Color','w','Name','Windowed particle numbers in two components');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
plot(tPlot,N1groundSpatialPlot,'r-','LineWidth',1.7, ...
    'DisplayName','$N_{1,g}^{\rm win}$');
plot(tPlot,N1excitedSpatialPlot,'b--','LineWidth',1.7, ...
    'DisplayName','$N_{1,e}^{\rm win}$');
plot(tPlot,N1total,'k:','LineWidth',1.2, ...
    'DisplayName','$N_1$');
xline(tPlot(iCollision),'Color',[0.35 0.35 0.35],'LineStyle',':', ...
    'HandleVisibility','off');
ylabel('particle number');
title('Component 1');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);

nexttile; hold on; box on;
plot(tPlot,N2groundSpatialPlot,'r-','LineWidth',1.7, ...
    'DisplayName','$N_{2,g}^{\rm win}$');
plot(tPlot,N2excitedSpatialPlot,'b--','LineWidth',1.7, ...
    'DisplayName','$N_{2,e}^{\rm win}$');
plot(tPlot,N2total,'k:','LineWidth',1.2, ...
    'DisplayName','$N_2$');
xline(tPlot(iCollision),'Color',[0.35 0.35 0.35],'LineStyle',':', ...
    'HandleVisibility','off');
xlabel(timeLabel,'Interpreter','latex');
ylabel('particle number');
title('Component 2');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);

% 固定参考态的投影只作为内部模诊断量。泵浦势随时间变化时它可以振荡，
% 因此这里记作 P^{proj}，不再把它误称为严格粒子数。
figure('Color','w','Name','Projected ground/excited modal weights');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
plot(tPlot,N1groundPlot,'r-','LineWidth',1.5, ...
    'DisplayName','$P_{1,g}^{\rm proj}$');
plot(tPlot,N1excitedPlot,'b--','LineWidth',1.5, ...
    'DisplayName','$P_{1,e}^{\rm proj}$');
xline(tPlot(iCollision),'k:','HandleVisibility','off');
ylabel('projected weight'); title('Component 1');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);
nexttile; hold on; box on;
plot(tPlot,N2groundPlot,'r-','LineWidth',1.5, ...
    'DisplayName','$P_{2,g}^{\rm proj}$');
plot(tPlot,N2excitedPlot,'b--','LineWidth',1.5, ...
    'DisplayName','$P_{2,e}^{\rm proj}$');
xline(tPlot(iCollision),'k:','HandleVisibility','off');
xlabel(timeLabel,'Interpreter','latex');
ylabel('projected weight'); title('Component 2');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);

figure('Color','w','Name','Collision centers and density overlap');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
if ana.useIsolatedReferenceEvolution
    plot(tPlot,XgroundIsolated,':','Color',[1 0.55 0.55],'LineWidth',1.1, ...
        'DisplayName','$X_g^{\rm iso}$');
    plot(tPlot,XexcitedIsolated,':','Color',[0.45 0.65 1],'LineWidth',1.1, ...
        'DisplayName','$X_e^{\rm iso}$');
end
plot(tPlot,XgroundDisplay,'r-','LineWidth',1.6,'DisplayName','$X_g$');
plot(tPlot,XexcitedDisplay,'b--','LineWidth',1.6,'DisplayName','$X_e$');
xline(tPlot(iCollision),'k:','HandleVisibility','off');
ylabel('center');
legend('Location','best','Interpreter','latex');
title(sprintf('%s: %s',ana.caseName,collisionClass),'Interpreter','none');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);

nexttile; hold on; box on;
trackBaselineIndex=1:max(3,round(0.1*Nt));
trackBaselineGround=median(trackMassGround(trackBaselineIndex));
trackBaselineExcited=median(trackMassExcited(trackBaselineIndex));
plot(tPlot,trackMassGround/(trackBaselineGround+eps),'r-', ...
    'LineWidth',1.4,'DisplayName','$S_g$');
plot(tPlot,trackMassExcited/(trackBaselineExcited+eps),'b--', ...
    'LineWidth',1.4,'DisplayName','$S_e$');
yline(ana.survivalThreshold,'k:','LineWidth',1, ...
    'DisplayName','survival threshold');
xline(tPlot(iCollision),'k:','HandleVisibility','off');
ylabel('survival ratio');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1);

nexttile; hold on; box on;
plot(tPlot,densityOverlap,'m-','LineWidth',1.6, ...
    'DisplayName','$\mathcal O_\rho^{\rm data}$');
if ana.useIsolatedReferenceEvolution
    plot(tPlot,referenceApproachOverlap,'Color',[0.2 0.55 0.2], ...
        'LineStyle','--','LineWidth',1.2, ...
        'DisplayName','$\mathcal O_\rho^{\rm iso}$');
else
    plot(tPlot,distanceGE/max(max(distanceGE),eps),'Color',[0.2 0.55 0.2], ...
        'LineStyle','--','LineWidth',1.2,'DisplayName','$D/D_{\max}$');
end
xline(tPlot(iCollision),'k:','HandleVisibility','off');
xlabel(timeLabel,'Interpreter','latex');
ylabel('overlap / distance');
ylim([0,1.05]);
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1);

%% ======================== 相位--间距静态扫描 ========================
phaseScan = struct();
if ana.doPhaseScan
    fprintf('\n========== Static phase/separation scan ==========\n');
    if isempty(ana.scanCenter)
        if ana.useIsolatedReferenceEvolution
            scanCenter = 0.5*(XgroundIsolated(iCollision)+ ...
                              XexcitedIsolated(iCollision));
        else
            scanCenter = 0.5*(Xground(iCollision)+Xexcited(iCollision));
        end
    else
        scanCenter = ana.scanCenter;
    end
    if isempty(ana.scanTime)
        scanTime = tdata(iCollision);
    else
        scanTime = ana.scanTime;
    end
    if isempty(ana.scanSelectedDistance)
        if ana.useIsolatedReferenceEvolution
            selectedDistance = abs(XexcitedIsolated(iCollision)- ...
                                   XgroundIsolated(iCollision));
        else
            selectedDistance = min(distanceGE(iStart:end));
        end
    else
        selectedDistance = ana.scanSelectedDistance;
    end

    phaseScan = scan_phase_separation_collision( ...
        x,k,dx,phiGround,phiExcited,xGroundRef,xExcitedRef, ...
        lat,S,scanCenter,scanTime,selectedDistance,ana);

    fprintf('scan center/time         = %.8g / %.8g\n',scanCenter,scanTime);
    fprintf('selected separation      = %.8g\n',phaseScan.selectedDistance);
    fprintf('U_int(pi)-U_int(0)       = %.12g\n',phaseScan.Upi-phaseScan.U0);
    fprintf('J_eff=[U(pi)-U(0)]/2     = %.12g\n',phaseScan.Jeffective);
    fprintf('static phase preference  = %s\n',phaseScan.preference);

    figure('Color','w','Name','Phase-dependent density overlap and interaction energy');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; box on;
    for jp = 1:numel(ana.phaseCases)
        plot(phaseScan.separation,phaseScan.interactionEnergy(jp,:), ...
            'LineWidth',1.6,'DisplayName',phase_label_collision(ana.phaseCases(jp)));
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('$D$','Interpreter','latex'); ylabel('$U_{\rm int}$','Interpreter','latex');
    legend('Location','best','Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1);

    nexttile; hold on; box on;
    for jp = 1:numel(ana.phaseCases)
        plot(phaseScan.separation,phaseScan.densityInterference(jp,:), ...
            'LineWidth',1.6,'DisplayName',phase_label_collision(ana.phaseCases(jp)));
    end
    yline(0,'k:','HandleVisibility','off');
    xlabel('$D$','Interpreter','latex');
    ylabel('$\mathcal I_\rho$','Interpreter','latex');
    legend('Location','best','Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1);

    nexttile; hold on; box on;
    plot(phaseScan.phaseGrid/pi,phaseScan.phaseInteractionEnergy, ...
        'k-','LineWidth',1.6);
    xlabel('$\Delta\phi/\pi$','Interpreter','latex');
    ylabel('$U_{\rm int}$','Interpreter','latex');
    title(sprintf('$D=%.3g$',phaseScan.selectedDistance),'Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1);

    nexttile; hold on; box on;
    plot(phaseScan.phaseGrid/pi,phaseScan.phaseDensityInterference, ...
        'Color',[0.75 0.1 0.55],'LineWidth',1.6);
    yline(0,'k:','HandleVisibility','off');
    xlabel('$\Delta\phi/\pi$','Interpreter','latex');
    ylabel('$\mathcal I_\rho$','Interpreter','latex');
    title('coherent density interference');
    set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1);
end

%% ================================ 输出 ===============================
CollisionTable = table(tdata(:),XgroundRaw(:),XexcitedRaw(:), ...
    Xground(:),Xexcited(:),XgroundIsolated(:),XexcitedIsolated(:), ...
    validGround(:),validExcited(:),trackMassGround(:),trackMassExcited(:), ...
    distanceGE(:),signedSeparation(:),densityOverlap(:), ...
    referenceApproachOverlap(:),scoreGround(:),scoreExcited(:), ...
    N1total(:),N1local(:),N1ground(:),N1excited(:), ...
    N2total(:),N2local(:),N2ground(:),N2excited(:), ...
    N1groundSpatial(:),N1excitedSpatial(:), ...
    N2groundSpatial(:),N2excitedSpatial(:), ...
    N1projectionResidual(:),N2projectionResidual(:), ...
    gramCondition1(:),gramCondition2(:), ...
    'VariableNames',{'Time','XgroundRaw','XexcitedRaw','Xground','Xexcited', ...
    'XgroundIsolated','XexcitedIsolated','GroundTrackValid','ExcitedTrackValid', ...
    'GroundTrackMass','ExcitedTrackMass','Distance','SignedSeparation', ...
    'DensityOverlap','ReferenceApproachOverlap','GroundTrackScore', ...
    'ExcitedTrackScore','N1Total','N1Local','N1GroundMode','N1ExcitedMode', ...
    'N2Total','N2Local','N2GroundMode','N2ExcitedMode', ...
    'N1GroundSpatial','N1ExcitedSpatial','N2GroundSpatial','N2ExcitedSpatial', ...
    'N1ProjectionResidual','N2ProjectionResidual','GramCondition1','GramCondition2'});

analysisInfo = struct();
analysisInfo.caseName = ana.caseName;
analysisInfo.collisionClass = collisionClass;
analysisInfo.collisionFrame = iCollision;
analysisInfo.collisionTime = tdata(iCollision);
analysisInfo.collisionTimeOverT = tdata(iCollision)/pumpPeriod;
analysisInfo.minimumDistance = min(distanceGE);
analysisInfo.collision = collisionInfo;
analysisInfo.groundReferenceCenter = xGroundRef;
analysisInfo.excitedReferenceCenter = xExcitedRef;
analysisInfo.usedIsolatedReferenceEvolution = ana.useIsolatedReferenceEvolution;
analysisInfo.referenceApproachOverlap = referenceApproachOverlap;
analysisInfo.groundMassInWindow = groundMassInWindow;
analysisInfo.excitedMassInWindow = excitedMassInWindow;
analysisInfo.totalNormRelativeDrift = [normDrift1,normDrift2];
analysisInfo.earlyWindowRelativeVariation = earlyVariation;
analysisInfo.phaseScan = phaseScan;

if ana.saveOutput
    outputStem = [ana.outputPrefix '_' safe_name_collision(ana.caseName)];
    writetable(CollisionTable,[outputStem '.csv']);
    save([outputStem '.mat'],'CollisionTable','analysisInfo','phaseScan','ana', ...
        'phiGround','phiExcited','N1ground','N1excited','N2ground','N2excited', ...
        'Xground','Xexcited','densityOverlap','-v7.3');
    fprintf('\nSaved: %s.csv and %s.mat\n',outputStem,outputStem);
end

%% ====================================================================
%% 局部函数
%% ====================================================================
function n = odd_positive_collision(n)
    n = max(1,round(n));
    if mod(n,2)==0, n=n+1; end
end

function value = relative_span_collision(y)
    y=real(y(:));
    value=(max(y)-min(y))/(abs(mean(y))+eps);
end

function V = external_potential_collision(x,t,lat)
    V = -lat.p1*cos(pi*x/lat.d1).^2 ...
        -lat.p2*cos(pi*x/lat.d2-lat.nu*t).^2;
end

function xc = density_center_collision(x,u,dx)
    rho = abs(u).^2;
    xc = real(dx*sum(x.*rho)/(dx*sum(rho)+eps));
end

function uShift = spectral_shift_collision(u,k,delta)
    uShift = ifft(fft(u(:)).*exp(-1i*k*delta));
end

function w = soft_window_collision(x,c,halfWidth,taperWidth)
    taperWidth = min(max(taperWidth,0),halfWidth);
    r = abs(x-c);
    w = zeros(size(x));
    inner = max(0,halfWidth-taperWidth);
    w(r<=inner)=1;
    if taperWidth>0
        idx = r>inner & r<halfWidth;
        w(idx)=0.5*(1+cos(pi*(r(idx)-inner)/taperWidth));
    end
end

function value = density_bhattacharyya_collision(rhoA,rhoB,dx)
    rhoA=max(0,real(rhoA));
    rhoB=max(0,real(rhoB));
    numerator=(dx*sum(sqrt(rhoA.*rhoB)))^2;
    denominator=(dx*sum(rhoA))*(dx*sum(rhoB));
    value=real(numerator/(denominator+eps));
end

function overlap = reference_approach_overlap_collision(rhoG,rhoE,dx)
    Nt=size(rhoG,2);
    overlap=zeros(1,Nt);
    for it=1:Nt
        overlap(it)=density_bhattacharyya_collision( ...
            rhoG(:,it),rhoE(:,it),dx);
    end
    peak=max(overlap);
    if peak>eps, overlap=overlap/peak; end
    overlap=min(1,max(0,overlap));
end

function xc = density_center_from_rho_collision(x,rho,dx)
    xc=real(dx*sum(x.*rho)/(dx*sum(rho)+eps));
end

function signature = reference_signature_collision( ...
        x,tdata,phiGround,phiExcited,lat,S,ana)
    if isempty(ana.referenceDt)
        if isfield(S,'pump') && isfield(S.pump,'dt')
            dt=S.pump.dt;
        else
            dt=1e-3;
        end
    else
        dt=ana.referenceDt;
    end
    if isempty(ana.referenceUseFourthOrder)
        if isfield(S,'pump') && isfield(S.pump,'useFourthOrder')
            fourth=double(S.pump.useFourthOrder);
        else
            fourth=1;
        end
    else
        fourth=double(ana.referenceUseFourthOrder);
    end
    weightG=1; weightE=1;
    spinG=[1;1]/sqrt(2); spinE=spinG;
    if isfield(S,'job')
        if isfield(S.job,'groundWeight'), weightG=S.job.groundWeight; end
        if isfield(S.job,'excitWeight'), weightE=S.job.excitWeight; end
        if isfield(S.job,'groundSpinor'), spinG=S.job.groundSpinor(:); end
        if isfield(S.job,'excitedSpinor'), spinE=S.job.excitedSpinor(:); end
    end
    dx=x(2)-x(1);
    signature=[numel(x),numel(tdata),dx,tdata(end),dt,fourth, ...
        lat.p1,lat.p2,lat.d1,lat.d2,lat.nu,lat.g11,lat.g22,lat.g12, ...
        real(dx*sum(abs(phiGround).^2)),real(dx*sum(abs(phiGround).^4)), ...
        real(dx*sum(abs(phiExcited).^2)),real(dx*sum(abs(phiExcited).^4)), ...
        weightG,weightE,abs(spinG(1)),abs(spinG(2)), ...
        abs(spinE(1)),abs(spinE(2))];
end

function ref = evolve_isolated_references_collision( ...
        x,k,tdata,phiGround,phiExcited,lat,S,ana)
    fprintf('\n========== Isolated dynamic-reference evolution ==========\n');
    if isfield(S,'job'), job=S.job; else, job=struct(); end
    if isfield(job,'groundSpinor')
        chiG=job.groundSpinor(:);
    else
        chiG=[1;1]/sqrt(2);
    end
    if isfield(job,'excitedSpinor')
        chiE=job.excitedSpinor(:);
    else
        chiE=[1;1]/sqrt(2);
    end
    if isfield(job,'groundPhase')
        chiG=chiG.*exp(1i*job.groundPhase(:));
    end
    if isfield(job,'excitedPhase')
        chiE=chiE.*exp(1i*job.excitedPhase(:));
    end
    chiG=chiG/sqrt(sum(abs(chiG).^2)+eps);
    chiE=chiE/sqrt(sum(abs(chiE).^2)+eps);
    if isfield(job,'groundWeight'), weightG=job.groundWeight; else, weightG=1; end
    if isfield(job,'excitWeight'), weightE=job.excitWeight; else, weightE=1; end

    if isempty(ana.actualGroundExtraPhase)
        if isfield(S,'vphi'), extraPhase=S.vphi; else, extraPhase=0; end
    else
        extraPhase=ana.actualGroundExtraPhase;
    end
    uG1=weightG*chiG(1)*phiGround*exp(1i*extraPhase);
    uG2=weightG*chiG(2)*phiGround;
    uE1=weightE*chiE(1)*phiExcited;
    uE2=weightE*chiE(2)*phiExcited;

    if isempty(ana.referenceDt)
        if isfield(S,'pump') && isfield(S.pump,'dt')
            dt=S.pump.dt;
        else
            dt=1e-3;
            warning('数据中没有 pump.dt，孤立参考演化采用 ana.referenceDt=1e-3。');
        end
    else
        dt=ana.referenceDt;
    end
    if isempty(ana.referenceUseFourthOrder)
        if isfield(S,'pump') && isfield(S.pump,'useFourthOrder')
            useFourthOrder=S.pump.useFourthOrder;
        else
            useFourthOrder=true;
        end
    else
        useFourthOrder=ana.referenceUseFourthOrder;
    end

    targetStep=round(tdata/dt);
    if max(abs(targetStep*dt-tdata))>max(1e-10,1e-6*dt)
        warning('tdata 不是参考时间步长的整数倍；参考态保存时刻采用最近时间步。');
    end
    if targetStep(1)~=0
        error('孤立参考演化要求 tdata(1)=0。');
    end
    Nt=numel(tdata);
    Nx=numel(x);
    g1=complex(zeros(Nx,Nt)); g2=complex(zeros(Nx,Nt));
    e1=complex(zeros(Nx,Nt)); e2=complex(zeros(Nx,Nt));
    g1(:,1)=uG1; g2(:,1)=uG2; e1(:,1)=uE1; e2(:,1)=uE2;

    if useFourthOrder
        c=1/(2-2^(1/3));
        a=[c/2,(1-c)/2,(1-c)/2,c/2];
        b=[c,1-2*c,c];
    else
        a=[1/2,1/2];
        b=1;
    end
    kinetic=cell(numel(a),1);
    for jj=1:numel(a)
        kinetic{jj}=exp(-1i*a(jj)*dt*(k.^2)/2);
    end

    nextFrame=2;
    nmax=max(targetStep);
    reportStep=max(1,round(nmax/10));
    for nn=1:nmax
        tMid=(nn-0.5)*dt;
        V=external_potential_collision(x,tMid,lat);
        if useFourthOrder
            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{1});
            [uG1,uG2]=nonlinear_pair_collision(uG1,uG2,V,b(1)*dt,lat);
            [uE1,uE2]=nonlinear_pair_collision(uE1,uE2,V,b(1)*dt,lat);

            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{2});
            [uG1,uG2]=nonlinear_pair_collision(uG1,uG2,V,b(2)*dt,lat);
            [uE1,uE2]=nonlinear_pair_collision(uE1,uE2,V,b(2)*dt,lat);

            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{3});
            [uG1,uG2]=nonlinear_pair_collision(uG1,uG2,V,b(3)*dt,lat);
            [uE1,uE2]=nonlinear_pair_collision(uE1,uE2,V,b(3)*dt,lat);

            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{4});
        else
            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{1});
            [uG1,uG2]=nonlinear_pair_collision(uG1,uG2,V,dt,lat);
            [uE1,uE2]=nonlinear_pair_collision(uE1,uE2,V,dt,lat);
            [uG1,uG2,uE1,uE2]=kinetic_four_collision( ...
                uG1,uG2,uE1,uE2,kinetic{2});
        end

        while nextFrame<=Nt && targetStep(nextFrame)==nn
            g1(:,nextFrame)=uG1; g2(:,nextFrame)=uG2;
            e1(:,nextFrame)=uE1; e2(:,nextFrame)=uE2;
            nextFrame=nextFrame+1;
        end
        if mod(nn,reportStep)==0 || nn==nmax
            fprintf('isolated references: %6.1f%%\n',100*nn/nmax);
        end
    end
    if nextFrame<=Nt
        error('孤立参考演化没有覆盖全部 tdata 保存帧。');
    end

    Xg=zeros(1,Nt); Xe=zeros(1,Nt);
    for it=1:Nt
        Xg(it)=density_center_from_rho_collision( ...
            x,abs(g1(:,it)).^2+abs(g2(:,it)).^2,x(2)-x(1));
        Xe(it)=density_center_from_rho_collision( ...
            x,abs(e1(:,it)).^2+abs(e2(:,it)).^2,x(2)-x(1));
    end
    dx=x(2)-x(1);
    ref=struct();
    ref.g1=g1; ref.g2=g2; ref.e1=e1; ref.e2=e2;
    ref.Xground=Xg; ref.Xexcited=Xe;
    ref.normGround=real(dx*sum(abs(g1(:,1)).^2+abs(g2(:,1)).^2));
    ref.normExcited=real(dx*sum(abs(e1(:,1)).^2+abs(e2(:,1)).^2));
    ref.dt=dt; ref.useFourthOrder=useFourthOrder;
end

function [uG1,uG2,uE1,uE2] = kinetic_four_collision( ...
        uG1,uG2,uE1,uE2,factor)
    uG1=ifft(fft(uG1).*factor);
    uG2=ifft(fft(uG2).*factor);
    uE1=ifft(fft(uE1).*factor);
    uE2=ifft(fft(uE2).*factor);
end

function [u1,u2] = nonlinear_pair_collision(u1,u2,V,dtStage,lat)
    rho1=abs(u1).^2;
    rho2=abs(u2).^2;
    u1=u1.*exp(-1i*dtStage*(V+lat.g11*rho1+lat.g12*rho2));
    u2=u2.*exp(-1i*dtStage*(V+lat.g22*rho2+lat.g12*rho1));
end

function [center,quality,localMass,valid] = track_anchored_reference_collision( ...
        x,k,rhoData,rhoReferenceData,referenceCenter,referenceNorm,ana)
    Nt=size(rhoData,2);
    center=zeros(1,Nt);
    quality=zeros(1,Nt);
    localMass=zeros(1,Nt);
    valid=false(1,Nt);
    dx=x(2)-x(1);

    for it=1:Nt
        anchor=referenceCenter(it);
        idxCandidate=find(abs(x-anchor)<=ana.referenceAnchorHalfWidth);
        if it>1 && valid(it-1)
            referenceStep=abs(referenceCenter(it)-referenceCenter(it-1));
            maximumStep=ana.trackMaxStep+referenceStep;
            idxCandidate=idxCandidate(abs(x(idxCandidate)-center(it-1))<=maximumStep);
        end
        if isempty(idxCandidate)
            [~,ic]=min(abs(x-anchor));
            idxCandidate=ic;
        end
        rhoNow=rhoData(:,it);
        rhoReference=rhoReferenceData(:,it);
        bestScore=-inf;
        bestCenter=anchor;
        bestMass=0;
        for jc=1:numel(idxCandidate)
            c=x(idxCandidate(jc));
            rhoRefShift=real(spectral_shift_collision( ...
                rhoReference,k,c-referenceCenter(it)));
            rhoRefShift=max(0,rhoRefShift);
            w=soft_window_collision(x,c,ana.trackWindowHalfWidth, ...
                                    ana.trackTaperWidth);
            massNow=real(dx*sum(w.*rhoNow));
            bc=real(dx*sum(w.*sqrt(max(0,rhoRefShift.*rhoNow))));
            nr=real(dx*sum(w.*rhoRefShift));
            shapeScore=bc^2/(nr*massNow+eps);
            anchorCost=ana.trackMotionPenalty*((c-anchor)/ ...
                       max(ana.referenceAnchorHalfWidth,eps))^2;
            thisScore=shapeScore-anchorCost;
            if thisScore>bestScore
                bestScore=thisScore;
                bestCenter=c;
                bestMass=massNow;
            end
        end
        quality(it)=bestScore;
        localMass(it)=bestMass;
        valid(it)=bestMass>=ana.trackMinMassFraction*referenceNorm && ...
                  bestScore>=ana.trackMinQuality;
        if valid(it)
            center(it)=bestCenter;
        else
            % 没有足够密度时回到无碰撞参考轨迹，绝不继续追踪噪声。
            center(it)=anchor;
        end
    end
end

function [center,quality] = track_template_center_collision( ...
        x,k,rhoData,phiRef,xRef,ana)
    Nt = size(rhoData,2);
    center = zeros(1,Nt);
    quality = zeros(1,Nt);
    center(1)=xRef;
    rhoTemplate0 = abs(phiRef).^2;

    for it=1:Nt
        if it==1
            prediction=xRef;
        elseif it==2
            prediction=center(it-1);
        else
            step=center(it-1)-center(it-2);
            step=max(-ana.trackMaxStep,min(ana.trackMaxStep,step));
            prediction=center(it-1)+step;
        end
        % 同时限制在预测点附近和上一帧中心附近。trackMaxStep 在这里是
        % 实际硬约束，防止某一帧突然从激发态峰跳到基态峰。
        if it==1
            allowedByPrevious=true(size(x));
        else
            allowedByPrevious=abs(x-center(it-1))<=ana.trackMaxStep;
        end
        idxCandidate=find(abs(x-prediction)<=ana.trackSearchHalfWidth & ...
                          allowedByPrevious);
        if isempty(idxCandidate)
            if it==1
                [~,ic]=min(abs(x-prediction));
            else
                target=max(center(it-1)-ana.trackMaxStep, ...
                           min(center(it-1)+ana.trackMaxStep,prediction));
                [~,ic]=min(abs(x-target));
            end
            idxCandidate=ic;
        end
        score=-inf(size(idxCandidate));
        rhoNow=rhoData(:,it);
        for jc=1:numel(idxCandidate)
            c=x(idxCandidate(jc));
            rhoRef=real(spectral_shift_collision(rhoTemplate0,k,c-xRef));
            rhoRef=max(0,rhoRef);
            w=soft_window_collision(x,c,ana.trackWindowHalfWidth, ...
                                    ana.trackTaperWidth);
            bc=real(sum(w.*sqrt(max(0,rhoRef.*rhoNow))));
            nr=real(sum(w.*rhoRef));
            nd=real(sum(w.*rhoNow));
            shapeScore=bc^2/(nr*nd+eps);
            motionCost=ana.trackMotionPenalty*((c-prediction)/ ...
                       max(ana.trackSearchHalfWidth,eps))^2;
            score(jc)=shapeScore-motionCost;
        end
        [bestScore,jbest]=max(score);
        center(it)=x(idxCandidate(jbest));
        quality(it)=bestScore;
    end
end

function q = local_wave_number_collision(u,k,dx,w)
    du=ifft(1i*k.*fft(u));
    den=real(dx*sum(w.*abs(u).^2));
    q=real(dx*sum(w.*imag(conj(u).*du)))/(den+eps);
end

function [Ng,Ne,conditionG] = mode_population_collision( ...
        u,x,k,dx,wUnion,wg,we,phiGround,phiExcited, ...
        xGroundRef,xExcitedRef,xg,xe,ana)
    phiG=spectral_shift_collision(phiGround,k,xg-xGroundRef);
    phiE=spectral_shift_collision(phiExcited,k,xe-xExcitedRef);
    if ana.removeTranslationalPhase
        qgData=local_wave_number_collision(u,k,dx,wg);
        qeData=local_wave_number_collision(u,k,dx,we);
        qgReference=local_wave_number_collision(phiG,k,dx,wg);
        qeReference=local_wave_number_collision(phiE,k,dx,we);
        phiG=phiG.*exp(1i*(qgData-qgReference)*(x-xg));
        phiE=phiE.*exp(1i*(qeData-qeReference)*(x-xe));
    end
    Phi=[phiG,phiE];
    G=dx*(Phi'*(wUnion.*Phi));
    G=0.5*(G+G');
    [Q,Lam]=eig(G);
    lam=real(diag(Lam));
    lamMax=max(max(lam),eps);
    lamSafe=max(lam,ana.gramEigenvalueFloor*lamMax);
    GminusHalf=Q*diag(1./sqrt(lamSafe))*Q';
    b=dx*(Phi'*(wUnion.*u));
    coeff=GminusHalf*b;
    Ng=real(abs(coeff(1))^2);
    Ne=real(abs(coeff(2))^2);
    conditionG=lamMax/max(min(lamSafe),eps);
end

function [className,info] = classify_collision_collision( ...
        t,Xg,Xe,overlap,iCollision,ana)
    Nt=numel(t);
    nPost=max(3,round(ana.postFraction*Nt));
    postIdx=max(iCollision+1,Nt-nPost+1):Nt;
    if isempty(postIdx), postIdx=max(1,Nt-nPost+1):Nt; end
    preEnd=max(1,iCollision-1);
    preStart=max(1,round(ana.ignoreInitialFraction*Nt));
    if preStart>preEnd, preStart=1; end
    preIdx=preStart:preEnd;
    signedSep=Xe-Xg;
    preSign=sign(median(signedSep(preIdx)));
    postSign=sign(median(signedSep(postIdx)));
    signChanged=(preSign~=0 && postSign~=0 && preSign~=postSign);
    postDistance=median(abs(signedSep(postIdx)));
    postOverlap=median(overlap(postIdx));

    if postDistance<ana.captureDistance && postOverlap>ana.captureOverlap
        className='CAPTURE / FUSION';
    elseif postDistance>ana.separatedDistance && signChanged
        className='TRANSMISSION / PASSING';
    elseif postDistance>ana.separatedDistance && ~signChanged
        className='REFLECTION / RE-SEPARATION';
    else
        className='PARTIAL CAPTURE / AMBIGUOUS';
    end
    info=struct('preSign',preSign,'postSign',postSign, ...
        'signChanged',signChanged,'postDistance',postDistance, ...
        'postOverlap',postOverlap,'postIndices',postIdx);
end

function [className,info] = classify_collision_survival_collision( ...
        t,Xg,Xe,massG,massE,validG,validE,iCollision,ana)
    Nt=numel(t);
    nPost=max(3,round(ana.postFraction*Nt));
    postIdx=max(iCollision+1,Nt-nPost+1):Nt;
    if isempty(postIdx), postIdx=max(1,Nt-nPost+1):Nt; end
    preStart=max(1,round(ana.ignoreInitialFraction*Nt));
    preEnd=max(preStart,iCollision-2);
    preIdx=preStart:preEnd;

    baseG=median(massG(preIdx));
    baseE=median(massE(preIdx));
    groundSurvival=median(massG(postIdx))/(baseG+eps);
    excitedSurvival=median(massE(postIdx))/(baseE+eps);
    groundValidFraction=mean(validG(postIdx));
    excitedValidFraction=mean(validE(postIdx));

    preBoth=preIdx(validG(preIdx)&validE(preIdx));
    postBoth=postIdx(validG(postIdx)&validE(postIdx));
    preSign=0; postSign=0; postDistance=nan;
    if ~isempty(preBoth)
        preSign=sign(median(Xe(preBoth)-Xg(preBoth)));
    end
    if ~isempty(postBoth)
        postSign=sign(median(Xe(postBoth)-Xg(postBoth)));
        postDistance=median(abs(Xe(postBoth)-Xg(postBoth)));
    end
    signChanged=preSign~=0 && postSign~=0 && preSign~=postSign;

    groundSurvives=groundSurvival>=ana.survivalThreshold && ...
                   groundValidFraction>=0.5;
    excitedSurvives=excitedSurvival>=ana.survivalThreshold && ...
                    excitedValidFraction>=0.5;
    if ~excitedSurvives && groundSurvival>=ana.captureReceiverThreshold
        className='CAPTURE: GROUND ABSORBS EXCITED';
    elseif ~groundSurvives && excitedSurvival>=ana.captureReceiverThreshold
        className='CAPTURE: EXCITED ABSORBS GROUND';
    elseif groundSurvives && excitedSurvives && signChanged
        className='TRANSMISSION / PASSING';
    elseif groundSurvives && excitedSurvives && ~signChanged
        className='REFLECTION / RE-SEPARATION';
    else
        className='PARTIAL CAPTURE / AMBIGUOUS';
    end
    info=struct('preSign',preSign,'postSign',postSign, ...
        'signChanged',signChanged,'postDistance',postDistance, ...
        'groundBaselineMass',baseG,'excitedBaselineMass',baseE, ...
        'groundSurvival',groundSurvival,'excitedSurvival',excitedSurvival, ...
        'groundValidFraction',groundValidFraction, ...
        'excitedValidFraction',excitedValidFraction, ...
        'preIndices',preIdx,'postIndices',postIdx);
end

function out = scan_phase_separation_collision( ...
        x,k,dx,phiG,phiE,xgRef,xeRef,lat,S,xc,tScan,Dselected,ana)
    if ~ismember(ana.phaseScanComponent,[1,2])
        error('ana.phaseScanComponent 必须为 1 或 2。');
    end
    if isfield(S,'job')
        job=S.job;
    else
        job=struct();
    end
    if isfield(job,'groundSpinor')
        chiG=job.groundSpinor(:);
    else
        chiG=[1;1]/sqrt(2);
    end
    if isfield(job,'excitedSpinor')
        chiE=job.excitedSpinor(:);
    else
        chiE=[1;1]/sqrt(2);
    end
    chiG=chiG/sqrt(sum(abs(chiG).^2)+eps);
    chiE=chiE/sqrt(sum(abs(chiE).^2)+eps);
    if isfield(job,'groundWeight'), weightG=job.groundWeight; else, weightG=1; end
    if isfield(job,'excitWeight'), weightE=job.excitWeight; else, weightE=1; end
    ampG=weightG*abs(chiG);
    ampE=weightE*abs(chiE);

    V=external_potential_collision(x,tScan,lat);
    wScan=soft_window_collision(x,xc,ana.scanWindowHalfWidth, ...
                                ana.scanWindowTaperWidth);
    Dgrid=ana.separationGrid(:).';
    phaseCases=ana.phaseCases(:).';
    U=zeros(numel(phaseCases),numel(Dgrid));
    Irho=zeros(size(U));
    peakRatio=zeros(size(U));

    for jd=1:numel(Dgrid)
        D=Dgrid(jd);
        [gD,eD]=placed_references_collision( ...
            phiG,phiE,k,xgRef,xeRef,xc,D);
        for jp=1:numel(phaseCases)
            [U(jp,jd),Irho(jp,jd),peakRatio(jp,jd)] = ...
                pair_metrics_collision(gD,eD,ampG,ampE,phaseCases(jp), ...
                ana.fixedRelativePhase,ana.phaseScanComponent, ...
                k,dx,V,lat,wScan);
        end
    end

    Dselected=max(min(Dselected,max(Dgrid)),min(Dgrid));
    phaseGrid=ana.phaseGrid(:).';
    Uphase=zeros(size(phaseGrid));
    Iphase=zeros(size(phaseGrid));
    Pphase=zeros(size(phaseGrid));
    [gSel,eSel]=placed_references_collision( ...
        phiG,phiE,k,xgRef,xeRef,xc,Dselected);
    for jp=1:numel(phaseGrid)
        [Uphase(jp),Iphase(jp),Pphase(jp)] = ...
            pair_metrics_collision(gSel,eSel,ampG,ampE,phaseGrid(jp), ...
            ana.fixedRelativePhase,ana.phaseScanComponent, ...
            k,dx,V,lat,wScan);
    end

    U0=interp1(phaseGrid,Uphase,0,'linear','extrap');
    Upi=interp1(phaseGrid,Uphase,pi,'linear','extrap');
    Jeff=0.5*(Upi-U0);
    if U0<Upi
        preference='Delta phi = 0 is energetically capture-favoring; pi is pass/repulsion-favoring';
    elseif Upi<U0
        preference='Delta phi = pi is energetically capture-favoring; 0 is pass/repulsion-favoring';
    else
        preference='no resolved phase preference';
    end

    out=struct();
    out.separation=Dgrid;
    out.phaseCases=phaseCases;
    out.interactionEnergy=U;
    out.densityInterference=Irho;
    out.peakDensityRatio=peakRatio;
    out.phaseGrid=phaseGrid;
    out.phaseInteractionEnergy=Uphase;
    out.phaseDensityInterference=Iphase;
    out.phasePeakDensityRatio=Pphase;
    out.selectedDistance=Dselected;
    out.scanCenter=xc;
    out.scanTime=tScan;
    out.U0=U0;
    out.Upi=Upi;
    out.Jeffective=Jeff;
    out.preference=preference;
end

function [gD,eD] = placed_references_collision( ...
        phiG,phiE,k,xgRef,xeRef,xc,D)
    xg=xc-D/2;
    xe=xc+D/2;
    gD=spectral_shift_collision(phiG,k,xg-xgRef);
    eD=spectral_shift_collision(phiE,k,xe-xeRef);
end

function [Uint,Irho,peakRatio] = pair_metrics_collision( ...
        gD,eD,ampG,ampE,phaseValue,fixedPhase,phaseComponent, ...
        k,dx,V,lat,w)
    delta=fixedPhase(:);
    if numel(delta)~=2, delta=[0;0]; end
    delta(phaseComponent)=phaseValue;

    uG1=ampG(1)*exp(1i*delta(1))*gD;
    uG2=ampG(2)*exp(1i*delta(2))*gD;
    uE1=ampE(1)*eD;
    uE2=ampE(2)*eD;
    u1=uG1+uE1;
    u2=uG2+uE2;

    Epair=two_component_energy_collision(u1,u2,k,dx,V,lat);
    Eg=two_component_energy_collision(uG1,uG2,k,dx,V,lat);
    Ee=two_component_energy_collision(uE1,uE2,k,dx,V,lat);
    Uint=real(Epair-Eg-Ee);

    rhoPair=abs(u1).^2+abs(u2).^2;
    rhoIncoherent=abs(uG1).^2+abs(uG2).^2+abs(uE1).^2+abs(uE2).^2;
    deltaRho=rhoPair-rhoIncoherent;
    % 用非相干密度加权，避免激发态正负波瓣的干涉贡献完全相消。
    Irho=real(dx*sum(w.*deltaRho.*rhoIncoherent))/ ...
         (real(dx*sum(w.*rhoIncoherent.^2))+eps);
    mask=w>0.5;
    if any(mask)
        peakRatio=max(rhoPair(mask))/(max(rhoIncoherent(mask))+eps);
    else
        peakRatio=max(rhoPair)/(max(rhoIncoherent)+eps);
    end
end

function E = two_component_energy_collision(u1,u2,k,dx,V,lat)
    du1=ifft(1i*k.*fft(u1));
    du2=ifft(1i*k.*fft(u2));
    rho1=abs(u1).^2;
    rho2=abs(u2).^2;
    density=0.5*abs(du1).^2+0.5*abs(du2).^2 ...
        +V.*(rho1+rho2)+0.5*lat.g11*rho1.^2 ...
        +0.5*lat.g22*rho2.^2+lat.g12*rho1.*rho2;
    E=real(dx*sum(density));
end

function label = phase_label_collision(phi)
    if abs(mod(phi,2*pi))<1e-10
        label='$\Delta\phi=0$';
    elseif abs(mod(phi,2*pi)-pi)<1e-10
        label='$\Delta\phi=\pi$';
    else
        label=sprintf('$\Delta\phi=%.3g\pi$',phi/pi);
    end
end

function s = safe_name_collision(s)
    s=regexprep(char(s),'[^a-zA-Z0-9_-]','_');
end
