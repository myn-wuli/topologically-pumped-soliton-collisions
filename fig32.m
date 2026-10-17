%% analyze_soliton_capture_energy.m
% 双组分泵浦碰撞的纯数据后处理：不重新演化波函数。
%
% 先运行原程序，使工作区中存在
%   x, tdata, udata1, udata2, lat
% 然后运行本程序。也可以在 cap.dataFile 中指定已保存的 MAT 文件。
%
% 本程序计算：
%   1. 小窗口内跟踪的两个孤子中心、速度和相对距离；
%   2. 局域密度重叠、交叉相互作用能和相互作用力；
%   3. 相对平动能、结合能、俘获参数和相对能量；
%   4. 完整 Gross--Pitaevskii 能量、泵浦做功及能量守恒误差；
%   5. 固定孤子形状时，结合能和吸引力随间距的变化。
load("g12=-5,G=10,E=10.mat");
%% ============================ 用户设置 ===============================
cap.dataFile = '';             % 空：直接使用当前工作区中的演化数据
cap.refFrame = 1;               % 构造 E_bind(D) 时使用的参考帧

% 中心位置跟踪参数。碰撞杂波较强时，优先减小 momentHalfWidth；若孤子
% 跃出跟踪区，再增大 searchHalfWidth。
cap.searchHalfWidth = 2.0;      % 在上一帧预测中心附近搜索孤子
cap.momentHalfWidth = 1.20;     % 真正用于中心矩积分的小窗口半宽
cap.centerIterations = 3;       % 小窗口中心迭代次数
cap.centerSmoothFrames = 9;     % 求速度前的移动平均窗口，建议取奇数

% 能量、重叠和相互作用力的局域窗口。采用余弦软边界，减小杂波和硬截断。
cap.coreHalfWidth = 1.60;
cap.coreTaperWidth = 0.25;

% 固定形状的距离扫描。
cap.distanceList = 0:0.05:8;
cap.collisionOverlapThreshold = 0.20; % 第一次超过此重叠度，定义为首次碰撞
cap.overlapThreshold = 0.50;
cap.incomingAverageFrames = 15;       % 碰撞前取中位数，定义入射相对动能
cap.maxCaptureRatioPlot = 1e8;        % 仅限制画图范围，不修改输出数据
cap.saveOutput = true;
cap.outputPrefix = 'soliton_capture_energy';

%% ============================== 读数据 ===============================
needLoad = ~(exist('x','var') == 1 && exist('tdata','var') == 1 && ...
             exist('udata1','var') == 1 && exist('udata2','var') == 1 && ...
             exist('lat','var') == 1);

if needLoad
    if isempty(cap.dataFile)
        error(['当前工作区缺少 x、tdata、udata1、udata2 或 lat。',newline, ...
               '请先运行原演化程序，或在 cap.dataFile 中填写结果 MAT 文件。']);
    end
    S = load(cap.dataFile,'x','tdata','udata1','udata2','lat');
    requiredNames = {'x','tdata','udata1','udata2','lat'};
    for jj = 1:numel(requiredNames)
        if ~isfield(S,requiredNames{jj})
            error('MAT 文件中缺少变量 %s。',requiredNames{jj});
        end
    end
    x = S.x;
    tdata = S.tdata;
    udata1 = S.udata1;
    udata2 = S.udata2;
    lat = S.lat;
end

x = x(:);
tdata = tdata(:).';
Nx = numel(x);
Nt = numel(tdata);
dx = x(2)-x(1);
Lx = Nx*dx;

if size(udata1,1) ~= Nx && size(udata1,2) == Nx
    udata1 = udata1.';
end
if size(udata2,1) ~= Nx && size(udata2,2) == Nx
    udata2 = udata2.';
end
if size(udata1,1) ~= Nx || size(udata2,1) ~= Nx
    error('udata1/udata2 的空间维必须等于 length(x)。');
end
if size(udata1,2) ~= Nt || size(udata2,2) ~= Nt
    error('udata1/udata2 的时间帧数必须等于 length(tdata)。');
end

requiredLat = {'p1','p2','d1','d2','nu','g11','g12','g22'};
for jj = 1:numel(requiredLat)
    if ~isfield(lat,requiredLat{jj})
        error('lat 中缺少参数 lat.%s。',requiredLat{jj});
    end
end

if mod(Nx,2)==0
    k = [0:Nx/2-1,-Nx/2:-1].'*(2*pi/Lx);
else
    k = [0:(Nx-1)/2,-(Nx-1)/2:-1].'*(2*pi/Lx);
end

cap.centerSmoothFrames = max(1,round(cap.centerSmoothFrames));
if mod(cap.centerSmoothFrames,2)==0
    cap.centerSmoothFrames = cap.centerSmoothFrames+1;
end

%% ======================== 小窗口跟踪孤子中心 =========================
rho10 = abs(udata1(:,1)).^2;
rho20 = abs(udata2(:,1)).^2;
x10 = sum(x.*rho10)/(sum(rho10)+eps);
x20 = sum(x.*rho20)/(sum(rho20)+eps);

[X1raw,~] = track_local_center(x,abs(udata1).^2,dx,x10,cap);
[X2raw,~] = track_local_center(x,abs(udata2).^2,dx,x20,cap);

X1 = movmean(X1raw,cap.centerSmoothFrames,'Endpoints','shrink');
X2 = movmean(X2raw,cap.centerSmoothFrames,'Endpoints','shrink');
v1 = gradient(X1,tdata);
v2 = gradient(X2,tdata);
Dactual = abs(X2-X1);

fprintf('\nLocal center tracking:\n');
fprintf('  center search half width = %.6g\n',cap.searchHalfWidth);
fprintf('  center moment half width = %.6g\n',cap.momentHalfWidth);
fprintf('  smoothing frames         = %d\n',cap.centerSmoothFrames);

%% ================= 相互作用、能量和泵浦做功 =========================
N1core = zeros(1,Nt);
N2core = zeros(1,Nt);
overlapCore = zeros(1,Nt);
E12core = zeros(1,Nt);
E12full = zeros(1,Nt);
F12attract = zeros(1,Nt);
F12full = zeros(1,Nt);

Ekin1 = zeros(1,Nt); Ekin2 = zeros(1,Nt);
Eself1 = zeros(1,Nt); Eself2 = zeros(1,Nt);
Eext1 = zeros(1,Nt); Eext2 = zeros(1,Nt);
Etot = zeros(1,Nt);
pumpPower = zeros(1,Nt);

for it = 1:Nt
    u1 = udata1(:,it);
    u2 = udata2(:,it);
    rho1 = abs(u1).^2;
    rho2 = abs(u2).^2;

    w1 = soft_core_window(x,X1(it),cap.coreHalfWidth,cap.coreTaperWidth);
    w2 = soft_core_window(x,X2(it),cap.coreHalfWidth,cap.coreTaperWidth);
    wUnion = max(w1,w2);

    N1core(it) = real(dx*sum(w1.*rho1));
    N2core(it) = real(dx*sum(w2.*rho2));

    numerator = real(dx*sum(wUnion.*rho1.*rho2));
    denominator = sqrt(real(dx*sum(wUnion.*rho1.^2))* ...
                       real(dx*sum(wUnion.*rho2.^2)) + eps);
    overlapCore(it) = min(1,max(0,numerator/(denominator+eps)));

    E12core(it) = real(lat.g12*dx*sum(wUnion.*rho1.*rho2));
    E12full(it) = real(lat.g12*dx*sum(rho1.*rho2));

    drho1 = real(ifft(1i*k.*fft(rho1)));
    drho2 = real(ifft(1i*k.*fft(rho2)));
    F1core = real(-lat.g12*dx*sum(wUnion.*rho1.*drho2));
    F2core = real(-lat.g12*dx*sum(wUnion.*rho2.*drho1));
    FpairCore = 0.5*(F1core-F2core);

    F1all = real(-lat.g12*dx*sum(rho1.*drho2));
    F2all = real(-lat.g12*dx*sum(rho2.*drho1));
    FpairAll = 0.5*(F1all-F2all);

    sepSign = sign(X2(it)-X1(it));
    if sepSign==0
        sepSign = 1;
    end
    % 正值表示使 D 减小的吸引力，负值表示瞬时排斥方向。
    F12attract(it) = sepSign*FpairCore;
    F12full(it) = sepSign*FpairAll;

    du1 = ifft(1i*k.*fft(u1));
    du2 = ifft(1i*k.*fft(u2));
    Vext = external_potential_capture(x,tdata(it),lat);
    dVdt = external_potential_time_derivative(x,tdata(it),lat);

    Ekin1(it) = real(0.5*dx*sum(abs(du1).^2));
    Ekin2(it) = real(0.5*dx*sum(abs(du2).^2));
    Eself1(it) = real(0.5*lat.g11*dx*sum(rho1.^2));
    Eself2(it) = real(0.5*lat.g22*dx*sum(rho2.^2));
    Eext1(it) = real(dx*sum(Vext.*rho1));
    Eext2(it) = real(dx*sum(Vext.*rho2));

    Etot(it) = Ekin1(it)+Ekin2(it)+Eself1(it)+Eself2(it)+ ...
               Eext1(it)+Eext2(it)+E12full(it);
    pumpPower(it) = real(dx*sum(dVdt.*(rho1+rho2)));
end

% 用局域孤子粒子数作为平动质量，避免辐射背景被计入相对动能。
reducedMass = N1core.*N2core./(N1core+N2core+eps);
Krel = 0.5*reducedMass.*(v2-v1).^2;
EbindDirect = -E12core;
ErelDirect = Krel+E12core;
etaCapture = EbindDirect./(Krel+1e-14);

pumpWork = cumtrapz(tdata,pumpPower);
energyBalanceError = (Etot-Etot(1))-pumpWork;

%% ================ 固定孤子形状：E_bind(D) 和 F(D) ==================
iref = min(max(round(cap.refFrame),1),Nt);
rho1ref = abs(udata1(:,iref)).^2;
rho2ref = abs(udata2(:,iref)).^2;
w1ref = soft_core_window(x,X1(iref),cap.coreHalfWidth,cap.coreTaperWidth);
w2ref = soft_core_window(x,X2(iref),cap.coreHalfWidth,cap.coreTaperWidth);
rho1ref = rho1ref.*w1ref;
rho2ref = rho2ref.*w2ref;

Dscan = cap.distanceList(:);
if any(Dscan<0) || numel(Dscan)<3
    error('cap.distanceList 必须包含至少三个非负、递增的距离。');
end
if any(diff(Dscan)<=0)
    error('cap.distanceList 必须严格递增。');
end

U12scan = zeros(size(Dscan));
for id = 1:numel(Dscan)
    targetX2 = X1(iref)+Dscan(id);
    rho2shift = shift_real_profile(rho2ref,k,targetX2-X2(iref));
    U12scan(id) = real(lat.g12*dx*sum(rho1ref.*rho2shift));
end

% 用最大距离处的相互作用能作为分离态基准，避免有限计算区中的微小尾部误差。
nBase = max(1,round(0.05*numel(Dscan)));
U12separated = mean(U12scan(end-nBase+1:end));
EbindScan = U12separated-U12scan;
generalizedForce = -gradient(U12scan,Dscan); % 对相对坐标 D 的力，吸引时为负
attractiveForceScan = abs(generalizedForce);

EbindModel = interp_inside(Dscan,EbindScan,Dactual);
ErelModel = Krel-EbindModel;

%% ============================= 碰撞指标 =============================
[Dmin,iclosest] = min(Dactual);
collisionCandidate = find(overlapCore>=cap.collisionOverlapThreshold,1,'first');
if isempty(collisionCandidate)
    icollision = iclosest;
else
    icollision = collisionCandidate;
end
tcollision = tdata(icollision);

fprintf('\nCollision/capture diagnostics:\n');
fprintf('  first collision (overlap >= %.3f): t = %.8g, D = %.8g\n', ...
        cap.collisionOverlapThreshold,tcollision,Dactual(icollision));
fprintf('  closest approach: Dmin = %.8g at t = %.8g\n', ...
        Dmin,tdata(iclosest));
fprintf('  localized density overlap at first collision = %.8g\n', ...
        overlapCore(icollision));
fprintf('  Ebind(core) at first collision = %.8g\n',EbindDirect(icollision));
fprintf('  Krel at first collision        = %.8g\n',Krel(icollision));
fprintf('  eta=Ebind/Krel                  = %.8g\n',etaCapture(icollision));
fprintf('  Erel=Krel+E12                   = %.8g\n',ErelDirect(icollision));

postStart = min(Nt,icollision+cap.centerSmoothFrames);
postMask = postStart:Nt;
if ~isempty(postMask)
    negativeFraction = mean(ErelDirect(postMask)<0);
    lockedFraction = mean(overlapCore(postMask)>=cap.overlapThreshold);
    fprintf('  post-collision fraction with Erel<0       = %.6f\n',negativeFraction);
    fprintf('  post-collision fraction with overlap>=%.2f = %.6f\n', ...
            cap.overlapThreshold,lockedFraction);
else
    negativeFraction = nan;
    lockedFraction = nan;
end

fprintf('  max |energy balance error| = %.6e\n',max(abs(energyBalanceError)));

%% ============== 俘获系数与孤子间距的关系 ===========================
% 固定形状模型中，以碰撞前的入射相对动能作为统一分母：
%       eta_model(D) = Ebind(D)/Krel_in.
% 这样 eta_model 是 D 的单值函数；eta_model=1 给出近似俘获临界距离。
preEnd = max(1,icollision-cap.centerSmoothFrames);
preStart = max(1,preEnd-round(cap.incomingAverageFrames)+1);
preRange = preStart:preEnd;
validIncoming = isfinite(Krel(preRange)) & Krel(preRange)>1e-14;
if any(validIncoming)
    KrelIncoming = median(Krel(preRange(validIncoming)));
else
    positiveK = Krel(isfinite(Krel) & Krel>1e-14);
    if isempty(positiveK)
        KrelIncoming = nan;
    else
        KrelIncoming = median(positiveK);
    end
end

if isfinite(KrelIncoming) && KrelIncoming>0
    etaDistanceModel = max(EbindScan,0)/KrelIncoming;
else
    etaDistanceModel = nan(size(EbindScan));
end

captureDistanceMask = isfinite(etaDistanceModel) & etaDistanceModel>=1;
if any(captureDistanceMask)
    criticalCaptureDistance = max(Dscan(captureDistanceMask));
else
    criticalCaptureDistance = nan;
end

fprintf('\nCapture coefficient versus distance:\n');
fprintf('  incoming Krel (median over frames %d:%d) = %.8g\n', ...
        preStart,preEnd,KrelIncoming);
if isfinite(criticalCaptureDistance)
    fprintf('  eta_model(D)=1 gives approximate critical distance Dc = %.8g\n', ...
            criticalCaptureDistance);
else
    fprintf('  eta_model(D) does not reach 1 in the scanned distance range.\n');
end

%% ================================ 作图 ===============================
if abs(lat.nu)>0
    pumpPeriod = pi/abs(lat.nu);
    tplot = tdata/pumpPeriod;
    timeLabel = '$t/T$';
else
    pumpPeriod = nan;
    tplot = tdata;
    timeLabel = '$t$';
end

figure('Color','w','Name','Quantitative soliton capture');
tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

nexttile; box on;
yyaxis left
plot(tplot,Dactual,'k-','LineWidth',1.5,'DisplayName','$D(t)$');
ylabel('$D(t)$','Interpreter','latex');
yyaxis right
plot(tplot,overlapCore,'Color',[0.55 0.10 0.70],'LineWidth',1.35, ...
     'DisplayName','$\Lambda_{\rm core}(t)$');
ylabel('$\Lambda_{\rm core}$','Interpreter','latex');
ylim([0,1.02]);
xline(tplot(icollision),'k:','LineWidth',1.1,'HandleVisibility','off');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);

nexttile; hold on; box on;
plot(tplot,EbindDirect,'LineWidth',1.5,'DisplayName','$E_{\rm bind}=-E_{12}$');
plot(tplot,Krel,'LineWidth',1.5,'DisplayName','$K_{\rm rel}$');
xline(tplot(icollision),'k:','LineWidth',1.1,'HandleVisibility','off');
ylabel('energy');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);

nexttile; box on;
etaPlot = min(max(etaCapture,1e-12),cap.maxCaptureRatioPlot);
semilogy(tplot,etaPlot,'Color',[0.1 0.55 0.25],'LineWidth',1.5);
yline(1,'k:','LineWidth',1.1);
xline(tplot(icollision),'k:','LineWidth',1.1,'HandleVisibility','off');
ylabel('$\eta=E_{\rm bind}/K_{\rm rel}$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);

nexttile; hold on; box on;
plot(tplot,ErelDirect,'LineWidth',1.5,'DisplayName','$K_{\rm rel}+E_{12}^{\rm core}$');
plot(tplot,ErelModel,'--','LineWidth',1.3,'DisplayName','$K_{\rm rel}-E_{\rm bind}(D)$');
yline(0,'k:','LineWidth',1.1,'HandleVisibility','off');
xline(tplot(icollision),'k:','LineWidth',1.1,'HandleVisibility','off');
ylabel('$E_{\rm rel}$','Interpreter','latex');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);

nexttile; hold on; box on;
plot(tplot,F12attract,'LineWidth',1.5,'DisplayName','core-window force');
plot(tplot,F12full,'--','LineWidth',1.1,'DisplayName','full-domain force');
yline(0,'k:','LineWidth',1.0,'HandleVisibility','off');
xline(tplot(icollision),'k:','LineWidth',1.1,'HandleVisibility','off');
xlabel(timeLabel,'Interpreter','latex');
ylabel('$F_{12}^{\rm attr}$','Interpreter','latex');
legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);

figure('Color','w','Name','Binding energy and attraction versus distance');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; box on;
plot(Dscan,EbindScan,'LineWidth',1.6);
xlabel('$D$','Interpreter','latex');
ylabel('$E_{\rm bind}(D)$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.0);
nexttile; box on;
plot(Dscan,attractiveForceScan,'LineWidth',1.6);
xlabel('$D$','Interpreter','latex');
ylabel('$|F_{12}(D)|$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.0);

figure('Color','w','Name','Capture coefficient versus soliton distance');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; box on; hold on;
etaModelPlot = min(max(etaDistanceModel,1e-12),cap.maxCaptureRatioPlot);
semilogy(Dscan,etaModelPlot,'Color',[0.10 0.45 0.80],'LineWidth',1.7);
yline(1,'k:','LineWidth',1.2,'DisplayName','$\eta_{\rm model}=1$');
if isfinite(criticalCaptureDistance)
    xline(criticalCaptureDistance,'--','Color',[0.75 0.15 0.10], ...
          'LineWidth',1.3,'DisplayName','$D_{\rm c}$');
end
xlabel('$D$','Interpreter','latex');
ylabel('$\eta_{\rm model}(D)=E_{\rm bind}(D)/K_{\rm rel}^{\rm in}$', ...
       'Interpreter','latex');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.0);

nexttile; box on; hold on;
etaActualPlot = min(max(etaCapture,1e-12),cap.maxCaptureRatioPlot);
approachIndex = 1:icollision;
capturedIndex = icollision:Nt;
semilogy(Dactual(approachIndex),etaActualPlot(approachIndex),'o-', ...
    'Color',[0.10 0.45 0.80],'MarkerSize',3.5,'LineWidth',1.1, ...
    'DisplayName','approach');
semilogy(Dactual(capturedIndex),etaActualPlot(capturedIndex),'.-', ...
    'Color',[0.85 0.25 0.10],'MarkerSize',7,'LineWidth',1.0, ...
    'DisplayName','after collision');
yline(1,'k:','LineWidth',1.2,'HandleVisibility','off');
xlabel('$D(t)=|X_2(t)-X_1(t)|$','Interpreter','latex');
ylabel('$\eta(t)=E_{\rm bind}(t)/K_{\rm rel}(t)$','Interpreter','latex');
legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.0);

figure('Color','w','Name','Energy conversion and pump work');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
plot(tplot,Ekin1+Ekin2-(Ekin1(1)+Ekin2(1)),'LineWidth',1.2, ...
     'DisplayName','$\Delta E_{\rm kin}$');
plot(tplot,Eself1+Eself2-(Eself1(1)+Eself2(1)),'LineWidth',1.2, ...
     'DisplayName','$\Delta E_{\rm self}$');
plot(tplot,Eext1+Eext2-(Eext1(1)+Eext2(1)),'LineWidth',1.2, ...
     'DisplayName','$\Delta E_{\rm ext}$');
plot(tplot,E12full-E12full(1),'LineWidth',1.4,'DisplayName','$\Delta E_{12}$');
ylabel('energy change');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1.0);

nexttile; hold on; box on;
plot(tplot,Etot-Etot(1),'LineWidth',1.5,'DisplayName','$E(t)-E(0)$');
plot(tplot,pumpWork,'--','LineWidth',1.5,'DisplayName','$W_{\rm pump}(t)$');
ylabel('total energy / work');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1.0);

nexttile; box on;
plot(tplot,energyBalanceError,'k','LineWidth',1.3);
xlabel(timeLabel,'Interpreter','latex');
ylabel('$E-E_0-W_{\rm pump}$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1.0);

figure('Color','w','Name','Local center tracking check');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
plot(tplot,X1raw,':','Color',[0.2 0.5 0.9],'DisplayName','$X_{1,\rm raw}$');
plot(tplot,X1,'-','Color',[0.0 0.25 0.75],'LineWidth',1.4,'DisplayName','$X_1$');
plot(tplot,X2raw,':','Color',[0.9 0.45 0.2],'DisplayName','$X_{2,\rm raw}$');
plot(tplot,X2,'-','Color',[0.75 0.15 0.0],'LineWidth',1.4,'DisplayName','$X_2$');
ylabel('local center');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1.0);
nexttile; hold on; box on;
plot(tplot,N1core,'LineWidth',1.3,'DisplayName','$N_{1,\rm core}$');
plot(tplot,N2core,'LineWidth',1.3,'DisplayName','$N_{2,\rm core}$');
xlabel(timeLabel,'Interpreter','latex');
ylabel('localized norm');
legend('Location','best','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',1.0);

%% =============================== 输出 ================================
TimeTable = table(tdata(:),X1raw(:),X2raw(:),X1(:),X2(:),v1(:),v2(:), ...
    Dactual(:),N1core(:),N2core(:),overlapCore(:),E12core(:),EbindDirect(:), ...
    Krel(:),ErelDirect(:),etaCapture(:),F12attract(:),Etot(:),pumpWork(:), ...
    energyBalanceError(:), ...
    'VariableNames',{'Time','X1raw','X2raw','X1','X2','Velocity1','Velocity2', ...
    'Distance','N1core','N2core','CoreOverlap','E12core','BindingEnergy', ...
    'RelativeKineticEnergy','RelativeEnergy','CaptureRatio','AttractiveForce', ...
    'TotalEnergy','PumpWork','EnergyBalanceError'});

DistanceTable = table(Dscan,U12scan,EbindScan,generalizedForce,attractiveForceScan, ...
    etaDistanceModel, ...
    'VariableNames',{'Distance','CrossInteractionEnergy','BindingEnergy', ...
    'GeneralizedForce','AttractiveForceMagnitude','CaptureRatioModel'});

captureInfo = struct();
captureInfo.collisionIndex = icollision;
captureInfo.collisionTime = tcollision;
captureInfo.closestApproachIndex = iclosest;
captureInfo.closestApproachTime = tdata(iclosest);
captureInfo.minimumDistance = Dmin;
captureInfo.captureRatioAtCollision = etaCapture(icollision);
captureInfo.incomingRelativeKineticEnergy = KrelIncoming;
captureInfo.criticalCaptureDistance = criticalCaptureDistance;
captureInfo.relativeEnergyAtCollision = ErelDirect(icollision);
captureInfo.negativeRelativeEnergyFractionAfterCollision = negativeFraction;
captureInfo.lockedOverlapFractionAfterCollision = lockedFraction;
captureInfo.referenceFrame = iref;
captureInfo.referenceTime = tdata(iref);
captureInfo.pumpPeriod = pumpPeriod;

% if cap.saveOutput
%     writetable(TimeTable,[cap.outputPrefix '_time.csv']);
%     writetable(DistanceTable,[cap.outputPrefix '_distance.csv']);
%     save([cap.outputPrefix '.mat'],'TimeTable','DistanceTable','captureInfo','cap', ...
%          'X1','X2','v1','v2','Dactual','overlapCore','E12core','EbindDirect', ...
%          'Krel','ErelDirect','etaCapture','F12attract','Etot','pumpWork', ...
%          'energyBalanceError','Dscan','U12scan','EbindScan','attractiveForceScan', ...
%          'etaDistanceModel','KrelIncoming','criticalCaptureDistance');
%     fprintf('\nSaved: %s_time.csv, %s_distance.csv and %s.mat\n', ...
%             cap.outputPrefix,cap.outputPrefix,cap.outputPrefix);
% end

%% ====================================================================
%% 局部函数
%% ====================================================================
function [center,normLocal] = track_local_center(x,rhoData,dx,xInitial,cap)
    Nt = size(rhoData,2);
    center = zeros(1,Nt);
    normLocal = zeros(1,Nt);
    center(1) = xInitial;

    for it = 1:Nt
        if it==1
            xPredict = center(1);
        elseif it==2
            xPredict = center(it-1);
        else
            stepPredict = center(it-1)-center(it-2);
            stepPredict = max(-0.5*cap.searchHalfWidth, ...
                              min(0.5*cap.searchHalfWidth,stepPredict));
            xPredict = center(it-1)+stepPredict;
        end

        rho = rhoData(:,it);
        searchMask = abs(x-xPredict)<=cap.searchHalfWidth;
        if sum(rho(searchMask))*dx<=1e-14
            center(it) = xPredict;
            normLocal(it) = 0;
            continue;
        end

        xNow = sum(x(searchMask).*rho(searchMask))/(sum(rho(searchMask))+eps);
        for jj = 1:cap.centerIterations
            momentMask = abs(x-xNow)<=cap.momentHalfWidth;
            localNorm = sum(rho(momentMask));
            if localNorm<=1e-14
                break;
            end
            xNow = sum(x(momentMask).*rho(momentMask))/(localNorm+eps);
        end
        center(it) = xNow;
        finalMask = abs(x-center(it))<=cap.momentHalfWidth;
        normLocal(it) = real(dx*sum(rho(finalMask)));
    end
end

function w = soft_core_window(x,c,halfWidth,taperWidth)
    taperWidth = min(max(taperWidth,0),halfWidth);
    r = abs(x-c);
    w = zeros(size(x));
    inner = max(0,halfWidth-taperWidth);
    w(r<=inner) = 1;
    if taperWidth>0
        trans = r>inner & r<halfWidth;
        w(trans) = 0.5*(1+cos(pi*(r(trans)-inner)/taperWidth));
    end
end

function V = external_potential_capture(x,t,lat)
    V = -lat.p1*cos(pi*x/lat.d1).^2 ...
        -lat.p2*cos(pi*x/lat.d2-lat.nu*t).^2;
end

function dVdt = external_potential_time_derivative(x,t,lat)
    phase = pi*x/lat.d2-lat.nu*t;
    dVdt = -lat.p2*lat.nu*sin(2*phase);
end

function shifted = shift_real_profile(profile,k,delta)
    shifted = real(ifft(fft(profile(:)).*exp(-1i*k*delta)));
    % 傅里叶平移在软截断边缘可能产生极小的负振铃，密度必须非负。
    shifted(shifted<0) = 0;
end

function yq = interp_inside(x,y,xq)
    yq = nan(size(xq));
    inside = xq>=min(x) & xq<=max(x);
    yq(inside) = interp1(x,y,xq(inside),'pchip');
end
