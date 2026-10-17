%% 相对相位 0：碰撞后单基态轨迹 -> 移动窗口 -> 直接统计粒子数
% 先 load 相对相位为 0 的原演化数据，再直接运行本脚本。
% 也可填写 ana.dataFile；本程序只做统计，不修改原始演化或相位。
% 只删去一次连续的碰撞过程，其他保存帧全部积分。
% 碰撞前分别统计 G/E；碰撞后按指定口径将剩余波包计入 G，E 记为 0。
% G 始终积分轨迹窗口内的实际密度，窗外辐射另计；G/E 是波包统计标签。
% 第一幅图显示实际轨迹（实线）、积分窗口（半透明色带）及左右边界（粗点线）。
% 第二幅图显示四条粒子数曲线。计算、保存工作区后统一绘图。

%% 设置
ana=struct();
ana.dataFile='';
ana.caseName='';
ana.windowG=[1.0,1.0];      % 基态轨迹左、右两侧的窗口宽度
ana.windowE=[1.5,1.0];      % 激发态轨迹左、右两侧的窗口宽度
ana.collisionRange=[1.4,2.0]; % 按当前轨迹图设定碰撞段 [开始,结束]，单位 t/T
ana.postSearchHalfWidth=1.0; % 碰撞后逐帧搜索主波包的半宽
ana.postCenterHalfWidth=0.75; % 主峰附近用于确定轨迹的局部质心半宽
ana.postStartX=[];          % 空：从碰撞后首帧最强主峰起追；也可指定起始 x
ana.densitySmoothSigma=0.15; % 仅用于定位；粒子数始终积分原始密度
ana.saveOutput=true;
ana.outputPrefix='N0';     % 相位 0 的输出，避免覆盖 pi 的分析结果

% 以下参数仅用于在真实密度中定位两条轨迹，不参与粒子数积分。
ana.trackReacquireHalfWidth=3.0;
ana.trackMotionPenalty=1e-4;
ana.trackWindowHalfWidth=1.50;
ana.trackTaperWidth=0.25;
ana.referenceDt=[];
ana.referenceUseFourthOrder=[];
ana.actualGroundExtraPhase=0;
ana.reuseReferenceCache=true;
ana.referenceCacheFile='ref_p0.mat';

%% 读取原始数据
if isempty(ana.dataFile)
    needed = {'x','tdata','udata1','udata2','lat','psiGround'};
    S = struct();
    for jj = 1:numel(needed)
        if exist(needed{jj},'var')~=1
            error('缺少 %s：请先 load 原演化数据，或设置 ana.dataFile。',needed{jj});
        end
        S.(needed{jj}) = eval(needed{jj});
    end
    if exist('psi','var')==1, S.psi=psi; end
    if exist('psiExcited','var')==1, S.psiExcited=psiExcited; end
    if exist('job','var')==1, S.job=job; end
    if exist('pump','var')==1, S.pump=pump; end
    if exist('k','var')==1, S.k=k; end
    if exist('vphi','var')==1, S.vphi=vphi; end
    if exist('isolatedReference','var')==1 && exist('referenceCacheSignature','var')==1
        S.isolatedReference=isolatedReference;
        S.referenceCacheSignature=referenceCacheSignature;
    end
else
    S = load(ana.dataFile);
end
% 核对保存的相位；不能把 pi 的演化数据通过后处理改成相位 0。
phaseUsed=0; phaseSource='user specified phase 0';
if isfield(S,'vphi')
    if ~isscalar(S.vphi) || ~isreal(S.vphi) || ~isfinite(S.vphi) || ...
            abs(angle(exp(1i*S.vphi)))>1e-8
        error('此程序用于相对相位 0，请加载对应的演化数据；当前 vphi 不是 0 (mod 2*pi)。');
    end
    phaseSource='saved vphi';
end
if ~isfield(S,'psi') && isfield(S,'psiExcited'), S.psi=S.psiExcited; end
required = {'x','tdata','udata1','udata2','lat','psiGround','psi'};
for jj=1:numel(required)
    if ~isfield(S,required{jj}), error('数据中缺少 %s。',required{jj}); end
end
x=S.x(:); tdata=S.tdata(:).'; lat=S.lat;
phiGround=S.psiGround(:); phiExcited=S.psi(:);
u1Data=S.udata1; u2Data=S.udata2;
Nx=numel(x); Nt=numel(tdata);
if Nx<3 || Nt<2, error('空间网格/保存帧太少。'); end
dx=x(2)-x(1); Lx=Nx*dx;
if any(~isfinite(x)) || any(~isfinite(tdata)) || dx<=0 || ...
        max(abs(diff(x)-dx))>1e-8*dx || any(diff(tdata)<=0)
    error('x 必须为递增等距网格，tdata 必须严格递增。');
end
if size(u1Data,1)~=Nx && size(u1Data,2)==Nx, u1Data=u1Data.'; end
if size(u2Data,1)~=Nx && size(u2Data,2)==Nx, u2Data=u2Data.'; end
if ~isequal(size(u1Data),[Nx,Nt]) || ~isequal(size(u2Data),[Nx,Nt]) || ...
        numel(phiGround)~=Nx || numel(phiExcited)~=Nx
    error('波函数尺寸应为 length(x) x length(tdata)，参考态长度应为 length(x)。');
end
if any(~isfinite(u1Data(:))) || any(~isfinite(u2Data(:)))
    error('原始波函数包含 NaN/Inf。');
end
if isfield(S,'k') && numel(S.k)==Nx
    k=S.k(:);
else
    k=[0:ceil(Nx/2)-1,-floor(Nx/2):-1].'*(2*pi/Lx);
end
if isempty(ana.caseName)
    ana.caseName=sprintf('g%g_G%g_E%g',lat.g12, ...
        max(abs(phiGround))^2,max(abs(phiExcited))^2);
end
xGroundRef=density_center_collision(x,phiGround,dx);
xExcitedRef=density_center_collision(x,phiExcited,dx);
rhoTotal=abs(u1Data).^2+abs(u2Data).^2;

%% 明确碰撞前、碰撞段与碰撞后，不用双峰判据删除整个碰撞后时段
if isfield(lat,'nu') && abs(lat.nu)>0
    pumpPeriod=pi/abs(lat.nu); tPlot=tdata/pumpPeriod; timeLabel='t/T';
else
    pumpPeriod=nan; tPlot=tdata; timeLabel='t';
end
if numel(ana.collisionRange)~=2 || any(~isfinite(ana.collisionRange)) || ...
        ana.collisionRange(1)>=ana.collisionRange(2)
    error('ana.collisionRange 应为有限递增的 [开始,结束]，单位与时间轴一致。');
end
beforeCollision=tPlot<ana.collisionRange(1);
afterCollision=tPlot>ana.collisionRange(2);
collisionMask=~beforeCollision & ~afterCollision;
keep=~collisionMask;
if ~any(beforeCollision) || ~any(afterCollision)
    error('数据需同时包含碰撞前、后的保存帧，请调整 ana.collisionRange。');
end
collisionRangeUsed=ana.collisionRange;
fprintf('碰撞段：%s = %g 到 %g；碰撞后仅追踪 G，E 记为 0。\n', ...
    timeLabel,collisionRangeUsed(1),collisionRangeUsed(2));

%% 动态参考态：严格核对网格、保存时刻、初态和参数，避免误用旧缓存
% 双波包参考只用于碰撞前；合并后直接搜索实际密度，不再演化参考态。
referenceTimes=tdata(beforeCollision);
referenceCacheSignatureNew=reference_signature_collision( ...
    x,k,referenceTimes,phiGround,phiExcited,lat,S,ana);
loadedReferenceCache=false;
    if ana.reuseReferenceCache && isfield(S,'isolatedReference') && ...
            isfield(S,'referenceCacheSignature') && ...
            isequaln(S.referenceCacheSignature,referenceCacheSignatureNew)
        isolatedReference=S.isolatedReference;
        loadedReferenceCache=true;
    elseif ana.reuseReferenceCache && isfile(ana.referenceCacheFile)
        cacheData=load(ana.referenceCacheFile,'isolatedReference','referenceCacheSignature');
        if isfield(cacheData,'isolatedReference') && isfield(cacheData,'referenceCacheSignature') && ...
                isequaln(cacheData.referenceCacheSignature,referenceCacheSignatureNew)
            isolatedReference=cacheData.isolatedReference;
            loadedReferenceCache=true;
        end
        clear cacheData
    end
    if ~loadedReferenceCache
        isolatedReference=evolve_isolated_references_collision( ...
            x,k,referenceTimes,phiGround,phiExcited,lat,S,ana);
    end
referenceCacheSignature=referenceCacheSignatureNew;
XgroundIsolated=isolatedReference.Xground(:).';
XexcitedIsolated=isolatedReference.Xexcited(:).';
rhoGroundReferenceData=abs(isolatedReference.g1).^2+abs(isolatedReference.g2).^2;
rhoExcitedReferenceData=abs(isolatedReference.e1).^2+abs(isolatedReference.e2).^2;
%% 碰撞前追踪两个波包，碰撞后只沿原始密度追踪一个主波包
Xground=nan(1,Nt); Xexcited=Xground; trackFitError=Xground;
[Xground(beforeCollision),Xexcited(beforeCollision),trackFitError(beforeCollision)]= ...
    track_pair_collision(x,k,rhoTotal(:,beforeCollision), ...
    rhoGroundReferenceData,rhoExcitedReferenceData, ...
    XgroundIsolated,XexcitedIsolated,ana);
radius=max(1,ceil(4*ana.densitySmoothSigma/dx));
z=(-radius:radius).'*dx;
kernel=exp(-z.^2/(2*ana.densitySmoothSigma^2)); kernel=kernel/sum(kernel);
rhoLocate=conv2(rhoTotal,kernel,'same');
Xground(afterCollision)=track_single_collision(x,tdata(afterCollision), ...
    rhoTotal(:,afterCollision),rhoLocate(:,afterCollision),ana);
XgroundRaw=Xground; XexcitedRaw=Xexcited;

%% 再设置轨迹左右窗口，并对原始密度直接积分
N1ground=nan(1,Nt); N2ground=N1ground; N1excited=N1ground; N2excited=N1ground;
XgroundLeft=N1ground; XgroundRight=N1ground;
XexcitedLeft=N1ground; XexcitedRight=N1ground;
windowG=zeros(Nx,Nt); windowE=zeros(Nx,Nt);
if numel(ana.windowG)~=2 || numel(ana.windowE)~=2 || ...
        any(~isfinite([ana.windowG(:);ana.windowE(:)])) || ...
        any([ana.windowG(:);ana.windowE(:)]<=0)
    error('ana.windowG、ana.windowE 必须各含两个正数：[左宽,右宽]。');
end
for it=find(keep(:)).'
    if afterCollision(it)
        % 碰撞后的唯一窗口随真实主波包移动；不再用虚构的 E 轨迹切割 G。
        left=max(x(1)-dx/2,Xground(it)-ana.windowG(1));
        right=min(x(end)+dx/2,Xground(it)+ana.windowG(2));
        XgroundLeft(it)=left; XgroundRight(it)=right;
        windowG(:,it)=max(0,min(x+dx/2,right)-max(x-dx/2,left))/dx;
        N1ground(it)=dx*sum(windowG(:,it).*abs(u1Data(:,it)).^2);
        N2ground(it)=dx*sum(windowG(:,it).*abs(u2Data(:,it)).^2);
        N1excited(it)=0; N2excited(it)=0;
        continue;
    end
    limits=[Xground(it)-ana.windowG(1),Xground(it)+ana.windowG(2); ...
            Xexcited(it)-ana.windowE(1),Xexcited(it)+ana.windowE(2)];
    limits(:,1)=max(limits(:,1),x(1)-dx/2);
    limits(:,2)=min(limits(:,2),x(end)+dx/2);
    % 宽窗口在靠近阶段若相交，在两条实际轨迹的中点切开，避免重复计数。
    left=max(limits(:,1)); right=min(limits(:,2));
    if left<right
        cut=min(max(0.5*(Xground(it)+Xexcited(it)),left),right);
        if Xground(it)<Xexcited(it)
            limits(1,2)=cut; limits(2,1)=cut;
        else
            limits(2,2)=cut; limits(1,1)=cut;
        end
    end
    if any(limits(:,2)<=limits(:,1))
        error('第 %d 帧无法分开积分窗口，请检查 ana.collisionRange。',it);
    end
    XgroundLeft(it)=limits(1,1); XgroundRight(it)=limits(1,2);
    XexcitedLeft(it)=limits(2,1); XexcitedRight(it)=limits(2,2);
    % 几何覆盖率只处理窗口端点落在网格单元内的情况；窗口内部权重恒为 1。
    windowG(:,it)=max(0,min(x+dx/2,limits(1,2))-max(x-dx/2,limits(1,1)))/dx;
    windowE(:,it)=max(0,min(x+dx/2,limits(2,2))-max(x-dx/2,limits(2,1)))/dx;
    N1ground(it)=dx*sum(windowG(:,it).*abs(u1Data(:,it)).^2);
    N2ground(it)=dx*sum(windowG(:,it).*abs(u2Data(:,it)).^2);
    N1excited(it)=dx*sum(windowE(:,it).*abs(u1Data(:,it)).^2);
    N2excited(it)=dx*sum(windowE(:,it).*abs(u2Data(:,it)).^2);
end
% 兼容原绘图变量名；G 为窗口积分，E 在碰撞前积分、碰撞后记为 0。
N1groundSpatial=N1ground; N2groundSpatial=N2ground;
N1excitedSpatial=N1excited; N2excitedSpatial=N2excited;
N1groundPlot=N1ground; N2groundPlot=N2ground;
N1excitedPlot=N1excited; N2excitedPlot=N2excited;
N1groundSpatialPlot=N1ground; N2groundSpatialPlot=N2ground;
N1excitedSpatialPlot=N1excited; N2excitedSpatialPlot=N2excited;
Nground=N1ground+N2ground; Nexcited=N1excited+N2excited;
NgroundSpatial=Nground; NexcitedSpatial=Nexcited;
N1total=dx*sum(abs(u1Data).^2,1); N2total=dx*sum(abs(u2Data).^2,1);
N1outside=N1total-N1ground-N1excited; N2outside=N2total-N2ground-N2excited;
XgroundPlot=Xground; XexcitedPlot=Xexcited;
XgroundPlot(collisionMask)=nan; XexcitedPlot(collisionMask)=nan;
assert(all(isfinite(N1ground(keep))) && all(isfinite(N2ground(keep))) && ...
       all(isfinite(N1excited(keep))) && all(isfinite(N2excited(keep))), ...
       '碰撞区间以外不得有缺失的粒子数。');
assert(max(windowG(:)+windowE(:))<=1+1e-10,'两个窗口重复计数。');
assert(all(N1excited(afterCollision)==0) && all(N2excited(afterCollision)==0), ...
    '碰撞后 E 必须按指定口径记为 0。');

%% 最后保存整个工作区，再统一绘图
CollisionTable=table(tdata(:),tPlot(:),collisionMask(:),afterCollision(:), ...
    Xground(:),Xexcited(:),XgroundLeft(:),XgroundRight(:),XexcitedLeft(:),XexcitedRight(:), ...
    N1ground(:),N2ground(:),N1excited(:),N2excited(:),N1total(:),N2total(:), ...
    'VariableNames',{'Time','TimeOverT','Collision','AfterCollision','Xground','Xexcited', ...
    'GroundLeft','GroundRight','ExcitedLeft','ExcitedRight', ...
    'N1ground','N2ground','N1excited','N2excited','N1total','N2total'});
analysisInfo=struct('definition','Direct integral of original densities inside moving windows', ...
    'collisionRange',collisionRangeUsed,'retainedFrames',nnz(keep), ...
    'postCollisionDefinition','Single tracked packet assigned to G; E set to zero by user convention', ...
    'phase',phaseUsed,'phaseSource',phaseSource, ...
    'windowG',ana.windowG,'windowE',ana.windowE);
if ana.reuseReferenceCache && ~loadedReferenceCache
    save(ana.referenceCacheFile,'isolatedReference','referenceCacheSignature','-v7.3');
end
outputStem=[ana.outputPrefix '_' safe_name_collision(ana.caseName)];
if ana.saveOutput
    save([outputStem '.mat'],'-v7.3');
    writetable(CollisionTable,[outputStem '.csv']);
    fprintf('已保存：%s.mat / %s.csv\n',outputStem,outputStem);
end

%% 绘图参数
fontSize=26; legendSize=22; lineWidth=3; axisWidth=2;
windowAlpha=.18; edgeWidth=3; pathWidth=2;
groundColor=[1,1,1]; excitedColor=[1,0,1];
visibleEdges=[XgroundLeft,XgroundRight,XexcitedLeft,XexcitedRight];
visibleEdges=visibleEdges(isfinite(visibleEdges));
plotLimits=[max(x(1),min(visibleEdges)-.5),min(x(end),max(visibleEdges)+.5)];
% 每段单独填色，不能跨过被删除的碰撞时段。
segmentChange=diff([false,keep(:).',false]);
segmentStart=find(segmentChange==1); segmentEnd=find(segmentChange==-1)-1;
figure('Color','w','Position',[60,60,1450,720],'Name','Actual trajectories and windows', ...
    'Visible','on','WindowStyle','normal');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for j=1:2
    nexttile;
    if j==1, imagesc(x,tPlot,abs(u1Data.').^2); else, imagesc(x,tPlot,abs(u2Data.').^2); end
    set(gca,'YDir','normal'); hold on; colormap(jet(256)); colorbar;
    for iSeg=1:numel(segmentStart)
        ids=segmentStart(iSeg):segmentEnd(iSeg);
        patch([XgroundLeft(ids),fliplr(XgroundRight(ids))], ...
            [tPlot(ids),fliplr(tPlot(ids))],groundColor, ...
            'FaceAlpha',windowAlpha,'EdgeColor','none','HandleVisibility','off');
        if all(isfinite(XexcitedLeft(ids))) && all(isfinite(XexcitedRight(ids)))
            patch([XexcitedLeft(ids),fliplr(XexcitedRight(ids))], ...
                [tPlot(ids),fliplr(tPlot(ids))],excitedColor, ...
                'FaceAlpha',windowAlpha,'EdgeColor','none','HandleVisibility','off');
        end
    end
    plot(XgroundPlot,tPlot,'-','Color',groundColor,'LineWidth',pathWidth,'HandleVisibility','off');
    plot(XexcitedPlot,tPlot,'-','Color',excitedColor,'LineWidth',pathWidth,'HandleVisibility','off');
    hg=plot(XgroundLeft,tPlot,':','Color',groundColor,'LineWidth',edgeWidth);
    plot(XgroundRight,tPlot,':','Color',groundColor,'LineWidth',edgeWidth,'HandleVisibility','off');
    he=plot(XexcitedLeft,tPlot,':','Color',excitedColor,'LineWidth',edgeWidth);
    plot(XexcitedRight,tPlot,':','Color',excitedColor,'LineWidth',edgeWidth,'HandleVisibility','off');
    xlim(plotLimits); ylim([tPlot(1),tPlot(end)]); xlabel('x'); ylabel(timeLabel);
    legend([hg,he],{'G window','E window'},'Location','southoutside', ...
        'Orientation','horizontal','FontSize',legendSize,'FontName','Times New Roman', ...
        'Color',[.1,.1,.1],'TextColor','w','EdgeColor','none');
    text(.03,.95,sprintf('Component %d',j),'Units','normalized','Color','w','FontSize',legendSize);
    set(gca,'FontName','Times New Roman','FontSize',fontSize,'LineWidth',axisWidth);
end
drawnow;

figure('Color','w','Position',[100,100,1200,700],'Name','Particle numbers', ...
    'Visible','on','WindowStyle','normal');
hold on; box on;
plot(tPlot,N1ground,'b-','LineWidth',lineWidth,'DisplayName','N_{1,g}');
plot(tPlot,N2ground,'b:','LineWidth',lineWidth,'DisplayName','N_{2,g}');
plot(tPlot,N1excited,'r-','LineWidth',lineWidth,'DisplayName','N_{1,e}');
plot(tPlot,N2excited,'r:','LineWidth',lineWidth,'DisplayName','N_{2,e}');
xlabel(timeLabel); ylabel('N'); xlim([tPlot(1),tPlot(end)]);
legend('Location','southoutside','FontSize',legendSize,'NumColumns',4);
set(gca,'FontName','Times New Roman','FontSize',fontSize,'LineWidth',axisWidth);
drawnow;

%% 局部计算函数（不含绘图）
function X=track_single_collision(x,t,rho,rhoLocate,ana)
    % 平滑密度定位主峰；局部原始密度质心给出实际轨迹。
    dx=x(2)-x(1); X=nan(1,numel(t));
    if ~isscalar(ana.postSearchHalfWidth) || ~isfinite(ana.postSearchHalfWidth) || ...
            ana.postSearchHalfWidth<=0 || ~isscalar(ana.postCenterHalfWidth) || ...
            ~isfinite(ana.postCenterHalfWidth) || ana.postCenterHalfWidth<=0
        error('碰撞后搜索半宽与质心半宽必须为正有限数。');
    end
    if ~isempty(ana.postStartX) && (~isscalar(ana.postStartX) || ~isfinite(ana.postStartX))
        error('ana.postStartX 应为空或一个有限位置。');
    end
    for it=1:numel(t)
        if it==1
            prediction=ana.postStartX;
        else
            prediction=X(it-1);
            if it>2
                step=(X(it-1)-X(it-2))*(t(it)-t(it-1))/(t(it-1)-t(it-2));
                prediction=prediction+max(-ana.postSearchHalfWidth/2, ...
                    min(ana.postSearchHalfWidth/2,step));
            end
            prediction=max(x(1),min(x(end),prediction));
        end
        if isempty(prediction), ix=(1:numel(x)).';
        else, ix=find(abs(x-prediction)<=ana.postSearchHalfWidth); end
        if isempty(ix), error('碰撞后第 %d 帧的轨迹搜索超出网格。',it); end
        [~,ip]=max(rhoLocate(ix,it)); p=ix(ip); peak=x(p);
        if p>1 && p<numel(x)
            v=rhoLocate(p-1:p+1,it); curvature=v(1)-2*v(2)+v(3);
            if curvature<0
                peak=peak+dx*max(-.5,min(.5,.5*(v(1)-v(3))/curvature));
            end
        end
        h=ana.postCenterHalfWidth;
        w=max(0,min(x+dx/2,peak+h)-max(x-dx/2,peak-h))/dx;
        mass=sum(w.*rho(:,it));
        if mass<=0 || ~isfinite(mass), error('碰撞后第 %d 帧没有可追踪的波包。',it); end
        X(it)=sum(x.*w.*rho(:,it))/mass;
    end
end

function signature=reference_signature_collision(x,k,t,phiG,phiE,lat,S,ana)
    if isfield(S,'job'), job=S.job; else, job=struct(); end
    if isempty(ana.referenceDt)
        if isfield(S,'pump') && isfield(S.pump,'dt'), dt=S.pump.dt; else, dt=1e-3; end
    else, dt=ana.referenceDt; end
    if isempty(ana.referenceUseFourthOrder)
        if isfield(S,'pump') && isfield(S.pump,'useFourthOrder')
            fourth=S.pump.useFourthOrder;
        else, fourth=true; end
    else, fourth=ana.referenceUseFourthOrder; end
    phase=ana.actualGroundExtraPhase;
    if isempty(phase)
        if isfield(S,'vphi'), phase=S.vphi; else, phase=0; end
    end
    signature=struct('version',2,'x',x,'k',k,'time',t, ...
        'ground',phiG,'excited',phiE,'lat',lat,'job',job, ...
        'dt',dt,'fourth',fourth,'groundExtraPhase',phase);
end

function [Xg,Xe,fitError]=track_pair_collision( ...
        x,k,rho,rhoG,rhoE,Xgr,Xer,ana)
    dx=x(2)-x(1); Nt=size(rho,2); radius=ana.trackReacquireHalfWidth;
    C=nan(2,Nt); fitError=nan(1,Nt); references=[Xgr;Xer];
    options=optimset('Display','off','TolX',dx/100,'TolFun',1e-8, ...
        'MaxIter',80,'MaxFunEvals',180);
    for it=1:Nt
        anchor=references(:,it);
        prediction=anchor;
        if it>1 && all(isfinite(C(:,it-1)))
            prediction=C(:,it-1)+references(:,it)-references(:,it-1);
        end
        cg=x(abs(x-prediction(1))<=radius | abs(x-anchor(1))<=radius).';
        ce=x(abs(x-prediction(2))<=radius | abs(x-anchor(2))<=radius).';
        if isempty(cg) || isempty(ce), error('第 %d 帧无法定位波包。',it); end
        fg=fft(rhoG(:,it)); fe=fft(rhoE(:,it)); y=rho(:,it);
        Tg=tracking_templates_collision(x,k,fg,cg,anchor(1),ana);
        Te=tracking_templates_collision(x,k,fe,ce,anchor(2),ana);
        [err,~,~]=positive_pair_fit_collision(Tg,Te,y,dx);
        denominator=dx*sum(y.^2)+eps;
        cost=err/denominator+ana.trackMotionPenalty* ...
            ((cg(:)-prediction(1)).^2+(ce-prediction(2)).^2)/radius^2;
        [best,jbest]=min(cost(:));
        if ~isfinite(best), error('第 %d 帧轨迹匹配失败。',it); end
        [ig,ie]=ind2sub(size(cost),jbest); seed=[cg(ig);ce(ie)];
        % 网格内细化，不需要 Optimization Toolbox。
        objective=@(z)tracking_objective_collision(z,seed,x,k,fg,fe,anchor,y,dx,ana);
        candidate=fminsearch(objective,seed,options);
        if objective(candidate)<=objective(seed), c=candidate; else, c=seed; end
        Tg=tracking_templates_collision(x,k,fg,c(1),anchor(1),ana);
        Te=tracking_templates_collision(x,k,fe,c(2),anchor(2),ana);
        err=positive_pair_fit_collision(Tg,Te,y,dx);
        C(:,it)=c; fitError(it)=err/denominator;
        if mod(it,max(1,round(Nt/10)))==0 || it==Nt
            fprintf('实际波包跟踪：%5.1f%%\n',100*it/Nt);
        end
    end
    Xg=C(1,:); Xe=C(2,:);
end

function T=tracking_templates_collision(x,k,f,c,anchor,ana)
    c=c(:).';
    T=max(0,real(ifft(f.*exp(-1i*k*(c-anchor)),[],1)));
    r=abs(x-c); h=ana.trackWindowHalfWidth; taper=ana.trackTaperWidth;
    if taper<=0
        W=double(r<=h);
    else
        inner=max(0,h-taper); W=zeros(size(r)); W(r<=inner)=1;
        mask=r>inner & r<h;
        W(mask)=0.5*(1+cos(pi*(r(mask)-inner)/(h-inner)));
    end
    T=T.*W;
end

function [err,ag,ae]=positive_pair_fit_collision(G,E,y,dx)
    % 两列非负最小二乘的解析解，包括其中一个波包振幅为零的边界解。
    gg=dx*sum(G.^2,1).'; ee=dx*sum(E.^2,1);
    ge=dx*(G.'*E); gy=dx*(G.'*y); ey=dx*(y.'*E);
    determinant=gg.*ee-ge.^2;
    good=determinant>1e-10*(gg.*ee) & determinant>eps^2;
    safe=determinant; safe(~good)=1;
    ag=(gy.*ee-ge.*ey)./safe; ae=(gg.*ey-ge.*gy)./safe;
    yy=dx*sum(y.^2);
    err=yy-ag.*gy-ae.*ey;
    err(~good | ag<0 | ae<0)=inf;
    onlyG=max(gy./max(gg,eps),0); onlyE=max(ey./max(ee,eps),0);
    errG=repmat(yy-onlyG.*gy,1,size(E,2));
    mask=errG<err; bg=repmat(onlyG,1,size(E,2));
    err(mask)=errG(mask); ag(mask)=bg(mask); ae(mask)=0;
    errE=repmat(yy-onlyE.*ey,size(G,2),1);
    mask=errE<err; be=repmat(onlyE,size(G,2),1);
    err(mask)=errE(mask); ag(mask)=0; ae(mask)=be(mask);
    err=max(0,err);       % 只去掉舍入导致的负平方残差
end

function value=tracking_objective_collision(c,seed,x,k,fg,fe,anchor,y,dx,ana)
    if any(abs(c-seed)>dx) || any(~isfinite(c)), value=inf; return; end
    G=tracking_templates_collision(x,k,fg,c(1),anchor(1),ana);
    E=tracking_templates_collision(x,k,fe,c(2),anchor(2),ana);
    value=positive_pair_fit_collision(G,E,y,dx)/(dx*sum(y.^2)+eps);
end

function xc = density_center_collision(x,u,dx)
    rho = abs(u).^2;
    xc = real(dx*sum(x.*rho)/(dx*sum(rho)+eps));
end

function xc = density_center_from_rho_collision(x,rho,dx)
    xc=real(dx*sum(x.*rho)/(dx*sum(rho)+eps));
end

function [uG1,uG2,uE1,uE2]=initial_vectors_collision(phiGround,phiExcited,S,ana)
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
end

function ref = evolve_isolated_references_collision( ...
        x,k,tdata,phiGround,phiExcited,lat,S,ana)
    fprintf('\n========== Isolated dynamic-reference evolution ==========\n');
    [uG1,uG2,uE1,uE2]=initial_vectors_collision(phiGround,phiExcited,S,ana);

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

function V = external_potential_collision(x,t,lat)
    V = -lat.p1*cos(pi*x/lat.d1).^2 ...
        -lat.p2*cos(pi*x/lat.d2-lat.nu*t).^2;
end

function s = safe_name_collision(s)
    s=regexprep(char(s),'[^a-zA-Z0-9_.-]','_');
end
