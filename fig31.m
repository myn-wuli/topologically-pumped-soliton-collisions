function out = fig31(resultFile,userCfg)
% 第二组分：跟踪实际主波包 -> 取窗 -> 瞬时裸晶格能带投影。
% out = fig31('g12=-5,G=10,E=10.mat');
% 已手动加载原始数据时：out = fig31;
% 可选设置：out = fig31(file,struct('window',[1,1]));
% Pband/PbandCore 按窗内粒子数归一化；PbandFull 按全空间粒子数归一化。
% 取窗作用于波函数：uCore=sqrt(windowWeight).*u2，而非直接乘已算出的占比。
% 不去除平均动量，不强制任一能带占比为 0 或 1。

%% 设置
if nargin<1, resultFile=''; end
if nargin<2, userCfg=struct(); end
cfg.window=[1.0,1.0];       % 主轨迹左、右两侧的积分宽度
cfg.searchHalfWidth=1.5;    % 相邻帧主峰搜索半宽
cfg.centerHalfWidth=.75;   % 主峰附近局部质心的半宽
cfg.reacquireRatio=.35;    % 局部峰过弱时重新定位第二组分最强峰
cfg.startX=[];             % 空：第二组分初帧最强峰；也可填初始中心
cfg.nBandsTrack=5;
cfg.savePrefix='Pc';       % 输出名自动附 g12 和初始峰值密度
cfg.saveOutput=true;
cfg.makePlots=true;
cfg.saveFigures=true;
cfg.plotXLim=[];           % 空：随完整主轨迹自动设置
cfg.fontSize=50; cfg.legendSize=30; cfg.lineWidth=6; cfg.axisWidth=2;
if ~isstruct(userCfg) || ~isscalar(userCfg), error('userCfg 必须为标量结构体。'); end
fields=fieldnames(userCfg);
for j=1:numel(fields)
    if ~isfield(cfg,fields{j}), error('未知设置：%s。',fields{j}); end
    cfg.(fields{j})=userCfg.(fields{j});
end
if numel(cfg.window)~=2 || any(~isfinite(cfg.window)) || any(cfg.window<=0)
    error('cfg.window 必须为两个正有限数：[左宽,右宽]。');
end
cfg.window=cfg.window(:).';
if ~isscalar(cfg.searchHalfWidth) || ~isfinite(cfg.searchHalfWidth) || cfg.searchHalfWidth<=0 || ...
        ~isscalar(cfg.centerHalfWidth) || ~isfinite(cfg.centerHalfWidth) || cfg.centerHalfWidth<=0 || ...
        ~isscalar(cfg.reacquireRatio) || ~isfinite(cfg.reacquireRatio) || ...
        cfg.reacquireRatio<=0 || cfg.reacquireRatio>1
    error('跟踪半宽须为正有限数，reacquireRatio 须在 (0,1]。');
end
if ~isempty(cfg.startX) && (~isscalar(cfg.startX) || ~isfinite(cfg.startX))
    error('cfg.startX 应为空或一个有限位置。');
end

%% 读取原始波函数
if isempty(resultFile)
    names={'x','tdata','udata1','udata2','lat','psi','psiGround'}; D=struct();
    for j=1:numel(names)
        nm=names{j};
        if evalin('base',['exist(''' nm ''',''var'')==1'])
            D.(nm)=evalin('base',nm);
        end
    end
else
    D=load(resultFile);
end
required={'x','tdata','udata1','udata2','lat'};
for j=1:numel(required)
    if ~isfield(D,required{j}), error('缺少原始变量 %s。',required{j}); end
end
x=D.x(:); t=D.tdata(:).'; lat=D.lat; N=numel(x); Nt=numel(t);
if N<4 || Nt<2 || any(~isfinite(x)) || any(~isfinite(t)) || any(diff(t)<=0)
    error('x/tdata 网格无效。');
end
dx=mean(diff(x)); L=N*dx;
if dx<=0 || max(abs(diff(x)-dx))>1e-8*dx, error('x 必须为递增等距网格。'); end
u1all=D.udata1; u2all=D.udata2;
if size(u1all,1)~=N && size(u1all,2)==N, u1all=u1all.'; end
if size(u2all,1)~=N && size(u2all,2)==N, u2all=u2all.'; end
if ~isequal(size(u1all),[N,Nt]) || ~isequal(size(u2all),[N,Nt]) || ...
        any(~isfinite(u1all(:))) || any(~isfinite(u2all(:)))
    error('udata1/udata2 应为 length(x) x length(tdata) 的有限复波函数。');
end
if ~isfield(lat,'a'), lat.a=lat.d2; end
requiredLat={'a','p1','p2','d1','d2','nu','g12'};
for j=1:numel(requiredLat)
    nm=requiredLat{j};
    if ~isfield(lat,nm) || ~isscalar(lat.(nm)) || ~isfinite(lat.(nm))
        error('lat.%s 必须为有限标量。',nm);
    end
end
if min([lat.a,lat.d1,lat.d2])<=0, error('晶格周期必须为正。'); end
M=round(L/lat.a); Npc=round(lat.a/dx);
if M<1 || Npc<3 || M*Npc~=N || abs(M*lat.a-L)>1e-8*L || ...
        any(abs(lat.a./[lat.d1,lat.d2]-round(lat.a./[lat.d1,lat.d2]))>1e-10)
    error('计算区域须含整数个晶胞，a/dx、a/d1、a/d2 须为整数。');
end
if ~isscalar(cfg.nBandsTrack) || cfg.nBandsTrack<2 || ...
        cfg.nBandsTrack~=round(cfg.nBandsTrack) || ~isfinite(cfg.nBandsTrack)
    error('nBandsTrack 必须为不小于 2 的整数。');
end
nb=min(cfg.nBandsTrack,Npc);
rho1=abs(u1all).^2; rho2=abs(u2all).^2;
N1=dx*sum(rho1,1); N2=dx*sum(rho2,1);
if any(N2<=0), error('第二组分存在零粒子数帧，无法定义占比。'); end
if lat.nu~=0
    T=pi/abs(lat.nu); tPlot=t/T; timeLabel='t/T';
else
    T=nan; tPlot=t; timeLabel='t';
end

%% 逐帧跟踪第二组分主波包：局部质心不受远处辐射拖动
[xcCore,reacquired]=track_core_band(x,t,rho2,cfg);
windowLeft=max(x(1)-dx/2,xcCore-cfg.window(1));
windowRight=min(x(end)+dx/2,xcCore+cfg.window(2));
windowWeight=max(0,min(x+dx/2,windowRight)-max(x-dx/2,windowLeft))/dx;
windowWeight=min(1,windowWeight);
Ncore=dx*sum(windowWeight.*rho2,1);
if any(Ncore<=0), error('主波包窗口为空，请检查跟踪参数。'); end
coreFraction=Ncore./N2;
Noutside=N2-Ncore;

%% 在完整 FFT 网格上按准动量分块，使用每帧瞬时裸势 H0(t)
% 每个块的全部正交本征态构成该准动量的完整基底。
% 全空间波函数与取窗波函数共享同一基底，不做 k/相位网格插值。
n=[0:ceil(N/2)-1,-floor(N/2):-1].'; k=2*pi*n/L;
blockIdx=cell(M,1); potentialIdx=cell(M,1); kinetic=cell(M,1);
for iq=1:M
    idx=find(mod(n,M)==iq-1);
    blockIdx{iq}=idx;
    potentialIdx{iq}=mod(n(idx)-n(idx).',N)+1;
    kinetic{iq}=diag(.5*k(idx).^2);
end
PallFull=zeros(Npc,Nt); PallCore=PallFull;
for it=1:Nt
    V=-lat.p1*cos(pi*x/lat.d1).^2-lat.p2*cos(pi*x/lat.d2-lat.nu*t(it)).^2;
    Vhat=fft(V)/N;
    Uhat=fft([u2all(:,it),sqrt(windowWeight(:,it)).*u2all(:,it)],[],1)/sqrt(N);
    bandMass=zeros(Npc,2);
    for iq=1:M
        idx=blockIdx{iq}; H=kinetic{iq}+Vhat(potentialIdx{iq}); H=(H+H')/2;
        [B,e]=eig(H,'vector'); [~,order]=sort(real(e),'ascend'); B=B(:,order);
        coefficients=B'*Uhat(idx,:);
        bandMass=bandMass+dx*abs(coefficients).^2;
    end
    PallFull(:,it)=bandMass(:,1)/N2(it);
    PallCore(:,it)=bandMass(:,2)/Ncore(it);
    if mod(it,max(1,round(Nt/10)))==0 || it==Nt
        fprintf('瞬时能带投影：%5.1f%%\n',100*it/Nt);
    end
end
closureFull=sum(PallFull,1); closureCore=sum(PallCore,1);
if max(abs([closureFull,closureCore]-1))>1e-9
    error('投影完备性检查失败：所有能带占比之和应为 1。');
end
PbandFull=PallFull(1:nb,:); PbandCore=PallCore(1:nb,:);
PotherFull=sum(PallFull(3:end,:),1); PotherCore=sum(PallCore(3:end,:),1);
E12=lat.g12*dx*sum(rho1.*rho2,1);

%% 保存整个计算工作区；所有绘图在保存之后
out=struct(); out.schema='tracked-core-bands-v1'; out.cfg=cfg; out.lat=lat;
out.x=x; out.t=t; out.T=T; out.tPlot=tPlot; out.timeLabel=timeLabel;
out.xcCore=xcCore; out.windowLeft=windowLeft; out.windowRight=windowRight;
out.windowWeight=windowWeight; out.reacquired=reacquired;
out.N1=N1; out.N2=N2; out.Ncore=Ncore; out.Noutside=Noutside;
out.coreFraction=coreFraction; out.rho2=rho2; out.E12=E12;
out.Pband=PbandCore; out.PbandCore=PbandCore; out.PbandFull=PbandFull;
out.PallCore=PallCore; out.PallFull=PallFull;
out.PotherCore=PotherCore; out.PotherFull=PotherFull;
out.Phigher=sum(PallCore(nb+1:end,:),1);
out.projectionClosureCore=closureCore; out.projectionClosureFull=closureFull;
out.definitionCore='||Pi_n(t)[sqrt(W(t))*psi2(t)]||^2 / integral W*|psi2|^2 dx';
out.definitionFull='||Pi_n(t)psi2(t)||^2 / integral |psi2|^2 dx';
out.basis='Instantaneous bare lattice H0=-0.5*d_xx+V(x,t), full periodic simulation grid';
out.note='Windowed spectral weights depend on the chosen window; nonlinear-mode populations are not imposed.';
Gpeak=max(rho1(:,1)); Epeak=max(rho2(:,1));
if isfield(D,'psiGround'), Gpeak=max(abs(D.psiGround(:)).^2); end
if isfield(D,'psi'), Epeak=max(abs(D.psi(:)).^2); end
outputStem=sprintf('%s_g%g_G%g_E%g',cfg.savePrefix,lat.g12,Gpeak,Epeak);
outputStem=regexprep(outputStem,'[^a-zA-Z0-9_.-]','_'); out.outputStem=outputStem;
Results=table(t(:),tPlot(:),xcCore(:),windowLeft(:),windowRight(:),N2(:),Ncore(:), ...
    coreFraction(:),PbandCore(1,:).',PbandCore(2,:).',PotherCore(:), ...
    PbandFull(1,:).',PbandFull(2,:).',PotherFull(:),E12(:), ...
    'VariableNames',{'Time','TimeAxis','Xcore','Left','Right','N2','Ncore','CoreFraction', ...
    'P1core','P2core','PotherCore','P1full','P2full','PotherFull','E12'});
fprintf('\n末帧 t=%g：窗口占比=%g，核心 P1=%g、P2=%g；全局 P1=%g、P2=%g。\n', ...
    tPlot(end),coreFraction(end),PbandCore(1,end),PbandCore(2,end),PbandFull(1,end),PbandFull(2,end));
if cfg.saveOutput
    save([outputStem '.mat'],'-v7.3');
    writetable(Results,[outputStem '.csv']);
    fprintf('已保存：%s.mat / %s.csv\n',outputStem,outputStem);
end
if cfg.makePlots
    plot_collision_conversion_results_a; % 单独绘图脚本，直接使用当前 out
end
end

function [xc,reacquired]=track_core_band(x,t,rho,cfg)
% 跟踪第二组分自身，避免初期误追第一组分的基态。
dx=x(2)-x(1); Nt=numel(t); xc=nan(1,Nt); reacquired=false(1,Nt);
for it=1:Nt
    [globalPeak,ig]=max(rho(:,it));
    if it==1
        prediction=cfg.startX;
    else
        prediction=xc(it-1);
        if it>2
            step=(xc(it-1)-xc(it-2))*(t(it)-t(it-1))/(t(it-1)-t(it-2));
            prediction=prediction+max(-cfg.searchHalfWidth/2,min(cfg.searchHalfWidth/2,step));
        end
        prediction=max(x(1),min(x(end),prediction));
    end
    if isempty(prediction), ids=(1:numel(x)).';
    else, ids=find(abs(x-prediction)<=cfg.searchHalfWidth); end
    if isempty(ids), error('第 %d 帧主波包搜索范围超出网格。',it); end
    [localPeak,j]=max(rho(ids,it)); p=ids(j);
    if it>1 && localPeak<cfg.reacquireRatio*globalPeak
        p=ig; reacquired(it)=true;
    end
    peak=x(p); h=cfg.centerHalfWidth;
    w=max(0,min(x+dx/2,peak+h)-max(x-dx/2,peak-h))/dx;
    mass=sum(w.*rho(:,it));
    if mass<=0, error('第 %d 帧主峰附近粒子数为零。',it); end
    xc(it)=sum(x.*w.*rho(:,it))/mass;
end
end
