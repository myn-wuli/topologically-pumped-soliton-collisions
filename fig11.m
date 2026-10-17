clc;
clear; 
close all;
% --------------------------------------------------------
% 伪谱法求解线性实空间能带, 计算Bloch能带, 计算超胞布洛赫态, 
% 计算陈数
% --------------------------------------------------------
L = 32;               % 区间长度 [-L/2, L/2]
N = 640;              % 网格数
Ncell = N/L;     % 单原胞离散点数
dx = L / N;       % 离散步长
I = eye(N);                   % 单位矩阵
T_full = circshift(I, -Ncell, 2);  % 沿列循环平移 Ncell， 向右平移了一个单位

xarr = linspace(-L/2, L/2 - dx, N).';   % 空间坐标列向量
k = (-N/2:N/2 - 1).' * (2*pi/L);   % 傅里叶频率（列向量）
k = fftshift(k);                  % 将频率中心化
K2 = k.^2;

u = 1; Nt = 200; T = pi/u; tarr = linspace(0, T, Nt); % u = 1; 没有涉及含时演化，这里全是解析计算所以u*t在这里相当于ky,空间周期是kx

%p1 = 25; p2 = 25; d1 = 1/2; d2 = 2/3; a = 2; % 晶格常数
p1 = 25; p2 = 25; d1 = 1/2; d2 = 1; a = 1; % 晶格常数
Allcell = L/a; % 总的原胞数
Vp = @(t) -(p1*cos(pi*xarr/d1).^2 + p2*cos(pi*xarr/d2 - u*t).^2); % 该周期是1

Parameters = struct('L', L, 'N', N, 'k', k, 'K2', K2, 'Vp', Vp, ...
    'tarr', tarr, 'Nt', Nt,'a', a, 'Allcell', Allcell, 'T_full', T_full);


% % 平移算符检查
% Gpsi = sech(xarr);
% Gpsi2 = T_full*Gpsi;
% hold on
% 
% plot(xarr, Gpsi)
% plot(xarr, Gpsi2,'--')


%%

[UPsi, ERealMatrix] = Real_Spectrum(Parameters); %  UPsi = zeros(Nt, N, N);
% 计算按照能量排序的布洛赫态, UPsi 是直接对角化的实空间所有特征态 Band 第几个能带
Band = 1;  % 选择能带
[SBlochV, EBloch] = Bloch_Psi(UPsi, Band, Parameters) ; % SBlochV = zeros(Nt, N, Allcell);


%%
figure  % 直接对角化特征态态在布洛赫态上的投影
hold on 
nn = 2;
RealV = squeeze(UPsi(1,:, (Band-1)*Allcell + nn ));

plot(xarr, real(RealV))
plot(xarr, real(SBlochV(1, :, nn)),'.' )
legend('实空间','Bloch')
hold off
% SBlochV(1, :, nn)
%%
figure % 实空间能谱
hold on 
plot(tarr/pi, ERealMatrix.','k-')
plot(tarr/pi, ERealMatrix(:,1:Allcell).','r-')
plot(tarr/pi, ERealMatrix(:,Allcell+1:2*Allcell).','g-')
plot(tarr/pi, ERealMatrix(:,2*Allcell+1:3*Allcell).','b-')
ylim([-30,-5])
%ylim([-40,20])
title('实空间能谱');
box on
ylabel('$\mu$','FontSize',40,'FontName','Times New Roman',Interpreter='latex');
xlabel('$t/T$','FontSize',40,'FontName','Times New Roman',Interpreter='latex');
zlabel('Density','FontSize',40,'FontName','Times New Roman',Interpreter='latex');
set(gca,'LineWidth',2);
set(gca,'FontSize',40,'FontName','Times New Roman');
h = gcf; % 获取当前图形句柄
print(h, '-djpeg', '-r600', 'fig2.jpg'); % 保
% plot(xarr, abs(psi).^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        下面是对角化布洛赫哈密顿量部分，这里可以计算扩展超胞布洛赫和陈数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p1 = 25; p2 = 25; d1 = 1; d2 = 1/2; a = 1; % 晶格常数
N = 64; dx = a / N;
xarr = linspace(-a/2, a/2 - dx, N).';
Vp = @(t) -(p1*cos(pi*xarr/d1).^2 + p2*cos(pi*xarr/d2 - u*t).^2); % 该周期是1
G = (-N/2: N/2 - 1).'*(2*pi/a);   % 倒格子基矢

BZkx = linspace(0, 2*pi/a, 200); % 晶格准动量维度
BZky = linspace(0, pi, 200); % 时间准动量维度
% BZky = 0; % 时间准动量维度

Parameters.Vp = Vp; Parameters.G = G; Parameters.BZkx = BZkx; Parameters.BZky = BZky;
Parameters.xarr = xarr; Parameters.N = N;
[EBlochMatrix, Bloch_U] = Bloch_Spectrum(Parameters); % 默认选择了前三条低能带
%%
figure
Bloch_U1 = squeeze(Bloch_U(1, 1, :, 1));
hold on 
plot(xarr, Vp(0))
plot(xarr, Vp(pi/4),'o')


%%  构建超胞态的的布洛赫态  假设总原胞数为16个；
Superpsi = zeros(16*N, 1);
superxarr = linspace(-8*a, 8*a - dx, 16*N).';
Vp = @(t) -(p1*cos(pi*superxarr/d1).^2 + p2*cos(pi*superxarr/d2 - u*t).^2); % 该周期是1
N_L = 8; N_R = 7;
% 循环叠加
nk = 1;
Cell = squeeze(Bloch_U(nk, 1, :, 1));
kk = BZkx(nk);
for n = -N_L : N_R
    idx_start = (n + N_L)*N + 1;  % 索引起点
    idx_end = idx_start + N - 1;  % 索引终点
    Superpsi(idx_start: idx_end) = Cell * exp(1i * kk * n);
end


%%
figure % 前几个布洛赫能级 % 默认选择了前三条低能带
hold on 
EBlochsqu1 = squeeze(EBlochMatrix(:, 100, :)); % 定下某个ky
EBlochsqu = squeeze(EBlochMatrix(:, 1, :)); % 定下某个ky
plot(BZkx, EBlochsqu1.','k-')
plot(BZkx, EBlochsqu.','k-')
%
figure % 前几个布洛赫能级 % 默认选择了前三条低能带
hold on
surf(BZkx, BZky,EBlochMatrix(:, :, 1).')
surf(BZkx, BZky,EBlochMatrix(:, :, 2).')
surf(BZkx, BZky,EBlochMatrix(:, :, 3).')
surf(BZkx, BZky,EBlochMatrix(:, :, 4).')
zlabel('E')
xlabel('Kx')
ylabel('Ky')
title('K空间能谱')
shading interp  %阴影插值函数
colormap('jet')
%%   计算陈数

Bloch_U_Band1 = squeeze(Bloch_U(:, :, :, 1)) ; % 选择一个带
[Chen_number1, Berry_curvature1] = Wilson_loop_2D(BZkx, BZky, Bloch_U_Band1);
Bloch_U_Band2 = squeeze(Bloch_U(:, :, :, 2)) ; % 选择一个带
[Chen_number2, Berry_curvature2] = Wilson_loop_2D(BZkx, BZky, Bloch_U_Band2);

hold on 
surf(BZkx(1:end-1), BZky(1:end-1), Berry_curvature1)
% surf(BZkx(1:end-1), BZky(1:end-1), Berry_curvature2)
shading interp  %阴影插值函数
colormap('jet')
xlabel('Kx')
ylabel('Ky')
title('磁场分布')


%%
% save('Linear Spectrum.mat');

function [Chen_number, Berry_curvature] = Wilson_loop_2D(BZkx, BZky, Bloch_U_Band) % wilson loop 计算陈数，无需规范依赖 
    lenkx = length(BZkx);
    lenky = length(BZky);
    Berry_curvature = zeros(lenkx-1,lenky-1);
    for kxj = 1 : lenkx - 1
        disp(kxj/lenkx)
        for kyj = 1 : lenky -1 
            F1 = dot(Bloch_U_Band(kxj, kyj, :) ,Bloch_U_Band(kxj + 1, kyj, :));
            F2 = dot(Bloch_U_Band(kxj+1, kyj, :) ,Bloch_U_Band(kxj+1, kyj+1, :));
            F3 = dot(Bloch_U_Band(kxj+1,kyj+1, :) ,Bloch_U_Band(kxj, kyj+1, :));
            F4 = dot(Bloch_U_Band(kxj, kyj+1, :) ,Bloch_U_Band(kxj, kyj, :));
            loop = angle(F1*F2*F3*F4);
            Berry_curvature(kxj, kyj) = loop;
        end  
    end
    Chen_number = sum(Berry_curvature,"all")/2/pi;
end


% 计算按照能量排序的布洛赫态, UPsi 是直接对角化的实空间所有特征态 Band 第几个能带
function [SBlochV, EBloch] = Bloch_Psi(UPsi, Band, Parameters) 
Nt = Parameters.Nt; Allcell = Parameters.Allcell; N = Parameters.N;
T_full = Parameters.T_full; tarr = Parameters.tarr;

SBlochV = zeros(Nt, N, Allcell);
EBloch = zeros(Nt, Allcell);
    for ti = 1 : Nt
        disp(ti/Nt)
        V = squeeze(UPsi(ti, :, (Band - 1)*Allcell + 1 : Band * Allcell)); % 子空间基矢 维度 N*L
        SubT = V' * T_full * V; % 构建子空间平移算子
        [ym, ~] = eig(SubT);
        BlochV = V * ym;
        Hami = R_Hami(tarr(1), Parameters);
        [E, idx] = sort(real(diag(BlochV'*Hami*BlochV)));
        EBloch(ti, :) = E;
        SBlochV(ti, :, :) = BlochV(:, idx);
    end
end

function [UPsi, EMatrix] = Real_Spectrum(Parameters)  % 计算实空间能谱
Nt = Parameters.Nt; N = Parameters.N; 
tarr = Parameters.tarr; 
EMatrix = zeros(Nt, N); UPsi = zeros(Nt, N, N);
    for ti = 1 : Nt
        disp(ti/Nt)
        t = tarr(ti);
        Hami = R_Hami(t, Parameters);
        [U, E] = eig(Hami);  % U: 本征函数, E: 对角本征值
        [Evalues, idx] = sort(real(diag(E)));   % 提取本征值 ε
        sortU = U(:,idx);
        UPsi(ti, :, :) = sortU;
        EMatrix(ti, :) = Evalues;
    end
end

function Hami = R_Hami(t, Parameters)  % 实空间哈密顿量
    N = Parameters.N; K2 = Parameters.K2; Vp = Parameters.Vp; 
    
    T_hat = 0.5 * diag(K2);           % 傅里叶空间动能算符
    F = fft(eye(N))/sqrt(N);          % 傅里叶变换矩阵
    F_inv = ifft(eye(N))*sqrt(N);     % 逆变换
    T_x = real(F_inv * T_hat * F);    % 动能在位置空间的表示（Hermitian）
    V = Vp(t); 
    V_op = diag(V);                   % 势能算符矩阵
    H = T_x + V_op;                   % 整体哈密顿矩阵（Hermitian）
    Hami = (H + H')/2; 
end

function [EMatrix, Bloch_U] = Bloch_Spectrum(Parameters)  % 计算Bloch空间能谱
N = Parameters.N;  BZkx = Parameters.BZkx; BZky = Parameters.BZky;
Vp = Parameters.Vp; G = Parameters.G;
xarr = Parameters.xarr;
% k 实际上就是G倒格子基矢
Nkx = length(BZkx); Nky = length(BZky);
Band = 4;   %默认能带选取
EMatrix = zeros(Nkx, Nky, Band);
Bloch_U = zeros(Nkx, Nky, N, Band); % 取前三个带的布洛赫态

    for kyi = 1 : Nky
        disp(kyi/Nky)
        ky = BZky(kyi);
        Vx = Vp(ky); % 特定时刻的势场分布
        V_G = fftshift(fft(Vx)/N);   % 归一化 + 中心化
        V_fftp = [V_G(N/2+1 : end); zeros(N/2, 1)];
        V_fftn = [V_G(N/2+1:-1:2); zeros(N/2, 1)];
        V_m = toeplitz(V_fftp, V_fftn); 
        
        for kxi = 1 : Nkx
            kx = BZkx(kxi);
            T_m =diag( 0.5 * (kx + G).^2); 
            H = T_m + V_m;                   % 整体哈密顿矩阵（Hermitian）
            H = (H + H')/2; 
            [coeffs, E] = eig(H);  % U: 本征函数, E: 对角本征值
            [Evalues, idx] = sort(real(diag(E)));   % 提取本征值 ε
            Scoeffs = coeffs(:, idx);
            EMatrix(kxi, kyi, :) = Evalues(1: Band);
            % Scoeffs: N x 3, 列对应三条带
            Scoeffs_shifted = fftshift(Scoeffs(:, 1 : Band), 1);   % 对第一维（行）做 shift
            u_x = ifft(Scoeffs_shifted, [], 1) * N;  % N x 3, 每列是一条带的 u(x)
            psi_x_mat = exp(1i*kx*xarr).* u_x;  %  布洛赫态
            Bloch_U(kxi, kyi, :, :) = psi_x_mat;% 最低带的布洛赫态

        end
    end
end


