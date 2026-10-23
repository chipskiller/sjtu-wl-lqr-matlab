% 计算K矩阵
clear;

%% 符号定义

% 变量
syms theta_w_l theta_w_r; % 驱动轮转角
syms ddot_theta_w_l ddot_theta_w_r; % 驱动轮转角加速度
syms theta_l_l theta_l_r; % 腿倾转角
syms dot_theta_l_l dot_theta_l_r; % 腿倾转角速度
syms ddot_theta_l_l ddot_theta_l_r; % 腿倾转角加速度
syms theta_b; % 机体倾斜角
syms dot_theta_b; % 机体倾斜角速度
syms ddot_theta_b; % 机体倾斜角加速度
syms phi; % 偏航角
syms dot_phi; % 偏航角速度
syms ddot_phi; % 偏航角加速度
syms s; % 自然坐标系下机器人水平方向移动距离
syms dot_s; % 自然坐标系下机器人水平方向移动速度
syms ddot_s; % 自然坐标系下机器人水平方向移动加速度
syms s_b; % 自然坐标系下机体质心水平方向移动距离
syms h_b; % 自然坐标系下机体质心垂直方向移动距离
syms s_l_l s_l_r; % 自然坐标系下腿部质心水平方向移动距离
syms h_l_l h_l_r; % 自然坐标系下腿部质心垂直方向移动距离
syms T_w_l T_w_r; % 驱动轮转矩（腿-轮）
syms T_b_l T_b_r; % 腿部转矩（机体-腿）
syms f_l f_r; % 地面对驱动轮摩擦力
syms F_ws_l F_ws_r; % 驱动轮对腿水平方向作用力
syms F_wh_l F_wh_r; % 驱动轮对腿竖直方向作用力
syms F_bs_l F_bs_r; % 腿对机体水平方向作用力
syms F_bh_l F_bh_r; % 腿对机体竖直方向作用力

% 参数 ------------------
syms g; % 重力加速度

syms R_w; % 驱动轮半径
syms R_l; % 驱动轮轮距/2

syms l_l l_r; % 腿部长度
syms l_w_l l_w_r; % 驱动轮到腿部质心距离*1 (1.3)
syms l_b_l l_b_r; % 机体质心到腿部质心距离*1 (1.2)
syms l_c; % 机体质心到腿部关节距离

syms m_w; % 驱动轮质量
syms m_l; % 腿部质量
syms m_b; % 机体质量

syms I_w; % 驱动轮转动惯量
syms I_l_l I_l_r; % 腿部转动惯量(1.1)
syms I_b; % 机体转动惯量
syms I_z; % 机体yaw转动惯量

%% 1、机器人参数
% ─── 车奉朝轮腿模型物理参数 ────────────────────────────
%  来源: robot.xml  +  env.py + 估算
% ─────────────────────────────────────────────────────────
g_val   = 9.81;      % 重力加速度 [m/s^2]

R_w_val = 0.077;     % 驱动轮半径 [m]              env.py wheel_r=0.77dm
R_l_val = 0.20;      % 驱动轮半轮距 [m]            XML髋关节 Y≈±0.20

l_c_val = 0.10;      % 机体CoM到髋关节距离 [m]      估计值, 需实测


m_w_val = 0.322;     % 单个驱动轮质量 [kg]         wheel_right/left mass="0.322"
m_l_val = 2.751;     % 单腿总质量 [kg]             AG+GH+AB+BE+EC+CF = 0.667+1.64+0.202+0.057+0.128+0.057
m_b_val = 8.4;       % 机体质量 [kg]              base mass="8.4"

% ── 惯量参数 (MuJoCo由mesh自动计算, 未显式给定, 以下为估计) ──
I_w_val = 0.0020;    % 驱动轮转动惯量 [kg·m^2]      1/2*m_w*R_w^2 ≈ 0.001, 略加裕量
% I_l_l / I_l_r 是随腿长变化的符号变量, 不在此赋值, 在循环内由 I_l = k_I*l + b_I 计算
I_b_val = 0.15;      % 机体俯仰转动惯量 [kg·m^2]     1/12*m_b*体长^2 量级
I_z_val = 0.10;      % 机体yaw转动惯量 [kg·m^2]       ≈ I_b * 0.7

%% 随腿长变化的腿部参数近似 (上交文档 WBR_control.html (1.1)~(1.3))
%   (1.1) 单侧腿平面转动惯量  I_l ≈ 0.4  *l + 0.07
%   (1.2) 腿部质心到髋距离    l_b ≈ 0.218*l + 0.075
%   (1.3) 髋到腿部质心距离    l_w = l - l_b ≈ k_lw*l + b_lw
% 斜率/截距用 y = k*x + b 形式表示, 系数需按本模型五连杆实测标定(此处沿用上交开源值)
k_Il  = 0.4;    b_Il  = 0.07; % (1.1)  I_l ≈ k_I*l + b_I
k_lb = 0.218;  b_lb = 0.075; % (1.2)  l_b ≈ k_lb*l + b_lb
k_lw = 1 - k_lb;  b_lw = -b_lb; % (1.3)  l_w = l - l_b ≈ k_lw*l + b_lw

%% 状态空间方程
% 通过建模方程，得到对于每个状态量的符号表达式，再分别对每个状态量方程求各个状态量的偏导，得到10*10的雅可比矩阵A和10*4的雅可比矩阵B，代入线性化点，得到A和B的数值矩阵，最后通过lqr函数计算K矩阵。

equ1 = (I_w*l_l/R_w + m_w*R_w*l_l + m_l*R_w*l_b_l)*ddot_theta_w_l + (m_l*l_w_l*l_b_l - I_l_l)*ddot_theta_l_l + (m_l*l_w_l + 0.5*m_b*l_l)*g*theta_l_l + T_b_l - T_w_l*(1 + l_l/R_w);
equ2 = (I_w*l_r/R_w + m_w*R_w*l_r + m_l*R_w*l_b_r)*ddot_theta_w_r + (m_l*l_w_r*l_b_r - I_l_r)*ddot_theta_l_r + (m_l*l_w_r + 0.5*m_b*l_r)*g*theta_l_r + T_b_r - T_w_r*(1 + l_r/R_w);
equ3 = -(m_w*R_w*R_w + I_w + m_l*R_w*R_w + 0.5*m_b*R_w*R_w)*ddot_theta_w_l - (m_w*R_w*R_w + I_w + m_l*R_w*R_w + 0.5*m_b*R_w*R_w)*ddot_theta_w_r - (m_l*R_w*l_w_l + 0.5*m_b*R_w*l_l)*ddot_theta_l_l - (m_l*R_w*l_w_r + 0.5*m_b*R_w*l_r)*ddot_theta_l_r + T_w_l + T_w_r;
% 第四条方程的“机体俯仰角转动惯量*俯仰角加速度”一项，手写稿是这样写的，html不是，目前按照手写稿写
equ4 = (m_w*R_w*l_c + I_w*l_c/R_w + m_l*R_w*l_c)*ddot_theta_w_l + (m_w*R_w*l_c + I_w*l_c/R_w + m_l*R_w*l_c)*ddot_theta_w_r + m_l*l_w_l*l_c*ddot_theta_l_l + m_l*l_w_r*l_c*ddot_theta_l_r - I_b*ddot_theta_b + m_b*g*l_c*theta_b - (T_w_l + T_w_r)*l_c/R_w - (T_b_l + T_b_r);
equ5 = (0.5*I_z*R_w/R_l + I_w*R_l/R_w)*ddot_theta_w_l - (0.5*I_z*R_w/R_l + I_w*R_l/R_w)*ddot_theta_w_r + 0.5*I_z*l_l/R_l*ddot_theta_l_l - 0.5*I_z*l_r/R_l*ddot_theta_l_r - T_w_l*R_l/R_w + T_w_r*R_l/R_w;

% 解状态量二阶表达式（s和phi的二阶表达式在运动学部分已转化为以下量）
[ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b] = solve([equ1, equ2, equ3, equ4, equ5], [ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b])

% s和phi的运动学二阶表达式
ddot_s = 0.5*R_w*(ddot_theta_w_l + ddot_theta_w_r);
ddot_phi = 0.5*R_w/R_l*(-ddot_theta_w_l + ddot_theta_w_r) - 0.5*l_l/R_l*ddot_theta_l_l + 0.5*l_r/R_l*ddot_theta_l_r;
dot_x = [dot_s, ddot_s, dot_phi, ddot_phi, dot_theta_l_l, ddot_theta_l_l, dot_theta_l_r, ddot_theta_l_r, dot_theta_b, ddot_theta_b];

% 状态向量和控制向量
x = [s, dot_s, phi, dot_phi, theta_l_l, dot_theta_l_l, theta_l_r, dot_theta_l_r, theta_b, dot_theta_b];
u = [T_w_l, T_w_r, T_b_l, T_b_r];

% 状态空间方程的雅可比矩阵
Ja = jacobian(dot_x, x)
Jb = jacobian(dot_x, u)

%% LQR计算与多项式拟合
% 依据《上海交通大学RoboMaster2023平衡步兵控制系统开源》WBR_control.html:
%   (2.1) A(l_l,l_r) = K_A1 + K_A2*l_l + K_A3*l_l^2 + K_A4*l_r + K_A5*l_r^2 + K_A6*l_l*l_r
%   (2.2) B(l_l,l_r) = K_B1 + K_B2*l_l + K_B3*l_l^2 + K_B4*l_r + K_B5*l_r^2 + K_B6*l_l*l_r
%   (3.5) K(l_l,l_r) = K_K1 + K_K2*l_l + K_K3*l_l^2 + K_K4*l_r + K_K5*l_r^2 + K_K6*l_l*l_r
% 三条多项式的系数统一按公式下标顺序存放 (第 1..6 项):
%   1:常数项  2:l_l  3:l_l^2  4:l_r  5:l_r^2  6:l_l*l_r
% 该顺序与上交文档公式下标一一对应, 便于把系数原样导出到嵌入式代码中按同一顺序求值。

%% Step 0: 代入固定物理参数
% 符号变量 g, m_b, ... 未被覆盖, 直接用符号名→数值代入
% 说明: 本模型方程已是小角度线性形式, A、B 与状态量/控制量的值无关, 无需代入 0。
Ja_fixed = subs(Ja, [g, m_b, m_w, m_l, R_w, R_l, l_c, I_w, I_b, I_z], ...
                    [g_val, m_b_val, m_w_val, m_l_val, R_w_val, R_l_val, l_c_val, I_w_val, I_b_val, I_z_val]);
Jb_fixed = subs(Jb, [g, m_b, m_w, m_l, R_w, R_l, l_c, I_w, I_b, I_z], ...
                    [g_val, m_b_val, m_w_val, m_l_val, R_w_val, R_l_val, l_c_val, I_w_val, I_b_val, I_z_val]);

%% LQR 代价函数权重
% J = ∫(x'Qx + u'Ru)dt, Q 为 10 个状态量权重, R 为 4 个控制量权重
% diag(v): 由向量 v 生成对角矩阵(对角元为 v, 其余元素为 0)
Q_cost = diag([1, 1, 1, 1, 10, 1, 10, 1, 100, 10]);
R_cost = diag([1, 1, 1, 1]);

%% 系数顺序转换
% fit(...,'poly22') 返回的曲面模型系数顺序为 [p00 p10 p01 p20 p11 p02],
% 即 [常数, l_l, l_r, l_l^2, l_l*l_r, l_r^2], 与上交文档顺序不同;
% 下面的匿名函数把 6 个系数重排成文档顺序 [1 2 4 3 6 5],
% 即 [常数, l_l, l_l^2, l_r, l_r^2, l_l*l_r]
% 匿名函数写法: f = @(入参) 表达式, 之后可像普通函数一样调用 f(...)
reorder22 = @(c) c([1, 2, 4, 3, 6, 5]);

%% 第一阶段: 采样 (A,B), 拟合 A(l_l,l_r)、B(l_l,l_r)   [文档 (2.1)(2.2)]
l_min = 0.08;   l_max = 0.30;      % 左右腿长工作区间 [m]
step_ab = 0.02;                    % 采样步长 [m]

% meshgrid: 生成所有 (l_l,l_r) 组合; (:): 按列优先拉直成列向量
[LL, LR] = meshgrid(l_min:step_ab:l_max);
LL = LL(:);  LR = LR(:);
n_ab = numel(LL);                  % 采样点总数

% zeros(...): 生成全零数组, 用于预分配内存
% 用三维数组存矩阵(第 3 维是"页"), 避免把矩阵拉直后算线性下标出错
A_s = zeros(n_ab,10,10);           % 第 k 页存一个 10×10 的 A
B_s = zeros(n_ab,10,4);            % 第 k 页存一个 10×4  的 B
for k = 1:n_ab
    li = LL(k);  ri = LR(k);       % 当前左右腿长
    % 由当前腿长算出随腿长变化的 4 个数值参数; 变量名加前缀以区别于符号变量
    % (lb_l≈l_b_l, lw_l≈l_w_l, Il_l≈I_l_l, 切勿写成 l_b_l... 否则会覆盖符号)
    lb_l = k_lb*li + b_lb;   lb_r = k_lb*ri + b_lb; % (1.2)  l_b ≈ k_lb*l + b_lb
    lw_l = li - lb_l;        lw_r = ri - lb_r; % (1.3)  l_w ≈ k_lw*l + b_lw
    Il_l = k_Il*li + b_Il;     Il_r = k_Il*ri + b_Il; % (1.1)  I_l ≈ k_I*l + b_I

    % subs(S,old,new): 用 new 替换符号表达式 S 中的 old(符号变量);
    % double(...): 把符号结果转换成双精度数值矩阵, 整个矩阵直接存入第 k 页
    A_s(k,:,:) = double(subs(Ja_fixed, ...
        [l_l l_r l_w_l l_w_r l_b_l l_b_r I_l_l I_l_r], ...
        [li  ri  lw_l  lw_r  lb_l  lb_r  Il_l  Il_r]));
    B_s(k,:,:) = double(subs(Jb_fixed, ...
        [l_l l_r l_w_l l_w_r l_b_l l_b_r I_l_l I_l_r], ...
        [li  ri  lw_l  lw_r  lb_l  lb_r  Il_l  Il_r]));
end

% 逐元素拟合: 对 A、B 的每个元素在 (l_l,l_r) 平面上做二元二次拟合, 取 6 个系数
% fit([X Y],Z,'poly22'): 曲线拟合工具箱的曲面拟合, X、Y 为两个自变量, Z 为函数值;
%   'poly22' 表示二元二次多项式 p00+p10*x+p01*y+p20*x^2+p11*x*y+p02*y^2
% coeffvalues(fitobj): 取出拟合对象的系数向量
K_A = zeros(10,10,6);    % A 每个元素(r,c)的 6 个多项式系数(文档顺序)
for r = 1:10
    for c = 1:10
        % A_s(:,r,c) 取所有采样点下 A 的第 (r,c) 个元素, 沿腿长做曲面拟合
        K_A(r,c,:) = reorder22(coeffvalues(fit([LL LR], A_s(:,r,c), 'poly22')));
    end
end

K_B = zeros(10,4,6);     % B 每个元素的 6 个多项式系数(文档顺序)
for r = 1:10
    for c = 1:4
        K_B(r,c,:) = reorder22(coeffvalues(fit([LL LR], B_s(:,r,c), 'poly22')));
    end
end
fprintf('阶段一完成: 已拟合 A(l_l,l_r)、B(l_l,l_r) (系数按文档顺序)\n');

%% 第二阶段: 用 A(l_l,l_r)、B(l_l,l_r) 求 K, 再拟合 K(l_l,l_r)   [文档 (3.5)]
step_k = 0.01;                          % 第二阶段网格步长(比第一阶段更密)
[LLk, LRk] = meshgrid(l_min:step_k:l_max);
LLk = LLk(:);  LRk = LRk(:);
n_k = numel(LLk);

K_s = zeros(n_k,4,10);                  % 第 k 页存一个 4×10 的 K
for k = 1:n_k
    li = LLk(k);  ri = LRk(k);

    % 按 (2.1)(2.2) 还原当前腿长下的 A、B (对 10×10 / 10×4 数组整体运算)
    A = K_A(:,:,1) + K_A(:,:,2)*li + K_A(:,:,3)*li^2 + K_A(:,:,4)*ri + K_A(:,:,5)*ri^2 + K_A(:,:,6)*li*ri;
    B = K_B(:,:,1) + K_B(:,:,2)*li + K_B(:,:,3)*li^2 + K_B(:,:,4)*ri + K_B(:,:,5)*ri^2 + K_B(:,:,6)*li*ri;

    % lqrd(A,B,Q,R,Ts): 离散线性二次调节器(discrete LQR)求解, Ts 为采样周期;
    %   返回状态反馈矩阵 K(4×10), 满足 u = -K*x
    K_s(k,:,:) = lqrd(A, B, Q_cost, R_cost, 0.0005);
end

K_K = zeros(4,10,6);                    % K 每个元素(r,c)的 6 个多项式系数(文档顺序)
for r = 1:4
    for c = 1:10
        K_K(r,c,:) = reorder22(coeffvalues(fit([LLk LRk], K_s(:,r,c), 'poly22')));
    end
end
fprintf('阶段二完成: 已拟合 K(l_l,l_r) (系数按文档顺序)\n');

%% 导出系数 (顺序: [常数, l_l, l_l^2, l_r, l_r^2, l_l*l_r])
save('LQR_K_coeff.mat', 'K_A', 'K_B', 'K_K');
fprintf('完成: K_A(10,10,6), K_B(10,4,6), K_K(4,10,6) 已保存到 LQR_K_coeff.mat\n');

% K矩阵
% val(:,:,1) =

%    -0.6995   -1.9929   -0.6809   -0.7028   -1.4450   -0.5642   -1.0666   -0.6765   -9.3259   -1.6775
%    -0.6995   -1.9929    0.6809    0.7028   -1.0666   -0.6765   -1.4450   -0.5642   -9.3259   -1.6775
%    -0.1848   -0.3801    0.1343    0.1574   -0.8057    0.6522    1.7385   -0.1065   -8.9580   -2.0584
%    -0.1848   -0.3801   -0.1343   -0.1574    1.7385   -0.1065   -0.8057    0.6522   -8.9580   -2.0584


% val(:,:,2) =

%    -0.1446    1.4406    0.0834    0.1090  -32.6612   -4.2682    7.5245    1.9569   16.8647    2.8012
%     0.1544    2.2523    0.2101    0.2513  -40.4784   -3.9767   15.2875    3.0059   18.9841    2.8702
%     1.7612    3.3493   -0.3380   -0.3642   63.8512    3.9458  -18.4760   -3.7884  -23.2010   -3.4169
%     0.3318    1.0803    0.0611    0.0962  -27.1526   -1.3697    5.6119    1.3387   11.6010    1.2225


% val(:,:,3) =

%    -0.0596   -2.7796    0.2079    0.1915    6.0403    1.2756   -7.1778   -1.2832   -3.8532   -0.7625
%     0.1654   -2.6716   -0.0218   -0.0640   30.3916    2.9116  -22.6991   -3.5354  -10.0286   -1.5298
%    -4.6490   -9.6031    0.3106    0.3344  -91.6800   -7.6208   52.2895    7.1238   11.8360    2.0601
%     0.3059   -0.0142   -0.1689   -0.2458   59.8841    4.0549  -31.9509   -4.6122  -13.8698   -1.4316


% val(:,:,4) =

%     0.1544    2.2523   -0.2101   -0.2513   15.2875    3.0059  -40.4784   -3.9767   18.9841    2.8702
%    -0.1446    1.4406   -0.0834   -0.1090    7.5245    1.9569  -32.6612   -4.2682   16.8647    2.8012
%     0.3318    1.0803   -0.0611   -0.0962    5.6119    1.3387  -27.1526   -1.3697   11.6010    1.2225
%     1.7612    3.3493    0.3380    0.3642  -18.4760   -3.7884   63.8512    3.9458  -23.2010   -3.4169


% val(:,:,5) =

%     0.1654   -2.6716    0.0218    0.0640  -22.6991   -3.5354   30.3916    2.9116  -10.0286   -1.5298
%    -0.0596   -2.7796   -0.2079   -0.1915   -7.1778   -1.2832    6.0403    1.2756   -3.8532   -0.7625
%     0.3059   -0.0142    0.1689    0.2458  -31.9509   -4.6122   59.8841    4.0549  -13.8698   -1.4316
%    -4.6490   -9.6031   -0.3106   -0.3344   52.2895    7.1238  -91.6800   -7.6208   11.8360    2.0601


% val(:,:,6) =

%    -0.1099   -1.5197    0.0258    0.0283   23.0424    0.3845   41.7791    2.2957  -35.3207   -4.6340
%    -0.1099   -1.5197   -0.0258   -0.0283   41.7791    2.2957   23.0424    0.3845  -35.3207   -4.6340
%    -0.6376   -1.1206   -1.1552   -1.3691  110.9887   13.5065 -144.7826  -14.9866   16.5149    2.3778
%    -0.6376   -1.1206    1.1552    1.3691 -144.7826  -14.9866  110.9887   13.5065   16.5149    2.3778

