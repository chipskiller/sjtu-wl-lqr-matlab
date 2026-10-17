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

% 参数
syms g; % 重力加速度
syms R_w; % 驱动轮半径
syms R_l; % 驱动轮轮距/2
syms l_l l_r; % 腿部长度
syms l_w_l l_w_r; % 驱动轮到腿部质心距离*1
syms l_b_l l_b_r; % 机体质心到腿部质心距离*1
syms l_c; % 机体质心到腿部关节距离
syms m_w; % 驱动轮质量
syms m_l; % 腿部质量
syms m_b; % 机体质量
syms I_w; % 驱动轮转动惯量
syms I_l_l I_l_r; % 腿部转动惯量
syms I_b; % 机体转动惯量
syms I_z; % 机体yaw转动惯量

%% 状态空间方程
% 通过建模方程，得到对于每个状态量的符号表达式，再分别对每个状态量方程求各个状态量的偏导，得到10*10的雅可比矩阵A和10*4的雅可比矩阵B，代入线性化点，得到A和B的数值矩阵，最后通过lqr函数计算K矩阵。

equ1 = (I_w*l_l/R_w + m_w*R_w*l_l + m_l*R_w*l_b_l)*ddot_theta_w_l + (m_l*l_w_l*l_b_l - I_l_l)*ddot_theta_l_l + (m_l*l_w_l + 0.5*m_b*l_l)*g*theta_l_l + T_b_l - T_w_l*(1 + l_l/R_w);
equ2 = (I_w*l_r/R_w + m_w*R_w*l_r + m_l*R_w*l_b_r)*ddot_theta_w_r + (m_l*l_w_r*l_b_r - I_l_r)*ddot_theta_l_r + (m_l*l_w_r + 0.5*m_b*l_r)*g*theta_l_r + T_b_r - T_w_r*(1 + l_r/R_w);
equ3 = -(m_w*R_w*R_w + I_w + m_l*R_w*R_w + 0.5*m_b*R_w*R_w)*ddot_theta_w_l - (m_w*R_w*R_w + I_w + m_l*R_w*R_w + 0.5*m_b*R_w*R_w)*ddot_theta_w_r - (m_l*R_w*l_w_l + 0.5*m_b*R_w*l_l)*ddot_theta_l_l - (m_l*R_w*l_w_r + 0.5*m_b*R_w*l_r)*ddot_theta_l_r + T_w_l + T_w_r;
% 第四条方程的“机体俯仰角转动惯量*俯仰角加速度”一项，手写稿是这样写的，html不是，目前按照手写稿写
equ4 = (m_w*R_w*l_c + I_w*l_c/R_w + m_l*R_w*l_c)*ddot_theta_w_l + (m_w*R_w*l_c + I_w*l_c/R_w + m_l*R_w*l_c)*ddot_theta_w_r + m_l*l_w_l*l_c*ddot_theta_l_l + m_l*l_w_r*l_c*ddot_theta_l_r - I_b*ddot_theta_b + m_b*g*l_c*theta_b - (T_w_l + T_w_r)*l_c/R_w - (T_b_l + T_b_r);
equ5 = (0.5*I_z*R_w/R_l + I_w*R_l/R_w)*ddot_theta_w_l - (0.5*I_z*R_w/R_l + I_w*R_l/R_w)*ddot_theta_w_r + 0.5*I_z*l_l/R_l*ddot_theta_l_l - 0.5*I_z*l_r/R_l*ddot_theta_l_r - T_w_l*R_l/R_w + T_w_r*R_l/R_w;

% s和phi的运动学二阶表达式
ddot_s = 0.5*R_w*(ddot_theta_w_l + ddot_theta_w_r);
ddot_phi = 0.5*R_w/R_l*(-ddot_theta_w_l + ddot_theta_w_r) - 0.5*l_l/R_l*ddot_theta_l_l + 0.5*l_r/R_l*ddot_theta_l_r;

% 解状态量二阶表达式（s和phi的二阶表达式在运动学部分已转化为以下量）
[ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b, ] = solve([equ1, equ2, equ3, equ4, equ5], [ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b])
dot_x = [dot_s, ddot_s, dot_phi, ddot_phi, dot_theta_l_l, ddot_theta_l_l, dot_theta_l_r, ddot_theta_l_r, dot_theta_b, ddot_theta_b];

% 状态向量和控制向量
x = [s, dot_s, phi, dot_phi, theta_l_l, dot_theta_l_l, theta_l_r, dot_theta_l_r, theta_b, dot_theta_b];
u = [T_w_l, T_w_r, T_b_l, T_b_r];

% 状态空间方程的雅可比矩阵
Ja = jacobian(dot_x, x)
Jb = jacobian(dot_x, u)

%% ─── 车奉朝轮腿模型物理参数 ────────────────────────────
%  来源: robot.xml  +  env.py + 估算
% ─────────────────────────────────────────────────────────

g   = 9.81;      % 重力加速度 [m/s^2]
m_b = 8.4;       % 机体质量 [kg]              base mass="8.4"
m_w = 0.322;     % 单个驱动轮质量 [kg]         wheel_right/left mass="0.322"
m_l = 2.751;     % 单腿总质量 [kg]             AG+GH+AB+BE+EC+CF = 0.667+1.64+0.202+0.057+0.128+0.057
R_w = 0.077;     % 驱动轮半径 [m]              env.py wheel_r=0.77dm
R_l = 0.20;      % 驱动轮半轮距 [m]            XML髋关节 Y≈±0.20
l_c = 0.10;      % 机体CoM到髋关节距离 [m]      估计值, 需实测

% ── 惯量参数 (MuJoCo由mesh自动计算, 未显式给定, 以下为估计) ──
I_w = 0.0020;    % 驱动轮转动惯量 [kg·m^2]      1/2*m_w*R_w^2 ≈ 0.001, 略加裕量
I_l_l = 0.015;   % 左腿绕髋转动惯量 [kg·m^2]     m_l*L^2/3 量级
I_l_r = 0.015;   % 右腿绕髋转动惯量 [kg·m^2]     m_l*L^2/3 量级
I_b = 0.15;      % 机体俯仰转动惯量 [kg·m^2]     1/12*m_b*体长^2 量级
I_z = 0.10;      % 机体yaw转动惯量 [kg·m^2]       ≈ I_b * 0.7

%% LQR计算

Q_cost = diag([0 0 0 0 0 0 0 0 0 0]); % 状态成本矩阵
R_cost = diag([0 0 0 0]); % 控制成本矩阵

leg_len = 0.096; % 腿长
K = zeros(40, 40); % 40组，每组10状态量*4控制量
leglen_arr = zeros(40, 2); % 腿长数组
for i = 1:40
    leg_var = leg_var 
end