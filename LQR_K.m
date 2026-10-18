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

% 解状态量二阶表达式（s和phi的二阶表达式在运动学部分已转化为以下量）
[ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b, ] = solve([equ1, equ2, equ3, equ4, equ5], [ddot_theta_w_l, ddot_theta_w_r, ddot_theta_l_l, ddot_theta_l_r, ddot_theta_b])

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

%% 物理参数赋值
% ─── 车奉朝轮腿模型物理参数 ────────────────────────────
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

%% LQR计算和拟合

% Step 0: 代入固定物理参数 (给除了和腿长相关的参数赋值）
% 剩余符号: l_l, l_r, l_w_l, l_w_r, l_b_l, l_b_r, I_l_l, I_l_r
params_fixed = [g, m_b, m_w, m_l, R_w, R_l, l_c, I_w, I_b, I_z];
Ja_fixed = subs(Ja, params_fixed);
Jb_fixed = subs(Jb, params_fixed);


% 腿长相关参数用 y = kx + b 线性近似 (上交 (1.1)~(1.3), 系数需根据本模型五连杆标定)
%   I_l(l) = k_I * l + b_I
%   l_b(l) = k_lb * l + b_lb
%   l_w(l) = l - l_b = (1 - k_lb)*l - b_lb
k_I  = 0.4;    b_I  = 0.07;
k_lb = 0.218;  b_lb = 0.075;

% Q/R 权重矩阵
Q_cost = diag([1, 1, 1, 1, 10, 1, 10, 1, 100, 10]);
R_cost = diag([1, 1, 1, 1]);

% ═══════════════════════════════════════════════════════════════
% 第一阶段：拟合 A(l_l,l_r) 和 B(l_l,l_r) 多项式 （上交 2.1）
%   多项式形式（每个矩阵元素 6 个系数）：
%     A(l_l,l_r) = K_A1 + K_A2*l_l + K_A3*l_l² + K_A4*l_r + K_A5*l_r² + K_A6*l_l*l_r
%     B(l_l,l_r) = K_B1 + K_B2*l_l + K_B3*l_l² + K_B4*l_r + K_B5*l_r² + K_B6*l_l*l_r
%   流程：左右腿长范围内按步长全排列 → 每点算 A,B 数值 → 最小二乘拟合 6 系数
% ═══════════════════════════════════════════════════════════════

% 腿长范围（左右独立，便于后期不对称扩展）
l_l_min = 0.08;   l_l_max = 0.30;   % 左腿长范围 [m]
l_r_min = 0.08;   l_r_max = 0.30;   % 右腿长范围 [m]
leg_step = 0.02;                      % 离散A/B采样步长 [m]

% 生成左右腿长的全排列网格 → 每个 (l_l, l_r) 组合都参与A/B拟合
l_l_vec = l_l_min:leg_step:l_l_max;
l_r_vec = l_r_min:leg_step:l_r_max;
[l_l_grid, l_r_grid] = meshgrid(l_l_vec, l_r_vec);   % 2D网格, 左右腿长所有排列组合
l_l_flat = l_l_grid(:);                              % 展开为列向量
l_r_flat = l_r_grid(:);
n_pts = length(l_l_flat);                            % 总采样点数

% 存储每点的 A(:) 和 B(:) —— 每行是一个 (l_l,l_r) 对应的矩阵展开
A_stack = zeros(n_pts, 100);   % 10×10 展开为 100 列
B_stack = zeros(n_pts, 40);    % 10×4 展开为 40 列

fprintf('=== 阶段一：离散采样 %d 个 (l_l,l_r) 点, 计算 A/B ===\n', n_pts);
for i = 1:n_pts
    li = l_l_flat(i);          % 当前左腿长
    ri = l_r_flat(i);          % 当前右腿长

    % 由腿长计算腿的质心位置、转动惯量 (y = kx + b 线性近似)
    l_b_l = k_lb*li + b_lb;    l_b_r = k_lb*ri + b_lb;   % 腿质心位置
    l_w_l = li - l_b_l;        l_w_r = ri - l_b_r;        % 轮心到腿质心
    I_l_l = k_I*li + b_I;      I_l_r = k_I*ri + b_I;      % 腿转动惯量

    % 代入腿长相关参数 → 符号矩阵变为纯数值
    A_num = double(subs(Ja_fixed, ...
        [l_l, l_r, l_w_l, l_w_r, l_b_l, l_b_r, I_l_l, I_l_r], ...
        [li, ri, l_w_l, l_w_r, l_b_l, l_b_r, I_l_l, I_l_r]));
    B_num = double(subs(Jb_fixed, ...
        [l_l, l_r, l_w_l, l_w_r, l_b_l, l_b_r, I_l_l, I_l_r], ...
        [li, ri, l_w_l, l_w_r, l_b_l, l_b_r, I_l_l, I_l_r]));

    A_stack(i, :) = A_num(:)';     % 展开存为一行
    B_stack(i, :) = B_num(:)';
end

% 构造 2D 二次多项式设计矩阵 X: 每行 = [1, l_l, l_l², l_r, l_r², l_l·l_r]
% 对矩阵每个元素做最小二乘: coeff = X \ y, 得到该元素的 6 个多项式系数
X_ab = [ones(n_pts,1), l_l_flat, l_l_flat.^2, l_r_flat, l_r_flat.^2, l_l_flat.*l_r_flat];

fprintf('  对 A, B 每个元素做最小二乘拟合 (每个元素 → 6 系数)...\n');
K_A = zeros(10, 10, 6);   % A 矩阵每个元素的 6 个多项式系数
K_B = zeros(10, 4, 6);    % B 矩阵每个元素的 6 个多项式系数

for row = 1:10
    for col = 1:10
        idx = (row-1)*10 + col;           % 矩阵元素 → 在展开向量中的列索引
        K_A(row, col, :) = X_ab \ A_stack(:, idx);  % 最小二乘解 6 系数
    end
end
for row = 1:10
    for col = 1:4
        idx = (row-1)*4 + col;
        K_B(row, col, :) = X_ab \ B_stack(:, idx);
    end
end
fprintf('  阶段一完成: K_A(10×10×6), K_B(10×4×6)\n');

% ═══════════════════════════════════════════════════════════════
% 第二阶段：用拟合的 A(l_l,l_r), B(l_l,l_r) 在稠密网格上计算 K（上交 3.1）
%   多项式形式（K 的每个元素 6 个系数）：
%     K(l_l,l_r) = K_K1 + K_K2*l_l + K_K3*l_l² + K_K4*l_r + K_K5*l_r² + K_K6*l_l*l_r
%   流程：左右腿长全排列（步长低于阶段一）→ 每点用A/B多项式快速求A,B → lqrd得K → 最小二乘拟合K
% ═══════════════════════════════════════════════════════════════

leg_step_k = 0.01;                                   % K采样步长 [m], 比阶段一的 leg_step(0.02) 更密
% 生成左右腿长的全排列网格（步长更密，点数更多）
l_l_k_vec = l_l_min:leg_step_k:l_l_max;
l_r_k_vec = l_r_min:leg_step_k:l_r_max;
[l_l_k_grid, l_r_k_grid] = meshgrid(l_l_k_vec, l_r_k_vec);  % 所有 (l_l,l_r) 排列组合
l_l_k_flat = l_l_k_grid(:);
l_r_k_flat = l_r_k_grid(:);
n_k_pts = length(l_l_k_flat);

K_stack = zeros(n_k_pts, 40);   % 4×10 展开为 40 列, 每行一个 (l_l,l_r)

fprintf('=== 阶段二：稠密采样 %d 个 (l_l,l_r) 点, 用 A,B 多项式快速计算 K ===\n', n_k_pts);
for i = 1:n_k_pts
    li = l_l_k_flat(i);         % 当前左腿长
    ri = l_r_k_flat(i);         % 当前右腿长

    % 用阶段一拟合的 K_A,K_B 多项式，快速计算当前腿长下的 A, B 矩阵
    % basis = [1, l_l, l_l², l_r, l_r², l_l·l_r], 点乘 6 系数即得该元素值
    basis = [1, li, li^2, ri, ri^2, li*ri];
    A_fit = zeros(10, 10);
    B_fit = zeros(10, 4);
    for row = 1:10
        for col = 1:10
            A_fit(row, col) = sum(squeeze(K_A(row, col, :))' .* basis);  % A(row,col) = Σ 系数_i * basis_i
        end
    end
    for row = 1:10
        for col = 1:4
            B_fit(row, col) = sum(squeeze(K_B(row, col, :))' .* basis);
        end
    end

    % 用重构的 A, B 计算离散 LQR 增益 (Ts=0.001)
    K_mat = lqrd(A_fit, B_fit, Q_cost, R_cost, 0.001);
    K_stack(i, :) = K_mat(:)';    % 4×10 → 1×40 行向量, 存入一行
end

% 对 K 的每个元素做最小二乘拟合：K(row,col) = f(l_l,l_r) 的 6 系数
% 设计矩阵与阶段一相同: [1, l_l, l_l², l_r, l_r², l_l·l_r]
X_k = [ones(n_k_pts,1), l_l_k_flat, l_l_k_flat.^2, l_r_k_flat, l_r_k_flat.^2, l_l_k_flat.*l_r_k_flat];
K_K = zeros(4, 10, 6);    % K 矩阵每个元素的 6 个多项式系数
for row = 1:4
    for col = 1:10
        idx = (row-1)*10 + col;           % K(row,col) → 在 K_stack 中的列号
        K_K(row, col, :) = X_k \ K_stack(:, idx);  % 最小二乘解 6 系数
    end
end
fprintf('  阶段二完成: K_K(4×10×6)\n');

% ═══════════════════════════════════════════════
% 验证：用 K_K 多项式重建某腿长下的 K 矩阵
% ═══════════════════════════════════════════════
l_test = 0.20;
basis_test = [1, l_test, l_test^2, l_test, l_test^2, l_test*l_test];
K_test = zeros(4, 10);
for row = 1:4
    for col = 1:10
        K_test(row, col) = sum(squeeze(K_K(row, col, :))' .* basis_test);
    end
end
fprintf('\n腿长 l=%.2f m 时的 LQR 反馈增益 K (4 控制 × 10 状态):\n', l_test);
disp(K_test);