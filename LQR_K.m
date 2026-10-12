% 计算K矩阵
clear;

%% 符号定义

% 变量
syms theta_w_l theta_w_r; % 驱动轮转角
syms theta_l_l theta_l_r; % 腿倾转角
syms theta_b; % 机体倾斜角
syms phi; % 偏航角
syms s; % 自然坐标系下机器人水平方向移动距离
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
syms I_i; % 腿部转动惯量
syms I_b; % 机体转动惯量
syms I_z; % 机体yaw转动惯量

g = 9.81;
