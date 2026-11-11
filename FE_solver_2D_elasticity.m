function [solution,P,T,Pb,Tb]=FE_solver_2D_elasticity(left,right,bottom,top,basis_type_trial,basis_type_test,N1,N2)
% 二维线弹性方程组  -div(sigma(u)) = f  的三角元有限元求解
% 求解区域：[left,right] × [bottom,top]，四条边均为 Dirichlet 边界 u1=u2=0
% 位移为二维向量 u=(u1,u2)^T，应力 sigma(u) = lambda*(div u)*I + 2*mu*eps(u)
% basis_type_trial / basis_type_test：
%   201 一次三角元（每单元 3 个局部基）；202 二次三角元（每单元 6 个局部基）
% N1,N2：x、y 方向上的均匀剖分段数
%
% 输出：
%   solution：2*Nb 维列向量，前 Nb 个分量为 u1 在有限元节点上的值，
%             后 Nb 个分量为 u2 在有限元节点上的值
%   P,T,Pb,Tb：网格数据（供 compute_error_2D 计算误差使用）
%                                                          % [修改1] 增加网格输出

% [修改5] 增加一致性检查：本求解器在同一有限元空间内做 Galerkin 离散，
%          试探基与检验基必须相同（否则节点表、自由度编号无法对应）。
if basis_type_trial ~= basis_type_test
    error('主求解器：basis_type_trial 与 basis_type_test 必须相同（Galerkin 离散使用同一空间）。');
end

%% 1. 生成网格与有限元节点
[P,T]=generate_PT(left,right,bottom,top,N1,N2);
[Pb,Tb]=generate_PbTb(left,right,bottom,top,N1,N2,basis_type_trial);

switch basis_type_trial
    case 201
        Pb_trial=Pb; Tb_trial=Tb; number_of_local_basis_fun_trial=3;
    case 202
        Pb_trial=Pb; Tb_trial=Tb; number_of_local_basis_fun_trial=6;
    otherwise
        error('主求解器：不支持的试探基类型（201 一次，202 二次）。');
end

switch basis_type_test
    case 201
        Pb_test=Pb; Tb_test=Tb; number_of_local_basis_fun_test=3;
    case 202
        Pb_test=Pb; Tb_test=Tb; number_of_local_basis_fun_test=6;
    otherwise
        error('主求解器：不支持的检验基类型（201 一次，202 二次）。');
end

%% 2. 生成边界节点（基于坐标，一次/二次元通用）
boundarynodes=generate_boundarynodes(Pb_test,left,right,bottom,top);

Nb=size(Pb_test,2);                     % 每个位移分量的自由度个数   % [修改2]
matrix_size=[Nb,Nb];
number_of_elements=size(T,2);

%% 3. 组装刚度矩阵
% 双线性形式： a(u,v) = ∫ lambda*(div u)*(div v) + 2*mu*eps(u):eps(v) dx
%              = ∫ lambda*(u1,x+u2,y)*(v1,x+v2,y)
%                + 2*mu*( u1,x*v1,x + u2,y*v2,y + (u1,y+u2,x)*(v1,y+v2,x)/2 )
% 记 Iab(i,j) = ∫ phi_j,a * psi_i,b （phi 为试探基，psi 为检验基，i 检验、j 试探），
% 则四个分块为：
%   A11 = (lambda+2*mu)*Ixx + mu*Iyy
%   A12 = lambda*Ixy  + mu*Iyx        (试探为 u1、检验为第 2 个方程)
%   A21 = lambda*Iyx  + mu*Ixy        (试探为 u2、检验为第 1 个方程)
%   A22 = (lambda+2*mu)*Iyy + mu*Ixx
A1=assemble_matrix_2D('function_lamda',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,1,0,basis_type_test,1,0);          % lambda*Ixx
A2=assemble_matrix_2D('function_mu',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,1,0,basis_type_test,1,0);          % mu*Ixx
A3=assemble_matrix_2D('function_mu',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,1,basis_type_test,0,1);          % mu*Iyy
A4=assemble_matrix_2D('function_lamda',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,1,basis_type_test,1,0);          % lambda*Iyx
A5=assemble_matrix_2D('function_mu',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,1,0,basis_type_test,0,1);          % mu*Ixy
A6=assemble_matrix_2D('function_lamda',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,1,0,basis_type_test,0,1);          % lambda*Ixy
A7=assemble_matrix_2D('function_mu',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,1,basis_type_test,1,0);          % mu*Iyx
A8=assemble_matrix_2D('function_lamda',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,1,basis_type_test,0,1);          % lambda*Iyy
% [修改3] 修正两个非对角块的导数配对：
%   原代码为 A=[A1+2*A2+A3  A4+A5; A6+A7  A8+2*A3+A2]，即 A12=lambda*Iyx+mu*Ixy、
%   A21=lambda*Ixy+mu*Iyx，把 lambda 与 mu 所乘的积分张冠李戴（lambda≠mu 时结果错误）。
%   正确配对：A12=lambda*Ixy+mu*Iyx=A6+A7，A21=lambda*Iyx+mu*Ixy=A4+A5。
A=[A1+2*A2+A3, A6+A7; A4+A5, A8+2*A3+A2];

%% 4. 组装载荷向量 b = ∫ f*psi
% 两个位移分量分别对应 f1、f2（f 为向量函数，必须按分量分别组装）
vector_size=[Nb,1];
b1=assemble_vector_2D('function_f1',vector_size,number_of_elements,P,T, ...
    Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);
b2=assemble_vector_2D('function_f2',vector_size,number_of_elements,P,T, ...
    Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);
b=[b1;b2];

%% 5. 处理 Dirichlet 边界条件（向量值：u1 与 u2 均为 0）
% [修改4] 原调用只传了一个标量边界函数 'function_g'，无法区分两个位移分量；
%         改为分别传入 u1、u2 的边界函数 function_g1 / function_g2。
[A,b]=treat_Dirichlet_boundary(A,b,boundarynodes,Pb_test,'function_g1','function_g2',-1);

%% 6. 求解整体线性方程组
solution=A\b;

end
