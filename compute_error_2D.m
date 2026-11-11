function [L2_err,H1_err,Linf_err,L2_rel,H1_rel]=compute_error_2D(solution,P,T,Pb,Tb,basis_type)
% 计算二维弹性有限元数值解相对解析解的误差（一次元 201、二次元 202 通用）
% 输入 solution：2*Nb 维列向量，前 Nb 个为 u1，后 Nb 个为 u2
%
%   L2_err  ：||u-u_h||_0  = sqrt(||u1-u1h||_0^2 + ||u2-u2h||_0^2)
%   H1_err  ： |u-u_h|_1   = sqrt(|u1-u1h|_1^2 + |u2-u2h|_1^2)
%   Linf_err：||u-u_h||_inf = max(||u1-u1h||_inf, ||u2-u2h||_inf)，
%             其中 ||e||_inf = sup_Omega |e|，在每个单元内部细分采样求 sup
%   L2_rel / H1_rel：相应的相对误差（分母取解析解的范数）
%
% L2、H1 的单元积分采用 7 点 5 阶三角高斯公式。
% [修改7] 原函数是标量 Poisson 版本（单个 solution、单个 exact_fun），
%         现改为向量版本：按分量取自由度、分别累加平方再开方；
%         解析解及其导数统一由 exact_u_2D(x,y,comp,der_x,der_y) 提供。
% [修改8] L∞ 由“仅在节点取最大”改为“在 Ω 上求真上确界”：按图片定义
%         ||u-u_h||_inf = max( sup|u1-u1h|, sup|u2-u2h| )，
%         每个单元内部用 (n_sub+1)(n_sub+2)/2 个采样点逼近 sup
%         （n_sub 取偶数可保证包含单元各边中点，线性元误差峰值恰在这些点上）。

switch basis_type
    case 201
        number_of_local_basis=3;
    case 202
        number_of_local_basis=6;
    otherwise
        error('compute_error_2D：不支持的 basis_type（201 一次，202 二次）。');
end

Nb=size(Pb,2);
if numel(solution) ~= 2*Nb
    error('compute_error_2D：solution 的长度应为 2*Nb（两个位移分量）。');
end
sol1=solution(1:Nb);            % u1 的节点值
sol2=solution(Nb+1:end);        % u2 的节点值

Gpn=5;
number_of_elements=size(T,2);

L2_1=0; L2_2=0; H1_1=0; H1_2=0;
L2_den=0; H1_den=0;

for n=1:number_of_elements

    vertices=P(:,T(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local_tri(vertices,Gpn);
    Ng=length(Gauss_weights);

    for k=1:Ng
        x=Gauss_nodes(k,1);
        y=Gauss_nodes(k,2);

        % ---- 数值解及其梯度（局部基线性组合）----
        uh1=0; uh1x=0; uh1y=0;
        uh2=0; uh2x=0; uh2y=0;
        for beta=1:number_of_local_basis
            phi  =FE_basis_local_fun_2D(x,y,vertices,basis_type,beta,0,0);
            phix =FE_basis_local_fun_2D(x,y,vertices,basis_type,beta,1,0);
            phiy =FE_basis_local_fun_2D(x,y,vertices,basis_type,beta,0,1);

            c1=sol1(Tb(beta,n));
            c2=sol2(Tb(beta,n));

            uh1 =uh1 +c1*phi;   uh1x=uh1x+c1*phix;   uh1y=uh1y+c1*phiy;
            uh2 =uh2 +c2*phi;   uh2x=uh2x+c2*phix;   uh2y=uh2y+c2*phiy;
        end

        % ---- 解析解及其梯度 ----
        u1 =exact_u_2D(x,y,1,0,0);   u1x=exact_u_2D(x,y,1,1,0);   u1y=exact_u_2D(x,y,1,0,1);
        u2 =exact_u_2D(x,y,2,0,0);   u2x=exact_u_2D(x,y,2,1,0);   u2y=exact_u_2D(x,y,2,0,1);

        L2_1=L2_1+Gauss_weights(k)*(uh1-u1)^2;
        L2_2=L2_2+Gauss_weights(k)*(uh2-u2)^2;
        H1_1=H1_1+Gauss_weights(k)*((uh1x-u1x)^2+(uh1y-u1y)^2);
        H1_2=H1_2+Gauss_weights(k)*((uh2x-u2x)^2+(uh2y-u2y)^2);

        L2_den=L2_den+Gauss_weights(k)*(u1^2+u2^2);
        H1_den=H1_den+Gauss_weights(k)*(u1x^2+u1y^2+u2x^2+u2y^2);
    end
end

L2_err =sqrt(L2_1+L2_2);
H1_err =sqrt(H1_1+H1_2);
L2_rel =L2_err/sqrt(L2_den);
H1_rel =H1_err/sqrt(H1_den);

%% ---- L∞ 误差：Ω 上的真上确界（逐单元细分采样）----
n_sub=20;                                   % 每个单元每条参考方向细分数
[II,JJ]=meshgrid(0:n_sub,0:n_sub);
mask=(II+JJ)<=n_sub;                        % 只保留参考单元内的采样点
xi_s=(II(mask)/n_sub).';
eta_s=(JJ(mask)/n_sub).';

sup1=0; sup2=0;
for n=1:number_of_elements

    vertices=P(:,T(:,n));

    % 参考坐标 -> 物理坐标（与 FE_basis_local_fun_2D 中的仿射映射一致）
    X=vertices(1,1)+(vertices(1,2)-vertices(1,1))*xi_s+(vertices(1,3)-vertices(1,1))*eta_s;
    Y=vertices(2,1)+(vertices(2,2)-vertices(2,1))*xi_s+(vertices(2,3)-vertices(2,1))*eta_s;

    uh1=0; uh2=0;
    for beta=1:number_of_local_basis
        phi=FE_basis_local_fun_2D(X,Y,vertices,basis_type,beta,0,0);
        uh1=uh1+sol1(Tb(beta,n))*phi;
        uh2=uh2+sol2(Tb(beta,n))*phi;
    end

    sup1=max(sup1,max(abs(uh1-exact_u_2D(X,Y,1,0,0))));
    sup2=max(sup2,max(abs(uh2-exact_u_2D(X,Y,2,0,0))));
end

Linf_err=max(sup1,sup2);

end
