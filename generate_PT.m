function [P,T]=generate_PT(left,right,bottom,top,N1,N2)
% 生成结构化矩形-三角网格
%   P：2×Nm 节点坐标矩阵，每列一个节点
%   T：3×N  单元顶点编号表，每列一个三角单元
% 编号规则（先 y 后 x）：节点全局编号 k = i*(N2+1)+j+1，
%   i=0..N1 为 x 方向索引，j=0..N2 为 y 方向索引。
% 每个小矩形沿“左下-右上”对角线剖分为两个逆时针三角形。

h1=(right-left)/N1;
h2=(top-bottom)/N2;
Nm=(N1+1)*(N2+1);
N=2*N1*N2;

%% 1. 节点坐标
P=zeros(2,Nm);
for i=0:N1
    for j=0:N2
        k=i*(N2+1)+j+1;
        P(1,k)=left+i*h1;
        P(2,k)=bottom+j*h2;
    end
end

%% 2. 三角单元连接表
T=zeros(3,N);
for i=0:N1-1
    for j=0:N2-1
        v1=i*(N2+1)+j+1;        % 矩形左下顶点
        v2=(i+1)*(N2+1)+j+1;    % 矩形右下顶点
        v3=i*(N2+1)+j+2;        % 矩形左上顶点
        v4=(i+1)*(N2+1)+j+2;    % 矩形右上顶点

        n=2*(i*N2+j)+1;
        T(:,n)  =[v1;v2;v3];    % 左下三角形（逆时针）
        T(:,n+1)=[v3;v2;v4];    % 右上三角形（逆时针）
    end
end

end