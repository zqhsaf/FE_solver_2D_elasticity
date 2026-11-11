function boundarynodes=generate_boundarynodes(Pb,left,right,bottom,top)
% 基于节点坐标识别矩形四条边上的全部边界节点（一次元、二次元通用）
% 输出 boundarynodes 为 2 行矩阵：
%   第 1 行：边界类型标记，-1 表示 Dirichlet
%   第 2 行：边界节点在 Pb 中的全局编号
% 二次元时边中点同样落在边界上，会被一并识别。

tol=1e-10;
idx=find( abs(Pb(1,:)-left)<tol  | abs(Pb(1,:)-right)<tol  | ...
         abs(Pb(2,:)-bottom)<tol| abs(Pb(2,:)-top)<tol );

boundarynodes=[-ones(1,length(idx)); idx];

end