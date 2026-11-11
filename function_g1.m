function g1=function_g1(x,y)
% 位移分量 u1 的 Dirichlet 边界值
% 本算例边界上 u1=sin(pi*x)*sin(pi*y)=0（四边均为齐次边界）
% [修改14] 新增：原 function_g 中两个分量的 if-else 完全重复，
%          且返回 [g1;g2] 向量，无法作为标量边界函数使用。

g1=exact_u_2D(x,y,1,0,0);

end
