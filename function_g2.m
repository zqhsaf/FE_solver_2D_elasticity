function g2=function_g2(x,y)
% 位移分量 u2 的 Dirichlet 边界值
% 本算例边界上 u2=x*(x-1)*y*(y-1)=0（四边均为齐次边界）
% [修改15] 新增：原因同 function_g1.m

g2=exact_u_2D(x,y,2,0,0);

end
