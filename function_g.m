function g=function_g(x,y)
% 二维弹性问题 Dirichlet 边界上的位移向量 g=[g1;g2]（本算例两者均为 0）
% 仅用于检查/画图；求解时由 function_g1、function_g2 分别给定两个分量。
% [修改16] 原代码两个分量的判断完全重复，且注释里的边界坐标（-1）与实际区域（0,1）不符。

g=[function_g1(x,y); function_g2(x,y)];

end
