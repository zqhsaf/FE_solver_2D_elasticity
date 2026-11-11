function f1=function_f1(x,y)
% 二维弹性问题右端项 f 的第一个分量 f1
% 对应解析解 u1=sin(pi*x)*sin(pi*y)、u2=x*(x-1)*y*(y-1)，lambda=1、mu=2
%   f1 = (lambda+3*mu)*pi^2*sin(pi*x)*sin(pi*y) - (lambda+mu)*(2*x-1)*(2*y-1)
% [修改10] 新增：原代码 function_f 用字符串做算术 ("function_lamda"+3*"function_mu")，
%         且 mu 未定义、返回的是 2×1 向量，无法作为标量被积函数传给高斯积分。

lam=function_lamda(x,y);
mu =function_mu(x,y);

f1=(lam+3*mu)*pi^2*sin(pi*x).*sin(pi*y) - (lam+mu)*(2*x-1).*(2*y-1);

end
