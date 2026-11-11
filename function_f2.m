function f2=function_f2(x,y)
% 二维弹性问题右端项 f 的第二个分量 f2
% 对应解析解 u1=sin(pi*x)*sin(pi*y)、u2=x*(x-1)*y*(y-1)，lambda=1、mu=2
%   f2 = -(lambda+2*mu)*2*x*(x-1) - (lambda+mu)*pi^2*cos(pi*x)*cos(pi*y)
%        - mu*2*y*(y-1)
% [修改11] 新增：原因同 function_f1.m

lam=function_lamda(x,y);
mu =function_mu(x,y);

f2=-(lam+2*mu).*(2*x.*(x-1)) ...
   -(lam+mu)*pi^2*cos(pi*x).*cos(pi*y) ...
   -mu.*(2*y.*(y-1));

end
