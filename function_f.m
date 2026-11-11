function f=function_f(x,y)
% 二维弹性问题右端项向量 f=[f1;f2]（仅用于检查/画图，求解器按分量分别组装）


f1=function_f1(x,y);
f2=function_f2(x,y);

f=[f1;f2];

end
