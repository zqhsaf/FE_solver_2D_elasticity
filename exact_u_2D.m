function value=exact_u_2D(x,y,comp,der_x,der_y)
% 二维弹性算例的解析解及其偏导数
%   comp=1： u1 = sin(pi*x)*sin(pi*y)
%   comp=2： u2 = x*(x-1)*y*(y-1)
% der_x,der_y：对 x、y 的偏导阶数（只支持 0/1 阶，二阶返回 0）
% 支持 x、y 为同维数组（逐元素运算），便于在节点上向量化求值。
% [修改13] 新增：解析解按分量/导数阶次统一取值，供 compute_error_2D 使用。

if comp==1
    if der_x==0 && der_y==0
        value=sin(pi*x).*sin(pi*y);
    elseif der_x==1 && der_y==0
        value=pi*cos(pi*x).*sin(pi*y);
    elseif der_x==0 && der_y==1
        value=pi*sin(pi*x).*cos(pi*y);
    else
        value=zeros(size(x));
    end

elseif comp==2
    if der_x==0 && der_y==0
        value=x.*(x-1).*y.*(y-1);
    elseif der_x==1 && der_y==0
        value=(2*x-1).*y.*(y-1);
    elseif der_x==0 && der_y==1
        value=x.*(x-1).*(2*y-1);
    else
        value=zeros(size(x));
    end

else
    error('exact_u_2D：comp 只能取 1（u1）或 2（u2）。');
end

end
