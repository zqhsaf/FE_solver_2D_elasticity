function result=reference_basis_2D(xh,yh,basis_type,basis_index,basis_x_der,basis_y_der)
% 参考三角（xi>=0, eta>=0, xi+eta<=1）上的有限元基函数
% basis_type==201：一次三角元（3 个局部基，节点为 3 个顶点）
% basis_type==202：二次三角元（6 个局部基，3 顶点 + 3 边中点）
%   局部节点 1:(0,0) 2:(1,0) 3:(0,1)
%             4:(1/2,0) 边1-2中点  5:(1/2,1/2) 边2-3中点  6:(0,1/2) 边1-3中点
% basis_x_der,basis_y_der：对参考坐标 xh、yh 的偏导阶数
%
% [修改9] 全部改用逐元素运算符（.^ .* ./），使 xh、yh 可以是标量，也可以是
%          同维数组（用于在单元内批量采样估计 L∞ 误差）。标量输入时结果不变。

if basis_type==201
    %% ===================== 一次三角元 =====================
    switch basis_index
        case 1      % 1-xh-yh
            if basis_x_der==0 && basis_y_der==0
                result=1-xh-yh;
            elseif basis_x_der==1 && basis_y_der==0
                result=-1;
            elseif basis_x_der==0 && basis_y_der==1
                result=-1;
            else
                result=0;
            end

        case 2      % xh
            if basis_x_der==0 && basis_y_der==0
                result=xh;
            elseif basis_x_der==1 && basis_y_der==0
                result=1;
            else
                result=0;
            end

        case 3      % yh
            if basis_x_der==0 && basis_y_der==0
                result=yh;
            elseif basis_x_der==0 && basis_y_der==1
                result=1;
            else
                result=0;
            end

        otherwise
            result=0;
    end

elseif basis_type==202
    %% ===================== 二次三角元 =====================
    switch basis_index
        case 1  % 2xh^2+2yh^2+4xh yh-3xh-3yh+1
            if basis_x_der==0 && basis_y_der==0
                result=2*xh.^2+2*yh.^2+4*xh.*yh-3*xh-3*yh+1;
            elseif basis_x_der==1 && basis_y_der==0
                result=4*xh+4*yh-3;
            elseif basis_x_der==0 && basis_y_der==1
                result=4*xh+4*yh-3;
            elseif basis_x_der==2 && basis_y_der==0
                result=4;
            elseif basis_x_der==0 && basis_y_der==2
                result=4;
            elseif basis_x_der==1 && basis_y_der==1
                result=4;
            else
                result=0;
            end

        case 2  % 2xh^2-xh
            if basis_x_der==0 && basis_y_der==0
                result=2*xh.^2-xh;
            elseif basis_x_der==1 && basis_y_der==0
                result=4*xh-1;
            elseif basis_x_der==2 && basis_y_der==0
                result=4;
            else
                result=0;
            end

        case 3  % 2yh^2-yh
            if basis_x_der==0 && basis_y_der==0
                result=2*yh.^2-yh;
            elseif basis_x_der==0 && basis_y_der==1
                result=4*yh-1;
            elseif basis_x_der==0 && basis_y_der==2
                result=4;
            else
                result=0;
            end

        case 4  % -4xh^2-4xh yh+4xh
            if basis_x_der==0 && basis_y_der==0
                result=-4*xh.^2-4*xh.*yh+4*xh;
            elseif basis_x_der==1 && basis_y_der==0
                result=-8*xh-4*yh+4;
            elseif basis_x_der==0 && basis_y_der==1
                result=-4*xh;
            elseif basis_x_der==2 && basis_y_der==0
                result=-8;
            elseif basis_x_der==1 && basis_y_der==1
                result=-4;
            else
                result=0;
            end

        case 5  % 4xh yh
            if basis_x_der==0 && basis_y_der==0
                result=4*xh.*yh;
            elseif basis_x_der==1 && basis_y_der==0
                result=4*yh;
            elseif basis_x_der==0 && basis_y_der==1
                result=4*xh;
            elseif basis_x_der==1 && basis_y_der==1
                result=4;
            else
                result=0;
            end

        case 6  % -4yh^2-4xh yh+4yh
            if basis_x_der==0 && basis_y_der==0
                result=-4*yh.^2-4*xh.*yh+4*yh;
            elseif basis_x_der==1 && basis_y_der==0
                result=-4*yh;
            elseif basis_x_der==0 && basis_y_der==1
                result=-8*yh-4*xh+4;
            elseif basis_x_der==0 && basis_y_der==2
                result=-8;
            elseif basis_x_der==1 && basis_y_der==1
                result=-4;
            else
                result=0;
            end

        otherwise
            result=0;
    end

else
    error('reference_basis_2D：不支持的 basis_type（201 一次，202 二次）。');
end

end
