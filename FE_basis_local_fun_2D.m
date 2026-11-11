function result=FE_basis_local_fun_2D(x,y,vertices,basis_type,basis_index,basis_x_der,basis_y_der)
% 物理三角单元上的局部有限元基函数（及其指定阶偏导数）
% basis_type==201：一次三角元（局部基 1/2/3，仅顶点）
% basis_type==202：二次三角元（局部基 1..6，顶点+边中点）
% vertices：2×3 顶点矩阵，每列为一个顶点 [x;y]
% basis_x_der,basis_y_der：对物理坐标 x、y 的偏导阶数
% 几何上一次、二次元均采用相同的顶点仿射映射，故坐标变换与链式法则通用。

if basis_type~=201 && basis_type~=202
    error('FE_basis_local_fun_2D：不支持的 basis_type（201 一次，202 二次）。');
end

xn1=vertices(1,1); xn2=vertices(1,2); xn3=vertices(1,3);
yn1=vertices(2,1); yn2=vertices(2,2); yn3=vertices(2,3);

% 仿射映射的雅可比行列式（逆时针顶点时为正，等于单元面积的 2 倍）
Jn=(xn2-xn1)*(yn3-yn1)-(xn3-xn1)*(yn2-yn1);

% 物理坐标 (x,y) -> 参考坐标 (xh,yh)
xh=((yn3-yn1)*(x-xn1)-(xn3-xn1)*(y-yn1))/Jn;
yh=(-(yn2-yn1)*(x-xn1)+(xn2-xn1)*(y-yn1))/Jn;

if basis_x_der==0 && basis_y_der==0
    % 基函数本身
    result=reference_basis_2D(xh,yh,basis_type,basis_index,0,0);

elseif basis_x_der==1 && basis_y_der==0
    % 对物理坐标 x 的偏导（链式法则）
    result=((yn3-yn1)/Jn)*reference_basis_2D(xh,yh,basis_type,basis_index,1,0) ...
          +((yn1-yn2)/Jn)*reference_basis_2D(xh,yh,basis_type,basis_index,0,1);

elseif basis_x_der==0 && basis_y_der==1
    % 对物理坐标 y 的偏导（链式法则）
    result=((xn1-xn3)/Jn)*reference_basis_2D(xh,yh,basis_type,basis_index,1,0) ...
          +((xn2-xn1)/Jn)*reference_basis_2D(xh,yh,basis_type,basis_index,0,1);
else
    result=0;
end

end