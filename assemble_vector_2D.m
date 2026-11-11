function b=assemble_vector_2D(f_fun,vector_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,basis_der_x_test,basis_der_y_test)
% 组装二维有限元整体载荷向量（逐单元累加）
% 被积形式：f_fun * 检验基（或其指定阶导数）

b=zeros(vector_size);

% 每个单元统一使用 7 点 5 阶三角高斯公式（Gpn=5）
Gpn=5;

for n=1:number_of_elements

    vertices=P(:, T(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local_tri(vertices,Gpn);

    for beta=1:number_of_local_basis_fun_test

        int_value=Gauss_quad_2D_trial_test_vector(f_fun,Gauss_weights,Gauss_nodes,vertices, ...
            basis_type_test,beta,basis_der_x_test,basis_der_y_test);

        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+int_value;

    end
end

end