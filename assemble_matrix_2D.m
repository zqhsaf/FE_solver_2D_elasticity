function A=assemble_matrix_2D(coe_fun,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_der_x_test,basis_der_y_test)
% 组装二维有限元整体矩阵（逐单元累加）
% 被积形式：coe_fun * (D_trial 试探基) * (D_test 检验基)
% 导数阶次由 basis_der_x/y_trial、basis_der_x/y_test 指定

A=sparse(matrix_size(1),matrix_size(2));

% 每个单元统一使用 7 点 5 阶三角高斯公式（Gpn=5）
Gpn=5;

for n=1:number_of_elements

    vertices=P(:, T(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local_tri(vertices,Gpn);

    for alpha=1:number_of_local_basis_fun_trial
        for beta=1:number_of_local_basis_fun_test

            int_value=Gauss_quad_2D_trial_test_matrix(coe_fun,Gauss_weights,Gauss_nodes,vertices, ...
                basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial, ...
                basis_type_test,beta,basis_der_x_test,basis_der_y_test);

            i=Tb_test(beta,n);
            j=Tb_trial(alpha,n);
            A(i,j)=A(i,j)+int_value;

        end
    end
end

end