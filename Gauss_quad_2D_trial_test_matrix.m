function int_value=Gauss_quad_2D_trial_test_matrix(coe_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_trial,basis_index_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_index_test,basis_der_x_test,basis_der_y_test)
% 在单个三角单元上对 系数函数×试探基×检验基 做高斯数值积分
% Gauss_nodes 为 Ng×2（每行一个高斯点），Gauss_weights 为 Ng×1

int_value=0;
Ng=length(Gauss_weights);

for k=1:Ng
    x=Gauss_nodes(k,1);
    y=Gauss_nodes(k,2);
    int_value=int_value+Gauss_weights(k) ...
        *feval(coe_fun,x,y) ...
        *FE_basis_local_fun_2D(x,y,vertices,basis_type_trial,basis_index_trial,basis_der_x_trial,basis_der_y_trial) ...
        *FE_basis_local_fun_2D(x,y,vertices,basis_type_test,basis_index_test,basis_der_x_test,basis_der_y_test);
end

end