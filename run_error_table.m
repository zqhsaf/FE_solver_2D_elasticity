function run_error_table
% 复现讲义中的两张收敛性表格（Example 1，lambda=1、mu=2）
%   表 1：一次三角元（basis_type=201）
%   表 2：二次三角元（basis_type=202）
% 网格步长 h = 1/8、1/16、1/32、1/64（即 N1=N2=8,16,32,64），区域 [0,1]×[0,1]
% [修改17] 新增：用于统一输出 ||u-uh||_inf、||u-uh||_0、|u-uh|_1 及收敛阶

clc;
N_all  = [8 16 32 64];
basis_all = [201 202];
names  = {'一次三角元 (basis_type = 201)','二次三角元 (basis_type = 202)'};

for ib = 1:2

    basis = basis_all(ib);
    fprintf('\n===== %s =====\n', names{ib});
    fprintf('%8s  %14s  %14s  %14s\n','h','||u-uh||_inf','||u-uh||_0','|u-uh|_1');

    E = zeros(numel(N_all),3);
    for k = 1:numel(N_all)
        N  = N_all(k);
        t0 = tic;
        [sol,P,T,Pb,Tb] = FE_solver_2D_elasticity(0,1,0,1,basis,basis,N,N);
        [L2,H1,Linf]    = compute_error_2D(sol,P,T,Pb,Tb,basis);
        E(k,:) = [Linf, L2, H1];
        fprintf('%8s  %14.4e  %14.4e  %14.4e   (%.1f s)\n', ...
            sprintf('1/%d',N), E(k,1), E(k,2), E(k,3), toc(t0));
    end

    fprintf('收敛阶 log2(e_{k-1}/e_k)：\n');
    for k = 2:numel(N_all)
        fprintf('  h=1/%-3d -> 1/%-3d :  %.2f   %.2f   %.2f\n', ...
            N_all(k-1), N_all(k), ...
            log2(E(k-1,1)/E(k,1)), log2(E(k-1,2)/E(k,2)), log2(E(k-1,3)/E(k,3)));
    end

end

end
