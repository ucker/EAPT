function [L2_err, cost_time] = demo_APT_synthetic(p, c, n, q, r)
    addpath('tensor_toolbox_2.6');
    addpath('toolbox')
    fprintf('Starting...\n');
    switch nargin
    case 4
        r = 5;
    case 3
        q = 5;
        r = 5;
    case 2
        n = 800;
        q = 5;
        r = 5;
    case 1
        c = 0;
        n = 800;
        q = 5;
        r= 5;
    case 0
        p = 0.65;
        c = 0;
        n = 800;
        q = 5;
        r = 5;
    end
    UU = dctmtx(q);
    %% Generate a Tensor RPCA problem
    m = n;
    A = randn(m, r);
    B = randn(r, n*q);
    L_unfold = A * B;
    L_true = reshape(L_unfold, [m, n, q]);

    S_supp_idx = randsample(m*n*q, round(p*m*n*q), false);
    S_range = 3.*mean(abs(L_true(:)));
    S_temp = 2*S_range*rand(m, n, q) - S_range; 
    S_true = zeros(m, n, q);
    S_true(S_supp_idx) = S_temp(S_supp_idx);                       

    D = L_true + S_true;
    Dn = D + randn(size(D)) * c;
    para.mu        = 1.1*get_mu_tensor(UU, L_true, r);
    para.tol       = 1e-4;
    para.gamma     = 0.7;
    para.max_iter  = 100;
    para.p = p;
    para.c = c;
    para.n = n;
    para.q = q;
    para.r = r;
    para.L = L_true;
    para.trimming = true;
    tic;
    [L2, ~] = APT_syn(UU, Dn, r, para);
    cost_time = toc;
    L2_err = norm(L2(:)-L_true(:), 2)/norm(L_true(:), 2);
    fprintf('cost time is %.2f, error is %f\n', cost_time, L2_err);
end
