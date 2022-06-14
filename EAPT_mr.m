function [L, S] = EAPT_mr(UU, D, r_list, r, para)
    [m, n, q] = size(D);
    norm_of_D = norm(D(:), 2);
    addpath tensor_toolbox_2.6;
    %% Default/Inputed parameters
    max_iter  = 100;
    tol       = 1e-5;
    if isfield(para, 'mu')
        mu = para.mu;
    else
        mu = 5;
    end
    beta = mu * r / (2*q*(m+n));
    beta_init = 4*beta;
    gamma     = 0.65;
    
    if isfield(para,'beta_init')
        beta_init = para.beta_init; 
        fprintf('beta_init = %f set.\n', beta_init);
    else
        fprintf('using default beta_init = %f.\n', beta_init);
    end
    if isfield(para,'beta') 
        beta = para.beta; 
        fprintf('beta = %f set.\n', beta);
    else
        fprintf('using default beta = %f.\n', beta);
    end
    if isfield(para,'gamma') 
        gamma = para.gamma; 
        fprintf('gamma = %f set.\n', tol);
    else
        fprintf('using default gamma = %f.\n', gamma);
    end
    if isfield(para,'mu') 
        mu = para.mu; 
        fprintf('mu = [%f,%f] set.\n', mu(1), mu(end));
    else
        fprintf('using default mu = [%f,%f].\n', mu, mu);
    end
    if isfield(para,'trimming') 
        trimming = para.trimming; 
        if  trimming
            fprintf('trimming = true set.\n');
        else 
            fprintf('trimming = false set.\n');
        end
    else
        fprintf('using default trimming = False.\n');
    end
    if isfield(para,'max_iter')   
        max_iter = para.max_iter; 
        fprintf('max_iter = %d set.\n', max_iter);
    else
        fprintf('using default max_iter = %d.\n', max_iter);
    end
    if isfield(para,'tol')        
        tol= para.tol; 
        fprintf('tol = %e set.\n', tol);
    else
        fprintf('using default tol = %e.\n', tol);
    end 
    
    err    = -1*ones(max_iter,1);
    timer  = -1*ones(max_iter+1,1);
    psnr_ls = -1 * ones(max_iter+1, 1);
    maxP = para.maxP;

    %% Initilization
    Unitary_fft = unitary_diag(UU, D);
    zeta = beta_init * max(abs(Unitary_fft(:)));
    
    S = wthresh(D ,'h',zeta);
    [U, Sigma, V] = ttsvd(UU, D - S, 'econ');
    Unitary_fft = unitary_diag(UU, D - S);
    L = ttprod(UU, ttprod(UU, U(:,1:r,:), Sigma(1:r, 1:r,:)), ttran(UU, V(:,1:r,:)));
    
    timer(1) = toc;
    L_hat = permute(L, [2, 3 ,1]);
    L_hat = max(L_hat, 0);
    L_hat = min(L_hat, maxP);
    psnr_ls(1) = psnr(L_hat, para.X, maxP);
    
    zeta = beta * max(abs(Unitary_fft(:))); 
    S = wthresh(D - L, 'h', zeta);
    
    init_timer = 0;
    init_err = norm(D(:)-L(:)-S(:),2)/norm_of_D;
    fprintf('Init: error: %e, timer: %f \n', init_err, init_timer);
    
    U = unitary_trans(UU, U(:,1:r,:));
    V = unitary_trans(UU, V(:,1:r,:));
    for i = 1 : max_iter
        Z = D - S;
        Z_unitary = unitary_trans(UU, Z);
        sigma_list = zeros(2);
        
        Sigma = zeros([2*r, 2*r, q]);
        for s = 1 : q
            [Q1, R1] = qr(Z_unitary(:,:,s)'*U(:,:,s)-V(:,:,s)*((Z_unitary(:,:,s)*V(:,:,s))'*U(:,:,s)), 0);
            [Q2, R2] = qr(Z_unitary(:,:,s)*V(:,:,s)-U(:,:,s)*(U(:,:,s)'*Z_unitary(:,:,s)*V(:,:,s)), 0);
            
            M = [ U(:,:,s)'*Z_unitary(:,:,s)*V(:,:,s), R1';
                  R2                             , zeros(size(R2)); ];
            [U_of_M, Sigma(:,:,s), V_of_M] = svd(M,'econ');
            U(:,:,s) = [U(:,:,s), Q2] * U_of_M(:,1:r);
            V(:,:,s) = [V(:,:,s), Q1] * V_of_M(:,1:r);
            L(:,:,s) = U(:,1:r_list(s),s) * Sigma(1:r_list(s),1:r_list(s),s) * V(:,1:r_list(s),s)'; 
        end
        for s = 1 : q
            sigma_list(1) = max([diag(Sigma(r_list(s)+1:end, r_list(s)+1:end, s)); sigma_list(1)]);
            sigma_list(2) = max([diag(Sigma(1:r_list(s), 1:r_list(s), s)); sigma_list(2)]);
        end
        L = inv_unitary_trans(UU, L);

        timer(i+1) = toc;
        L_hat = permute(L, [2, 3 ,1]);
        L_hat = max(L_hat, 0);
        L_hat = min(L_hat, maxP);
        psnr_ls(i+1) = psnr(L_hat, para.X, maxP);
    
        %% Update S
        zeta = beta * (sigma_list(1) + (gamma^i) * sigma_list(2));
        S = wthresh(D - L, 'h', zeta);
        
        %% Stop Condition
        err(i) = norm(D(:) - L(:) - S(:), 2)/norm_of_D;
        if err(i) < tol  
            fprintf('Total %d iteration, final error: %e, total time without init: %f , with init: %f\n======================================\n', i, err(i), sum(timer(timer>0)),sum(timer(timer>0))+init_timer);
            break;
        else
            fprintf('Iteration %d: error: %e, timer: %f \n', i, err(i), timer(i));
        end
    end
    if i == max_iter
        fprintf('Maximum iterations reached, final error: %e.\n======================================\n', err(i));
    end
    end