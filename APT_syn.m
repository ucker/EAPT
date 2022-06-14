function [L, S] = APT_syn(UU, D, r, para)
    [m, n, q] = size(D);
    norm_of_D = norm(D(:), 2);
    addpath tensor_toolbox_2.6;
    %% Default/Inputed parameters
    max_iter  = 200;
    tol       = 1e-5;
    if isfield(para, 'mu')
        mu = para.mu;
    else
        mu = 5;
    end
    beta = 2 * mu * r / (m+n);
    beta_init = 4*beta;
    gamma     = 0.85;    
    
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
    error_L = -1 * ones(max_iter+1, 1);
    norm_of_L = norm(para.L(:), 2);
    %% Initilization
    Unitary_fft = unitary_diag(UU, D);
    zeta = beta_init * max(abs(Unitary_fft(:)));
    
    S = wthresh(D ,'h',zeta);
    [U, Sigma, V] = ttsvd(UU, D - S, 'econ');
    Unitary_fft = unitary_diag(UU, D - S);
    L = ttprod(UU, ttprod(UU, U(:,1:r,:), Sigma(1:r, 1:r,:)), ttran(UU, V(:,1:r,:)));

    error_L(1) = 1;
    timer(1) = toc;
    zeta = beta * max(abs(Unitary_fft(:))); 
    S = wthresh(D - L, 'h', zeta);

    for i = 1 : max_iter
        %% Update L
        LL = D - S;
        LL_t = unitary_trans(UU, LL);
        L_t = zeros(size(LL));
        min_mn = min(m, n);
        Sigma = zeros(min_mn, min_mn, q);
        for s = 1 : q
            [U, Sigma(:,:,s), V] = svd(LL_t(:,:,s), 'econ');
            L_t(:,:,s) = U(:, 1:r) * Sigma(1:r, 1:r, s) * V(:, 1:r)';
        end
        sigma_max = 0;
        for s = 1 : q
            sigma_max = max([diag(Sigma(1:r, 1:r, s)); sigma_max]);
        end
        L = inv_unitary_trans(UU, L_t);
        timer(i+1) = toc;
        error_L(i+1) = norm(L(:)-para.L(:), 2) / norm_of_L;
        %% Update S
        zeta = beta * (gamma^(i)) * sigma_max;
        S = wthresh(D - L, 'h', zeta);

        %% Stop Condition
        err(i) = norm(D(:) - L(:) - S(:), 2)/norm_of_D;
        if err(i) < tol
            fprintf('Total %d iteration, final error: %e, total time without init: %f , with init: %f\n======================================\n', i, err(i), sum(timer(timer>0)),sum(timer(timer>0)));
            break;
        else
            fprintf('Iteration %d: error: %e, timer: %f \n', i, err(i), timer(i));
        end
    end
    if i == max_iter
        fprintf('Maximum iterations reached, final error: %e.\n======================================\n', err(i));
    end
    end