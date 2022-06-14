function [L, S] = EAPT_syn(UU, D, r, para)
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
    beta = mu * r / (2 * (m+n));
    beta_init = 4*beta;
    gamma     = 0.65;    
    trimming  = false;
    
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
    
    err    = -1*ones(max_iter+1,1);
    timer  = -1*ones(max_iter+1,1);
    error_L = -1 * ones(max_iter+1, 1);
    norm_of_L = norm(para.L(:), 2);
    %% Initilization
    Unitary_fft = unitary_diag(UU, D);
    zeta = beta_init * max(abs(Unitary_fft(:)));
    
    S = wthresh(D ,'h',zeta);
    [U, Sigma, V] = ttsvd(UU, D - S, 'econ');
    Sigma_tran = unitary_trans(UU, Sigma(1:r, 1:r, :));
    Unitary_fft = unitary_diag(UU, D - S);

    L = ttprod(UU, ttprod(UU, U(:,1:r,:), Sigma(1:r, 1:r,:)), ttran(UU, V(:,1:r,:)));
    
    error_L(1) = 1;
    timer(1) = toc;
    zeta = beta * max(abs(Unitary_fft(:))); 
    S = wthresh(D - L, 'h', zeta);
    
    U = unitary_trans(UU, U(:,1:r,:));
    V = unitary_trans(UU, V(:,1:r,:));

    for i = 1 : max_iter
        U = inv_unitary_trans(UU, U(:,1:r,:));
        V = inv_unitary_trans(UU, V(:,1:r,:));
        if trimming
            [U, V] = ttrim(UU, U, Sigma_tran(1:r,1:r,:), V, mu(1), mu(end), r * q);
        end
        %% Update L
        U = unitary_trans(UU, U(:,1:r,:));
        V = unitary_trans(UU, V(:,1:r,:));
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
            L(:,:,s) = U(:,:,s) * Sigma(1:r,1:r,s) * V(:,:,s)'; 
        end
        % Sigma_ifft = ifft(Sigma, [], 3);
        Sigma_vector = vec(abs(Sigma(1:r,1:r,:)));
        Sigma_vector_extra = vec(abs(Sigma(r+1:end, r+1:end, :)));
        sorted_Sigma = sort(Sigma_vector, 'descend');
        sorted_Sigma_extra = sort(Sigma_vector_extra, 'descend');
        sigma_list(1) = sorted_Sigma_extra(1);
        sigma_list(2) = sorted_Sigma(1);
        Sigma_tran = Sigma;
        % Sigma = inv_unitary_trans(UU, Sigma);
        L = inv_unitary_trans(UU, L);
        timer(i+1) = toc;
        error_L(i+1) = norm(L(:)-para.L(:), 2) / norm_of_L;

        %% Update S
        zeta = beta * (sigma_list(1) + (gamma^i) * sigma_list(2));
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