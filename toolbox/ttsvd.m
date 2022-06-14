function [U, S, V] = ttsvd(UU, X, opt)
    if ~exist('opt', 'var')
        opt = 'full';
    end
    
    [n1,n2,n3] = size(X);
    X = unitary_trans(UU, X);
    if strcmp(opt,'skinny') == 1 || strcmp(opt,'econ') == 1 
        min12 = min(n1,n2);
        U = zeros(n1,min12,n3);
        S = zeros(min12,min12,n3);
        V = zeros(n2,min12,n3);
        
        for i = 1 : n3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i),'econ');
        end
        if strcmp(opt,'skinny') == 1
            s1 = diag(sum(S,3))/n3;
            tol = max(n1,n2)*eps(max(s1));
            trank = sum(s1 > tol); % tensor tubal rank
            U = U(:,1:trank,:);
            V = V(:,1:trank,:);
            S = S(1:trank,1:trank,:);        
        end
        
    elseif strcmp(opt,'full') == 1
        U = zeros(n1,n1,n3);
        S = zeros(n1,n2,n3);
        V = zeros(n2,n2,n3);
        
        for i = 1 : n3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i));
        end
    end
    
    U = inv_unitary_trans(UU, U);
    S = inv_unitary_trans(UU, S);
    V = inv_unitary_trans(UU, V);
    
    end