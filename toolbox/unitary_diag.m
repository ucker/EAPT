function [S] = unitary_diag(UU, X)

    [n1,n2,n3] = size(X);
    X = unitary_trans(UU, X);
    min_n = min(n1, n2);
    S = zeros(min_n, n3);
    
    for i = 1 : n3
        S(:, i) = svds(X(:,:,i), min_n);
    end
end