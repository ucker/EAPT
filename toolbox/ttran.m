function Xt = ttran(UU, X)
    [n1,n2,n3] = size(X);
    X = unitary_trans(UU, X);
    Xt = zeros(n2,n1,n3);
    for i = 1 : n3
        Xt(:,:,i) = X(:,:,i)';
    end
    Xt = inv_unitary_trans(UU, Xt);