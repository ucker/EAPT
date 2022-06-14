function mu = get_mu_tensor(UU, A, r)
    [U, S, V] = ttsvd(UU, A, 'skinny');
    [m, ~, ~] = size(U);
    mu = -1;
    for i = 1:m
        mu = max(mu, norm(vec(U(i,:,:)), 2)^2 * m / r);
    end
    [n,~,~] = size(V);
    for i = 1:n 
        mu = max(mu, norm(vec(V(i,:,:)), 2)^2 * n / r);
    end
end