function [ U_out, V_out ] = ttrim(UU, U, Sig, V, mu_U, mu_V, s_r)
    % [m, r, q] = size(U);
    % U_mat = reshape(U, [m, r*q]);
    % row_norm_square_U = sum(U_mat.^2,2);
    % big_rows_U = row_norm_square_U > (mu_U*s_r/m/q);
    % U_mat(big_rows_U,:) = bsxfun(@times,U_mat(big_rows_U,:),(sqrt(mu_U*s_r/m/q) ./ sqrt(row_norm_square_U(big_rows_U))));   %This is A in paper
    % U = reshape(U_mat, [m ,r, q]);

    % [n, r, q] = size(V);
    % V_mat = reshape(V, [n, r*q]);
    % row_norm_square_V = sum(V_mat.^2,2);
    % big_rows_V = row_norm_square_V > (mu_V*s_r/n/q);
    % V_mat(big_rows_V,:) = bsxfun(@times,V_mat(big_rows_V,:),(sqrt(mu_V*s_r/n/q) ./ sqrt(row_norm_square_V(big_rows_V))));   %This is B in paper
    % V = reshape(V_mat, [n, r, q]);

    [m, r, q] = size(U);
    for i = 1:m
        U_i_F = norm(reshape(U(i,:,:), [], 1), 2);
        coef = min(1, sqrt(mu_U*s_r/m/q)/U_i_F);
        U(i,:,:) = coef .* U(i,:,:);
    end

    [n, r, q] = size(V);
    for j = 1:n
        V_j_F = norm(reshape(V(j,:,:), [], 1), 2);
        coef = min(1, sqrt(mu_V*r/n)/V_j_F);
        V(j,:,:) = coef .* V(j,:,:);
    end



    % These 2 QR can be computed parallelly
    U_tran = unitary_trans(UU, U);
    V_tran = unitary_trans(UU, V);

    Q1 = zeros(m, r, q);
    R1 = zeros(r, r, q);
    Q2 = zeros(n, r, q);
    R2 = zeros(r, r, q);
    for i = 1:q
        [Q1(:,:,i), R1(:,:,i)] = qr(U_tran(:,:,i), 0);
        [Q2(:,:,i), R2(:,:,i)] = qr(V_tran(:,:,i), 0);
    end
    U_temp = zeros(r, r, q);
    Sig_out = zeros(r, r, q);
    V_temp = zeros(r, r, q);
    for i = 1:q
        [U_temp(:,:,i), Sig_out(:,:,i), V_temp(:,:,i)] = svd(R1(:,:,i)*Sig(:,:,i)*R2(:,:,i)', 'econ');
    end

    U_out = zeros(m, r, q);
    V_out = zeros(n, r, q);
    for i = 1:q
        U_out(:,:,i) = Q1(:,:,i) * U_temp(:,:,i);
        V_out(:,:,i) = Q2(:,:,i) * V_temp(:,:,i);
    end
    U_out = inv_unitary_trans(UU, U_out);
    V_out = inv_unitary_trans(UU, V_out);
    end