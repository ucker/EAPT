function C = ttprod(UU, A, B)
    [n1, n2, n3] = size(A);
    [m1, m2, m3] = size(B);
    
    if n2 ~= m1 || n3 ~= m3 
        error('Inner tensor dimensions must agree.');
    end
    
    A = unitary_trans(UU, A);
    B = unitary_trans(UU, B);
    C = zeros([n1, m2, n3]);
    
    for i = 1 : n3
        C(:, :, i) = A(:, :, i) * B(:, :, i);
    end
    C = inv_unitary_trans(UU, C);
end
    