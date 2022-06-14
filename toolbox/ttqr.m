function [Q, R] = ttqr(UU, X, opt)
[n1,n2,n3] = size(X);
X = unitary_trans(UU, X);

if n1>n2 && exist('opt', 'var') && strcmp(opt,'econ') == 1    
    Q = zeros(n1,n2,n3);
    R = zeros(n2,n2,n3);    
    for i = 1 : n3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
    end
else    
    Q = zeros(n1,n1,n3);
    R = zeros(n1,n2,n3);
    for i = 1 : n3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i));
    end
end

Q = inv_unitary_trans(UU, Q);
R = inv_unitary_trans(UU, R);
end