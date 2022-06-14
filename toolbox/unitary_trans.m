function YY = unitary_trans(UU, X)
    O = tenmat(X, [3]);
    SS = O.data;
    Y = UU'*SS;
    Y = tensor(tenmat(Y, O.rdims, O.cdims, O.tsize));
    YY = Y.data;
end