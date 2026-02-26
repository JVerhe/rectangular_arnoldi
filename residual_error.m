function res = residual_error(A,B,pencil_size,V,D)
    n = pencil_size(2);
    A = A(:,1:n);
    B = B(:,1:n);
    res = NaN(size(V,2),1);
    for i = 1:size(V,2)
        lam = D(i);
        x = V(1:n,i);
        o = V(n+1:end,i);

        res(i) = (norm(A*x-lam*B*x))/ ((norm(A,"fro") + abs(lam)*norm(B,"fro")) * norm(x)) + norm(o);
    end
end