function res = residual_error(A,B,V,D)
    m = size(A,1); n = size(A,2);

    if m ~= n
        A = A(:,1:n);
        B = B(:,1:n);
        res = NaN(size(V,2),1);
        for i = 1:size(V,2)
            lam = D(i);
            x = V(1:n,i);
            o = V(n+1:end,i);
    
            res(i) = (norm(A*x-lam*B*x))/ ((norm(A,"fro") + abs(lam)*norm(B,"fro")) * norm(x)) + norm(o);
        end

    else      
        A = A(1:m,1:m);
        B = B(1:m,1:m);
        res = NaN(size(V,2),1);
        for i = 1:size(V,2)
            lam = D(i);
            x = V(1:n,i);
            o = V(m+1:end,i);
    
            res(i) = (norm(A*x-lam*B*x))/ ((norm(A,"fro") + abs(lam)*norm(B,"fro")) * norm(x)) + norm(o);
        end

    end
end