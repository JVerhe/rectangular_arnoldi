function [A,B] = make_pencil_square(A,B)
    [n,m] = size(A);
    if n == m
        return
    elseif n > m
        A = [A, zeros(n, n-m)];
        B = [B, zeros(n, n-m)];
    else
        A = [A; zeros(m-n, m)];
        B = [B; zeros(m-n, m)];
    end
end