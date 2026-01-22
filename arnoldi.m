function [Q,H] = arnoldi(Op,b,k)
    % Creates an Arnoldi factorization A*Qk = Qk+1 * H_
    % Uses double GS orthogonalization

    n = length(b);
    Q = zeros(n,k+1);
    H = zeros(k+1,k);
    Q(:,1) = b / norm(b);

    for i = 1:k
        w = Op(Q(:,i));

        for j = 1:i
            h = Q(:,j)'*w;
            w = w - Q(:,j)*h;
            H(j,i) = h;
        end

        for j = 1:i
            g = Q(:,j)'*w;
            w = w - Q(:,j)*g;
            H(j,i) = H(j,i) + g;
        end
        
        beta = norm(w);

        if beta < 1e-12
            fprintf("Breakdown of the Arnoldi method occured at iteration %d\n", i-1)
            Q = Q(:,1:i);
            H = H(1:i,1:i);
            return
        end

        H(i+1,i) = beta;
        Q(:,i+1) = w / beta;
    end
end