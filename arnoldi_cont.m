function [Q,H] = arnoldi_cont(Op,max_dim,Q,H)
    % Given an Arnoldi factorization of the form A*Qm = Qm+1H_
    % continue the algorithm until a dimension of max_dim is reached
    m = size(Q,2);

    for i = m:max_dim
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

        if beta < 1e-15
            fprintf("Breakdown of the Arnoldi-cont method occured at iteration %d\n", i-1)
            Q = Q(:,1:i);
            H = H(1:i,1:i);
            return
        end

        H(i+1,i) = beta;
        Q(:,i+1) = w / beta;
    end
end