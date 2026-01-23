function [Q,H] = restarted_arnoldi(Op,b,min_dim,max_dim,restarts)
    % Creates an Arnoldi factorization A*Qk = Qk+1 * H_
    %
    %   Op: defines matrix vector product
    %   b: starting vector
    %   min_dim: final dimension of Krylov subspace
    %   max_dim: maximal dimension of Krylov subspace before restart
    %   restarts: number of restarts
    [Q,H] = arnoldi(Op,b,max_dim);
    for r = 1:restarts
        
        Qm = Q(:,1:max_dim);
        qres = Q(:,end);
        Hm = H(1:max_dim,1:max_dim);
        
        % Compute ordered schur
        [U,T] = schur(Hm);
        e = diag(T) ;
        [~,ind] = sort(-abs(e)) ;
        select = ones(1,max_dim) ;
        select(ind(min_dim+1:end)) = 0 ;
        [U,T] = ordschur(U,T,select) ;
        
        % Truncate
        U1 = U(:,1:min_dim); 
        T1 = T(1:min_dim, 1:min_dim);
        
        % 'Restart' Matrices
        Q = [Qm*U1, qres];
        hm = U1'*H(end,:)';
        H = [T1; hm'];

        if norm(hm) < 1e-15
            fprintf("Breakdown of restarted Arnoldi after %d restarts\n",r);
            return;
        end
        
        [Q,H] = arnoldi_cont(Op,max_dim,Q,H);
    end
end