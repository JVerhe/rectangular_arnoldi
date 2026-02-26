function [Q,H,Qhist,Hhist] = restarted_arnoldi(Op,b,min_dim,max_dim,restarts,pencil_size)
    % Creates an Arnoldi factorization A*Qk = Qk+1 * H_
    %
    %   Op: defines matrix vector product
    %   b: starting vector
    %   min_dim: final dimension of Krylov subspace
    %   max_dim: maximal dimension of Krylov subspace before restart
    %   restarts: number of restarts
    m = pencil_size(1); n = pencil_size(2);
    Qhist = zeros(m,min_dim+1,restarts);
    Hhist = zeros(min_dim+1,min_dim,restarts);

    [Q,H] = arnoldi(Op,b,min_dim);

    for r = 1:restarts
        
        [Q,H] = arnoldi_cont(Op,max_dim,Q,H);

        Qm = Q(:,1:max_dim);
        qres = Q(:,end);
        Hm = H(1:max_dim,1:max_dim);

        % % Compute ordered schur based on norm
        % [U,T] = schur(Hm);
        % e = diag(T);
        % [~,ind] = sort(-abs(e));
        % select = ones(1,max_dim);
        % select(ind(min_dim+1:end)) = 0;
        % [U,T] = ordschur(U,T,select);
        
        % Compute ordered schur based on residual for rect pencil
        [U,T] = schur(Hm);
        v = Qm*U;
        res = vecnorm(v(n+1:end,:));
        [~,ind] = sort(abs(res));
        select = ones(1,max_dim);
        select(ind(min_dim+1:end)) = 0;
        [U,T] = ordschur(U,T,select);
        
        % Truncate
        U1 = U(:,1:min_dim); 
        T1 = T(1:min_dim, 1:min_dim);
        
        % 'Restart' Matrices
        Q = [Qm*U1, qres];
        hm = U1'*H(end,:)';
        H = [T1; hm'];

        Qhist(:,:,r) = Q;
        Hhist(:,:,r) = H;

        if norm(hm) < 1e-15
            fprintf("Breakdown of restarted Arnoldi after %d restarts\n",r);
            return;
        end
    end
end