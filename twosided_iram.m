function [V,W,H,K] = twosided_iram(A,v1,w1,min,max,restarts,pencil_size)
    
    m = pencil_size(1); n = pencil_size(2);
    Op1 = @(v) A*v; Op2 = @(v) A'*v;
    [V,H] = arnoldi(Op1,v1,min);
    [W,K] = arnoldi(Op2,w1,min);

    for r = 1:restarts
        
        [V,H] = arnoldi_cont(Op1,max,V,H);
        [W,K] = arnoldi_cont(Op2,max,W,K);

        ritz = eig(H(1:max,1:max));
        [~, idx] = sort(abs(ritz), 'ascend');  % Unwanted eigenvalues
        shifts = ritz(idx(1:(max - min)));
        
        vres = V(:, max+1);
        h_m = H(max+1, max);
        wres = W(:, max+1);
        k_m = K(max+1, max);
        
        % Initialize accumulators for Q
        Q_v = eye(max);
        Q_w = eye(max);
        
        for i = 1:length(shifts)
            [Qm1,~] = qr(H(1:max,1:max) - shifts(i)*eye(max));
            H(1:max,1:max) = Qm1'* H(1:max,1:max) * Qm1;
            Q_v = Q_v * Qm1;
            
            [Qm2,~] = qr(K(1:max,1:max) - shifts(i)*eye(max));
            K(1:max,1:max) = Qm2' * K(1:max,1:max) * Qm2;
            Q_w = Q_w * Qm2;
        end
        
        % Update bases
        V(:,1:max) = V(:,1:max) * Q_v;
        W(:,1:max) = W(:,1:max) * Q_w;
        
        % Calculate new residuals
        f_v = V(:, min+1) * H(min+1, min) + vres * h_m * Q_v(max, min);
        f_w = W(:, min+1) * K(min+1, min) + wres * k_m * Q_w(max, min);
        
        beta_v = norm(f_v);
        vnew = f_v / beta_v;
        
        beta_w = norm(f_w);
        wnew = f_w / beta_w;
        
        % Truncate
        V = [V(:, 1:min), vnew]; 
        W = [W(:, 1:min), wnew];
        
        H = H(1:min+1, 1:min);
        H(min+1, min) = beta_v;
        
        K = K(1:min+1, 1:min);
        K(min+1, min) = beta_w;
    end
end