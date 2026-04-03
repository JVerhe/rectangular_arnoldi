close all, clear all
rng(1);

m = 90;
n = 90;
pencil_size = [m,n];
sigma = -1; % shift
rest = 6; % number of restarts
min_dim = 4; % #approx. eigenvalues
max_dim = 2*min_dim; % restart threshold
filt_method = 2;    % 0: radius filtering, 
                    % 1: radius filtering with inf eigenvalue filter
                    % 2: residue filtering
                    % 3: residue filtering with inf eigenvalue filter
pencil_type = 1;

if pencil_type == 0 % Regular square pencil
    
    k = min(m,n);

    A = diag(diag(randn(m,n)));
    B = eye(m,n);
    exact_eig = diag(A);
    
    [Q1,~] = qr(randn(m,m));
    [Q2,~] = qr(randn(m,m));
    A = Q2'*A*Q1; B = Q2'*B*Q1;

    S = A-sigma*B;
    
    v1 = randn(n,1); w1 = randn(m,1);
    [L1,U1,P1] = lu(A-sigma*B);
    Op1 = @(v) U1\(L1\(P1*B*v));
    [L2,U2,P2] = lu((A-sigma*B)');
    Op2 = @(v) U2\(L2\(P2*B'*v));

elseif pencil_type == 1 % (2m-r)x(2m-r) regularized square pencil
    
    k = floor(0.8*m); % rank of singular pencil

    Ap = zeros(m,m); Bp = zeros(m,m); 
    Ap(1:k,1:k) = diag(diag(randn(k)));
    Bp(1:k,1:k) = eye(k);
    exact_eig = diag(Ap)./diag(Bp);

    [Q1,~] = qr(randn(m,m));
    [Q2,~] = qr(randn(m,m));
    Ap = Q2'*Ap*Q1; Bp = Q2'*Bp*Q1;   
    
    A = [Ap, randn(m, m-k); randn(m-k, m), randn(m-k,m-k)];
    B = [Bp, zeros(m, m-k); zeros(m-k, m), zeros(m-k,m-k)];

    S = A-sigma*B;

    v1 = randn(size(S,2),1); w1 = randn(size(S,1),1);

elseif pencil_type == 2 % mxm singular square pencil
    
    k = floor(0.8*m); % rank of singular pencil

    A = zeros(m,m); B = zeros(m,m); 
    A(1:k,1:k) = diag(diag(randn(k)));
    B(1:k,1:k) = eye(k);
    exact_eig = diag(A)./diag(B);

    [Q1,~] = qr(randn(m,m));
    [Q2,~] = qr(randn(m,m));
    A = Q2'*A*Q1; B = Q2'*B*Q1;   

    S = A-sigma*B;

    v1 = randn(size(S,2),1); w1 = randn(size(S,1),1);

end

fac_error_1 = NaN(rest,1);
fac_error_2 = NaN(rest,1);
LRE = NaN(min_dim,rest);
RES = NaN(min_dim,rest);
LAM = NaN(min_dim,rest);

for r = 1:rest
    [V,W,H,K] = twosided_iram(A,B,pencil_size,sigma,v1,w1,min_dim,max_dim,r,filt_method);

    fac_error_1(r) = norm(S \ (B * V(:,1:end-1)) - V*H);
    fac_error_2(r) = norm(S' \ (B' * W(:,1:end-1)) - W*K);
    
    [Rev, ritz_mat] = eig(H(1:min_dim, 1:min_dim));
    [Lev, leftritz_mat] = eig(K(1:min_dim, 1:min_dim));

    theta = diag(ritz_mat);
    left_theta = diag(leftritz_mat);
    approx_eig = (1 ./ theta) + sigma;
    l_approx_eig = (1 ./ left_theta) + sigma;

    X = V(:,1:end-1) * Rev;
    Y = W(:,1:end-1) * Lev;

    % Vr = V(:,1:end-1); Wr = W(:,1:end-1);
    % [R,ritz,L] = eig(Wr'*A*Vr,Wr'*B*Vr);
    % X = Vr * R;
    % Y = Wr * L;
    % theta = diag(ritz);
    % approx_eig = (1 ./ theta) + sigma;
    % l_approx_eig = approx_eig;

    [approx_eig, idx_R] = sort(approx_eig, 'ComparisonMethod', 'abs');
    X = X(:, idx_R);
    [l_approx_eig, idx_L] = sort(l_approx_eig, 'ComparisonMethod', 'abs');
    Y = Y(:, idx_L);

       
    LRE(:,r) = log_relative_error(approx_eig,exact_eig);
    RES(:,r) = residual_error(A,B,X,approx_eig);
    LAM(:,r) = approx_eig;
end

figure
semilogy(fac_error_1)
hold on
semilogy(fac_error_2)
hold off
axis on, grid on
xlabel("#Restarts");ylabel("Factorization Error");
legend("SV-VH","S'W-WK");

figure
for i = 1:min_dim
    plot(real(LAM(i,:)))
    hold on;
end
xlim([1 rest]); ylim([sigma-1.5 sigma+1.5]);
labels = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k),1:min_dim,'UniformOutput',false);
legend(labels, 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("Value of $\mathcal{R}(\lambda_i)$",'Interpreter','latex');
title("Eigenvalue Approximations");
axis on; grid on;
hold off;

figure
for i = 1:min_dim
    plot(LRE(i,:))
    hold on;
end
xlim([1 rest]); ylim([0 17]);
labels = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k),1:min_dim,'UniformOutput',false);
legend(labels, 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("$-\log_{10}[(\lambda-\tilde{\lambda})/\lambda]$",'Interpreter','latex');
title("Log Relative Error for computed eigenvalues");
axis on; grid on;
hold off;

figure
for i = 1:min_dim
    semilogy(RES(i,:))
    hold on;
end
xlim([1 rest]);
labels = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k),1:min_dim,'UniformOutput',false);
legend(labels, 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("Residual");
title("Residual Error for computed eigenvalues");
axis on; grid on;
hold off;

figure
title('Eigenvalues')
scatter(real(exact_eig),imag(exact_eig),'ro');
hold on
scatter(real(sigma),imag(sigma),'kx');
scatter(real(approx_eig),imag(approx_eig),'b*');
axis on; grid on;
xlabel("Real Axis");ylabel("Imag axis");
legend("Exact eigvals","Shift","Approx eigvals");
hold off

figure
tiledlayout(1,3)
nexttile
heatmap(real(approx_eig))
title("Re(eig)")
nexttile
heatmap(abs(Y))
title('Left eigenvectors')
nexttile
heatmap(abs(X))
title('Right eigenvectors')