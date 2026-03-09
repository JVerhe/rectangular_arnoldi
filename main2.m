close all, clear all
rng(1);

m = 25;
n = 25;
pencil_size = [m,n];
sigma = -1; % shift
rest = 5; % number of restarts
min_dim = 2; % #approx. eigenvalues
max_dim = 2*min_dim; % restart threshold
filt_method = 1;
pencil_type = 1;

if pencil_type == 0 % Regular square pencil
    
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
    
    [V,W,H,K] = twosided_iram(Op1,Op2,v1,w1,min_dim,max_dim,rest,filt_method,pencil_size);
    
    norm(inv(S)*B*V(:,1:end-1)-V*H)
    norm(inv(S')*B'*W(:,1:end-1)-W*K)
    ritz = eig(H(1:min_dim,1:min_dim));

elseif pencil_type == 1 % mxm singular square pencil
    
    r = floor(0.8*m); % rank of singular pencil

    Ap = zeros(m,m); Bp = zeros(m,m); 
    Ap(1:r,1:r) = diag(diag(randn(r)));
    Bp(1:r,1:r) = eye(r);
    exact_eig = diag(Ap)./diag(Bp);

    [Q1,~] = qr(randn(m,m));
    [Q2,~] = qr(randn(m,m));
    Ap = Q2'*Ap*Q1; Bp = Q2'*Bp*Q1;   
    
    A = [Ap, randn(m, m-r); randn(m-r, m), randn(m-r,m-r)];
    B = [Bp, zeros(m, m-r); zeros(m-r, m), zeros(m-r,m-r)];

    S = A-sigma*B;

    [L1,U1,P1] = lu(A-sigma*B);
    Op1 = @(v) U1\(L1\(P1*B*v));
    [L2,U2,P2] = lu((A-sigma*B)');
    Op2 = @(v) U2\(L2\(P2*B'*v));

    v1 = randn(size(S,2),1); w1 = randn(size(S,1),1);

    [V,W,H,K] = twosided_iram(Op1,Op2,v1,w1,min_dim,max_dim,rest,filt_method,pencil_size);
    
    fac_error = norm(inv(S)*B*V(:,1:end-1)-V*H)
    fac_error = norm(inv(S')*B'*W(:,1:end-1)-W*K)
    [Rev, ritz] = eig(H(1:min_dim,1:min_dim));
    [Lev, leftritz] = eig(K(1:min_dim,1:min_dim));

end

approx_eig = diag((1./ritz) + sigma);
l_approx_eig = diag((1./leftritz) + sigma);
R = V(:,1:end-1) * Rev;
L = W(:,1:end-1) * Lev;

log_relative_error(approx_eig,exact_eig)
log_relative_error(l_approx_eig,exact_eig)

residual_error(A,B,pencil_size,R,approx_eig)
residual_error(A',B',pencil_size,L,l_approx_eig)

figure
title('Eigenvalues')
scatter(real(exact_eig),imag(exact_eig),'ro');
hold on
scatter(real(sigma),imag(sigma),'kx');
scatter(real(approx_eig),imag(approx_eig),'b*');
scatter(real(l_approx_eig),imag(l_approx_eig),'g*');
axis on; grid on;
xlabel("Real Axis");ylabel("Imag axis");
legend("Exact eigvals","Shift","Right approximated eigvals","Left approximated eigvals");
hold off

figure
tiledlayout(1,2)
nexttile
heatmap(abs(L))
title('Left eigenvectors')
nexttile
heatmap(abs(R))
title('Right eigenvectors')