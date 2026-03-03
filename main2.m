close all, clear all
rng(1);

m = 30;
n = 30;
pencil_size = [m,n];
sigma = -0.5; % shift
rest = 10; % number of restarts
min_dim = 3; % #approx. eigenvalues
max_dim = 2*min_dim; % restart threshold
pencil_type = -1;

if pencil_type == 0

    num_eig = floor(0.5*n);
    
    A = zeros(m,n); B = zeros(m,n);
    A(1:num_eig,1:num_eig) = diag(diag(randn(num_eig)));
    B(1:num_eig,1:num_eig) = diag(diag(randn(num_eig)));
    exact_eig = diag(A)./diag(B);

    [Q1,~] = qr(randn(n,n));
    [Q2,~] = qr(randn(m,m));
    A = Q2'*A*Q1; B = Q2'*B*Q1;

    approx_eig = eig(A,B);

elseif pencil_type == -1
    
    A = randn(10,10);
    exact_eig = eig(A);
    
    v1 = randn(10,1); w1 = randn(10,1);
    Op1 = @(v) A*v;
    Op2 = @(v) A'*v;

    [V,W,H,K] = twosided_iram(A,v1,w1,min_dim,max_dim,rest,pencil_size);
    
    % Check residual = 0?
    size(V), size(H)
    size(W), size(K)
    norm(A*V(:,1:end-1)-V*H)
    norm(A'*W(:,1:end-1)-W*K)

    approx_eig = eig(H(1:min_dim,1:min_dim));
end

figure
scatter(real(exact_eig),imag(exact_eig),'ro');
hold on
scatter(real(approx_eig),imag(approx_eig),'b*');
axis on; grid on;
xlabel("Real Axis");ylabel("Imag axis");
legend("Exact eigvals","approx");
hold off