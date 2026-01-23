close all, clear all
rng(1);

n = 120;
m = 30;
num_eigs = floor(0.4*m); % for the singular pencil
sigma = -0.5; % shift
rest = 6; % number of restarts
min_dim = 5; % #approx. eigenvalues
max_dim = 2*min_dim; % restart threshold
pencil_type = 2; % 0: regular, 1: singular, 2: singular bordered

if pencil_type == 0 %Regular pencil
    A = randn(m,m);
    B = randn(m,m);
    [A,B] = make_pencil_square(A,B);

    b = randn(m,1);
    [L,U,P] = lu(A-sigma*B);
    Op = @(v) U\(L\(P*B*v));

    exact_eigval = eig(A,B);
    qz_eig = exact_eigval;

elseif pencil_type == 1 % Singular pencil with orthogonal transformation

    A = zeros(n,m);
    A(1:num_eigs,1:num_eigs) = diag(2*rand(num_eigs,1)-1);

    B = zeros(n,m);
    B(1:num_eigs,1:num_eigs) = eye(num_eigs);
    exact_eigval = diag(A)./diag(B);

    [Q1,~] = qr(randn(m,m));
    [~,R2] = qr(randn(n,m));
    A = R2'*A*Q1; B = R2'*B*Q1;

    b = randn(m,1);
    [L,U,P] = lu(A-sigma*B);
    Op = @(v) U\(L\(P*B*v));
    
    qz_eig = eig(A,B);


elseif pencil_type == 2 % Singular bordered pencil
    
    
    A = zeros(n,m);
    A(1:num_eigs,1:num_eigs) = diag(2*rand(num_eigs,1)-1);

    B = zeros(n,m);
    B(1:num_eigs,1:num_eigs) = eye(num_eigs);
    exact_eigval = diag(A)./diag(B);

    % Transform so the problem is not trivial
    [Q,~] = qr(randn(n,n));
    A = Q'*A; 
    B = Q'*B;

    th = 1e-12;
    S = A-sigma*B;
    alpha = norm(S(:),inf);
    [L,U,P,C,D,~,~]=LU_border_partial(S,th,alpha);
    
    Ap = [A , D ; C', zeros(size(C,2),size(D,2))] ;
    Bp = [B , zeros(size(D,1),size(D,2)) ; zeros(size(C,2),size(C,1)), zeros(size(C,2),size(D,2))] ;

    Op = @(x) U\(L\(P*x)) ;
    b = randn(size(Ap,2),1);

    [Aqz, Bqz] = make_pencil_square(A,B);
    qz_eig = eig(Aqz,Bqz);

end

% Compute eigenvalues and plot
figure

[Q,H] = restarted_arnoldi(Op,b,min_dim,max_dim,rest);
ritzval = eig(H(1:min_dim,1:min_dim));
iram_eig = (1./ritzval) + sigma;

iram_err = eigapprox_error(exact_eigval,iram_eig,sigma)
qz_err = eigapprox_error(exact_eigval,qz_eig,sigma)

str = sprintf("iram %d restarts",rest);

scatter(real(iram_eig),imag(iram_eig),'b*','LineWidth',1,'DisplayName',str); hold on;
qz_eig = qz_eig(~isinf(qz_eig));
% scatter(real(qz_eig),imag(qz_eig),'g+','LineWidth',1,'DisplayName','qz'); hold on;
scatter(real(exact_eigval),imag(exact_eigval),'LineWidth',1,'DisplayName',"Exact"); hold on;
scatter(real(sigma),imag(sigma),'kx','LineWidth',1.5,'DisplayName',"Shift");
legend; xlim([-2 3]); ylim([-1 1]); xlabel("Real axis"); ylabel("Imag axis");
title("Eigenvalue approximation with Iram")
axis on; grid on;
hold off