close all, clear all
rng(1);

m = 60;
n = 24;
pencil_size = [m,n];
sigma = -0.5; % shift
rest = 8; % number of restarts
min_dim = 6; % #approx. eigenvalues
max_dim = 2*min_dim; % restart threshold
pencil_type = 3; % 0: regular, 1: singular, 2: singular bordered

if pencil_type == 0 % Regular pencil
    A = randn(m,m);
    B = randn(m,m);
    [A,B] = make_pencil_square(A,B);

    b = randn(m,1);
    [L,U,P] = lu(A-sigma*B);
    Op = @(v) U\(L\(P*B*v));

    exact_eigval = eig(A,B);
    qz_eig = exact_eigval;

elseif pencil_type == 1 % Singular pencil with orthogonal transformation
    
    num_eigs = floor(0.4*n); % for the singular pencil
    A = zeros(m,n);
    A(1:num_eigs,1:num_eigs) = diag(2*rand(num_eigs,1)-1);

    B = zeros(m,n);
    B(1:num_eigs,1:num_eigs) = eye(num_eigs);
    exact_eigval = diag(A)./diag(B);

    [Q1,~] = qr(randn(n,n));
    [~,R2] = qr(randn(m,n));
    A = R2'*A*Q1; B = R2'*B*Q1;

    b = randn(n,1);
    [L,U,P] = lu(A-sigma*B);
    Op = @(v) U\(L\(P*B*v));
    
    qz_eig = eig(A,B);


elseif pencil_type == 2 % Singular bordered pencil
    
    num_eigs = floor(0.4*n); % for the singular pencil
    A = zeros(m,n);
    A(1:num_eigs,1:num_eigs) = diag(2*rand(num_eigs,1)-1);

    B = zeros(m,n);
    B(1:num_eigs,1:num_eigs) = eye(num_eigs);
    exact_eigval = diag(A)./diag(B);

    % Transform so the problem is not trivial
    [Q,~] = qr(randn(m,m));
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

elseif pencil_type == 3
    
    A = zeros(m,n);
    A(1:n,1:n) = diag(2*rand(n,1)-1);

    B = zeros(m,n);
    B(1:n,1:n) = eye(n);
    exact_eigval = diag(A)./diag(B);

    [Q1,~] = qr(randn(n,n));
    A = A*Q1; B = B*Q1;
    
    EA = randn(m,m-n);
    EB = zeros(m,m-n);
    A = [A, EA];
    B = [B, EB];

    b = randn(m,1);
    [L,U,P] = lu(A-sigma*B);
    Op = @(v) U\(L\(P*B*v));
    
    qz_eig = eig(A,B);
    
end


% Compute eigenvalues and plot
[Q,H,Qhist,Hhist] = restarted_arnoldi(Op,b,min_dim,max_dim,rest,pencil_size);


LRE = NaN(min_dim,rest);
RES = NaN(min_dim,rest);
LAM = NaN(min_dim,rest);
for r = 1:rest  
    [V, ritzval] = eig(Hhist(1:min_dim,1:min_dim,r));
    V = Qhist(:,1:end-1,r) * V;
    
    [iram_eig, idx] = sort((1./diag(ritzval)) + sigma);
    
    V = V(:, idx);

    LAM(:,r) = iram_eig;
    LRE(:,r) = log_relative_error(iram_eig,exact_eigval);
    RES(:,r) = residual_error(A,B,pencil_size,V,iram_eig);
end

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
for i = 1:min_dim
    plot(LAM(i,:))
    hold on;
end
xlim([1 rest]);
labels = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k),1:min_dim,'UniformOutput',false);
legend(labels, 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("Eigenvalue",'Interpreter','latex');
title("Computed eigenvalues");
axis on; grid on;
hold off;

figure
tiledlayout(1,2)
nexttile
heatmap(diag(iram_eig));
nexttile
heatmap(V);

figure
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