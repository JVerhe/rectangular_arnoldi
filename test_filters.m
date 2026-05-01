clear all
rng(1); warning('off','all')


% sizes = [64, 128, 256, 384, 512, 768, 1024];
% n_approximations = [1, 2, 4, 8, 16, 32];
% true_eig_frac = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95];

sizes = [512, 1024];
n_approximations = [4, 8, 16, 32];
true_eig_frac = [0.05, 0.25, 0.75, 0.95];
runs = 1e2;
rest = 12; % number of restarts

LRE_data = zeros(4,rest,runs);
RRES_data = zeros(4,rest,runs);
% LRES_data = zeros(4,rest,runs);

tic
for run = 1:runs
    
    for sz_idx = 1:length(sizes)
        % m = sizes(sz_idx); n = sizes(sz_idx);
        m = 90; n = 90;
    
    
    for frac_idx = 1:length(true_eig_frac)
        % frac_true = true_eig_frac(frac_idx);
        frac_true = 0.1;
    
    
    for n_idx = 1:length(n_approximations)
        % min_dim = n_approximations(n_idx);
        min_dim = 4;
        
    end
    end
    end
    
    mess = sprintf("Now computing: %d x %d, frac %.2f, n %d", m,n,frac_true,min_dim);
    disp(mess)
    
    sigma = rand; % shift
    pencil_size = [m,n];
    max_dim = 2*min_dim; % restart threshold
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
        
        k = floor(frac_true*m); % rank of singular pencil
    
        Ap = zeros(m,m); Bp = zeros(m,m); 
        Ap(1:k,1:k) = diag(randn(k,1));
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
        
        k = floor(frac_true*m); % rank of singular pencil
    
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
    
    for fm = 1:1:4
        for r = 1:rest
            [V,W,H,K] = twosided_iram(A,B,pencil_size,sigma,v1,w1,min_dim,max_dim,r,fm);
            % Check for breakdown
            mdim = size(H,2);
            mdim2 = size(K,2);
            
            [Rev, ritz_mat] = eig(H(1:mdim,1:mdim));
            [Lev, leftritz_mat] = eig(K(1:mdim2,1:mdim2));
        
            theta = diag(ritz_mat);
            left_theta = diag(leftritz_mat);
            approx_eig = (1 ./ theta) + sigma;
            l_approx_eig = (1 ./ left_theta) + sigma;
        
            X = V(:,1:mdim) * Rev;
            Y = W(:,1:mdim2) * Lev;
        
            LRE = log_relative_error(approx_eig,exact_eig);
            LRE(isinf(LRE)) = 16; % Machine precision if infinite precision is reached
            LRE_data(fm,r,run) = mean(LRE(~isnan(LRE)));
            RRES = residual_error(A,B,pencil_size,X,approx_eig);
            RRES_data(fm,r,run) = mean(RRES(~isnan(RRES)));
            % LRES = residual_error(A',B',pencil_size,Y,l_approx_eig);
            % LRES_data(fm,r,run) = mean(LRES(~isnan(LRES)));
        end
    end

end
toc

close all

median_LRE = NaN(4,rest);
% median_LRES = NaN(4,rest);
median_RRES = NaN(4,rest);

for fm = 1:4
    for r = 1:rest
        median_LRE(fm,r) = median(LRE_data(fm,r,:),"omitmissing");
        % median_LRES(fm,r) = median(LRES_data(fm,r,:),"omitmissing");
        median_RRES(fm,r) = median(RRES_data(fm,r,:),"omitmissing");
    end
end

txt_num_eigs = sprintf("Fraction of true eigenvalues: %.2f", frac_true);
f2=figure;
for i = 1:4
    plot(median_LRE(i,:),'LineWidth',1)
    hold on;
end
xlim([1 rest]); ylim([-0.2 17]);
legend("Radius filterings", "Radius-based + inf. shift", "Residual filtering", "Residual filtering + inf. shift", 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("$-\log_{10}[(\lambda-\tilde{\lambda})/\lambda]$",'Interpreter','latex');
str = sprintf("Mean Log Relative Error for computed eigenvalues (n=%d) of pencil size %dx%d", min_dim, m, n);
title(str);
annotation("textbox", [0.1397 0.1217 0.3198 0.05201], "String", txt_num_eigs, "FontSize", 13)
axis on; grid on;
hold off;
fname = sprintf("meanLRE_size%d_approx%d_frac%.2f.pdf", m, min_dim,frac_true);
% exportgraphics(f2,fname);

f3=figure;
for i = 1:4
    semilogy(median_RRES(i,:),'LineWidth',1)
    hold on;
end
xlim([1 rest]);
legend("Radius filterings", "Radius-based + inf. shift", "Residual filtering", "Residual filtering + inf. shift", 'Interpreter', 'latex')
xlabel("Restart number"); ylabel("Residual");
str = sprintf("Mean Residual Error for (n=%d) Right computed eigenvalues of pencil size %dx%d", min_dim, m, n);
title(str);
annotation("textbox", [0.1397 0.1217 0.3198 0.05201], "String", txt_num_eigs, "FontSize", 13)
axis on; grid on;
hold off;
fname = sprintf("meanRRES_size%d_approx%d_frac%.2f.pdf", m, min_dim,frac_true);
% exportgraphics(f3,fname);

% f4=figure;
% for i = 1:4
%     semilogy(median_LRES(i,:),'LineWidth',1)
%     hold on;
% end
% xlim([1 rest]);
% legend("Radius filterings", "Radius-based + inf. shift", "Residual filtering", "Residual filtering + inf. shift", 'Interpreter', 'latex')
% xlabel("Restart number"); ylabel("Residual");
% str = sprintf("Mean Residual Error for (n=%d) Left computed eigenvalues of pencil size %dx%d", min_dim, m, n);
% title(str);
% annotation("textbox", [0.1397 0.1217 0.3198 0.05201], "String", txt_num_eigs, "FontSize", 13)
% axis on; grid on;
% hold off;
% fname = sprintf("meanLRES_size%d_approx%d_frac%.2f.pdf", m, min_dim,frac_true);
% exportgraphics(f4,fname);
