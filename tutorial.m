close all, clear all
rng(1);

n1 = 50;
eigs1 = 0.3*randn(n1,1) + 0.3*1i*randn(n1,1);

% n2 = 5;
% eigs2 = 0.3*randn(n2,1) + 0.3*1i*randn(n2,1) + 2 + 2*1i;

n2 = 6;
eigs2 = randn(n2,1) + 1i*randn(n2,1);
r = 2;
for i = 1:n2
    eigs2(i) = r*(eigs2(i) / norm(eigs2(i)));
end


n = n1+n2;
A = blkdiag(diag(eigs1),diag(eigs2));
[Q,~] = qr(rand(n));
A = Q'*A*Q;

b = rand(n,1);
b = b / norm(b);
Op = @(v) A*v;

f0=figure;
scatter(real(eigs1),imag(eigs1),'LineWidth',1);
    hold on;
    scatter(real(eigs2),imag(eigs2),'LineWidth',1);
    xlim([-2.5 2.5]);ylim([-2.5 2.5]);
    xlabel("Real Axis");ylabel("Imag Axis");
    legend("Group 1","Group 2",'Location','nw');
    title("True Eigenvalues")
    axis on; grid on; axis square;

    exportgraphics(f0,"arnoldi_ex_0.pdf");

for k = 4:4:16
    [~,H] = arnoldi(Op,b,k);
    
    ritz = eig(H(1:k,1:k));
    
    f1 = figure;
    scatter(real(eigs1),imag(eigs1),'LineWidth',1);
    hold on;
    scatter(real(eigs2),imag(eigs2),'LineWidth',1);
    hold on;
    scatter(real(ritz),imag(ritz),'kx','linewidth',1.5);
    xlim([-2.5 2.5]);ylim([-2.5 2.5]);
    xlabel("Real Axis");ylabel("Imag Axis");
    legend("Group 1","Group 2","Approximations",'Location','nw');
    title(['Eigenvalue Approximations at Iteration k = ', num2str(k)]);
    axis on; grid on; axis square;

    str = sprintf("arnoldi_ex_%d.pdf",k/4);
    exportgraphics(f1,str);
end