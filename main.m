close all, clear all
rng(1);

n = 10;
m = 25;
sigma = 0; % shift
k = 5; % size of Krylov Subspace

A = randn(m,m);
B = randn(m,m);
[A,B] = make_pencil_square(A,B);

b = randn(max(m,n),1);
[L,U,P] = lu(A-sigma*B);
Op = @(v) U\(L\(P*B*v));

rest = 15;
computed_eigval = zeros(rest,k);
figure
for r = 1:rest
    [Q,H] = restarted_arnoldi(Op,b,k,2*k,r);
    
    eigval = sort(eig(H(1:k,1:k)));
    computed_eigval(r,:) = (1./eigval) + sigma;
    str = sprintf("Restart %d",r);

    if r == rest
        scatter(real(computed_eigval(r,:)),imag(computed_eigval(r,:)),'LineWidth',1,'DisplayName',str); hold on;
    end
end

exact_eigval = eig(A,B);

scatter(real(exact_eigval),imag(exact_eigval),'LineWidth',1,'DisplayName',"Exact"); hold on;
scatter(real(sigma),imag(sigma),'kx','LineWidth',1.5,'DisplayName',"Shift");
legend; xlim([-2 3]); ylim([-1 1]);
axis on; grid on;
hold off