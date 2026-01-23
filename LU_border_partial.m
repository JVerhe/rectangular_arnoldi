function [L,U,P,C,D,QL,QR]=LU_border_partial(A,threshold,alpha)
% LU decom for bordered GEP
% [L,U,P,C,D,QL,QR]=LU_border_partial(A,threshold,alpha)
% ****** P*Ap=L\U ******
% where Ap=[A,D;C',zeros(size(C,2))];
[n1,n2]=size(A);
m=min(n1,n2);
P=speye(n1,n1);
C=sparse(n2,0);
D=sparse(n1,0);
QL = sparse(n1,0);
QR = sparse(n2,0);

for i=1:n2
%    disp(sprintf("i=%d", i))
    [pivot,I]=max(abs(A(i:end,i)));
    I = I + i - 1 ;
    if abs(pivot)<threshold || i>size(A,1)
        A(size(A,1)+1,i) = alpha ; % C
        C = [C , A(end,1:n2)'] ;
%        if size(A,1)>size(A,2)
%          A(i,end+1)=alpha ;
%        end
        P=blkdiag(P,1);
        I = size(A,1); % set pivot
    end

    % Pivot rows
    tmp = A(i,:) ; A(i,:) = A(I,:) ; A(I,:) = tmp ;
    tmp = P(i,:) ; P(i,:) = P(I,:) ; P(I,:) = tmp ;

    % Update
    A(i+1:end,i) = A(i+1:end,i) / A(i,i) ;
    A(i+1:end,i+1:end) = A(i+1:end,i+1:end) - A(i+1:end,i) * A(i,i+1:end) ;

    QR(i,end+1) = 1 ;
end

% Expand for columns between n2 and n1 --> D
A(n2+1:size(A,1),n2+1:size(A,1)) = alpha*speye(size(A,1)-n2) ;
if size(A,2)~=size(A,1)
  disp('A is not square: bug')
end

L = tril(A,-1) + speye(size(A));
U = triu(A) ;

D=P(n2+1:end,1:n1)'*alpha ;

QL = P(1:size(A,2),:)' ;
