function [X]=QRsolveX(A,b)
[n,m]=size(A) ; U=[];
for k=1:m
    X=A(k:n,k); tau=sign(X(1))*norm(X); gamma=1+X(1)/tau;
    uk=[1;X(2:end)/(tau+X(1))]; U=[U, [zeros(k-1,1);uk]];
    for j=k:m
        A(k:n,j)=A(k:n,j)-gamma*uk*(uk'*A(k:n,j));             
    end
    b(k:n)=b(k:n)-gamma*uk*(uk'*b(k:n));
end
c=b;R=A;
if n==m
    X=BS(R,c);
else
    X=BS(R(1:m,1:m),c);
end