%QR-factorization
%Input: matrix A nxm (n>=m)
%Output(1): U -- array of vectors uk to generate reflectors
%Output(2): R -- upper triangular factor of A=Q*R
%Note: Q is not explicitly built

function [U,R,c,x]=QRsolve(A,b)
[n,m]=size(A) ; U=[];
for k=1:m
    x=A(k:n,k); tau=sign(x(1))*norm(x); gamma=1+x(1)/tau;
    uk=[1;x(2:end)/(tau+x(1))]; U=[U, [zeros(k-1,1);uk]];
    for j=k:m
        A(k:n,j)=A(k:n,j)-gamma*uk*(uk'*A(k:n,j));             
    end
    b(k:n)=b(k:n)-gamma*uk*(uk'*b(k:n));
end
c=b;R=A;
if n==m
    x=BS(R,c);
else
    x=BS(R(1:m,1:m),c);
end