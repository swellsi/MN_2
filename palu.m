% PA=LU factorization (partial pivoting)
function [P,L,U]=palu(A)
[m,n]=size(A);

if m~=n; error('not square matrix'); end
U=A; L=eye(n); P=eye(n);

for k=1:n-1
    [foo,imax]=max(abs(U(k:end,k)));    %Find max pivot
    imax=imax+k-1;
    i1=[k,imax]; i2=[imax,k];          %exchange row index(es)
    U(i1,:)=U(i2,:); P(i1,:)=P(i2,:);  %Row exch-U / Update Perm-P
    L(i1,1:k-1)=L(i2,1:k-1);          %LowTriang Row exchange of L
    for j=[k+1:n]                     %Row index for next multipliers
        L(j,k)=U(j,k)/U(k,k);         %j-th Multiplier
        U(j,k:n)=U(j,k:n)-L(j,k)*U(k,k:n);  %Updating U
    end
end