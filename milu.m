%A=LU factorization (no pivoting)
function [L,U]=milu(A)
[m,n]=size(A);

if m~=n; error('not square matrix'); end
U=A; L=eye(n);
for k=1:n-1
    for j=[k+1:n]
        ['*******************k=' int2str(k) '*******************']
        'U before'
        U
        L(j,k)=U(j,k)/U(k,k);           %j-th Multiplier m_jk
        U(j,k:n)=U(j,k:n)-L(j,k)*U(k,k:n);    %Updating U
        'U after'
        U
        ['************************************']
        pause
    end
end