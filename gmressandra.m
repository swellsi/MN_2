function [x,k]=gmressandra(Afun,b,tol,dimkryl)
% GMRES function using Gram-Schmidt-Arnoldi reorthogolanization
%Matrix-free version
%Afun: you have to provide an M-file capable of computing the linear action
%x -----> A*x
%b: the rhs
%tol: the iteration stops if |A*x_n-b|_2 < tol
%x: approximation found
%k: number of iterations/Krylov dimension used

hkp1=1; resd=1; k=1;
Q(:,1)=b/norm(b); H=[];              %q_1=b/|b|

while resd>tol & hkp1>eps & k<=dimkryl
    disp(['****************** k=' int2str(k) '*********************'])
    Q(:,k+1)=feval(Afun,Q(:,k));    %q_k+1=Aq_k
    h=Q(:,1:k)'*Q(:,k+1);           %[h_1k ... h_kk]=q_{1:k}'*q_{k+1}
    Q(:,k+1)=Q(:,k+1)-Q(:,1:k)*h;   %Sum
        %Reorthogonalization:
        S=Q(:,1:k)'*Q(:,k+1); Q(:,k+1)=Q(:,k+1)-Q(:,1:k)*S;
        h=h+S;
        disp(['-----------h:-------------'])
        disp(h)
        disp(['--------------------------'])
    hkp1=norm(Q(:,k+1));            %h_{k+1,k}=|q_{k+1}|
    H(1:k+1,k)=[h;hkp1];            %[h_1k ... h_kk h_k+1,1]'
    Q(:,k+1)=Q(:,k+1)/hkp1;         %q_{k+1}=q_{k+1}/h_{k+1,k}
    [c,r]=qr(H(1:k+1,1:k));         %QR(-Hn) cheating!
    y=r\(c'*norm(b)*[1;zeros(k,1)]); x=Q(:,1:k)*y; %x=QK*y
    resd=norm(feval(Afun,x)-b);
    disp(['----------- -H: ------------'])
    disp(H)
    disp(['----------------------------'])
    disp(['RESIDUAL=' num2str(resd)])
    disp(['----------- x: -------------'])
    disp(x)
    disp(['----------------------------'])
    disp(['*************************************************'])
    %pause
    k=k+1;
end
    
    
    