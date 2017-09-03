function [x,k] = mygmres(Afun,b,tol,dimkryl)
% GMRES function using Gram-Schimdt-Arnoldi reotho

% GOAL: min(c es R^n) de norm(A*Kn*c-b) on xn=Kn*c;
% Kn es la matriu de la base de xn, Kn=[b Ab A^2b ... A^(n-1)b]
% REORTHOGONALIZATION OF Kn matrix with ARNOLDI'S ITERATION
% Q matriu de la nova base reortogonalitzada, vectors en columnes
H=[]; % PARTIAL REDUCTION TO THE UPPER HESSEMBERG FORM
      % A*Q_n=Q_n+1*H_n
hkp1=1; resd=1; k=1;
Q(:,1)=b/norm(b);
while resd>tol & hkp1>eps & k<=dimkryl
    % ANEM CONSTRUINT LA Q (nova base de xn) i H (matriu Hessiana)
    Q(:,k+1)=feval(Afun,Q(:,k));
    h=Q(:,1:k)'*Q(:,k+1);
    Q(:,k+1)=Q(:,k+1)-Q(:,1:k)*h;
    % reorthogonalization
    s=Q(:,1:k)'*Q(:,k+1);
    Q(:,k+1)=Q(:,k+1)-Q(:,1:k)*s;
    h=h+s; % updating h used with s used from reorth.; qk=vk-(r+s)q(1:1-k)
    % if modul 0... no afegeix dimensio, 
    hkp1=norm(Q(:,k+1));
    H(1:k+1,k)=[h;hkp1];
    Q(:,k+1)=Q(:,k+1)/hkp1;
    
    % Kn*c = Qn*y
    % NEW GOAL: min(y es R^n) de norm(Hn*y-Q_n+1'*b) on Hn=Q_n+1'*A*Qn
    [c,r]=qr(H(1:k+1,1:k)); % C es la Q i r es la R de la part de H
                            % Q_n+1'*b=Q_n+1'*q1*norm(b)=norm(b)*e1
                            % Hn*y-Q_n+1'*b=c*r*y-norm(b)*e1
                            % si ho volem minimitzar, busquem y tq faci =0
                             % c*r*y=norm(b)*e1 ---> r*y=c'*norm(b)*e1
                            % y=r\c'*norm(b)*e1
    y=r\(c'*norm(b)*[1;zeros(k,1)]);
    x=Q(:,1:k)*y; %=Kn*c
    resd=norm(feval(Afun,x)-b);
    k=k+1;    
end


    
    