
function [ xfinal,xk,resk,it ] = newtonraphson(x,fun,tol,itmax,m,n)
it=1; xk=x; resk=max(abs(feval(fun,x))); tolk=1;
while it<itmax & tolk>tol;
    x=xk(:,end);
    J=jacobian(fun,x,m,n);
    b=-feval(fun,x);
    [P,L,U]=PALU(J); % PJ=LU ---> P-1*L*U*Axk=-F(xk)=b ---> LUAxk=Pb=c
    c=P*b; d=L\c; Axk=U\d;
    xk=[xk Axk+x];
    it=it+1; tolk=max(abs((xk(:,end)-xk(:,end-1))));
    resk=[resk; max(abs(feval(fun,xk(:,end))))];
end
xfinal=xk(:,end);
end

function DF = jacobian(F,x,m,n)
% x vector de variables vertical
% m = numero de variables (x1,...,xm)
% n = numero de funcions (F1,...,Fn)
h=sqrt(eps); DF=zeros(n,m);
I=eye(m);
for i=1:n
    DF(:,i)=1/h*(feval(F,x+I(:,i)*h)-feval(F,x));
end
end

function [P,L,U]=PALU(A) % factorization (WTIH pivoting)
% PA=LU; L = lowertriangular; U = uppertriangular; P = exhanges matrix
[m,n]=size(A);
if m~=n; error('not square matriz'); end
L=eye(n); U=A; %U sera la A que es va transformant, A(2), A(3),.. fins A(n)
P=eye(n);
for k=1:(n-1) %iteracio
    % CANVI DE PIVOT
    [foo imax]=max(abs(A(k:end,k))); %trobem el pivot m�xim
    iimax=imax+k-1; %index de la fila en A que te maxima a_jk
    if iimax~=k %si la a_jk maxima (a_iimax,k) no �s la primera (a_kk)
        %canviem la fila de la a_iimax,k per la de a_kk
        i1=[k iimax]; i2=[iimax k];
        U(i1,:)=U(i2,:); %row-exchange U
        P(i1,:)=P(i2,:); %update P
        L(i1,1:(k-1))=L(i2,1:(k-1)); %low triangle row exchange for L,
                                    %no volem fer row-exch de la identitat!
    end %no caldria if, ja que si iimax=k les operacions no canvien U,P :)
    %GAUSS SUBSTITUTION
    for j=(k+1):n %fila
        L(j,k)=U(j,k)/U(k,k); %�s la m_jk, diferent per cada fila
        U(j,k:n)=U(j,k:n)-L(j,k)*U(k,k:n); %updating U
    end
end
end