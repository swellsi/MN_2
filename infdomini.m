function [aDy,aDy2,Dy,Dy2]=infdomini(L,N)
% per canvis del tipus y=Lx/sqrt(1-x.^2);
% dom(y)=(-inf,inf); dom(x)=[-1,1];
% x i D de chebyshev
[x,w,lambda,D]=legcheby(N);
Q=1-x.^2;

aDy=diag((Q.^(3/2))./L)*D;
aDy2=diag((Q.^3)./L^2)*D*D-diag((3*x.*Q.^2)/L^2)*D;
Dy=aDy([2:end-1],[2:end-1]);
Dy2=aDy2([2:end-1],[2:end-1]);
end