%Matrix of a F field: F:R^n ---> R^n
%On input:
%  1) @name of the function F
%  2) number of independent variables of F: m
%  3) number of components of F: n
%  4) point xo (column m-vector) in R^m at which Newton's method initiates

function [x,e,it]=Newtonvec(F,m,n,xo,emax,itmax)
x=xo; it=0; e=1;
while e>emax &it<itmax
    A=jac(F,m,n,x(:,end));  deltaxk=QRsolveX(A,-F(x(:,end)));
    x=[x, x(:,end)+deltaxk];
    e=norm(x(:,end)-x(:,end-1)); it=it+1;
end
end