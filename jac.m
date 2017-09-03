%Computation of the Jacobian
%Matrix of a F field: F:R^n ---> R^n
%Warning!: h takken identical to ALL variables
%On input:
%  1) @name of the function F
%  2) number of independent variables of F: m
%  3) number of components of F: n
%  4) the point x (column m-vector) in R^m at which DF is computed
%On output:
%  DF(x) n x m - Jacobian matrix at x

function [DF]=jac(F,m,n,x)
I=eye(m); DF=zeros(n,m); h=sqrt(eps);
for j=1:m
    f2=feval(F,x+I(:,j)*h); f1=feval(F,x-I(:,j)*h);
    DF(:,j)=(f2-f1)/(2*h);
end
