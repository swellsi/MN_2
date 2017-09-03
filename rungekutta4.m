function [v]=rungekutta4(k,T,v1,funvt)
%k: step size
t = [0:k:T]';
v = zeros(length(v1),length(t)); 
v(:,1)=v1; %v1: column initial condition
%funvt: 1st variable is v, then t

for n=1:(length(t)-1)
    a = k*funvt(v(:,n),t(n));
    b = k*funvt(v(:,n)+0.5*a,t(n)+0.5*k);
    c = k*funvt((v(:,n)+0.5*b),(t(n)+0.5*k));
    d = k*funvt((v(:,n)+c),(t(n)+k));

    v(:,n+1) = v(:,n) + (1/6)*(a+2*b+2*c+d); 
end
end