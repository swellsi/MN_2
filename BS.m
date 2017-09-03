function [x]=BS(U,b)
m=size(U); n=m(1); x=zeros(n,1); x(n)=b(n)/U(n,n); 
for i= n-1:-1:1
    x(i)=(1/U(i,i)).*(b(i)-(U(i,i+1:n)*x(i+1:n)));
end
  