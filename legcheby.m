function [x,w,lambda,D]=legcheby(N)
w=zeros(N,1); lambda=zeros(N,1); D=zeros(N,N);

j=[0:N-1]'; x=cos((pi*j)./(N-1));
v=ones(N-2,1)*(pi/(2*(N-1))); w=[pi/(N-1);v;pi/(N-1)];

lambda(1)=(2^(N-1-2))/(N-1);
lambda(N)=(((-1)^(N-1))*2^(N-1-2))/(N-1);
for j=1:N-2; lambda(j+1)=(((-1)^j)*2^(N-1-1))/(N-1); end

for j=1:N for i=1:N  
    if i~=j; D(i,j)=lambda(j)/(lambda(i)*(x(i)-x(j)));
    else for k=1:N if k~=j; D(j,j)=D(j,j)+((x(j)-x(k))^(-1)); end; end; end 
end; end
end