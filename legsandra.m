function [x,w,lambda,D]=legsandra(N)
beta=zeros(N,1); x=zeros(N,1); w=zeros(N,1); lambda=zeros(N,1); D=zeros(N,N);

for j=1:N; beta(j)=j./((4*(j.^2)-1).^.5); end

diag=[beta,zeros(N,1),[0;beta(1:N-1)]];
T=full(spdiags(diag,[-1:1],N,N))
[v,d]=eig(T);

for j=1:N; x(j)=d(j,j); w(j)=2*((v(1,j))^2); end
for k=0:N-1; lambda(k+1)=((-1)^k)*sqrt((1-(x(k+1))^2)*w(k+1)); end

for j=1:N for i=1:N  
    if i~=j; D(i,j)=lambda(j)/(lambda(i)*(x(i)-x(j)));
    else for k=1:N if k~=j; D(j,j)=D(j,j)+((x(j)-x(k))^(-1)); end; end; end 
end; end
end
