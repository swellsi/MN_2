function [D]=Dreal(N)
h=2*pi/N; D=zeros(N);
for j=1:N; for l=1:N;
        for k=0:N-1
        D(j,l)=D(j,l)+(1/N)*i*(k-N/2)*exp(i*(k-N/2)*h*(j-l));
        end
end; end
end
