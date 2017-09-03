function [C,fj]=IDFT(Fk)
%N is an even number
N=length(Fk); h=2*pi/N; T=[];
WN=exp(-i*h);
for j=0:N-1
for k=0:N-1
C(j+1,k+1)=WN^(-((k-(N/2))*j));
end
end
fj=C*Fk;
end