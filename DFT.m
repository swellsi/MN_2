function [T,Fk]=DFT(fj)
%N is an even number
N=length(fj); h=2*pi/N; T=[];
WN=exp(-i*h);
for k=0:N-1
for j=0:N-1
T(k+1,j+1)=(1/N)*WN^((k-(N/2))*j);
end
end
Fk=T*fj;
end