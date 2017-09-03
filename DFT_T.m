function [T,Fk]=DFT_T(N)
%N is an even number
h=2*pi/N; T=[];
WN=exp(-i*h);
for k=0:N-1
for j=0:N-1
T(k+1,j+1)=(1/N)*WN^((k-(N/2))*j);
end
end
end