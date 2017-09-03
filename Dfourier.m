function [D]=Dfourier(N)
k=[-N/2:N/2-1];
D=diag(i*k);
end