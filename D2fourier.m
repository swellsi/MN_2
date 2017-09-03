function [D2]=D2fourier(N)
k=[-N/2:N/2-1];
D2=diag(-k.^2);
end