function [x]=resoldresistemaambPALU(A,b)
[P,L,U]=palu(A);           %Ax=b; P^(-1)LUx=b; LUx=Pb; Ux=y; 
c=P*b;                     %Ly=Pb=c
[y]=FS(L,c);     %-->y     %Ux=y
[x]=BS(U,y);     %-->x
end