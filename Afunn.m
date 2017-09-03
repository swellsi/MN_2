function Ax=Afun(x)
% li entra un vector x i li calcula l'ACCIO de la matriu al vector; no cal
% construir la matriu:))
Ax=[];
R=3; V=1; m=15; Ax=[];
Ax=[Ax; 2*R*x(1) + R*x(2)];
Ax=[Ax; x(1)-x(2)-x(3)];
for j=3:m
    Ax = [Ax ; 2*R*x(2*j-3) + R*x(2*j-2) - R*x(2*j-4)];
    Ax = [Ax ; x(2*j-3) - x(2*j-2) - x(2*j-1)];
end
Ax = [Ax ; 3*R*x(2*m-1) - R*x(2*m-2)];
