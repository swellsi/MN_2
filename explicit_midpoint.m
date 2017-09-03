function [z]=explicit_midpoint(initial_values,k,it,fd)
% z matriu de posicions (totes les x; totes les y)
% size(z)=(2,it)
% initial_values son la primera posicio necessària (es d'ordre 1) = [x0 x1;y0 x1]
% k es el pas (aprox. 0.02)
% it es el nombre d'iteracions, el nombre de passos
% fd es la funcio derivadora; fd([x;y])=[x';y']
% fd HA DE RETORNAR I REBRE UN VECTOR DIAGONAAAAAL!!!!!!! 
z=initial_values; 
for i=3:it
    z(:,i)=z(:,i-2)+2*k*fd(z(:,i-1));
end
