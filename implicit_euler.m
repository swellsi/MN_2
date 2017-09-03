function [z]=implicit_euler(z0,k,it)
% z0 es vector vertical (x;y)
% z es la solucio [totes les x; totes les y];
% cal canviar la funcio derivada de dins de IEULER, fd(z) per la funcio que
% et calcula la derivada del teu problema en concret
z=z0;
for i=2:it
    [xfinal,xk,resk,it]=newtonraphson_impliciteuler(z(:,i-1),@ieuler,1e-8,50,2,2,z(:,i-1),k);
    z(:,i)=xfinal;
end