function[x,e,it]=Newton(f,xo,emax,itmax)
format long; format compact
it=0; e=1; x=xo; dx=1e-7;
while e>emax & it<itmax
    x=[x;x(end)-(feval(f,x(end))/((feval(f,x(end)+dx)-feval(f,x(end)))/dx))];
    e=abs(x(end)-x(end-1));
    it=it+1;
end