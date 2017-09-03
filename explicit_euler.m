function z=explicit_euler(z0,deltat,it)
% z0 es vector vertical (x;y)
% z es la solucio [totes les x; totes les y];
% fd es la funcio que deriva fd([x;y])=[x';y']
global Re; global N;
z=z0; k=[0:N-1]'; 
dreal=Dreal(N); dft=DFT_T(N); idft=IDFT_C(N);  
global gn
for i=2:it
    g(:,i-1)=dft*(dreal*(.5*((idft*z(:,i-1))).^2));
    gn=g(:,i-1);
    z(:,i)=(z(:,i-1)+deltat*gn)./(1+(deltat*(k-N/2).^2)./Re);
end

