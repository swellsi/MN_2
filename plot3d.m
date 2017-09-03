function [XX,YY,W] = plot3d(xx,u,n)
Nt = 100 ; tt = [0:2*pi/100:2*(Nt)*pi/Nt];
[RR,TT] = meshgrid(xx,tt); [XX,YY] = pol2cart(TT,RR);
U = repmat(u,[Nt+1 1]);
v = cos(n*tt).'; V = repmat(v,[1 length(u)]);
W = U.*V ;  W=sign(W(1,end))*W/max(max(abs(W)));
mesh(XX,YY,W) ; view(26,15); colormap([0 0 0]) ; daspect([1 1 1.5]);
axis([-1.25 1.25 -1.25 1.25 -.5 .5])
axis off