function z_der=fd(z)
% z vector (x,y)
% z vector (x',y')
z_der(1)=1./z(2)-1;
z_der(2)=1-1./z(1);
z_der=z_der';
end
