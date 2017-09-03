function [F]=ieuler(z,z0,k)
%fder=[(1/z(2))-1 ; 1-(1/z(1))]
%F=-z+z0+k*[(1/z(2))-1 ; 1-(1/z(1))];
F=-z+z0+k*fd(z);
end

function z_der=fd(z)
% z vector (x,y)
% z vector (x',y')
z_der(1)=1./z(2)-1;
z_der(2)=1-1./z(1);
z_der=z_der';
end