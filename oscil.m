function [F]=oscil(phi,alfa)
phi1=phi(1); phi2=phi(2);
F1=tan(phi1)-alfa*(2*sin(phi1)+sin(phi2));
F2=tan(phi2)-2*alfa*(sin(phi1)+sin(phi2));
F=[F1;F2];
end

