function [H]=hamiltonian(x,y)
gamma=1; H0=-2.601;
H=log(y*x^gamma)-x-y-H0;
end