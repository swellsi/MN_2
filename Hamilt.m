function [H]=Hamilt(vector_posicio)
gamma=1;
H0=-2.601;
x=vector_posicio(1); y=vector_posicio(2);
H=log(y.*(x).^gamma)-x-y-H0;
end
