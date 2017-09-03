function [t]=funciotonta(x)
a=x(1,1);b=x(2,1);c=x(3,1);
t1=3*a-sin(b)/c; t2=exp(c/a)+a*b*c-tan(sin(exp(a*b*c)));
% t1=a+1; t2=b-1; 
t3=a+b; t4=c-a;
t=[t1;t2;t3;t4];
end