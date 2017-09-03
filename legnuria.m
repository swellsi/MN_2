% GAUSS-LEGENDRE APROXIMATION THEORY
% Provides:
% 1) gl_nodes x=(x0,...,xn)
% 2) quadrature weights w=(w0,...,wn)
% 3) barycentric weights (lambda0,...,lambdan)
% 4) differentiation matrix Dij
% for an arbitrary number N=n+1 nodes

function [x,w,lamb,D] = legnuria(N)

% NODES AND WEIGHTS (using the companion matrix technique)
j=[1:N]';
beta=j./(4.*j.^2-1).^0.5;
diag=[beta,zeros(N,1),[0;beta(1:N-1)]];
T=full(spdiags(diag,[-1:1],N,N));
[eigenvectors,eigenvalues]=eig(T);
x=spdiags(eigenvalues);       % vertical vector of nodes   (N,1)
w=2.*eigenvectors(1,:).^2';   % vertical vector of weights (N,1)

% BARYCENTRIC WEIGHTS (using the magical formula)
j=[0:N-1]';
lamb=(-1).^j.*((1-x.^2).*w).^0.5;

% DIFFERENTIATION MATRIX
% Atencio! els nodes van de [0:n], pero els vectors de nodes, pesos,... van
% de [1:N] en matlab!! (N=n+1)
D=[];
for i=1:N
    for j=1:N
        if i~=j
            D(i,j)=lamb(j)/(lamb(i)*(x(i)-x(j)));
        else
            xnew=[x(1:j-1);x(j+1:N)];
            D(i,j)=sum(1./(x(j)-xnew));
        end
    end
end

end

            