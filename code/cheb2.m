% SECOND ORDER CHEBYSHEV DIFFERENTIATION MATRIX (Formula (3.4.37))

% cheb2 compute D^2 = second-order differentiation matrix
%               x = Chebyshev grid

% Based on El-Baghdady, G. I., & El-Azab, M. S. (2016). 
% Chebyshev-Gauss-Lobatto Pseudo-Spectral Method for One-Dimensional 
% Advection-Diffusion Equation with Variable Coefficients. Math, 3(1), 1-8.

function [D2,x] = cheb2(N)
[D,x] = cheb(N);   % First order Cheb. differentiation matrix and Cheb. pts
for i = 1:N+1;     % Formula (3.4.37)
    D2(i,i) = 0;
    for j = 1:N+1;
        if i~=j;                                        
            D2(i,j) = 2*D(i,j)*(D(i,i)-(1/(x(i,1)-x(j,1))));
            D2(i,i) = D2(i,i) - D2(i,j);
        end
    end
end
end

