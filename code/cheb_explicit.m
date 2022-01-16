% FIRST ORDER CHEBYSHEV DIFFERENTIATION MATRIX (Theorem 14)

% cheb_explicit compute D = differentiation matrix, 
%                       x = Chebyshev grid

% Compute first order Chebyshev differentiation matrix using explicit
% formulas in Theorem 14


function [D, Z] = cheb_explicit(N)

Z = cos(pi*(0:N)/N)';                                    % Chebyshev pts.

for i = 1:N+1
    for j = 1:N+1
        
        if i == 1 && i == j;                             % Formula (3.2.34)                             
            D(i,j) = (1/6)*(2*(N^2)+1);
        
        elseif i == N+1 && i == j;                       % Formula (3.4.34)
            D(i,j) = -(1/6)*(2*(N^2)+1);
        
        elseif i == j;                                   % Formula (3.4.43)
            D(i,j) = -Z(j)/(2*(1-Z(j)^2));
        
        else                                             % Formula (3.4.42)
            
            if i == 1 || i == N + 1;
                c_i = 2;
            else
                c_i = 1;
            end
            
            if j == 1 || j == N + 1;
                c_j = 2;
            else
                c_j = 1;
            end     
            D(i,j) = (c_i / c_j)*((-1)^(i+j)/(Z(i)-Z(j)));
        
        end
    end
end
end

