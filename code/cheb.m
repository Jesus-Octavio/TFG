% FIRST ORDER CHEBYSHEV DIFFERENTIATION MATRIX (Theorem 14)

% cheb compute D = differentiation matrix, 
%              x = Chebyshev grid

% Based on Trefethen, L. N. (2000). Spectral methods in MATLAB. Society 
% for industrial and applied mathematics. 

function [D,x] = cheb(N)
    if N == 0, D = 0; 
        return, 
    end                              
    x = cos(pi*(0:N)/N)';                            % Chebyshev pts.
    c = [2; ones(N-1,1);2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';
    D = (c*(1./c)')./(dX+(eye(N+1)));                % Off-diagonal entries
    D = D - diag(sum(D'));                           % Diagonal entries
end

