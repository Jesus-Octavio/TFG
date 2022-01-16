%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.

function [Yz, xis, T] = slm(N,n,iters)

% INPUT:
% N        Natural.            N+1 is the No. of Chebyshev collocation pts.
% n:       Real in (0,5).      Polytropic index
% iters:   Natural.            No. of iterations


% OUTPUT:
% Yz:            (N+1)x(N+1)  double array.         y_{i} value at grid.
% xis:           (iters+1)x 1 double array.         xi_{k}.
% T:             (N+1)x1      double array          gridpoints t_{i}.

format longg

Yz = zeros(N+1, iters+1);
xis = zeros(iters+1, 1);

[D,Z] = cheb(N);      % Chebyshev differentiation matrix and Chebyshev grid
D_2 = cheb2(N);       % Square of Chebyshev differentiation matriz
T = slm_transform(Z,N);  % Transform region [-1,1] into [0,1]: grid -> grid
Diag = diag(4./T);    % Diagonal matrix of elements 4/T(i)


y1 = @(z) 1-((z+1)/2)^2;                          % Initial guess function.
for k=1:N+1 
      Yz(k,1) = y1(Z(k));
end        
xis(1) = sqrt(6*((4/3)^(n)));            % Initial guess xi_{0} root.

for i = 1:iters    
    
    % Matrix r_{i-1}=[r_{i-1}(z_{0}), r_{i-1}(z_{1}) ..., r_{i-1}(z_{N})]'
    % Formula (4.2.13)
    R_i1 = slm_coefsR(N,n,i+1,Z,D,D_2,xis,Yz);
    
    % Matrix a_{1,i-1}
    % Formula (4.2.13)
    a_i1 = slm_coefsA(N,n,i+1,xis,Yz);
            
    % Matrix a_{2,i-1}
    % Formula (4.2.13)
    B_i1 = slm_coefsB(N,n,i+1,xis,Yz);
       
    % Matrix A_{i-1}
    % Formula (4.2.13)
    A_i1 = 4.*D_2 + Diag*D + a_i1;
        
    
    % M_1 coefficient matrix
    % M_2 vector
    % Ignore 1st and (N+1)st rows of the coefficient matrix and vector
    % Ignore 1sr and (N+1)st columns of the coefficient matrix
    % System (4.2.15)
    M_1 = [A_i1(2:N,2:N) B_i1(2:N); D(N+1,2:N) 0];
    M_2 = [R_i1(2:N); 0];
    
    % Solve system (4.2.15)
    M_3 = M_1\M_2;
        
    % Separate results and update matrices with results
    Yz(2:N,i+1) = M_3(1:N-1);
    xis(i+1) = M_3(N);
end
end