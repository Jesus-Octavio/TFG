%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.

% coefsA compute A = matrix a_{1,i-1} (Formula (4.2.13))

function [A] = slm_coefsA(N,n,i,xis,Yz)
format longg

Aaux = zeros(N+1,1);

factor_1 = sum(xis(1:i-1));
factor_1 = factor_1^(2);
    
for j = 1:N+1
    
    factor_2 = sum(Yz(j,1:i-1));
    factor_2 = factor_2^(n-1);
      
    Aaux(j) = n*factor_1*factor_2;

end
A = diag(Aaux);
end