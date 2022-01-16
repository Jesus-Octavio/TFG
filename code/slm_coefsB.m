%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.

% coefsB compute B = matrix B_{1,i-1} (Formula (4.2.13))

function [B] = slm_coefsB(N,n,i,xis,Yz)
format longg

B = zeros(N+1,1);

factor_1 = sum(xis(1:i-1)); 

for j = 1:N+1 %grid points z_{j}
      
    factor_2 = sum(Yz(j,1:i-1));
    factor_2 = factor_2 ^(n);

    B(j) = 2*factor_1*factor_2;
end
end