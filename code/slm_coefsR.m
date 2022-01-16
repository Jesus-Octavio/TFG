%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.

% coefsR compute A = matrix R_{1,i-1} (Formula (4.2.13))

function [R] = slm_coefsR(N,n,i,Z,D,D_2,xis,Yz)
format longg

R = zeros(N+1,1);

factor_3 = sum(xis(1:i-1)); % Sum of roots xis(1)+...+xis(i-1)
factor_3 = factor_3^(2);

for j = 1:N+1  % Grid point z_{j}
    
    factor_1 = 0;
    for m = 1:i-1 % Calculate second derivative at z_{j]
        for k = 1:N+1
            factor_1 = factor_1 + D_2(j,k)*Yz(k,m);
        end
    end
    
    
    factor_2 = 0;
    for m=1:i-1 %Calculate first derivative at z_{j]
        for k = 1:N+1
            factor_2 = factor_2 + D(j,k)*Yz(k,m);
        end
    end
    %factor_2 = (2/Z(j))*factor_2; %what if passing Z?
    factor_2 = (2/((Z(j)+1)/2))*factor_2; %Trnasform z_{j} into t_{j}
    
    
    factor_4 = sum(Yz(j, 1:i-1)); 
    factor_4 = factor_4^(n);
        
    R(j) = -(factor_1 + factor_2 + ((factor_3)*(factor_4)));    
end
end
    