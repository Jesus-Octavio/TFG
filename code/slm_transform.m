%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.

% Transform the chebyshev collocation point
% on physical region [-1,1] into colloction points on [0,1].
% This can be achieved using t=(z+1)/2

% Transform Chebyshev grid z into grid t.

function [T] = transform(Z,N)
format longg
T = zeros(N+1,1);
for i = 1:N+1
    T(i) = (Z(i) + 1)/2;
end
end