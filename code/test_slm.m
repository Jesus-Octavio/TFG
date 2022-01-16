%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUCCESSIVE LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Motsa, S. S., & Shateyi, S. (2012). A successive Linearization 
% Method Approach to Solve Lane-Emden Type of Equations. 
% Mathematical Problems in Engineering, 2012.


N = 100;                                             % No. of Chebyshev pts
n = 1;                                               % Polytropic index
iters = 10;                                          % No. of iterations

[Yz, xis, T] = slm(N,n,iters);

fprintf('First zero (xi) is %s', sum(xis))

Y = sum(Yz,2);
x_values = sum(xis).*T;
plot(x_values, Y, 'r-*');
met = 'SLM';
s=sprintf('Solve Lane-Emden with ,met=%s, chebyshev nodes N=%g,polytropic index n=%g,iterations=%g',met,N,n,iters);
title({' ','Successive linearization method (SLM)',s,});
grid on

