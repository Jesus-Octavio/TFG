%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% QUASI - LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.

% Calculate dy(1)/dx 

% According to Boyd's article, derived quantities such as the first and
% second derivatives can be avaluated at an arbitrary point by using the
% trigonometric functions for D0Mat, D1Mat, D2Mat (2.10 in article) with
% t_{i} replaced by t=arccos(x). There's another way: using the Cjhebyshev
% coefficients and calculating D1Mat*a. The first derivative grid points
% values of y are given by D0Mat*a.

% (Table 4.5)

format longg                                  % Variable precision
n = 2;                                        % Polytropic index
N = 100;                                      % Number of collocation pts
iters = 10;                                   % Number of Newton iterations
 
[xi, Y, XCheb, a, D] = qlm_first_derivative(N,n,iters);
R = xi*XCheb;
deriv = D*a;
disp('dy/dx(1) =  ');
disp(deriv(1))
