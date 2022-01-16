% Comparison:
% Fourth order convergence of the finite difference differentiation process
% vs
% Spectral accuracy (case of Fourier methods)

% Based on Trefethen, L. N. (2000). Spectral methods in MATLAB. Society 
% for industrial and applied mathematics. 

% (Figure 3.1)

% Convergence of fourth-order finite differences for exp(sin(x))

% For various N, set up grid in [-pi,pi] and function u(x):
Nvec = 2.^(3:12);
clf, subplot('position',[.1 .4 .8 .5])

for N = Nvec
    h = 2*pi/N;  
    x = -pi + (1:N)'*h;   
    u = exp(sin(x));                                % Function exp(sin(x))
    uprime = cos(x).*u;                             % Analytic derivative

    % Construct sparse fourth-order differentiation matrix
    e = ones(N,1);
    D =   sparse(1:N,[2:N 1],2*e/3,N,N)...
        - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf);
    loglog(N,error,'r.','markersize',15), hold on
end

semilogy(Nvec,Nvec.^(-4),'--') 
text(105,5e-8,'N^{-4}','fontsize',18)

% Convergence of periodic spectral method for exp(sin(x))

% Function exp(sin(x)) gives periodic data on the domain [-pi, pi], so the
% natural choice is a trigonometric polynomial on an equispaced grid to
% implement 'Fourier' methods'. However, for non-periodic domains,
% algebraic polynomials on irregular grids are the right choice. 

% For various N (even), set up grid:
for N = 2:2:100;
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
    u = exp(sin(x));                                % Function exp(sin(x))
    uprime = cos(x).*u;                             % Analytic derivative
    
    % Construct spectral differentiation matrix.
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    D = toeplitz(column,column([1 N:-1:2]));

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf);
    loglog(N,error,'b.','markersize',15), hold on
end

grid on, xlabel N, ylabel error
title({'Convergence of finite diference method (red) vs',
    'convergence of spectral differentiation (blue).'})