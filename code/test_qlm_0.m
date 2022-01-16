%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% QUASI - LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.

% Chebyshev coefficients for n = 1/10000 [thick solid]. The asymptote for
% large degree m to roughly 0.001/m^5. The exact solution for n=0 is a
% parabola whose three Chebyshev coefficients are marked by the black disk.
% The perturbation of changing n from 0 to 1/10000 has modified the
% Chebyshev coefficients by rendering nonzero all coefficients of degree
% three and greater. However, because of the sum of a series of Chebyshev
% coefficients is bounded by the absolute values of the includad
% coefficients, it follows just from inspection of the tiny magtinutde of
% a_{3}, a_{4},..., that the algebraically converging sum will modify y(x)
% by less than 1e-4. Although the perturbation has drastically changed the
% quantitative appearance of a graph of Chebyshev chefficients is still
% very well approximated by the parabola. 

% (Figure 4.1)

format longg

N = 100;                                      % Collocation points
iters = 30;                                   % Number of Newton iterations
n = 1/10000;                                  % Polytropic index

f = @(x) 0.001./x.^(5);                       % Convergence
x = linspace(1,N,N);
y = f(x);


[xi, Y, XCheb, a] = qlm(N,n,iters);           % Function call to QSLM
loglog(1:N, abs(a),'k-')                      % log scale in both axles 
hold on
plot(1:N, y, 'r--')                           % Plot convergence lineguide 

hold on
plot([1 2 3],abs(a(1:3)),'-o',...             % Parabola coefficients
    'MarkerFaceColor','black',...
    'MarkerSize',8)

hold on                                       % Display coeffs values
xt = [1 2 3];
yt = [abs(a(1)) abs(a(2)) abs(a(3))];
str = {num2str(a(1)), num2str(a(2)), num2str(a(3))};
text(xt, yt, str);

xlabel('Degree');
ylabel('Absolute value of coefficients.'); 
text(6, 10e-10,'$\frac{0.001}{m^{5}}$','Interpreter','latex',...,
    'fontsize',18, 'Color', 'red')

title('Chebyshev coefficients for n = 1/10000.');

