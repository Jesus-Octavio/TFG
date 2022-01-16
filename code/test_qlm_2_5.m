%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% QUASI - LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.

% Chebishev coefficients fo polytropic index n = 5/2. The left plot, which
% uses a linear sclae in degree but a logarithmic scale for the
% coefficients, shows that the leading coefficients fall geometrically.
% The dashed line is the graph of 2.5exp(-1.12N). The right graph is the
% same butwith both scales logarithmic, and with the range of coefficients
% extened. The coefficients of degree 27 to 400 are well-fit by the dashed
% line, indistinguishable from the graph of the coefficients, 75/N^10. The
% other dashed guidelines, clearly not oof the same slope as the thickk
% curve of the coefficients, have slopes proportional to N^-9, N^-11

% (Figure 4.2)

format longg

N = 400;                                      % Collocation points
iters = 50;                                   % Number of Newton iterations
n = 2.5;                                      % Polytropic index

f = @(x) 2.5*exp(-1.12*x); 
f_2 = @(x) 75./x.^(10); 
f_3 = @(x) x.^(-9); 
f_4 = @(x) x.^(-10); 


[xi, Y, XCheb, a] = qlm(N,n,iters);          % Function call to QSLM                            

subplot(1,2,1)
semilogy(1:N,abs(vpa(a)),'r-')   
hold on

x = linspace(1,100,400);
y = f(x); 

subplot(1,2,1)
semilogy(x,y,'k--');

xlim([0 40]);
xlabel('Degree.')
ylim([10e-15 10e0]);
ylabel('Chebyshev coefficients');
title('`n=2/5: semilog plot')
hold on



subplot(1,2,2);
loglog(1:N,abs(a),'r-');
hold on

x = linspace(27,200,400);
y_2 = f_2(x); 
y_3 = f_3(x); 
y_4 = f_4(x); 
subplot(1,2,2)
loglog(x,y_2,'k--');
hold on
subplot(1,2,2)
loglog(x,y_3,'b--');
hold on
subplot(1,2,2)
loglog(x,y_4,'g--');
hold on
xlim([0 100]);
xlabel('Degree.')
ylim([10e-25 10e0]);
ylabel('Chebyshev coefficients')
title('loglog plot')