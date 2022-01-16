%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% QUASI - LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.

% Calculate first derivative


function [xi, Y, XCheb, a, D1Mat] = qlm(N,n,iters)


% [xi, Y, XCheb, a] = qlm(N,n,iters) computes the Quasi-linearization
% method for the Lane-Emden equation.

% Input:
% N:        Natural.                  No.r of nodes.
% n:        Real in [0,5].            Polytropic index.
% iters:    Natural.                  No. of Newton-Kantorovich iterations.

% Output:
% xi:       Double.                   Eigenvalue.
% Y:        Nx1 array of double       Value at grid.
% XCheb:    Nx1 array of double       Gridpoints.
% a:        Nx1 array of double       Chebyshev serie coefficients.
% D1Mat :   N+N array of double       1st order differentiation matrix.

    % Initialize matrix
    Y0 = zeros(N,1);

    % Compute grid points, XCheb 
    i = 1:N;
    ta = pi*(i-1)/(N-1);   
    XCheb = 0.5.*(1+cos(ta))';
    
    % Compute differentiation matrices
    t = ta(2:N-1);
    ss = sin(t);
    cc = cos(t);
    
    D0Mat(2:N-1,1:N) = cos(t.*(0:N-1)')';
    pt = -(0:N-1).*sin(t.*(0:N-1)')';
    ptt = -((0:N-1).^2).*D0Mat(2:N-1,1:N);
    D1Mat(2:N-1,1:N) = (-2.*pt./ss');
    D2Mat(2:N-1,1:N) = 4.*(ptt./(ss'.^2)-cc'.*pt./(ss'.^3));
          
    % Apply non-trig formulas at the endpoints
    j = 1:N;
    D0Mat(1,j) = 1;
    D0Mat(N,j) = (-1).^(j+1);
    D1Mat(1,j) = 2.*(j-1).^2;
    D1Mat(N,j) = 2.*(-1).^(j).*(j-1).^2;
    D2Mat(1,j) = (j-1).^(2).*((j-1).^(2)-1).*(4/3);
    D2Mat(N,j) = -(-1).^(j).*(j-1).^(2).*((j-1).^(2)-1).*(4/3);
           
    %Initial guess
    Y0 = cos((pi/2)*XCheb);
    xi0 = 3;
    a0 = D0Mat\Y0;
    a = a0;
    xi = xi0;

    % Newton-Kantorovich iterations
    for iter=1:iters

        Y = D0Mat*a;

        % Compute residual          
        Res = -xi^(2)*Y(2:N-1).^n;
        Res =  Res - D2Mat(2:end-1,:)*a - (2./XCheb(2:end-1)).*...
            D1Mat(2:end-1,:)*a;
                        
        Ytzero = 0;
        Ytzero = Ytzero + D0Mat(N,:)*a;
        Res(N+1) = -(Ytzero-1);


        % Compute Jacobian matrix       
        Jacobian = D2Mat(2:end-1,:)+(2./XCheb(2:end-1)).*D1Mat(2:end-1,:)+...
            ((xi^2)*n*Y(2:end-1).^(n-1)).*D0Mat(2:end-1,:);
                    
        j = 1:N;
        Jacobian(N-1,j) = D0Mat(1,j);
        Jacobian(N,j) = D1Mat(N,j);
        Jacobian(N+1,j) = D0Mat(N,j);
                
        i = 1:N-2;
        Jacobian(i,N+1) = 2*xi*Y(i+1).^n;
        
        % Solve system
        delta_a_and_xi = Jacobian\Res;

        % Update Chebyshev coefficients
        a = a + delta_a_and_xi(1:N);

        % Update eigenvalue
        xi = xi + delta_a_and_xi(N+1);
    
    end
    
    Y = D0Mat*a;
  
    %disp(xi) % Display eigenvalue
    %disp(a)  % Display Chebyshev coefficients
    
end
