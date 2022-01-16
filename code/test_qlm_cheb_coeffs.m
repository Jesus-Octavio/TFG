%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% QUASI - LINEARIZATION METHOD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.
% Absolute values of Chebyshev coefficientes for ten different values of
% the polytropic index n. The coefficientes for m=n, m=n+1/2 (same colorm
% dashed) are similar for small degree but then differ framatically for
% large degree. When n is an integer, there is no singularity and the
% exponential decay continues to all degrees. When n is not an integer, the
% coefficients have a 'long tail of decay as an inverse power law
% proportional to n^(-(2m+5))

% (Figure 4.3)

format longg

% Color palette
myColor = [1 0 1; 0 0 1; 0 1 0; 0.9290, 0.6940, 0.1250; 1 0 0; 0 0 0];

N = 100;                                      % Collocation points
iters = 30;                                   % Number of Newton iterations


for n = 0:0.5:4.5
    [xi, Y, XCheb, a] = qlm(N,n,iters);       % Function call to QSLM
    
    
    % If integer
    if floor(n) == n                              
        % Plot absolute values of Chebyshev coefficientes
        semilogy(1:N, abs(a),'-','LineWidth', 1,'Color',myColor(floor(n)+1,:))                  
    
    % If not integer (dashed)
    else
        % Plot absolute values of Chebyshev coefficientes
        semilogy(1:N, abs(a),'--','LineWidth', 1.5,...,
            'Color',myColor(floor(n)+1,:))               
    end
    hold on
end

ylim([10e-19 10e0]);
xticks([0 20 40 60 80 100]);
xlabel('Degree'), ylabel('Chebyshev coefficients.'); 
title('Absolute values of Chebyshev coeffficients for different values of the index n.');
legend('n=0','n=0.5','n=1','n=1.5','n=2','n=2.5','n=3','n=3.5','n=4',...,
    'n=4.5','Location','NorthEast','NumColumns',5)
