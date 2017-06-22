function ta = my_acf(y,maxLag,doGraph)
% ACF - Compute Autocorrelations Through p Lags
% myacf = my_acf(y,maxLag) 
% myacf = my)_acf(y,maxLag,true) % for producing the graph
%
% Inputs:
% y - matrix to compute acf for, nxk, operates on columns
% maxLag - total number of lags, 1x1 integer, optional
%
% Output:
% myacf - maxLag x k matrix containing autocorrelations
%        (First lag computed is lag 1. Lag 0 not computed)
%
%
% Optionally: A bar graph of the autocorrelations is also produced, with
% rejection region bands for testing individual autocorrelations = 0.
%
% Note that lag 0 autocorelation is not computed, 
% and is not shown on this graph.
%
% Example:
% >> acf(randn(100,1), 10)
%


% --------------------------
% USER INPUT CHECKS
% --------------------------

[N, k] = size(y) ;

if ~exist('doGraph', 'var') && exist('maxLag', 'var') && islogical(maxLag)
    doGraph = maxLag;
    maxLag = [];
end;
if ~exist('doGraph', 'var') || isempty(doGraph)
    doGraph = false;
end;

if ~exist('maxLag', 'var') || isempty(maxLag)
    maxLag = N-1;
end;
[a1, a2] = size(maxLag) ;
if ~((a1==1 && a2==1) && (maxLag<N))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end



% -------------
% BEGIN CODE
% -------------

% Collect ACFs at each lag i
yBar = mean(y);
yZero = y - repmat(yBar, N, 1);
yvar = zeros(1,k);
for k1 = 1:k
    yvar(k1) = yZero(:,k1)'*yZero(:,k1);
end;

ta = zeros(maxLag,k) ;
cross_sum = zeros(1,k);
for lag = 1:maxLag
    for k1 = 1:k
        cross_sum(k1) = yZero(1:(N-lag),k1)'*yZero(lag+1:N,k1);
    end;
    ta(lag,:) = cross_sum ./ yvar ;
end;

if doGraph
    % Plot ACF
    % Plot rejection region lines for test of individual autocorrelations
    % H_0: rho(tau) = 0 at alpha=.05
    for k1 = 1:k
        ta1 = ta(:,k1);
        bar(ta1)
        line([0 maxLag+.5], (1.96)*(1/sqrt(N))*ones(1,2))
        line([0 maxLag+.5], (-1.96)*(1/sqrt(N))*ones(1,2))
    
        % Some figure properties
        line_hi = (1.96)*(1/sqrt(N))+.05;
        line_lo = -(1.96)*(1/sqrt(N))-.05;
        bar_hi = max(ta1)+.05 ;
        bar_lo = -max(ta1)-.05 ;
    
        if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
            axis([0 maxLag+.60 line_lo line_hi])
        else
            axis([0 maxLag+.60 bar_lo bar_hi])
        end
        title({['Vector ' num2str(k1)],'Sample Autocorrelations',' '})
        xlabel('Lag Length')
        set(gca,'YTick',-1:.20:1)
        % set number of lag labels shown
        if (maxLag<28 && maxLag>4)
            set(gca,'XTick',floor(linspace(1,maxLag,4)))
        elseif (maxLag>=28)
            set(gca,'XTick',floor(linspace(1,maxLag,8)))
        end
        set(gca,'TickLength',[0 0])
    end;
end;

