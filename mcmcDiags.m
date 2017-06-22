function [diags, figs] = mcmcDiags(p, pName)
% diags = mcmcDiags(p, pName)
% diags = mcmcDiags(p)
% [diags, figs] = mcmcDiags(...)
%
% Calculate diagnostics for the MCMC samples in p. Displays the sample
% traces, the autocorrelations, the Gelman-Rubin shrinkage, and the density
% histogram (with HDI) for each chain.
%
% p (k x N x m) k chains of N samples each. If the parameter is a vector
% parameter, there will be m elements of the vector. In this case, a
% separate figure is displayed for each element.
%
% pName : a character array or string giving the parameter name
%
% Output:
% 
% diags (optional): returns a struct with the ESS, MCSE, and final
% Shrinkage
% figs (optional): returns the figure handles for all the ifgures
%

if ~exist('pName') || isempty(pName)
    pName = 'Param';
end;

[numChains, numSamples, numElements] = size(p);
figHandles = zeros(1, numElements);
for elementNum = 1:numElements
    thisP = p(:,:,elementNum);

    figHandles(elementNum) = figure('Units', 'Pixels', 'Position', [403    85   731   581]);
    
    axes('Position',[0.13 0.55 0.33 0.30]);
    plot(thisP');
    xlabel('Iterations');
    ylabel('Param. Value');
    yLims = prctile(thisP(:), [1 99]);
    ylim(yLims);
    
    axes('Position',[0.57 0.55 0.33 0.30]);
    
    xc1 = my_acf(thisP');
    maxLag = find(max(xc1,[],2) < 0.05, 1);
    maxLag = max(maxLag, 20);
    xc = my_acf(thisP', maxLag);
    xc0 = [ones(1, size(thisP,1)); xc];
    plot(0:maxLag, xc0);
    line([0 maxLag], [0 0], 'LineStyle', '--', 'Color', 'k');
    ylabel('Autocorrelation');
    xlabel('Lag');
    maxAcf = prctile(xc0(:), 95);
    ylim([-0.05 maxAcf*1.1]);
    chainESS = numSamples./(1+2*sum(xc));
    ESS = sum(chainESS);
    text(maxLag-2, maxAcf*0.9, ['ESS: ' sprintf('%.0f',ESS)], 'HorizontalAlignment', 'right');

    if numChains > 1
        axes('Position',[0.13 0.11 0.33 0.30]);
        R_stat = zeros(numSamples,1);
        for thisS = 1:numSamples
            n = thisS;
            m = numChains;
            B = n*var(mean(thisP(:,1:n)'));
            W = mean(var(thisP(:,1:n)'));
            sigma2 = ((n - 1)/n) * W + (1/n) * B;
            R_stat(thisS) = sqrt((m + 1)/m * sigma2 ./ W - (n-1)/m/n);
        end;
        plot(R_stat);
        yMax = prctile(R_stat, 90);
        ylim([0.99 yMax+(yMax-1)*1.1]);
        ylabel('Shrink factor');
        xlabel('Num iterations');
    end;
    
    axes('Position',[0.57 0.11 0.33 0.30]);
    hHandle = zeros(numChains, 1);
    colors = get(gca, 'ColorOrder');
    for histNum = 1:numChains
        hHandle(histNum) = histogram(thisP(histNum,:), 'DisplayStyle', 'stairs', 'EdgeColor', colors(histNum,:));
        hold on;
    end;
    hold off;
    xlabel('Value');
    ylabel('Counts');
    MCSE = std(thisP(:))./ESS;
    t = text(0,0, ['MCSE: ' sprintf('%.4f', MCSE)]);
    set(t, 'Units', 'normalized', 'Position', [0.95 0.95], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    
    HDI = prctile(thisP', [2.5 97.5]);
    yLim = get(gca, 'YLim');
    yEnd = yLim(2)/5;
    for HDInum = 1:numChains
        thisH = HDI(:, HDInum);
        xVals = [thisH(1) thisH(1) nan thisH(2) thisH(2)];
        yVals = [0 yEnd nan 0 yEnd];
        line(xVals, yVals, 'Color', get(hHandle(HDInum), 'EdgeColor'));
    end;
    text(mean(HDI(:)), yEnd/2, 'HDI', 'FontSize', 8, 'HorizontalAlignment', 'center');
    
    thisName = pName;
    if numElements > 1
        thisName = [thisName '(' num2str(elementNum) ')'];
    end;
    annotation('textbox', [0.3 0.92 0.4 0.047], 'String',thisName, ...
        'Interpreter', 'none', ...
        'FontSize', get(gca, 'FontSize')+4, 'FontWeight', 'bold', ...
        'LineStyle', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');


    if nargout > 0
        diags(elementNum).ESS = ESS;
        diags(elementNum).MCSE = MCSE;
        diags(elementNum).HDI = HDI;
        if numChains > 1
            diags(elementNum).Shrinkage = R_stat(end);
        end;
    end;


end;


if nargout > 1
    figs = figHandles;
end;