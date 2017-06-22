function hHandle = plotSamples(samples, varargin)
% plotSamples(samples)
% plotSamples(samples, 'Name', value)
% h = plotSamples(...)
%
% Plots the MCMC samples in the variable samples as a histogram.
%
% Name/Value parameters:
%
% 'Name': ['Samples'] String to title the graph
% 'HDI':  ['true'] true/false whether to add the HDI to the graph
% 'ROPE': [''] If non-empty must be a two element vector indicating the
% beginning and end of the rope.
%

samples = samples(:);

p = inputParser;

checkString = @(x) assert(ischar(x) || isstring(x) || iscellstr(x), ...
    'Must be a string or character');
checkBoolean = @(x) assert(islogical(x) || x==0 || x==1, 'Must be boolean');
checkTwo = @(x) assert(isnumeric(x) && isvector(x) && length(x)==2, 'Must be a two element vector');

addParameter(p, 'Name','Samples', checkString);
addParameter(p, 'HDI', true, checkBoolean);
addParameter(p, 'ROPE', [], checkTwo);

parse(p, varargin{:});

doHDI = p.Results.HDI;
ROPE = p.Results.ROPE;
sampName = p.Results.Name;

dataColor = 'b';
hdiColor = 'r';
ropeColor = 'm';
fontSize = 13;

textGap = 0.01;
textTop = 1-textGap;
textBegin = 0.03;
boxTop = textTop+0.01;
boxLeft = textBegin-0.01;

h = histogram(samples, 'FaceColor', dataColor);

[dsamples, dX] = ksdensity(samples);
[~, modeIndex] = max(dsamples);
densityMode = dX(modeIndex);
tString = sprintf('Mode: %.3g', densityMode);
t = text(0,0,tString);
set(t, 'Units', 'normalized', 'Position', [textBegin textTop], ...
    'Color', dataColor, 'FontSize', fontSize, 'VerticalAlignment', 'top');
textExtent = get(t, 'Extent');
textTop = textTop - textExtent(4) - textGap;
boxString = {tString};
boxRight = textExtent(3)+0.01;

yLim = get(gca, 'YLim');
yEnd = yLim(2)*0.6;

HDI = prctile(samples, [2.5 97.5]);
if doHDI
    xVals = [HDI(1) HDI(1) nan HDI(2) HDI(2)];
    yVals = [0 yEnd nan 0 yEnd];
    line(xVals, yVals, 'Color', hdiColor, 'LineStyle', '--', 'LineWidth', 2);
    tString = sprintf('HDI: [%.3g %.3g]', HDI);
    t = text(0,0,tString);
    set(t, 'Units', 'normalized', 'Position', [textBegin textTop], ...
        'Color', hdiColor, 'FontSize', fontSize);
    
    boxString = [boxString; tString];
    textExtent = get(t, 'Extent');
    textTop = textTop - textExtent(4) - textGap;
    boxRight = max(boxRight, textExtent(3)+0.01);
end;

if ~isempty(ROPE)
    inRope = mean(samples > ROPE(1) & samples < ROPE(2));
    xVals = [ROPE(1) ROPE(1) nan ROPE(2) ROPE(2)];
    yVals = [0 yEnd nan 0 yEnd];
    line(xVals, yVals, 'Color', ropeColor, 'LineStyle', '-.', 'LineWidth', 2);
    tString = sprintf('In [%.3g %.3g]: %.3g%%', ROPE(1), ROPE(2), inRope*100);
    t = text(0,0,tString);
    set(t, 'Units', 'normalized', 'Position', [textBegin textTop], ...
        'Color', ropeColor, 'FontSize', fontSize);
    
    boxString = [boxString; tString];
    textExtent = get(t, 'Extent');
    textTop = textTop - textExtent(4) - textGap;
    boxRight = max(boxRight, textExtent(3)+0.01);

end;    

boxBottom = textTop + 0.01;

% p = patch([0 0 0 0], [0 0 0 0]);
% set(p, 'Units', 'Normalized', 'Position', [BoxLeft BotBottom BoxTop]);
title(sampName);


if nargout > 0
    hHandle = h;
end;

