function [AUC,OPT, outname ] = plotROC( output, dataf, color, plottext, display, name )
%PLOTROC Summary of this function goes here
%   Detailed explanation goes here

set( 0 , 'DefaultAxesFontSize', 20 );
set( 0, 'DefaultTextFontSize', 20 );
set( 0, 'DefaultLineLinewidth', 2.5 );

res = '-r900';
format = '-depsc';

if nargin < 3
    color = {};
end

if nargin < 6
    name = '';
end

if nargin < 5
    display = true;
end

if nargin < 4
    plottext = false;
end

%output = 'models/';

if isempty( name )
    if iscell( dataf )
        [ pathstr, name, ~ ] = fileparts( dataf{1} );
        for i = 2 : length( dataf )
            [ ~, n, ~ ] = fileparts( dataf{i} );
            name = [ name '_and_' n ];
        end
    else
        [ pathstr, name, ~ ] = fileparts( dataf );
    end
    name = [ name '_roccurve.eps' ];
end
if isempty( output )
    output = [ pathstr '/' ];
end

%figures = [ output 'figures/' ];
figures = output;

if ~exist( figures, 'dir' )
    mkdir( figures );
end

outname = [ figures name ];

figargs = {};
if ~display
    figargs = { 'Visible', 'off' };
end

f = figure( figargs{:} );
hold on;

if iscell( dataf )
    AUC = cell( length( dataf ) );
    OPT = cell( length( dataf ) );
    plotargs = {};
    for i = 1 : length(dataf)
        [a,o,x,y] = readAndROC( dataf{i} );
        AUC{i} = a;
        OPT{i} = o;
        plotargs = [ plotargs, {x}, {y} ];
        if ~isempty( color )
            plotargs = [ plotargs, { color{i} } ];
        end
    end
else
    [AUC,OPT, x, y] = readAndROC( dataf );
    plotargs = { x, y };
    if ~isempty( color )
        plotargs{3} = color;
    end
end

plot( plotargs{:} );
xlabel( 'False Positives' );
ylabel( 'True Positives' );
if plottext
    text( 0.65, 0.1, [ 'AUC = ' num2str(AUC) ], 'HorizontalAlignment', 'left' );
end

hold off

print( f, outname, res, format );

if ~display
    close( f );
end

    function [AUC,OPT, X, Y] = readAndROC( dataf )
        data = importdata( dataf );
        labels = data( :, 1);
        predictions = data(:,2);
        
        [X,Y,T,AUC,OPT,~,~] = perfcurve( labels, predictions, 1 );
        if length(X) > 3000
            I = randsample( 1:length(X), 3000 );
            I = sort(I);
            X = X(I);
            Y = Y(I);
        end
        
        
    end

end

