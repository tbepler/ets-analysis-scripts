function f = plotParameterSelection( params, err, selected, xlab, ylab, ttl, fname, res, format, supress  )
%PLOTPARAMETERSELECTION Summary of this function goes here
%   Detailed explanation goes here

MARKER_SIZE = 9;
TITLE_SIZE = 24;

if nargin < 10
    supress = false;
end

if nargin < 9
    format = '-depsc';
end

if nargin < 8
    res = '-r900';
end

if nargin < 7
    fname = [];
end

if nargin < 6
    ttl = [];
end

if nargin < 5
    ylab = [];
end

if nargin < 4
    xlab = [];
end

if nargin < 3
    selected = [];
end

figargs = {};
if supress
    figargs = { 'Visible', 'off' };
end

f = figure( figargs{:} );

hold on;
axis square;

plot( params, err, '-o', 'MarkerSize', MARKER_SIZE );

if ~isempty( selected )
    plot( params(selected), err(selected), 'o' ...
        , 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', MARKER_SIZE );
end

if ~isempty( xlab )
    xlabel( xlab );
end

if ~isempty( ylab )
    ylabel( ylab );
end

if ~isempty( ttl )
    title( ttl, 'FontSize', TITLE_SIZE );
end

hold off;

if ~isempty( fname )
    print( f, fname, res, format );
end

if supress
    close(f);
end

end

