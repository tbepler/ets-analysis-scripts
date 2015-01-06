function [AUC,OPT] = plotROC( output, dataf )
%PLOTROC Summary of this function goes here
%   Detailed explanation goes here

set( 0 , 'DefaultAxesFontSize', 20 );
set( 0, 'DefaultTextFontSize', 20 );
set( 0, 'DefaultLineLinewidth', 2.5 );

res = '-r900';
format = '-depsc';

%output = 'models/';

[ pathstr, name, ~ ] = fileparts( dataf );
if isempty( output )
    output = [ pathstr '/' ];
end

figures = [ output 'figures/' ];

if ~exist( figures, 'dir' )
    mkdir( figures );
end

outname = [ figures name '_roccurve' ];

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

f = figure();
hold on
plot( X, Y );
xlabel( 'False Positives' );
ylabel( 'True Positives' );
text( 0.65, 0.1, [ 'AUC = ' num2str(AUC) ], 'HorizontalAlignment', 'left' );
hold off

print( f, outname, res, format );

end

