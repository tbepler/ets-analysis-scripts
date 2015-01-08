function report  = setup( output, id, bp, use_bias, time )
%SETUP Summary of this function goes here
%   Detailed explanation goes here

set( 0 , 'DefaultAxesFontSize', 20 );
set( 0, 'DefaultTextFontSize', 20 );
set( 0, 'DefaultLineLinewidth', 2.5 );

if ~isempty( bp )
    bp = [ '_' bp 'bpcore' ];
end

name = [ id bp ];

if ~use_bias
    name = [ name '_nobias' ];
end

reportfile = [ output name '_report.txt' ];

report = Report( reportfile, name, '', '    ', 0, '', output, time );

end

