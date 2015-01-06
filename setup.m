function [ figures, res, format, name, title, report ]  = setup( output, id, bp, use_bias, time )
%SETUP Summary of this function goes here
%   Detailed explanation goes here

set( 0 , 'DefaultAxesFontSize', 20 );
set( 0, 'DefaultTextFontSize', 20 );
set( 0, 'DefaultLineLinewidth', 2.5 );

res = '-r900';
format = '-depsc';

%output = 'models/';
figures = [ output 'figures/' ];

%time_format = 'yyyy_mm_dd_HH.MM.SS';
%timestr = datestr( time, time_format );

%name = [ tf '_' num2str(c) 'nM_bound_log_' o '_' timestr ];
if ~isempty( bp )
    bp = [ '_' bp 'bpcore' ];
end

name = [ id bp ];

if ~use_bias
    name = [ name '_nobias' ];
end
title = name;

if ~exist( output, 'dir' )
    mkdir( output );
end

if ~exist( figures, 'dir' )
    mkdir( figures );
end

report = [ output name '_report.txt' ];
report = fopen( report, 'w' );

fprintf( report, '%s - %s\n\n', name, datestr( time ) );

end

