function writeMatrix( fname, M, colheaders, rowheaders )
%WRITEMATRIX Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    rowheaders = [];
end

if nargin < 3
    colheaders = [];
end

f = fopen( fname, 'w' );

if ~isempty( rowheaders ) && ~isempty( colheaders )
    colheaders = [ {''}, colheaders ];
end

if ~isempty( colheaders )
    for j = 1 : length( colheaders )
        fprintf( f, '%s ', colheaders{j} );
    end
    fprintf( f, '\n' );
end
for j = 1 : size( M, 1 )
    if ~isempty( rowheaders )
        fprintf( f, '%s ', rowheaders{j} );
    end
    fprintf( f, '%f ', M(j,:) );
    fprintf( f, '\n');
end
fclose( f );

end

