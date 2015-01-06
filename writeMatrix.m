function writeMatrix( fname, M, colheaders )
%WRITEMATRIX Summary of this function goes here
%   Detailed explanation goes here

f = fopen( fname, 'w' );
if nargin > 2
    for j = 1 : length( colheaders )
        fprintf( f, '%s ', colheaders{j} );
    end
    fprintf( f, '\n' );
end
for j = 1 : size( M, 1 )
    fprintf( f, '%f ', M(j,:) );
    fprintf( f, '\n');
end
fclose( f );

end

