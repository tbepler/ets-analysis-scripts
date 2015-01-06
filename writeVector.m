function writeVector( fname, W, colheaders )
%WRITEWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

f = fopen( fname, 'w' );
if nargin > 2
    for j = 1 : length( colheaders )
        fprintf( f, '%s ', colheaders{j} );
    end
    fprintf( f, '\n' );
end
fprintf( f, '%f ', W );
fprintf( f, '\n' );
fclose( f );


end

