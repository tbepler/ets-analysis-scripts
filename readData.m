function [Xtrain, Ytrain, Xtest, Ytest, colheaders, name] = readData( id, bp, kmer )
if ~isempty( kmer )
    kmer =['_' kmer];
end

if ~isempty( bp )
    bp = [ '_' bp 'bpcore' ];
end

% name = [ tf '_' num2str(c) 'nM_bound_negctrl_log_' o bp kmer ];
% test = [ 'test/' tf '_' num2str(c) 'nM_bound_negctrl_log_' o '_test' bp kmer '.txt' ];
% train = [ 'train/' tf '_' num2str(c) 'nM_bound_negctrl_log_' o '_train' bp kmer '.txt' ];

name = [ id bp kmer ];
test = [ 'test/' id '_test' bp kmer '.txt' ];
train = [ 'train/' id '_train' bp kmer '.txt' ];


if ~isempty( kmer )
    
    delim = ' ';
    headerlines = 1;
    
    fprintf( 'Reading file: %s\n', test );
    testData = importdata( test, delim, headerlines );
    Xtest = testData.data( : , 1 : end - 1 );
    Ytest = testData.data( : , end );
    
    fprintf( 'Reading file: %s\n', train );
    trainData = importdata( train, delim, headerlines );
    Xtrain = trainData.data( : , 1 : end - 1 );
    Ytrain = trainData.data( : , end );
    
    colheaders = trainData.colheaders(1:end-1);
    
else
    
    fprintf( 'Reading file: %s\n', test );
    file = fopen( test, 'r' );
    cols = textscan( file, '%s %f' );
    Xtest = cols{1};
    Xtest = char( Xtest );
    Ytest = cols{2};
    fclose( file );
    
    fprintf( 'Reading file: %s\n', train );
    file = fopen( train, 'r' );
    cols = textscan( file, '%s %f' );
    Xtrain = cols{1};
    Xtrain = char( Xtrain );
    Ytrain = cols{2};
    fclose( file );
    
    colheaders = [];
    
end