function [auc, opt, figname] = generateROC( seqs, weightsfile, modeltype, outdir, display )

if isempty( seqs )
    return;
end

if nargin < 5
    display = false;
end

figdir = [ outdir 'figures/' ];

%first score seqs with the weights
[ ~, seqsname, ~ ] = fileparts( seqs );
[ ~, weightsname, ~ ] = fileparts( weightsfile );
weightsname = regexprep( weightsname, '_[^_]*$', '' );
scoresname = [ outdir seqsname '_' weightsname '_' modeltype 'Scored.txt' ];

if display
    fprintf( 'Performing ROC analysis on %s...\n', seqsname );
end

cmd = [ './scoreSeq.out ' weightsfile ' ' seqs ' > ' scoresname ];
if display
    fprintf( [ 'Running: ' cmd '\n' ] );
end
system( cmd );

if display
    fprintf( 'Plotting ROC curve\n' );
end
[auc, opt, figname] = plotROC( figdir, scoresname, false, false );


end
