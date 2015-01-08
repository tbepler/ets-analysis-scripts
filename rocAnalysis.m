function rocAnalysis( )

COLORS = { 'b', 'r', 'g', 'k', 'c', 'y' };

time_format = 'yyyy_mm_dd_HH.MM.SS';
timestr = datestr( now, time_format );

peaksdir = 'idr_peaks/';
rocdir = [ 'roc/' timestr '/' ];
peaks = {
    'wgEncodeSydhTfbsK562Elk112771IggrabAlnPoolRep1-2.MACS2-IDR0.05_peaks_0_and_wgEncodeHaibTfbsK562Ets1V0416101AlnPoolRep1-2.MACS2-IDR0.05_peaks_1.seqs'
};

models = {
    'models/ridge-regression/ELK1_100nM_boundNeg_log_combinedReps_fwd_2015_01_08_00.23.14/ELK1_100nM_boundNeg_log_combinedReps_fwd_1mer_weights.txt'
    'models/ridge-regression/ELK1_100nM_boundNeg_log_combinedReps_fwd_2015_01_08_00.23.14/ELK1_100nM_boundNeg_log_combinedReps_fwd_1-2mer_weights.txt'
    'models/ridge-regression/ELK1_100nM_boundNeg_log_combinedReps_fwd_2015_01_08_00.23.14/ELK1_100nM_boundNeg_log_combinedReps_fwd_1-3mer_weights.txt'
    'kmer_pwm/ELK1_pwm.txt'
};

suffix = 'ELK1_allROCs';

reverse_labels = {
    true
};

aucmatrix = zeros( length( models ), length( peaks ) );

for i = 1 : length( peaks )
    peak = [ peaksdir peaks{i} ];
    fprintf( 'Running ROC analysis on seqs: %s\n\n', peak );
    [ ~, peak_name, ~ ] = fileparts( peak );
    outdir = [ rocdir peak_name '_' suffix '/' ];
    if ~exist( outdir, 'dir' )
        mkdir( outdir );
    end
    figdir = [ outdir 'figures/' ];
    if ~exist( figdir, 'dir' )
        mkdir( figdir );
    end
    report = fopen( [ outdir peak_name '_' suffix '_report.txt' ], 'w' );
    scorefiles = cell( length( models ) );
    for j = 1 : length( models )
        model = models{ j };
        [ ~, model_name, ~ ] = fileparts( model );
        output = [ outdir peak_name '_' model_name 'Scored.txt' ];
        if ~exist( output, 'file' )
            if isempty( strfind( model, 'pwm' ) )
                cmd = [ './scoreSeq.out ' model ' ' peak ' > ' output ];
            else
                cmd = [ 'python pwmScoreSeq.py ' model ' ' peak ' > ' output ];
            end
            fprintf( [ 'Running: ' cmd '\n' ] );
            system( cmd );
        else
            fprintf( '%s already exists.\n', output );
        end
        scorefiles{ j } = output;
    end
    [ aucs, opts, figfile ] = plotROC( figdir, scorefiles, COLORS( 1 : length( scorefiles ) ), false, false, [ peak_name '_' suffix '.eps' ], reverse_labels{i} );
    aucmatrix(:,i)=aucs;
    fprintf( report, 'Peaks: %s\n', peak_name );
    for j = 1 : length( models )
        [ ~, model_name, ~ ] = fileparts( models{j} );
        fprintf( report, 'Model: %s\n', model_name );
        fprintf( '\n%s\n', model_name );
        fprintf( report, '    AUC: %f\n', aucs(j) );
        fprintf( 'AUC: %f\n', aucs(j) );
        fprintf( report, '    OPT: %f, %f\n', opts(j,1), opts(j,2) );
        fprintf( 'OPT: %f, %f\n', opts(j,1), opts(j,2) );
        fprintf( report, '    Peak scores: %s\n', scorefiles{j} );
        fprintf( report, '    Color: %s\n', COLORS{j} );
    end
    fprintf( report, 'Plot: %s\n', figfile );
    fprintf( '\nPlot: %s\n\n', figfile );
    
    fclose( report );
end
    
peak_names = cell( length( peaks ) );
for i = 1 : length( peaks );
    [ ~, name, ~ ] = fileparts( peaks{i} );
    peak_names{i} = name;
end
model_names = cell( length( models ) );
for i = 1 : length( models );
    [ ~, name, ~ ] = fileparts( models{i} );
    model_names{i} = name;
end
writeMatrix( [ rocdir 'auc_matrix.txt' ], aucmatrix, peak_names, model_names );


end