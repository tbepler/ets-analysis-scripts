function runAnalysis( ids )

IDs = { 'bayes', 'ridge', 'feature-ridge', 'pkk-ridge' };

if nargin < 1
    ids = IDs;
end

if ischar( ids )
    ids = { ids };
end

bayes = false;
ridge = false;
featureRidge = false;
pkk_ridge = false;

for i = 1 : length(ids)
    id = ids{i};
    switch id
        case 'bayes'
            bayes = true;
        case 'ridge'
            ridge = true;
        case 'feature-ridge'
            featureRidge = true;
        case 'pkk-ridge'
            pkk_ridge = true;
    end
end

args = {
    { 'ELK1_10nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/elk1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ELK1_50nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/elk1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ELK1_100nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/elk1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ELK1_uPBM_log', true, 'idr_peaks/elk1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ETS1_100nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/ets1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ETS1_10nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/ets1AndNonOverlappingDnaseLabeled.seqs' }
    { 'ETS1_uPBM_log', true, 'idr_peaks/ets1AndNonOverlappingDnaseLabeled.seqs' }
    { 'GABPA_100nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/gabpaAndNonOverlappingDnaseLabeled.seqs' }
    { 'GABPA_50nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/gabpaAndNonOverlappingDnaseLabeled.seqs' }
    { 'GABPA_25nM_boundNeg_log_combinedReps_fwd', true, 'idr_peaks/gabpaAndNonOverlappingDnaseLabeled.seqs' }
    { 'GABPA_uPBM_log', true, 'idr_peaks/gabpaAndNonOverlappingDnaseLabeled.seqs' }
    };

% args = {
%      { 'ETS1', 100, 'fwd' } ...
%      { 'GABPA', 100, 'fwd' } ...
%      { 'GABPA', 50, 'fwd' } ...
%      { 'ETS1', 10, 'fwd' } ...
%      { 'GABPA', 25, 'fwd' } ...
% };

% args = {
% %     { 'ELK1', 100, 'fwd' } ...
%      { 'ETS1', 100, 'fwd' } ...
%      { 'GABPA', 100, 'fwd' } ...
% %     { 'ELK1', 50, 'fwd' } ...
%      { 'GABPA', 50, 'fwd' } ...
% %     { 'ELK1', 10, 'fwd' } ...
%      { 'ETS1', 10, 'fwd' } ...
%      { 'GABPA', 25, 'fwd' } ...
% %     { 'ELK1', 100, 'rvs' }, ...
% %     { 'ETS1', 100, 'rvs' }, ...
% %     { 'GABPA', 100, 'rvs' }, ...
% %     { 'ELK1', 50, 'rvs' }, ...
% %     { 'GABPA', 50, 'rvs' }, ...
% %     { 'ELK1', 10, 'rvs' }, ...
% %     { 'ETS1', 10, 'rvs' }, ...
% %     { 'GABPA', 25, 'rvs' }
%     };

def_use_bias = true;
kmers = { '1mer', '1-2mer', '1-3mer' };
bpcores = { '' ...
%     '12', ...
%     '24' ...
    };
calcBetaOn = 1;
kfold = 5;
lambdas = 10.^(-5:5);
samples = 2;
sample_size = 10000;

basedir = 'models/';

time = now;
time_format = 'yyyy_mm_dd_HH.MM.SS';
timestr = datestr( time, time_format );

if bayes
    if use_bias
        dir = [ basedir 'bayesian-regression/' ];
    else
        dir = [ basedir 'bayesian-regression-nobias/' ];
    end
    fprintf( 'Starting Bayesian regression analysis...\n' );
    for i = 1 : length( args )
        for j = 1 : length( bpcores )
            core = bpcores{j};
            arg = args{i};
            odir = toDir( dir, arg{:}, core );
            buildModelsBayesianRegression( odir, arg, core, kmers, use_bias, calcBetaOn, kfold, lambdas, time );
        end
    end
end

if ridge
    %if use_bias
        dir = [ basedir 'ridge-regression/' ];
    %else
    %    dir = [ basedir 'ridge-regression-nobias/' ];
    %end
    fprintf( 'Starting Ridge regression analysis...\n' );
    for i = 1 : length( args )
        for j = 1 : length( bpcores )
            core = bpcores{j};
            arg = args{i};
            if length( arg ) >= 2
                use_bias = arg{2};
            else
                use_bias = def_use_bias;
            end
            if length( arg ) >= 3
                roc = arg{3};
            else
                roc = '';
            end
            if length( arg ) >= 4
                this_kmers = arg{4};
            else
                this_kmers = kmers;
            end
            arg = arg{1};
            odir = toDir( dir, arg, core, use_bias );
            buildModelsRidgeRegression( odir, arg, core, this_kmers, use_bias, kfold, lambdas, roc, time );
        end
    end
end

if featureRidge
    %if use_bias
        dir = [ basedir 'feature-selected-ridge/' ];
    %else
    %    dir = [ basedir 'feature-selected-ridge-nobias/' ];
    %end
    fprintf( 'Starting feature selection - ridge regression analysis...\n' );
    for i = 1 : length( args )
        for j = 1 : length( bpcores )
            core = bpcores{j};
            arg = args{i};
            use_bias = arg{2};
            arg = arg{1};
            odir = toDir( dir, arg, core, use_bias );
            buildModelsFeatureRidge( odir, arg, core, kmers, lambdas, kfold, use_bias, time );
        end
    end
end

if pkk_ridge
    dir = [ basedir 'pkk-ridge/' ];
    fprintf( 'Starting positional k-mer kernel ridge regression analysis...\n' );
    for i = 1 : length( args )
        for j = 1 : length( bpcores )
            core = bpcores{j};
            arg = args{i};
            odir = toDir( dir, arg{:}, core );
            buildModelsPKKRidge( odir, arg{:}, core, time, lambdas, kfold, samples, sample_size );
        end
    end
end


    function dir = toDir( basedir, id, bpcores, use_bias )
        dir = [ basedir toStr( id, bpcores, use_bias ) '_' timestr '/' ];
    end

    function str = toStr( id, bpcores, use_bias )
        if ~isempty( bpcores )
            bpcores = [ '_' bpcores 'bpcore' ];
        end
        str = [ id bpcores ];
        if ~use_bias
            str = [ str '_nobias' ];
        end
    end

end