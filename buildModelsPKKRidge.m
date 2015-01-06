function models = buildModelsPKKRidge(output, tf, c, o, core, time, lambdas, kfold, samples, sample_size)

INDENT = '    ';
LEVEL = 0;

titleSize = 24;
scatterMarker = 50;

if nargin < 9
    sample_size = 0;
end

if nargin < 8
    samples = 0;
end

if samples < 0
    samples = 0;
end

assert( samples == 0 || sample_size > 0 );

if nargin < 7 || isempty( kfold )
    kfold = 5;
end

if nargin < 6
    lambdas = [];
end

[figures, res, format, name, mTitle, report ] = setup( output, tf, c, o, core, time );


fprintf( 'Performing positional k-mer kernel ridge regression analysis on %s...\n', name );
fprintf( 'Output dir = %s\n', output );

[ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( tf, c, o, core, [] );
%Xtest = Xtest( 1:1000, : );
%Ytest = Ytest( 1:1000, : );
[~, centered] = empiricalBeta( [ Xtrain ; Xtest ], [Ytrain ; Ytest ], false );
centered = gather( centered );
centeredTrain = centered( 1 : length(Ytrain) );
centeredTest = centered( length(Ytrain)+1 : end );
cmax = 3*std( centered );

if samples > 0
    
    fprintf( '\nBuilding models on %d subsamples of size %d of %s...\n', samples, sample_size, fname );
    recordLn( 'Samples: %d', samples );
    recordLn( 'Sample size: %d', sample_size );
    
    d = size( Xtrain, 1 );
    models = cell( samples, 1 );
%     errs = zeros( samples, 1 );
%     hold_out = randsample( 1:d, sample_size );
%     Xval = Xtrain( hold_out );
%     Yval = Ytrain( hold_out );
%     Xtrain( hold_out ) = [];
%     Ytrain( hold_out ) = [];
%     d = d - sample_size;
    
    for i = 1 : samples
        trainI = randsample( 1:d, sample_size );
        X = Xtrain( trainI, : );
        Y = Ytrain( trainI );
        fprintf( 'Training model: %d\n', i );
        recordLn( 'Sample %d:', i );
        DOWN_LEVEL();
        models{i} = buildModel( X, Y, [ 'rep' num2str(i) ], centeredTrain( trainI ) );
        UP_LEVEL();
    end
    
    
else
    
    fprintf( '\nBuilding model on data %s...\n', fname );
    models = { buildModel( Xtrain, Ytrain, [] ) };
    
end

fclose( report );

    function M = buildModel( Xtrain, Ytrain, rep, centeredTrain )
        
        [ M, lambda, xvalerr, ~ ] = ridgeRegressionKernel( Xtrain, Ytrain, lambdas, @posKmerKernel, [], false, kfold );

        xlab = 'log_{10}\lambda';
        ylab = 'Mean Squared Error';
        plot_ridgexval = [ figures name '_' rep '_pkkRidgeCrossValidation.eps' ];
        plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
            xlab, ylab, [ mTitle ' PKK Ridge Model Selection' ], plot_ridgexval, res, format, true );
        
        recordCrossValidation( kfold, lambdas, xvalerr, lambda, plot_ridgexval);
        
        fprintf( '\nTesting model on training set...\n' );
        Yh = predictKernel( Xtrain, M );
        err = mse( Ytrain, Yh );
        r2 = corr( Ytrain, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
        Ydif = Ytrain - Yh;
        drm_ydif_r2 = corr( centeredTrain, Ydif )^2;
        fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( gather(Ytrain), gather(Yh), scatterMarker, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [mTitle ' Train Actual vs Predicted Intensity'], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' rep '_trainingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        hold on;
        axis square;
        %colormap jet;
        caxis( [-cmax cmax] );
        scatter( gather(Ytrain), gather(Yh), scatterMarker, centeredTrain, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [mTitle ' Train Actual vs Predicted Intensity'], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        cb = colorbar;
        ylabel( cb, 'Difference from replicate mean' );
        hold off;
        plot_col_scatter = [ figures name '_' rep '_trainingPredictColored.eps' ];
        print( f, plot_col_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( gather( centeredTrain ), gather( Ydif ), scatterMarker, '.');
        xlabel( 'Difference from replicate mean' );
        ylabel( 'Prediction error' );
        title( [ mTitle ' Train DRM vs Prediction Error' ], 'FontSize', titleSize );
        plot_drm_err_scatter = [ figures name '_' rep '_trainingDifRepMeanVsPredErr.eps' ];
        print( f, plot_drm_err_scatter, res, format );
        close(f);
        
        recordPrediction( 'Training', err, r2, drm_ydif_r2, plot_scatter, plot_col_scatter, plot_drm_err_scatter );
        
        fprintf( '\nTesting model on testing set...\n' );
        Yh = predictKernel( Xtest, M );
        err = mse( Ytest, Yh );
        r2 = corr( Ytest, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
        Ydif = Ytest - Yh;
        drm_ydif_r2 = corr( centeredTest, Ydif )^2;
        fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        scatter( gather( Ytest ), gather( Yh ), scatterMarker, '.' );
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [mTitle ' Test Actual vs Predicted Intensity'], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' rep '_testingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        hold on;
        axis square;
        %colormap jet;
        caxis( [-cmax cmax] );
        scatter( gather( Ytest ), gather( Yh ), scatterMarker, centeredTest, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [mTitle ' Test Actual vs Predicted Intensity'], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        cb = colorbar;
        ylabel( cb, 'Difference from replicate mean' );
        hold off;
        plot_col_scatter = [ figures name '_' rep '_testingPredictColored.eps' ];
        print( f, plot_col_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( gather( centeredTest) , gather( Ydif ), scatterMarker, '.');
        xlabel( 'Difference from replicate mean' );
        ylabel( 'Prediction error' );
        title( [ mTitle ' Test DRM vs Prediction Error' ], 'FontSize', titleSize );
        plot_drm_err_scatter = [ figures name '_' rep '_testingDifRepMeanVsPredErr.eps' ];
        print( f, plot_drm_err_scatter, res, format );
        close(f);
        
        recordPrediction( 'Testing', err, r2, drm_ydif_r2, plot_scatter, plot_col_scatter, plot_drm_err_scatter );
        
        UP_LEVEL();
        recordLn('');
        fprintf( '\n' );
    end


    function DOWN_LEVEL()
        LEVEL = LEVEL + 1;
    end

    function UP_LEVEL()
        LEVEL = LEVEL - 1;
        if LEVEL < 0
            LEVEL = 0;
        end
    end

    function record( str, varargin )
        indt = LEVEL;
        while indt > 0
            fprintf( report, INDENT );
            indt = indt - 1;
        end
        fprintf( report, str, varargin{:} );
    end

    function recordLn( str, varargin )
        record( str, varargin{:} );
        fprintf( report, '\n' );
    end

    function recordCrossValidation( kfold, lambdas, xvalerr, lambda, plotf )
        recordLn( '%d-fold cross validation:', kfold );
        DOWN_LEVEL();
        
        record( 'Lambdas: ' );
        fprintf( report, '%f ', lambdas );
        fprintf( report, '\n' );
        
        record( 'MSE: ' );
        fprintf( report, '%f ', xvalerr );
        fprintf( report, '\n');
        
        recordLn( 'Lambda: %f', lambda );
        recordLn( 'Plot: %s', plotf );
        
        UP_LEVEL();
    end

    function recordPrediction( str, err, r2, drm_ydif_r2, r2plot, r2colplot, drm_ydif_r2_plot )
        recordLn( '%s:', str );
        DOWN_LEVEL();
        recordLn( 'MSE: %f', err );
        recordLn( 'Actual vs Predicted R^2: %f', r2 );
        recordLn( 'Difference from replicate mean vs prediction error R^2: %f', drm_ydif_r2 );
        recordLn( 'Plots:');
        DOWN_LEVEL();
        recordLn( '%s', r2plot);
        recordLn( '%s', r2colplot);
        recordLn( '%s', drm_ydif_r2_plot);
        UP_LEVEL();
        UP_LEVEL();
    end





end