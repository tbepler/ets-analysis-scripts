function models = buildModelsBayesianRegression(output, tf, c, o, bp, kmers, use_bias, calcBetaOn, kfold, lambdas, time )
    
    titleSize = 24;
    scatterMarker = 50;
    
    if nargin < 10
        lambdas = [];
    end
    
    if nargin < 9 || isempty( kfold )
        kfold = 5;
    end
    
    if nargin < 8 || isempty( calcBetaOn )
        calcBetaOn = 1;
    end
    
    if isempty( kmers )
        msgId = 'buildModelsBayesianRegression:kmers';
        msg = 'kmers must not be empty';
        exc = MException( msgId, msg );
        throw( exc );
    end
    
    if ischar( kmers )
        kmers = { kmers };
    end
    
    [figures, res, format, name, mTitle, report ] = setup( output, tf, c, o, bp, use_bias, time );
    
    fprintf( 'Performing Bayesian regression analysis on %s...\n', name );
    fprintf( 'Output dir = %s\n', output );
    
    %calc beta
    [ Xtrain, Ytrain, Xtest, Ytest, ~, fname ] = readData( tf, c, o, bp, kmers{calcBetaOn} );
    fprintf( 'Estimating beta from %s...\n', fname );
    f = figure( 'Visible', 'off' );
    hold on;
    axis square;
    [beta, centered] = empiricalBeta( [ Xtrain ; Xtest ], [Ytrain ; Ytest ], true );
    centeredTrain = centered( 1 : length(Ytrain) );
    centeredTest = centered( length(Ytrain)+1 : end );
    xlabel( 'Replicate Normalized Intensity' );
    ylabel( 'Counts' );
    title( [ mTitle ' Measurement Noise' ], 'FontSize', titleSize );
    hold off;
    plot_scatter = [ figures name '_replicateHistogram.eps' ];
    print( f, plot_scatter, res, format ); 
    close(f);
    
    fprintf( 'Beta = %f\n', beta );
    fprintf( report, 'Beta: %f\nHistogram: %s\n\n', beta, plot_scatter );
    
    models = cell( length( kmers ), 1 );
    
    alphas = beta * lambdas;
    
    cmax = 3*std( centered );
    %cmax = max( abs( centered ) );
    
    for i = 1 : length( kmers )
        kmer = kmers{i};
        [ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( tf, c, o, bp, kmer );
        fprintf( '\nTraining model on %s...\n', fname );
        fprintf( report, '%s model:\n', kmer );
        [ M, S, b, xvalerr, a, alphagrid, betagrid ] = bayesRegression( Xtrain, Ytrain, alphas, beta, kfold, use_bias);
        models{i} = struct( 'mean', M, 'cov', S );
        
        fprintf( report, '    Weights:\n' );
        
        colheaders = [ {'Bias'}, colheaders ] ;
        fname = [ output name '_' kmer '_weightsMean.txt' ];
        writeVector( fname, M, colheaders );

        fprintf( report, '        Mean: %s\n', fname );
        
        fname = [ output name '_' kmer '_weightsCovariance.txt' ];
        writeMatrix( fname, S );
        
        fprintf( report, '        Covariance: %s\n', fname );
        
        xlab = 'log_{10}\alpha';
        ylab = 'Mean Squared Error';
        plot_scatter = [ figures name '_' kmer '_crossValidation.eps' ];
        
        plotParameterSelection( log10(alphagrid), xvalerr, alphagrid == a, ...
            xlab, ylab, [ mTitle ' Model Selection' ], plot_scatter, res, format, true );
        
        fprintf( report, '    %d-fold cross validation:\n', kfold );
        fprintf( report, '        Alphas: ');
        fprintf( report, '%f ', alphagrid );
        fprintf( report, '\n' );
        fprintf( report, '        MSE: ');
        fprintf( report, '%f ', xvalerr );
        fprintf( report, '\n');
        fprintf( report, '        Alpha: %f\n', a );
        fprintf( report, '        Plot: %s\n', plot_scatter );
        
        fprintf( '\nTesting model on training set...\n' );
        Yh = bayesPredict( Xtrain, M, S, b );
        err = mse( Ytrain, Yh );
        r2 = corr( Ytrain, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
        Ydif = Ytrain - Yh;
        drm_ydif_r2 = corr( centeredTrain, Ydif )^2;
        fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( Ytrain, Yh, scatterMarker, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [ mTitle ' Train Actual vs Predicted Intensity' ], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' kmer '_trainingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        hold on;
        axis square;
        %colormap jet;
        caxis( [-cmax cmax] );
        scatter( Ytrain, Yh, scatterMarker, centeredTrain, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [ mTitle ' Train Actual vs Predicted Intensity' ], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        cb = colorbar;
        ylabel( cb, 'Difference from replicate mean' );
        hold off;
        plot_col_scatter = [ figures name '_' kmer '_trainingPredictColored.eps' ];
        print( f, plot_col_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( centeredTrain, Ydif, scatterMarker, '.');
        xlabel( 'Difference from replicate mean' );
        ylabel( 'Prediction error' );
        title( [ mTitle ' Train DRM vs Prediction Error' ], 'FontSize', titleSize );
        plot_drm_err_scatter = [ figures name '_' kmer '_trainingDifRepMeanVsPredErr.eps' ];
        print( f, plot_drm_err_scatter, res, format );
        close(f);
        
        fprintf( report, '    Training:\n' );
        fprintf( report, '        MSE: %f\n', err );
        fprintf( report, '        Actual vs Predicted R^2: %f\n', r2 );
        fprintf( report, '        Difference from replicate mean vs prediction error R^2: %f\n', drm_ydif_r2 );
        fprintf( report, '        Plots:\n');
        fprintf( report, '            %s\n', plot_scatter);
        fprintf( report, '            %s\n', plot_col_scatter);
        fprintf( report, '            %s\n', plot_drm_err_scatter);
        
        fprintf( '\nTesting model on testing set...\n' );
        Yh = bayesPredict( Xtest, M, S, b );
        err = mse( Ytest, Yh );
        r2 = corr( Ytest, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
        Ydif = Ytest - Yh;
        drm_ydif_r2 = corr( centeredTest, Ydif )^2;
        fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        scatter( Ytest, Yh, scatterMarker, '.' );
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [ mTitle ' Test Actual vs Predicted Intensity' ], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' kmer '_testingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        hold on;
        axis square;
        %colormap jet;
        caxis( [-cmax cmax] );
        scatter( Ytest, Yh, scatterMarker, centeredTest, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [ mTitle ' Test Actual vs Predicted Intensity' ], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        cb = colorbar;
        ylabel( cb, 'Difference from replicate mean' );
        hold off;
        plot_col_scatter = [ figures name '_' kmer '_testingPredictColored.eps' ];
        print( f, plot_col_scatter, res, format );
        close(f);
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( centeredTest, Ydif, scatterMarker, '.');
        xlabel( 'Difference from replicate mean' );
        ylabel( 'Prediction error' );
        title( [ mTitle ' Train DRM vs Prediction Error' ], 'FontSize', titleSize );
        plot_drm_err_scatter = [ figures name '_' kmer '_testingDifRepMeanVsPredErr.eps' ];
        print( f, plot_drm_err_scatter, res, format );
        close(f);
        
        fprintf( report, '    Testing:\n' );
        fprintf( report, '        MSE: %f\n', err );
        fprintf( report, '        Actual vs Predicted R^2: %f\n', r2 );
        fprintf( report, '        Difference from replicate mean vs prediction error R^2: %f\n', drm_ydif_r2 );
        fprintf( report, '        Plots:\n');
        fprintf( report, '            %s\n', plot_scatter);
        fprintf( report, '            %s\n', plot_col_scatter);
        fprintf( report, '            %s\n', plot_drm_err_scatter);
        fprintf( report, '\n' );
        
        fprintf( '\n' );
        
    end
        

    fclose( report );









end