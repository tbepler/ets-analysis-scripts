function models = buildModelsFeatureRidge(output, id, bp, kmers, lambdas, kfold, use_bias, time )
    
    INDENT = '    ';
    LEVEL = 0;

    titleSize = 24;
    scatterMarker = 50;
    
    if nargin < 7
        use_bias = true;
    end
    
    if nargin < 6 || isempty( kfold )
        kfold = 5;
    end
    
    if nargin < 5
        lambdas = [];
    end
    
    if isempty( kmers )
        msgId = 'buildModelsFeatureRidge:kmers';
        msg = 'kmers must not be empty';
        exc = MException( msgId, msg );
        throw( exc );
    end
    
    if ischar( kmers )
        kmers = { kmers };
    end
    
    [figures, res, format, name, mTitle, report ] = setup( output, id, bp, use_bias, time );
    
    
    fprintf( 'Performing feature selection - ridge regression analysis on %s...\n', name );
    fprintf( 'Output dir = %s\n', output );
    
    %centered = [];
    %cmax = 0;
    
    models = cell( length( kmers ), 1 );
    
    for i = 1 : length( kmers )
        kmer = kmers{i};
        [ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( id, bp, kmer );
        %Xtrain = gpuArray( Xtrain );
        %Ytrain = gpuArray( Ytrain );
        %Xtest = gpuArray( Xtest );
        %Ytest = gpuArray( Ytest );
        if use_bias
            colheaders = [ {'Bias'}, colheaders ];
        end
        D = size( Xtrain, 2 );
        
        %get the centered data if haven't already for use later
%         if isempty( centered )
%             [~, centered] = empiricalBeta( [ Xtrain ; Xtest ], [Ytrain ; Ytest ], false );
%             centeredTrain = centered( 1 : length(Ytrain) );
%             centeredTest = centered( length(Ytrain)+1 : end );
%             cmax = 3*std( centered );
%         end
        
        recordLn( '%s model:', kmer );
        DOWN_LEVEL();
        
        fprintf( '\nPerforming feature selection on %s...\n', fname );
        recordLn( 'Feature selection:' );
        DOWN_LEVEL();
        %perform feature selection on the training data
        [ I, W, lambda, ~, xvalerr ] = featureSelection( Xtrain, Ytrain, lambdas, kfold, ~use_bias );
        df = sum( I );
        recordLn( 'Original features: %d', D );
        recordLn( 'Selected features: %d', df );
        fprintf( 'Original features: %d\n', D );
        fprintf( 'Selected features: %d\n', df );
        
        fsfname = [ output name '_' kmer '_lassoWeights.txt' ];
        writeVector( fsfname, W , colheaders );
        recordLn( 'Lasso weights: %s', fsfname );
        
        xlab = 'log_{10}\lambda';
        ylab = 'Mean Squared Error';
        plot_lassoxval = [ figures name '_' kmer '_lassoCrossValidation.eps' ];
        plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
            xlab, ylab, [ mTitle ' Lasso Model Selection' ], plot_lassoxval, res, format, true );
        
        recordCrossValidation( kfold, lambdas, xvalerr, lambda, plot_lassoxval );
        UP_LEVEL();
        
        fprintf( '\nTraining model on %s...\n', fname );
        %train ridge regression model on training data subsetted by the I
        %feature selection vector
        [ Ws, lambda, xvalerr, ~ ] = ridgeRegression( Xtrain( : , I ), Ytrain, lambdas, ~use_bias, kfold );
        models{i} = struct( 'weights', Ws, 'features', I );
        
        if use_bias
            W = W(2:end);
            W(I) = Ws(2:end);
            W = [ Ws(1); W ];
        else
            W(I) = Ws;
        end
        
        fname = [ output name '_' kmer '_ridgeWeights.txt' ];
        writeVector( fname, W, colheaders );
        recordLn( 'Ridge weights: %s', fname );
        
        xlab = 'log_{10}\lambda';
        ylab = 'Mean Squared Error';
        plot_ridgexval = [ figures name '_' kmer '_ridgeCrossValidation.eps' ];
        plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
            xlab, ylab, [ mTitle ' Ridge Model Selection' ], plot_ridgexval, res, format, true );
        
        recordCrossValidation( kfold, lambdas, xvalerr, lambda, plot_ridgexval);
        
        fprintf( '\nTesting model on training set...\n' );
        Yh = predict( Xtrain( : , I), Ws );
        err = mse( Ytrain, Yh );
        r2 = corr( Ytrain, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
%         Ydif = Ytrain - Yh;
%         drm_ydif_r2 = corr( centeredTrain, Ydif )^2;
%         fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        axis square;
        scatter( Ytrain, Yh, scatterMarker, '.');
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [mTitle ' Training Predictions'], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' kmer '_trainingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
%         f = figure( 'Visible', 'off' );
%         hold on;
%         axis square;
%         %colormap jet;
%         caxis( [-cmax cmax] );
%         scatter( Ytrain, Yh, scatterMarker, centeredTrain, '.');
%         xlabel( 'Intensity(actual)' );
%         ylabel( 'Intensity(predicted)' );
%         title( [ mTitle ' Training Predictions' ], 'FontSize', titleSize );
%         xl = xlim();
%         yl = ylim();
%         mi = min( [xl(1) yl(1)] );
%         ma = max( [xl(2) yl(2)] );
%         xlim( [mi ma] );
%         ylim( [mi ma] );
%         cb = colorbar;
%         ylabel( cb, 'Difference from replicate mean' );
%         hold off;
%         plot_col_scatter = [ figures name '_' kmer '_trainingPredictColored.eps' ];
%         print( f, plot_col_scatter, res, format );
%         close(f);
        
%         f = figure( 'Visible', 'off' );
%         axis square;
%         scatter( centeredTrain, Ydif, scatterMarker, '.');
%         xlabel( 'Difference from replicate mean' );
%         ylabel( 'Prediction error' );
%         title( [ mTitle ' Training Error vs DRM' ], 'FontSize', titleSize );
%         plot_drm_err_scatter = [ figures name '_' kmer '_trainingDifRepMeanVsPredErr.eps' ];
%         print( f, plot_drm_err_scatter, res, format );
%         close(f);
        
%        recordPrediction( 'Training', err, r2, drm_ydif_r2, plot_scatter, plot_col_scatter, plot_drm_err_scatter );
        recordPrediction( 'Training', err, r2, plot_scatter );
        
        fprintf( '\nTesting model on testing set...\n' );
        Yh = predict( Xtest( : , I), Ws );
        err = mse( Ytest, Yh );
        r2 = corr( Ytest, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
%         Ydif = Ytest - Yh;
%         drm_ydif_r2 = corr( centeredTest, Ydif )^2;
%         fprintf( 'Difference from replicate mean vs prediction error R^2 = %f\n', drm_ydif_r2 );
        
        f = figure( 'Visible', 'off' );
        scatter( Ytest, Yh, scatterMarker, '.' );
        xlabel( 'Intensity(actual)' );
        ylabel( 'Intensity(predicted)' );
        title( [ mTitle ' Testing Predictions' ], 'FontSize', titleSize );
        xl = xlim();
        yl = ylim();
        mi = min( [xl(1) yl(1)] );
        ma = max( [xl(2) yl(2)] );
        xlim( [mi ma] );
        ylim( [mi ma] );
        plot_scatter = [ figures name '_' kmer '_testingPredict.eps' ];
        print( f, plot_scatter, res, format );
        close(f);
        
%         f = figure( 'Visible', 'off' );
%         hold on;
%         axis square;
%         %colormap jet;
%         caxis( [-cmax cmax] );
%         scatter( Ytest, Yh, scatterMarker, centeredTest, '.');
%         xlabel( 'Intensity(actual)' );
%         ylabel( 'Intensity(predicted)' );
%         title( [ mTitle ' Testing Predictions'], 'FontSize', titleSize );
%         xl = xlim();
%         yl = ylim();
%         mi = min( [xl(1) yl(1)] );
%         ma = max( [xl(2) yl(2)] );
%         xlim( [mi ma] );
%         ylim( [mi ma] );
%         cb = colorbar;
%         ylabel( cb, 'Difference from replicate mean' );
%         hold off;
%         plot_col_scatter = [ figures name '_' kmer '_testingPredictColored.eps' ];
%         print( f, plot_col_scatter, res, format );
%         close(f);
        
%         f = figure( 'Visible', 'off' );
%         axis square;
%         scatter( centeredTest, Ydif, scatterMarker, '.');
%         xlabel( 'Difference from replicate mean' );
%         ylabel( 'Prediction error' );
%         title( [ mTitle ' Training Error vs DRM' ], 'FontSize', titleSize );
%         plot_drm_err_scatter = [ figures name '_' kmer '_testingDifRepMeanVsPredErr.eps' ];
%         print( f, plot_drm_err_scatter, res, format );
%         close(f);
        
%         recordPrediction( 'Testing', err, r2, drm_ydif_r2, plot_scatter, plot_col_scatter, plot_drm_err_scatter );
        recordPrediction( 'Testing', err, r2, plot_scatter );
        
        UP_LEVEL();
        recordLn('');
        fprintf( '\n' );
        
    end
        

    fclose( report );


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
    
%     function recordPrediction( str, err, r2, drm_ydif_r2, r2plot, r2colplot, drm_ydif_r2_plot )
%         recordLn( '%s:', str );
%         DOWN_LEVEL();
%         recordLn( 'MSE: %f', err );
%         recordLn( 'Actual vs Predicted R^2: %f', r2 );
%         recordLn( 'Difference from replicate mean vs prediction error R^2: %f', drm_ydif_r2 );
%         recordLn( 'Plots:');
%         DOWN_LEVEL();
%         recordLn( '%s', r2plot);
%         recordLn( '%s', r2colplot);
%         recordLn( '%s', drm_ydif_r2_plot);
%         UP_LEVEL();
%         UP_LEVEL();
%     end

    function recordPrediction( str, err, r2, r2plot )
        recordLn( '%s:', str );
        DOWN_LEVEL();
        recordLn( 'MSE: %f', err );
        recordLn( 'Actual vs Predicted R^2: %f', r2 );
        recordLn( 'Plots:');
        DOWN_LEVEL();
        recordLn( '%s', r2plot);
        UP_LEVEL();
        UP_LEVEL();
    end    





end