function models = buildModelsRidgeRegression(output, id, bp, kmers, use_bias, kfold, lambdas, roc, time )
    
    titleSize = 24;
    scatterMarker = 50;
    
    LEVEL = 0;
    INDENT = '    ';
    ID = '';
    
    if nargin < 8
        roc = '';
    end
    
    if nargin < 7
        lambdas = [];
    end
    
    if nargin < 6 || isempty( kfold )
        kfold = 5;
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
    
    [figures, res, format, name, mTitle, report ] = setup( output, id, bp, use_bias, time );
    
    fprintf( 'Performing Ridge regression analysis on %s...\n', name );
    fprintf( 'Output dir = %s\n', output );
    
    %cmax = max( abs( centered ) );
    
    for i = 1 : length( kmers )
        kmer = kmers{i};
        ID = kmer;
        [ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( id, bp, kmer );
        fprintf( '\nTraining model on %s...\n', fname );
        record( 'Model:\n' );
        DOWN_LEVEL();
        [ M, lambda, xvalerr, ~ ] = ridgeRegression( Xtrain, Ytrain, lambdas, ~use_bias, kfold );
        models{i} = struct( 'weights', M );
        
        if use_bias
            colheaders = [ {'Bias'}, colheaders ] ;
        end
        weightsfile = [ output name '_' kmer '_weights.txt' ];
        writeVector( weightsfile, M, colheaders );

        record( 'Weights: %s\n', weightsfile );
        
        xlab = 'log_{10}\lambda';
        ylab = 'Mean Squared Error';
        plot_scatter = [ figures name '_' kmer '_crossValidation.eps' ];
        
        plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
            xlab, ylab, [ mTitle ' Model Selection' ], plot_scatter, res, format, true );
        
        record( 'Cross validation:\n' );
        DOWN_LEVEL();
        ID = [ kmer '-xval' ];
        record( 'k-fold: %f\n', kfold );
        record( 'Lambdas: ');
        fprintf( report, '%f ', lambdas );
        fprintf( report, '\n' );
        record( 'MSE: ');
        fprintf( report, '%f ', xvalerr );
        fprintf( report, '\n');
        record( 'Lambda: %f\n', lambda );
        record( 'Plot: %s\n', plot_scatter );
        UP_LEVEL();
        ID = kmer;
        
        fprintf( '\nTesting model on training set...\n' );
        Yh = predict( Xtrain, M );
        err = mse( Ytrain, Yh );
        r2 = corr( Ytrain, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
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
        
        record( 'Training:\n' );
        DOWN_LEVEL();
        ID = [ kmer '-train' ];
        record( 'MSE: %f\n', err );
        record( 'r2: %f\n', r2 );
        record( 'Plot: %s\n', plot_scatter);
        UP_LEVEL();
        ID = kmer;
        
        fprintf( '\nTesting model on testing set...\n' );
        Yh = predict( Xtest, M );
        err = mse( Ytest, Yh );
        r2 = corr( Ytest, Yh )^2;
        fprintf( 'MSE = %f\n', err );
        fprintf( 'Actual vs Predicted R^2 = %f\n', r2 );
        
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
        
        record( 'Testing:\n' );
        DOWN_LEVEL();
        ID = [ kmer '-test' ];
        record( 'MSE: %f\n', err );
        record( 'r2: %f\n', r2 );
        record( 'Plot: %s\n', plot_scatter);
        UP_LEVEL();
        ID = kmer;
        
        fprintf( '\n' );
        
        [ auc, opt, rocfig ] = generateROC( roc, weightsfile, 'ridge', output );
        record( 'ROC:\n' );
        DOWN_LEVEL();
        ID = [ kmer '-roc' ];
        record( 'AUC: %f\n', auc );
        record( 'OPT: %f, %f\n', opt(1), opt(2) );
        record( 'Plot: %s\n', rocfig );
        UP_LEVEL();
        ID = kmer;
        
        UP_LEVEL();
        fprintf( report, '\n' );
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
        fprintf( report,  '[%s]', ID );
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






end