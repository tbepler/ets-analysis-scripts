function models = buildModelsRidgeRegression(output, id, bp, kmers, use_bias, kfold, lambdas, time )
    
    titleSize = 24;
    scatterMarker = 50;
    
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
        [ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( id, bp, kmer );
        fprintf( '\nTraining model on %s...\n', fname );
        fprintf( report, '%s model:\n', kmer );
        [ M, lambda, xvalerr, ~ ] = ridgeRegression( Xtrain, Ytrain, lambdas, ~use_bias, kfold );
        models{i} = struct( 'weights', M );
        
        if use_bias
            colheaders = [ {'Bias'}, colheaders ] ;
        end
        fname = [ output name '_' kmer '_weights.txt' ];
        writeVector( fname, M, colheaders );

        fprintf( report, '    Weights: %s\n', fname );
        
        xlab = 'log_{10}\lambda';
        ylab = 'Mean Squared Error';
        plot_scatter = [ figures name '_' kmer '_crossValidation.eps' ];
        
        plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
            xlab, ylab, [ mTitle ' Model Selection' ], plot_scatter, res, format, true );
        
        fprintf( report, '    %d-fold cross validation:\n', kfold );
        fprintf( report, '        Lambdas: ');
        fprintf( report, '%f ', lambdas );
        fprintf( report, '\n' );
        fprintf( report, '        MSE: ');
        fprintf( report, '%f ', xvalerr );
        fprintf( report, '\n');
        fprintf( report, '        Lambda: %f\n', lambda );
        fprintf( report, '        Plot: %s\n', plot_scatter );
        
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
        
        fprintf( report, '    Training:\n' );
        fprintf( report, '        MSE: %f\n', err );
        fprintf( report, '        Actual vs Predicted R^2: %f\n', r2 );
        fprintf( report, '        Plots:\n');
        fprintf( report, '            %s\n', plot_scatter);
        
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
        
        fprintf( report, '    Testing:\n' );
        fprintf( report, '        MSE: %f\n', err );
        fprintf( report, '        Actual vs Predicted R^2: %f\n', r2 );
        fprintf( report, '        Plots:\n');
        fprintf( report, '            %s\n', plot_scatter);
        fprintf( report, '\n' );
        
        fprintf( '\n' );
        
    end
        

    fclose( report );









end