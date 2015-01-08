function models = buildModelsRidgeRegression(output, id, bp, kmers, use_bias, kfold, lambdas, roc, time )
    
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
    
    report = setup( output, id, bp, use_bias, time );
    
    fprintf( 'Performing Ridge regression analysis on %s...\n', report.name );
    fprintf( 'Output dir = %s\n', report.outdir );
    
    model_peak_scores = cell( length( kmers ) );
    names = cell( length( kmers ) );
    
    %cmax = max( abs( centered ) );
    
    for i = 1 : length( kmers )
        kmer = kmers{i};
        if isempty( kmer )
            kmer = 'pkk';
        end
        names{i} = [ report.name '_' kmer ];
        report.setKmer( kmer );
        if strcomp( kmer, 'pkk' )
            if ~isempty( bp )
                bp = [ '_' bp 'bpcore' ];
            end
            fname = [ id bp ];
            test = [ 'test/' id '_test' bp kmer '.txt' ];
            train = [ 'train/' id '_train' bp kmer '.txt' ];
            
            fprintf( '\nTraining model on %s...\n', fname );
            report.printf( 'Model:\n' );
            report.down_level();
            
            %determine lambda using cross validation on the 1-3mer feature
            %set
            [ Xtrain, Ytrain, ~, ~, ~, ~ ] = readData( id, bp, '1-3mer' );
            train = @( x, y, l ) ridgeRegression( x, y, l, false );
            [ lambda, ~ ] = crossValidate( Xtrain, Ytrain, kfold, lambdas, train, @predict, @mse );
            report.printAndDisplayf( true, 'Lambda: %f\n', lambda );
            Xtrain = [];
            Ytrain = [];
            
            modelfile = [ report.outdir report.name '_pkkRidge.txt' ];
            time_format = 'yyyy_mm_dd_HH.MM.SS';
            timestr = datestr( time, time_format );
            predfile = [ '.tmp.preds.' timestr ];
            cmd = [ 'pkkridge train -d ' train ' -l ' num2str( lambda ) ' -m ' modelfile ' -p ' predfile ];
            
                
        else
            [ Xtrain, Ytrain, Xtest, Ytest, colheaders, fname ] = readData( id, bp, kmer );
            fprintf( '\nTraining model on %s...\n', fname );
            report.printf( 'Model:\n' );
            report.down_level();
            [ M, lambda, xvalerr, ~ ] = ridgeRegression( Xtrain, Ytrain, lambdas, ~use_bias, kfold );
            models{i} = struct( 'weights', M );
            
            if use_bias
                colheaders = [ {'Bias'}, colheaders ] ;
            end
            weightsfile = report.weights( M, colheaders, '' );
            
            report.crossValidation( kfold, lambdas, xvalerr, lambda, true );
            
            fprintf( '\nTesting model on training set...\n' );
            Yh = predict( Xtrain, M );
            report.predictions( true, Ytrain, Yh, 'train' );
            
            fprintf( '\nTesting model on testing set...\n' );
            Yh = predict( Xtest, M );
            report.predictions( true, Ytest, Yh, 'test' );
            
            fprintf( '\n' );
            
            model_peak_scores{i} = report.roc( roc, weightsfile, 'ridge', true );
            
            report.up_level();
            report.setId( '' );
            report.printf( '\n' );
            fprintf( '\n' );
            
        end
        
    end
    
    if ~isempty( roc )
        report.setId( 'roc-all' );
        report.overlayROCs( model_peak_scores, names );
    end
        

    report.close();




end