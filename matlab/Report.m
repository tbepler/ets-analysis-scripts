classdef Report < handle
    
    properties( Constant = true )
        TITLE_SIZE = 24;
        SCATTER_MARKER = 50;
        RES = '-r900';
        FORMAT = '-depsc';
        COLORS = { 'b', 'r', 'g', 'k', 'c', 'y' };
    end
    
    properties
        file;
        indent;
        level;
        id;
        name;
        kmer;
        titlename;
        outdir;
        figdir;
    end
    
    methods
        
        function obj = Report( filepath, name, kmer, indent, level, id, outdir, time )
            
            if nargin < 8
                time = now;
            end
            
            if nargin < 7
                [ outdir, ~, ~ ] = fileparts( filepath );
                outdir = [ outdir '/' ];
            end
            
            if nargin < 6
                id = '';
            end
            
            if nargin < 5
                level = 0;
            end
            
            if nargin < 4
                indent = '    ';
            end
            
            obj.outdir = outdir;
            obj.figdir = [ outdir 'figures/' ];
            
            if ~exist( obj.outdir, 'dir' )
                mkdir( obj.outdir );
            end
            
            if ~exist( obj.figdir, 'dir' )
                mkdir( obj.figdir );
            end
            
            obj.file = fopen( filepath, 'w' );
            obj.indent = indent;
            obj.level = level;
            obj.id = id;
            obj.name = name;
            obj.kmer = kmer;
            obj.titlename = [ obj.name '_' obj.kmer ];
            
            obj.printf( '%s - %s\n\n', name, datestr( time ) );
            
        end
        
        function self = down_level( self )
            self.level = self.level + 1;
        end
        
        function self = up_level( self )
            self.level = self.level - 1;
            if self.level < 0
                self.level = 0;
            end
        end
        
        function self = printf( self, str, varargin )
            indt = self.level;
            if ~isempty( self.id )
                fprintf( self.file,  '[%s]', self.id );
            end
            while indt > 0
                fprintf( self.file, self.indent );
                indt = indt - 1;
            end
            fprintf( self.file, str, varargin{:} );
        end
        
        function self = printAndDisplayf( self, display, str, varargin )
            self.printf( str, varargin{:} );
            if display
                fprintf( str, varargin{:} );
            end
        end
        
        function self = setId( self, id )
            self.id = id;
        end
        
        function self = setKmer( self, kmer )
            self.kmer = kmer;
            self.id = kmer;
        end
        
        function scoresfile = roc( self, seqspath, weightspath, mtype, display )
            if nargin < 5
                display = false;
            end
            [ auc, opt, rocfig, scoresfile ] = generateROC( seqspath, weightspath, mtype, self.outdir, display );
            self.printf( 'ROC:\n' );
            self.down_level();
            previd = self.id;
            self.setId( [ previd '-roc' ] );
            self.printAndDisplayf( display, 'AUC: %f\n', auc );
            self.printAndDisplayf( display, 'OPT: %f, %f\n', opt(1), opt(2) );
            self.printf( 'Plot: %s\n', rocfig );
            self.up_level();
            self.setId( previd );
        end
        
        function self = overlayROCs( self, scorefiles, names )
            [~,~,outname] = plotROC( self.figdir, scorefiles, self.COLORS( 1:length(scorefiles) ), false, false, [ self.name '_rocAll.eps' ] );
            self.printf( 'ROC-all:\n' );
            self.down_level();
            self.printf( 'Plot: %s\n', outname );
            for i = 1 : length( scorefiles )
                self.printf( '%s: %s\n', self.COLORS{i}, names{i} );
            end
            self.up_level();
        end
            
            
        
        function self = predictions( self, display, Yactual, Ypred, id )
            
            err = mse( Yactual, Ypred );
            r2 = corr( Yactual, Ypred )^2;
            
            plotfile = self.plotScatter( Ypred, Yactual, id );
          
            self.printf( '%s:\n', id );
            self.down_level();
            previd = self.id;
            self.setId( [ previd '-' id ] );
            self.printAndDisplayf( display, 'MSE: %f\n', err );
            self.printAndDisplayf( display, 'r2: %f\n', r2 );
            self.printf( 'Plot: %s\n', plotfile);
            self.up_level();
            self.setId( previd );
        end
        
        function plotfile = plotScatter( self, Ypred, Yactual, id )
            f = figure( 'Visible', 'off' );
            scatter( Yactual, Ypred, self.SCATTER_MARKER, '.' );
            xlabel( 'Intensity(actual)' );
            ylabel( 'Intensity(predicted)' );
            title( [ self.titlename ' ' id ' Actual vs Predicted Intensity' ], 'FontSize', self.TITLE_SIZE );
            xl = xlim();
            yl = ylim();
            mi = min( [xl(1) yl(1)] );
            ma = max( [xl(2) yl(2)] );
            xlim( [mi ma] );
            ylim( [mi ma] );
            plotfile = [ self.figdir self.name '_' self.kmer '_' id 'Predict.eps' ];
            print( f, plotfile, self.RES, self.FORMAT );
            close(f);
        end
        
        function self = crossValidation( self, kfold, lambdas, xvalerr, lambda, display )
            xlab = 'log_{10}\lambda';
            ylab = 'Mean Squared Error';
            plotfile = [ self.figdir self.name '_' self.kmer '_crossValidation.eps' ];
            
            plotParameterSelection( log10(lambdas), xvalerr, lambdas == lambda, ...
                xlab, ylab, [ self.titlename ' Model Selection' ], plotfile, self.RES, self.FORMAT, true );
            
            self.printf( 'Cross validation:\n' );
            self.down_level();
            previd = self.id;
            self.setId( [ previd '-xval' ] );
            self.printf( 'k-fold: %f\n', kfold );
            self.printf( 'Lambdas: ');
            fprintf( self.file, '%f ', lambdas );
            fprintf( self.file, '\n' );
            self.printf( 'MSE: ');
            fprintf( self.file, '%f ', xvalerr );
            fprintf( self.file, '\n');
            self.printAndDisplayf( display, 'Lambda: %f\n', lambda );
            self.printf( 'Plot: %s\n', plotfile );
            self.up_level();
            self.setId( previd );
        end
        
        function weightsfile = weights( self, W, headers, suffix )
            weightsfile = [ self.outdir self.name '_' self.kmer '_' suffix 'weights.txt' ];
            writeVector( weightsfile, W, headers );
            self.printf( 'Weights: %s\n', weightsfile );
        end
        
        function self = close( self )
            fclose( self.file );
        end
    
    end
    
end