classdef Report < handle
    
    properties
        file;
        indent;
        level;
        id;
        outdir;
    end
    
    methods
        
        function obj = Report( filepath, indent, level, id, outdir )
            
            if nargin < 5
                [ outdir, ~, ~ ] = fileparts( filepath );
                outdir = [ outdir '/' ];
            end
            
            if nargin < 4
                id = '';
            end
            
            if nargin < 3
                level = 0;
            end
            
            if nargin < 2
                indent = '    ';
            end
            
            obj.file = fopen( filepath, 'w' );
            obj.indent = indent;
            obj.level = level;
            obj.id = id;
            obj.outdir = outdir;
            
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
        
        function self = setID( self, id )
            self.id = id;
        end
        
        function self = roc( self, seqspath, weightspath, mtype, display )
            if nargin < 5
                display = false;
            end
            [ auc, opt, rocfig ] = generateROC( seqspath, weightspath, mtype, self.outdir, display );
            record( 'ROC:\n' );
            self.down_level();
            previd = self.id;
            self.setID( [ previd '-roc' ] );
            self.printf( 'AUC: %f\n', auc );
            self.printf( 'OPT: %f, %f\n', opt(1), opt(2) );
            self.printf( 'Plot: %s\n', rocfig );
            self.up_level();
            self.setID( previd );
        end
    
    end
    
end