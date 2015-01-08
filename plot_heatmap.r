require("gplots")
require("RColorBrewer")

DPI <- 300
WIDTH <- 10
HEIGHT <- 10
FONT <- 12

args <- commandArgs( trailingOnly = TRUE )

palette <- colorRampPalette( c( "red", "yellow", "green" ) )( n = 299 )

for( file in args ){

    filename <- strsplit( file, "\\." )[[1]][1]
    outname <- paste( c( filename, "_heatmap.png" ), collapse="" )

    data <- read.table( file, header = TRUE, row.names = 1 )
    m <- data.matrix( data )
    colnames( m ) <- colnames( data )
    rownames( m ) <- rownames( data )

    png( outname, width = WIDTH*DPI, height = HEIGHT*DPI, res = DPI, pointsize = FONT )

    heatmap.2(
        m,
        density.info = "none",
        trace = "none",
        col = palette,
        dendrogram = "both"
        )

    def.off()

}
