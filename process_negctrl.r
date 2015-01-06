args <- commandArgs( trailingOnly = TRUE )

for ( file in args ){

    print( sprintf( "Processing file: %s", file ) )

    out <- gsub( "_extracted.txt", ".txt", file )
    out <- gsub( ".txt", "_negctrl.txt", out )

    data <- read.table( file, stringsAsFactors = FALSE )
    idx = grep( "NegControl", data[,1] )
    data <- data[ idx, -1 ]
    data[ , 2 ] <- log( as.numeric( data[ , 2 ] ) )

    print( sprintf( "Writing to file: %s", out ) )
    write.table( data, file = out, quote = FALSE, row.names = FALSE, col.names = FALSE )

}
