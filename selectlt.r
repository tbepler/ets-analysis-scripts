args <- commandArgs( trailingOnly = TRUE )

lt = as.numeric( args[1] )

if( length( args ) > 1 ){
    data <- read.table( args[2] )
}else{
    data <- read.table( stdin() )
}

data <- data[ as.numeric( data[ , 2 ] ) < lt , ]

write.table( data , quote = FALSE, row.names = FALSE, col.names = FALSE )
