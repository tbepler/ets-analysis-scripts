function load_sequences( file_path )
    seqs = ASCIIString[]
    scores = Float64[]
    open( file_path ) do f
        for line in eachline( f )
            spl = split( line )
            push!( seqs, spl[1] )
            push!( scores, float64( spl[2] ) )
        end
    end
    return seqs,scores
end
