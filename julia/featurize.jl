immutable KmerSpace
    alphabet
    dict
    ks
    function KmerSpace( alphabet, ks )
        dict = [ alphabet[i]=>i-1 for i in 1:Base.length(alphabet) ]
        new( alphabet, dict, ks )
    end
end

nkmers( k, alphabet ) = Base.length( alphabet ) .^ k

Base.length( ks::KmerSpace ) = sum( nkmers( ks.ks, ks.alphabet ) )

function make_kmer( k, alphabet, idx )
    n = Base.length( alphabet )
    kmer = Array( eltype(alphabet), k )
    #use 1 as starting index for consistency
    idx -= 1
    for i = k:-1:1
        offset = idx % n
        kmer[i] = alphabet[ offset + 1 ]
        idx = floor( idx / n )
    end
    return kmer
end

function kmer_index( string, dict )
    idx = 0
    n = Base.length( dict )
    for i = 1:Base.length(string)
        idx *= n
        idx += dict[ string[i] ]
    end
    return idx + 1
end

type KmerSpaceState
    i
    ki
    ke
end

Base.start( ks::KmerSpace ) = KmerSpaceState( 1, 1, nkmers( ks.ks[1], ks.alphabet ) )
function Base.next( ks::KmerSpace, s::KmerSpaceState )
    i = s.i
    ki = s.ki
    ke = s.ke
    if i > ke
        i = 1
        ki += 1
        ke = nkmers( ks.ks[ki], ks.alphabet )
    end
    kmer = make_kmer( ks.ks[ki], ks.alphabet, i )
    next_state = KmerSpaceState( i+1, ki, ke )
    return kmer, next_state
end
Base.done( ks::KmerSpace, s::KmerSpaceState ) = s.i > s.ke && s.ki >= Base.length( ks.ks )

function nfeatures( len::Int, kmer_space )
    n = 0
    for k in kmer_space.ks
        n += ( len - k + 1 ) * nkmers( k, kmer_space.alphabet )
    end
    return n    
end

nfeatures( string, kmer_space ) = nfeatures( Base.length(string), kmer_space )

function featurize( string, kmer_space )
    fs = falses( nfeatures( string, kmer_space ) )
    offset = 0
    for k = kmer_space.ks
        n = nkmers( k, kmer_space.alphabet )
        for p = 1:Base.length(string)-k+1
            kmer = string[p:p+k-1]
            i = kmer_index( kmer, kmer_space.dict )
            fs[ offset + (p-1)*n + i ] = true
        end
        offset += n*(Base.length(string)-k+1)
    end
    return fs
end

function feature_names( len, kmer_space )
    names = Array( ASCIIString, nfeatures( len, kmer_space ) )
    offset = 0
    for k = kmer_space.ks
        n = nkmers( k, kmer_space.alphabet )
        for p = 1:len-k+1
            for i = 1:n
                kmer = join( make_kmer( k, kmer_space.alphabet, i ) )
                names[ offset + (p-1)*n + i ] = string( "[",p,"]",kmer )
            end
        end
        offset += n*(len-k+1)
    end
    return names
end

function featurize_all( xs, kmer_space )
    m = size( xs, 1 )
    if m > 0
        n = nfeatures( xs[1], kmer_space )
        fs = Array( Bool, m, n )
        @parallel for i = 1:m
            fs[i,:] = featurize( xs[i], kmer_space )
        end
        return fs, feature_names( Base.length( xs[1] ), kmer_space )
    end
    return Array( bool, 0, 0 ), Array( ASCIIString, 0 )
end


