function positional_kmer_kernel( str1, str2 )
    s = 0
    prev = 0
    for i = 1:min( Base.length(str1), Base.length(str2) )
        if str1[i] == str2[i]
            prev += 1
            s += prev
        else
            prev = 0
        end
    end
    return s
end

function positional_kmer_kernel( n, str1, str2 )
    s = 0
    prev= 0
    for i = 1:min( Base.length(str1), Base.length(str2) )
        if str1[i] == str2[i]
            prev = min( n, prev + 1 )
            s += prev
        else
            prev = 0
        end
    end
    return s
end

function positional_kmer_kernel( n )
    kernel(a,b) = positional_kmer_kernel(n,a,b)
    return kernel
end
