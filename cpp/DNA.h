#ifndef DNA_INCLUDED
#define DNA_INCLUDED

inline char complement( const char c ){
    switch(c){
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return c;
}

inline std::string reverseCompliment( const std::string& seq ){
    std::string rc;
    unsigned long n = seq.size();
    unsigned long last = n - 1;
    rc.resize( n );
    for( unsigned long i = 0 ; i < n ; ++i ){
        rc[ last - i ] = complement( seq[i] );
    }
    return rc;
}


#endif
