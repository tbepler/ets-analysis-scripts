#include"Motif.h"
#include"StringUtil.h"
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sstream>

#define M Motif

using namespace std;

typedef unsigned long ulong;

inline double safeLog( double n ){
    if( n == 0 ) return -INFINITY;
    return log( n );
}

inline void log( vector<double>& vec ){
    for( unsigned long i = 0 ; i < vec.size() ; ++i ){
        vec[i] = safeLog( vec[i] );
    }
}

inline double log2( double n )
{
    return log( n ) / log ( 2 );
}

inline unsigned long nextPowerOf2( unsigned long n )
{
    double exp = ceil( log2( n ) );
    n = (unsigned long) pow( 2, exp );
    return n;
}

typedef struct M::Data{
    double* kmerScores;
    ulong bits;
    ulong mask;
} Data;

//TODO implement data initialization
inline Data* initializeData( const double* scoreMatrix, ulong k, ulong n,
                                const char* alphabet, ulong maxBytes ){
    ulong a = nextPowerOf2( n );
    ulong aBytes = log2( a );
    ulong kmers = ( unsigned long ) pow( a, k );
    ulong kmerScoresBytes = kmers * sizeof( double );
    if( kmerScoresBytes <= MOTIF_MAX_BYTES ){

    }

    return NULL;
}

M::M( const double* scoreMatrix, ulong k, ulong n, const char* alphabet, ulong maxBytes ){
     
    this->data = initializeData( scoreMatrix, k, n, alphabet, maxBytes );
    this->k = k;
    this->n = n;
    ulong alphBytes = n * sizeof( *alphabet );
    this->alphabet = (char*) malloc( alphBytes );
    memcpy( this->alphabet, alphabet, alphBytes );
    ulong scoreMBytes = n * k * sizeof( *(this->scoreMatrix) );
    this->scoreMatrix = (double*) malloc( scoreMBytes );
    for( ulong i = 0 ; i < k ; ++i ){
        for( ulong j = 0 ; j < n ; ++j ){
            this->scoreMatrix[ n*i + j ] = scoreMatrix[i*n + j];
        }
    }

}

M::~M(){
    if( this->data ) free( this->data );
    this->data = NULL;
    if( this->alphabet ) free( this->alphabet );
    this->alphabet = NULL;
    if( this->scoreMatrix ) free( this->scoreMatrix );
    this->scoreMatrix = NULL;
    this->k = 0;
    this->n = 0;
}

inline double sum( double* start, double* end ){
    double s = 0;
    for( double* i = start ; i < end ; ++i ){
        s += *i;
    }
    return s;
}

inline void toProbabilities( vector<double>& occurences, const ulong k, const ulong n){
    for( ulong i = 0 ; i < k ; ++i ){
        double s = sum( &occurences[i*n], &occurences[i*n + n] );
        for( ulong j = 0 ; j < n ; ++j ){
            occurences[ i*n + j ] = occurences[ i*n + j ] / s;
        }
    }
}

#define PROB 0
#define LOG 1
#define OCCURENCES 2

M* M::parse( istream& in, ulong maxBytes ){
    
    string line;
    //parse alphabet from first line
    getline( in, line );
    vector<string> tokens;
    splitAny( line, tokens );
    unsigned long n = tokens.size();
    char* alph = (char*) malloc( n * sizeof(*alph) );
    for( unsigned long i = 0 ; i < n ; ++i ){
        alph[i] = tokens[i][0];
    }

    //parse score matrix
    unsigned long k = 0;
    vector<double> scoreMatrix;
    stringstream ss;
    double entry;
    int type = PROB;
    while( getline( in, line ) && !line.empty() ){
        ++k;
        ss.clear();
        ss.str(line);
        for( int i = 0 ; i < n ; ++i ){
            ss >> entry;
            scoreMatrix.push_back( entry );
            if( entry < 0 ) type = LOG;
            if( entry > 1 ) type = OCCURENCES;
        }
    }

    switch(type){
        case PROB:
            log( scoreMatrix );
            break;
        case OCCURENCES:
            toProbabilities( scoreMatrix, k, n );
            log( scoreMatrix );
            break;
    }

    //create motif
    return new Motif( &scoreMatrix[0], k, n, alph, maxBytes );

}

istream& operator>>( istream& in, M& motif ){
    //clear motif
    motif.~M();

    string line;
    //parse alphabet from first line
    if( !getline( in, line ) || line.empty() ){
        return in;
    }
    vector<string> tokens;
    splitAny( line, tokens );
    unsigned long n = tokens.size();
    char* alph = (char*) malloc( n * sizeof(*alph) );
    for( unsigned long i = 0 ; i < n ; ++i ){
        alph[i] = tokens[i][0];
    }

    //parse score matrix
    unsigned long k = 0;
    vector<double> scoreMatrix;
    stringstream ss;
    double entry;
    int type = PROB;
    while( getline( in, line ) && !line.empty() ){
        ++k;
        ss.clear();
        ss.str(line);
        for( int i = 0 ; i < n ; ++i ){
            ss >> entry;
            scoreMatrix.push_back( entry );
            if( entry < 0 ) type = LOG;
            if( entry > 1 ) type = OCCURENCES;
        }
    }
   
    switch(type){
        case PROB:
            log( scoreMatrix );
            break;
        case OCCURENCES:
            toProbabilities( scoreMatrix, k, n );
            log( scoreMatrix );
            break;
    }

    //allocate new motif where motif is
    new(&motif) M( &scoreMatrix[0], k, n, alph );
    
    return in;
}

template< typename T >
inline ostream& report( ostream& o, T* start, T* end, const string& delim = " " ){
    for( T* i = start ; i < end ; ++i ){
        if( i != start ) o << delim;
        o << *i;
    }
    return o;
}

ostream& operator<<( ostream& o, const M& motif ){
    double* scoreMatrix = motif.scoreMatrix;
    ulong k = motif.k;
    ulong n = motif.n;
    char* alph = motif.alphabet;
    //first output the alphabet header
    report( o, alph, alph + n );
    //report the score matrix
    for( ulong i = 0 ; i < k ; ++i ){
        o << endl;
        report( o, scoreMatrix + i*n, scoreMatrix + (i+1)*n );
    }
    return o;
}

string M::toString() const{
    stringstream ss;
    ss << this;
    return ss.str();
}

inline double scoreUnsafe( const double* scoreMatrix, const ulong n, ulong* start, ulong* end ){
    double score = 0;
    ulong k = end - start;
    for( ulong i = 0 ; i < k ; ++i ){
        score += scoreMatrix[ i * n + *(start + i ) ];
    }
    return score;
}

//TODO include usage of hashed kmer scores
double M::score( ulong* start, ulong* end ) const{
    double* scoreMatrix = this->scoreMatrix;
    ulong n = this->n;
    ulong k = this->k;
    //check that input is valid
    assert( k == end - start );
    for( ulong* i = start ; i < end ; ++i ){
        assert( *i < n );
    }
    
    return scoreUnsafe( scoreMatrix, n, start, end );
    
}

inline ulong endIndex( const ulong k, const ulong* start, const ulong* end ){
    return ( end - start ) - k + 1;
}

inline void scoreAllUnsafe( const double* scoreMatrix, const ulong k, const ulong n,
                        ulong* start, ulong* end, double* buffer ){
    
    ulong e = endIndex( k, start, end );
    for( ulong i = 0 ; i < e ; ++i ){
        buffer[i] = scoreUnsafe( scoreMatrix, n, start + i, start + i + k );
    }

}

inline void scoreAllUnsafe( const double* scoreMatrix, const ulong k, const ulong n,
                        ulong* start, ulong* end, vector<double>& buffer ){
    ulong e = endIndex( k, start, end );
    for( ulong i = 0 ; i < e ; ++i ){
        buffer.push_back( scoreUnsafe( scoreMatrix, n, start + i, start + i + k ) );
    }
}

inline double maxScoreUnsafe( const double* scoreMatrix, const ulong k, const ulong n,
                        ulong* start, ulong* end ){
    end = end - k + 1;
    double max = -INFINITY;
    double s;
    for( ulong* i = start ; i < end ; ++i ){
        s = scoreUnsafe( scoreMatrix, n, i, i + k );
        if( s > max ) max = s;
    }
    return max;
}

inline bool scoreAllConditionsMet( const double* scoreMatrix, const ulong k, const ulong n,
                                    ulong* start, ulong* end ){
    if( k > end - start ) return false;
    for( ulong* i = start ; i < end ; ++i ){
        if( *i >= n ) return false;
    }
    return true;
}

//TODO include usage of hashed kmer scores
void M::scoreAll( ulong* start, ulong* end, double* buffer ) const{
    double* scoreMatrix = this->scoreMatrix;
    ulong n = this->n;
    ulong k = this->k;
    //check that input is valid
    assert( scoreAllConditionsMet( scoreMatrix, k, n, start, end ) );
    
    scoreAllUnsafe( scoreMatrix, k, n, start, end, buffer );

}

//TODO include usage of hashed kmer scores
void M::scoreAll( ulong* start, ulong* end, vector<double>& buffer ) const{
    double* scoreMatrix = this->scoreMatrix;
    ulong n = this->n;
    ulong k = this->k;
    //check that input is valid
    assert( scoreAllConditionsMet( scoreMatrix, k, n, start, end ) );
    
    scoreAllUnsafe( scoreMatrix, k, n, start, end, buffer );

}

inline bool maxScoreConditionsMet( const double* scoreMatrix, const ulong k, const ulong n,
                                ulong* start, ulong* end ){
    return scoreAllConditionsMet( scoreMatrix, k, n, start, end );
}

//TODO include usage of hashed kmer scores
double M::maxScore( ulong* start, ulong* end ) const{
    double* scoreMatrix = this->scoreMatrix;
    ulong n = this->n;
    ulong k = this->k;
    //check that input is valid
    assert( maxScoreConditionsMet( scoreMatrix, k, n, start, end ) );
    
    return maxScoreUnsafe( scoreMatrix, k, n, start, end );

}

