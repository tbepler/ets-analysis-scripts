#include<string>
#include<vector>
#include<unordered_map>
#include<set>
#include<functional>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<locale>

using namespace std;

struct Feature{

    size_t pos;
    string kmer;

    inline bool containedIn( const string& str ){
        return str.substr( pos, kmer.size() ) == kmer;
    }

    bool operator== ( const Feature & f ) const{
        return pos == f.pos && kmer == f.kmer;
    }

};

namespace std{

    template<>
    struct hash<Feature>{

        size_t operator() ( const Feature& f ) const{
            hash<string> str_hash;
            size_t h = 5;
            h = h*31 + f.pos;
            h = h*31 + str_hash( f.kmer );
            return h;
        }

    };

}

static size_t featureHash( const Feature& f ){

    hash<string> str_hash;
    size_t h = 5;
    h = h*31 + f.pos;
    h = h*31 + str_hash( f.kmer );
    return h;

}

istream& operator>> ( istream& in, Feature& f ){

    string pos_kmer;
    in >> pos_kmer;
    //parse pos and kmer from string "[pos]kmer"
    replace( pos_kmer.begin(), pos_kmer.end(), '[', ' ' );
    replace( pos_kmer.begin(), pos_kmer.end(), ']', ' ' );
    stringstream ss( pos_kmer );
    ss >> f.pos;
    ss >> f.kmer;
    return in;

}

struct Model{

    double bias;
    unordered_map< Feature, double > features;
    set<size_t> sizes;
    size_t length;

    inline double score( const string& str ) const {
        double s = bias;
        Feature f;
        for( auto it = sizes.begin() ; it != sizes.end() ; ++it ){
            size_t size = *it;
            //expand all features of size from str
            for( size_t i = 0 ; i < str.size() - size + 1 ; ++i ){
                f.pos = i;
                f.kmer = str.substr( i, size );
                auto item = features.find( f );
                if( item != features.end() ){
                    s += (*item).second;
                }
            }
        }

        return s;
    }

    inline double maxScore( const string& str ) const {
        double max = -INFINITY;
        for( size_t i = 0 ; i < str.size() - length + 1 ; ++i ){
            double s = score( str.substr( i, length ) );
            if( s > max ) max = s;
        }
        return max;
    }

    inline double averageScore( const string& str ) const {
        double s = 0;
        for( size_t i = 0 ; i < str.size() - length + 1 ; ++i ){
            s += score( str.substr( i, length ) );
        }
        s = s / ( double ) ( str.size() - length + 1 );
        return s;
    }

};

istream& operator>> ( istream& in, Model& m ){

    m.features.clear();
    m.sizes.clear();
    size_t length = 0;
    string line;
    while( getline( in, line ) ){
        if( !line.empty() ){
            stringstream ss( line );
            if( line.substr( 0, 4 ) == "Bias" ){
                string ignore;
                ss >> ignore;
                ss >> m.bias;
            }else{
                Feature f;
                double w;
                ss >> f;
                ss >> w;
                if( w > 0 ){
                    m.features[ f ] = w;
                    m.sizes.insert( f.kmer.size() );
                }
                if( f.pos + f.kmer.size() > length ){
                    length = f.pos + f.kmer.size();
                }
            }
        }
    }
    m.length = length;
    return in;

}

inline string& upper( string& str ){
    for( size_t i = 0 ; i < str.size() ; ++i ){
        str[i] = toupper( str[i] );
    }
    return str;
}

inline char comp( char c ){

    switch( c ){
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return c;

}

inline string rvscomp( const string& str ){

    string rc;
    for( long i = str.size() - 1 ; i >= 0 ; --i ){
        rc += comp( str[i] );
    }
    return rc;

}

void scoreSeqs( const Model& m, istream& seqs, ostream& out ){
    string line;
    while( getline( seqs, line ) ){
        stringstream ss( line );
        string seq;
        ss >> seq;
        string ignore;
        ss >> ignore;
        string label;
        ss >> label;
        //upper case the sequence
        upper( seq );
        //ignore sequences containing anything other than A, C, G, T
        bool seqValid = true;
        char c;
        for( size_t i = 0 ; i < seq.size() ; ++i ){
            c = seq[i];
            if( c != 'A' && c != 'C' && c != 'G' && c != 'T' ){
                seqValid = false;
                break;
            }
        }
        if( !seqValid ){
            continue;
        }
        //get reverse compliment of seq
        string rvs = rvscomp( seq );

        //assign score to be average of scores on fwd and rvs
        //double fwdScore = m.averageScore( seq );
        //double rvsScore = m.averageScore( rvs );
        //double score = ( fwdScore + rvsScore ) / 2.0;
        
        //assign score to be max of scores on fwd and rvs
        double fwdScore = m.maxScore( seq );
        double rvsScore = m.maxScore( rvs );
        double score = max( fwdScore, rvsScore );

        out << label << " " << score << endl;

    }
}

int main( int argc, const char* argv[] ){

    ifstream m_in( argv[1] );
    Model m;
    m_in >> m;
    m_in.close();

    ifstream seq_in( argv[2] );
    scoreSeqs( m, seq_in, cout );
    seq_in.close();

    return 0;

}


