#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<unistd.h>
#include<cmath>
#include"Parser.h"
#include"Motif.h"
#include"DNA.h"

using namespace std;

//void processStringDouble( istream& in, 

inline double maxScoreCenter( const Motif& m, const std::string& seq ){
    std::vector<unsigned long> trnsl;
    m.translate( seq, trnsl );
    //only score the center motif of the seq
    unsigned long start_pos = ( trnsl.size() - m.size() ) / 2;
    unsigned long* start = &trnsl[0] + start_pos;
    return m.score( start , start + m.size() );
}

inline double maxScoreAll( const Motif& m, const std::string& seq ){
    std::vector<unsigned long> trnsl;
    m.translate( seq, trnsl );
    return m.maxScore( &trnsl[0], &trnsl[0] + seq.size() );
}

class Function{

    public:
        ostream& fOut;
        ostream& rOut;
        bool labels;
        const Motif& m;
        bool eliminate;

        Function( const Motif& motif, ostream& fOut = cout, ostream& rOut = cout, bool labels = true, bool elim = false )
            : m(motif), fOut( fOut ), rOut( rOut ), labels( labels ), eliminate( elim ) {};
        
        void operator()(const StringDouble& sd ){
            std::string seq = sd.first;
            std::string rvscomp = reverseCompliment( seq );
            double d = sd.second;
            double fwdMax = maxScoreCenter( m, seq );
            double rvsMax = maxScoreCenter( m, rvscomp );
            bool elim = false;
            if( eliminate ){
                double s = max( maxScoreAll( m, seq ), maxScoreAll( m, rvscomp ) );
                elim = s > fwdMax && s > rvsMax;
            }
            if( !elim ){
                char c = fwdMax >= rvsMax ? 'F' : 'R' ;
                ostream& out = fwdMax >= rvsMax ? fOut : rOut;
                out << seq << " " << d;
                if( labels ) out << " " << c;
                out << endl;
            }
        }

};

void usage(){
    cerr << "orient -p PWM [-s SEQS] [-f FWD OUTPUT] [-r RVS OUTPUT] [-l] [-e]" << endl;
}

int main( int argc, char* argv[] ){

    if( argc <= 1 ){
        usage();
        return 0;
    }

    ifstream pwmFile;
    ifstream seqFile;
    ofstream fwdFile;
    ofstream rvsFile;
    bool labels = true;
    bool elim = false;

    char* val;
    int c;
    while( ( c = getopt( argc, argv, "p:s:f:r:le" ) ) != -1 ){
        switch(c){
            case 'p':
                val = optarg;
                pwmFile.open( val );
                if( !pwmFile.is_open() ){
                    cerr << "Unable to open file: " << val << endl;
                    return 1;
                }
                break;
            case 's':
                val = optarg;
                seqFile.open( val );
                if( !seqFile.is_open() ){
                    cerr << "Unable to open file: " << val << endl;
                    return 1;
                }
                break;
            case 'f':
                val = optarg;
                fwdFile.open( val );
                if( !fwdFile.is_open() ){
                    cerr << "Unable to open file: " << val << endl;
                    return 1;
                }
                break;
            case 'r':
                val = optarg;
                rvsFile.open( val );
                if( !rvsFile.is_open() ){
                    cerr << "Unable to open file: " << val << endl;
                    return 1;
                }
                break;
            case 'l':
                labels = false;
                break;
            case 'e':
                elim = true;
                break;
        }
    }

    if( !pwmFile.is_open() ){
        cerr << "Error: a PWM must be specified" << endl;
        usage();
        return 1;
    }

    istream& seqIn = seqFile.is_open() ? seqFile : cin;
    ostream& fOut = fwdFile.is_open() ? fwdFile : cout;
    ostream& rOut = rvsFile.is_open() ? rvsFile : cout;

    Motif m;
    pwmFile >> m;
    pwmFile.close();

    Function f( m, fOut, rOut, labels, elim );
    processStringDouble( seqIn, f );

    return 0;

    /*
    switch(argc){

        case 2:{
                ifstream pwmIn( argv[1] );
                if( pwmIn.is_open() ){
                    Motif m;
                    pwmIn >> m;
                    pwmIn.close();
                    Function f( m );
                    processStringDouble( cin, f );
                }else{
                    cerr << "Error opening file: " << argv[1] << endl;
                }
            }
            break;

        case 3:{
                ifstream pwmIn( argv[1] );
                if( pwmIn.is_open() ){
                    Motif m;
                    pwmIn >> m;
                    pwmIn.close();
                    Function f( m );
                    ifstream in( argv[2] );
                    if( in.is_open() ){
                        processStringDouble( in, f );
                        in.close();
                    }else{
                        cerr << "Error opening file: " << argv[2] << endl;
                    }
                }else{
                    cerr << "Error opening file: " << argv[1] << endl;
                }
            }
            break;

        default:
            cerr << "Usage: orient PWM_FILE [SEQS_FILE] [ < SEQS_FILE ]" << endl;
            break;
    }
            
    */

}
