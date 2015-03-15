#include<stdio.h>
#include<iostream>
#include"GenomicRegion.h"
#include"Chromosome.h"
#include"Motif.h"

//#define CHROM
//#define GR
#define MOTIF

#define K 8
#define N 4

using namespace std;

int main( int argc, const char* argv[] ){
    
#ifdef CHROM
    cout << "Input 2 chromosomes:" << endl;
    Chromosome c1;
    Chromosome c2;
    cin >> c1 >> c2;
    hash<Chromosome> chrHash;
    cout << "c1 = " << c1 << ", hash = " << chrHash(c1) << endl;
    cout << "c2 = " << c2 << ", hash = " << chrHash(c2) << endl;
    if (c1 == c2) cout << c1 << " == " << c2 << endl;
    if (c1 != c2) cout << c1 << " != " << c2 << endl;
    if (c1 < c2) cout << c1 << " < " << c2 << endl;
    if (c1 <= c2) cout << c1 << " <= " << c2 << endl;
    if (c1 > c2) cout << c1 << " > " << c2 << endl;
    if (c1 >= c2) cout << c1 << " >= " << c2 << endl;
#endif

#ifdef GR
    GenomicRegion r1;
    GenomicRegion r2;
    cout << "Input 2 genomic regions:" << endl;
    cin >> r1 >> r2;
    hash<GenomicRegion> grHash;
    cout << "r1 = " << r1 << ", hash = " << grHash( r1 ) << endl;
        inline std::string toString() const{
            std::stringstream ss;
            ss << this;
            return ss.str();
        }
    cout << "r2 = " << r2 << ", hash = " << grHash( r2 ) << endl;
    if ( r1 == r2 ) cout << r1 << " == " << r2 << endl;
    if ( r1 != r2 ) cout << r1 << " != " << r2 << endl;
    if ( r1 < r2 ) cout << r1 << " < " << r2 << endl;   
    if ( r1 <= r2 ) cout << r1 << " <= " << r2 << endl;
    if ( r1 > r2 ) cout << r1 << " > " << r2 << endl;
    if ( r1 >= r2 ) cout << r1 << " >= " << r2 << endl;
#endif

#ifdef MOTIF
/*
    char A[] = { 'A', 'C', 'G', 'T' };
    double M[][N] = {
        { 1, 0, 0, 0 },
        { 1, 1, 2, 5 },
        { 1, 0, 0, 1 },
        { 1, 5, 2, 1 },
        { 0, 0, 10, 0 },
        { 5, 1, 1, 1 },
        { 1, 5, 5, 3 },
        { 1, 1, 1, 1 }
    };
    Motif m( M[0], K, N, A );
    */
    cout << "Enter PWM:" << endl;
    //Motif* mp = Motif::parse( cin );
    //if( !mp ) return 0;
    //Motif m = *mp;
    Motif m;
    cin >> m;
    if( m.empty() ) return 0;
    cout << m << endl;
    string line;
    vector<unsigned long> seq;
    cout << "Input sequences to score: " << endl;
    while( getline( cin, line ) ){
        if( line == "q" || line == "Q" ) break;
        seq.clear();
        m.translate( line, seq );
        for( int i = 0 ; i < seq.size() ; ++i ){
            cout << seq[i] << " ";
        }
        cout << endl;
        vector<double> scores;
        unsigned long* start = &seq[0];
        m.scoreAll( start, start + seq.size(), scores );
        for( int i = 0 ; i < scores.size() ; ++i ){
            cout << i << " - " << i + m.size() -1 << ": " << scores[i] << endl;
        }
        cout << "Max: " << m.maxScore( start, start + seq.size() ) << endl;
    }
#endif

}
