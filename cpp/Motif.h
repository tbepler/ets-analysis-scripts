#ifndef MOTIF_INCLUDED
#define MOTIF_INCLUDED

#include<iostream>
#include<vector>
#include<string>

#define MOTIF_MAX_BYTES 200e6


class Motif{
    
    public:
        Motif () : data(NULL), k(0), n(0), scoreMatrix(NULL), alphabet(NULL)  {};
        Motif( const double* scoreMatrix, unsigned long k, unsigned long n,
                const char* alphabet, unsigned long maxBytes = MOTIF_MAX_BYTES );

        ~Motif();

        static Motif* parse( std::istream& i, unsigned long maxBytes = MOTIF_MAX_BYTES );

        inline unsigned long size() const { return k; }; 
        inline bool empty() const { return k <= 0 && n <= 0; };
        
        double score( unsigned long* start, unsigned long* end ) const;
        void scoreAll( unsigned long* start, unsigned long* end, double* buffer ) const;
        void scoreAll( unsigned long* start, unsigned long* end, std::vector<double>& buffer ) const;
        double maxScore( unsigned long* start, unsigned long* end ) const; 

        inline unsigned long translate( char c ) const{
            char* alph = this->alphabet;
            unsigned long n = this->n;
            for( unsigned long i = 0 ; i < n ; ++i ){
                if( alph[i] == c ) return i;
            }
            return n;
        }

        inline void translate( char* str, std::vector<unsigned long>& buffer ) const {
            char c;
            unsigned long i = 0;
            while( (c = str[i++]) != '\0' ){
                buffer.push_back( translate( c ) );
            }
        };
        
        inline void translate( std::string str, std::vector<unsigned long>& buffer ) const {
            char c;
            for( unsigned long i = 0 ; i < str.size() ; ++i ){
                c = str[i];
                buffer.push_back( translate( c ) );
            }
        }

        friend std::istream& operator>>( std::istream& i, Motif& m );
        friend std::ostream& operator<<( std::ostream& o, const Motif& m );
        std::string toString() const;


    private:
        struct Data;
        Data* data;

        unsigned long k;
        unsigned long n;
        double* scoreMatrix;
        char* alphabet;


};






#endif
