#ifndef GENOMICREGION_INCLUDED
#define GENOMICREGION_INCLUDED

#include<iostream>
#include<string>
#include"Chromosome.h"

class GenomicRegion{

    private:
        unsigned long mStart;
        unsigned long mEnd;
        Chromosome mChrom;

        friend inline std::ostream& operator<<( std::ostream& o, const GenomicRegion& r ){
            return r.report( o );
        }
        friend inline std::istream& operator>>( std::istream& i, GenomicRegion& r ){
            return i >> r.mChrom >> r.mStart >> r.mEnd ;
        }

    public:
        
        GenomicRegion( ) : mChrom(), mStart(0), mEnd(0) {};
        GenomicRegion( const Chromosome&, const unsigned long, const unsigned long );
        
        inline unsigned long start() const{ return mStart; };
        inline unsigned long end() const{ return mEnd; };
        inline const Chromosome& chromosome() const{ return mChrom; };
        
        inline std::ostream& report( std::ostream& o, const std::string& delim = " " ) const{ 
            return o << mChrom << delim << mStart << delim << mEnd;
        };

        inline std::string toString(const std::string& delim = " ") const{
            return ((mChrom.toString() += delim)
                += std::to_string( mStart )
                += delim) += std::to_string( mEnd );
        };
       
        inline std::size_t hash() const{
            std::size_t h = 5;
            h += 11 * mChrom.hash();
            h += 11 * mStart;
            h += 11 * mEnd;
            return h;
        }

        inline int compareTo( const GenomicRegion& r ) const{
            if( this == &r ) return 0;
            int comp = mChrom.compareTo( r.mChrom );
            if( comp != 0 ) return comp;
            comp = mStart - r.mStart;
            if( comp != 0 ) return comp;
            comp = mEnd - r.mEnd;
            return comp;
        }

        inline bool operator==( const GenomicRegion& r ) const{
            return mChrom == r.mChrom && mStart == r.mStart && mEnd == r.mEnd;
        };
        inline bool operator!=( const GenomicRegion& r ) const{
            return !this->operator==(r);
        };
        inline bool operator< ( const GenomicRegion& r ) const{
            return this->compareTo( r ) < 0;
        };
        inline bool operator<= ( const GenomicRegion& r ) const{
            return this->compareTo( r ) <= 0;
        };
        inline bool operator> ( const GenomicRegion& r ) const{
            return this->compareTo( r ) > 0;
        }
        inline bool operator>= ( const GenomicRegion& r ) const{
            return this->compareTo( r ) >= 0;
        }

};

namespace std{
    template<>
    class hash<GenomicRegion> {
        public:
            size_t operator()( const GenomicRegion& r ) const{
                return r.hash();
            }
    };
}

#endif
