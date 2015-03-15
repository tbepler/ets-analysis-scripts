#ifndef CHROMOSOME_INCLUDED
#define CHROMOSOME_INCLUDED

#include<iostream>
#include<string>

class Chromosome{

    private:
        std::string mName;
        friend inline std::ostream& operator<<( std::ostream& o, const Chromosome& chrom ){
            return o << chrom.mName;
        }
        friend inline std::istream& operator>>( std::istream& i, Chromosome& chrom ){
            return i >> chrom.mName;
        }

    public:
        Chromosome( ) : mName("") {};
        Chromosome( const std::string& name );
        inline std::string toString() const{ return mName; };

        inline std::size_t hash() const{
            return std::hash<std::string>()( mName );
        }

        inline int compareTo( const Chromosome& c ) const{
            return mName.compare( c.mName );
        }
        inline bool operator==( const Chromosome& rhs ) const{
            return mName == rhs.mName;
        }
        inline bool operator!=( const Chromosome& rhs ) const{
            return mName != rhs.mName;
        }
        inline bool operator<( const Chromosome& rhs ) const{
            return mName < rhs.mName;
        }
        inline bool operator<=( const Chromosome& rhs ) const{
            return mName <= rhs.mName;
        }
        inline bool operator>( const Chromosome& rhs ) const{
            return mName > rhs.mName;
        }
        inline bool operator>=( const Chromosome& rhs ) const{
            return mName >= rhs.mName;
        }
};

namespace std{
    template<>
    class hash<Chromosome> {
         public:
             size_t operator()(const Chromosome& c ) const{
                return c.hash();
            }
    };
}

#endif
