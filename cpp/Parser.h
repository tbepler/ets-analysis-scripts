#ifndef PARSER_INCLUDED
#define PARSER_INCLUDED

#include<iostream>
#include<string>
#include<sstream>

template< typename T1, typename T2 >
struct Tuple2{
    T1 first;
    T2 second;

    friend std::istream& operator>>( std::istream& i, Tuple2<T1,T2>& t ){
        return i >> t.first >> t.second;
    }

    friend std::ostream& operator<<( std::ostream& o, Tuple2<T1,T2>& t ){
        return o << t.first << " " << t.second;
    }

};

typedef struct Tuple2<std::string,double> StringDouble;

template< typename Function >
void processStringDouble( std::istream& in, Function f ){
    StringDouble obj;
    std::string line;
    while( std::getline( in, line ) ){
        if( !line.empty() ){
            std::stringstream ss( line );
            ss >> obj;
            f( obj );
        }
    }
}




#endif
