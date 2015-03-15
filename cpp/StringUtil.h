#ifndef STRINGUTIL_INCLUDED
#define STRINGUTIL_INCLUDED

#include<string>
#include<vector>
#include<algorithm>

#define WHITESPACE " \t\n\v\f\r"

inline bool startsWith( const std::string& str, const std::string& prefix ){
    return str.substr( 0, prefix.size() ) == prefix;
}

inline std::string toUpper( const std::string& str ){
    std::string upper;
    upper.resize( str.size() );
    std::transform( str.begin(), str.end(), upper.begin(), ::toupper );
    return upper;
}

inline std::string toLower( const std::string& str ){
    std::string lower;
    lower.resize( str.size() );
    std::transform( str.begin(), str.end(), lower.begin(), ::tolower );
    return lower;
}

inline std::vector<std::string> split( const std::string& str, const std::string& delim ){
    std::vector<std::string> tokens;
    std::size_t pos = 0;
    std::size_t prev = pos;
    while( ( pos = str.find( delim, pos ) ) < std::string::npos ) {
        tokens.push_back( str.substr( prev, pos ) );
        pos += delim.size();
        prev = pos;
    }
    //add last token
    tokens.push_back( str.substr( prev, str.size() ) );
    
    return tokens;
}

inline std::vector<std::string>& split( const std::string& str, const std::string& delim,
                                        std::vector<std::string>& tokens ){
    std::size_t pos = 0;
    std::size_t prev = pos;
    while( ( pos = str.find( delim, pos ) ) < std::string::npos ) {
        tokens.push_back( str.substr( prev, pos ) );
        pos += delim.size();
        prev = pos;
    }
    //add last token
    tokens.push_back( str.substr( prev, str.size() ) );
    
    return tokens;
}

inline std::vector<std::string> splitAny( const std::string& str, const std::string& delims = WHITESPACE ){
    std::vector<std::string> tokens;
    std::size_t pos = 0;
    std::size_t prev = pos;
    while( ( pos = str.find_first_of( delims, pos ) ) < std::string::npos ) {
        tokens.push_back( str.substr( prev, pos ) );
        ++pos;
        prev = pos;
    }
    //add last token
    tokens.push_back( str.substr( prev, str.size() ) );
    
    return tokens;
}

inline std::vector<std::string>& splitAny( const std::string& str, std::vector<std::string>& tokens,
                                            const std::string& delims = WHITESPACE){
    std::size_t pos = 0;
    std::size_t prev = pos;
    while( ( pos = str.find_first_of( delims, pos ) ) < std::string::npos ) {
        tokens.push_back( str.substr( prev, pos ) );
        ++pos;
        prev = pos;
    }
    //add last token
    tokens.push_back( str.substr( prev, str.size() ) );
    return tokens;
}



#endif


