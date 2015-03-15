#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <getopt.h>

using namespace std;

static const char ALPH[] = {'A', 'C', 'G', 'T'};
static const int ALPH_SIZE = 4;

static const int KMERS[] = {1,2,3};
static const int KMERS_SIZE = 3;

void usage ( ostream& out )
{
	out << "Usage: featurize [--input/-i FILE] [--output/-o FILE] -k INT1 INT2 ... " << endl;
}

int pow(int base, int exp)
{
	int result = 1;
	while (exp)
	{
		if ( exp & 1 )
		{
			result *= base;
		}
		exp >>= 1;
		base *= base;
	}
	return result;

}

void permutation( const vector<char>& alphabet, int v, int lkmer, char* kmer )
{
	int k = alphabet.size();
	for ( int i = 0 ; i < lkmer ; ++i )
	{
		kmer[lkmer - i - 1] = alphabet[ ( v / pow(k,i) ) % k ];
	}
}

bool matches( const string& str, int lkmer, const char* kmer, int start )
{
	if ( str.size() - start < lkmer ) return false;
	for( int i = 0 ; i < lkmer ; ++i )
	{
		//cerr << str[i+start] << " " << kmer[i] << endl;
		if( str[i+start] != kmer[i] ) return false;
	}
	return true;
}

vector<bool>& featurize( const string& str, const vector<char>& alphabet, int k, vector<bool>& features)
{
	char* kmer = new char[k];
	int max = pow(alphabet.size(), k);
	int flen = (str.size() - k + 1) * max;
	bool* feats = new bool[ flen ];
	for ( int i = 0 ; i < max ; ++i )
	{
		permutation( alphabet, i, k, kmer);
		for ( int j = 0 ; j < str.size() - k + 1 ; ++j )
		{
			feats[ j * max + i ] = matches( str, k, kmer, j );
			//cerr << str.substr(j,k) << " " << kmer << " " << feats[ j * max + i ] << endl;
		}
	}
	for ( int i = 0 ; i < flen ; ++i )
	{
		features.push_back(feats[i]);
	}
	delete kmer;
	delete feats;
	
	return features;
		
}

vector<bool>& featurize( const string& str, const vector<char>& alphabet, const int* kmers, int nkmers, vector<bool>& features)
{
	for ( int i = 0 ; i < nkmers ; ++i )
	{
		featurize( str, alphabet, kmers[i], features );
	}
	return features;
}

vector<string>& split ( const string& str, const string& delim, vector<string>& tokens )
{
	int pos = 0 ;
	while ( pos < str.size() )
	{
		//cerr << "Pos = " << pos << endl;
		int next = str.find( delim, pos );
		if( next < 0 ) next = str.size();
		//cerr << "Next = " << next << endl;
		string sub = str.substr(pos,next);
		//cerr << sub << endl;
		if (sub.size() > 0 )
		{
			tokens.push_back(sub);
		}
		pos = next + 1;
		//cerr << pos << "  " << str.size() << endl;
	}
	return tokens;
}

string& toUpper( string& str )
{
	for ( int i = 0 ; i < str.size() ; ++i )
	{
		str[i] = toupper( str[i] );
	}
	return str;
}

void printFeatures( const vector<bool>& features, double intensity, ostream& out )
{
	for ( int i = 0 ; i < features.size() ; ++i )
	{
		if ( features[i] )
		{
			out << "1 ";
		}
		else
		{
			out << "0 ";
		}
	}
	out << intensity << endl;
}

void printFeatureNames ( int strlen, const vector<char>& alphabet, const int* kmers, int nkmers, ostream& out )
{
	for ( int i = 0 ; i < nkmers ; ++i )
	{
		int k = kmers[i];
		char* kmer = new char[k];
		//cerr << "Calling pow" << endl;
		int max = pow(alphabet.size(),k);
		//cerr << "Pow done" << endl;
		for ( int j = 0 ; j < strlen - k + 1 ; ++j )
		{
			for ( int v = 0 ; v < max ; ++v )
			{
				permutation( alphabet, v, k, kmer );
				out << "[" << (j+1) << "]" << kmer << " ";
			}
		}
		delete kmer;
	}
	out << "Intensity" << endl;
}

/*
double stod ( const string& str )
{
	double d;
	stringstream s(str);
	s >> d;
	return d;
}
*/

void featurize ( istream& in, ostream& out , const vector<char>& alphabet, const int* kmers, int nkmers, bool printHeader = true )
{
	string line;
	while ( !in.eof() )
	{
		getline(in, line);
		if ( line.size() > 0 ){
			vector<string> tokens;
			//cerr << "Calling split" << endl;
			split( line, " ", tokens );
			if ( tokens.size() < 2 )
			{
				//try splitting around tab
				tokens.clear();
				split( line, "\t", tokens);
				if ( tokens.size() < 2 )
				{
					//now report an error
					cerr << "Error parsing line: " << endl << line << endl;
					continue;
				}
			}
			//cerr << "Calling toUpper" << endl;
			string seq = toUpper(tokens[0]);
			//cerr << "toUpper done" << endl;
			if (printHeader)
			{
				printFeatureNames( seq.size(), alphabet, kmers, nkmers, out );
				printHeader = false;
			}
			vector<bool> features;
			featurize( seq, alphabet, kmers, nkmers, features );
			printFeatures( features, stod(tokens[1]), out);
		}
	}
}

string& filename( const string& filepath, string& name )
{
	string pathdelim;
	if (filepath.find("/") != filepath.npos )
	{
		pathdelim = "/";
	}
	else
	{
		pathdelim = "\\";
	}

	int start = filepath.find_last_of(pathdelim) + 1;
	if ( start >= filepath.npos ) start = 0;
	
	int end = filepath.find(".",start);
	
	name = filepath.substr(start,end-start);
	return name;
	
}

int parseOutputDir( int argc, const char* argv[], string& outputDir, int* ignore )
{
	for ( int i = 1 ; i < argc ; ++i )
	{	
		//cerr << argv[i] << " " << (argv[i] == "-o") << endl;
		if ( strcmp(argv[i],"-o") == 0 )
		{
			++i;
			if ( i >= argc )
			{
				cerr << "Error: -o requires an argument" << endl;
				return 1;
			}
			outputDir = string(argv[i]);
			*ignore = i - 1;
			return 0;
		}
	}
	outputDir = "";
	return 0;
}

typedef struct Args{

    vector<int> kmers;
    ofstream out;
    ifstream in;

} Args;

bool parse( int argc, char* argv[], Args& args ){
    static struct option long_options[] =
    {
        { "input", required_argument, 0, 'i' },
        { "output", required_argument, 0, 'o' },
        { "kmers", required_argument, 0, 'k' },
        { 0, 0, 0, 0}
    };

    int option_index;
    int c;
    while( ( c = getopt_long( argc, argv, "i:o:k:", long_options, &option_index) ) != -1 ){
        switch(c){
            case 'i':
                args.in.open(optarg);
                if( !args.in.is_open() ){
                    cerr << "Unable to open file: " << optarg << endl;
                    return false;
                }
                break;
            case 'o':
                args.out.open(optarg);
                if( !args.out.is_open() ){
                    cerr << "Unable to open file: " << optarg << endl;
                    return false;
                }
                break;
            case 'k':
                --optind;
                while( optind < argc && *argv[optind] != '-' ){
                    args.kmers.push_back( atoi( argv[optind] ) );
                    ++optind;
                }
                break;
        }
    }

    return true;
}

int main( int argc, char* argv[] )
{

    Args args;
    if( !parse( argc, argv, args ) ){
        usage(cerr);
        return 1;
    }

    if( args.kmers.empty() ){
        cerr << "Error: must specify k-mer lengths for featurizing" << endl;
        usage(cerr);
        return 1;
    }

    ostream& out = args.out.is_open() ? args.out : cout;
    istream& in = args.in.is_open() ? args.in : cin;

	const vector<char> alphabet(ALPH, ALPH + ALPH_SIZE);

    featurize( in, out, alphabet, &args.kmers[0], args.kmers.size() );
    
    /*
	for ( int i = 1 ; i < argc ; ++i )
	{
		if ( i != ignore && i != ignore + 1 )
		{
			ifstream in;
			in.open( argv[i] );
			if ( in.good() )
			{
				string name;
				filename( argv[i], name );
				ofstream out;
				string outname = outputDir + name + "_features.txt";
				out.open( outname.c_str() );
				cerr << "Featurizing " << argv[i] << " to " << outname <<  endl;
				featurize( in, out, alphabet, KMERS, KMERS_SIZE );
			}
			else
			{
				cerr << "Unable to open file: " << argv[i] << endl;
			}
		}
	}
    */
}
