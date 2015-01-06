import sys
import os

class Model:

    def __init__( self, bias, features, length ):
        self.bias = bias
        self.features = features
        self.length = length

    def score( self, seq ):
        s = self.bias
        for ( i, kmer, w ) in self.features:
            if seq[i:i+len(kmer)] == kmer:
                s += w
        return s

    def __len__( self ):
        return self.length

def readModel( fpath ):
    features = []
    length = 0
    bias = 0
    with open( fpath, 'r' ) as f:
        for line in f:
            if line[0:4] == 'Bias':
                bias = float( line.strip().split()[1] )
            else:
                line = line.replace( '[', ' ' )
                line =line.replace( ']', ' ' )
                [ i, kmer, w ] = line.strip().split()
                i = int(i) - 1
                w = float(w)
                if w > 0:
                    features.append( ( i, kmer, w ) )
                if i + len( kmer ) > length:
                    length = i + len( kmer )
    return Model( bias, features, length )

def scoreSeqs( M, fpath ):
    with open( fpath, 'r' ) as f:
        for line in f:
            [ seq, _, label ] = line.strip().split()
            seq = seq.upper()
            for base in seq:
                if base != 'A' or base != 'C' or base != 'G' or base != 'T':
                    continue
            smax = 0
            for i in range( len( seq ) - len( M ) + 1 ):
                s = M.score( seq[i:i+len(M)] )
                if s > smax:
                    smax = s
            print label, smax
    return


def main( args ):
    M = readModel( args[0] )
    scoreSeqs( M, args[1] )


if __name__ == "__main__":
    main( sys.argv[1:] )
