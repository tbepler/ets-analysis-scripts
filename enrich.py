import sys
import Enrichment


def parse( rin ):
    weights = []
    labels = []
    for line in rin:
        [ lab, w ] = line.strip().split()
        lab = int( lab )
        w = float( w )
        if lab == 0:
            lab = -1
        weights.append( w )
        labels.append( lab )
    return ( weights, labels )

def main( args ):
    if len( args ) == 0:
        rin = sys.stdin
    else:
        rin = open( args[0], 'r' )
    ( weights, labels ) = parse( rin )
    rin.close()
    ( s, es, scores ) = Enrichment.enrichment( weights, labels )
    print s, es


if __name__ == '__main__':
    main( sys.argv[1:] )
