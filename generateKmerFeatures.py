import subprocess
import sys
import os

KMERS = [ "1", "1 2", "1 2 3" ]
NAMES = [ "1mer", "1-2mer", "1-3mer" ]
FEATURIZE = "featurize -k "

def featurize( path ):
    fname = os.path.basename( path )
    root, ext = os.path.splitext( fname )
    for ( kmer, name ) in zip( KMERS, NAMES ):
        ofname = root + "_" + name + ext
        cmd = FEATURIZE + kmer + " < " + path + " > " + ofname
        print >> sys.stderr, "Executing: " + cmd
        subprocess.call( [ "featurize", "-k" ] + kmer.split(), stdout = open( ofname, 'w' ), stdin = open( path, 'r' ) )

def main( args ):
    for arg in args:
        featurize( arg )

if __name__ == "__main__":
    main( sys.argv[1:] )
