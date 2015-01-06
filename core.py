import sys
import os

def core( bp, fpath ):
    root, ext = os.path.splitext( fpath )
    wname = root + '_' + str(bp) + 'bpcore' + ext;
    with open( wname, 'w' ) as out:
        with open( fpath, 'r' ) as fin:
            for line in fin:
                [ seq, inten ] = line.strip().split()
                start = ( len( seq ) - bp ) / 2
                out.write( seq[start:start+bp] + " " + inten + "\n" )


def main( args ):
    bp = int( args[0] )
    for arg in args[1:]:
        core( bp, arg )

if __name__ == "__main__":
    main( sys.argv[1:] )
