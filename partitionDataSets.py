import sys
import random
import os

N_TEST = 5000

def partition( path ):
    lines = open( path, 'r' ).readlines()
    random.shuffle( lines )
    fname = os.path.basename( path )
    root, ext = os.path.splitext( fname )
    testName = root + "_test" + ext
    trainName = root + "_train" + ext
    open( testName, 'w' ).writelines( lines[0:N_TEST] )
    open( trainName, 'w' ).writelines( lines[N_TEST:] )
    

def main( args ):
    for arg in args:
        partition( arg )






if __name__ == '__main__':
    main( sys.argv[1:] )
