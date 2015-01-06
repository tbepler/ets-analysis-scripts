import PyPositionWeightMatrix as PWM
import sys
import os

comp = { 'A' : 'T', 'C' : 'G', 'G' : 'C' , 'T' : 'A' }

def rvscomp( seq ):
    return [ comp[c] for c in seq[::-1] ]

def scoreSeqs( pwm, fpath ):
    with open( fpath, 'r' ) as f:
        for line in f:
            [ seq, _, label ] = line.strip().split()
            seq = seq.upper()
            seqValid = True
            for base in seq:
                if base != 'A' and base != 'C' and base != 'G' and base != 'T':
                    seqValid = False
                    break
            if not seqValid:
                continue
            rvs = rvscomp( seq )
            fwdScore = max( pwm.scoreAll( pwm.translate( seq ) ) )
            rvsScore = max( pwm.scoreAll( pwm.translate( rvs ) ) )
            
            print label, max( [ fwdScore, rvsScore ] )
    return


def main( args ):
    f = open( args[0], 'r' )
    pwm = PWM.parsePFM( f )
    f.close()
    scoreSeqs( pwm, args[1] )


if __name__ == "__main__":
    main( sys.argv[1:] )
