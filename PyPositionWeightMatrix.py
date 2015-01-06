import PositionWeightMatrix as PWM
import math
import random
import multiprocessing
import itertools
from collections import Counter


class PositionWeightMatrix( PWM.PositionWeightMatrix ):

    def __init__(self,scores,alphabet):
        super( PositionWeightMatrix, self ).__init__(scores)
        self._alphabet = alphabet
    
    def __eq__(self,other):
        #TODO update because broken
        return pwmEq(self._struct, other._struct, 0)

    def __hash__(self):
        #TODO update and stuff
        h = 5
        for c in self._alphabet:
            for j in range(len(self._struct[i])):
                h = h*11 + hash(self._struct[c][j])
        return h
    
    def __reduce__(self):
        return  self.__class__, (self.toList(), self.alphabet())

    def occurs( self, seq, cutoff = 0 ):
        assert len(seq) >= self.length()
        for i in range( len(seq) - self.length() + 1 ):
            if self.score( seq[i:i+self.length()] ) > cutoff:
                return True
        return False

    def occurrences( self, seq, cutoff = 0 ):
        return sum( 1 for s in self.scoreAll( seq ) if s > cutoff )

    def alphabetSize(self):
        return len( self._alphabet )
    
    def alphabet(self):
        return self._alphabet

    def adjustToBG( self, bg ):
        return adjustToBG( self, bg )

    def translate( self, seq ):
        return [ self._alphabet[c] for c in seq ]

    def __str__(self):
        s = ' '.join(self._alphabet)
        for row in self:
            s += '\n' + ' '.join( map( str, row ) )
        return s
        
#TODO BROKE
def pwmEq(pwm1,pwm2,tolerance):
    for i in pwm1.keys():
        if len(pwm1[i]) == len(pwm2[i]):
            for j in range(len(pwm1[i])):
                if abs(exp(pwm1[i][j]) - exp(pwm2[i][j])) > tolerance: return False
        else:
            return False
    return True

def adjustToBG(pwm,bg):
    pwm = PositionWeightMatrix(
        [ [ pwm[i][j] - log(bg[j]) for j in range( pwm.alphabetSize() ) ] for i in range( len( pwm ) ) ],
        pwm.alphabet() 
        )
    return pwm

def log(x):
    if x == 0: return float("-inf")
    return math.log(x)

def parseBG( bgF ):
    length = 0
    counts = Counter()
    for line in bgF:
        if line[0] != ">":
            length += len(line)
            counts.update(line.upper().strip())
    return [ float(counts[x])/float(length) for x in alphabet ]

def nucleotideFrequency( seq, pseudocount = 0 ):
    elems = max( seq )
    bg = [ 0 for _ in xrange( elems + 1 ) ]
    for i in seq:
        bg[i] += 1
    addedCounts = len( bg ) * pseudocount;
    denom = float( len(seq) ) + addedCounts
    for i in xrange( len( bg ) ):
        bg[i] = float( bg[i] + pseudocount ) / denom
    return bg

def nullScoreDistribution( pwm, seq, reps = 1000, tPool = None ):
    return pwm.sampleNullDistribution( seq, reps )

def scoreCutoff( nullDist, pval ):
    i = float( len( nullDist ) ) -  pval * len( nullDist )
    upper = int( math.ceil( i ) ) 
    lower = int( math.floor( i ) )
    if upper == lower:
        return nullDist[upper]
    wl = float( upper ) - i
    wu = i - float( lower )
    return wu * nullDist[upper] + wl * nullDist[lower]

def parsePFM( f, pseudocount = 0 ):
    bases = f.readline().strip().split()
    bases = { c : i for ( i, c ) in enumerate( bases )  }
    pfm = [ [ float(s)+pseudocount for s in line.strip().split() ] for line in f ]
    pfm = [ [ log(elem) - log( sum( row ) ) for elem in row ] for row in pfm ]
    return PositionWeightMatrix(pfm,bases)
