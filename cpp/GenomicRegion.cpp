#include<stdio.h>
#include<iostream>
#include"GenomicRegion.h"

using namespace std;

typedef GenomicRegion GR;

GR::GenomicRegion( const Chromosome& chrom, const unsigned long start, const unsigned long end )
    :mStart( start ), mEnd( end ), mChrom( chrom ) {}


