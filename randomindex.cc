/*
 * Random indexing of n-dimensional arrays, v 1.0
 *
 * URL: github.com/fresan/nri (released under The MIT License).
 *
 * Copyright Fredrik Sandin (2010-2016), fredrik.sandin@gmail.com
 */

#include <iostream>
using namespace std;

#include <cstdlib>
#include "randomindex.h"

MTRand RandomIndex::mtrand_(0x12345);   // Constant seed so that experiments can be repeated

RandomIndex::RandomIndex(unsigned dims, unsigned *datarange, unsigned *numrind)
: dims_(dims), datasat_(0), datanumel_(1), distnumel_(1), ritables_(NULL), riunroll_(NULL) {
    
    assert(__MIN(RIIndex) == 0);        // Better be unsigned
	assert(__MIN(RIData) < 0);          // Must be signed
	assert(dims > 0);
    
    unsigned i,j;
    ritables_ = new RITable* [dims];
    riunroll_ = new unsigned [dims];
    for(i=0; i<dims; i++) {
        ritables_[i] = new RITable(datarange[i], numrind[i]);
        datanumel_ *= datarange[i];
        distnumel_ *= numrind[i];
        riunroll_[i] = 1;
        
        // Lookup table for loop unrolling, row-major ordering
        for(j=i+1; j<dims_; j++) {
            riunroll_[i] *= numrind[j];
        }
    }
}

RandomIndex::~RandomIndex() {
    if(ritables_ != NULL) {
        for(unsigned i=0; i<dims(); i++) {
            delete ritables_[i];
        }
        delete[] ritables_;
        delete[] riunroll_;
    }
}

unsigned RandomIndex::setrange(unsigned dim, unsigned range) {
    assert(dim < dims_);
    return ritables_[dim]->rows(range);
}

unsigned RandomIndex::indexsize(void) const {
    unsigned i, s=0;
    for(i=0; i<dims_; i++) {
        s += ritables_[i]->size();
    }
    return s;
}

void RandomIndex::encode(RIData *data, const unsigned *ind, RIData weight) {
    if(weight==0) return; // Nothing to do, return early
    unsigned i, *coli = new unsigned[dims_];
    for(unsigned j=0; j<distnumel_; j++) { // Unrolled loop, row-major ordering
        RIData w = weight;
		unsigned temp=j, dataindex=0;
        for(i=0; i<dims_; i++) {
			coli[i] = temp / riunroll_[i];
			temp -= coli[i] * riunroll_[i];
			if(coli[i] >= ritables_[i]->numpositive()) w = -w;
			coli[i] = ritables_[i]->item(ind[i], coli[i]);
            dataindex = coli[i] + ritables_[i]->datarange()*dataindex;
        }
        long v = data[dataindex] + w;
        if(v >= __MIN(RIData) && v <= __MAX(RIData)) {
            data[dataindex] = v;
        } else { // Saturation logic
            datasat_ ++;
            if(v < __MIN(RIData)) {
                data[dataindex] = __MIN(RIData);
            } else {
                data[dataindex] = __MAX(RIData);
            }
        }
    }
  	delete[] coli;
}

double RandomIndex::decode(const RIData *data, const unsigned *ind) const {
    double weight = 0;
    unsigned i, *coli = new unsigned[dims_];
	for(unsigned j=0; j<distnumel_; j++) { // Unrolled loop, row-major ordering
        RIData s = 1;
		unsigned temp=j, dataindex=0;
        for(i=0; i<dims_; i++) {
			coli[i] = temp / riunroll_[i];
			temp -= coli[i] * riunroll_[i];
			if(coli[i] >= ritables_[i]->numpositive()) s = -s;
			coli[i] = ritables_[i]->item(ind[i], coli[i]);
            dataindex = coli[i] + ritables_[i]->datarange()*dataindex;
        }
        weight += s*data[dataindex];
    }
  	delete[] coli;
    return weight / distnumel_;
}

double RandomIndex::cosa(const RIData *d1, const unsigned *i1,
                         const RIData *d2, const unsigned *i2) const {
    unsigned i, j, dnumel = 1;
    unsigned *coli1 = new unsigned[dims_];
    unsigned *coli2 = new unsigned[dims_];
    unsigned *unroll = new unsigned[dims_];
    double weight = 0, norm1 = 0, norm2 = 0;
    for(i=0; i<dims_; i++) {
        if(i1[i]==Average || i2[i]==Average) {
            assert(i1[i] == i2[i]);
            dnumel *= ritables_[i]->datarange();
        } else {
            dnumel *= ritables_[i]->cols();
        }
        // Unrolling lookup table, row-major order
        unroll[i] = 1;
        for(j=i+1; j<dims_; j++) {
            if(i1[j] == Average) {
                unroll[i] *= ritables_[j]->datarange();
            } else {
                unroll[i] *= ritables_[j]->cols();
            }
        }
    }
    for(j=0; j<dnumel; j++) { // Unrolled loop, row-major order
        unsigned temp=j, dind1=0, dind2=0;
        for(i=0; i<dims_; i++) {
            coli1[i] = coli2[i] = temp / unroll[i];
            temp -= coli1[i] * unroll[i];
            if(i1[i] != Average) { // Get random index
                coli1[i] = ritables_[i]->item(i1[i], coli1[i]);
                coli2[i] = ritables_[i]->item(i2[i], coli2[i]);
            }
            dind1 = coli1[i] + ritables_[i]->datarange()*dind1;
            dind2 = coli2[i] + ritables_[i]->datarange()*dind2;
        }
        weight += d1[dind1]*d2[dind2]; // RI signs cancel
        norm1 += d1[dind1]*d1[dind1];
        norm2 += d2[dind2]*d2[dind2];
    }
    delete[] coli1;
    delete[] coli2;
    delete[] unroll;
    return weight / sqrt(norm1*norm2);
}

unsigned long RandomIndex::saturation(void) const {
    return datasat_;
}

RandomIndex::RITable::RITable(unsigned datarange, unsigned cols)
: datarange_(datarange), cols_(cols), rows_(0), table_(NULL) {
    assert(datarange_ <= __MAX(RIIndex));   // Range of RIIndex type
    assert(!(cols_ & 0x1));                 // Even number of cols
    assert(cols_ < datarange_);             // Necessary for unique indices
}

RandomIndex::RITable::~RITable() {
    if(table_ != NULL) {
        free(table_);
        table_ = NULL;
    }
}

unsigned RandomIndex::RITable::rows(unsigned n) {
    unsigned i0=rows_, i, j, k, upper=datarange_-1;
    if(n > rows_) {                     // Add rows to table?
        if(table_ != NULL) {
            table_ = (RIIndex *)realloc(table_, n*cols_*sizeof(RIIndex));
        } else {
            table_ = (RIIndex *)malloc(n*cols_*sizeof(RIIndex));
        }
        assert(table_ != NULL);
        rows_ = n;
        for(i=i0; i<rows_; i++) {       // For each new row
            for(j=0; j<cols_; j++) {    // ... and for each random index
                RIIndex rnd;
                while(true) {           // Generate unique index in [0,upper]
                    rnd = RandomIndex::mtrand_.randInt(upper);
                    for(k=0; k<j; k++) {
                        if(item(i,k) == rnd) k=cols_;
                    }
                    if(k<cols_) break;
                }                       // Linear search, short list
                item(i,j) = rnd;
            }
        }
    }
    return rows_;
}

/*

Changelog

v1.0 First release on http://github.com/fresan/nri

v0.4.1
- Corrected a bug in RandomIndex::RITable::rows()
 
v0.4.0
- Rewrite of code.
  A RandomIndex class is defined that only specifies the random indexing
  tables and methods needed to encode, decode and approximate inner products
  given data arrays that are defined by the user. This way it becomes possible
  to combine explicit and random indexing without hardcoding of special cases.
  This design simplifies the code and makes it more flexible and easy to use.
- File I/O has not been ported/implemented in this version, nor templates
  or the Matlab/MEX interface.
 
v0.3.1
 - Updated an assert statement in stateSize() that was incompatible
 with one-way indexing.

v0.3.0
- Made randinit() compatible with one-way indexing. This special
  case was omitted in the initial implementation.
 
Early developments not documented.
 
EOF */
