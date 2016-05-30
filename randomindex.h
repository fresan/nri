/*
 * Random indexing of n-dimensional arrays, v 1.0
 *
 * URL: github.com/fresan/nri (released under The MIT License).
 *
 * Copyright Fredrik Sandin (2010-2016), fredrik.sandin@gmail.com
 *
 * Requires the Mersenne Twister algorithm implemented in mtrand.h
 * at http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 */

#ifndef RANDOMINDEX_H
#define RANDOMINDEX_H

#include <cassert>
#include "mtrand.h"

#define __HALF_MAX_SIGNED(type) ((type)1 << (sizeof(type)*8-2))
#define __MAX_SIGNED(type) (__HALF_MAX_SIGNED(type) - 1 + __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type) (-1 - __MAX_SIGNED(type))
#define __MIN(type) ((type)-1 < 1?__MIN_SIGNED(type):(type)0)
#define __MAX(type) ((type)~__MIN(type))

/* Random indexing data types */
typedef signed short    RIData;  // Must be signed
typedef unsigned short  RIIndex; // Sign not needed

class RandomIndex {
public:
    
    RandomIndex(unsigned dims, unsigned *datarange, unsigned *numrind);
    virtual ~RandomIndex();
    
    //! Number of random indices and number of dimensions of distributional array
    unsigned dims(void) const;
    
    //! Get range of random index
    unsigned range(unsigned dim) const;
    
    //! Set and return range of random index (does not decrease the range)
    unsigned setrange(unsigned dim, unsigned range);
    
    //! Number of random indices used for distributional coding
    unsigned numrind(unsigned dim) const;
    
    //! Range of distributional array index in each dimension
    unsigned datarange(unsigned dim) const;
    
    //! Size of distributional representation in bytes
    unsigned datasize(void) const;
    
    //! Size of random index tables in bytes
    unsigned indexsize(void) const;
    
    //! Encode weight at location ind[] in array
    void encode(RIData *data, const unsigned *ind, RIData weight);
    
    //! Decode approximated data at location ind[] in array
    double decode(const RIData *data, const unsigned *ind) const;
    
    //! Cos(alpha) estimation for two distributional representations
    double cosa(const RIData *d1, const unsigned *i1, const RIData *d2, const unsigned *i2) const;
    
    //! Read saturation counter
    unsigned long saturation(void) const;
    
    //! Averaging directive for cosa()
    enum { Average = __MAX(unsigned) };
    
protected:
    
    /**
     * Encapsulates functionality needed for random indexing of one dimension.
     * Multiple RITable objects are used for RI of multidimensional arrays.
     */
    class RITable {
    public:
        RITable(unsigned datarange, unsigned cols);
        virtual ~RITable();
        unsigned rows(void) const;
        unsigned cols(void) const;
        unsigned numpositive(void) const;
        unsigned size(void) const;
        unsigned datarange(void) const;
        RIIndex& item(unsigned row, unsigned col) const;
        unsigned rows(unsigned n);
    protected:
        unsigned datarange_;
        unsigned rows_;
        unsigned cols_;
        RIIndex *table_;
    };
    
    unsigned dims_;         // Dimensionality, number or random index tables
    unsigned datanumel_;    // Number of elements in distributional array
    unsigned long datasat_; // Indicates how many times data elements have been saturated
    unsigned distnumel_;    // Size of distributed representation (number of elements)
    unsigned *riunroll_;    // Index transformation vector for loop unrolling
    RITable **ritables_;    // Array of pointers to random index tables
    
    static MTRand mtrand_;  // Pseudorandom number generator
};

inline unsigned RandomIndex::RITable::rows(void) const {
    return rows_;
}

inline unsigned RandomIndex::RITable::cols(void) const {
    return cols_;
}

inline unsigned RandomIndex::RITable::numpositive(void) const {
    return cols_>>1;
}

inline unsigned RandomIndex::RITable::size(void) const {
    return sizeof(RIIndex)*cols_*rows_;
}

inline unsigned RandomIndex::RITable::datarange(void) const {
    return datarange_;
}

inline RIIndex& RandomIndex::RITable::item(unsigned row, unsigned col) const {
    assert(row < rows_);
    assert(col < cols_);
    assert(table_[row*cols_ + col] < datarange_);
    return table_[row*cols_ + col];
}

inline unsigned RandomIndex::dims(void) const {
    return dims_;
}

inline unsigned RandomIndex::range(unsigned dim) const {
    assert(dim < dims_);
    return ritables_[dim]->rows();
}

inline unsigned RandomIndex::numrind(unsigned dim) const {
    assert(dim < dims_);
    return ritables_[dim]->cols();
}

inline unsigned RandomIndex::datarange(unsigned dim) const {
    assert(dim < dims_);
    return ritables_[dim]->datarange();
}

inline unsigned RandomIndex::datasize(void) const {
    return datanumel_*sizeof(RIData);
}

#endif // RANDOMINDEX_H
