/*
 * Random indexing of n-dimensional arrays, v 1.0
 *
 * URL: github.com/fresan/nri (released under The MIT License).
 *
 * Copyright Fredrik Sandin (2010-2016), fredrik.sandin@gmail.com
 */

#include "randomindex.h"

#include <iostream>
using namespace std;

int main(void) {

    /*
     * In this example the concept of distributional vectors
     * representing terms is extended to 2D distributional arrays,
     * which for example enables analysis of context dependence.
     */

    // Distributional datatype with 64k elements
    const unsigned D1 = 2048;       // Co-occurrence dimensions
    const unsigned D2 = 32;         // Context dimensions
    typedef RIData dtype[D1*D2];
    
    // Set of 1k distributional arrays
    const unsigned D3 = 1024;
    dtype *terms = new dtype[D3];
    
    cout << "Created " << D3
        << " distributional arrays of size "
        << D1 << "x" << D2 << endl;
    
    // Define 2D random index for dtype
    unsigned dr[] = {D1,D2};
    unsigned nr[] = {8,4};
    RandomIndex ri(2, dr, nr);
    
    // Check that parameters are consistent
    assert(ri.datasize() == sizeof(dtype));
    
    // Set range of first random index to [0,9999]
    ri.setrange(0, 10000);
    
    // Set range of second random index to [0,999]
    ri.setrange(1, 1000);
    
    cout << "Encoding 1M random co-occurence weights..." << endl;

    unsigned l;
    MTRand rnd(0x2345);
    for(l=0; l<1000000; l++) {
        
        unsigned i = (unsigned)rnd.randInt(D3-1); // Random term
        unsigned j = (unsigned)rnd.randInt(9999); // Random co-occurrence
        unsigned k = (unsigned)rnd.randInt(999);  // Random context
        unsigned w = (unsigned)rnd.randInt(10);   // Random weight

        unsigned rind[] = {j,k};
        ri.encode(terms[i], rind, w);
    }
    
    // Identify max cos(angle) between one term and all others
    // by averaging over all random indices, which implies
    // averaging over co-occurrences and contexts
   
    double cmax = -1;
    unsigned lmax = 0;
    for(l=1; l<D3; l++) {
        unsigned ind[] = {RandomIndex::Average, RandomIndex::Average};
        double cosa = ri.cosa(terms[0],ind,terms[l],ind);
        if(cosa > cmax) {
            lmax = l;
            cmax = cosa;
        }
    }
    cout << "Term 0 is most similar to term "
        << lmax << " with cos(alpha) " << cmax << endl;

    cout << "Maximizing the context-specific cos(alpha) for these terms..." << endl;
    
    cmax = -1;
    unsigned lcmax = 0;
    for(l=0; l<999; l++) {
        unsigned ind[] = {RandomIndex::Average, l};
        double cosa = ri.cosa(terms[0],ind,terms[lmax],ind);
        if(cosa > cmax) {
            lcmax = l;
            cmax = cosa;
        }
    }

    cout << "Term 0 is most similar to term "
        << lmax << " in context " << lcmax
        << " with cos(alpha) " << cmax << endl;
    
    delete[] terms;
    return 0;
}

/* EOF */