
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Histogram_h
#define MEDYAN_Histogram_h

#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>

#include "common.h"

namespace medyan {
/// A class to hold frequency of occurences in a general set of data.

/*!
 *  The Histogram class holds frequencies of events in
 *  certain ranges (bins). It contains functions to
 *  increment and store frequencies in a data structure.
 */
class Histogram {
   
private:
    int _numBins; ///< Number of bins in histogram
    
    floatingpoint _histMin; ///< Minimum histogram value
    floatingpoint _histMax; ///< Maximum histogram value
    
    floatingpoint _range; ///< range of bin

    vector<int> _frequencies; ///< The histogram data
    
public:
    Histogram(int numBins, floatingpoint histMin = 0, floatingpoint histMax = 0)
        : _numBins(numBins), _histMin(histMin), _histMax(histMax),
          _frequencies(numBins, 0) {

        //calculate range
        _range = (_histMax + _histMin) / _numBins;
    }
    
    ///Adds a value to the histogram by updating the correct bin
    void addValue(floatingpoint value) {
        
        int bin = max(0, (int) (value / _range));
        assert(bin >= 0 && bin < _numBins && "Histogram error - trying to add value outside of bin range.");
        
        _frequencies[bin] += 1;
    }
    
    ///Clears all values from histogram
    void clearValues() {
        _frequencies.assign(_numBins, 0);
    }
    
    ///Print the histogram
    void print(ofstream& outputFile) {
        
        int bin = 0;
        for(auto freq : _frequencies) {
            
            float binMin = bin * _range;
            float binMax = (bin + 1) * _range;
            
            outputFile << binMin << " " << binMax << " " << freq << " ";
            
            bin++;
        }
    }
    floatingpoint getMin() {return _histMin;}
    floatingpoint getMax() {return _histMax;}
    
};

} // namespace medyan

#endif
