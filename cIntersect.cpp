// mex
#include "mex.h"

// stdlib
#include <math.h>
#include <vector>

// our headers
#include "mexUtil.h"

using namespace std;


// -- Type definitions --
typedef Array1D_t<double> Array1D;
typedef Array2D_t<double> Array2D;
typedef Array3D_t<double> Array3D;


// -- Prototypes -- //

/* Intersection of vectors y and y+lag
    Input: y - array of input values
 *         lag - lag value
 *         values_out - intersecting values
 *         idx1_out - indices of intersections for array y
 *         idx2_out - indices of intersections for array y+lag
 */
void intersect(Array1D& y, double lag, std::vector<double>& values_out, std::vector<unsigned int>& idx1_out, std::vector<unsigned int>& idx2_out);

// Fuzzy comparison needed for doubles
// return a<b    (fuzzy: a-b<-1e-10 )
bool fuzzyLess(double a, double b);


// --  Global Variables  -- //
mxArray** plhs;  // Access to left hand side is global, so we can output to it in functions



/* Usage:  [values idx1 idx2] = cIntersect( y, lag)
* Compute the intersection of arrays (y) and (y+lag)
* Input:
*       y - 1D array
*       lag - lag value
* Ouptut:
*          values - intersecting values
*          idx1 - indices of intersections for array y
*          idx2 - indices of intersections for array y+lag */
void mexFunction(int nlhs, mxArray* plhs_arg[], int nrhs, const mxArray* prhs[])
{	    
    // Access to left hand side is global, so we can output to it in functions
    plhs = plhs_arg;

	 /* Check for proper number of arguments. */
 	if(nrhs!=2) {
    	mexErrMsgTxt("Two inputs required.");
  	} 
    if(nlhs!=3) {
    	mexErrMsgTxt("Three outputs required.");
  	}

	/* Map access to input data */
	Array1D y( prhs[0] );
    double lag = *mxGetPr(prhs[1]);
    
    // Set up containers for output (allocate memory with guess)
    std:vector<double> values;        values.reserve(y.nElements/4);
    std::vector<unsigned int> idx1;   idx1.reserve(y.nElements/4);
    std::vector<unsigned int> idx2;   idx2.reserve(y.nElements/4);

    
    /* ---- Computing the intersection ---- */
    intersect(y, lag, values, idx1, idx2);
    unsigned int num_intersections = values.size();
    
            
    // -- Copy Result to Matlab output -- //  
     mwSize dim_out[1] = { num_intersections };
     plhs[0] = mxCreateNumericArray( 1, dim_out , mxDOUBLE_CLASS, mxREAL);
     plhs[1] = mxCreateNumericArray( 1, dim_out , mxDOUBLE_CLASS, mxREAL);
     plhs[2] = mxCreateNumericArray( 1, dim_out , mxDOUBLE_CLASS, mxREAL);
     
     // Map access to output
     Array1D mx_val_out ( plhs[0] );
     Array1D mx_idx1_out( plhs[1] );
     Array1D mx_idx2_out( plhs[2] );
     
     // Copy to output
     for(unsigned int i=0; i< num_intersections; ++i)
     {
         mx_val_out[i]  = values[i];
         mx_idx1_out[i] = idx1[i];
         mx_idx2_out[i] = idx2[i];
     }
}


// Fuzzy comparison needed for doubles
// return a<b    (fuzzy: a-b<-1e-10 )
bool fuzzyLess(double a, double b)
{
  if ( (a-b) > -1e-10)        
      return false;
  
  return true;
}


/* Intersection of vectors y and y+lag
    Input: y - array of input values
 *         lag - lag value
 *  Output:
 *         values_out - intersecting values
 *         idx1_out - indices of intersections for array y
 *         idx2_out - indices of intersections for array y+lag
 */
void intersect(Array1D& y, double lag, std::vector<double>& values_out, std::vector<unsigned int>& idx1_out, std::vector<unsigned int>& idx2_out)
{        
    unsigned int pos0 = 0;
    unsigned int pos1 = 0;
    
    // Search for elements in y equal to their lagged value
     while(pos0 < (unsigned int) y.nElements  && pos1 < (unsigned int) y.nElements)
     {
         if( fuzzyLess( y[pos1]+lag, y[pos0]) ) ++pos1;
         else if( fuzzyLess( y[pos0], y[pos1]+lag) ) ++pos0;
         else {
              values_out.push_back( y[pos0] );
              idx1_out.push_back( pos0+1 );
              idx2_out.push_back( pos1+1 );   
              ++pos0;
              ++pos1;
         }
     }    
}


