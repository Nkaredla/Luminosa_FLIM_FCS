/*
 * IMMUNE_CELL_MIET_WINDOW_ACCUMULATE  Scatter photons into overlapping windows.
 *
 * [contributionCount, contributed, overflow] = ...
 *     immune_cell_MIET_window_accumulate(pixelCube, rowBin, columnBin, ...
 *         timeBin, selectedLookup, dims)
 *
 * dims = [imageHeight, maximumAnchorRow, maximumAnchorColumn, windowHeight,
 *         windowWidth, selectedPixelCount, gateLength]
 *
 * Each source photon is added to every selected upper-left anchor whose
 * window contains it. This is the same reduction the MATLAB path performs
 * with accumarray, but as a direct pointer increment: no index array is
 * materialised and nothing is sorted, so the cost is O(photons * window
 * area) with a very small constant.
 *
 * pixelCube is a uint16 array of [selectedPixelCount 1 gateLength] and is
 * modified IN PLACE. The caller must own it outright - it is created by
 * zeros() inside immune_cell_MIET_reassigned_sliding_tcspc and never shared
 * or aliased - because MATLAB copy-on-write is bypassed here.
 *
 * Saturation is not silently tolerated: if any bin would exceed 65535 the
 * overflow flag is returned and the caller raises the same error the MATLAB
 * path raises.
 *
 * Build:
 *   mex -O immune_cell_MIET_window_accumulate.c
 */

#include "mex.h"
#include <stdint.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:nrhs",
                          "Six inputs are required.");
    }
    if (!mxIsUint16(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:cube",
                          "pixelCube must be a full uint16 array.");
    }
    if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3])) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:bins",
                          "rowBin, columnBin and timeBin must be double.");
    }
    if (!mxIsUint32(prhs[4])) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:lookup",
                          "selectedLookup must be uint32.");
    }
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 7) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:dims",
                          "dims must be a double vector with seven entries.");
    }

    uint16_t *cube = (uint16_t *) mxGetData(prhs[0]);
    const double *rowBin = mxGetPr(prhs[1]);
    const double *columnBin = mxGetPr(prhs[2]);
    const double *timeBin = mxGetPr(prhs[3]);
    const uint32_t *lookup = (const uint32_t *) mxGetData(prhs[4]);
    const double *dims = mxGetPr(prhs[5]);

    const mwSize photonCount = mxGetNumberOfElements(prhs[1]);
    if (mxGetNumberOfElements(prhs[2]) != photonCount ||
        mxGetNumberOfElements(prhs[3]) != photonCount) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:binLength",
                          "rowBin, columnBin and timeBin must be the same length.");
    }

    const long imageHeight        = (long) dims[0];
    const long maximumAnchorRow   = (long) dims[1];
    const long maximumAnchorCol   = (long) dims[2];
    const long windowHeight       = (long) dims[3];
    const long windowWidth        = (long) dims[4];
    const long selectedPixelCount = (long) dims[5];
    const long gateLength         = (long) dims[6];

    const mwSize lookupCount = mxGetNumberOfElements(prhs[4]);
    const mwSize cubeCount = mxGetNumberOfElements(prhs[0]);
    if (cubeCount != (mwSize)(selectedPixelCount * gateLength)) {
        mexErrMsgIdAndTxt("immune_cell_MIET_window_accumulate:cubeSize",
                          "pixelCube size does not match selectedPixelCount*gateLength.");
    }

    plhs[0] = mxCreateDoubleScalar(0.0);
    plhs[1] = mxCreateLogicalMatrix(photonCount, 1);
    plhs[2] = mxCreateLogicalScalar(false);
    double *contributionCount = mxGetPr(plhs[0]);
    mxLogical *contributed = mxGetLogicals(plhs[1]);
    mxLogical *overflow = mxGetLogicals(plhs[2]);

    double contributions = 0.0;
    int sawOverflow = 0;

    for (mwSize i = 0; i < photonCount; ++i) {
        const long sourceRow = (long) rowBin[i];
        const long sourceCol = (long) columnBin[i];
        const long t = (long) timeBin[i];
        if (t < 1 || t > gateLength) {
            continue;
        }
        const long timeOffset = (t - 1) * selectedPixelCount;
        int any = 0;

        for (long rowOffset = 0; rowOffset < windowHeight; ++rowOffset) {
            const long anchorRow = sourceRow - rowOffset;
            if (anchorRow < 1 || anchorRow > maximumAnchorRow) {
                continue;
            }
            for (long colOffset = 0; colOffset < windowWidth; ++colOffset) {
                const long anchorCol = sourceCol - colOffset;
                if (anchorCol < 1 || anchorCol > maximumAnchorCol) {
                    continue;
                }
                /* 1-based full-image linear index, as in the MATLAB path. */
                const long full = anchorRow + (anchorCol - 1) * imageHeight;
                if (full < 1 || (mwSize) full > lookupCount) {
                    continue;
                }
                const uint32_t compactPixel = lookup[full - 1];
                if (compactPixel == 0u) {
                    continue;
                }
                const long index = (long) compactPixel - 1 + timeOffset;
                if (index < 0 || (mwSize) index >= cubeCount) {
                    continue;
                }
                if (cube[index] == UINT16_MAX) {
                    sawOverflow = 1;
                } else {
                    cube[index] += 1u;
                }
                contributions += 1.0;
                any = 1;
            }
        }
        if (any) {
            contributed[i] = 1;
        }
    }

    *contributionCount = contributions;
    *overflow = sawOverflow ? 1 : 0;
}
