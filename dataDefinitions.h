#pragma once
#include "mat.h"
#include <array>
#include <vector>
#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <math.h>
#include <string>
#include "mex.h"
#include "matrix.h"
#include <stdint.h>
#include <ctime>
#include <algorithm>
#include "vectorfunctions.h"
#include <random>
#include <iterator>
#include <functional>
#include <valarray>

#define NR_FIELDS_AGG (sizeof(field_names_sAggregate)/sizeof(*field_names_sAggregate)) // number of fields in structure sAggregate

// define namespace
using namespace std;

//// define data types
// define fixed-size matrix class
class Matrix
{
private:
	int64_t m_SizeCols;
	int64_t m_SizeRows;
	int64_t m_NumElements;
	vector<double> m_Data;

public:
	// default constructor
	Matrix()
		: m_SizeCols(0), m_SizeRows(0), m_NumElements(0)
	{
		;
	}
	
	// constructor of empty matrix; inputs: columns, rows
	Matrix(int64_t intSizeCols, int64_t intSizeRows)
		: m_SizeCols(intSizeCols), m_SizeRows(intSizeRows), m_NumElements(intSizeCols * intSizeRows)
	{
		m_Data.resize(m_NumElements); // allocate matrix
	}
	// constructor of empty matrix; inputs: ptrMatlabVector
	Matrix(int64_t intSizeCols, int64_t intSizeRows, const mxArray *ptrMatlabVector)
		: m_SizeCols(intSizeCols), m_SizeRows(intSizeRows), m_NumElements(intSizeCols * intSizeRows)
	{
		m_Data.resize(m_NumElements); // allocate matrix
		fillVals(ptrMatlabVector); // fill with values
	}

	// destructor
	~Matrix()
	{
		m_Data.clear(); // Clear the data
		m_SizeCols = 0;
		m_SizeRows = 0;
		m_NumElements = 0;
	}

	// functions for setting size; NB: this clears all data!
	void setSizeCols(int64_t intSizeCols)
	{
		m_SizeCols = intSizeCols;
		m_NumElements = intSizeCols * m_SizeRows;
		m_Data.resize(m_NumElements); // allocate matrix
	}
	void setSizeRows(int64_t intSizeRows)
	{
		m_SizeRows = intSizeRows; 
		m_NumElements = intSizeRows * m_SizeCols;
		m_Data.resize(m_NumElements); // allocate matrix
	}
	
	// functions for returning size
	int64_t getSizeCols() { return m_SizeCols; }
	int64_t getSizeRows() { return m_SizeRows; }
	int64_t getSize() { return m_NumElements; }

	// retrieve value from matrix
	double getVal(const int64_t &intCol, const int64_t &intRow)
	{
		return m_Data[intRow * m_SizeCols + intCol];
	}
	double getValAtIdx(const int64_t &intIdx)
	{
		return m_Data[intIdx];
	}
	vector<double> getAllVals()
	{
		return m_Data;
	}

	// set value to matrix
	void setVal(const int64_t &intCol, const int64_t &intRow, const double &dblVal)
	{
		m_Data[intRow * m_SizeCols + intCol] = dblVal;
	}
	void setValAtIdx(const int64_t &intIdx, const double &dblVal)
	{
		m_Data[intIdx] = dblVal;
	}

	// fill matrix with values from mxArray
	void fillVals(const mxArray *ptrMatlabVector)
	{
		// define variables
		double *ptrElement;	// pointer to elements in matlab array
		mwSize intNumElements, intIdx; // integer for number of elements and index counter
		intNumElements = mxGetNumberOfElements(ptrMatlabVector); // get number of elements in vector
		
		 // get starting memory location
		ptrElement = mxGetPr(ptrMatlabVector);

		// assign all values to output
		for (intIdx = 0; intIdx < intNumElements; intIdx++)
		{
			m_Data[intIdx] = *ptrElement++;
		}
	}

	// fill mxArray with values
	mxArray* exportToMx()
	{
		// define variables
		mxArray *mxMat;
		int64_t intNrElements = getSize();
		int64_t intRows = getSizeRows();
		int64_t intCols = getSizeCols();
		double  *ptrElement;
		mwSize intIdx;
		
		// assign content
		mxMat = mxCreateDoubleMatrix((mwSize)intRows, (mwSize)intCols, mxREAL); // Create an m-by-n mxArray; will copy existing data into it
		ptrElement = mxGetPr(mxMat);

		// loop through all elements
		for (intIdx = 0; intIdx < intNrElements; intIdx++)
		{
			ptrElement[intIdx] = m_Data[intIdx]; // Copy data into the mxArray
		}

		// return filled matlab matrix
		return mxMat;
	}
};

// define cell array alias
typedef int64_t* ptrMem_t; // pointer is of size [int64_t] because we're using an x64 system
typedef vector<vector<Matrix>> cellArrayStim; // vector with pointers to constituent vectors
typedef vector<vector<double>> cellArraySpikes; // vector with pointers to constituent vectors

//define the fields of input/output data structure
struct sAggregate // aggregate data structure
{
	vector<double>	vecOverallT;
	double			dblDeltaT;
	Matrix			matCortConn;					// size: [intCortexSynapses x 2] (from <cell1> - to <cell2>)
	double			dblSynSpikeMem;
	vector<int64_t>	vecCortSynType;					// size: [intCortexSynapses x 1]
	int				intCortexCells;
	vector<double>	vecCortDelay;					// size: [intCortexSynapses x 1]
	vector<double>	vecCortConductance;				// size: [intCortexSynapses x 1]
	vector<double>	vecCellThresh;					// size: [intCortexCells x 1]
	vector<double>	vecTauPeakByType;				// size: [intCellTypes x 1]
	vector<double>	vecCellV_E;						// size: [intCortexCells x 1]
	vector<double>	vecCellV_I;						// size: [intCortexCells x 1]
	vector<double>	vecCellV_AHP;					// size: [intCortexCells x 1]
	vector<double>	vecCellV_Leak;					// size: [intCortexCells x 1]
	vector<double>	vecCellCm; 						// size: [intCortexCells x 1]
	vector<double>	vecCellG_Leak;					// size: [intCortexCells x 1]
	vector<double>	vecCellG_AHP;					// size: [intCortexCells x 1]
	vector<double>	vecSynConductanceON_to_Cort;	// size: [intLGNSynapses x 1]
	vector<double>	vecSynConductanceOFF_to_Cort;	// size: [intLGNSynapses x 1]
	vector<double>	vecSynWeightON_to_Cort;			// size: [intLGNSynapses x 1]
	vector<double>	vecSynWeightOFF_to_Cort;		// size: [intLGNSynapses x 1]
	vector<double>	vecSynDelayON_to_Cort;			// size: [intLGNSynapses x 1]
	vector<double>	vecSynDelayOFF_to_Cort;			// size: [intLGNSynapses x 1]
	Matrix			matSynConnON_to_Cort;			// size: [intLGNSynapses x 2] (from <cell1> - to <cell2>)
	Matrix			matSynConnOFF_to_Cort;			// size: [intLGNSynapses x 2] (from <cell1> - to <cell2>)
	Matrix			matBlankLGN_ON;					// size: [intGridSizeLGN x intGridSizeLGN] (21 x 21)
	Matrix			matBlankLGN_OFF;				// size: [intGridSizeLGN x intGridSizeLGN] (21 x 21)
	cellArrayStim	cellLGN_ON;
	cellArrayStim	cellLGN_OFF;
	vector<double>	vecTrialOris;
	vector<int64_t>	vecTrialOriIdx;
	vector<double>	vecStimStartSecs;
	vector<double>	vecTrialEndSecs;
	vector<double>	vecThisV;
	int				boolStimPresent;
	int				intPrevTrial;
	int				intTrialT;
	int64_t			intIter;
	cellArraySpikes	cellSpikeTimesLGN_ON;		// use vector.reserve() to check reserved versus used elements
	cellArraySpikes	cellSpikeTimesLGN_OFF;
	cellArraySpikes	cellSpikeTimesCortex;
	vector<int64_t> vecSpikeCounterLGN_ON;		// keeps track of the number of assigned spikes across vectors cellSpikeTimesLGN_ON >> compare with vector.size() to check when vector.resize() is needed
	vector<int64_t> vecSpikeCounterLGN_OFF;	// keeps track of the number of assigned spikes across vectors cellSpikeTimesLGN_OFF 
	vector<int64_t> vecSpikeCounterCortex;		// keeps track of the number of assigned spikes across vectors cellSpikeTimesCortex 
	int64_t			intPreAllocationSize;		// size of pre-allocation blocks
};

// set field names for MATLAB structure
const char *field_names_sAggregate[] =
{
	"vecOverallT",
	"dblDeltaT",
	"matCortConn",
	"dblSynSpikeMem",
	"vecCortSynType",
	"intCortexCells",
	"vecCortDelay",
	"vecCortConductance",
	"vecCellThresh",
	"vecTauPeakByType",
	"vecCellV_E",
	"vecCellV_I",
	"vecCellV_AHP",
	"vecCellV_Leak",
	"vecCellCm",
	"vecCellG_Leak",
	"vecCellG_AHP",
	"vecSynConductanceON_to_Cort",
	"vecSynConductanceOFF_to_Cort",
	"vecSynWeightON_to_Cort",
	"vecSynWeightOFF_to_Cort",
	"vecSynDelayON_to_Cort",
	"vecSynDelayOFF_to_Cort",
	"matSynConnON_to_Cort",
	"matSynConnOFF_to_Cort",
	"matBlankLGN_ON",
	"matBlankLGN_OFF",
	"cellLGN_ON",
	"cellLGN_OFF",
	"vecTrialOris",
	"vecTrialOriIdx",
	"vecStimStartSecs",
	"vecTrialEndSecs",
	"vecThisV",
	"boolStimPresent",
	"intPrevTrial",
	"intTrialT",
	"intIter",
	"cellSpikeTimesLGN_ON", // use vector.reserve() or vector.size() to check allocated versus used elements
	"cellSpikeTimesLGN_OFF",
	"cellSpikeTimesCortex",
	"vecSpikeCounterLGN_ON", // keeps track of the number of assigned spikes across vectors cellSpikeTimesLGN_ON >> compare with vector.size() to check when vector.resize() is needed
	"vecSpikeCounterLGN_OFF", // keeps track of the number of assigned spikes across vectors cellSpikeTimesLGN_OFF 
	"vecSpikeCounterCortex", // keeps track of the number of assigned spikes across vectors cellSpikeTimesCortex 
	"intPreAllocationSize" // size of pre-allocation blocks
};