// runSimEXE.cpp : Defines the entry point for the console application.
//
/*
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime> // for time()
#include <random> // for std::random_device and std::mt19937
#include <stdint.h>
#include <array>
#include <vector>
#include <stdio.h>
#include <windows.h>
*/

// define headers ==> all other includes are in stdafx.h
//#pragma once
#include "stdafx.h"
#include "dataDefinitions.h"
#include <chrono>
#pragma warning(disable : 4996)

// declare functions
//int runSimLoop(sAggregate);
//int loadData(sAggregate, const mxArray*);
//int saveData();


//// FUNCTIONS
// PSP function
double fPSP(const double &dblT, const vector<double>&vecSpikeTimes, const double &dblTauPeak)
{
	//fPSP = @(dblT,vecSpikeTimes,dblTauPeak) sum(max(0,dblT - vecSpikeTimes).*(exp(1)/dblTauPeak).*exp(-(dblT - vecSpikeTimes)/dblTauPeak));

	// vecExp = exp(-(vecTempOffset) / dblTauPeak)
	const vector<double> vecTempOffset = vMult(vMinus(dblT, vecSpikeTimes), (double)-1);

	// return values; sum across: get values of [get element-wise max[(dblT - vecSpikeTimes),0]] times [vecExp]
	return (exp((double)1) / dblTauPeak) * vSum(vMult(vMax(vecTempOffset, (double)0), vExp(vDiv(vecTempOffset, dblTauPeak))));
}
vector<double> getSumPSPsPerCell(Matrix &matC, const vector<double> &vecPSPs, const int64_t &intCortexCells)
{
	// allocate output vector
	vector<double> vecSummedPSPs(intCortexCells);
	std::fill(vecSummedPSPs.begin(), vecSummedPSPs.end(), 0);
	
	// sum all PSPs
	for (int64_t intS = 0; intS < vecPSPs.size(); intS++)
	{
		vecSummedPSPs[(int64_t)matC.getVal(intS, 2) - 1] += vecPSPs[intS];
	}
	return vecSummedPSPs;
}

// generate random numbers
void vRand(vector<double> &vecRand)
{
	// generate random numbers
	random_device rd;
	mt19937_64 mersenne_engine(rd());
	uniform_real_distribution<double> dist(0, 1);
	auto gen = bind(dist, mersenne_engine);
	generate(begin(vecRand), end(vecRand), gen);
}

// run actual simulation loop
int runSimLoop(sAggregate &sData)
{
	//// run simulation
	// define variables
	const auto sResult = minmax_element(sData.vecOverallT.begin(), sData.vecOverallT.end());
	const double dblStartT = *sResult.first;
	const double dblStopT = *sResult.second;
	int64_t intThisTrial;
	double	dblOri;
	int64_t	intOri;
	int64_t intStopTrial;
	const int64_t intStimFrames = sData.cellLGN_ON[0].size();
	const vector<double> vecBlankLGN_ON = sData.matBlankLGN_ON.getAllVals();
	const vector<double> vecBlankLGN_OFF = sData.matBlankLGN_OFF.getAllVals();
	vector<double> vecRotLGN_ON (vecBlankLGN_ON.size());
	vector<double> vecRotLGN_OFF(vecBlankLGN_ON.size());
	vector<double> vecThisIterON(vecBlankLGN_ON.size());
	vector<double> vecThisIterOFF(vecBlankLGN_ON.size());
	time_t now = time(0);
	int intTempSubType;
	int64_t intNumSynapsesLGN = sData.vecSynWeightOFF_to_Cort.size();
	int64_t intNumSynapsesCortex = sData.vecCortSynType.size();

	/*
	To create a valarray from a vector:
	valarray<double> vav(vec.data(), vec.size());

	To write the data back into a vector:
	vec.assign(begin(vav), end(vav));
	*/

	// define aliases for sData
	auto &intPrevTrial = sData.intPrevTrial;
	auto &vecThisV = sData.vecThisV;
	auto &boolStimPresent = sData.boolStimPresent;
	auto &intTrialT = sData.intTrialT;
	auto &intIter = sData.intIter;
	auto &cellLGN_ON = sData.cellLGN_ON;
	auto &cellLGN_OFF = sData.cellLGN_OFF;
	auto &vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	auto &vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	auto &vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
	auto &cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	auto &cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	auto &cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	auto &intPreAllocationSize = sData.intPreAllocationSize;
	auto &dblSynSpikeMem = sData.dblSynSpikeMem;
	auto &vecTauPeakByType = sData.vecTauPeakByType;
	auto &intCortexCells = sData.intCortexCells;
	auto &matCortConn = sData.matCortConn;
	auto &dblDeltaT = sData.dblDeltaT;

	// create temporary storage variables
	vector<double> vecRandVals (vecBlankLGN_ON.size());
	vector<double> vecSW;
	vector<double> vecSC;
	vector<double> vecSD;
	Matrix matC;
	cellArraySpikes cellSpikes;
	vector<double> vecTempLGNCurrents (sData.intCortexCells);
	vector<double> vecTempExcCurrents(sData.intCortexCells);
	vector<double> vecTempInhCurrents(sData.intCortexCells);
	vector<double> vecTempAHPCurrents(sData.intCortexCells);
	vector<double> vecEPSPs(intNumSynapsesCortex);
	vector<double> vecIPSPs(intNumSynapsesCortex);
	vector<double> vecPrevV = vecThisV;
	vector<double> vecTempSpiking(sData.intCortexCells);
	vector<double> vecAllSTs;
	vector<double> vecSelect;
	vector<double> vecSpikeTimes;

	// run
	for (const double &dblCurT : sData.vecOverallT)
	{
		// check which trial
		intIter++;
		intThisTrial = vSum(vLT(dblCurT, sData.vecStimStartSecs));

		// prepare trial
		if (intThisTrial > intPrevTrial) //start new trial
		{
			//// get trial data
			intPrevTrial = intThisTrial;
			dblOri = sData.vecTrialOris[intThisTrial];
			intOri = sData.vecTrialOriIdx[intThisTrial];

			//set counter for this trial
			intTrialT = 0;
			boolStimPresent = true;
		}

		//check to stop trial
		intStopTrial = vSum(vLT(dblCurT,sData.vecTrialEndSecs));
		if (intStopTrial >= intThisTrial)
		{
			boolStimPresent = false;
		}
		
		
		// check if stimulus should be presented
		intTrialT += intTrialT;
		if ((boolStimPresent==1) & (intTrialT < intStimFrames))
		{
			//rotate template LGN response
			vecRotLGN_ON = cellLGN_ON[intOri][intTrialT].getAllVals();
			vecRotLGN_OFF = cellLGN_OFF[intOri][intTrialT].getAllVals();

			//create stochastic spiking
			vRand(vecRandVals); // generate random numbers
			vecThisIterON = vLT(vecRotLGN_ON, vecRandVals);
			vRand(vecRandVals); // generate more random numbers
			vecThisIterOFF = vLT(vecRotLGN_OFF, vecRandVals);


			//mexPrintf("Stim! size cellLGN_ON[intOri][intTrialT]: %d; vecRotLGN_ON: %d; vecRand: %d\n", cellLGN_ON[intOri][intTrialT].getSize(), vecRotLGN_ON.size(), vecRandVals.size());
			//mexEvalString("pause(eps);");
		}
		else
		{
			//create stochastic spiking
			vRand(vecRandVals); // generate random numbers
			vecThisIterON = vLT(vecBlankLGN_ON, vecRandVals);
			vRand(vecRandVals); // generate more random numbers
			vecThisIterOFF = vLT(vecBlankLGN_OFF, vecRandVals);

			const vector<double> vecIn1 = { 1,2,3 };
			const vector<double> vecIn2 = { 3,2,1 };
			vector<double> vecOut = vLT(vecIn1, vecIn2);
			//mexPrintf("  <runSimMEX> vLT vals: [%f %f %f] (size: %d)\n", vecOut[0], vecOut[1], vecOut[2], vecOut.size());

			//mexPrintf("No stim! values vecBlankLGN_ON: [%f %f %f] (size: %d)\n", vecBlankLGN_ON[0], vecBlankLGN_ON[1], vecBlankLGN_ON[2], vecBlankLGN_ON.size());
			//mexPrintf("No stim! values vecRand:        [%f %f %f] (size: %d)\n", vecRandVals[0], vecRandVals[1], vecRandVals[2], vecRandVals.size());
			//mexPrintf("No stim! values vecThisIterON:  [%f %f %f] (size: %d)\n", vecThisIterON[0], vecThisIterON[1], vecThisIterON[2], vecThisIterON.size());

			//mexPrintf("No stim! size vecBlankLGN_ON: %d; vecRand: %d\n", vecBlankLGN_ON.size(), vecRandVals.size());
			//mexEvalString("pause(eps);");
		}
		//mexPrintf("vecThisIterON: %d; vecThisIterOFF: %d\n", vecThisIterON.size(), vecThisIterOFF.size());
		//mexEvalString("pause(eps);");

		// assign new LGN spikes to cell arrays; skip if no spikes
		if (any_of(vecThisIterON.begin(), vecThisIterON.end(), [](int i) {return i == 1; }))
		{
			// loop through all elements
			for (int64_t intN = 0; intN < vecThisIterON.size(); intN++)
			{
				// check if this cell was active, skip otherwise
				if (vecThisIterON[intN] == 1)
				{
					//check number of spikes and pre - allocated size
					if (cellSpikeTimesLGN_ON[intN].capacity() <= cellSpikeTimesLGN_ON[intN].size());
					{
						// increase vector size
						cellSpikeTimesLGN_ON[intN].reserve(cellSpikeTimesLGN_ON[intN].capacity() + intPreAllocationSize);
					}

					// add spike time
					cellSpikeTimesLGN_ON[intN].push_back(dblCurT);
				}
			}
		}
		
		if (any_of(vecThisIterOFF.begin(), vecThisIterOFF.end(), [](int i) {return i == 1; }))
		{
			// loop through all elements
			for (int64_t intN = 0; intN < vecThisIterOFF.size(); intN++)
			{
				// check if this cell was active, skip otherwise
				if (vecThisIterOFF[intN] == 1)
				{
					//check number of spikes and pre - allocated size
					if (cellSpikeTimesLGN_OFF[intN].capacity() <= cellSpikeTimesLGN_OFF[intN].size());
					{
						// increase vector size
						cellSpikeTimesLGN_OFF[intN].reserve(cellSpikeTimesLGN_OFF[intN].capacity() + intPreAllocationSize);
					}

					// add spike time
					cellSpikeTimesLGN_OFF[intN].push_back(dblCurT);
				}
			}
		}
		
		// calculate all excitatory currents due to spikes from LGN
		intTempSubType = 1;
		std::fill(vecTempLGNCurrents.begin(), vecTempLGNCurrents.end(), 0); // clear temp vector
		for (int intField = 1; intField < 3; intField++)
		{
			if (intField == 1)
			{
				vecSW = sData.vecSynWeightON_to_Cort;
				vecSC = sData.vecSynConductanceON_to_Cort;
				vecSD = sData.vecSynDelayON_to_Cort;
				matC = sData.matSynConnON_to_Cort;
				cellSpikes = cellSpikeTimesLGN_ON;
			}
			else
			{
				vecSW = sData.vecSynWeightOFF_to_Cort;
				vecSC = sData.vecSynConductanceOFF_to_Cort;
				vecSD = sData.vecSynDelayOFF_to_Cort;
				matC = sData.matSynConnOFF_to_Cort;
				cellSpikes = cellSpikeTimesLGN_OFF;
			}
			
			// create temp variables
			vector<double> vecPSPs (intNumSynapsesLGN);
			for (int64_t intSynapse = 0; intSynapse < intNumSynapsesLGN; intSynapse++)
			{
				// remove all spikes older than[dblSynSpikeMem] seconds to avoid
				//infinite memory requirements as simulation gets longer
				int64_t intNeuronLGN = matC.getVal(intSynapse, 1) - 1;
				vecAllSTs = cellSpikes[intNeuronLGN]; // get spike times for this synapse's cell
				vecSelect = vST(vMinus(vecAllSTs, dblCurT), dblSynSpikeMem);
				vecSpikeTimes = vSelect(vecAllSTs, vecSelect); // create index for spike times newer than dblCurT - dblSynSpikeMem and select subset
				if (vecSpikeTimes.size() > 0) // skip if no spikes
				{
					//get synapse data
					const double &dblWeight = vecSW[intSynapse];
					const double &dblDelay = vecSD[intSynapse];
					const double &dblTauPeak = vecTauPeakByType[intTempSubType];
					const double &dblCond = vecSC[intSynapse];
					const double dblT = dblCurT - dblDelay;

					//calculate PSPs
					vecPSPs[intSynapse] = dblWeight * dblCond * fPSP(dblT, vecSpikeTimes, dblTauPeak);
				}
			}
			// sum all transmissions incoming to each neuron
			vecTempLGNCurrents = vPlus(vecTempLGNCurrents, getSumPSPsPerCell(matC, vecPSPs, intCortexCells));
		}
		
		
		//// calculate currents from intra - cortical connections
		// create temp variables
		std::fill(vecTempExcCurrents.begin(), vecTempExcCurrents.end(), 0); // clear temp vector
		std::fill(vecTempInhCurrents.begin(), vecTempInhCurrents.end(), 0); // clear temp vector
		std::fill(vecEPSPs.begin(), vecEPSPs.end(), 0); // clear temp vector
		std::fill(vecIPSPs.begin(), vecIPSPs.end(), 0); // clear temp vector
		
		// run through synapses
		for (int64_t intSynapse = 0; intSynapse < intNumSynapsesCortex; intSynapse++)
		{
			// remove all spikes older than[dblSynSpikeMem] seconds to avoid
			//infinite memory requirements as simulation gets longer
			int64_t intNeuron = matCortConn.getVal(intSynapse, 1) - 1;
			int64_t intNeuronTarget = matCortConn.getVal(intSynapse, 2) - 1;
			//mexPrintf("Synapse [%d/%d], originating from neuron %d, projecting to neuron %d \n", intSynapse, intNumSynapsesCortex, intNeuron, intNeuronTarget);
			//mexEvalString("pause(eps);");
			
			vecAllSTs = cellSpikeTimesCortex[intNeuron]; // get spike times for this synapse's cell
			vecSelect = vST(vMinus(vecAllSTs, dblCurT), dblSynSpikeMem); // create index for spike times newer than dblCurT - dblSynSpikeMem
			vecSpikeTimes = vSelect(vecAllSTs, vecSelect); // select subset
			//mexPrintf("size vecAllSTS: %d; vecSelect: %d; vecSpikeTimes: %d\n", vecAllSTs.size(), vecSelect.size(), vecSpikeTimes.size());
			//mexEvalString("pause(eps);");
			if (vecSpikeTimes.size() > 0) // skip if no spikes
			{
				//mexPrintf("vecSpikeTime vals: [%f %f %f]\n", vecSpikeTimes[0], vecSpikeTimes[1], vecSpikeTimes[2]);
				//mexEvalString("pause(eps);");

				//get synapse data
				intTempSubType = sData.vecCortSynType[intSynapse];
				const double dblDelay = sData.vecCortDelay[intSynapse];
				const double dblTauPeak = vecTauPeakByType[intTempSubType];
				const double dblCond = sData.vecCortConductance[intSynapse];
				const double dblT = dblCurT - dblDelay;

				//calculate PSPs
				if (intTempSubType == 1)
				{
					vecEPSPs[intSynapse] = dblCond * fPSP(dblT, vecSpikeTimes, dblTauPeak);
				}
				else
				{
					vecIPSPs[intSynapse] = dblCond * fPSP(dblT, vecSpikeTimes, dblTauPeak);
				}
			}
		}
		
		
		// sum all transmissions incoming to each neuron
		vecTempExcCurrents = vPlus(vecTempExcCurrents, getSumPSPsPerCell(matCortConn, vecEPSPs, intCortexCells));
		vecTempInhCurrents = vPlus(vecTempInhCurrents, getSumPSPsPerCell(matCortConn, vecIPSPs, intCortexCells));


		//// calculate after hyper polarizations
		intTempSubType = 2;
		std::fill(vecTempAHPCurrents.begin(), vecTempAHPCurrents.end(), 0); // clear temp vector
		
		for (int64_t intNeuron = 0 ; intNeuron < intCortexCells; intNeuron++)
		{
			// remove all spikes older than[dblSynSpikeMem] seconds to avoid
			//infinite memory requirements as simulation gets longer
			const vector<double> &vecAllSTs = cellSpikeTimesCortex[intNeuron]; // get spike times for this cell
			const vector<double> vecSelectSTs = vST(vMinus(vecAllSTs, dblCurT), dblSynSpikeMem); // create index for spike times newer than dblCurT - dblSynSpikeMem
			const vector<double> vecSpikeTimes = vSelect(vecAllSTs, vecSelectSTs); // select subset
			if (vecSpikeTimes.size() > 0) // skip if no spikes
			{
				//get synapse data
				constexpr double dblDelay = 0;
				const double dblTauPeak = vecTauPeakByType[intTempSubType];
				const double dblCond = sData.vecCellG_AHP[intNeuron];
				const double dblT = dblCurT - dblDelay;

				//calculate AHPs
				vecTempAHPCurrents[intNeuron] = dblCond * fPSP(dblT, vecSpikeTimes, dblTauPeak);
			}
		}

		//// integrate all inputs
		vecPrevV = vecThisV;
		vecThisV = vecPrevV;
		vecThisV = vPlus(vecThisV,vMult(dblDeltaT,
			vMinus(vMinus(vMinus(vMinus(
				vMult(-1, vMult(vecTempLGNCurrents, (vMinus(vecPrevV, sData.vecCellV_E))))
				, vMult(vecTempExcCurrents, (vMinus(vecPrevV, sData.vecCellV_E))))
				, vMult(vecTempInhCurrents, (vMinus(vecPrevV, sData.vecCellV_I))))
				, vMult(sData.vecCellG_Leak, (vMinus(vecPrevV, sData.vecCellV_Leak))))
				, vMult(vecTempAHPCurrents, (vMinus(vecPrevV, sData.vecCellV_AHP))))
			));
		vecThisV = vDiv(vecThisV, sData.vecCellCm);
		
		//// spiking
		vecTempSpiking = vLT(vecThisV,sData.vecCellThresh);
		if (any_of(vecTempSpiking.begin(), vecTempSpiking.end(), [](int i) {return i == 1; }))
		{
			// loop through all elements
			for (int64_t intN = 0; intN < vecTempSpiking.size(); intN++)
			{
				if (vecTempSpiking[intN] == 1)
				{
					//check number of spikes and pre - allocated size
					vecSpikeCounterCortex[intN]++;
					if (cellSpikeTimesCortex[intN].size() < vecSpikeCounterCortex[intN]);
					{
						// increase vector size
						cellSpikeTimesCortex[intN].resize(cellSpikeTimesCortex[intN].size() + intPreAllocationSize);
					}

					// add spike time
					cellSpikeTimesCortex[intN][vecSpikeCounterCortex[intN]] = dblCurT;
				}
			}
		}
		
		
		////message
		if (intIter%10 == 0) { now = time(0); mexPrintf("\t  <runSimMEX> Now at t=%.3fs / %.3fs [%s\b]\n", dblCurT, dblStopT, ctime(&now)); mexEvalString("pause(eps);"); }
	}

	// return 0-flag
	return 0;
}

//read vector from matlab matrix
vector<double> getDblVectorFromMx(const mxArray *ptrMatlabVector)
{
	// define variables
	double *ptrElement;	// pointer to elements in matlab vector
	mwSize intNumElements, intIdx; // integer for number of elements and index counter
	intNumElements = mxGetNumberOfElements(ptrMatlabVector); // get number of elements in vector
	vector<double> vecOut; // create output vector
	vecOut.resize(intNumElements); // resize to number of elements in matlab vector

	// get starting memory location
	ptrElement = mxGetPr(ptrMatlabVector);
	
	// assign all values to output
	for (intIdx = 0; intIdx<intNumElements; intIdx++)
	{
		vecOut[intIdx] = *ptrElement++;
	}
	return vecOut;
}

// assign vector to matlab vector matrix
mxArray* putDblVectorToMx(const vector<double> &vecIn)
{
	// define variables
	mxArray *mxVec; 
	int intNrElements = vecIn.size();
	mxVec = mxCreateDoubleMatrix((mwSize)intNrElements, 1, mxREAL); // Create an m-by-n mxArray; will copy existing data into it
	double  *ptrElement;
	mwSize intIdx;
	ptrElement = mxGetPr(mxVec);

	// loop through all elements
	for (intIdx = 0; intIdx < intNrElements; intIdx++)
	{
		ptrElement[intIdx] = vecIn[intIdx]; // Copy data into the mxArray
	}

	// return filled matlab vector
	return mxVec; 
}

//read int vector from matlab matrix
vector<int64_t> getIntVectorFromMx(const mxArray *ptrMatlabVector)
{
	// define variables
	double *ptrElement;	// pointer to elements in matlab vector
	mwSize intNumElements, intIdx; // integer for number of elements and index counter
	intNumElements = mxGetNumberOfElements(ptrMatlabVector); // get number of elements in vector
	vector<int64_t> vecOut; // create output vector
	vecOut.resize(intNumElements); // resize to number of elements in matlab vector

	// get starting memory location
	ptrElement = mxGetPr(ptrMatlabVector);

	// assign all values to output
	for (intIdx = 0; intIdx<intNumElements; intIdx++)
	{
		vecOut[intIdx] = *ptrElement++;
	}
	return vecOut;
}

// assign int vector to matlab vector matrix
mxArray* putIntVectorToMx(const vector<int64_t> &vecIn)
{
	// define variables
	mxArray *mxVec;
	int64_t intNrElements = vecIn.size();
	mxVec = mxCreateDoubleMatrix((mwSize)intNrElements, 1, mxREAL); // Create an m-by-n mxArray; will copy existing data into it
	double  *ptrElement;
	mwSize intIdx;
	ptrElement = mxGetPr(mxVec);

	// loop through all elements
	for (intIdx = 0; intIdx < intNrElements; intIdx++)
	{
		ptrElement[intIdx] = vecIn[intIdx]; // Copy data into the mxArray
	}

	// return filled matlab vector
	return mxVec;
}

// assign int vector to matlab vector matrix
cellArraySpikes putMxCellToVecVec(const mxArray *ptrArray)
{
	// define variables
	mwIndex cellIdx = 0;			// index of cell array access
	const mxArray *ptrCell;			// pointer to a single cell
	cellArraySpikes cellOut;
	int64_t intNumEl = mxGetNumberOfElements(ptrArray);
	cellOut.resize(intNumEl);

	// loop through cell array and assign values to vector
	for (cellIdx = 0; cellIdx<intNumEl; cellIdx++)
	{
		ptrCell = mxGetCell(ptrArray, cellIdx);
		//intNumCellSize = mxGetNumberOfElements(ptrCell);
		cellOut[cellIdx] = getDblVectorFromMx(ptrCell);
	}

	// return
	return cellOut;
}

// assign int vector to matlab vector matrix
mxArray* putVecVecToMxCell(const cellArraySpikes &cellIn)
{
	// define variables
	mwIndex cellIdx = 0;			// index of cell array access
	mxArray *mxCellArray;			// pointer to mx cell array
	mxArray *mxVec;					// pointer to mx vector
	int64_t intNumEl = cellIn.size();
	
	// create cell array
	mxCellArray = mxCreateCellMatrix((mwSize)intNumEl, 1);

	// loop through cell array and assign values to vector
	for (cellIdx = 0; cellIdx<intNumEl; cellIdx++)
	{
		//transform vector contents to mxArray
		mxVec = putDblVectorToMx(cellIn[cellIdx]);

		//put vector in cell;
		mxSetCell(mxCellArray, cellIdx, mxVec);
	}

	// return
	return mxCellArray;
}

// assign matlab cell array to vector<vector<Matrix>>
cellArrayStim putMxCellStimToVecMat(const mxArray *ptrArray)
{
	// define variables
	mwIndex cellIdx = 0;			// index of cell array access
	const mxArray *ptrCell;			// pointer to a single cell
	cellArrayStim cellOut;
	int64_t intNumStims = mxGetNumberOfElements(ptrArray);
	cellOut.resize(intNumStims);

	// define cell variables
	const mwSize *arrDimSize;
	int intFrameSize;
	int intTimePoints;
	int intT;
	Matrix MatTemp;
	double *ptrElement;	// pointer to elements in matlab array
	mwSize intIdx; // for indexing within single 2D frame

	// loop through cell array and assign values to vector
	for (cellIdx = 0; cellIdx<intNumStims; cellIdx++)
	{
		// points to 3D matrix; [X x Y x T]
		ptrCell = mxGetCell(ptrArray, cellIdx);
		arrDimSize = mxGetDimensions(ptrCell); // [X x Y x T]

		// calculate size of [X x Y] matrix
		intFrameSize = arrDimSize[0] * arrDimSize[1];
		intTimePoints = arrDimSize[2];
		cellOut[cellIdx].resize(intTimePoints);

		// set pointer to cell array start
		ptrElement = mxGetPr(ptrCell);

		// loop through time points
		for (intT = 0; intT < intTimePoints; intT++)
		{
			// create matrix size
			if (MatTemp.getSize() != intFrameSize)
			{
				MatTemp.setSizeCols(arrDimSize[1]);
				MatTemp.setSizeRows(arrDimSize[0]);
			}

			// fill matrix with values
			for (intIdx = 0; intIdx < intFrameSize; intIdx++)
			{
				MatTemp.setValAtIdx(intIdx, *ptrElement++);
			}

			// assign matrix object to C-structure field
			cellOut[cellIdx][intT] = MatTemp;
		}
	}

	// return
	return cellOut;
}

// assign vector<vector<Matrix>> to matlab cell array
mxArray* putVecMatToMxCellStim(cellArrayStim &cellIn)
{
	// define variables
	mwIndex cellIdx = 0;			// index of cell array access
	mxArray *mxCellArray;			// pointer to mx cell array
	mxArray *mxMat;					// pointer to mx matrix
	int64_t intNumStims = cellIn.size();
	
	// create cell array
	mxCellArray = mxCreateCellMatrix((mwSize)intNumStims, 1);

	// loop through cell array and assign values to vector
	// define cell variables
	const mwSize arrDimSize[3]{ (mwSize)cellIn[0][0].getSizeRows(), (mwSize)cellIn[0][0].getSizeCols(), (mwSize)cellIn[0].size() };
	int64_t intFrameSize;
	int64_t intTimePoints;
	int64_t intT;
	Matrix MatTemp;
	double *ptrElement;	// pointer to elements in matlab array
	int64_t intFrameEl; // for indexing within single 2D frame
	mwSize intIdx; // for indexing across entire 3D matrix
	vector<int64_t> vecSize;

	// loop through cell array and assign values to vector
	for (cellIdx = 0; cellIdx<intNumStims; cellIdx++)
	{
		// points to 3D matrix; [X x Y x T]
		intTimePoints = cellIn[cellIdx].size();
		MatTemp = cellIn[cellIdx][0];
		vecSize = { MatTemp.getSizeRows(), MatTemp.getSizeCols(), intTimePoints};
		intFrameSize = MatTemp.getSizeRows() * MatTemp.getSizeCols();// calculate size of [X x Y] matrix		   
		
		// create 3D matrix
		mxMat = mxCreateNumericArray((mwSize)3, arrDimSize, mxDOUBLE_CLASS, mxREAL);
		ptrElement = mxGetPr(mxMat);
		intIdx = 0;

		// loop through time points
		for (intT = 0; intT < intTimePoints; intT++)
		{
			// get data
			MatTemp = cellIn[cellIdx][intT];

			// loop through all elements
			for (intFrameEl = 0; intFrameEl < intFrameSize; intFrameEl++)
			{
				ptrElement[intIdx++] = MatTemp.getValAtIdx(intFrameEl); // Copy data into the mxArray
			}
		}

		//put matrix in cell;
		mxSetCell(mxCellArray, cellIdx, mxMat);
	}

	// return
	return mxCellArray;
}

// transform matlab structure into c structure
int loadData(sAggregate &sData, const mxArray *sMxData)
{
	//read field values from matlab structure and assign to c-structure
	// define re-useable variables
	const mxArray *fPtr;			// field pointer
	mwIndex mIdx = 0;				// index of structure array access
	double *dblPtr;					// pointer to double data
	const mwIndex *arrMatSize;
	
	// read dblDeltaT
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "dblDeltaT"));
	sData.dblDeltaT = dblPtr[0];

	// read vecOverallT
	sData.vecOverallT = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecOverallT"));

	// read matCortConn
	fPtr = mxGetField(sMxData, mIdx, "matCortConn");
	arrMatSize = mxGetDimensions(fPtr);
	sData.matCortConn.setSizeCols(arrMatSize[1]);
	sData.matCortConn.setSizeRows(arrMatSize[0]);
	sData.matCortConn.fillVals(fPtr);

	// read dblSynSpikeMem;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "dblSynSpikeMem"));
	sData.dblSynSpikeMem = dblPtr[0];
	
	// read vecCortSynType;
	sData.vecCortSynType = getIntVectorFromMx(mxGetField(sMxData, mIdx, "vecCortSynType"));

	// read intCortexCells;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "intCortexCells"));
	sData.intCortexCells = (int)dblPtr[0];
	
	// read vecCortDelay;
	sData.vecCortDelay = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCortDelay"));
	
	// read vecCortConductance;
	sData.vecCortConductance = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCortConductance"));
	
	// read vecCellThresh;
	sData.vecCellThresh = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellThresh"));
	
	// read vecTauPeakByType;
	sData.vecTauPeakByType = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecTauPeakByType"));
	
	// read vecCellV_E;
	sData.vecCellV_E = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellV_E"));
	
	// read vecCellV_I;
	sData.vecCellV_I = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellV_I"));
	
	// read vecCellV_AHP;
	sData.vecCellV_AHP = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellV_AHP"));
	
	// read vecCellV_Leak;
	sData.vecCellV_Leak = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellV_Leak"));
	
	// read vecCellCm;
	sData.vecCellCm = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellCm"));
	
	// read vecCellG_Leak;
	sData.vecCellG_Leak = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellG_Leak"));
	
	// read vecCellG_AHP;
	sData.vecCellG_AHP = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecCellG_AHP"));
	
	// read vecSynConductanceON_to_Cort;
	sData.vecSynConductanceON_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynConductanceON_to_Cort"));
	
	// read vecSynConductanceOFF_to_Cort;
	sData.vecSynConductanceOFF_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynConductanceOFF_to_Cort"));
	
	// read vecSynWeightON_to_Cort;
	sData.vecSynWeightON_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynWeightON_to_Cort"));
	
	// read vecSynWeightOFF_to_Cort;
	sData.vecSynWeightOFF_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynWeightOFF_to_Cort"));
	
	// read vecSynDelayON_to_Cort;
	sData.vecSynDelayON_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynDelayON_to_Cort"));
	
	// read vecSynDelayOFF_to_Cort;
	sData.vecSynDelayOFF_to_Cort = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecSynDelayOFF_to_Cort"));
	
	// read matSynConnON_to_Cort;
	fPtr = mxGetField(sMxData, mIdx, "matSynConnON_to_Cort");
	arrMatSize = mxGetDimensions(fPtr);
	sData.matSynConnON_to_Cort.setSizeCols(arrMatSize[1]);
	sData.matSynConnON_to_Cort.setSizeRows(arrMatSize[0]);
	sData.matSynConnON_to_Cort.fillVals(fPtr);

	// read matSynConnOFF_to_Cort;
	fPtr = mxGetField(sMxData, mIdx, "matSynConnOFF_to_Cort");
	arrMatSize = mxGetDimensions(fPtr);
	sData.matSynConnOFF_to_Cort.setSizeCols(arrMatSize[1]);
	sData.matSynConnOFF_to_Cort.setSizeRows(arrMatSize[0]);
	sData.matSynConnOFF_to_Cort.fillVals(fPtr);

	// read matBlankLGN_ON;
	fPtr = mxGetField(sMxData, mIdx, "matBlankLGN_ON");
	arrMatSize = mxGetDimensions(fPtr);
	sData.matBlankLGN_ON.setSizeCols(arrMatSize[1]);
	sData.matBlankLGN_ON.setSizeRows(arrMatSize[0]);
	sData.matBlankLGN_ON.fillVals(fPtr);

	// read matBlankLGN_OFF;
	fPtr = mxGetField(sMxData, mIdx, "matBlankLGN_OFF");
	arrMatSize = mxGetDimensions(fPtr);
	sData.matBlankLGN_OFF.setSizeCols(arrMatSize[1]);
	sData.matBlankLGN_OFF.setSizeRows(arrMatSize[0]);
	sData.matBlankLGN_OFF.fillVals(fPtr);

	// read cellLGN_ON;
	sData.cellLGN_ON = putMxCellStimToVecMat(mxGetField(sMxData, mIdx, "cellLGN_ON"));

	// read cellLGN_OFF;
	sData.cellLGN_OFF = putMxCellStimToVecMat(mxGetField(sMxData, mIdx, "cellLGN_OFF"));

	// read vecTrialOris;
	sData.vecTrialOris = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecTrialOris"));
	
	// read vecTrialOriIdx;
	sData.vecTrialOriIdx = getIntVectorFromMx(mxGetField(sMxData, mIdx, "vecTrialOriIdx"));
	
	// read vecStimStartSecs;
	sData.vecStimStartSecs = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecStimStartSecs"));
	
	// read vecTrialEndSecs;
	sData.vecTrialEndSecs = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecTrialEndSecs"));
	
	// read vecThisV;
	sData.vecThisV = getDblVectorFromMx(mxGetField(sMxData, mIdx, "vecThisV"));
	
	// read boolStimPresent;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "boolStimPresent"));
	sData.boolStimPresent = (int)dblPtr[0];

	// read intPrevTrial;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "intPrevTrial"));
	sData.intPrevTrial = (int)dblPtr[0];
	
	// read intTrialT;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "intTrialT"));
	sData.intTrialT = (int)dblPtr[0];
	
	// read intIter;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "intIter"));
	sData.intIter = (int64_t)dblPtr[0];
	
	// read cellSpikeTimesLGN_ON;
	fPtr = mxGetField(sMxData, mIdx, "cellSpikeTimesLGN_ON");
	sData.cellSpikeTimesLGN_ON = putMxCellToVecVec(fPtr);

	// read cellSpikeTimesLGN_OFF;
	fPtr = mxGetField(sMxData, mIdx, "cellSpikeTimesLGN_OFF");
	sData.cellSpikeTimesLGN_OFF = putMxCellToVecVec(fPtr);
	
	// read cellSpikeTimesCortex;
	fPtr = mxGetField(sMxData, mIdx, "cellSpikeTimesCortex");
	sData.cellSpikeTimesCortex = putMxCellToVecVec(fPtr);
	
	// read vecSpikeCounterLGN_ON;
	sData.vecSpikeCounterLGN_ON = getIntVectorFromMx(mxGetField(sMxData, mIdx, "vecSpikeCounterLGN_ON"));
	
	// read vecSpikeCounterLGN_OFF;
	sData.vecSpikeCounterLGN_OFF = getIntVectorFromMx(mxGetField(sMxData, mIdx, "vecSpikeCounterLGN_OFF"));
	
	// read vecSpikeCounterCortex;
	sData.vecSpikeCounterCortex = getIntVectorFromMx(mxGetField(sMxData, mIdx, "vecSpikeCounterCortex"));
	
	// read intPreAllocationSize;
	dblPtr = mxGetPr(mxGetField(sMxData, mIdx, "intPreAllocationSize"));
	sData.intPreAllocationSize = (int64_t)dblPtr[0];


	// return
	return 0;
}


// The gateway function //
// calls input function, simulation function, and assigns output
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	// starting message
	time_t now = time(0);
	mexPrintf("\n  <runSimMEX> \tStarting simulation run [%s\b]\n\n", ctime(&now));
	mexEvalString("pause(eps);");

	// create data structure
	struct sAggregate sData;
	
	//get data
	int intLoadFlag = loadData(sData, prhs[0]);
	now = time(0);
	mexPrintf("  <runSimMEX> \tMATLAB mxArray input processed successfully; loading flag was <%d>\n\t\t\t\t\t... Proceeding with simulation run [%s\b]\n\n", intLoadFlag, ctime(&now));
	mexEvalString("pause(eps);");


	const vector<double> vecIn1 = { 1,2,3 };
	const vector<double> vecIn2 = { 3,2,1 };
	vector<double> vecOut = vPlus(vecIn1, vecIn2);
	/*
	mexPrintf("  <runSimMEX> Plus vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	vecOut = vMinus(vecIn1, vecIn2);
	mexPrintf("  <runSimMEX> Minus vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	vecOut = vDiv(vecIn1, vecIn2);
	mexPrintf("  <runSimMEX> Div vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	vecOut = vMult(vecIn1, vecIn2);
	mexPrintf("  <runSimMEX> Mult vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	vecOut = vLT(vecIn1, vecIn2);
	mexPrintf("  <runSimMEX> LT vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	vecOut = vST(vecIn1, vecIn2);
	mexPrintf("  <runSimMEX> ST vals: [%f %f %f]\n", vecOut[0], vecOut[1], vecOut[2]);
	mexEvalString("pause(eps);");
	*/

	//run simulation
	int intSimFlag = runSimLoop(sData);
	now = time(0);
	mexPrintf("  <runSimMEX> \tSimulation run completed; output flag was <%d>\n\t\t\t\t\t... Transforming data to MATLAB mxArray structure [%s\b]\n\n", intSimFlag, ctime(&now));
	mexEvalString("pause(eps);");
	
	//create output MATLAB structure
	mwIndex mIdx = 0; 
	mwSize dims[2] = { 1, 1 };
	plhs[0] = mxCreateStructArray(2, dims, NR_FIELDS_AGG, field_names_sAggregate);
	mxArray **mxOut = &plhs[0]; //create alias


	//// assign data to fields
	mxArray *mD_TempFieldValue;
	// set dblDeltaT
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = sData.dblDeltaT;
	mxSetField(*mxOut, mIdx, "dblDeltaT", mD_TempFieldValue);

	// set vecOverallT
	mxSetField(*mxOut, mIdx, "vecOverallT", putDblVectorToMx(sData.vecOverallT));
	
	// set matCortConn
	mxSetField(*mxOut, mIdx, "matCortConn", sData.matCortConn.exportToMx());

	// set dblSynSpikeMem;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = sData.dblSynSpikeMem;
	mxSetField(*mxOut, mIdx, "dblSynSpikeMem", mD_TempFieldValue);

	// set vecCortSynType;
	mxSetField(*mxOut, mIdx, "vecCortSynType", putIntVectorToMx(sData.vecCortSynType)); // assign mxArray to field

	// set intCortexCells;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.intCortexCells;
	mxSetField(*mxOut, mIdx, "intCortexCells", mD_TempFieldValue);
	
	// set vecCortDelay;
	mxSetField(*mxOut, mIdx, "vecCortDelay", putDblVectorToMx(sData.vecCortDelay));
	
	// set vecCortConductance;
	mxSetField(*mxOut, mIdx, "vecCortConductance", putDblVectorToMx(sData.vecCortConductance));
	
	// set vecCellThresh;
	mxSetField(*mxOut, mIdx, "vecCellThresh", putDblVectorToMx(sData.vecCellThresh));
	
	// set vecTauPeakByType;
	mxSetField(*mxOut, mIdx, "vecTauPeakByType", putDblVectorToMx(sData.vecTauPeakByType));
	
	// set vecCellV_E;
	mxSetField(*mxOut, mIdx, "vecCellV_E", putDblVectorToMx(sData.vecCellV_E));
	
	// set vecCellV_I;
	mxSetField(*mxOut, mIdx, "vecCellV_I", putDblVectorToMx(sData.vecCellV_I));
	
	// set vecCellV_AHP;
	mxSetField(*mxOut, mIdx, "vecCellV_AHP", putDblVectorToMx(sData.vecCellV_AHP));
	
	// set vecCellV_Leak;
	mxSetField(*mxOut, mIdx, "vecCellV_Leak", putDblVectorToMx(sData.vecCellV_Leak));
	
	// set vecCellCm;
	mxSetField(*mxOut, mIdx, "vecCellCm", putDblVectorToMx(sData.vecCellCm));
	
	// set vecCellG_Leak;
	mxSetField(*mxOut, mIdx, "vecCellG_Leak", putDblVectorToMx(sData.vecCellG_Leak));
	
	// set vecCellG_AHP;
	mxSetField(*mxOut, mIdx, "vecCellG_AHP", putDblVectorToMx(sData.vecCellG_AHP));
	
	// set vecSynConductanceON_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynConductanceON_to_Cort", putDblVectorToMx(sData.vecSynConductanceON_to_Cort));
	
	// set vecSynConductanceOFF_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynConductanceOFF_to_Cort", putDblVectorToMx(sData.vecSynConductanceOFF_to_Cort));
	
	// set vecSynWeightON_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynWeightON_to_Cort", putDblVectorToMx(sData.vecSynWeightON_to_Cort));
	
	// set vecSynWeightOFF_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynWeightOFF_to_Cort", putDblVectorToMx(sData.vecSynWeightOFF_to_Cort));
	
	// set vecSynDelayON_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynDelayON_to_Cort", putDblVectorToMx(sData.vecSynDelayON_to_Cort));
	
	// set vecSynDelayOFF_to_Cort;
	mxSetField(*mxOut, mIdx, "vecSynDelayOFF_to_Cort", putDblVectorToMx(sData.vecSynDelayOFF_to_Cort));
	
	// set matSynConnON_to_Cort;
	mxSetField(*mxOut, mIdx, "matSynConnON_to_Cort", sData.matSynConnON_to_Cort.exportToMx());

	// set matSynConnOFF_to_Cort;
	mxSetField(*mxOut, mIdx, "matSynConnOFF_to_Cort", sData.matSynConnOFF_to_Cort.exportToMx());

	// set matBlankLGN_ON;
	mxSetField(*mxOut, mIdx, "matBlankLGN_ON", sData.matBlankLGN_ON.exportToMx());

	// set matBlankLGN_OFF;
	mxSetField(*mxOut, mIdx, "matBlankLGN_OFF", sData.matBlankLGN_OFF.exportToMx());

	// set cellLGN_ON;
	// matrix (x-y) in vector (t) in vector (stim)
	mD_TempFieldValue = putVecMatToMxCellStim(sData.cellLGN_ON);
	mxSetField(*mxOut, mIdx, "cellLGN_ON", mD_TempFieldValue); // assign mxArray to field


	// set cellLGN_OFF;
	// matrix (x-y) in vector (t) in vector (stim)
	mD_TempFieldValue = putVecMatToMxCellStim(sData.cellLGN_OFF);
	mxSetField(*mxOut, mIdx, "cellLGN_OFF", mD_TempFieldValue); // assign mxArray to field

	
	// set vecTrialOris;
	mxSetField(*mxOut, mIdx, "vecTrialOris", putDblVectorToMx(sData.vecTrialOris));
	
	// set vecTrialOriIdx;
	mxSetField(*mxOut, mIdx, "vecTrialOriIdx", putIntVectorToMx(sData.vecTrialOriIdx));
	
	// set vecStimStartSecs;
	mxSetField(*mxOut, mIdx, "vecStimStartSecs", putDblVectorToMx(sData.vecStimStartSecs));
	
	// set vecTrialEndSecs;
	mxSetField(*mxOut, mIdx, "vecTrialEndSecs", putDblVectorToMx(sData.vecTrialEndSecs));
	
	// set vecThisV;
	mxSetField(*mxOut, mIdx, "vecThisV", putDblVectorToMx(sData.vecThisV));
	
	// set boolStimPresent;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.boolStimPresent;
	mxSetField(*mxOut, mIdx, "boolStimPresent", mD_TempFieldValue);

	// set intPrevTrial;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.intPrevTrial;
	mxSetField(*mxOut, mIdx, "intPrevTrial", mD_TempFieldValue);
	
	// set intTrialT;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.intTrialT;
	mxSetField(*mxOut, mIdx, "intTrialT", mD_TempFieldValue);
	
	// set intIter;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.intIter;
	mxSetField(*mxOut, mIdx, "intIter", mD_TempFieldValue);
	
	// set cellSpikeTimesLGN_ON;
	mD_TempFieldValue = putVecVecToMxCell(sData.cellSpikeTimesLGN_ON);
	mxSetField(*mxOut, mIdx, "cellSpikeTimesLGN_ON", mD_TempFieldValue); // assign mxArray to field

	// set cellSpikeTimesLGN_OFF;
	mD_TempFieldValue = putVecVecToMxCell(sData.cellSpikeTimesLGN_OFF);
	mxSetField(*mxOut, mIdx, "cellSpikeTimesLGN_OFF", mD_TempFieldValue); // assign mxArray to field
	
	// set cellSpikeTimesCortex;
	mD_TempFieldValue = putVecVecToMxCell(sData.cellSpikeTimesCortex);
	mxSetField(*mxOut, mIdx, "cellSpikeTimesCortex", mD_TempFieldValue); // assign mxArray to field

	// set vecSpikeCounterLGN_ON;
	mxSetField(*mxOut, mIdx, "vecSpikeCounterLGN_ON", putIntVectorToMx(sData.vecSpikeCounterLGN_ON));
	
	// set vecSpikeCounterLGN_OFF;
	mxSetField(*mxOut, mIdx, "vecSpikeCounterLGN_OFF", putIntVectorToMx(sData.vecSpikeCounterLGN_OFF));
	
	// set vecSpikeCounterCortex;
	mxSetField(*mxOut, mIdx, "vecSpikeCounterCortex", putIntVectorToMx(sData.vecSpikeCounterCortex));
	
	// set intPreAllocationSize;
	mD_TempFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mD_TempFieldValue) = (double)sData.intPreAllocationSize;
	mxSetField(*mxOut, mIdx, "intPreAllocationSize", mD_TempFieldValue);

	// exit message
	now = time(0);
	mexPrintf("  <runSimMEX> \tData transformation complete! [%s\b]\n\n", ctime(&now));
	mexEvalString("pause(eps);");

	//no return flag
	return;
}

/*
// define aggregate data structure used for input and output
struct sAggregate
{
short intID;
double dblTest;
};


struct sAggregate main(struct sAggregate sIn)
{
// define namespace
using namespace std;

// define variables
constexpr double dblMaxRand = static_cast<double>(UINT32_MAX);

// random values
random_device rd; // Use a hardware entropy source if available, otherwise use PRNG
mt19937 mersenne(rd()); // initialize our mersenne twister with a random seed

// transform to double with [0 - 1] range
uint32_t uintRand = mersenne();
double dblRand = static_cast<double>(uintRand);
dblRand = dblRand / dblMaxRand;

// print random number
cout << dblRand << "\t";

// prepare output
sAggregate sOut;
sOut = sIn;
sOut.dblTest = dblRand;
return sOut;
}

*/

/*
__declspec(dllexport) int __stdcall add2(int num) {
return num + 2;
}

__declspec(dllexport) int __stdcall mul(int num1, int num2) {
return num1 * num2;
}

struct getSimulationRun(struct sIn)
{
// define namespace
using namespace std;

// define variables
constexpr double dblMaxRand = static_cast<double>(UINT32_MAX);

// random values
random_device rd; // Use a hardware entropy source if available, otherwise use PRNG
mt19937 mersenne(rd()); // initialize our mersenne twister with a random seed

// transform to double with [0 - 1] range
uint32_t uintRand = mersenne();
double dblRand = static_cast<double>(uintRand);
dblRand = dblRand / dblMaxRand;

// print output
cout << dblRand << "\t";
}
*/


/*
test: F7 > ctrl-F5

// fixed arrays; remember, it's [0 - N-1] !
int arrArray[5] = { 7, 4, 5 } // initialize fixed array
int arrArray[] = { 7, 4, 5, 9, 11 } // initialize fixed array


// built-in array package, uses size()
std::array<int, 3> myarray = { 9, 7, 5, 3, 1 }; // initialization list
void printLength(const std::array<double, 5> &myarray) // when taking array as input


// for loop with array
const int intArraySize = sizeof(arrArray) / sizeof(arrArray[0]);
for (int intIdx = 0; intIdx < intArraySize; ++intIdx)
{
arrArray[intIdx];
}


// multidimensional array
int arrArray[3][5] = // a 3-element array of 5-element arrays [row col]
{
{ 1, 2, 3, 4, 5 }, // row 0
{ 6, 7, 8, 9, 10 }, // row 1
{ 11, 12, 13, 14, 15 } // row 2
};


// dynamic memory allocation
int *ptr1 = new int (5);
delete ptr1; // return the memory pointed to by ptr to the operating system
ptr1 = nullptr; // set ptr to be a null pointer


// dynamic memory allocation for array
int intLength = 10;
int *arrArray = new int[intLength];
delete[] arrArray; // return the memory pointed to by vecArray to the operating system
arrArray = nullptr; // set ptr to be a null pointer


// array from array package
for (int intVal : arrArray) // returns all values in arrArray iteratively
for (auto &element : arrArray) // returns all values in arrArray iteratively through direct referencing


// vector from vector package
std::vector<int> array2 = { 9, 7, 5, 3, 1 };
array.resize(5); // to resize
stack.reserve(5); // to reserve (but not fill) elements

// structs
struct Employee
{
short id;
int age;
double wage;
};

// references
int &ref = other.something.value1; // ref can now be used in place of other.something.value1


// inputs
int main(int argc, char *argv[])
*/

//DLL in/output
/*
//In .h:
#ifdef BUILD_DLL
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __declspec(dllimport)
#endif

extern "C" // Only if you are using C++ rather than C
{
EXPORT struct __stdcall getSimulationRun(struct sIn);
}


//in .cpp:
extern "C"
{
EXPORT struct __stdcall add2(struct sIn)
{
return sOut;
}
}

*/