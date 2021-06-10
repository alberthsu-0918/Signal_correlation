/*
	++ Source:
			CrossCorrelation-BC.cpp

	++ Description:
			A class that computes the cross-correlation for a vector of time stamps and a continuous signal (Binary - Continuous)
			Cross-Correlation is computed in the classical way by counting
			For pairs of spikes and LFP signals in neuroscience, this measure is also known as "spike-triggered LFP average"

	++ History:
			22-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/

#define div /


#include <stdio.h>
#include <math.h>
#include "MyMath.h"
#include "Constants.h"
#include "CrossCorrelation-BC.h"

//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
//Constructor parameters are as follows:
/// @param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use the same units as the sampling units
/// @param piTrialLength		- The size of the trial in original sampling units
/// @param piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL;
/// @warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
CCrossCorrelationBC::CCrossCorrelationBC(int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faCrossCorrelogram = NULL;
	piIsError = 0;
	if(CheckWindowsConsistency(piCorrelationWindow,piTrialLength))
	{
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units
		AllocateData();
	}
	else piIsError = 1;
}

//Destructor
CCrossCorrelationBC::~CCrossCorrelationBC()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

//Checks if the sizes of the windows are consistent
int	CCrossCorrelationBC::CheckWindowsConsistency(int piCorrW,int piTrialLength)
{
	if((piCorrW > piTrialLength) || (1.0*piCorrW*piTrialLength <= 0.0)) return 0;
	return 1;
}


//Allocates data for the cross correlogram
void CCrossCorrelationBC::AllocateData()
{
	if(faCrossCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faCrossCorrelogram = new float[2*iCorrelationWindow+1];
	iaNrSummedBins = new int[2*iCorrelationWindow+1];
}

//Releases the data structures used by the cross-correlation
void CCrossCorrelationBC::CleanupData()
{
	if(faCrossCorrelogram != NULL)		//safety check
	{
		delete[] faCrossCorrelogram;
		faCrossCorrelogram = NULL;
		delete[] iaNrSummedBins;
		iaNrSummedBins = NULL;
	}
}

//Clears the bytes in the cross correlogram
void CCrossCorrelationBC::ResetCrossCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faCrossCorrelogram[i] = 0;
		iaNrSummedBins[i] = 0;
	}
	iNumberOfSummedCorrelograms = 0;
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CCrossCorrelationBC::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++) 
	{
		faCrossCorrelogram[i] = NAN;
		iaNrSummedBins[i] = 0;
	}
	iNumberOfSummedCorrelograms = 0;
}

//The function that computes one pass cross-correlation for a pair of trials
void CCrossCorrelationBC::ComputeTrialCrossCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB)
{
	int i;
	int iOffset;							//The current offset for computing the cross correlogram
	int idxA;								//Index to move in the binary channel
	int iLeft,iRight;
	float fNorm;

	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	//For each time stamp in A
	for(idxA=piIdxTrialStartA;idxA<=piIdxTrialEndA;idxA++)
	{
		iLeft = max(piIdxTrialStartB,iaTimeStampsA[idxA] - iCorrelationWindow);
		iRight = min(piIdxTrialEndB,iaTimeStampsA[idxA] + iCorrelationWindow);

		//Add the values to the average
		for(i=iLeft;i<=iRight;i++)
		{
			iOffset = i - iaTimeStampsA[idxA] + iCorrelationWindow;
			faCrossCorrelogram[iOffset] += faSignalB[i];
			iaNrSummedBins[iOffset]++;
		}
		iNumberOfSummedCorrelograms++;
	}
}

//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		

//Compute the cross-correlation of one binary signal with one continuos signal
//Parameters for ComputeCrossCorrelation
///	 @param	piaTimeStampsA			- The vector holding the time stamps of the binary variable
///  @param	piNrTimeStampsInA		- Number of time stamps in the first vector
///  @param	pfaSamplesB				- The vector holding samples of the continuous signal
///  @param	piNrSamplesInB			- Number of samples in the second vector
///	 @param	pbNormalizeCorrelogram	- Set to 1 to normalize the correlogram; ; normalization is made by dividing to the number of time stamps used to compute the correlogram
///  @warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match! \n
///			  The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
void CCrossCorrelationBC::ComputeCrossCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbNormalizeCorrelogram)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA;					//The number of the current trial on channel A
	int iPrevTrialA;						//The number of the previous trial on channel A
	int isEndOfTrialA;					//Flag to determine whether we are at the end of a trial in channel A

	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrSamplesInB <= 0) || (piNrTimeStampsInA <= 0) || (pfaSamplesB == NULL) || (piaTimeStampsA == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	iaTimeStampsA	= piaTimeStampsA;
	faSignalB		= pfaSamplesB;
	idxB			= 0;
	idxA			= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= 0;
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	isEndOfTrialA	= 0;

	iCurrentTrialA	= piaTimeStampsA[idxA] div iTrialLength;;
	iPrevTrialA		= iCurrentTrialA;
	
	//Main loop that crosses over the data; we parse the time stamps and find matching windows in the continuous signal
	while(idxA < piNrTimeStampsInA)
	{
		//Compute the current trial number
		iCurrentTrialA = piaTimeStampsA[idxA] div iTrialLength;

		//If the trial has just changed then we have parsed a whole trial so we can compute its bounds
		if(iCurrentTrialA != iPrevTrialA) 
		{
			idxEndTrialA  = idxA - 1;
			isEndOfTrialA = 1;
		}

		//If we have a trial detected then we have to compute cross-correlation for the trial
		if(isEndOfTrialA)
		{
			//We identify the trial in the continuous signal
			idxStartTrialB = iPrevTrialA * iTrialLength;
			idxEndTrialB   = min((iPrevTrialA + 1) * iTrialLength-1,piNrSamplesInB);

			//Check the end of the continuous signal
			if((idxStartTrialB > idxEndTrialB) || ((idxEndTrialB - idxStartTrialB < iCorrelationWindow))) break;  

			//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
			ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);

			//Last operation is to start a new trial
			idxStartTrialA = idxA;
			isEndOfTrialA  = 0;
			iPrevTrialA	   = iCurrentTrialA;
		}


		//Advance the index
		if(isEndOfTrialA)			//If we were at the and of the trial, now after the increment we are no more!
		{	
			iPrevTrialA = iCurrentTrialA;
			isEndOfTrialA = 0;
			idxStartTrialA = idxA;
		}
		idxA++;						//Advance index in the time stamp buffer
	}


	//Process the last trial if we need to
	if(idxA == piNrTimeStampsInA)
	{
		idxEndTrialA = idxA - 1;

		//We identify the trial in the continuous signal
		idxStartTrialB = iCurrentTrialA * iTrialLength;
		idxEndTrialB   = min((iCurrentTrialA + 1) * iTrialLength-1,piNrSamplesInB - 1);

		//Check the end of the continuous signal
		if((idxStartTrialB < idxEndTrialB) || (idxEndTrialB - idxStartTrialB >= iCorrelationWindow))
		{
			//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
			ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);
		}
	}

	//Normalize the cross correlogram
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		if(iaNrSummedBins[i] <= 0) faCrossCorrelogram[i] = NAN;
		else
		{
			if(pbNormalizeCorrelogram == 1) faCrossCorrelogram[i] /= (float)iaNrSummedBins[i];
			else faCrossCorrelogram[i] = faCrossCorrelogram[i] = faCrossCorrelogram[i] * (float)iNumberOfSummedCorrelograms / (float)iaNrSummedBins[i]; //correct for border anomalies
		}
	}
}


//Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow.
/// @return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///			The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
float* CCrossCorrelationBC::GetCrossCorrelogram()
{
	return faCrossCorrelogram;
}

//Set the size of the correlation window; for example for a window of -100..+100 pass 100 as a parameter; 
///@param piNewCorrelationWindow - Specifies the new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
int CCrossCorrelationBC::ModifyCorrelationWindow(int piNewCorrelationWindow)
{
	SetCorrelogramToNAN();
	if(piNewCorrelationWindow <= 0) return -1;
	if(piNewCorrelationWindow > iTrialLength) return -2;
	CleanupData();
	iCorrelationWindow = piNewCorrelationWindow;
	AllocateData();
	ResetCrossCorrelogram();
	return 0;
}

//Set the length of the trial in original sampling units
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the correlation window
int CCrossCorrelationBC::ModifyTrialLength(int piNewTrialLength)
{
	SetCorrelogramToNAN();
	if(piNewTrialLength < iCorrelationWindow) return -1;
	iTrialLength = piNewTrialLength;
	ResetCrossCorrelogram();
	return 0;
}

//Sets all parameters at once
///@param piNewCorrelationWindow - The new correlation window
///@param piNewTrialLength - The new length of the trial 
///@return Returns 0 on success, -1 if the new correlation window is too small; -2 if the trial length is smaller than the correlation window
int CCrossCorrelationBC::ModifyAllParameters(int piNewCorrelationWindow,int piNewTrialLength)
{
	SetCorrelogramToNAN();
	if(piNewCorrelationWindow <= 0) return -1;
	if(piNewCorrelationWindow > piNewTrialLength) return -2;

	CleanupData();
	iCorrelationWindow = piNewCorrelationWindow;
	AllocateData();
	iTrialLength = piNewTrialLength;	
	ResetCrossCorrelogram();

	return 0;
}

//Get the size of the correlation window
///@return The size of the correlation window
int CCrossCorrelationBC::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
///@return The length of the trial in sampling units
int CCrossCorrelationBC::GetTrialLength()
{
	return iTrialLength;
}