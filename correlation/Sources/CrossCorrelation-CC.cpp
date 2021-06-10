/*
	++ Source:
			CrossCorrelation-CC.cpp

	++ Description:
			A class that computes the cross-correlation for two continuous signals (Continuous - Continuous)
			Cross-Correlation is computed in the classical way by counting

	++ History:
			19-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/

#define div /

#include <stdio.h>
#include <math.h>
#include "MyMath.h"
#include "Constants.h"
#include "CrossCorrelation-CC.h"


//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
// Constructor parameters are as follows:
///	@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with laggs from -80 to +80); use the same units as the sampling of your samples
///	@param piTrialLength			- The size of the trial in units of your samples
///	@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL;
CCrossCorrelationCC::CCrossCorrelationCC(int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faCrossCorrelogram   = NULL;
	piIsError			  = 0;
	if(CheckWindowsConsistency(piCorrelationWindow,piTrialLength))
	{
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units, same as the samples in the input vectors
		AllocateData();
	}
	else piIsError = 1;
}

///Destructor
CCrossCorrelationCC::~CCrossCorrelationCC()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

///Checks if the sizes of the windows are consistent
int	CCrossCorrelationCC::CheckWindowsConsistency(int piCorrW,int piTrialLength)
{
	if((piCorrW > piTrialLength) || (1.0*piCorrW*piTrialLength <= 0.0)) return 0;
	return 1;
}


///Allocates data for the cross correlogram
void CCrossCorrelationCC::AllocateData()
{
	if(faCrossCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faCrossCorrelogram = new float[2*iCorrelationWindow+1];
	iaNumberOfSummedCorrelations = new int[2*iCorrelationWindow+1];
}

///Releases the data structures used by the cross-correlation
void CCrossCorrelationCC::CleanupData()
{
	if(faCrossCorrelogram != NULL)		//safety check
	{
		delete[] faCrossCorrelogram;
		faCrossCorrelogram = NULL;
	}
	if(iaNumberOfSummedCorrelations != NULL)
	{
		delete[] iaNumberOfSummedCorrelations;
		iaNumberOfSummedCorrelations = NULL;
	}
}

///Clears the bytes in the cross correlogram
void CCrossCorrelationCC::ResetCrossCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++) 
	{
		faCrossCorrelogram[i] = 0;
		iaNumberOfSummedCorrelations[i] = 0;
	}
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CCrossCorrelationCC::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++) 
	{
		faCrossCorrelogram[i] = NAN;
		iaNumberOfSummedCorrelations[i] = 0;
	}
}

///The function that computes one pass cross-correlation for a pair of trials
void CCrossCorrelationCC::ComputeTrialCrossCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB)
{
	int   i,j;
	int   iLeft,iRight;						//Indexes for the valid windows
	float fSignalACache;
	int   iOffset;
	
	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0) || (piIdxTrialStartA != piIdxTrialStartB)) return;


	//Method with throwing out borders - can create symetry problems in autocorrelograms!!!
	/*
	iLeft  = max(piIdxTrialStartA,piIdxTrialStartB) + iCorrelationWindow;
	iRight = min(piIdxTrialEndA,piIdxTrialEndB) - iCorrelationWindow;
	iTotalWindowSize = 2*iCorrelationWindow+1;


	for(i=iLeft;i<=iRight;i++)
	{
		fSignalACache = faSignalA[i];
		iOffset = i - iCorrelationWindow;		
		for(j=0;j<iTotalWindowSize;j++) faCrossCorrelogram[j] += fSignalACache * faSignalB[j+iOffset];
		iNumberOfSummedCorrelograms++;
	}
	*/

	//Method with checking border values
	for(i=piIdxTrialStartA;i<=piIdxTrialEndA;i++)
	{
		iLeft = max(piIdxTrialStartB,i-iCorrelationWindow); 
		iRight = min(piIdxTrialEndB,i+iCorrelationWindow);
		iOffset = i - iCorrelationWindow;
		fSignalACache = faSignalA[i];

		for(j=iLeft;j<=iRight;j++) 
		{
			faCrossCorrelogram[j - iOffset] += fSignalACache * faSignalB[j];
			iaNumberOfSummedCorrelations[j - iOffset]++;
		}
	}
}

///The function that computes one pass cross-correlation for a pair of trials without checking for borders
void CCrossCorrelationCC::ComputeTrialCrossCorrelationNoBorderChecking(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB)
{
	int   i,j;
	int   iLeft,iRight;						//Indexes for the valid windows
	float fSignalACache;
	int   iOffset;
	
	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0) || (piIdxTrialStartA != piIdxTrialStartB)) return;

	//Method without checking border values
	for(i=piIdxTrialStartA;i<=piIdxTrialEndA;i++)
	{
		iLeft = i-iCorrelationWindow; 
		iRight = i+iCorrelationWindow;
		iOffset = i - iCorrelationWindow;
		fSignalACache = faSignalA[i];

		for(j=iLeft;j<=iRight;j++) 
		{
			faCrossCorrelogram[j - iOffset] += fSignalACache * faSignalB[j];
			iaNumberOfSummedCorrelations[j - iOffset]++;
		}
	}
}

//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		
//Compute the cross-correlation of two vectors of samples
///If the number of samples passed is larger than the size of one trial, then the signals are divided into trials and the results are eventually averaged over trials
//Parameters for ComputeCrossCorrelation
/// @param pfaSamplesA				- The vector holding the samples of the first signal
/// @param piNrSamplesInA			- Number of samples in the first vector
/// @param pfaSamplesB				- The vector holding the samples of the second signal
/// @param piNrSamplesInB			- Number of samples in the second vector
/// @param piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
void CCrossCorrelationCC::ComputeCrossCorrelation(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag)
{
	int i;
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iTrialCountA,iTrialCountB;			//Number of trials for the two Sampless
	int iCentralPeakCount;					//Number of correlation bina added up for lag 0
	
	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL))
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= 0;
	faSignalA		= pfaSamplesA;
	faSignalB		= pfaSamplesB;

	//Compute the trials and call the handlers
	iTrialCountA = (int) ceilf((float)piNrSamplesInA / iTrialLength);
	iTrialCountB = (int) ceilf((float)piNrSamplesInB / iTrialLength);

	for(i=0;i<min(iTrialCountA,iTrialCountB);i++)
	{
		idxStartTrialA = i * iTrialLength;
		idxEndTrialA   = min(piNrSamplesInA,(i+1)*iTrialLength) - 1;
		idxStartTrialB = i * iTrialLength;
		idxEndTrialB   = min(piNrSamplesInB,(i+1)*iTrialLength) - 1;

		ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);
	}

	//Normalize the cross correlogram and deal with bias
	iCentralPeakCount = iaNumberOfSummedCorrelations[iCorrelationWindow];
	switch(piNormalizationFlag)
	{
		case NORM_BIASED: 
				if(iCentralPeakCount)
					for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] /= (float)iCentralPeakCount;
				break;

		case UNNORM_UNBIASED:
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] = faCrossCorrelogram[i] / (float)iaNumberOfSummedCorrelations[i] * (float)iCentralPeakCount;
				break;

		case NORM_UNBIASED:	
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] /= (float)iaNumberOfSummedCorrelations[i];
	}	
}

//Compute the cross-correlation of two vectors of samples without checking for borders
//WARNING: for each trial, the buffer must contain before and after the trial a number of samples equal to the correlation window!
///If the number of samples passed is larger than the size of one trial, then the signals are divided into trials and the results are eventually averaged over trials
//Parameters for ComputeCrossCorrelation
/// @param pfaSamplesA				- The vector holding the samples of the first signal
/// @param piNrSamplesInA			- Number of samples in the first vector
/// @param pfaSamplesB				- The vector holding the samples of the second signal
/// @param piNrSamplesInB			- Number of samples in the second vector
/// @param piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
void CCrossCorrelationCC::ComputeCrossCorrelationNoBorderChecking(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag)
{
	int i;
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iTrialCountA,iTrialCountB;			//Number of trials for the two Sampless
	int iCentralPeakCount;					//Number of correlation bina added up for lag 0
	
	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL))
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= 0;
	faSignalA		= pfaSamplesA;
	faSignalB		= pfaSamplesB;

	//Compute the trials and call the handlers
	iTrialCountA = (int) ceilf((float)piNrSamplesInA / iTrialLength);
	iTrialCountB = (int) ceilf((float)piNrSamplesInB / iTrialLength);

	for(i=0;i<min(iTrialCountA,iTrialCountB);i++)
	{
		idxStartTrialA = i * iTrialLength;
		idxEndTrialA   = min(piNrSamplesInA,(i+1)*iTrialLength) - 1;
		idxStartTrialB = i * iTrialLength;
		idxEndTrialB   = min(piNrSamplesInB,(i+1)*iTrialLength) - 1;

		ComputeTrialCrossCorrelationNoBorderChecking(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);
	}

	//Normalize the cross correlogram and deal with bias
	iCentralPeakCount = iaNumberOfSummedCorrelations[iCorrelationWindow];
	switch(piNormalizationFlag)
	{
		case NORM_BIASED: 
				if(iCentralPeakCount)
					for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] /= (float)iCentralPeakCount;
				break;

		case UNNORM_UNBIASED:
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] = faCrossCorrelogram[i] / (float)iaNumberOfSummedCorrelations[i] * (float)iCentralPeakCount;
				break;

		case NORM_UNBIASED:	
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] /= (float)iaNumberOfSummedCorrelations[i];
	}	
}

//Compute the cross-correlation of two digitized continuous signals, on a partial window of the trial; 
/// @warning ONLY ACCEPTS ONE TRIAL so do not pass signals longer than the size of a single trial! Subsequent trials will be discarded.
//Parameters for ComputeWindowedCrossCorrelationPerTrial
/// @param	pfaSamplesA				- The vector holding the samples of the first signal 
/// @param	piNrSamplesInA			- Number of samples in the first vector
/// @param	pfaSamplesB				- The vector holding the samples of the second signal
/// @param	piNrSamplesInB			- Number of samples in the second vector
/// @param	piFromOffsetInTrial		- The start offset in the trial where the desired window starts
/// @param	piToOffsetInTrial		- The end offset in the trial where the desired window stops
/// @param	piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
void CCrossCorrelationCC::ComputeWindowedCrossCorrelationPerTrial(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int piNormalizationFlag)
{

	int i,j;
	int iTrialCountA,iTrialCountB;		//Number of trials for the two Sampless
	int iLeft,iRight;					//Indexes for the valid windows
	float fSignalACache;
	int iOffset;
	int iCentralPeakCount;				//Number of correlation bina added up for lag 0

	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	if((piFromOffsetInTrial < 0) || (piToOffsetInTrial >= piNrSamplesInA)||(piFromOffsetInTrial >= piNrSamplesInA) || (piToOffsetInTrial < 0)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Compute the trials and call the handlers
	iTrialCountA = (int) ceilf((float)piNrSamplesInA / iTrialLength);
	iTrialCountB = (int) ceilf((float)piNrSamplesInB / iTrialLength);

	if((iTrialCountA <= 0)||(iTrialCountB <= 0)) return;

	//Method with checking border values
	for(i=piFromOffsetInTrial;i<=piToOffsetInTrial;i++)
	{
		iLeft = max(0,i-iCorrelationWindow); 
		iRight = min(piNrSamplesInB-1,i+iCorrelationWindow);
		iOffset = i - iCorrelationWindow;
		fSignalACache = pfaSamplesA[i];

		for(j=iLeft;j<=iRight;j++)
		{
			faCrossCorrelogram[j - iOffset] += fSignalACache * pfaSamplesB[j];
			iaNumberOfSummedCorrelations[j - iOffset]++;
		}
	}

	//Correct for border anomalies
	/*iATrialLength = piToOffsetInTrial-piFromOffsetInTrial+1;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		fNorm = (float)max(0,iATrialLength-abs(i-iCorrelationWindow)) / (float)iATrialLength;
		if(fNorm) faCrossCorrelogram[i] = faCrossCorrelogram[i] / fNorm;
	}*/

	//Normalize the cross correlogram and deal with bias
	iCentralPeakCount = iaNumberOfSummedCorrelations[iCorrelationWindow];
	switch(piNormalizationFlag)
	{
		case NORM_BIASED: 
				if(iCentralPeakCount)
					for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] /= (float)iCentralPeakCount;
				break;

		case UNNORM_UNBIASED:
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] = faCrossCorrelogram[i] / (float)iaNumberOfSummedCorrelations[i] * (float)iCentralPeakCount;
				break;

		case NORM_UNBIASED:	
				for(i=0;i<2*iCorrelationWindow+1;i++) 
					if(iaNumberOfSummedCorrelations[i]) faCrossCorrelogram[i] /= (float)iaNumberOfSummedCorrelations[i];
	}
}


//Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the cross-correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
/// @return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///			The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
float* CCrossCorrelationCC::GetCrossCorrelogram()
{
	return faCrossCorrelogram;
}

//Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter;
///@param piNewCorrelationWindow - Specifies the new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
int CCrossCorrelationCC::ModifyCorrelationWindow(int piNewCorrelationWindow)
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
int CCrossCorrelationCC::ModifyTrialLength(int piNewTrialLength)
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
int CCrossCorrelationCC::ModifyAllParameters(int piNewCorrelationWindow,int piNewTrialLength)
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
/// @return The size of the correlation window
int CCrossCorrelationCC::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
/// @return The trial length in sampling units
int CCrossCorrelationCC::GetTrialLength()
{
	return iTrialLength;
}