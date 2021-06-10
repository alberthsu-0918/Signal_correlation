/*
	++ Source:
			ScaledCorrelation-BC.cpp

	++ Description:
			A class that computes the scaled correlation for a vector of time stamps and a continuous signal (Binary - Continuous)
			For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer

	++ History:
			21-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


#define div /


#include <stdio.h>
#include <math.h>
#include "MyMath.h"
#include "ScaledCorrelation-BC.h"
#include "Constants.h"

//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
//Constructor parameters are as follows:
///@param piScaleWindow			- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
///@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
///@param piTrialLength			- The size of the trial in original sampling units
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
///@warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
CScaledCorrelationBC::CScaledCorrelationBC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faScaledCorrelogram = NULL;
	faPBsCoefficientsSum   = NULL;
	iaPBsCoefficientsCount = NULL;
	iaPBsCoefficientsDistribution = NULL;
	piIsError = 0;
	if(CheckWindowsConsistency(piScaleWindow,piCorrelationWindow,piTrialLength))
	{
		iScaleWindow		 = piScaleWindow;					//The size of the scale window
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units
		bUseFisherZTransform = 0;
		fBinSizeOfCorrCoeffDistribution = (float) BIN_SIZE_OF_COEFF_DISTRIBUTION;					//The size of a bin in the distribution of correlation coefficients
		iPBsCorrCoeffDistributionBins = (int)(2*ceil(1.0 / fBinSizeOfCorrCoeffDistribution)+1);	//Number of bins in the distribution of fi correlation coefficients; they span the interval of [-1..1]
		AllocateData();
	}
	else piIsError = 1;
}

//Destructor
CScaledCorrelationBC::~CScaledCorrelationBC()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

//Checks if the sizes of the windows are consistent
int	CScaledCorrelationBC::CheckWindowsConsistency(int piScaleW,int piCorrW,int piTrialLength)
{
	if((piScaleW > piTrialLength) || (piCorrW > piTrialLength) || (1.0*piCorrW*piScaleW*piTrialLength <= 0.0)) return 0;
	return 1;
}


//Allocates data for the cross correlogram
void CScaledCorrelationBC::AllocateData()
{
	if(faScaledCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faScaledCorrelogram = new float[2*iCorrelationWindow+1];
	faPBsCoefficientsSum = new float[2*iCorrelationWindow+1];
	iaPBsCoefficientsCount = new int[2*iCorrelationWindow+1];
	iaPBsCoefficientsDistribution = new int[iPBsCorrCoeffDistributionBins];
}

//Releases the data structures used by the cross-correlation
void CScaledCorrelationBC::CleanupData()
{
	if(faScaledCorrelogram != NULL)		//safety check
	{
		delete[] faScaledCorrelogram;
		faScaledCorrelogram = NULL;
	}

	if(faPBsCoefficientsSum != NULL)
	{
		delete[] faPBsCoefficientsSum;
		faPBsCoefficientsSum = NULL;
	}

	if(iaPBsCoefficientsCount != NULL)
	{
		delete[] iaPBsCoefficientsCount;
		iaPBsCoefficientsCount = NULL;
	}
	
	if(iaPBsCoefficientsDistribution != NULL)
	{
		delete[] iaPBsCoefficientsDistribution;
		iaPBsCoefficientsDistribution = NULL;
	}
}

//Clears the bytes in the cross correlogram
void CScaledCorrelationBC::ResetScaledCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = 0;
		faPBsCoefficientsSum[i] = 0;
		iaPBsCoefficientsCount[i] = 0;
	}
	for(i=0;i<iPBsCorrCoeffDistributionBins;i++) iaPBsCoefficientsDistribution[i] = 0;
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CScaledCorrelationBC::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = NAN;
		faPBsCoefficientsSum[i] = NAN;
		iaPBsCoefficientsCount[i] = 0;
	}
	for(i=0;i<iPBsCorrCoeffDistributionBins;i++) iaPBsCoefficientsDistribution[i] = 0;
}

//Determines the Point Biserial coefficient of correlation for a pair of segments;
float CScaledCorrelationBC::ComputePBsCorrelation(int* piaBufferA,int piBufferASize,float* pfaBufferB,int piBufferBSize,int piWindowStartTime)					
{
	
	float x0,x1,p,Sx,sumA,meanA;
	int idxB,idxA;

	if((piBufferBSize*piBufferASize == 0) || (piBufferASize > piBufferBSize)) return NAN;

	idxB = 0;
	idxA = 0;

	//Compute x1 unnormalized
	x1 = 0;
	for(idxA=0;idxA<piBufferASize;idxA++) x1 += pfaBufferB[piaBufferA[idxA] - piWindowStartTime];


	//Compute the sum over the signal
	sumA = 0;
	for(idxB=0;idxB<piBufferBSize;idxB++) sumA += pfaBufferB[idxB];

	//Compute x0
	if((piBufferBSize - piBufferASize) > 0) x0 = (sumA - x1) / (float)(piBufferBSize - piBufferASize);
	else x0 = 0;


	//Normalize x1
	x1 /= (float) piBufferASize;

	
	//Compute the mean of the continuous signal
	meanA = sumA / (float)piBufferBSize;


	//Compute the standard deviation of A
	Sx = 0;
	for(idxB=0;idxB<piBufferBSize;idxB++) Sx += (pfaBufferB[idxB] - meanA) * (pfaBufferB[idxB] - meanA);
	if(piBufferBSize > 1) Sx = (float)(sqrt(Sx / (float)(piBufferBSize - 1)));
	else Sx = (float)sqrt(Sx);

	
	//Compute the probability that there is a binary event in the window
	p =  (float)piBufferASize / (float)iScaleWindow;	


	//Compute the correlation coeff and ret
	if(Sx != 0)	return (float) (((x1 - x0) * sqrt(p * (1 - p))) / Sx);
	else return NAN;	
}

//The function that computes one pass scale correlation for a pair of trials
void CScaledCorrelationBC::ComputeTrialScaledCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i;
	int iOffset;							//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;	//The indexes where the current windows start in the first and in the second channel
	int idxWindowEndA,idxWindowEndB;		//The indexes where the current windows end in the first and in the second channel
	int iTrialStartTime,iTrialEndTime;		//The times when the trial starts and when the trial ends
	int iWindowA;							//The number of the current window in the first channel
	int iNewWindowA;						//The number of the current window in the first channel in the parsing of the end of the windows
	int idxA;								//Index to move in the binary channel
	float pbsCoefficientOfCorrelation;		//The Point Biserial coefficient of correlation for a given window
	int iIdxTrialStartB;					//The index where trial B should start; this is adjusted to an index such that adding the offset does not exceed the trial beginning or its end
	int iIdxTrialEndB;						//The index where trial B should end; this is adjusted to an index such that adding the offset does not exceed the trial beginning or its end

	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	iTrialStartTime = piTrialIndex * iTrialLength;
	iTrialEndTime   = (piTrialIndex + 1) * iTrialLength - 1;
	
	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= i - iCorrelationWindow;	//The offset for comparing the two trials

		//Compute the true bounds on trial B
		if((piIdxTrialStartB + iOffset > iTrialEndTime) || (piIdxTrialEndB + iOffset < iTrialStartTime)) continue;	//trial B is not within bounds when we add the offset
		else
		{
			iIdxTrialStartB = min(max(piIdxTrialStartB + iOffset,piIdxTrialStartB),piIdxTrialEndB);
			iIdxTrialEndB = max(min(piIdxTrialEndB + iOffset,piIdxTrialEndB),piIdxTrialStartB);
		}

		if(iIdxTrialStartB > iIdxTrialEndB) continue;

		idxWindowStartA = piIdxTrialStartA;;
		idxWindowStartB = iIdxTrialStartB;
		idxWindowEndA	= piIdxTrialStartA;
		idxWindowEndB	= iIdxTrialStartB;


		//Loop through the trial of time stamps in A
		while(idxWindowStartA <= piIdxTrialEndA)
		{		
			iWindowA = (iaTimeStampsA[idxWindowStartA] - iTrialStartTime) div iScaleWindow;			

			if((iaTimeStampsA[idxWindowStartA] >= iIdxTrialStartB) && (iaTimeStampsA[idxWindowStartA] <= (iIdxTrialEndB - iScaleWindow))) //safe check if we are aligned
			{
				//Now we have to determine the end of the scale windows
				idxA = idxWindowStartA;
				iNewWindowA = iWindowA;

				//Try to advance window A and find it's end
				while((iNewWindowA == iWindowA)&&(idxA <= piIdxTrialEndA))
				{
					idxA++;
					iNewWindowA = (iaTimeStampsA[idxA] - iTrialStartTime) div iScaleWindow;
				}
				idxWindowEndA = idxA - 1; //found the last index in the scale window

				//Compute the bounds of the scale window in the continuous signal
				idxWindowStartB = min(max(iWindowA * iScaleWindow + iOffset + iTrialStartTime,iIdxTrialStartB),iIdxTrialEndB);
				idxWindowEndB	= max(min((iWindowA + 1) * iScaleWindow - 1 + iOffset + iTrialStartTime,iIdxTrialEndB),iIdxTrialStartB);

				//Check if we have a full scale segment for the continuous signal
				if((idxWindowEndB - idxWindowStartB + 1) != iScaleWindow) 
				{
					idxWindowStartA = idxA;
					continue;
				}
			
				//Whe got the windows so we can compute the correlation
				pbsCoefficientOfCorrelation = ComputePBsCorrelation((iaTimeStampsA+idxWindowStartA),idxWindowEndA-idxWindowStartA+1,(faSignalB+idxWindowStartB),idxWindowEndB-idxWindowStartB+1,iWindowA * iScaleWindow + iTrialStartTime);
				
				if(pbsCoefficientOfCorrelation != NAN)
				{					
					if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += pbsCoefficientOfCorrelation;
					else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+pbsCoefficientOfCorrelation/1.12)/(1.0-pbsCoefficientOfCorrelation/1.12)));
					
					if((pbsCoefficientOfCorrelation >= -1)&&(pbsCoefficientOfCorrelation <= 1)) iaPBsCoefficientsDistribution[round((float)((pbsCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;
					iaPBsCoefficientsCount[i]++;
				}
				
				//Now start some new window
				idxWindowStartA = idxA;
			} 
			else idxWindowStartA++;
		}
	}
}

//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		

//Compute the scaled correlation of one binary signal with one continuos signal
//Parameters for ComputeScaledCorrelation
///@param piaTimeStampsA			- The vector holding the time stamps of the binary variable
///@param piNrTimeStampsInA			- Number of time stamps in the first vector
///@param pfaSamplesB				- The vector holding samples of the continuous signal
///@param piNrSamplesInB			- Number of samples in the second vector
///@param pbUseFisherZTransform		- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
///@warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match! \n
///			The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
void CScaledCorrelationBC::ComputeScaledCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA;					//The number of the current trial on channel A
	int iPrevTrialA;						//The number of the previous trial on channel A
	int isEndOfTrialA;					//Flag to determine whether we are at the end of a trial in channel A

	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrSamplesInB <= 0) || (piNrTimeStampsInA <= 0) || (pfaSamplesB == NULL) || (piaTimeStampsA == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
	faSignalB		= pfaSamplesB;
	iaTimeStampsA	= piaTimeStampsA;
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

		//If we have a trial detected then we have to compute the Point Biserial correlation coefficient scores
		if(isEndOfTrialA)
		{
			//We identify the trial in the continuous signal
			idxStartTrialB = iPrevTrialA * iTrialLength;
			idxEndTrialB   = min((iPrevTrialA + 1) * iTrialLength-1,piNrSamplesInB);

			//Check the end of the continuous signal
			if((idxStartTrialB > idxEndTrialB) || ((idxEndTrialB - idxStartTrialB < iScaleWindow))) break;  

			//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
			ComputeTrialScaledCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iPrevTrialA);

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
		if((idxStartTrialB < idxEndTrialB) || (idxEndTrialB - idxStartTrialB >= iScaleWindow))
		{
			//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
			ComputeTrialScaledCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iCurrentTrialA);
		}
	}


	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPBsCoefficientsCount[i] > 0)
			{
				faPBsCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaPBsCoefficientsCount[i]);
			}
			else
			{
				faPBsCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPBsCoefficientsCount[i] > 0)
			{
				faPBsCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaPBsCoefficientsCount[i]));
			}
			else
			{
				faPBsCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
}


//Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
///@return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///		   The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the original signals.
float* CScaledCorrelationBC::GetScaledCrossCorrelogram()
{
	return faScaledCorrelogram;
}

//Returns the sum of valid Point Biserial coefficients of correlation for each bin of the correlogram
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PBsCoefficientSums[i] / PBsCoefficientCounts[i].
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
float* CScaledCorrelationBC::GetPBsCoefficientSums()
{
	return faPBsCoefficientsSum;
}

//Returns how many Point Biserial coefficients of correlation have been averaged for each bin of the scaled cross correlogram
///@return A pointer to the internal buffer containing the counts of Point Biserial correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of valid Point Biserial coefficients, that were averaged to compute the scaled correlogram at lag 0, is at position iCorrelationWindow
int* CScaledCorrelationBC::GetPBsCoefficientCounts()
{
	return iaPBsCoefficientsCount;
}

//Returns the distribution of coefficients of correlation
///@return A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
int* CScaledCorrelationBC::GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize)
{
	piNumberOfBins = iPBsCorrCoeffDistributionBins;
	pfBinSize = fBinSizeOfCorrCoeffDistribution;
	return iaPBsCoefficientsDistribution;
}

//Changes the scale window for the scaled correlation; 
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
int CScaledCorrelationBC::ModifyScaleWindow(int piNewScale)
{
	SetCorrelogramToNAN();
	if(piNewScale <= 0) return -1;
	if(piNewScale > iTrialLength) return -2;
	iScaleWindow = piNewScale;
	ResetScaledCorrelogram();
	return 0;
}

//Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter; 
///@param piNewCorrelationWindow - The new size of the correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial
int CScaledCorrelationBC::ModifyCorrelationWindow(int piNewCorrelationWindow)
{
	SetCorrelogramToNAN();
	if(piNewCorrelationWindow <= 0) return -1;
	if(piNewCorrelationWindow > iTrialLength) return -2;
	CleanupData();
	iCorrelationWindow = piNewCorrelationWindow;
	AllocateData();
	ResetScaledCorrelogram();
	return 0;
}

//Set the length of the trial in original sampling units
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the scale window; -2 if the trial length is smaller than the correlation window
int CScaledCorrelationBC::ModifyTrialLength(int piNewTrialLength)
{
	SetCorrelogramToNAN();
	if(piNewTrialLength < iScaleWindow) return -1;
	if(piNewTrialLength < iCorrelationWindow) return -2;
	iTrialLength = piNewTrialLength;
	ResetScaledCorrelogram();
	return 0;
}

//Sets all parameters at once
///@param piNewScale - The new size of the scale segment
///@param piNewCorrelationWindow - The new size of the correlation window
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success, -1 if the new scale window is too small; -2 if the new correlation window is too small; -3 if the new scale window is larger than the trial size; -4 if the trial length is smaller than the correlation window
int CScaledCorrelationBC::ModifyAllParameters(int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
{
	SetCorrelogramToNAN();
	if(piNewScale <= 0) return -1;
	if(piNewCorrelationWindow <= 0) return -2;
	if(piNewScale > piNewTrialLength) return -3;
	if(piNewCorrelationWindow > piNewTrialLength) return -4;

	iScaleWindow = piNewScale;
	CleanupData();
	iCorrelationWindow = piNewCorrelationWindow;
	AllocateData();
	iTrialLength = piNewTrialLength;	
	ResetScaledCorrelogram();

	return 0;
}

//Get the size of the current scale window
///@return Returns the size of the scale segment window
int CScaledCorrelationBC::GetScaleWindow()
{
	return iScaleWindow;
}

//Get the size of the correlation window
///@return Returns the size of the correlation window
int CScaledCorrelationBC::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
///@return Returns the length of the trial in original sampling units
int CScaledCorrelationBC::GetTrialLength()
{
	return iTrialLength;
}

//Get the number of bins of the distribution of coefficients of correlation
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
int CScaledCorrelationBC::GetDistributionOfCorrelationCoefficientsBinNr()
{
	return iPBsCorrCoeffDistributionBins;
}

//Get the size of a bin of the distribution of coefficients of correlation
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
float CScaledCorrelationBC::GetDistributionOfCorrelationCoefficientsBinSize()
{
	return fBinSizeOfCorrCoeffDistribution;
}

