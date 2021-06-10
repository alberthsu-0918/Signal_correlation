/*
	++ Source:
			ScaledCorrelation-CC.cpp

	++ Description:
			A class that computes the scaled correlation for two continuous signals (Continuous - Continuous)
			For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer

	++ History:
			17-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/

#define div /

#include <stdio.h>
#include <math.h>
#include "Constants.h"
#include "MyMath.h"
#include "ScaledCorrelation-CC.h"


//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
//Constructor parameters are as follows:
///@param piScaleWindow			- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your signals
///@param piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use again the same units as the sampling of your signals
///@param piTrialLength			- The size of the trial in sampling units of your signals
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
CScaledCorrelationCC::CScaledCorrelationCC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faScaledCorrelogram = NULL;
	faPearsonCoefficientsSum   = NULL;
	iaPearsonCoefficientsCount = NULL;
	iaPearsonCoefficientsDistribution = NULL;
	piIsError = 0;
	if(CheckWindowsConsistency(piScaleWindow,piCorrelationWindow,piTrialLength))
	{
		iScaleWindow		 = piScaleWindow;					//The size of the scale window
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units, same as the samples in the input vectors
		bUseFisherZTransform = 0;
		fBinSizeOfCorrCoeffDistribution = (float) BIN_SIZE_OF_COEFF_DISTRIBUTION;					//The size of a bin in the distribution of correlation coefficients
		iPearsonCorrCoeffDistributionBins = (int)(2*ceil(1.0 / fBinSizeOfCorrCoeffDistribution)+1);	//Number of bins in the distribution of Pearson correlation coefficients; they span the interval of [-1..1]
		AllocateData();
	}
	else piIsError = 1;
}

//Destructor
CScaledCorrelationCC::~CScaledCorrelationCC()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

//Checks if the sizes of the windows are consistent
int	CScaledCorrelationCC::CheckWindowsConsistency(int piScaleW,int piCorrW,int piTrialLength)
{
	if((piScaleW > piTrialLength) || (piCorrW > piTrialLength) || (1.0*piCorrW*piScaleW*piTrialLength <= 0.0)) return 0;
	return 1;
}


//Allocates data for the cross correlogram
void CScaledCorrelationCC::AllocateData()
{
	if(faScaledCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faScaledCorrelogram = new float[2*iCorrelationWindow+1];
	faPearsonCoefficientsSum = new float[2*iCorrelationWindow+1];
	iaPearsonCoefficientsCount = new int[2*iCorrelationWindow+1];
	iaPearsonCoefficientsDistribution = new int[iPearsonCorrCoeffDistributionBins];
}

//Releases the data structures used by the cross correlation
void CScaledCorrelationCC::CleanupData()
{
	if(faScaledCorrelogram != NULL)		//safety check
	{
		delete[] faScaledCorrelogram;
		faScaledCorrelogram = NULL;
	}
	if(faPearsonCoefficientsSum != NULL)
	{
		delete[] faPearsonCoefficientsSum;
		faPearsonCoefficientsSum = NULL;
	}
	if(iaPearsonCoefficientsCount != NULL)
	{
		delete[] iaPearsonCoefficientsCount;
		iaPearsonCoefficientsCount = NULL;
	}
	if(iaPearsonCoefficientsDistribution != NULL)
	{
		delete[] iaPearsonCoefficientsDistribution;
		iaPearsonCoefficientsDistribution = NULL;
	}
}

//Clears the bytes in the cross correlogram
void CScaledCorrelationCC::ResetScaledCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = 0;
		faPearsonCoefficientsSum[i] = 0;
		iaPearsonCoefficientsCount[i] = 0;
	}
	for(i=0;i<iPearsonCorrCoeffDistributionBins;i++) iaPearsonCoefficientsDistribution[i] = 0;
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CScaledCorrelationCC::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = NAN;
		faPearsonCoefficientsSum[i] = NAN;
		iaPearsonCoefficientsCount[i] = 0;
	}
	for(i=0;i<iPearsonCorrCoeffDistributionBins;i++) iaPearsonCoefficientsDistribution[i] = 0;
}

//Determines the Pearson coefficient of correlation for a pair of segments;
float CScaledCorrelationCC::ComputePearsonCorrelation(float* pfaBufferA,int piBufferASize,float* pfaBufferB,int piBufferBSize)
{
	//Classic way
	int i;
	double sumA,sumB,sumSqrA,sumSqrB,sumAB;
	double fDenominator;
	int iWindowSize;

	if((piBufferASize <= 0) || (piBufferBSize <= 0)) return NAN;

	sumA	= 0;
	sumB	= 0;
	sumSqrA = 0;
	sumSqrB = 0;
	sumAB	= 0;

	//Compute the correlation coeff and ret
	iWindowSize = min(piBufferASize,piBufferBSize);
	for(i=0;i<iWindowSize;i++)
	{
		sumA	+= pfaBufferA[i];
		sumSqrA += pfaBufferA[i] * pfaBufferA[i];
		sumB	+= pfaBufferB[i];
		sumSqrB += pfaBufferB[i] * pfaBufferB[i];
		sumAB   += pfaBufferA[i] * pfaBufferB[i];
	}

	fDenominator = ((double)iWindowSize * sumSqrA - sumA * sumA) * ((double)iWindowSize * sumSqrB - sumB * sumB);

	if(fDenominator <= 0) return NAN;
	fDenominator = sqrt(fDenominator);

	return (float) (((double)iWindowSize * sumAB - sumA*sumB) / fDenominator);
}

//The function that computes one pass scale correlation for a pair of trials - always segments relative to the beginning of signal A (reference); the autocorrelation is NOT symmetric; may discard first and last windows when partially filled with signals
void CScaledCorrelationCC::ComputeTrialScaledCorrelation_Asymmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i,j;
	int iOffset;									//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;			//The indexes where the current windows start in the first and in the second channel
	float fCoefficientOfCorrelation;				//The Pearson coeffcient of correlation for a given window
	int iWindowCountA;								//Number of windows for the first signal
	int idxLeftBoundA,idxLeftBoundB;				//The left bound on the common part of the two trials
	int idxRightBoundA,idxRightBoundB;				//The right bound on the common part of the two trials


	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		//Compute the real bounds of the trials, taking the offset into account
		idxLeftBoundA  = piIdxTrialStartA;
		idxRightBoundA = piIdxTrialEndA;
		idxLeftBoundB  = min(max(piIdxTrialStartB,piIdxTrialStartB - iOffset),piIdxTrialEndB);
		idxRightBoundB = max(min(piIdxTrialEndB,piIdxTrialEndB - iOffset),piIdxTrialStartB);
		
		if((idxLeftBoundA > idxRightBoundA)||(idxLeftBoundB > idxRightBoundB)) continue;

		//Compute the count of windows for the two trials
		iWindowCountA = (int) floorf((float)(idxRightBoundA - idxLeftBoundA + 1) / (float)iScaleWindow);
		
		//Compute the scaled correlation for each window
		for(j=0;j<iWindowCountA;j++)
		{
			idxWindowStartA = idxLeftBoundA + j*iScaleWindow;
			idxWindowStartB = piIdxTrialStartB + j*iScaleWindow;

			//Find matching window in signal B
			if((idxWindowStartB >= idxLeftBoundB) && ((idxWindowStartB + iScaleWindow) < idxRightBoundB))
			{						
				fCoefficientOfCorrelation = ComputePearsonCorrelation(faSignalA + idxWindowStartA,iScaleWindow,faSignalB + idxWindowStartB - iOffset,iScaleWindow);

				if((fCoefficientOfCorrelation != NAN) && (fCoefficientOfCorrelation >= -1) && (fCoefficientOfCorrelation <= 1)) //consistency check
				{
					if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fCoefficientOfCorrelation;
					else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fCoefficientOfCorrelation/1.12)/(1.0-fCoefficientOfCorrelation/1.12)));
				
					iaPearsonCoefficientsDistribution[round((float)((fCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;					
					iaPearsonCoefficientsCount[i]++;
				}
			}
		}
	}
}

//The function that computes one pass scale correlation for a pair of trials - segments the signals such that segments are aligned to signal B for right shift and signal A for left shift; the autocorrelation is symmetric
void CScaledCorrelationCC::ComputeTrialScaledCorrelation_Symmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i,j;
	int iOffset;										//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;				//The indexes where the current windows start in the first and in the second channel
	float fCoefficientOfCorrelation;					//The Pearson coeffcient of correlation for a given window
	int iWindowCountA,iWindowCountB,iWindowCount;	//Number of windows for the two trials at a given offset
	int idxLeftBoundA,idxLeftBoundB;					//The left bound on the common part of the two trials
	int idxRightBoundA,idxRightBoundB;				//The right bound on the common part of the two trials


	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		//Compute the real bounds of the trials, taking the offset into account
		idxLeftBoundA  = min(max(piIdxTrialStartA,piIdxTrialStartB + iOffset),piIdxTrialEndA);
		idxRightBoundA = max(min(piIdxTrialEndA,piIdxTrialEndB + iOffset),piIdxTrialStartA);
		idxLeftBoundB  = min(max(piIdxTrialStartB,piIdxTrialStartB - iOffset),piIdxTrialEndB);
		idxRightBoundB = max(min(piIdxTrialEndB,piIdxTrialEndB - iOffset),piIdxTrialStartB);
		
		if((idxLeftBoundA > idxRightBoundA)||(idxLeftBoundB > idxRightBoundB)) continue;

		//Compute the count of windows for the two trials
		iWindowCountA = (int) floorf((float)(idxRightBoundA - idxLeftBoundA + 1) / (float)iScaleWindow);
		iWindowCountB = (int) floorf((float)(idxRightBoundB - idxLeftBoundB + 1) / (float)iScaleWindow);
		iWindowCount  = min(iWindowCountA,iWindowCountB);
		
		//Compute the scaled correlation for each window
		for(j=0;j<iWindowCount;j++)
		{
			idxWindowStartA = idxLeftBoundA + j*iScaleWindow;
			idxWindowStartB = idxLeftBoundB + j*iScaleWindow;

			fCoefficientOfCorrelation = ComputePearsonCorrelation(faSignalA + idxWindowStartA,iScaleWindow,faSignalB + idxWindowStartB,iScaleWindow);

			if((fCoefficientOfCorrelation != NAN) && (fCoefficientOfCorrelation >= -1) && (fCoefficientOfCorrelation <= 1)) //consistency check
			{
				if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fCoefficientOfCorrelation;
				else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fCoefficientOfCorrelation/1.12)/(1.0-fCoefficientOfCorrelation/1.12)));
				
				iaPearsonCoefficientsDistribution[round((float)((fCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;					
				iaPearsonCoefficientsCount[i]++;
			}
		}
	}
}


//The function that computes one pass scale correlation for a pair of trials; there is no border cheking on the start and end of the trial - one must provide longer buffers by at least the size of the correlation window; the function provides symmetric autocorrelation
void CScaledCorrelationCC::ComputeTrialScaledCorrelationNoBorderChecking(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i,j;
	int iOffset;									//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;			//The indexes where the current windows start in the first and in the second channel
	float fCoefficientOfCorrelation;				//The Pearson coeffcient of correlation for a given window
	int iWindowCount;								//Number of windows for the two trials
	int idxLeftBoundA,idxLeftBoundB;				//The left bound on the common part of the two trials
	int idxRightBoundA,idxRightBoundB;				//The right bound on the common part of the two trials


	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	//Compute the matching parts of the trials
	idxLeftBoundA  = max(piIdxTrialStartA,piIdxTrialStartB);
	idxRightBoundA = min(piIdxTrialEndA,piIdxTrialEndB);
	idxLeftBoundB  = idxLeftBoundA;
	idxRightBoundB = idxLeftBoundB;

		
	if((idxLeftBoundA > idxRightBoundA)||(idxLeftBoundB > idxRightBoundB)) return;

	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		//Compute the count of windows for the two trials
		iWindowCount = (int) floorf((float)(idxRightBoundA - idxLeftBoundA + 1) / (float)iScaleWindow);
		
		//Compute the scaled correlation for each window
		for(j=0;j<iWindowCount;j++)
		{
			idxWindowStartA = idxLeftBoundA + j*iScaleWindow;
			idxWindowStartB = idxLeftBoundB + j*iScaleWindow + iOffset;

			fCoefficientOfCorrelation = ComputePearsonCorrelation(faSignalA + idxWindowStartA,iScaleWindow,faSignalB + idxWindowStartB,iScaleWindow);

			if((fCoefficientOfCorrelation != NAN) && (fCoefficientOfCorrelation >= -1) && (fCoefficientOfCorrelation <= 1)) //consistency check
			{
				if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fCoefficientOfCorrelation;
				else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fCoefficientOfCorrelation/1.12)/(1.0-fCoefficientOfCorrelation/1.12)));
				
				iaPearsonCoefficientsDistribution[round((float)((fCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;					
				iaPearsonCoefficientsCount[i]++;
			}
		}
	}
}


//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		

//Compute the scaled correlation of two digitized continuous signals
//Parameters for ComputeScaledCorrelation
///@param	pfaSamplesA				- The vector holding the first signal's input samples
///@param	piNrSamplesInA			- Number of samples in the first vector
///@param	pfaSamplesB				- The vector holding the second signal's input samples
///@param	piNrSamplesInB			- Number of samples in the second vector
///@param	pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
void CScaledCorrelationCC::ComputeScaledCorrelation(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	int i;
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iTrialCountA,iTrialCountB;		//Number of trials for the two Sampless
	
	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
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

		ComputeTrialScaledCorrelation_Symmetric(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,i);
	}



	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaPearsonCoefficientsCount[i]);
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaPearsonCoefficientsCount[i]));
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
}

///Compute the scaled correlation of two digitized continuous signals without checking for borders at the end of the trials
//WARNING: for each trial, the buffer must contain before and after the trial a number of samples equal to the correlation window!
//Parameters for ComputeScaledCorrelation
///@param	pfaSamplesA				- The vector holding the first signal's input samples
///@param	piNrSamplesInA			- Number of samples in the first vector
///@param	pfaSamplesB				- The vector holding the second signal's input samples
///@param	piNrSamplesInB			- Number of samples in the second vector
///@param	pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
void CScaledCorrelationCC::ComputeScaledCorrelationNoBorderChecking(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	int i;
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iTrialCountA,iTrialCountB;		//Number of trials for the two Sampless
	
	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
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

		ComputeTrialScaledCorrelationNoBorderChecking(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,i);
	}



	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaPearsonCoefficientsCount[i]);
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaPearsonCoefficientsCount[i]));
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
}

//Compute the scaled correlation of two digitized continuous signals, on a partial window of the trial; 
///@warning ONLY ACCEPTS ONE TRIAL!!! so do not pass signals longer than the size of a single trial! Subsequent trials will be discarded.
//Parameters for ComputeWindowedScaledCorrelationPerTrial
///@param	pfaSamplesA				- The vector holding the first signal's input samples
///@param	piNrSamplesInA			- Number of samples in the first vector
///@param	pfaSamplesB				- The vector holding the second signal's input samples
///@param	piNrSamplesInB			- Number of samples in the second vector
///@param	piFromOffsetInTrial		- The start offset in the trial where the desired window starts
///@param	piToOffsetInTrial		- The end offset in the trial where the desired window stops
///@param	pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
void CScaledCorrelationCC::ComputeWindowedScaledCorrelationPerTrial(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform)
{
	int i,j;
	int idxStartTrialA,idxEndTrialA;					//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;					//Index of the start of the trial and end of the trial on channel B
	int iTrialCountA,iTrialCountB;						//Number of trials for the two Sampless
	int iOffset;										//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;				//The indexes where the current windows start in the first and in the second channel
	float fCoefficientOfCorrelation;					//The Pearson coeffcient of correlation for a given window
	int iWindowCountA,iWindowCountB,iWindowCount;		//Number of windows for the two trials at a given offset
	int idxLeftBoundA,idxLeftBoundB;					//The left bound on the common part of the two trials
	int idxRightBoundA,idxRightBoundB;					//The right bound on the common part of the two trials
	
	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrSamplesInA <= 0) || (piNrSamplesInB <= 0) || (pfaSamplesA == NULL) || (pfaSamplesB == NULL) || ((piToOffsetInTrial - piFromOffsetInTrial + 1) < iScaleWindow)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= 0;
	faSignalA		= pfaSamplesA;
	faSignalB		= pfaSamplesB;

	//Compute the trial counts and check for consistency
	iTrialCountA = (int) ceilf((float)piNrSamplesInA / iTrialLength);
	iTrialCountB = (int) ceilf((float)piNrSamplesInB / iTrialLength);

	if((iTrialCountA <= 0)||(iTrialCountB <= 0)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	idxStartTrialA = 0;
	idxStartTrialB = 0;
	idxEndTrialA = iTrialLength - 1;
	idxEndTrialB = iTrialLength - 1;

	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		idxLeftBoundA  = min(max(piFromOffsetInTrial,idxStartTrialB + iOffset),piToOffsetInTrial);
		idxRightBoundA = max(min(piToOffsetInTrial,idxEndTrialB + iOffset),piFromOffsetInTrial);
		idxLeftBoundB  = min(max(idxLeftBoundA - iOffset,idxStartTrialB),idxEndTrialB);
		idxRightBoundB = max(min(idxRightBoundA - iOffset,idxEndTrialB),idxStartTrialB);

		if((idxLeftBoundA > idxRightBoundA)||(idxLeftBoundB > idxRightBoundB)) continue;

		//Compute the count of windows for the two trials
		iWindowCountA = (int) floorf((float)(idxRightBoundA - idxLeftBoundA + 1) / (float)iScaleWindow);
		iWindowCountB = (int) floorf((float)(idxRightBoundB - idxLeftBoundB + 1) / (float)iScaleWindow);
		iWindowCount  = min(iWindowCountA,iWindowCountB);
		
		//Compute the scaled correlation for each window
		for(j=0;j<iWindowCount;j++)
		{
			idxWindowStartA = idxLeftBoundA + j*iScaleWindow;
			idxWindowStartB = idxLeftBoundB + j*iScaleWindow;

			fCoefficientOfCorrelation = ComputePearsonCorrelation(faSignalA + idxWindowStartA,iScaleWindow,faSignalB + idxWindowStartB,iScaleWindow);

			if((fCoefficientOfCorrelation != NAN) && (fCoefficientOfCorrelation >= -1) && (fCoefficientOfCorrelation <= 1)) //consistency check
			{
				if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fCoefficientOfCorrelation;
				else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fCoefficientOfCorrelation/1.12)/(1.0-fCoefficientOfCorrelation/1.12)));
				
				iaPearsonCoefficientsDistribution[round((float)((fCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;					
				iaPearsonCoefficientsCount[i]++;
			}
		}
	}


	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaPearsonCoefficientsCount[i]);
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaPearsonCoefficientsCount[i] > 0)
			{
				faPearsonCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaPearsonCoefficientsCount[i]));
			}
			else
			{
				faPearsonCoefficientsSum[i] = NAN;
				faScaledCorrelogram[i] = NAN;
			}
		}
	}
}



//Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
/// @return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///			The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the input signals
float* CScaledCorrelationCC::GetScaledCrossCorrelogram()
{
	return faScaledCorrelogram;
}

//Returns the sum of valid Pearson coefficients of correlation for each bin of the correlogram
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PearsonCoefficientSums[i] / PearsonCoefficientCounts[i].
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
float* CScaledCorrelationCC::GetPearsonCoefficientSums()
{
	return faPearsonCoefficientsSum;
}

//Returns how many Pearson coefficients of correlation have been averaged for each bin of the scaled cross correlogram
///@return A pointer to the internal buffer containing the counts of correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of coefficients that were averaged for lag 0 is at position iCorrelationWindow.
int* CScaledCorrelationCC::GetPearsonCoefficientCounts()
{
	return iaPearsonCoefficientsCount;
}

//Returns the distribution of coefficients of correlation
///@return A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
int* CScaledCorrelationCC::GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize)
{
	piNumberOfBins = iPearsonCorrCoeffDistributionBins;
	pfBinSize = fBinSizeOfCorrCoeffDistribution;
	return iaPearsonCoefficientsDistribution;
}

//Changes the scale segment size for the scaled correlation; 
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
int CScaledCorrelationCC::ModifyScaleWindow(int piNewScale)
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
int CScaledCorrelationCC::ModifyCorrelationWindow(int piNewCorrelationWindow)
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
int CScaledCorrelationCC::ModifyTrialLength(int piNewTrialLength)
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
int CScaledCorrelationCC::ModifyAllParameters(int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
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
int CScaledCorrelationCC::GetScaleWindow()
{
	return iScaleWindow;
}

//Get the size of the correlation window
///@return Returns the size of the correlation window
int CScaledCorrelationCC::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
///@return Returns the length of the trial in original sampling units
int CScaledCorrelationCC::GetTrialLength()
{
	return iTrialLength;
}

//Get the number of bins of the distribution of coefficients of correlation
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
int CScaledCorrelationCC::GetDistributionOfCorrelationCoefficientsBinNr()
{
	return iPearsonCorrCoeffDistributionBins;
}

//Get the size of a bin of the distribution of coefficients of correlation
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
float CScaledCorrelationCC::GetDistributionOfCorrelationCoefficientsBinSize()
{
	return fBinSizeOfCorrCoeffDistribution;
}

