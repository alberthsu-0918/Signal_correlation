/*
	++ Source:
			ScaledCorrelation-BB.cpp

	++ Description:
			A class that computes the scaled correlation for two vectors of time stamps (Binary - Binary)
			For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer

	++ History:
			05-03-2007 - Raul C. Muresan : Created file and added main declarations

  	++Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk!
*/


#define div /


#include <stdio.h>
#include <math.h>
#include "Constants.h"
#include "MyMath.h"
#include "ScaledCorrelation-BB.h"

//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
//Constructor parameters are as follows:
///@param piScaleWindow			- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
///@param piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
///@param piTrialLength			- The size of the trial in units of your time stamps
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
CScaledCorrelationBB::CScaledCorrelationBB(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faScaledCorrelogram   = NULL;
	faFiCoefficientsSum   = NULL;
	iaFiCoefficientsCount = NULL;
	iaFiCoefficientsDistribution = NULL;
	piIsError			  = 0;
	if(CheckWindowsConsistency(piScaleWindow,piCorrelationWindow,piTrialLength))
	{
		iScaleWindow		 = piScaleWindow;					//The size of the scale window
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units, same as the time stamps in the input vectors
		bUseFisherZTransform = 0;
		fBinSizeOfCorrCoeffDistribution = (float) BIN_SIZE_OF_FI_DISTRIBUTION;					//The size of a bin in the distribution of correlation coefficients
		iFiCorrCoeffDistributionBins = (int)(2*ceil(1.0 / fBinSizeOfCorrCoeffDistribution)+1);	//Number of bins in the distribution of fi correlation coefficients; they span the interval of [-1..1]
		AllocateData();
	}
	else piIsError = 1;
}

//Destructor
CScaledCorrelationBB::~CScaledCorrelationBB()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

//Checks if the sizes of the windows are consistent
int	CScaledCorrelationBB::CheckWindowsConsistency(int piScaleW,int piCorrW,int piTrialLength)
{
	if((piScaleW > piTrialLength) || (piCorrW > piTrialLength) || (1.0*piCorrW*piScaleW*piTrialLength <= 0.0)) return 0;
	return 1;
}


//Allocates data for the cross correlogram
void CScaledCorrelationBB::AllocateData()
{
	if(faScaledCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faScaledCorrelogram = new float[2*iCorrelationWindow+1];
	faFiCoefficientsSum = new float[2*iCorrelationWindow+1];
	iaFiCoefficientsCount = new int[2*iCorrelationWindow+1];
	iaFiCoefficientsDistribution = new int[iFiCorrCoeffDistributionBins];
}

//Releases the data structures used by the cross correlation
void CScaledCorrelationBB::CleanupData()
{
	if(faScaledCorrelogram != NULL)		//safety check
	{
		delete[] faScaledCorrelogram;
		faScaledCorrelogram = NULL;
	}

	if(faFiCoefficientsSum != NULL)
	{
		delete[] faFiCoefficientsSum;
		faFiCoefficientsSum = NULL;
	}

	if(iaFiCoefficientsCount != NULL)
	{
		delete[] iaFiCoefficientsCount;
		iaFiCoefficientsCount = NULL;
	}
	if(iaFiCoefficientsDistribution != NULL)
	{
		delete[] iaFiCoefficientsDistribution;
		iaFiCoefficientsDistribution = NULL;
	}
}

//Clears the bytes in the cross correlogram
void CScaledCorrelationBB::ResetScaledCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = 0;
		faFiCoefficientsSum[i] = 0;
		iaFiCoefficientsCount[i] = 0;
	}
	for(i=0;i<iFiCorrCoeffDistributionBins;i++) iaFiCoefficientsDistribution[i] = 0;
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CScaledCorrelationBB::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		faScaledCorrelogram[i] = NAN;
		faFiCoefficientsSum[i] = NAN;
		iaFiCoefficientsCount[i] = 0;
	}
	for(i=0;i<iFiCorrCoeffDistributionBins;i++) iaFiCoefficientsDistribution[i] = 0;
}

//Determines the "Fi" (Pearson's) coefficient of correlation for a pair of segments;
float CScaledCorrelationBB::ComputeFiCorrelation(int* piaBufferA,int piBufferASize,int* piaBufferB,int piBufferBSize,int iOffsetB)
{
	int a,b,c,d;
	int idxA,idxB,iCompare;
	double dDenominator;

	if(piBufferASize*piBufferBSize == 0) return NAN; //normally this should never happen

	//compute the local a,b,c,d and Fi's ----- c = count(0,0); a = count(0,1); d = count(1,0); b = count(1,1)
	a = 0;
	b = 0;
	d = 0;

	idxA = 0;
	idxB = 0;

	while((idxA < piBufferASize) && (idxB < piBufferBSize))
	{
		iCompare = piaBufferA[idxA] - (piaBufferB[idxB] + iOffsetB);

		if(iCompare < 0)	//A < B
		{
			d++;
			idxA++;
		}
		else				
		{
			if(iCompare > 0) //A > B
			{
				a++;
				idxB++;
			}
			else             //A == B
			{
				b++;
				idxA++;
				idxB++;
			}
		}
	}

	if(idxA < piBufferASize) d += piBufferASize - idxA;
	if(idxB < piBufferBSize) a += piBufferBSize - idxB;

	c = iScaleWindow - a - b - d;

	dDenominator = ((double)(a + b)*(double)(c + d)*(double)(a + c)*(double)(b + d));

	if(dDenominator <= 0) return NAN;

	dDenominator = sqrt(dDenominator);

	//Compute the correlation coeff and ret
	return (float) (((double)b*(double)c - (double)a*(double)d) / dDenominator);
}

//The function that computes one pass scale correlation for a pair of trials - always segments relative to the beginning of signal A (reference); the autocorrelation is NOT symmetric
void CScaledCorrelationBB::ComputeTrialScaledCorrelation_Asymmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i,j;
	int iOffset;							//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;	//The indexes where the current windows start in the first and in the second channel
	int idxWindowEndA,idxWindowEndB;		//The indexes where the current windows end in the first and in the second channel
	int iTrialStartTime,iTrialEndTime;		//The times when the trial starts and when the trial ends
	int iWindowA,iWindowB;					//The number of the current window in the first and the second channel
	int iNewWindowA,iNewWindowB;			//The number of the current window in the first and the second channel in the parsing of the end of the windows
	int idxA,idxB;							//Indexes to move in the two vector channels
	float fiCoefficientOfCorrelation;		//The fi coefficient of correlation for a given window
	int iIdxTrialStartB;					//The index where trial B should start; this is adjusted to an index such that adding the offset does not exceed the trial beginning or its end
	int iIdxTrialEndB;						//The index where trial B should end; this is adjusted to an index such that adding the offset does not exceed the trial beginning or its end

	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	iTrialStartTime = piTrialIndex * iTrialLength;
	iTrialEndTime   = (piTrialIndex + 1) * iTrialLength - 1;
	
	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		//Compute the true bounds on trial B
		if((iaTimeStampsB[piIdxTrialStartB] + iOffset > iTrialEndTime) || (iaTimeStampsB[piIdxTrialEndB] + iOffset < iTrialStartTime)) continue;	//trial B is not within bounds when we add the offset
		else
		{
			j = piIdxTrialStartB;
			while(iaTimeStampsB[j] + iOffset < iTrialStartTime) j++;
			iIdxTrialStartB = j;
			j = piIdxTrialEndB;
			while(iaTimeStampsB[j] + iOffset > iTrialEndTime) j--;
			iIdxTrialEndB = j;
		}

		idxWindowStartA = piIdxTrialStartA;
		idxWindowStartB = iIdxTrialStartB;
		idxWindowEndA	= piIdxTrialStartA;
		idxWindowEndB	= iIdxTrialStartB;

		//Loop through the two trials' time stamps
		while((idxWindowStartA <= piIdxTrialEndA)&&(idxWindowStartB <= iIdxTrialEndB))
		{		
			iWindowA = (iaTimeStampsA[idxWindowStartA] - iTrialStartTime) div iScaleWindow;
			iWindowB = (iaTimeStampsB[idxWindowStartB] + iOffset - iTrialStartTime) div iScaleWindow;			

			//Bring the two scale windows to the same location in the two channels
			//First advance the first window if necessary
			while((iWindowA < iWindowB)&&(idxWindowStartA < piIdxTrialEndA)) 
			{
				idxWindowStartA++;
				iWindowA = (iaTimeStampsA[idxWindowStartA] - iTrialStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartA == piIdxTrialEndA) && (iWindowA < iWindowB)) break;

			//Next advance the second window if necessary
			while((iWindowA > iWindowB)&&(idxWindowStartB < iIdxTrialEndB)) 
			{
				idxWindowStartB++;
				iWindowB = (iaTimeStampsB[idxWindowStartB] + iOffset - iTrialStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartB == iIdxTrialEndB) && (iWindowA > iWindowB)) break;

			if(iWindowA == iWindowB) //safe check if we are aligned
			{
				//Now we have to determine the end of the scale windows
				idxA = idxWindowStartA;
				idxB = idxWindowStartB;
				iNewWindowA = iWindowA;
				iNewWindowB = iWindowB;

				//Try to advance window A and find it's end
				while((iNewWindowA == iWindowA)&&(idxA < piIdxTrialEndA))
				{
					idxA++;
					iNewWindowA = (iaTimeStampsA[idxA] - iTrialStartTime) div iScaleWindow;
				}
				if(iNewWindowA != iWindowA) idxWindowEndA = idxA - 1; //found the last index in the scale window
				else idxWindowEndA = idxA;

				//Try to advance window B and find it's end
				while((iNewWindowB == iWindowB)&&(idxB < piIdxTrialEndB))
				{
					idxB++;
					iNewWindowB = (iaTimeStampsB[idxB] + iOffset - iTrialStartTime) div iScaleWindow;
				}
				if(iNewWindowB != iWindowB) idxWindowEndB = idxB - 1; //found the last index in the scale window
				else idxWindowEndB = idxB;
			
				//We got the windows so we can compute the "Fi" correlation
				fiCoefficientOfCorrelation = ComputeFiCorrelation((iaTimeStampsA+idxWindowStartA),idxWindowEndA-idxWindowStartA+1,(iaTimeStampsB+idxWindowStartB),idxWindowEndB-idxWindowStartB+1,iOffset);

				if((fiCoefficientOfCorrelation >= -1) && (fiCoefficientOfCorrelation <= 1)) //consistency check
				{					
					if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fiCoefficientOfCorrelation;
					else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fiCoefficientOfCorrelation/1.12)/(1.0-fiCoefficientOfCorrelation/1.12)));					
					iaFiCoefficientsDistribution[round((float)((fiCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;
					iaFiCoefficientsCount[i]++;
				}

				//Now start some new windows
				idxWindowStartA = idxWindowEndA + 1;
				idxWindowStartB = idxWindowEndB + 1;
			}
		}
	}
}

//The function that computes one pass scale correlation for a pair of trials - segments the signals such that segments are aligned to signal B for right shift and signal A for left shift; the autocorrelation is symmetric
void CScaledCorrelationBB::ComputeTrialScaledCorrelation_Symmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i;
	int iOffset;							//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;	//The indexes where the current windows start in the first and in the second channel
	int idxWindowEndA,idxWindowEndB;		//The indexes where the current windows end in the first and in the second channel
	int iTrialStartTime,iTrialEndTime;		//The times when the trial starts and when the trial ends
	int iSegmentationStartTime;				//The time point that is taken as a reference for the segmentation of the common portion of the two trials after shifting
	int iWindowA,iWindowB;					//The number of the current window in the first and the second channel
	int iNewWindowA,iNewWindowB;			//The number of the current window in the first and the second channel in the parsing of the end of the windows
	int idxA,idxB;							//Indexes to move in the two vector channels
	float fiCoefficientOfCorrelation;		//The fi coefficient of correlation for a given window

	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	iTrialStartTime = piTrialIndex * iTrialLength;
	iTrialEndTime   = (piTrialIndex + 1) * iTrialLength - 1;
	
	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
 		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials
		if(iOffset > 0) iSegmentationStartTime = iTrialStartTime + iOffset;	//This setting is necessary in order to have a symmetric autocorrelation
		else iSegmentationStartTime = iTrialStartTime;

		idxWindowStartA = piIdxTrialStartA;
		idxWindowStartB = piIdxTrialStartB;
		idxWindowEndA	= piIdxTrialStartA;
		idxWindowEndB	= piIdxTrialStartB;

		//Loop through the two trials' time stamps
		while((idxWindowStartA <= piIdxTrialEndA)&&(idxWindowStartB <= piIdxTrialEndB))
		{		
			iWindowA = (iaTimeStampsA[idxWindowStartA] - iSegmentationStartTime) div iScaleWindow;
			iWindowB = (iaTimeStampsB[idxWindowStartB] + iOffset - iSegmentationStartTime) div iScaleWindow;			

			//Bring the two scale windows to the same location in the two channels
			//First advance the first window if necessary
			while((iWindowA < iWindowB)&&(idxWindowStartA < piIdxTrialEndA)) 
			{
				idxWindowStartA++;
				iWindowA = (iaTimeStampsA[idxWindowStartA] - iSegmentationStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartA == piIdxTrialEndA) && (iWindowA < iWindowB)) break;

			//Next advance the second window if necessary
			while((iWindowA > iWindowB)&&(idxWindowStartB < piIdxTrialEndB)) 
			{
				idxWindowStartB++;
				iWindowB = (iaTimeStampsB[idxWindowStartB] + iOffset - iSegmentationStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartB == piIdxTrialEndB) && (iWindowA > iWindowB)) break;

			if(iWindowA == iWindowB) //safe check if we are aligned
			{
				//Now we have to determine the end of the scale windows
				idxA = idxWindowStartA;
				idxB = idxWindowStartB;
				iNewWindowA = iWindowA;
				iNewWindowB = iWindowB;

				//Try to advance window A and find it's end
				while((iNewWindowA == iWindowA)&&(idxA < piIdxTrialEndA))
				{
					idxA++;
					iNewWindowA = (iaTimeStampsA[idxA] - iSegmentationStartTime) div iScaleWindow;
				}
				if(iNewWindowA != iWindowA) idxWindowEndA = idxA - 1; //found the last index in the scale window
				else idxWindowEndA = idxA;

				//Try to advance window B and find it's end
				while((iNewWindowB == iWindowB)&&(idxB < piIdxTrialEndB))
				{
					idxB++;
					iNewWindowB = (iaTimeStampsB[idxB] + iOffset - iSegmentationStartTime) div iScaleWindow;
				}
				if(iNewWindowB != iWindowB) idxWindowEndB = idxB - 1; //found the last index in the scale window
				else idxWindowEndB = idxB;
			
				//We got the windows so we can compute the "Fi" correlation
				fiCoefficientOfCorrelation = ComputeFiCorrelation((iaTimeStampsA+idxWindowStartA),idxWindowEndA-idxWindowStartA+1,(iaTimeStampsB+idxWindowStartB),idxWindowEndB-idxWindowStartB+1,iOffset);

				if((fiCoefficientOfCorrelation >= -1) && (fiCoefficientOfCorrelation <= 1)) //consistency check
				{					
					if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fiCoefficientOfCorrelation;
					else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fiCoefficientOfCorrelation/1.12)/(1.0-fiCoefficientOfCorrelation/1.12)));					
					iaFiCoefficientsDistribution[round((float)((fiCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;
					iaFiCoefficientsCount[i]++;
				}

				//Now start some new windows
				idxWindowStartA = idxWindowEndA + 1;
				idxWindowStartB = idxWindowEndB + 1;
			}
		}
	}	
}


//Same as ComputeTrialScaledCorrelation_Symmetric but faster version with local expansion and use of a memory buffer
void CScaledCorrelationBB::ComputeTrialScaledCorrelation_Symmetric_Fast(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex)
{
	int i,j;
	int iOffset;							//The current offset for computing the cross correlogram
	int idxWindowStartA,idxWindowStartB;	//The indexes where the current windows start in the first and in the second channel
	int idxWindowEndA,idxWindowEndB;		//The indexes where the current windows end in the first and in the second channel
	int iTrialStartTime,iTrialEndTime;		//The times when the trial starts and when the trial ends
	int iSegmentationStartTime;				//The time point that is taken as a reference for the segmentation of the common portion of the two trials after shifting
	int iWindowA,iWindowB;					//The number of the current window in the first and the second channel
	int iNewWindowA,iNewWindowB;			//The number of the current window in the first and the second channel in the parsing of the end of the windows
	int idxA,idxB;							//Indexes to move in the two vector channels
	int iSpikesInA,iSpikesInB;				//The number of spikes in buffers A and B
	float fiCoefficientOfCorrelation;		//The fi coefficient of correlation for a given window
	int	*iaSpikeBufferA,*iaSpikeBufferB;	//Temporary buffer for spikes in A and B (B is allocated, A is referenced only)
	int iIdxTrialEndA,iIdxTrialEndB;		//Indexes of last spikes in buffers A and B of spikes
	int a,b,c,d;							//Variables used to compute the fi correlation coefficient
	int idxWA,idxWB,iCompare;				//Indexing and auxiliary variables
	double dDenominator;					//Denominator used for fi correlation formula
	int *iaSpikesWindowA,*iaSpikesWindowB;	//Spikes in the matching windows
	int iSpikesInWindowA,iSpikesInWindowB;	//Count of spikes in the matched windows

	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	iTrialStartTime = 0;
	iTrialEndTime   = iTrialLength - 1;
	
	iSpikesInA = piIdxTrialEndA-piIdxTrialStartA+1;
	iSpikesInB = piIdxTrialEndB-piIdxTrialStartB+1;
	
	//Consistency check
	if((iSpikesInA <= 0) || (iSpikesInB <= 0)) return;
	
	iIdxTrialEndA = piIdxTrialEndA-piIdxTrialStartA;		
	iaSpikeBufferA = &iaTimeStampsA[piIdxTrialStartA];		
	iIdxTrialEndB = piIdxTrialEndB-piIdxTrialStartB;
	iaSpikeBufferB = new int[iSpikesInB];
	
	//For each offset
	for(i=0;i<2*iCorrelationWindow+1;i++)
	{
 		iOffset	= iCorrelationWindow - i;	//The offset for comparing the two trials

		//Create an already offsetted buffer for B
		for(j=0;j<iSpikesInB;j++) iaSpikeBufferB[j] = iaTimeStampsB[piIdxTrialStartB+j] + iOffset;
		
		if(iOffset > 0) iSegmentationStartTime = iOffset;	//This setting is necessary in order to have a symmetric autocorrelation
		else iSegmentationStartTime = 0;

		idxWindowStartA = 0;
		idxWindowStartB = 0;
		idxWindowEndA	= iIdxTrialEndA;
		idxWindowEndB	= iIdxTrialEndB;

		//Loop through the two trials' time stamps
		while((idxWindowStartA <= iIdxTrialEndA)&&(idxWindowStartB <= iIdxTrialEndB))
		{		
			iWindowA = (iaSpikeBufferA[idxWindowStartA] - iSegmentationStartTime) div iScaleWindow;
			iWindowB = (iaSpikeBufferB[idxWindowStartB] - iSegmentationStartTime) div iScaleWindow;			

			//Bring the two scale windows to the same location in the two channels
			//First advance the first window if necessary
			while((iWindowA < iWindowB)&&(idxWindowStartA < iIdxTrialEndA)) 
			{
				idxWindowStartA++;
				iWindowA = (iaSpikeBufferA[idxWindowStartA] - iSegmentationStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartA == iIdxTrialEndA) && (iWindowA < iWindowB)) break;

			//Next advance the second window if necessary
			while((iWindowA > iWindowB)&&(idxWindowStartB < iIdxTrialEndB)) 
			{
				idxWindowStartB++;
				iWindowB = (iaSpikeBufferB[idxWindowStartB] - iSegmentationStartTime) div iScaleWindow;
			}
			//If we reached the last index it means there are no more matching scale windows
			if((idxWindowStartB == iIdxTrialEndB) && (iWindowA > iWindowB)) break;

			if(iWindowA == iWindowB) //safe check if we are aligned
			{
				//Now we have to determine the end of the scale windows
				idxA = idxWindowStartA;
				idxB = idxWindowStartB;
				iNewWindowA = iWindowA;
				iNewWindowB = iWindowB;

				//Try to advance window A and find it's end
				while((iNewWindowA == iWindowA)&&(idxA < iIdxTrialEndA))
				{
					idxA++;
					iNewWindowA = (iaSpikeBufferA[idxA] - iSegmentationStartTime) div iScaleWindow;
				}
				if(iNewWindowA != iWindowA) idxWindowEndA = idxA - 1; //found the last index in the scale window
				else idxWindowEndA = idxA;

				//Try to advance window B and find it's end
				while((iNewWindowB == iWindowB)&&(idxB < iIdxTrialEndB))
				{
					idxB++;
					iNewWindowB = (iaSpikeBufferB[idxB] - iSegmentationStartTime) div iScaleWindow;
				}
				if(iNewWindowB != iWindowB) idxWindowEndB = idxB - 1; //found the last index in the scale window
				else idxWindowEndB = idxB;
	
				//We got the windows so we can compute the "Fi" correlation	
				iaSpikesWindowA = &iaSpikeBufferA[idxWindowStartA];
				iaSpikesWindowB = &iaSpikeBufferB[idxWindowStartB];
				iSpikesInWindowA = idxWindowEndA - idxWindowStartA + 1;
				iSpikesInWindowB = idxWindowEndB - idxWindowStartB + 1;

				//compute the local a,b,c,d and Fi's ----- c = count(0,0); a = count(0,1); d = count(1,0); b = count(1,1)
				a = 0;
				b = 0;
				d = 0;

				idxWA = 0;
				idxWB = 0;

				while((idxWA < iSpikesInWindowA) && (idxWB < iSpikesInWindowB))
				{
					iCompare = iaSpikesWindowA[idxWA] - iaSpikesWindowB[idxWB];

					if(iCompare < 0)	//A < B
					{
						d++;
						idxWA++;
					}
					else				//A > B
					{
						if(iCompare > 0)
						{
							a++;
							idxWB++;
						}
						else            //A == B
						{
							b++;
							idxWA++;
							idxWB++;
						}
					}					
				}

				if(idxWA < iSpikesInWindowA) d += iSpikesInWindowA - idxWA;
				if(idxWB < iSpikesInWindowB) a += iSpikesInWindowB - idxWB;

				c = iScaleWindow - a - b - d;

				dDenominator = ((double)(a + b)*(double)(c + d)*(double)(a + c)*(double)(b + d));

				if(dDenominator > 0)
				{
					dDenominator = sqrt(dDenominator);

					//Compute the correlation coeff and ret
					fiCoefficientOfCorrelation = (float) (((double)b*(double)c - (double)a*(double)d) / dDenominator);

					if((fiCoefficientOfCorrelation >= -1) && (fiCoefficientOfCorrelation <= 1)) //consistency check
					{					
						if(bUseFisherZTransform != 1) faScaledCorrelogram[i] += fiCoefficientOfCorrelation;
						else faScaledCorrelogram[i] = faScaledCorrelogram[i] + (float)(0.5*log((1.0+fiCoefficientOfCorrelation/1.12)/(1.0-fiCoefficientOfCorrelation/1.12)));					
						iaFiCoefficientsDistribution[round((float)((fiCoefficientOfCorrelation + 1.0) / fBinSizeOfCorrCoeffDistribution))]++;
						iaFiCoefficientsCount[i]++;
					}
				}

				//Now start some new windows
				idxWindowStartA = idxWindowEndA+1;
				idxWindowStartB = idxWindowEndB+1;
			}
		}
	}

	delete[] iaSpikeBufferB;
}


//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		
//Compute the scaled correlation of two vectors of time stamps
//Parameters for ComputeScaledCorrelation
///@param piaTimeStampsA		- The vector holding the time stamps of the first neuron
///@param piNrTimeStampsInA		- Number of time stamps in the first vector
///@param piaTimeStampsB		- The vector holding the time stamps of the second neuron
///@param piNrTimeStampsInB		- Number of time stamps in the second vector
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
void CScaledCorrelationBB::ComputeScaledCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int pbUseFisherZTransform)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA,iCurrentTrialB;	//The number of the current trial on channel A and channel B
	int iPrevTrialA,iPrevTrialB;			//The number of the previous trials on channel A and channel B
	int isEndOfTrialA,isEndOfTrialB;		//Flags to determine whether we are at the end of a trial in channel A or in channel B

	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrTimeStampsInA <= 0) || (piNrTimeStampsInB <= 0) || (piaTimeStampsA == NULL) || (piaTimeStampsB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
	iaTimeStampsA	= piaTimeStampsA;
	iaTimeStampsB	= piaTimeStampsB;
	idxA			= 0;
	idxB			= 0;
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= 0;
	iCurrentTrialA	= 0;
	iCurrentTrialB	= 0;
	iPrevTrialA		= 0;
	iPrevTrialB		= 0;
	isEndOfTrialA	= 0;
	isEndOfTrialB	= 0;
	
	//Main loop that crosses over the data
	while((idxA < piNrTimeStampsInA) && (idxB < piNrTimeStampsInB))
	{
		//Compute the current trial number
		iCurrentTrialA = piaTimeStampsA[idxA] div iTrialLength;
		iCurrentTrialB = piaTimeStampsB[idxB] div iTrialLength;

		//If the trial has just changed then we have parsed a whole trial so we can compute its bounds
		if(iCurrentTrialA != iPrevTrialA) 
		{
			idxEndTrialA  = idxA - 1;
			isEndOfTrialA = 1;
		}
		if(iCurrentTrialB != iPrevTrialB) 
		{
			idxEndTrialB  = idxB - 1;
			isEndOfTrialB = 1;
		}

		//If we have both trials detected then we have to compute the "fi" correlation coefficient scores
		if((isEndOfTrialA)&&(isEndOfTrialB)&&(iPrevTrialA == iPrevTrialB))
		{
			//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
			//ComputeTrialScaledCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iPrevTrialA);
			ComputeTrialScaledCorrelation_Symmetric_Fast(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iPrevTrialA);

			//Last operation is to start a new trial
			idxStartTrialA = idxA;
			idxStartTrialB = idxB;
			isEndOfTrialA  = 0;
			isEndOfTrialB  = 0;
			iPrevTrialA	   = iCurrentTrialA;
			iPrevTrialB	   = iCurrentTrialB;
		}


		//Advance the indexes, depending on which has to be moved to get the trials synchronously
		if(iCurrentTrialA < iCurrentTrialB)
		{
			if(isEndOfTrialA)				//If we were at the and of the trial, now after the increment we are no more!
			{	
				iPrevTrialA = iCurrentTrialA;
				isEndOfTrialA = 0;
				idxStartTrialA = idxA;
			}
			idxA++;							//Try to catch up with the second channel
		}
		else
		{
			if(iCurrentTrialA > iCurrentTrialB) 
			{
				if(isEndOfTrialB)			//If we were at the and of the trial, now after the increment we are no more!
				{	
					iPrevTrialB = iCurrentTrialB;
					isEndOfTrialB = 0;
					idxStartTrialB = idxB;
				}
				idxB++;						//Try to catch up with the first channel
			}
			else
			{
				if(isEndOfTrialA)				//If we were at the and of the trial, now after the increment we are no more!
				{	
					iPrevTrialA = iCurrentTrialA;
					isEndOfTrialA = 0;
					idxStartTrialA = idxA;
				}
				if(isEndOfTrialB)				//If we were at the and of the trial, now after the increment we are no more!
				{		
					iPrevTrialB = iCurrentTrialB;
					isEndOfTrialB = 0;
					idxStartTrialB = idxB;
				}

				idxA++;						//Simultaneous increment
				idxB++;
			}
		}
	}


	//Process the last trial, but only if they match... because we exited the while from the last branch "Simultaneous increment"
	if(iCurrentTrialA == iCurrentTrialB)
	{
		//Compute the last indexes for channels that did not end but where the index is whithin the current trial; if the while exited because of a buffer end then we anyway can compute the index of the end of the trials as idxA - 1 or idxB - 1
		while((iCurrentTrialA == iCurrentTrialB) && (idxA < piNrTimeStampsInA))
		{
			idxA++;
			if(idxA < piNrTimeStampsInA) iCurrentTrialA = piaTimeStampsA[idxA] div iTrialLength;	
		}
		idxEndTrialA = idxA - 1;
		while((iCurrentTrialA == iCurrentTrialB) && (idxB < piNrTimeStampsInB))
		{
			idxB++;
			if(idxB < piNrTimeStampsInB) iCurrentTrialB = piaTimeStampsB[idxB] div iTrialLength;	
		}
		idxEndTrialB = idxB - 1;

		//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
		//ComputeTrialScaledCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iCurrentTrialA);
		ComputeTrialScaledCorrelation_Symmetric_Fast(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iCurrentTrialA);
	}

	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaFiCoefficientsCount[i] > 0)
			{
				faFiCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaFiCoefficientsCount[i]);
			}
			else
			{
				faScaledCorrelogram[i] = NAN;
				faFiCoefficientsSum[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaFiCoefficientsCount[i] > 0)
			{
				faFiCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaFiCoefficientsCount[i]));
			}
			else 
			{
				faScaledCorrelogram[i] = NAN;
				faFiCoefficientsSum[i] = NAN;
			}
		}
	}
}


///Compute the scaled correlation of two vectors of time stamps, on a partial window of the trial; 
///@warning ONLY ACCEPTS ONE TRIAL!!! so do not pass signals longer than the size of a single trial! Subsequent trials will be discarded.
//Parameters for ComputeWindowedScaledCorrelationPerTrial
///@param piaTimeStampsA		- The vector holding the time stamps of the first neuron
///@param piNrTimeStampsInA		- Number of time stamps in the first vector
///@param piaTimeStampsB		- The vector holding the time stamps of the second neuron
///@param piNrTimeStampsInB		- Number of time stamps in the second vector
///@param piFromOffsetInTrial	- The start offset in the trial where the desired window starts (time stamp)
///@param piToOffsetInTrial		- The end offset in the trial where the desired window stops (time stamp)
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
void CScaledCorrelationBB::ComputeWindowedScaledCorrelationPerTrial(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA,iCurrentTrialB;	//The number of the current trial on channel A and channel B
	int isEndOfTrialA,isEndOfTrialB;		//Flags to determine whether we are at the end of a trial in channel A or in channel B

	//Clear up the cross correlogram
	ResetScaledCorrelogram();

	//Test for consistency
	if((faScaledCorrelogram == NULL) || (piNrTimeStampsInA <= 0) || (piNrTimeStampsInB <= 0) || (piaTimeStampsA == NULL) || (piaTimeStampsB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
	bUseFisherZTransform = pbUseFisherZTransform;
	iaTimeStampsA	= piaTimeStampsA;
	iaTimeStampsB	= piaTimeStampsB;
	idxA			= 0;
	idxB			= 0;
	idxStartTrialA	= 0;
	idxEndTrialA	= 0;
	idxStartTrialB	= 0;
	idxEndTrialB	= piNrTimeStampsInB - 1; //take the whole trial B
	iCurrentTrialA	= 0;
	iCurrentTrialB	= 0;
	isEndOfTrialA	= 0;
	isEndOfTrialB	= 0;

	//Compute the current trial number
	iCurrentTrialA = piaTimeStampsA[idxA] div iTrialLength;
	iCurrentTrialB = piaTimeStampsB[idxB] div iTrialLength;

	//Check if the trials actually match
	if(iCurrentTrialA != iCurrentTrialB) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Find the region of interest in the first vector
	while((idxA < piNrTimeStampsInA)&&(piaTimeStampsA[idxA] < piFromOffsetInTrial)) idxA++;
	if(idxA < piNrTimeStampsInA)
	{
		idxStartTrialA = idxA;
		
		while((idxA < piNrTimeStampsInA)&&(piaTimeStampsA[idxA] <= piToOffsetInTrial)) idxA++;
		idxEndTrialA = idxA - 1;
		
		if(piaTimeStampsA[idxEndTrialA] <= piToOffsetInTrial)
		{
			if(idxStartTrialA <= idxEndTrialA) //Valid window found verification
			{
				//Now we have both trials isolated between idxTrialStartA,B and idxTrialEndA,B respective
				ComputeTrialScaledCorrelation_Symmetric_Fast(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB,iCurrentTrialA);		
			}
		}
	}	

	//Compute the average of the scaled correlogram
	if(bUseFisherZTransform != 1)
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaFiCoefficientsCount[i] > 0)
			{
				faFiCoefficientsSum[i] = faScaledCorrelogram[i];
				faScaledCorrelogram[i] = faScaledCorrelogram[i] / (float)(iaFiCoefficientsCount[i]);
			}
			else
			{
				faScaledCorrelogram[i] = NAN;
				faFiCoefficientsSum[i] = NAN;
			}
		}
	}
	else
	{
		for(i=0;i<2*iCorrelationWindow+1;i++)
		{
			if(iaFiCoefficientsCount[i] > 0)
			{
				faFiCoefficientsSum[i] = (float) (1.12*tanh(faScaledCorrelogram[i]));
				faScaledCorrelogram[i] = (float) (1.12*tanh(faScaledCorrelogram[i] / (float) iaFiCoefficientsCount[i]));
			}
			else 
			{
				faScaledCorrelogram[i] = NAN;
				faFiCoefficientsSum[i] = NAN;
			}
		}
	}
}


//Returns the buffer with the computed scaled-cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
/// @return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///			The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the input vectors
float* CScaledCorrelationBB::GetScaledCrossCorrelogram()
{
	return faScaledCorrelogram;
}

//Returns the sum of valid Fi coefficients of correlation for each bin of the correlogram
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = FiCoefficientSums[i] / FiCoefficientCounts[i].
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
float* CScaledCorrelationBB::GetFiCoefficientSums()
{
	return faFiCoefficientsSum;
}

//Returns how many Fi coefficients of correlation have been averaged for each bin of the scaled cross correlogram
///@return A pointer to the internal buffer containing the counts of "Fi" correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of valid Fi coefficients, that were averaged to compute the scaled cross correlogram at lag 0, is at position iCorrelationWindow
int* CScaledCorrelationBB::GetFiCoefficientCounts()
{
	return iaFiCoefficientsCount;
}

//Returns the distribution of coefficients of correlation
///@return A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
int* CScaledCorrelationBB::GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize)
{
	piNumberOfBins = iFiCorrCoeffDistributionBins;
	pfBinSize = fBinSizeOfCorrCoeffDistribution;
	return iaFiCoefficientsDistribution;
}

//Changes the scale window for the scaled correlation; 
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
int CScaledCorrelationBB::ModifyScaleWindow(int piNewScale)
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
int CScaledCorrelationBB::ModifyCorrelationWindow(int piNewCorrelationWindow)
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
int CScaledCorrelationBB::ModifyTrialLength(int piNewTrialLength)
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
int CScaledCorrelationBB::ModifyAllParameters(int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
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
int CScaledCorrelationBB::GetScaleWindow()
{
	return iScaleWindow;
}

//Get the size of the correlation window
///@return Returns the size of the correlation window
int CScaledCorrelationBB::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
///@return Returns the length of the trial in original sampling units
int CScaledCorrelationBB::GetTrialLength()
{
	return iTrialLength;
}

//Get the number of bins of the distribution of coefficients of correlation
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
int CScaledCorrelationBB::GetDistributionOfCorrelationCoefficientsBinNr()
{
	return iFiCorrCoeffDistributionBins;
}

//Get the size of a bin of the distribution of coefficients of correlation
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
float CScaledCorrelationBB::GetDistributionOfCorrelationCoefficientsBinSize()
{
	return fBinSizeOfCorrCoeffDistribution;
}

