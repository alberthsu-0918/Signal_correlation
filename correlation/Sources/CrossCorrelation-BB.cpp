/*
	++ Source:
			CrossCorrelation-BB.cpp

	++ Description:
			A class that computes the cross-correlation for two vectors of time stamps (Binary - Binary)
			Cross-Correlation is computed in the classical way by counting

	++ History:
			16-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


#define div /

#include <stdio.h>
#include "Constants.h"
#include <math.h>
#include "MyMath.h"
#include "CrossCorrelation-BB.h"

//Construction interface -----------------------------------------------------------------------------------------------------------------------------------
		
//Constructor parameters are as follows:
///@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use the same units as the sampling of your time stamps
///@param piTrialLength			- The size of the trial in sampling units of the time stamps
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL;
CCrossCorrelationBB::CCrossCorrelationBB(int piCorrelationWindow,int piTrialLength,int& piIsError)
{
	faCrossCorrelogram   = NULL;
	piIsError			  = 0;
	if(CheckWindowsConsistency(piCorrelationWindow,piTrialLength))
	{
		iCorrelationWindow	 = piCorrelationWindow;				//The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		iTrialLength		 = piTrialLength;					//The size of the trial in original sampling time units, same as the time stamps in the input vectors
		AllocateData();
	}
	else piIsError = 1;
}

//Destructor
CCrossCorrelationBB::~CCrossCorrelationBB()
{
	CleanupData();
}

//Operational interface ------------------------------------------------------------------------------------------------------------------------------------

//Checks if the sizes of the windows are consistent
int	CCrossCorrelationBB::CheckWindowsConsistency(int piCorrW,int piTrialLength)
{
	if((piCorrW > piTrialLength) || (1.0*piCorrW*piTrialLength <= 0.0)) return 0;
	return 1;
}

//Allocates data for the cross correlogram
void CCrossCorrelationBB::AllocateData()
{
	if(faCrossCorrelogram != NULL)	CleanupData(); //do not reallocate; safety measure;
	faCrossCorrelogram = new float[2*iCorrelationWindow+1];
}

//Releases the data structures used by the cross-correlation
void CCrossCorrelationBB::CleanupData()
{
	if(faCrossCorrelogram != NULL)		//safety check
	{
		delete[] faCrossCorrelogram;
		faCrossCorrelogram = NULL;
	}
}

//Clears the bytes in the cross correlogram
void CCrossCorrelationBB::ResetCrossCorrelogram()																	
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] = 0;
	iNumberOfSummedCorrelograms = 0;
}

//Sets the correlogram to Not A Number values to signal an error or undefined correlation
void CCrossCorrelationBB::SetCorrelogramToNAN()
{
	int i;
	for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] = NAN;
	iNumberOfSummedCorrelograms = 0;
}

//The function that computes one pass cross-correlation for a pair of trials
void CCrossCorrelationBB::ComputeTrialCrossCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB)
{
	int iOffset;							//The current offset for computing the cross correlogram
	int idxA,idxB;							//Indexes to move in the two vector channels
	
	//Consistency check
	if((piIdxTrialStartA < 0) || (piIdxTrialStartB < 0) || (piIdxTrialEndA < 0) || (piIdxTrialEndB < 0)) return;

	
	idxA = piIdxTrialStartA;

	while((idxA <= piIdxTrialEndA))
	{
		idxB = piIdxTrialStartB;
		while(idxB <= piIdxTrialEndB)
		{
			iOffset = iaTimeStampsB[idxB] - iaTimeStampsA[idxA];
			if(iOffset > iCorrelationWindow) break;	//we are out of bounds so no computation is needed any more
			if(abs(iOffset) <= iCorrelationWindow)
			{
				faCrossCorrelogram[iOffset+iCorrelationWindow] = faCrossCorrelogram[iOffset+iCorrelationWindow] + 1;
			}
			idxB++;
		}
		idxA++;
		iNumberOfSummedCorrelograms++;
	}
}

//Functional interface -------------------------------------------------------------------------------------------------------------------------------------
		
//Compute the cross-correlation of two vectors of time stamps
//Parameters for ComputeCrossCorrelation
///@param	piaTimeStampsA			- The vector holding the time stamps of the first binary variable
///@param	piNrTimeStampsInA		- Number of time stamps in the first vector
///@param	piaTimeStampsB			- The vector holding the time stamps of the second binary variable
///@param	piNrTimeStampsInB		- Number of time stamps in the second vector
///@param	pbNormalizeCorrelogram	- Set to 1 to normalize the correlogram; normalization is made by dividing the correlogram to the number of time stamps in the first buffer (rate-normalized correlogram)
void CCrossCorrelationBB::ComputeCrossCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int pbNormalizeCorrelogram)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA,iCurrentTrialB;	//The number of the current trial on channel A and channel B
	int iPrevTrialA,iPrevTrialB;			//The number of the previous trials on channel A and channel B
	int isEndOfTrialA,isEndOfTrialB;		//Flags to determine whether we are at the end of a trial in channel A or in channel B

	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrTimeStampsInA <= 0) || (piNrTimeStampsInB <= 0) || (piaTimeStampsA == NULL) || (piaTimeStampsB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
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
			ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);


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
		ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);
	}
	//Normalize the cross correlogram
	if((iNumberOfSummedCorrelograms > 0)&&(pbNormalizeCorrelogram == 1)) for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] /= (float)iNumberOfSummedCorrelograms;
}


//Compute the cross-correlation of two vectors of time stamps, on a partial window of the trial. ONLY ACCEPTS ONE TRIAL!!!
//Parameters for ComputeWindowedCrossCorrelationPerTrial
///@param piaTimeStampsA		 - The vector holding the time stamps of the first binary variable
///@param piNrTimeStampsInA		 - Number of time stamps in the first vector
///@param piaTimeStampsB		 - The vector holding the time stamps of the second binary variable
///@param piNrTimeStampsInB		 - Number of time stamps in the second vector		
///@param piFromOffsetInTrial	 - The start offset in the trial where the desired window starts (time stamp)
///@param piToOffsetInTrial		 - The end offset in the trial where the desired window stops (time stamp)
///@param pbNormalizeCorrelogram - Set to 1 to normalize the correlogram
void CCrossCorrelationBB::ComputeWindowedCrossCorrelationPerTrial(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbNormalizeCorrelogram)
{
	int i;
	int idxA,idxB;						//Current position in the buffers, parsed so far
	int idxStartTrialA,idxEndTrialA;		//Index of the start of the trial and end of the trial on channel A
	int idxStartTrialB,idxEndTrialB;		//Index of the start of the trial and end of the trial on channel B
	int iCurrentTrialA,iCurrentTrialB;	//The number of the current trial on channel A and channel B
	int iPrevTrialA,iPrevTrialB;			//The number of the previous trials on channel A and channel B
	int isEndOfTrialA,isEndOfTrialB;		//Flags to determine whether we are at the end of a trial in channel A or in channel B

	//Clear up the cross correlogram
	ResetCrossCorrelogram();

	//Test for consistency
	if((faCrossCorrelogram == NULL) || (piNrTimeStampsInA <= 0) || (piNrTimeStampsInB <= 0) || (piaTimeStampsA == NULL) || (piaTimeStampsB == NULL)) 
	{
		SetCorrelogramToNAN();
		return;
	}

	//Initializing
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
	iPrevTrialA		= 0;
	iPrevTrialB		= 0;
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
				ComputeTrialCrossCorrelation(idxStartTrialA,idxEndTrialA,idxStartTrialB,idxEndTrialB);		
			}
		}
	}	

	//Normalize the cross correlogram
	if((iNumberOfSummedCorrelograms > 0)&&(pbNormalizeCorrelogram == 1)) for(i=0;i<2*iCorrelationWindow+1;i++) faCrossCorrelogram[i] /= (float)iNumberOfSummedCorrelograms;
}


//Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
///The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow.
/// @return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///			The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
float* CCrossCorrelationBB::GetCrossCorrelogram()
{
	return faCrossCorrelogram;
}

//Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter; 
///@param piNewCorrelationWindow - The new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
int CCrossCorrelationBB::ModifyCorrelationWindow(int piNewCorrelationWindow)
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
int CCrossCorrelationBB::ModifyTrialLength(int piNewTrialLength)
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
int CCrossCorrelationBB::ModifyAllParameters(int piNewCorrelationWindow,int piNewTrialLength)
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
int CCrossCorrelationBB::GetCorrelationWindow()
{
	return iCorrelationWindow;
}

//Get the trial length in original sampling units
///@return The length of the trial in sampling units
int CCrossCorrelationBB::GetTrialLength()
{
	return iTrialLength;
}