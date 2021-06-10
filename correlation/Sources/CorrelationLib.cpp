// CORRELATIONLIB.cpp : Defines the entry point for the DLL application.
//

#include "StdAfx.h"
#include "CorrelationLib.h"
#include "ScaledCorrelation-BB.h"
#include "ScaledCorrelation-CC.h"
#include "ScaledCorrelation-BC.h"
#include "CrossCorrelation-BB.h"
#include "CrossCorrelation-CC.h"
#include "CrossCorrelation-BC.h"

///Entry into the Dll space
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}


#define ERR_NULL_OBJECT -5


//Exports go here ----------------------------------------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//Interface for Binary - Binary scaled corr ======================================================================================================================================

///Creates an object that computes scaled correlation on a pair of binary - binary signals
// Constructor parameters are as follows:
///@param piScaleWindow			- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
///@param piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use the same units as the sampling of your time stamps
///@param piTrialLength			- The size of the trial in units of your time stamps
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL
///@return An untyped pointer to the created object
CORRELATIONLIB_API void* CreateScaledCorrelationComputerBB(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError)
{
	return new CScaledCorrelationBB(piScaleWindow,piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for scaled correlation objects on binary-binary pairs; always call this to cleanup your data space
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
CORRELATIONLIB_API void FreeScaledCorrelationComputerBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CScaledCorrelationBB*) (pvObject);
	}
}

///Compute the scaled correlation of two vectors of time stamps
//Parameters for ComputeScaledCorrelation
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piaTimeStampsA		- The vector holding the time stamps of the first binary variable
///@param piNrTimeStampsInA		- Number of time stamps in the first vector
///@param piaTimeStampsB		- The vector holding the time stamps of the second binary variable
///@param piNrTimeStampsInB		- Number of time stamps in the second vector
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
CORRELATIONLIB_API void ComputeScaledCorrelationBB(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationBB*) pvObject)->ComputeScaledCorrelation(piaTimeStampsA,piNrTimeStampsInA,piaTimeStampsB,piNrTimeStampsInB,pbUseFisherZTransform);
	}
}		

///Compute the scaled correlation of two vectors of time stamps, on a partial window of the trial; ONLY ACCEPTS ONE TRIAL!!!
//
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piaTimeStampsA		- The vector holding the time stamps of the first binary variable
///@param piNrTimeStampsInA		- Number of time stamps in the first vector
///@param piaTimeStampsB		- The vector holding the time stamps of the second binary variable
///@param piNrTimeStampsInB		- Number of time stamps in the second vector
///@param piFromOffsetInTrial	- The start offset in the trial where the desired window starts (time stamp)
///@param piToOffsetInTrial		- The end offset in the trial where the desired window stops (time stamp)
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation
CORRELATIONLIB_API void ComputeWindowedScaledCorrelationPerTrialBB(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationBB*) pvObject)->ComputeWindowedScaledCorrelationPerTrial(piaTimeStampsA,piNrTimeStampsInA,piaTimeStampsB,piNrTimeStampsInB,piFromOffsetInTrial,piToOffsetInTrial,pbUseFisherZTransform);
	}
}


///Returns the buffer with the computed scaled-cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///		   The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the input vectors
CORRELATIONLIB_API float* GetScaledCrossCorrelogramBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetScaledCrossCorrelogram();
	}
	else return NULL;
}


///Returns the sum of valid Fi coefficients of correlation for each bin of the correlogram
//
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = FiCoefficientSums[i] / FiCoefficientCounts[i].
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
CORRELATIONLIB_API float* GetFiCoefficientSumsBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetFiCoefficientSums();
	}
	else return NULL;
}	


///Returns how many Fi coefficients of correlation have been averaged for each bin of the scaled cross correlogram
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return A pointer to the internal buffer containing the counts of "Fi" correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of valid Fi coefficients, that were averaged to compute the scaled cross correlogram at lag 0, is at position iCorrelationWindow
CORRELATIONLIB_API int* GetFiCoefficientCountsBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetFiCoefficientCounts();
	}
	else return NULL;
}

///Returns the distribution of coefficients of correlation
//
///@param  pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piNumberOfBins - The number of bins in the distribution (output parameter)
///@param pfBinSize - The size of one bin of the distribution (output parameter) 
///@return - A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return - The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return - The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
CORRELATIONLIB_API int* GetDistributionOfCorrelationCoefficientsBB(void* pvObject,int& piNumberOfBins,float& pfBinSize)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetDistributionOfCorrelationCoefficients(piNumberOfBins,pfBinSize);
	}
	else return NULL;
}

///Changes the scale window for the scaled correlation;
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
CORRELATIONLIB_API int ScaledCorrelationModifyScaleWindowBB(void* pvObject,int piNewScale)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->ModifyScaleWindow(piNewScale);
	}
	else return ERR_NULL_OBJECT;
}

///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piNewCorrelationWindow - The new size of the correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial
CORRELATIONLIB_API int ScaledCorrelationModifyCorrelationWindowBB(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}

///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the scale window; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyTrialLengthBB(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}

///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@param piNewScale - The new size of the scale segment
///@param piNewCorrelationWindow - The new size of the correlation window
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success, -1 if the new scale window is too small; -2 if the new correlation window is too small; -3 if the new scale window is larger than the trial size; -4 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyAllParametersBB(void* pvObject,int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->ModifyAllParameters(piNewScale,piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}

///Get the size of the current scale window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns the size of the scale segment window
CORRELATIONLIB_API int ScaledCorrelationGetScaleWindowBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetScaleWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns the size of the correlation window
CORRELATIONLIB_API int ScaledCorrelationGetCorrelationWindowBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns the length of the trial in original sampling units
CORRELATIONLIB_API int ScaledCorrelationGetTrialLengthBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}

///Get the number of bins of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficientsBB)
CORRELATIONLIB_API int GetDistributionOfCorrelationCoefficientsBinNrBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetDistributionOfCorrelationCoefficientsBinNr();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of a bin of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBB function
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficientsBB)
CORRELATIONLIB_API float GetDistributionOfCorrelationCoefficientsBinSizeBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBB*) pvObject)->GetDistributionOfCorrelationCoefficientsBinSize();
	}
	else return ERR_NULL_OBJECT;
}

//====================================================================================================================================================================================





//Interface for Continuous - Continuous scaled corr ==================================================================================================================================
///Creates an object that computes scaled correlation on a pair of continuous - continuous signals
//Constructor parameters are as follows:
///@param piScaleWindow		- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your signals
///@param piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use again the same units as the sampling of your signals
///@param piTrialLength		- The size of the trial in sampling units of your signals
///@param piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL
///@return An untyped pointer to the created object
CORRELATIONLIB_API void* CreateScaledCorrelationComputerCC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError)
{
	return new CScaledCorrelationCC(piScaleWindow,piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for scaled correlation objects on continuous-continuous input pairs; always call this to cleanup your data space
//
///@param pvObject - pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
CORRELATIONLIB_API void FreeScaledCorrelationComputerCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CScaledCorrelationCC*) (pvObject);
	}
}


///Compute the scaled correlation of two digitized continuous signals
//Parameters for ComputeScaledCorrelation
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param pfaSamplesA			- The vector holding the first signal's input samples
///@param piNrSamplesInA		- Number of samples in the first vector
///@param pfaSamplesB			- The vector holding the second signal's input samples
///@param piNrSamplesInB		- Number of samples in the second vector
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
CORRELATIONLIB_API void ComputeScaledCorrelationCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationCC*) pvObject)->ComputeScaledCorrelation(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,pbUseFisherZTransform);
	}
}	

///Compute the scaled correlation of two digitized continuous signals without checking for borders at the end of the trials
//WARNING: for each trial, the buffer must contain before and after the trial a number of samples equal to the correlation window!	
//Parameters for ComputeScaledCorrelationNoBorderChecking
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param pfaSamplesA			- The vector holding the first signal's input samples
///@param piNrSamplesInA		- Number of samples in the first vector
///@param pfaSamplesB			- The vector holding the second signal's input samples
///@param piNrSamplesInB		- Number of samples in the second vector
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
CORRELATIONLIB_API void ComputeScaledCorrelationNoBorderCheckingCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationCC*) pvObject)->ComputeScaledCorrelationNoBorderChecking(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,pbUseFisherZTransform);
	}
}	

///Compute the scaled correlation of two digitized continuous signals, on a partial window of the trial; ONLY ACCEPTS ONE TRIAL!!!
//Parameters for ComputeWindowedScaledCorrelationPerTrial
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param pfaSamplesA			- The vector holding the first signal's input samples
///@param piNrSamplesInA		- Number of samples in the first vector
///@param pfaSamplesB			- The vector holding the second signal's input samples
///@param piNrSamplesInB		- Number of samples in the second vector
///@param piFromOffsetInTrial	- The start offset in the trial where the desired window starts
///@param piToOffsetInTrial		- The end offset in the trial where the desired window stops
///@param pbUseFisherZTransform	- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
CORRELATIONLIB_API void ComputeWindowedScaledCorrelationPerTrialCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationCC*) pvObject)->ComputeWindowedScaledCorrelationPerTrial(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,piFromOffsetInTrial,piToOffsetInTrial,pbUseFisherZTransform);
	}
}


///Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
///@param pvObject				- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///			The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the input signals
CORRELATIONLIB_API float* GetScaledCrossCorrelogramCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetScaledCrossCorrelogram();
	}
	else return NULL;
}

///Returns the sum of valid Pearson coefficients of correlation for each bin of the correlogram
//
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PearsonCoefficientSums[i] / PearsonCoefficientCounts[i].
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
CORRELATIONLIB_API float* GetPearsonCoefficientSumsCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetPearsonCoefficientSums();
	}
	else return NULL;
}

//Returns how many Pearson coefficients of correlation have been averaged for each bin of the scaled cross correlogram
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return A pointer to the internal buffer containing the counts of correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of coefficients that were averaged for lag 0 is at position iCorrelationWindow.
CORRELATIONLIB_API int* GetPearsonCoefficientCountsCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetPearsonCoefficientCounts();
	}
	else return NULL;
}

///Returns the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param piNumberOfBins - The number of bins in the distribution (output parameter)
///@param pfBinSize - The size of one bin of the distribution (output parameter) 
///@return A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
CORRELATIONLIB_API int* GetDistributionOfCorrelationCoefficientsCC(void* pvObject,int& piNumberOfBins,float& pfBinSize)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetDistributionOfCorrelationCoefficients(piNumberOfBins,pfBinSize);
	}
	else return NULL;
}

///Changes the scale segment size for the scaled correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
CORRELATIONLIB_API int ScaledCorrelationModifyScaleWindowCC(void* pvObject,int piNewScale)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->ModifyScaleWindow(piNewScale);
	}
	else return ERR_NULL_OBJECT;
}

///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param piNewCorrelationWindow - The new size of the correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial
CORRELATIONLIB_API int ScaledCorrelationModifyCorrelationWindowCC(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}

///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the scale window; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyTrialLengthCC(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}

///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@param piNewScale - The new size of the scale segment
///@param piNewCorrelationWindow - The new size of the correlation window
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success, -1 if the new scale window is too small; -2 if the new correlation window is too small; -3 if the new scale window is larger than the trial size; -4 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyAllParametersCC(void* pvObject,int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->ModifyAllParameters(piNewScale,piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}

///Get the size of the current scale window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns the size of the scale segment window
CORRELATIONLIB_API int ScaledCorrelationGetScaleWindowCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetScaleWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns the size of the correlation window
CORRELATIONLIB_API int ScaledCorrelationGetCorrelationWindowCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns the length of the trial in original sampling units
CORRELATIONLIB_API int ScaledCorrelationGetTrialLengthCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}

///Get the number of bins of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficientsCC)
CORRELATIONLIB_API int GetDistributionOfCorrelationCoefficientsBinNrCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetDistributionOfCorrelationCoefficientsBinNr();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of a bin of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerCC function
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficientsCC)
CORRELATIONLIB_API float GetDistributionOfCorrelationCoefficientsBinSizeCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationCC*) pvObject)->GetDistributionOfCorrelationCoefficientsBinSize();
	}
	else return ERR_NULL_OBJECT;
}

//====================================================================================================================================================================================





//Interface for Binary - Continuous scaled corr ====================================================================================================================================

//Construction interface

//Constructor parameters are as follows:
///Creates an object that computes scaled correlation on a pair of binary - continuous signals
//Constructor parameters are as follows:
///@param piScaleWindow			- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
///@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
///@param piTrialLength			- The size of the trial in original sampling units
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL
///@warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
///@return An untyped pointer to the created object!
CORRELATIONLIB_API void* CreateScaledCorrelationComputerBC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError)
{
	return new CScaledCorrelationBC(piScaleWindow,piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for scaled correlation objects on binary-continuous input pairs; always call this to cleanup your data space
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
CORRELATIONLIB_API void FreeScaledCorrelationComputerBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CScaledCorrelationBC*) (pvObject);
	}
}

//Functional interface

///Compute the scaled correlation of one binary signal with one continuos signal
//Parameters for ComputeScaledCorrelation
///@param pvObject					- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piaTimeStampsA			- The vector holding the time stamps of the binary variable
///@param piNrTimeStampsInA			- Number of time stamps in the first vector
///@param pfaSamplesB				- The vector holding samples of the continuous signal
///@param piNrSamplesInB			- Number of samples in the second vector
///@param pbUseFisherZTransform		- Set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
///@warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match! \n
///			The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
CORRELATIONLIB_API void ComputeScaledCorrelationBC(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform)
{
	if(pvObject != NULL)
	{
		((CScaledCorrelationBC*) pvObject)->ComputeScaledCorrelation(piaTimeStampsA,piNrTimeStampsInA,pfaSamplesB,piNrSamplesInB,pbUseFisherZTransform);
	}
}	

///Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the scaled correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
///@param pvObject					- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns a pointer to the internal scaled correlation buffer. Do not write into the buffer! \n
///		   The buffer has a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the original signals.
CORRELATIONLIB_API float* GetScaledCrossCorrelogramBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetScaledCrossCorrelogram();
	}
	else return NULL;
}

///Returns the sum of valid Point Biserial coefficients of correlation for each bin of the correlogram
//
///The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow.
///Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PBsCoefficientSums[i] / PBsCoefficientCounts[i].
///@param pvObject					- Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return A pointer to a buffer holding the sum of correlation coefficients for each bin (without dividing them to the counts; the scaled correlogram = the sum of coefficients / the count of coefficients).
///@warning If no valid coefficient was found for a given bin, the buffer containes a value of NAN at that position. Please check for this in the code!
CORRELATIONLIB_API float* GetPBsCoefficientSumsBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetPBsCoefficientSums();
	}
	else return NULL;
}

///Returns how many Point Biserial coefficients of correlation have been averaged for each bin of the scaled cross correlogram
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return A pointer to the internal buffer containing the counts of Point Biserial correlation coefficients averaged for each bin.
///The size of the buffer is (2*iCorrelationWindow+1), and the count of valid Point Biserial coefficients, that were averaged to compute the scaled correlogram at lag 0, is at position iCorrelationWindow
CORRELATIONLIB_API int* GetPBsCoefficientCountsBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetPBsCoefficientCounts();
	}
	else return NULL;
}

///Returns the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piNumberOfBins - The number of bins in the distribution (output parameter)
///@param pfBinSize - The size of one bin of the distribution (output parameter)
///@return A pointer to the internal buffer that holds the distribution of correlation coefficients
///@return The number of bins in the distribution (piNumberOfBins); this is equal to the size of the buffer
///@return The size of one bin of the distribution (pfBinSize); the size of the bin is in sampling units of the original signals
CORRELATIONLIB_API int* GetDistributionOfCorrelationCoefficientsBC(void* pvObject,int& piNumberOfBins,float& pfBinSize)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetDistributionOfCorrelationCoefficients(piNumberOfBins,pfBinSize);
	}
	else return NULL;
}

///Changes the scale window for the scaled correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piNewScale - The new size of the scale segment
///@return Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size
CORRELATIONLIB_API int ScaledCorrelationModifyScaleWindowBC(void* pvObject,int piNewScale)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->ModifyScaleWindow(piNewScale);
	}
	else return ERR_NULL_OBJECT;
}

///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter; 
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piNewCorrelationWindow - The new size of the correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial
CORRELATIONLIB_API int ScaledCorrelationModifyCorrelationWindowBC(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}

///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the scale window; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyTrialLengthBC(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}

///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@param piNewScale - The new size of the scale segment
///@param piNewCorrelationWindow - The new size of the correlation window
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success, -1 if the new scale window is too small; -2 if the new correlation window is too small; -3 if the new scale window is larger than the trial size; -4 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int ScaledCorrelationModifyAllParametersBC(void* pvObject,int piNewScale,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->ModifyAllParameters(piNewScale,piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}

///Get the size of the current scale window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns the size of the scale segment window
CORRELATIONLIB_API int ScaledCorrelationGetScaleWindowBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetScaleWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns the size of the correlation window
CORRELATIONLIB_API int ScaledCorrelationGetCorrelationWindowBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns the length of the trial in original sampling units
CORRELATIONLIB_API int ScaledCorrelationGetTrialLengthBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}

///Get the number of bins of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns the number of bins used to compute the distribution of correlation coefficients; the value equals the size of the internal buffer that stores the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
CORRELATIONLIB_API int GetDistributionOfCorrelationCoefficientsBinNrBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetDistributionOfCorrelationCoefficientsBinNr();
	}
	else return ERR_NULL_OBJECT;
}

///Get the size of a bin of the distribution of coefficients of correlation
//
///@param pvObject - Pointer to the object that computes the scaled correlation; please instantiate first using the CreateScaledCorrelationComputerBC function
///@return Returns the size of the bins used to compute the distribution of correlation coefficients (see also function: GetDistributionOfCorrelationCoefficients)
CORRELATIONLIB_API float GetDistributionOfCorrelationCoefficientsBinSizeBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CScaledCorrelationBC*) pvObject)->GetDistributionOfCorrelationCoefficientsBinSize();
	}
	else return ERR_NULL_OBJECT;
}

//====================================================================================================================================================================================





//Interface for Binary - Binary Cross corr ======================================================================================================================================


///Creates an object that computes cross-correlation on a pair of binary - binary signals
//Constructor parameters are as follows:
///@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with laggs from -80 to +80); use the same units as the sampling of your samples
///@param piTrialLength			- The size of the trial in units of your samples
///@param piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogramBB() function will return NULL
///@return An untyped pointer to the created object!
CORRELATIONLIB_API void* CreateCrossCorrelationComputerBB(int piCorrelationWindow,int piTrialLength,int &piIsError)
{
	return new CCrossCorrelationBB(piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for cross-correlation objects on binary-binary input pairs; always call this to cleanup your data space
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
CORRELATIONLIB_API void FreeCrossCorrelationComputerBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CCrossCorrelationBB*) (pvObject);
	}
}


///Compute the cross-correlation of two vectors of time stamps
//Parameters for ComputeCrossCorrelation
///@param	pvObject				- Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@param	piaTimeStampsA			- The vector holding the time stamps of the first binary variable
///@param	piNrTimeStampsInA		- Number of time stamps in the first vector
///@param	piaTimeStampsB			- The vector holding the time stamps of the second binary variable
///@param	piNrTimeStampsInB		- Number of time stamps in the second vector
///@param	pbNormalizeCorrelogram	- Set to 1 to normalize the correlogram; normalization is made by dividing the correlogram to the number of time stamps in the first buffer (rate-normalized correlogram)
CORRELATIONLIB_API void ComputeCrossCorrelationBB(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int pbNormalizeCorrelogram)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationBB*) pvObject)->ComputeCrossCorrelation(piaTimeStampsA,piNrTimeStampsInA,piaTimeStampsB,piNrTimeStampsInB,pbNormalizeCorrelogram);
	}
}


///Compute the cross-correlation of two vectors of time stamps, on a partial window of the trial. ONLY ACCEPTS ONE TRIAL!!!
//Parameters for ComputeWindowedCrossCorrelationPerTrial
///@param pvObject				- Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@param piaTimeStampsA		 - The vector holding the time stamps of the first binary variable
///@param piNrTimeStampsInA		 - Number of time stamps in the first vector
///@param piaTimeStampsB		 - The vector holding the time stamps of the second binary variable
///@param piNrTimeStampsInB		 - Number of time stamps in the second vector		
///@param piFromOffsetInTrial	 - The start offset in the trial where the desired window starts (time stamp)
///@param piToOffsetInTrial		 - The end offset in the trial where the desired window stops (time stamp)
///@param pbNormalizeCorrelogram - Set to 1 to normalize the correlogram
CORRELATIONLIB_API void ComputeWindowedCrossCorrelationPerTrialBB(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbNormalizeCorrelogram)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationBB*) pvObject)->ComputeWindowedCrossCorrelationPerTrial(piaTimeStampsA,piNrTimeStampsInA,piaTimeStampsB,piNrTimeStampsInB,piFromOffsetInTrial,piToOffsetInTrial,pbNormalizeCorrelogram);
	}
}
		

///Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow.
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///		   The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
CORRELATIONLIB_API float* GetCrossCorrelogramBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->GetCrossCorrelogram();
	}
	else return NULL;
}

///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
//
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@param piNewCorrelationWindow - The new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial
CORRELATIONLIB_API int CrossCorrelationModifyCorrelationWindowBB(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}

///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyTrialLengthBB(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}

///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@param piNewCorrelationWindow - The new correlation window
///@param piNewTrialLength - The new length of the trial 
///@return Returns 0 on success, -1 if the new correlation window is too small; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyAllParametersBB(void* pvObject,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->ModifyAllParameters(piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}

///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@return The size of the correlation window
CORRELATIONLIB_API int CrossCorrelationGetCorrelationWindowBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross correlation; please instantiate first using the CreateCrossCorrelationComputerBB function
///@return The length of the trial in sampling units
CORRELATIONLIB_API int CrossCorrelationGetTrialLengthBB(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBB*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}

//====================================================================================================================================================================================



//Interface for Continuous - Continuous cross correlation ============================================================================================================================
///Creates an object that computes cross-correlation on a pair of continuous - continuous signals
//Constructor parameters are as follows:
///	@param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with laggs from -80 to +80); use the same units as the sampling of your samples
///	@param piTrialLength		- The size of the trial in units of your samples
///	@param piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL
/// @return An untyped pointer to the created object.
CORRELATIONLIB_API void* CreateCrossCorrelationComputerCC(int piCorrelationWindow,int piTrialLength,int &piIsError)

{
	return new CCrossCorrelationCC(piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for cross-correlation objects on continuous-continuous input pairs; always call this to cleanup your data space
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
CORRELATIONLIB_API void FreeCrossCorrelationComputerCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CCrossCorrelationCC*) (pvObject);
	}
}


///Compute the cross-correlation of two vectors of samples
//Parameters for ComputeCrossCorrelation
///If the number of samples passed is larger than the size of one trial, then the signals are divided into trials and the results are eventually averaged over trials
/// @param pvObject					- Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
/// @param pfaSamplesA				- The vector holding the samples of the first signal
/// @param piNrSamplesInA			- Number of samples in the first vector
/// @param pfaSamplesB				- The vector holding the samples of the second signal
/// @param piNrSamplesInB			- Number of samples in the second vector
/// @param	piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
CORRELATIONLIB_API void ComputeCrossCorrelationCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationCC*) pvObject)->ComputeCrossCorrelation(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,piNormalizationFlag);
	}
}

///Compute the cross-correlation of two vectors of samples without checking for border limits (the two signals must be longer than 2*Correlation window)
//Parameters for ComputeCrossCorrelationNoBorderChecking
///If the number of samples passed is larger than the size of one trial, then the signals are divided into trials and the results are eventually averaged over trials
///WARNING: for each trial, the buffer must contain before and after the trial an additional number of samples equal to the correlation window!
/// @param pvObject					- Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
/// @param pfaSamplesA				- The vector holding the samples of the first signal
/// @param piNrSamplesInA			- Number of samples in the first vector
/// @param pfaSamplesB				- The vector holding the samples of the second signal
/// @param piNrSamplesInB			- Number of samples in the second vector
/// @param	piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
CORRELATIONLIB_API void ComputeCrossCorrelationNoBorderCheckingCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationCC*) pvObject)->ComputeCrossCorrelationNoBorderChecking(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,piNormalizationFlag);
	}
}

///Compute the cross-correlation of two digitized continuous signals, on a partial window of the trial
//Parameters for ComputeWindowedCrossCorrelationPerTrial
/// @warning ONLY ACCEPTS ONE TRIAL so do not pass signals longer than the size of a single trial! Subsequent trials will be discarded.
/// @param pvObject					- Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
/// @param	pfaSamplesA				- The vector holding the samples of the first signal 
/// @param	piNrSamplesInA			- Number of samples in the first vector
/// @param	pfaSamplesB				- The vector holding the samples of the second signal
/// @param	piNrSamplesInB			- Number of samples in the second vector
/// @param	piFromOffsetInTrial		- The start offset in the trial where the desired window starts
/// @param	piToOffsetInTrial		- The end offset in the trial where the desired window stops
/// @param	piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
CORRELATIONLIB_API void ComputeWindowedCrossCorrelationPerTrialCC(void* pvObject,float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int piNormalizationFlag)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationCC*) pvObject)->ComputeWindowedCrossCorrelationPerTrial(pfaSamplesA,piNrSamplesInA,pfaSamplesB,piNrSamplesInB,piFromOffsetInTrial,piToOffsetInTrial,piNormalizationFlag);
	}
}
		

///Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the cross-correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
/// @param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
/// @return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///			The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
CORRELATIONLIB_API float* GetCrossCorrelogramCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->GetCrossCorrelogram();
	}
	else return NULL;
}

///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
///@param piNewCorrelationWindow - Specifies the new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
CORRELATIONLIB_API int CrossCorrelationModifyCorrelationWindowCC(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}


///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyTrialLengthCC(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}


///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
///@param piNewCorrelationWindow - The new correlation window
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success, -1 if the new correlation window is too small; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyAllParametersCC(void* pvObject,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->ModifyAllParameters(piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}


///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
///@return The size of the correlation window
CORRELATIONLIB_API int CrossCorrelationGetCorrelationWindowCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerCC function
///@return The trial length in sampling units
CORRELATIONLIB_API int CrossCorrelationGetTrialLengthCC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationCC*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}
//===================================================================================================================================================================================




//Interface for Binary - Continuous cross correlation ======================================================================================================================================

///Creates an object that computes cross-correlation on a pair of binary - continuous signals
//Constructor parameters are as follows:
/// @param piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use the same units as the sampling units
/// @param piTrialLength		- The size of the trial in original sampling units
/// @param piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL
/// @warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
/// @return An untyped pointer to the created object.
CORRELATIONLIB_API void* CreateCrossCorrelationComputerBC(int piCorrelationWindow,int piTrialLength,int &piIsError)
{
	return new CCrossCorrelationBC(piCorrelationWindow,piTrialLength,piIsError);
}

///Destructor for cross-correlation objects on binary-continuous input pairs; always call this to cleanup your data space
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
CORRELATIONLIB_API void FreeCrossCorrelationComputerBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		delete (CCrossCorrelationBC*) (pvObject);
	}
}

///Compute the cross-correlation of one binary signal with one continuos signal
//Parameters for ComputeCrossCorrelation
///  @param pvObject				- Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///	 @param	piaTimeStampsA			- The vector holding the time stamps of the binary variable
///  @param	piNrTimeStampsInA		- Number of time stamps in the first vector
///  @param	pfaSamplesB				- The vector holding samples of the continuous signal
///  @param	piNrSamplesInB			- Number of samples in the second vector
///	 @param	pbNormalizeCorrelogram	- Set to 1 to normalize the correlogram; ; normalization is made by dividing to the number of time stamps used to compute the correlogram
///  @warning The sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match! \n
///			  The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
CORRELATIONLIB_API void ComputeCrossCorrelationBC(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbNormalizeCorrelogram)
{
	if(pvObject != NULL)
	{
		((CCrossCorrelationBC*) pvObject)->ComputeCrossCorrelation(piaTimeStampsA,piNrTimeStampsInA,pfaSamplesB,piNrSamplesInB,pbNormalizeCorrelogram);
	}
}		

///Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
//
///The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow.
/// @param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
/// @return Returns a pointer to the internal cross-correlation buffer. Do not write into the buffer! \n
///			The buffer is filled with NAN (Not A Number) = -100000000.0 if one of the input buffers was empty and thus the correlation was not defined; please check in your code for this!!!!
CORRELATIONLIB_API float* GetCrossCorrelogramBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->GetCrossCorrelogram();
	}
	else return NULL;
}

///Set the size of the correlation window; for example for a window of -100..+100 pass 100 as a parameter 
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///@param piNewCorrelationWindow - Specifies the new correlation window
///@return Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
CORRELATIONLIB_API int CrossCorrelationModifyCorrelationWindowBC(void* pvObject,int piNewCorrelationWindow)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->ModifyCorrelationWindow(piNewCorrelationWindow);
	}
	else return ERR_NULL_OBJECT;
}

///Set the length of the trial in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///@param piNewTrialLength - The new length of the trial
///@return Returns 0 on success and -1 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyTrialLengthBC(void* pvObject,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->ModifyTrialLength(piNewTrialLength);
	}
	else return ERR_NULL_OBJECT;
}

///Sets all parameters at once
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///@param piNewCorrelationWindow - The new correlation window
///@param piNewTrialLength - The new length of the trial 
///@return Returns 0 on success, -1 if the new correlation window is too small; -2 if the trial length is smaller than the correlation window
CORRELATIONLIB_API int CrossCorrelationModifyAllParametersBC(void* pvObject,int piNewCorrelationWindow,int piNewTrialLength)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->ModifyAllParameters(piNewCorrelationWindow,piNewTrialLength);
	}
	return ERR_NULL_OBJECT;
}

///Get the size of the correlation window
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///@return The size of the correlation window
CORRELATIONLIB_API int CrossCorrelationGetCorrelationWindowBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->GetCorrelationWindow();
	}
	else return ERR_NULL_OBJECT;
}

///Get the trial length in original sampling units
//
///@param pvObject - Pointer to the object that computes the cross-correlation; please instantiate first using the CreateCrossCorrelationComputerBC function
///@return The length of the trial in sampling units
CORRELATIONLIB_API int CrossCorrelationGetTrialLengthBC(void* pvObject)
{
	if(pvObject != NULL)
	{
		return ((CCrossCorrelationBC*) pvObject)->GetTrialLength();
	}
	else return ERR_NULL_OBJECT;
}

//===================================================================================================================================================================================

//Doxygen help comments
/**
* \mainpage
* \section introduction Introduction
* Correlation Library is a C++ library implemented by Raul C. Muresan. The code is free for non-commercial purposes. \n
* For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer.
* \n
* \section description Description
* The Correlation Library offers support for computing cross and scaled correlation between a pair of continuous valued, binary valued, or mixed signals. \n
* For each method of computing the correlation and each combination of signals, there is one specialized class that can handle the respective operations.\n
* In addition to C++ classes, the library also compiles into a Windows DLL that exports a number of functions as wrappers over the original classes.
* \section classe C++ classes
* \subsection naming Naming conventions
* Each class is named as: C_descriptor_sufix, where
* - descriptor - stands for CrossCorrelation or ScaledCorrelation\n
* - sufix - denotes the combination type of data that the class accepts: B = Binary valued discrete signal (represented as a vector of time stamps corresponding to events with value "1"), C = Continuous valued discrete signal (represented as a vector of samples) => CC = Continuous-Continuos; BC = Binary-Continuous; BB = Binary-Binary.
* \subsection classlist List of classes
* - CCrossCorrelationCC:	A class that computes the cross-correlation for two continuous valued signals (Continuous - Continuous)
* - CCrossCorrelationBC:	A class that computes the cross-correlation for a vector of time stamps corresponding to events in a binary signal and a continuous valued signal (Binary - Continuous)
* - CCrossCorrelationBB:	A class that computes the cross-correlation for two vectors of time stamps corresponding to events in binary signals (Binary - Binary)
* - CScaledCorrelationCC:	A class that computes the scaled correlation for two continuous valued signals (Continuous - Continuous)
* - CScaledCorrelationBC:	A class that computes the scaled correlation for a vector of time stamps corresponding to events in a binary signal and a continuous valued signal (Binary - Continuous)
* - CScaledCorrelationBB:	A class that computes the scaled correlation for two vectors of time stamps corresponding to events in binary signals (Binary - Binary)
* \subsection input Input data formatting
* All classes respect the same standard for accepting the input data. For a binary input, a vector of time stamps is expected as input while for a continuous valued input, a digitized vector of samples is expected as input. \n\n
* When multiple measurements (trials) are available for the data, the classes can handle multiple trials. The length of the trial, in original sampling units, has to be specified in the constructor.
* Next, when the input signals are provided in the method call to ComputeCrossCorrelation or ComputeScaledCorrelation, the method determines how many trials are available, computes the correlogram on each trial and finally averages the correlograms into one single correlogram.\n\n
* When multiple trials are available, the input should be structured as follows:\n\n
* |.........Trial 1...........||.........Trial 2...........||.........Trial 3...........||.........Trial 4...........||.........Trial 5...........|\n\n
* - For continuous valued signals, there must be as many samples in each trial as the previously specified trial length. In any case, the method only considers as many samples as specified in the call to ComputeCrossCorrelation or ComputeScaledCorrelation.\n\n
* - For binary signals, the number of available trials is decided by taking into account the time stamp with the largest value. The method then searches, for each trial, the presence of valid time stamps. Time stamps for consecutive trials have to have increasing values. The trial to which a time stamp belongs is decided by taking into account the length of a trial in sampling units.\n\n
*
* If the user cannot structure the trials as a contiguous series of datapoints, then the best option is to pass a single trial to the class, get the correlogram for each trial, and then the user should compute the average correlogram.
* \note For the case of scaled correlation, the user can use the sum of coefficients of correlation for each bin of the correlogram and the count of coefficients to efficiently compute the average, smooth, scaled correlogram for all trials (see section Example). Averaging the sums of correlation coefficients gives a better and smoother estimate of the final scaled correlogram than averaging the individual scaled correlograms for each trial (see Example: Scaled correlogram example with only one trial passed at a time).\n
* \subsection methodcall Method calls
* For the two types of correlations the classes expose an interface to compute the CrossCorrelation and the ScaledCorrelation respectively. Depending on the type of the input data, the methods accept different parameter types. However, they have a general format:\n
* - ComputeCrossCorrelation(VectorA, SizeOfVectorA, VectorB, SizeOfVectorB, pbNormalizeCorrelogram)\n
* - ComputeScaledCorrelation(VectorA, SizeOfVectorA, VectorB, SizeOfVectorB, pbUseFisherZTransform)\n
*\n
* The first four parameters specify the input vectors and their sizes (in number of elements). Each vector has a different type depending if the input signal is binary or continuous valued.\n
* For binary data, the vectors are arrays of integers holding the time stamps. For continuous valued data, the vectors are arrays of floats, holding the samples.\n
* The last parameter of each function model is a boolean. pbNormalizeCorrelogram - determines whether the cross-correlogram should be normalized. pbUseFisherZTransform - determines if the scaled correlogram should use the Fisher Z Transform to average the coefficients of correlation.\n
* 
* \subsection output Output data
* Each class exposes an interface that allows access to the computed cross or scaled correlogram. The correlogram is stored in the object as a vector of floats. The vector has a size of 2*CorrelationWindow+1.\n
* The correlogram vector is stored such that the correlation at lag -CorrelationWindow is at position 0, the correlation at lag 0 is stored at position CorrelationWindow, and the correlation at lag +CorrelationWindow is stored at position 2*CorrelationWindow.\n\n
* To retrieve the correlogram, call:
* - GetCrossCorrelogram() - retrieves a pointer to the vector holding the cross-correlogram \n
* - GetScaledCrossCorrelogram() - retrieves a pointer to the vector holding the scaled correlogram \n
* The user should not write into the vectors! \n
*
* \attention Very important: the cross and scaled correlation buffers might contain values of Not A Number (NAN = -100000000.0f) if there was no variance in the input signals and correlation was, as a consequence, not defined! Always check the output buffers for NAN values before using them! Otherwise results may be compromised.
* \subsection example Examples of usage
* Let us consider we have two continuous valued signals that we want to correlate. The signals, we suppose, have been sampled at 1 kHz. Each signal has 2 trials recorded, with a length of 1000 ms each. We want to compute the cross-correlation and the scaled correlation.\n\n
* \code
* //Cross correlogram example:
* #include "CrossCorrelation-CC.h"
* float SignalA[2000];
* float SignalB[2000];
* float* CrossCorrelogram;
* int iTrialLength = 1000;		//the legth of the trial is 1000 so we have signal vectors holding 2 trials (2000 entries)
* int iCorrelationWindow = 100;		//-100..+100 correlation lags
* 
* //Here you should fill the SignalA and SignalB arrays with values
*
* //And then compute the correlogram
* int isError = 0;
* CCrossCorrelationCC ComputeClassicCC(iCorrelationWindow,iTrialLength,isError);
* if(isError == 0)
* {
*	ComputeClassicCC.ComputeCrossCorrelation(SignalA,2000,SignalB,2000,0);
*	CrossCorrelogram = ComputeClassicCC.GetCrossCorrelogram();
*	for(int j=0;j<2*iCorrelationWindow+1;j++)
*	{
*		if(CrossCorrelogram[j] != NAN)
*		{
*			//use the CrossCorrelogram only if the value is not NAN (only if the correlation was defined)
*		}
*	}
* }
* \endcode
* \n\n
* \code
* //Scaled correlogram example:
* #include "ScaledCorrelation-CC.h"
* float SignalA[2000];
* float SignalB[2000];
* float* ScaledCorrelogram;
* int iTrialLength = 1000;		//the legth of the trial is 1000 so we have signal vectors holding 2 trials (2000 entries)
* int iCorrelationWindow = 100;		//-100..+100 correlation lags
* int iScaleWindow = 25;			//25 ms - if the sampling frequency of the signals is 1 kHz => we are looking to components faster than 40 Hz
* 
* //Here you should fill the SignalA and SignalB arrays with values
*
* //And then compute the correlogram
* int isError = 0;
* CScaledCorrelationCC ComputeScaledCC(iScaleWindow,iCorrelationWindow,iTrialLength,isError);
* if(isError == 0)
* {
*	ComputeScaledCC.ComputeScaledCorrelation(SignalA,2000,SignalB,2000,0);
*	ScaledCorrelogram = ComputeScaledCC.GetScaledCrossCorrelogram();
*	for(int j=0;j<2*iCorrelationWindow+1;j++)
*	{
*		if(ScaledCorrelogram[j] != NAN)
*		{
*			//use the ScaledCorrelogram only if the value is not NAN (only if the correlation was defined)
*		}
*	}
* }
* \endcode
* \n\n
* Let us now assume we have again two signals each with two trials. This time however, we want to compute the correlogram manually by computing the sum of correlation coefficients for each trial and then computing the average in the end.\n\n
* \code
* //Scaled correlogram example with only one trial passed at a time:
* //This is useful when the user wants to compute the correlogram per trial and in the end average the correlation coefficients himself.
* //Averaging the sums of correlation coefficients gives a better and smoother estimate of the final scaled correlogram than averaging the individual scaled correlograms for each trial!!!
* #include "ScaledCorrelation-CC.h"
* #include "Constants.h"
* float SignalA[2000];
* float SignalB[2000];
* float ScaledCorrelogram[201];
* int   SummedCoefficientsCount[201];
* float *coeff_sum_buffer;
* int *coeff_count_buffer;
* int iTrialLength = 1000;		//the legth of the trial is 1000 so we have signal vectors holding 2 trials (2000 entries)
* int iCorrelationWindow = 100;		//-100..+100 correlation lags
* int iScaleWindow = 25;			//25 ms - if the sampling frequency of the signals is 1 kHz => we are looking to components faster than 40 Hz
* int iTrialNr = 2;			//Two trials
* 
* //Here you should fill the SignalA and SignalB arrays with values
*
* //And then compute the correlogram
* int isError = 0;
* CScaledCorrelationCC ComputeScaledCC(iScaleWindow,iCorrelationWindow,iTrialLength,isError);
* 
* if(isError == 0)
* {
*	//Cleanup the final buffer and the counting buffer
*	for(int j=0;j<2*iCorrelationWindow+1;j++)
*	{
*		ScaledCorrelogram[j] = 0;
*		SummedCoefficientsCount[j] = 0;
*	}
*
*	//For each trial, compute the sum of correlation coefficients
*	for(int i=0;i<iTrialNr;i++)
*	{
*		ComputeScaledCC.ComputeScaledCorrelation(&(SignalA[i*iTrialLength]),iTrialLength,&(SignalB[i*iTrialLength]),iTrialLength,0);
*		coeff_sum_buffer = ComputeScaledCC.GetPearsonCoefficientSums();
*		coeff_count_buffer = ComputeScaledCC.GetPearsonCoefficientCounts();
*		for(j=0;j<2*iCorrelationWindow+1;j++)
*		{
*			if(coeff_sum_buffer[j] != NAN) ScaledCorrelogram[j] += coeff_sum_buffer[j];
*			SummedCoefficientsCount[j] += coeff_count_buffer[j];
*		}		
*	}
*
*	//Now compute the final scaled correlogram by computing the average
*	for(j=0;j<2*iCorrelationWindow+1;j++) if(SummedCoefficientsCount[j]) ScaledCorrelogram[j] /= SummedCoefficientsCount[j];
* }
* \endcode
*
* \section dllfunctions Dll functions exported
* The correlation library provides flexible means to interact with the classes that compute correlation.
* It exports a set of wrapper functions that can be used by an external caller to access the full class functionality.
* \subsection interfacingconvention Interfacing convention
* To interface a given class method, the library exports a function with a similar name that has an additional parameter, which is an untyped pointer to an already created object of that class.\n \n
* For example, let ClassA be the class of interest, and Method1 be a method of ClassA. To create an interface to ClassA functionality, one must first create an object of that class and then call a wrapper function to access the class method.\n\n
* \code
* void* Wrapper_ObjectConstructor_ClassA()
* {
*	return new ClassA();
* }
*
* void Wrapper_ObjectDestructor_ClassA(void* pObj)
* {
*	delete ((ClassA*) pObj);
* }
*
* void Wrapper_Method1_Of_ClassA(void* pObj,int iSomeParam)
* {
*	((ClassA*) pObj)->Method1(iSomeParam);
* }	
*
* void main()
* {
*	void *pObject;
*	pObject = Wrapper_ObjectConstructor_ClassA();	//A wrapper over the constructor of class A
*	Wrapper_Method1_Of_ClassA(pObject,1);		//A wrapper over Method1 of class A
*	Wrapper_ObjectDestructor_ClassA(pObject);	//A wrapper over the destructor of class A
* }
* \endcode
* \note The wrapper constructors create untyped pointers to objects. The DLL caller can use these pointers as handles for the objects of a given class. All class operations have, as an additional parameter, the pointer to the object being handled. Depending on the wrapper function that is called, a corresponding method will be called from a corresponding class of the object whose handle is passed to the wrapper function.
* \subsection naming Naming conventions
* Each function is named as: prefix_name_sufix, where
* - prefix - is optional and is usually added to avoid confusion, when a function with the same name is present in both Cross-Correlation and Scaled-Correlation classes. The prefix could be either CrossCorrelation or ScaledCorrelation
* - name - is usually the name of the class method that is being interfaced \n
* - sufix - denotes the combination type of data that the class accepts: B = Binary, C = Continuous => CC = Continuous-Continuos; BC = Binary-Continuous; BB = Binary-Binary.
* \subsection wrapperfunctionlist List of wrapper functions
*
* Cross-correlation for Continuous-Continuous data:
* - void* 	CreateCrossCorrelationComputerCC(int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void 	FreeCrossCorrelationComputerCC(void *pvObject)
* - void 	ComputeCrossCorrelationCC(void *pvObject, float *pfaSamplesA, int piNrSamplesInA, float *pfaSamplesB, int piNrSamplesInB, int pbNormalizeCorrelogram)
* - void 	ComputeWindowedCrossCorrelationPerTrialCC(void *pvObject, float *pfaSamplesA, int piNrSamplesInA, float *pfaSamplesB, int piNrSamplesInB, int piFromOffsetInTrial, int piToOffsetInTrial, int pbNormalizeCorrelogram)
* - float* 	GetCrossCorrelogramCC(void *pvObject)
* - int 	CrossCorrelationModifyCorrelationWindowCC(void *pvObject, int piNewCorrelationWindow)
* - int 	CrossCorrelationModifyTrialLengthCC(void *pvObject, int piNewTrialLength)
* - int 	CrossCorrelationModifyAllParametersCC(void *pvObject, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	CrossCorrelationGetCorrelationWindowCC(void *pvObject)
* - int 	CrossCorrelationGetTrialLengthCC(void *pvObject)
*
* Cross-correlation for Binary-Continuous data:
* - void* 	CreateCrossCorrelationComputerBC(int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void 	FreeCrossCorrelationComputerBC(void *pvObject)
* - void 	ComputeCrossCorrelationBC(void *pvObject, int *piaTimeStampsA, int piNrTimeStampsInA, float *pfaSamplesB, int piNrSamplesInB, int pbNormalizeCorrelogram)
* - float* 	GetCrossCorrelogramBC(void *pvObject)
* - int 	CrossCorrelationModifyCorrelationWindowBC(void *pvObject, int piNewCorrelationWindow)
* - int 	CrossCorrelationModifyTrialLengthBC(void *pvObject, int piNewTrialLength)
* - int 	CrossCorrelationModifyAllParametersBC(void *pvObject, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	CrossCorrelationGetCorrelationWindowBC(void *pvObject)
* - int 	CrossCorrelationGetTrialLengthBC(void *pvObject)
*
* Cross-correlation for Binary-Binary data:
* - void* 	CreateCrossCorrelationComputerBB(int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void 	FreeCrossCorrelationComputerBB(void *pvObject)
* - void 	ComputeCrossCorrelationBB(void *pvObject, int *piaTimeStampsA, int piNrTimeStampsInA, int *piaTimeStampsB, int piNrTimeStampsInB, int pbNormalizeCorrelogram)
* - void	ComputeWindowedCrossCorrelationPerTrialBB(void* pvObject,int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbNormalizeCorrelogram)
* - float* 	GetCrossCorrelogramBB(void *pvObject)
* - int 	CrossCorrelationModifyCorrelationWindowBB(void *pvObject, int piNewCorrelationWindow)
* - int 	CrossCorrelationModifyTrialLengthBB(void *pvObject, int piNewTrialLength)
* - int 	CrossCorrelationModifyAllParametersBB(void *pvObject, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	CrossCorrelationGetCorrelationWindowBB(void *pvObject)
* - int 	CrossCorrelationGetTrialLengthBB(void *pvObject)
*
* Scaled correlation for Continuous-Continuous data:
* - void* 	CreateScaledCorrelationComputerCC(int piScaleWindow, int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void 	FreeScaledCorrelationComputerCC(void *pvObject)
* - void 	ComputeScaledCorrelationCC(void *pvObject, float *pfaSamplesA, int piNrSamplesInA, float *pfaSamplesB, int piNrSamplesInB, int pbUseFisherZTransform)
* - void 	ComputeWindowedScaledCorrelationPerTrialCC(void *pvObject, float *pfaSamplesA, int piNrSamplesInA, float *pfaSamplesB, int piNrSamplesInB, int piFromOffsetInTrial, int piToOffsetInTrial, int pbUseFisherZTransform)
* - float* GetScaledCrossCorrelogramCC(void *pvObject)
* - float* GetPearsonCoefficientSumsCC(void *pvObject)
* - int* 	GetPearsonCoefficientCountsCC(void *pvObject)
* - int* 	GetDistributionOfCorrelationCoefficientsCC(void *pvObject, int &piNumberOfBins, float &pfBinSize)
* - int 	ScaledCorrelationModifyScaleWindowCC(void *pvObject, int piNewScale)
* - int 	ScaledCorrelationModifyCorrelationWindowCC(void *pvObject, int piNewCorrelationWindow)
* - int 	ScaledCorrelationModifyTrialLengthCC(void *pvObject, int piNewTrialLength)
* - int 	ScaledCorrelationModifyAllParametersCC(void *pvObject, int piNewScale, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	ScaledCorrelationGetScaleWindowCC(void *pvObject)
* - int 	ScaledCorrelationGetCorrelationWindowCC(void *pvObject)
* - int 	ScaledCorrelationGetTrialLengthCC(void *pvObject)
* - int 	GetDistributionOfCorrelationCoefficientsBinNrCC(void *pvObject)
* - float 	GetDistributionOfCorrelationCoefficientsBinSizeCC(void *pvObject)
*
* Scaled correlation for Binary-Continuous data:
* - void* 	CreateScaledCorrelationComputerBC(int piScaleWindow, int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void 	FreeScaledCorrelationComputerBC(void *pvObject)
* - void 	ComputeScaledCorrelationBC(void *pvObject, int *piaTimeStampsA, int piNrTimeStampsInA, float *pfaSamplesB, int piNrSamplesInB, int pbUseFisherZTransform)
* - float* 	GetScaledCrossCorrelogramBC(void *pvObject)
* - float* 	GetPBsCoefficientSumsBC(void *pvObject)
* - int* 	GetPBsCoefficientCountsBC(void *pvObject)
* - int* 	GetDistributionOfCorrelationCoefficientsBC(void *pvObject, int &piNumberOfBins, float &pfBinSize)
* - int		ScaledCorrelationModifyScaleWindowBC(void *pvObject, int piNewScale)
* - int 	ScaledCorrelationModifyCorrelationWindowBC(void *pvObject, int piNewCorrelationWindow)
* - int 	ScaledCorrelationModifyTrialLengthBC(void *pvObject, int piNewTrialLength)
* - int 	ScaledCorrelationModifyAllParametersBC(void *pvObject, int piNewScale, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	ScaledCorrelationGetScaleWindowBC(void *pvObject)
* - int 	ScaledCorrelationGetCorrelationWindowBC(void *pvObject)
* - int 	ScaledCorrelationGetTrialLengthBC(void *pvObject)
* - int 	GetDistributionOfCorrelationCoefficientsBinNrBC(void *pvObject)
* - float 	GetDistributionOfCorrelationCoefficientsBinSizeBC(void *pvObject)
*
* Scaled correlation for Binary-Binary data:
* - void*	CreateScaledCorrelationComputerBB(int piScaleWindow, int piCorrelationWindow, int piTrialLength, int &piIsError)
* - void	FreeScaledCorrelationComputerBB(void *pvObject)
* - void	ComputeScaledCorrelationBB(void *pvObject, int *piaTimeStampsA, int piNrTimeStampsInA, int *piaTimeStampsB, int piNrTimeStampsInB, int pbUseFisherZTransform)
* - void	ComputeWindowedScaledCorrelationPerTrialBB(void *pvObject, int *piaTimeStampsA, int piNrTimeStampsInA, int *piaTimeStampsB, int piNrTimeStampsInB, int piFromOffsetInTrial, int piToOffsetInTrial, int pbUseFisherZTransform)
* - float*	GetScaledCrossCorrelogramBB(void *pvObject)
* - float*	GetFiCoefficientSumsBB(void *pvObject)
* - int*	GetFiCoefficientCountsBB(void *pvObject)
* - int* 	GetDistributionOfCorrelationCoefficientsBB(void *pvObject, int &piNumberOfBins, float &pfBinSize)
* - int 	ScaledCorrelationModifyScaleWindowBB(void *pvObject, int piNewScale)
* - int 	ScaledCorrelationModifyCorrelationWindowBB(void *pvObject, int piNewCorrelationWindow)
* - int 	ScaledCorrelationModifyTrialLengthBB(void *pvObject, int piNewTrialLength)
* - int 	ScaledCorrelationModifyAllParametersBB(void *pvObject, int piNewScale, int piNewCorrelationWindow, int piNewTrialLength)
* - int 	ScaledCorrelationGetScaleWindowBB(void *pvObject)
* - int 	ScaledCorrelationGetCorrelationWindowBB(void *pvObject)
* - int 	ScaledCorrelationGetTrialLengthBB(void *pvObject)
* - int 	GetDistributionOfCorrelationCoefficientsBinNrBB(void *pvObject)
* - float 	GetDistributionOfCorrelationCoefficientsBinSizeBB(void *pvObject)
*/



