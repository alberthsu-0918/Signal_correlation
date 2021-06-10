/*
	++ Source:
			CrossCorrelation-CC.h

	++ Description:
			A class that computes the cross-correlation for two continuous signals (Continuous - Continuous)
			Cross-Correlation is computed in the classical way by counting

	++ History:
			19-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


//Doxygen comment
/**
* \brief
* A class that computes the cross-correlation for two continuous signals (Continuous - Continuous).
* \details
* Cross-correlation is computed in the classical way.\n
* Two continuous signals of equal length must be provided. The signals should be digitized with the same sampling frequency.
*/


#ifndef USED_CROSSCORRELATION_CC
#define USED_CROSSCORRELATION_CC

class CCrossCorrelationCC
{
	private:
		int		iCorrelationWindow;				///The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		int		iTrialLength;					///The size of the trial in original sampling time units, same as the samples in the input vectors
		float	*faCrossCorrelogram;			///The cross correlogram having a size of (2*iCorrelationWindow + 1), with lag 0 at faCrossCorrelogram[iCorrelationWindow]
		int		*iaNumberOfSummedCorrelations;	///The total number of correlations summed for each correation bin; used to normalize the cross correlogram;
		float	*faSignalA, *faSignalB;			///Internal pointers to the current sources provided in the call to ComputeCrossCorrelations

		//Operational interface
		///Checks if the sizes of the windows are consistent
		int	  CheckWindowsConsistency(int piCorrW,int piTrialLength);											
		///Allocates the space for the cross correlogram data
		void  AllocateData();																								
		///Cleans up the memory occupied by the correlogram data
		void  CleanupData();																								
		///Clears the bytes in the cross correlogram
		void  ResetCrossCorrelogram();	
		///Sets the correlogram to Not A Number (NAN) values to signal an error or undefined correlation
		void  SetCorrelogramToNAN();
		///The function that computes cross-correlation for a pair of trials
		void  ComputeTrialCrossCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB);
		///The function that computes one pass cross-correlation for a pair of trials without checking for borders
		void ComputeTrialCrossCorrelationNoBorderChecking(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB);

	public:
		//Construction interface
		
		///Constructor
		//Constructor parameters are as follows:
		//	piCorrelationWindow		- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use again the same units as the sampling of your samples
		//	piTrialLength			- The size of the trial in units of your samples
		//	piIsError				- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL.
		CCrossCorrelationCC(int piCorrelationWindow,int piTrialLength,int &piIsError);
		
		///Destructor
		~CCrossCorrelationCC();

		//Functional interface
		
		///Compute the cross-correlation of two vectors of samples
		//Parameters for ComputeCrossCorrelation
		//	pfaSamplesA				- the vector holding the samples of the first discrete variable
		//  piNrSamplesInA			- number of samples in the first vector
		//  pfaSamplesB				- the vector holding the samples of the second continuous variable
		//  piNrSamplesInB			- number of samples in the second vector
		//	piNormalizationFlag		- flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased
		void ComputeCrossCorrelation(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag);

		//Compute the cross-correlation of two vectors of samples without checking for borders
		//WARNING: for each trial, the buffer must contain before and after the trial a number of samples equal to the correlation window!
		///If the number of samples passed is larger than the size of one trial, then the signals are divided into trials and the results are eventually averaged over trials
		//Parameters for ComputeCrossCorrelation
		/// @param pfaSamplesA				- The vector holding the samples of the first signal
		/// @param piNrSamplesInA			- Number of samples in the first vector
		/// @param pfaSamplesB				- The vector holding the samples of the second signal
		/// @param piNrSamplesInB			- Number of samples in the second vector
		/// @param piNormalizationFlag		- Flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased; normalization removes the dependency on the length of the signals and gives values in cross-product units of the continuous variables; unbias removes the dependency on the finite length of signal windows
		void ComputeCrossCorrelationNoBorderChecking(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piNormalizationFlag);

		///Compute the cross-correlation of two digitized continuous signals, on a partial window of the trial. ONLY ACCEPTS ONE TRIAL!!!
		//Parameters for ComputeWindowedCrossCorrelationPerTrial
		//	pfaSamplesA				- the vector holding the first signal's input samples
		//  piNrSamplesInA			- number of samples in the first vector
		//  pfaSamplesB				- the vector holding the second signal's input samples
		//  piNrSamplesInB			- number of samples in the second vector
		//  piFromOffsetInTrial		- the start offset in the trial where the desired window starts
		//	piToOffsetInTrial		- the end offset in the trial where the desired window stops
		//	piNormalizationFlag		- flag that decides how to normalize the correlogram: 0 - unnormalized, biased; 1 - normalized, biased; 2 - unnormalized, unbiased; 3 - normalized, unbiased
		void ComputeWindowedCrossCorrelationPerTrial(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbNormalizeCorrelogram);
		
		///Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
		//The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
		//Returns NAN (Not A Number) = -100000000.0 if one of the buffers is empty and thus the correlation is not defined; please check in your code for this!!!!
		float* GetCrossCorrelogram();	
		
		///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
		//Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
		int ModifyCorrelationWindow(int piNewCorrelationWindow);

		///Set the length of the trial in original sampling units
		//Returns 0 on success and -1 if the trial length is smaller than the correlation window
		int ModifyTrialLength(int piNewTrialLength);

		///Sets all parameters at once
		//Returns 0 on success, -1 if the new correlation window is too small; -2 if the trial length is smaller than the correlation window
		int ModifyAllParameters(int piNewCorrelationWindow,int piNewTrialLength);

		///Get the size of the correlation window
		int GetCorrelationWindow();

		///Get the trial length in original sampling units
		int GetTrialLength();
};
#endif