/*
	++ Source:
			ScaledCorrelation-CC.h

	++ Description:
			A class that computes the scaled correlation for two continuous signals (Continuous - Continuous)
			For details on the Scaled Correlation algorithm see paper "Scaled Correlation Analysis" by Danko Nikolic, Wolf Singer, and Raul C. Muresan

	++ History:
			17-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


#ifndef USED_SCALEDCORRELATION_CC
#define USED_SCALEDCORRELATION_CC

#define BIN_SIZE_OF_COEFF_DISTRIBUTION 0.01			//The size of a bin in the distribution of correlation coefficients

//Doxygen comment
/**
* \brief
* A class that computes the scaled correlation for two continuous signals (Continuous - Continuous).
* \details
* For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer.\n
* Two continuous signals of equal length must be provided. The signals should be digitized with the same sampling frequency.
*/

class CScaledCorrelationCC
{
	private:
		int		iScaleWindow;						///The size of the scale window
		int		iCorrelationWindow;					///The size of the correlation window (for example to compute cross correlation in a window of -100..100, iCorrelationWindow = 100;
		int		iTrialLength;						///The size of the trial in original sampling time units
		float	*faScaledCorrelogram;				///The scale correlogram having a size of (2*iCorrelationWindow + 1), with lag 0 at faScaledCorrelogram[iCorrelationWindow]
		float	*faPearsonCoefficientsSum;			///For each lag in the correlogram it summs up the Pearson coefficients that were used to compute the correlogram (faScaledCorrelogram[i] = faPearsonCoefficientsSum[i] / iaPearsonCoefficientsCount[i])
		int		*iaPearsonCoefficientsCount;		///Counts how many Pearson coefficients have been added to each lag of the scale correlogram
		int		*iaPearsonCoefficientsDistribution;	///Distribution of correlation coefficients that is also calculated with the cross correlation measure; at index 0 we have the count for coefficients == -1, while at iPearsonCorrCoeffDistributionBins we have the count for coefficients == +1
		float	*faSignalA, *faSignalB;				///Internal pointers to the current sources provided in the call to ComputeScaledCorrelations
		int		bUseFisherZTransform;				///Flag that specifies whether Fisher's Z Transform should be applied to average the correlation coefficients
		int     iPearsonCorrCoeffDistributionBins;	///Number of bins in the distribution of Pearson correlation coefficients; they span the interval of [-1..1]
		float	fBinSizeOfCorrCoeffDistribution;	///The size of a bin in the distribution of (Pearson's) correlation coefficients

		//Operational interface
		///Checks if the sizes of the windows are consistent
		int	  CheckWindowsConsistency(int piScaleW,int piCorrW,int piTrialLength);											
		///Allocates the space for the cross correlogram data
		void  AllocateData();																								
		///Cleans up the memory occupied by the correlogram data
		void  CleanupData();																								
		///Clears the bytes in the cross correlogram
		void  ResetScaledCorrelogram();			
		///Sets the correlogram to Not A Number values to signal an error or undefined correlation
		void  SetCorrelogramToNAN();
		///Determines the Pearson's coefficient of correlation for a pair of segments;
		float ComputePearsonCorrelation(float* pfaBufferA,int piBufferASize,float* pfaBufferB,int piBufferBSize);					
		///The function that computes one pass scale correlation for a pair of trials - always segments relative to the beginning of signal A (reference); the autocorrelation is NOT symmetric
		void  ComputeTrialScaledCorrelation_Asymmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);
		///The function that computes one pass scale correlation for a pair of trials - segments the signals such that segments are aligned to signal B for right shift and signal A for left shift; the autocorrelation is symmetric
		void  ComputeTrialScaledCorrelation_Symmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);	
		///The function that computes one pass scale correlation for a pair of trials; there is no border cheking on the start and end of the trial - one must provide longer buffers by at least the size of the correlation window; the function provides symmetric autocorrelation
		void  ComputeTrialScaledCorrelationNoBorderChecking(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);

	public:
		//Construction interface
		
		///Constructor
		//Constructor parameters are as follows:
		//  piScaleWindow		- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your signals
		//  piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use again the same units as the sampling of your signals
		//	piTrialLength		- The size of the trial in sampling units of your signals
		//	piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
		CScaledCorrelationCC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError);
		
		///Destructor
		~CScaledCorrelationCC();

		//Functional interface
		
		///Compute the scaled correlation of two digitized continuous signals
		//Parameters for ComputeScaledCorrelation
		//	pfaSamplesA				- the vector holding the first signal's input samples
		//  piNrSamplesInA			- number of samples in the first vector
		//  pfaSamplesB				- the vector holding the second signal's input samples
		//  piNrSamplesInB			- number of samples in the second vector
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		void ComputeScaledCorrelation(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform);

		///Compute the scaled correlation of two digitized continuous signals without checking for borders at the end of the trials
		//WARNING: for each trial, the buffer must contain before and after the trial a number of samples equal to the correlation window!
		//Parameters for ComputeScaledCorrelationNoBorderChecking
		//	pfaSamplesA				- the vector holding the first signal's input samples
		//  piNrSamplesInA			- number of samples in the first vector
		//  pfaSamplesB				- the vector holding the second signal's input samples
		//  piNrSamplesInB			- number of samples in the second vector
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		void ComputeScaledCorrelationNoBorderChecking(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform);
		
		///Compute the scaled correlation of two digitized continuous signals, on a partial window of the trial; ONLY ACCEPTS ONE TRIAL!!!
		//Parameters for ComputeWindowedScaledCorrelationPerTrial
		//	pfaSamplesA				- the vector holding the first signal's input samples
		//  piNrSamplesInA			- number of samples in the first vector
		//  pfaSamplesB				- the vector holding the second signal's input samples
		//  piNrSamplesInB			- number of samples in the second vector
		//  piFromOffsetInTrial		- the start offset in the trial where the desired window starts
		//	piToOffsetInTrial		- the end offset in the trial where the desired window stops
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		void ComputeWindowedScaledCorrelationPerTrial(float *pfaSamplesA,int piNrSamplesInA,float *pfaSamplesB,int piNrSamplesInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform);

		///Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
		//The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
		//Returns a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the original signals
		float* GetScaledCrossCorrelogram();	

		///Returns the sum of valid Pearson coefficients of correlation for each bin of the correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow;
		//Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PearsonCoefficientSums[i] / PearsonCoefficientCounts[i]
		float* GetPearsonCoefficientSums();		
		
		///Returns how many Pearson's coefficients of correlation have been averaged for each bin of the scaled cross correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the count of coefficients averaged for lag 0 is at position iCorrelationWindow;
		int* GetPearsonCoefficientCounts();

		///Returns the distribution of coefficients of correlation
		//Outputs are: the buffer with the distribution of correlation coefficients, the number of bins in the distribution and the size of one bin
		int* GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize);
		
		///Changes the size of the scale segment used to compute the scaled correlation
		//Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size;
		int ModifyScaleWindow(int piNewScale);
		
		///Set the size of the correlation window; the parameter specifies the correlation window; for example for lags between -100 to +100 pass 100 as a parameter
		//Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
		int ModifyCorrelationWindow(int piNewCorrelationWindow);

		///Set the length of the trial in original sampling units
		//Returns 0 on success and -1 if the trial length is smaller than the scale window; -2 if the trial length is smaller than the correlation window
		int ModifyTrialLength(int piNewTrialLength);

		///Sets all parameters at once
		//Returns 0 on success, -1 if the new scale window is too small; -2 if the new correlation window is too small; -3 if the new scale window is larger than the trial size; -4 if the trial length is smaller than the correlation window
		int ModifyAllParameters(int piNewScale,int piNewCorrelationWindow,int piNewTrialLength);

		///Get the size of the current scale segment
		int GetScaleWindow();

		///Get the size of the correlation window
		int GetCorrelationWindow();

		///Get the trial length in original sampling units
		int GetTrialLength();

		///Get the number of bins of the distribution of coefficients of correlation
		int GetDistributionOfCorrelationCoefficientsBinNr();

		///Get the size of a bin of the distribution of coefficients of correlation
		float GetDistributionOfCorrelationCoefficientsBinSize();
};
#endif
