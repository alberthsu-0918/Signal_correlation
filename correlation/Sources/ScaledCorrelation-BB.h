/*
	++ Source:
			ScaledCorrelation-BB.h

	++ Description:
			A class that computes the scaled correlation for two vectors of time stamps (Binary - Binary)
			For details on the Scaled Correlation algorithm see paper "Scaled Correlation Analysis" by Danko Nikolic, Wolf Singer, and Raul C. Muresan

	++ History:
			05-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


#ifndef USED_SCALEDCORRELATION_BB
#define USED_SCALEDCORRELATION_BB

#define BIN_SIZE_OF_FI_DISTRIBUTION 0.01		//The size of a bin in the distribution of correlation coefficients

//Doxygen comment
/**
* \brief
* A class that computes the scaled correlation for two vectors of time stamps (Binary - Binary).
* \details
* For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer.\n
*/

class CScaledCorrelationBB
{
	private:
		int		iScaleWindow;					///The size of the scale window
		int		iCorrelationWindow;				///The size of the correlation window (for example to compute cross correlation in a window of -100..100, iCorrelationWindow = 100;
		int		iTrialLength;					///The size of the trial in original sampling time units, same as the time stamps in the input vectors
		float	*faScaledCorrelogram;			///The scale correlogram having a size of (2*iCorrelationWindow + 1), with lag 0 at faScaledCorrelogram[iCorrelationWindow]
		float	*faFiCoefficientsSum;			///For each lag in the correlogram it summs up the fi coefficients that were used to compute the correlogram (faScaledCorrelogram[i] = faFiCoefficientsSum[i] / iaFiCoefficientsCount[i])
		int		*iaFiCoefficientsCount;			///Counts how many fi coefficients have been added to each lag of the scale correlogram
		int		*iaFiCoefficientsDistribution;	///Distribution of correlation coefficients that is also calculated with the cross correlation measure; at index 0 we have the count for coefficients == -1, while at iFiCorrCoeffDistributionBins we have the count for coefficients == +1
		int		*iaTimeStampsA, *iaTimeStampsB;	///Internal pointers to the current sources provided in the call to ComputeScaledCorrelation
		int		bUseFisherZTransform;			///Flag that specifies whether Fisher Z Transform should be applied to average the correlation coefficients
		int     iFiCorrCoeffDistributionBins;	///Number of bins in the distribution of fi correlation coefficients; they span the interval of [-1..1]
		float	fBinSizeOfCorrCoeffDistribution;///The size of a bin in the distribution of correlation coefficients

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
		///Determines the "Fi" (Pearson's) coefficient of correlation for a pair of scale segments;
		float ComputeFiCorrelation(int* piaBufferA,int piBufferASize,int* piaBufferB,int piBufferBSize,int iOffsetB);					
		///The function that computes one pass scale correlation for a pair of trials - always segments relative to the beginning of signal A (reference); the autocorrelation is NOT symmetric
		void  ComputeTrialScaledCorrelation_Asymmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);	
		///The function that computes one pass scale correlation for a pair of trials - segments the signals such that segments are aligned to signal B for right shift and signal A for left shift; the autocorrelation is symmetric
		void  ComputeTrialScaledCorrelation_Symmetric(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);	
		///Same as ComputeTrialScaledCorrelation_Symmetric but faster version with local expansion and use of a memory buffer
		void  ComputeTrialScaledCorrelation_Symmetric_Fast(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);	

	public:
		//Construction interface
		
		///Constructor
		//Constructor parameters are as follows:
		//  piScaleWindow		- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
		//  piCorrelationWindow	- The size of the cross correlation window (eg. 80 for a cross correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
		//	piTrialLength		- The size of the trial in units of your time stamps
		//	piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
		CScaledCorrelationBB(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError);
		
		///Destructor
		~CScaledCorrelationBB();

		//Functional interface
		
		///Compute the scaled correlation of two vectors of time stamps
		//Parameters for ComputeScaledCorrelation
		//	piaTimeStampsA			- the vector holding the time stamps of the first neuron
		//  piNrTimeStampsInA		- number of time stamps in the first vector
		//  piaTimeStampsB			- the vector holding the time stamps of the second neuron
		//  piNrTimeStampsInB		- number of time stamps in the second vector
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		void ComputeScaledCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int pbUseFisherZTransform);

		///Compute the scaled correlation of two vectors of time stamps, on a partial window of the trial
		//ONLY ACCEPTS ONE TRIAL!!!
		//Parameters for ComputeWindowedScaledCorrelationPerTrial
		//	piaTimeStampsA			- the vector holding the time stamps of the first neuron
		//  piNrTimeStampsInA		- number of time stamps in the first vector
		//  piaTimeStampsB			- the vector holding the time stamps of the second neuron
		//  piNrTimeStampsInB		- number of time stamps in the second vector
		//  piFromOffsetInTrial		- the start offset in the trial where the desired window starts
		//	piToOffsetInTrial		- the end offset in the trial where the desired window stops
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		void ComputeWindowedScaledCorrelationPerTrial(int *piaTimeStampsA,int piNrTimeStampsInA,int *piaTimeStampsB,int piNrTimeStampsInB,int piFromOffsetInTrial,int piToOffsetInTrial,int pbUseFisherZTransform);
	
		///Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
		//The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
		//Returns a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the original signals
		float* GetScaledCrossCorrelogram();	

		///Returns the sum of valid Fi coefficients of correlation for each bin of the correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow;
		//Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = FiCoefficientSums[i] / FiCoefficientCounts[i]
		float* GetFiCoefficientSums();		

		///Returns how many Fi coefficients of correlation have been averaged for each bin of the scaled cross correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the count of coefficients averaged for lag 0 is at position iCorrelationWindow;
		int* GetFiCoefficientCounts();

		///Returns the distribution of coefficients of correlation
		//Outputs are: the buffer with the distribution of correlation coefficients, the number of bins in the distribution and the size of one bin
		int* GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize);
		
		///Changes the scale segment size for the scaled correlation
		//Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size;
		int ModifyScaleWindow(int piNewScale);
		
		///Set the size of the correlation window; the parameter specifies the of the correlation window; for example for a window of -100..+100 pass 100 as a parameter
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