/*
	++ Source:
			ScaledCorrelation-BC.h

	++ Description:
			A class that computes the scaled correlation for a vector of time stamps and a continuous signal (Binary - Continuous)
			For details on the Scaled Correlation algorithm see paper "Scaled Correlation Analysis" by Danko Nikolic, Wolf Singer, and Raul C. Muresan

	++ History:
			21-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/


#ifndef USED_SCALEBCORRELATION_BC
#define USED_SCALEBCORRELATION_BC

#define BIN_SIZE_OF_COEFF_DISTRIBUTION 0.01		//The size of a bin in the distribution of correlation coefficients

//Doxygen comment
/**
* \brief
* A class that computes the scaled correlation for a vector of time stamps and a continuous signal (Binary - Continuous).
* \details
* For details on the Scaled Correlation algorithm see paper on "Scaled Correlation Analysis" by Danko Nikolic, Raul C. Muresan, Weijia Feng, and Wolf Singer.\n
* For pairs of spikes and LFP signals in neuroscience, this measure is also known as "spike-triggered LFP average".\n
* Timestamps corresponding to events of one binary signal and one continuous signal must be provided. The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
*/

class CScaledCorrelationBC
{
	private:
		int		iScaleWindow;					///The size of the scale window
		int		iCorrelationWindow;				///The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		int		iTrialLength;					///The size of the trial in original sampling time units, same as the time stamps and samples in the input vectors
		float	*faScaledCorrelogram;			///The scale correlogram having a size of (2*iCorrelationWindow + 1), with lag 0 at faScaledCorrelogram[iCorrelationWindow]
		float	*faPBsCoefficientsSum;			///For each lag in the correlogram it summs up the Point Biserial coefficients that were used to compute the correlogram (faScaledCorrelogram[i] = faPBsCoefficientsSum[i] / iaPBsCoefficientsCount[i])
		int		*iaPBsCoefficientsCount;		///Counts how many Point Biserial coefficients have been added to each lag of the scale correlogram
		int		*iaPBsCoefficientsDistribution;	///Distribution of correlation coefficients that is also calculated with the cross-correlation measure; at index 0 we have the count for coefficients == -1, while at iPBsCorrCoeffDistributionBins we have the count for coefficients == +1
		int		*iaTimeStampsA;					///Internal pointer to the current source of binary signal provided in the call to ComputeScaledCorrelation
		float   *faSignalB;						///Internal pointer to the current source of continuous signal provided in the call to ComputeScaledCorrelation
		int		bUseFisherZTransform;			///Flag that specifies whether Fisher Z Transform should be applied to average the correlation coefficients
		int     iPBsCorrCoeffDistributionBins;	///Number of bins in the distribution of fi correlation coefficients; they span the interval of [-1..1]
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
		///Sets the correlogram to Not A Number (NAN) values to signal an error or undefined correlation
		void  SetCorrelogramToNAN();
		///Determines the Point Biserial coefficient of correlation a pair of segments;
		float ComputePBsCorrelation(int* piaBufferA,int piBufferASize,float* pfaBufferB,int piBufferBSize,int piWindowStartTime);					
		///The function that computes one pass scale correlation for a pair of trials
		void  ComputeTrialScaledCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB,int piTrialIndex);	

	public:
		//Construction interface
		
		///Constructor
		//Constructor parameters are as follows:
		//  piScaleWindow		- The scale on which you compute the scaled correlation; use the same units as the sampling unit of your time stamps
		//  piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
		//	piTrialLength		- The size of the trial in original sampling units
		//	piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetScaledCorrelogram() function will return NULL;
		//  Warning: the sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
		CScaledCorrelationBC(int piScaleWindow,int piCorrelationWindow,int piTrialLength,int &piIsError);
		
		///Destructor
		~CScaledCorrelationBC();

		//Functional interface
		
		///Compute the scaled correlation of one binary signal with one continuos signal
		//Parameters for ComputeScaledCorrelation
		//	piaTimeStampsA			- the vector holding the time stamps of the binary variable
		//  piNrTimeStampsInA		- number of time stamps in the first vector
		//  pfaSamplesB				- the vector holding samples of the continuous signal
		//  piNrSamplesInB			- number of samples in the second vector
		//  pbUseFisherZTransform	- set to 1 to use the Fisher Z transform to average the coefficients of correlation; set to 0 for normal computation;
		//  Warning: the sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
		void ComputeScaledCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbUseFisherZTransform);
		
		///Returns the buffer with the computed scaled correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
		//The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
		//Returns a NotANumber (NAN = -100000000.0) value for bins where correlation could not be defined because of lack of variance in the original signals
		float* GetScaledCrossCorrelogram();

		///Returns the sum of valid Point Biserial coefficients of correlation for each bin of the correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the sum of coefficients for lag 0 is at position iCorrelationWindow;
		//Use this function to get the sum of all valid coefficients that were used to compute the average correlation coefficient for each bin of the ScaledCrossCorrelogram: ScaledCrossCorrelogram[i] = PBsCoefficientSums[i] / PBsCoefficientCounts[i]
		float* GetPBsCoefficientSums();
		
		///Returns how many Point Biserial coefficients of correlation have been averaged for each bin of the scaled correlogram
		//The size of the buffer is (2*iCorrelationWindow+1), and the count of coefficients averaged for lag 0 is at position iCorrelationWindow;
		int* GetPBsCoefficientCounts();

		///Returns the distribution of coefficients of correlation
		//Outputs are: the buffer with the distribution of correlation coefficients, the number of bins in the distribution and the size of one bin
		int* GetDistributionOfCorrelationCoefficients(int& piNumberOfBins,float& pfBinSize);
		
		///Changes the scale segment size for the scaled correlation
		//Returns 0 on success, -1 if the scale window is too small; -2 if the new scale window is larger than the trial size;
		int ModifyScaleWindow(int piNewScale);
		
		///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter
		//Returns 0 on success, -1 if the window is too small; -2 if the window is larger than the length of the trial;
		int ModifyCorrelationWindow(int piNewCorrelationWindow);

		///Set the length of the trial
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