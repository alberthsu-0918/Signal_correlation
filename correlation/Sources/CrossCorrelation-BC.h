/*
	++ Source:
			CrossCorrelation-BC.h

	++ Description:
			A class that computes the cross-correlation for a vector of time stamps and a continuous signal (Binary - Continuous)
			Cross-Correlation is computed in the classical way by counting
			For pairs of spikes and LFP signals in neuroscience, this measure is also known as "spike-triggered LFP average"

	++ History:
			22-03-2007 - Raul C. Muresan : Created file and added main declarations

	++ Disclaimer:
			The code is free for non commercial purposes. You may use it freely, with the sole restriction that you may not claim that you wrote it. 
			I do not warrant that the code is 100% bug free. Use at your own risk! 
*/

//Doxygen comment
/**
* \brief
* A class that computes the cross-correlation for a vector of time stamps and a continuous signal (Binary - Continuous).
* \details
* Cross-correlation is computed in the classical way.\n
* For pairs of spikes and LFP signals in neuroscience, this measure is also known as "spike-triggered LFP average".\n
* Timestamps corresponding to events of one binary signal and one continuous signal must be provided. The time stamps of the binary signal must match the samples of the continuous signal, i.e. the first sample of the continuous signal corresponds to time stamp 0.
*/


#ifndef USED_CROSSCORRELATION_BC
#define USED_CROSSCORRELATION_BC

class CCrossCorrelationBC
{
	private:
		int		iCorrelationWindow;				///The size of the correlation window (for example to compute cross-correlation in a window of -100..100, iCorrelationWindow = 100;
		int		iTrialLength;					///The size of the trial in original sampling time units, same as the time stamps and samples in the input vectors
		float	*faCrossCorrelogram;			///The scale correlogram having a size of (2*iCorrelationWindow + 1), with lag 0 at faCrossCorrelogram[iCorrelationWindow]
		int		iNumberOfSummedCorrelograms;	///The total number of correlograms that were summed up when computing the cross correlogram; used to normalize the cross correlogram;
		int		*iaTimeStampsA;					///Internal pointer to the current source of binary signal provided in the call to ComputeCrossCorrelation
		float   *faSignalB;						///Internal pointer to the current source of continuous signal provided in the call to ComputeCrossCorrelation
		int		*iaNrSummedBins;				///An array that holds the number of correlations that have been summed for each correlation bin

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
		///The function that computes one pass scale correlation for a pair of trials
		void  ComputeTrialCrossCorrelation(int piIdxTrialStartA,int piIdxTrialEndA,int piIdxTrialStartB,int piIdxTrialEndB);	

	public:
		//Construction interface
		
		///Constructor
		//Constructor parameters are as follows:
		//  piCorrelationWindow	- The size of the cross-correlation window (eg. 80 for a cross-correlation with lags of -80..+80); use again the same units as the sampling of your time stamps
		//	piTrialLength		- The size of the trial in original sampling units
		//	piIsError			- A return parameter, set to 0 on success. It is set to -1 by the constructor if there is a mismatch of window sizes and the GetCrossCorrelogram() function will return NULL;
		//  Warning: the sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
		CCrossCorrelationBC(int piCorrelationWindow,int piTrialLength,int &piIsError);
		
		///Destructor
		~CCrossCorrelationBC();

		//Functional interface
		
		///Compute the cross-correlation of one binary signal with one continuos signal
		//Parameters for ComputeCrossCorrelation
		//	piaTimeStampsA			- the vector holding the time stamps of the binary variable
		//  piNrTimeStampsInA		- number of time stamps in the first vector
		//  pfaSamplesB				- the vector holding samples of the continuous signal
		//  piNrSamplesInB			- number of samples in the second vector
		//	pbNormalizeCorrelogram	- set to 1 to normalize the correlogram
		//  Warning: the sampling frequency of the continuous and binary signals has to be the same! Also the trial lengths have to match!
		void ComputeCrossCorrelation(int *piaTimeStampsA,int piNrTimeStampsInA,float *pfaSamplesB,int piNrSamplesInB,int pbNormalizeCorrelogram);
		
		///Returns the buffer with the computed cross-correlogram; returns NULL if there was an error; please check the buffer for NULL before using it!!
		//The size of the cross correlogram buffer is (2*iCorrelationWindow+1), and element with lag 0 is at position iCorrelationWindow;
		//Returns NAN (Not A Number) = -100000000.0 if one of the buffers is empty and thus the correlation is not defined; please check in your code for this!!!!
		float* GetCrossCorrelogram();	
		
		///Set the size of the correlation window; the parameter specifies the correlation window; for example for a window of -100..+100 pass 100 as a parameter; 
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