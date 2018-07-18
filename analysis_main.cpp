#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"

#define numScrMaps 10000

using namespace std;
//g++ analysis_main.cpp `root-config --cflags` `root-config --glibs` -o analysisMain

//sort array in increasing values
void sort_array(int numOfEvents, double array[]) 
{
    for(int j = 0; j < numOfEvents - 1; j++)
    {
        double currentMin = array[j];
        int currentMinIndex = j;

        for(int k = j+1; k < numOfEvents; k++)
        {
            if(currentMin > array[k])
            {
                currentMin = array[k];
                currentMinIndex = k;
            }
        }

        if(currentMinIndex != j)
        {
            array[currentMinIndex] = array[j];
            array[j] = currentMin;
        }
    }
}

//calculate time difference between consecutive events, given the number of consecutive events one wishes to look at (mutiplet)
double timeDiff_funct(const int numOfEvents, double timeStamp[], int firstEvent, int multiplet)
{

	int j = 1,k = 0;
	double deltaT = 0;
	while (j < multiplet)
	{
		deltaT += (multiplet - j)*(timeStamp[firstEvent + (multiplet - k - 1)] - timeStamp[firstEvent + k]);
		j += 2;
		k++;
	}
	return deltaT;
}

//calculate the test statistic of a given data set for a given multiplet value
void TS_funct(const int numOfEvents, int multiplet, double timeStamp[], double testStatistic[], int i)
{
	double timeDiff;

    int j = 0; 

    while (j <= numOfEvents - multiplet)
    {
        timeDiff = timeDiff_funct(numOfEvents,timeStamp,j,multiplet);

        if (timeDiff < testStatistic[i])
        {
            testStatistic[i] = timeDiff; // test statistic is defined as the smallest time difference 
        }

        j++;
    }
}

//calculate how many times two events are found within a 5 minutes time window
void Doublet_funct(const int numOfEvents, int multiplet, double timeStamp[], double counterNumOfDoublets_SC[], int i)
{
    int j = 0; 

    while (j <= numOfEvents - multiplet)
    {
        if (timeStamp[j+1] - timeStamp[j] <= 300)
        {
            counterNumOfDoublets_SC[i]++; // if 2 events are within 300 sec (5 min), counter is incremented by one
        }

        j++;
    }
}

//finds the smallest value in an array
double min_array(double array[], int numOfScrMaps)
{
    double min = array[0];

    for (int i = 0; i < numOfScrMaps; i++)
    {
        if (array[i] < min && array[i] != 0)
        {
            min = array[i];
        }
    }

    return min;
}

//finds the greatest value in an array
double max_array(double array[], int numOfScrMaps)
{
    double max = array[0];

    for (int i = 0; i < numOfScrMaps; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }

    return max;
}

/*fit the left tail of the test statistic distribution of scrambled maps to obtain
3 and 5 sigma values as well as the p-value and the sigma of test statistic of data*/
void fit_histo_TS(TH1D* histo, double testStatistic_data, double parameters_TS[])
{
    double startPoint = histo->GetBinLowEdge(2);   
    int maxBinLoc = histo->GetMaximumBin();
    double endPoint = histo->GetBinLowEdge(maxBinLoc-3);
    double binWidth = histo->GetBinWidth(maxBinLoc);

    double realStartPoint, realEndPoint;

    //--------- Chi-square minimization to obtain best fit--------------
    double minChiSquare = 100000, chiSquare;
    double tempEndPoint = startPoint+8*binWidth;

    while(tempEndPoint <= endPoint)
    {
        TF1 *fit1 = new TF1("fit1","expo",startPoint, endPoint);
        histo->Fit(fit1,"R Q");
        chiSquare = fit1->GetChisquare();

        if (chiSquare < minChiSquare)
        {
            realStartPoint = startPoint;
            realEndPoint = endPoint;
        }

        startPoint = startPoint + binWidth/2;
        tempEndPoint = startPoint + 8*binWidth; 
    }

    TF1 *fit1 = new TF1("fit1","expo",realStartPoint, realEndPoint);
    histo->Fit(fit1,"R Q");

    double expoConst = fit1->GetParameter(0);
    double expoSlope = fit1->GetParameter(1);

    //----------calculating integral on the right of the fit-------------
    int startBin = histo->FindBin(realEndPoint);
    int endBin = 100;
    double binRightTotalContent = 0;
    double rightIntegral;

    for (int i = startBin; i <= endBin; i++)
    {
        binRightTotalContent = binRightTotalContent + histo->GetBinContent(i);
    }

    rightIntegral = binRightTotalContent*binWidth;

    //------------- Matching the fit function constant with the bin content of startBin--------------
    double constant = log(histo->GetBinContent(startBin)) - expoSlope*realEndPoint;

    //--------calculating integral under fit-----------------------------
    double leftIntegral = (1./expoSlope)*(exp(expoSlope*realEndPoint+constant));

    //-------------normalizing total integral----------------
    double normLeftInt = leftIntegral/(leftIntegral+rightIntegral);
    double normRightInt = rightIntegral/(leftIntegral+rightIntegral);

    //----------- TS value to get 3 and 5 sigma ---------------
    double ts_3sig = (1./expoSlope)*(log(expoSlope*(leftIntegral+rightIntegral)*(1-0.998650102))-constant);
    double ts_5sig = (1./expoSlope)*(log(expoSlope*(leftIntegral+rightIntegral)*(1-0.999999713485))-constant);

    //--------------- Calculating p-value and sigma of data test statistic value -------------------
    double pvalue, dataIntegral, binDataContent = 0;
    int dataBin;

    if(log10(testStatistic_data) < realEndPoint)
    {
        pvalue = (1./expoSlope*exp(expoSlope*log10(testStatistic_data)+constant))/(leftIntegral+rightIntegral);
    }

    else
    {
        dataBin = histo->FindBin(log10(testStatistic_data));
        for (int i = dataBin; i <= endBin; i++)
        {
            binDataContent = binDataContent + histo->GetBinContent(i);
        }
        dataIntegral = binDataContent*binWidth;

        pvalue = 1.-dataIntegral/(leftIntegral+rightIntegral);
    }

    parameters_TS[0] = pvalue;
    parameters_TS[1] = ts_3sig;
    parameters_TS[2] = ts_5sig;
}

/*obtain parameters of doublet distribution of scrambled maps and calculate significance of number 
of doublets found in data*/
void fit_histo_doublet(TH1D* histo, double counterNumOfDoublets_data, double parameters_doublet[])
{    
    double mean = histo->GetMean();
    double std = histo->GetRMS();

    parameters_doublet[0] = (counterNumOfDoublets_data - mean)/std;
    parameters_doublet[1] = mean + 3*std;
    parameters_doublet[2] = mean + 5*std;
}
 

//int main(int multiplet, string user_name, int user_id, int currentDate)
int main(int nargv, char **argv)
{
	// extracting arguments
	int multiplet = atoi(argv[1]);
    string user_name = argv[2];
    int user_id = atoi(argv[3]);
    int currentDate = atoi(argv[4]);
        
    // Input file - timestamps
	ifstream input;
  	char data [256];
  	sprintf(data,"timestamp_%d.txt",user_id); //user id is defined in the bash script

  	cout << "Processing data from user " << user_name.c_str() << "..." << endl;

  	// Output file containing all the parameters needed for the plots (test statistic, number of doublets, sigma values, etc...)
  	FILE *pfile;
  	char output[256];
 	sprintf(output,"output_parameters_%s.txt",user_name.c_str()); //output file
 	pfile = fopen(output,"a");

 	string line;
 	int numOfEvents = 0;
  	double time_of_evt;
  	vector<double> timeStamp_Data_Vec_TS; //vector storing the timestamps used for the test statistic analysis
  	vector<double> timeStamp_Data_Vec_Doublet; //vector storing the timestamps used for the doublet analysis

 	input.open(data);
  	//Checking if the data file exists or/and can be opened
 	if (!input)
 	{
   		cout << "Unable to open file." << endl;
   		return 0;
 	}
 	
 	while (input >> time_of_evt)
 	{
 		timeStamp_Data_Vec_TS.push_back(time_of_evt);
    	timeStamp_Data_Vec_Doublet.push_back(time_of_evt);
   		
   		numOfEvents++; // Counting the number of events in data file
 	}
 	input.close();

 	cout << "Number of events in data (used for doublet analysis): " << numOfEvents << endl;

 	//removing all the events with the same timestamp for the test statistic analysis
 	for (int i = 1; i < timeStamp_Data_Vec_TS.size(); i++)
 	{
 		if (timeStamp_Data_Vec_TS[i] == timeStamp_Data_Vec_TS[i-1])
 		{
 			timeStamp_Data_Vec_TS.erase(timeStamp_Data_Vec_TS.begin() + i);
 			i--;
 		}
 	}

 	const int numOfEvents_TS = timeStamp_Data_Vec_TS.size();
 	const int numOfEvents_Doublet = numOfEvents;

 	cout << "Number of events used for the test statistic: " << numOfEvents_TS << endl;

 	// storing the timestamps in an array and coverting them to seconds
 	double timeStamp_Data_TS[numOfEvents_TS];
  	for(int i = 0; i < numOfEvents_TS; i++)
  	{
    	timeStamp_Data_TS[i] = timeStamp_Data_Vec_TS[i]/1000;
  	}

  	double timeStamp_Data_Doublet[numOfEvents_Doublet];
  	for(int i = 0; i < numOfEvents_Doublet; i++)
  	{
    	timeStamp_Data_Doublet[i] = timeStamp_Data_Vec_Doublet[i]/1000;
  	}

  	double firstTime_TS = timeStamp_Data_TS[0], lastTime_TS = timeStamp_Data_TS[numOfEvents_TS-1]; // saving first and last timestamp of data
	double timeSpan_TS = lastTime_TS - firstTime_TS; // time duration between first and last detection

  	double firstTime_Doublet = timeStamp_Data_Doublet[0], lastTime_Doublet = timeStamp_Data_Doublet[numOfEvents_Doublet-1]; // saving first and last timestamp of data
	double timeSpan_Doublet = lastTime_Doublet - firstTime_Doublet; // time duration between first and last detection
	
	// Scrambled maps (SM) variables
	TRandom *r1 = new TRandom3();
	double rand1; 
	double timeStamp_SC_TS[numOfEvents_TS], timeStamp_SC_Doublet[numOfEvents_Doublet]; // 2 dimensional array storing randomly generated timestamps for each SM
	double testStatistic_SC[numScrMaps]; //array storing the TS value for each SM
	double counterNumOfDoublets_SC[numScrMaps]; //array storing the number of doublets found in each SM

	//*********** Generating the scrambled maps in a 2-D array (numScrMaps by numOfEvents)**********
	for (int i = 0; i < numScrMaps; i++)
	{
		testStatistic_SC[i] = lastTime_TS;
		counterNumOfDoublets_SC[i] = 0;

		for (int j = 0; j < numOfEvents_Doublet; j++)
		{
			rand1 = r1->Rndm(j);

			if (j < numOfEvents_TS)
			{
				timeStamp_SC_TS[j] = firstTime_TS + rand1*timeSpan_TS;
			}

			timeStamp_SC_Doublet[j] = firstTime_Doublet + rand1*timeSpan_Doublet;		
		}

		// sorting arrays rows in increasing timestamps
		sort_array(numOfEvents_TS,timeStamp_SC_TS); 
		sort_array(numOfEvents_Doublet,timeStamp_SC_Doublet);

		// performing both TS and doublet analysis
		TS_funct(numOfEvents_TS,multiplet,timeStamp_SC_TS,testStatistic_SC,i); 
		Doublet_funct(numOfEvents_Doublet,multiplet,timeStamp_SC_Doublet,counterNumOfDoublets_SC,i);
				
	}	
	//**********************************************************************************************

	double maxBin_TS = max_array(testStatistic_SC,numScrMaps); // variable to store the highest TS for histogram purpose
	double minBin_TS = min_array(testStatistic_SC,numScrMaps); // variable to store the lowest TS for histogram purpose (arbitrary high enough value given)

	double maxBin_DOUBLET = max_array(counterNumOfDoublets_SC,numScrMaps); // variable to store the highest number of doublets for histogram purpose
	double minBin_DOUBLET = min_array(counterNumOfDoublets_SC,numScrMaps); // variable to store the lowest number of doublets for histogram purpose (arbitrary high enough value given)

	//******************Test statistics AND doublet histograms**************************************

	TCanvas *c1 = new TCanvas("c1","histogram_testStatistic");

	TString histo_file = TString::Format("TS_doublet_distributions_%d_user-%s.root",currentDate,user_name.c_str());
	TFile file1(histo_file, "RECREATE");

	// Test statistic
	TString histoTS = TString::Format("histogram_testStatistic");
	TString histoTS_title = TString::Format("Test statistic distribution for %d-plet - User %s",multiplet,user_name.c_str());
	TH1D *h1 = new TH1D(histoTS, histoTS_title, 100, log10(minBin_TS)-1, log10(maxBin_TS)+1);

	// Doublet
	TString histoDOUBLET = TString::Format("histogram_doublets");
	TString histoDOUBLET_title = TString::Format("Distribution of number of doublets in background - User %s",user_name.c_str());
	TH1D *h2 = new TH1D(histoDOUBLET, histoDOUBLET_title, 100, minBin_DOUBLET-1, maxBin_DOUBLET+1);

 	for (int i = 0; i < numScrMaps; i++)
	{
		h1->Fill(log10(testStatistic_SC[i]));
		h2->Fill(counterNumOfDoublets_SC[i]);
	}

	h1->GetXaxis()->SetTitle("log10(TS[sec])");
	h2->GetXaxis()->SetTitle("Number of doublets");

	h1->Write();
	h2->Write();
	//**********************************************************************************************

	//***********Unblinding data: calculating test statistic and number of doublets*****************
	double testStatistic_Data_array[1];
	double counterNumOfDoublets_Data_array[1];
	
	testStatistic_Data_array[0] = lastTime_TS;
	counterNumOfDoublets_Data_array[0] = 0;
	
	TS_funct(numOfEvents_TS,multiplet,timeStamp_Data_TS,testStatistic_Data_array,0);
	Doublet_funct(numOfEvents_Doublet,multiplet,timeStamp_Data_Doublet,counterNumOfDoublets_Data_array,0);
	//**********************************************************************************************

	double testStatistic_Data = testStatistic_Data_array[0];
	double counterNumOfDoublets_Data = counterNumOfDoublets_Data_array[0];

	//*******************Getting 3 and 5 sigma values for TS and Doublet histograms*****************
	double parameters_TS[3];
	double parameters_doublet[3];

	// fitting histograms to obtain sigma and p values
	fit_histo_TS(h1,testStatistic_Data,parameters_TS);
	fit_histo_doublet(h2,counterNumOfDoublets_Data,parameters_doublet);

	double pvalue_TS = parameters_TS[0], threeSig_TS = parameters_TS[1], fiveSig_TS = parameters_TS[2];
	double sig_TS = sqrt(2.)*TMath::ErfInverse(pvalue_TS);
	double meanTS_SC = h1->GetMean();

	double sig_doublet = parameters_doublet[0], threeSig_doublet = parameters_doublet[1], fiveSig_doublet = parameters_doublet[2];
	double pvalue_doublet = TMath::Erf(sig_doublet/sqrt(2));
	double meanNumOfDoublets_SC = h2->GetMean();

	cout << "*** User " << user_name.c_str() << " ***" << endl << endl;
	cout << "----- Test statistic analysis -----" << endl;
	cout << "log10(Expected test statistic [sec]) = " << meanTS_SC << " (" << pow(10,meanTS_SC) << " sec) || log10(Test statistic of data [sec]) = " << log10(testStatistic_Data) << " (" << testStatistic_Data << " sec) || pvalue = " << pvalue_TS << " || sigma = " << sig_TS << endl;
	cout << "3 sigma at " << threeSig_TS << " || 5 sigma at " << fiveSig_TS << endl << endl;

	cout << "----- Doublet analysis -----" << endl;
	cout << "Expected number of doublets = " << meanNumOfDoublets_SC << " || Number of doublets in data = " << counterNumOfDoublets_Data << " || pvalue = " << pvalue_doublet << " || sigma = " << sig_doublet << endl;
	cout << "3 sigma at " << threeSig_doublet << " || 5 sigma at " << fiveSig_doublet << endl;

	// saving all parameters needed for plots in output file
	fprintf(pfile, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",currentDate,meanTS_SC,log10(testStatistic_Data),threeSig_TS,fiveSig_TS,meanNumOfDoublets_SC,counterNumOfDoublets_Data,threeSig_doublet,fiveSig_doublet);
}
