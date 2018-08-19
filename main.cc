#include <MSTMixtureFitter.h>
#include <Logger.h>
#include <ConfigParser.h>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>
#include <TObjArray.h>
#include <TMatrixD.h>
#include <TColor.h>
#include <TApplication.h>
#include <TVector3.h>
#include <TView3D.h>
#include <TMarker.h>
#include <TPaletteAxis.h>
#include <TApplication.h>
#include <TGraph2D.h>

#include <RInside.h> // for the embedded R via RInside

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <getopt.h>


using namespace std;
using namespace MSTMixFitter_ns;


void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	cout<<"-c, --config=[CONFIG_FILENAME] \t Configuration file with options (NB: override all command line options if provided)"<<endl;
	cout<<"=============================="<<endl;
}

std::string GetStringLogLevel(int verbosity)
{
	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "config", required_argument, 0, 'c' },
	{(char*)0, (int)0, (int*)0, (int)0}
};


//Options
int verbosity= 4;//INFO level
std::string configFileName= "";

//Functions
//...

int main(int argc, char **argv)
{
	//====================================================
	//==         PARSE ARGS
	//=====================================================
	//## Check args
	if(argc<2){
		cout<<endl;
		cerr<< "ERROR: Incorrect number of arguments...see program usage!"<<endl;
		Usage(argv[0]);		
		exit(1);
	}

	//## Get command args
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hv::c:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;
			}
    	default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
	
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	//====================================================
	//==         PARSE CONFIG FILE
	//=====================================================
	ConfigParser parser;
	if(parser.ReadConfig(configFileName)<0){
		ERROR_LOG("Failed to read and parse config file "<<configFileName<<"!");
		return -1;
	}

	//====================================================
	//==         RUN FITTER
	//=====================================================
	INFO_LOG("Running fitter...");
	MSTMixtureFitter* fitter= new MSTMixtureFitter;
	if(fitter->Run()<0){
		cerr<<"ERROR: Failed to compute imputed data!"<<endl;
		return -1;
	}

	//Clear data
	INFO_LOG("Clearing data...");
	if(fitter){
		delete fitter;
		fitter= 0;
	}

	return 0;
	
}//close macro




