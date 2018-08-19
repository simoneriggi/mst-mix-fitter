/**
* @file ConfigParser.cc
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#include <ConfigParser.h>
#include <Logger.h>

#include <TMath.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

namespace MSTMixFitter_ns {


//### FILE INFO ###
std::string ConfigParser::fConfigFileName;
std::string ConfigParser::fInputFileName;
std::string ConfigParser::fOutputFileName;

//### DATA READ OPTIONS
std::string ConfigParser::fDataReadDelimiter;

//## ALGORITHM OPTIONS
int ConfigParser::fNDim;
int ConfigParser::fNComponents;
int ConfigParser::fNIterations;
bool ConfigParser::fUseStoppingCriteria;
double ConfigParser::fEpsilon;
bool ConfigParser::fRunMinuitFitterAtConvergence;
int ConfigParser::fFitterChoice;	
		

//## STARTING VALUES OPTIONS
int ConfigParser::fParInitMethod;
bool ConfigParser::fRandomizeStartPars;
bool ConfigParser::fRandomizeStartCovariancePars;
bool ConfigParser::fRandomizeStartMeanPars;
std::vector<TMatrixD> ConfigParser::fMeanStartPars;
std::vector<TMatrixD> ConfigParser::fDeltaStartPars;
std::vector<TMatrixD> ConfigParser::fCovarianceStartPars;
std::vector<double> ConfigParser::fFractionStartPars;
std::vector<double> ConfigParser::fNuStartPars;

//## CONSTRAINT OPTIONS
bool ConfigParser::fFixMeanPars;
bool ConfigParser::fFixCovariancePars;
bool ConfigParser::fFixFractionPars;
bool ConfigParser::fFixDeltaPars;
bool ConfigParser::fFixNuPars;
bool ConfigParser::fForceDiagonalCovariance;
bool ConfigParser::fUseRandomRegenerationAfterStuck;
bool ConfigParser::fUseConstraints;
double ConfigParser::fConstraintAlphaScale;
double ConfigParser::fConstraintAlphaTolerance;
bool ConfigParser::fUseCovarianceEigenBoundConstraint;
std::vector<TMatrixD> ConfigParser::fCovarianceEigenMinBound;
std::vector<TMatrixD> ConfigParser::fCovarianceEigenMaxBound;
bool ConfigParser::fUseMeanBoundConstraint;
std::vector<TMatrixD> ConfigParser::fMeanMinBound;
std::vector<TMatrixD> ConfigParser::fMeanMaxBound;
bool ConfigParser::fUseCovarianceBoundConstraint;
std::vector<TMatrixD> ConfigParser::fCovarianceMinBound;
std::vector<TMatrixD> ConfigParser::fCovarianceMaxBound;
bool ConfigParser::fUseDeltaBoundConstraint;
std::vector<TMatrixD> ConfigParser::fDeltaMinBound;
std::vector<TMatrixD> ConfigParser::fDeltaMaxBound;
bool ConfigParser::fUseNuBoundConstraint;
std::vector<double> ConfigParser::fNuMinBound;
std::vector<double> ConfigParser::fNuMaxBound;



ConfigParser::ConfigParser()
{
	//### FILE INFO ###
	fConfigFileName= "";
	fInputFileName= "";
	fOutputFileName= "FitOutput.root";
	//fConstraintDataFileName= "";
	//fIsInteractiveRun= false;
	
	fDataReadDelimiter= "";
	
	//### MAIN OPTIONS ###
	fNDim= 2;
	fNComponents= 3;
	fNIterations= 10;
	fUseStoppingCriteria= false;
	fEpsilon= 1.e-6;
	fRunMinuitFitterAtConvergence= false;
	fFitterChoice= 1;

	
	//### START VALUES OPTIONS
	fParInitMethod= 1;
	fRandomizeStartPars= false;
	fRandomizeStartCovariancePars= false;
	fRandomizeStartMeanPars= false;
	fMeanStartPars.clear();
	fDeltaStartPars.clear();
	fCovarianceStartPars.clear();
	fFractionStartPars.clear();
	fNuStartPars.clear();
	
	//### CONSTRAINT OPTIONS
	fFixMeanPars= false;
	fFixCovariancePars= false;
	fFixFractionPars= false;
	fFixDeltaPars= false;
	fFixNuPars= false;
	fForceDiagonalCovariance= false;
	fUseRandomRegenerationAfterStuck= true;
	fUseConstraints= true;
	fConstraintAlphaScale= 2;
	fConstraintAlphaTolerance= 1.e-4;
	fUseCovarianceEigenBoundConstraint= true;
	fCovarianceEigenMinBound.clear();
	fCovarianceEigenMaxBound.clear();
	fUseMeanBoundConstraint= false;
	fMeanMinBound.clear();
	fMeanMaxBound.clear();
	fUseCovarianceBoundConstraint= false;
	fCovarianceMinBound.clear();
	fCovarianceMaxBound.clear();
	fUseDeltaBoundConstraint= false;
	fDeltaMinBound.clear();
	fDeltaMaxBound.clear();
	fUseNuBoundConstraint= false;
	fNuMinBound.clear();
	fNuMaxBound.clear();

}//close costructor

ConfigParser::~ConfigParser(){

}//close destructor


int ConfigParser::ReadConfig(std::string filename) 
{
	//Set config filename
	fConfigFileName= filename;

	//## Read configuration from file
	INFO_LOG("Reading and parsing file "<<fConfigFileName.c_str()<<" ...");
	
  ifstream in;  
  in.open(fConfigFileName.c_str());
  if(!in.good()) {
    ERROR_LOG("Cannot read config file "<< fConfigFileName.c_str()<<"!");
    return -1;
  }

	//Start parsing the config file
	char buffer[1000];//container for a full line
  std::string descriptor;//container for config descriptor 
  in.getline(buffer,1000);//get the full line
	std::string parsedline;
	
	while(std::getline(in,parsedline)) {
		
		char first_char= *(parsedline.c_str());
		
		if(first_char!='#' && first_char!='\n' && first_char!=' '){
			stringstream line(parsedline);
			stringstream line_copy(parsedline);
     
			line_copy >> descriptor;
      
      if(descriptor!="\n" && descriptor!=""){
				//get all config parameters

				//########################
				//####  RUN CONFIG
				//########################	
				if(descriptor.compare("inputFile")==0){
					line >> descriptor >> fInputFileName;			
		  	}
				else if(descriptor.compare("outputFile")==0){
					line >> descriptor >> fOutputFileName;			
		  	}
				else if(descriptor.compare("dataDelimiter")==0){
					line >> descriptor >> fDataReadDelimiter;			
					cout<<"fDataReadDelimiter="<<fDataReadDelimiter<<" (l="<<fDataReadDelimiter.size()<<")"<<endl;
		  	}
				
				else if(descriptor.compare("nDim")==0){
					line >> descriptor >> fNDim;
				}
				else if(descriptor.compare("nComponents")==0){
					line >> descriptor >> fNComponents;
				}
				else if(descriptor.compare("nIterations")==0){
					line >> descriptor >> fNIterations;
				}
				else if(descriptor.compare("useStoppingCriteria")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >>fEpsilon;			
					if(thisFlagValue.compare("T")==0) fUseStoppingCriteria= true;
					else if(thisFlagValue.compare("F")==0) fUseStoppingCriteria= false;
					else{
						ERROR_LOG("Invalid setting for useStoppingCriteria, use T or F!");
    				return -1;
					}
		  	}//close else if

				else if(descriptor.compare("fitter")==0){
					line >> descriptor >> fFitterChoice;
				}
				else if(descriptor.compare("runMinuitAtConvergence")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fRunMinuitFitterAtConvergence= true;
					else if(thisFlagValue.compare("F")==0) fRunMinuitFitterAtConvergence= false;
					else{
						ERROR_LOG("Invalid setting for runMinuitAtConvergence, use T or F!");
    				return -1;
					}
		  	}//close else if
				
				//====  START PAR OPTIONS =====
				else if(descriptor.compare("parInitMethod")==0){
					line >> descriptor >> fParInitMethod;
				}
				else if(descriptor.compare("randomizeStartPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fRandomizeStartPars= true;
					else if(thisFlagValue.compare("F")==0) fRandomizeStartPars= false;
					else{
						ERROR_LOG("Invalid setting for randomizeStartPars, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("randomizeStartCovariancePars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fRandomizeStartCovariancePars= true;
					else if(thisFlagValue.compare("F")==0) fRandomizeStartCovariancePars= false;
					else{
						ERROR_LOG("Invalid setting for randomizeStartCovariancePars, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("randomizeStartMeanPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fRandomizeStartMeanPars= true;
					else if(thisFlagValue.compare("F")==0) fRandomizeStartMeanPars= false;
					else{
						ERROR_LOG("Invalid setting for randomizeStartMeanPars, use T or F!");
    				return -1;
					}
		  	}//close else if
				
				else if(descriptor.compare("meanStartPars")==0){
					TMatrixD meanPar(1,fNDim);
					meanPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						meanPar(0,j)= thisEntry;
					}
					fMeanStartPars.push_back(meanPar);	
				}
				else if(descriptor.compare("deltaStartPars")==0){
					TMatrixD deltaPar(1,fNDim);
					deltaPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						deltaPar(0,j)= thisEntry;
					}
					fDeltaStartPars.push_back(deltaPar);	
				}
				else if(descriptor.compare("covarianceStartPars")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					fCovarianceStartPars.push_back(covariancePar);	
				}
				else if(descriptor.compare("fractionStartPars")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					fFractionStartPars.push_back(thisEntry);
				}	
				else if(descriptor.compare("nuStartPars")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					fNuStartPars.push_back(thisEntry);
				}


				//====  CONSTRAINTS OPTIONS ====
				else if(descriptor.compare("fixMeanPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFixMeanPars= true;
					else if(thisFlagValue.compare("F")==0) fFixMeanPars= false;
					else{
						ERROR_LOG("Invalid setting for fixMeanPars, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("fixCovariancePars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFixCovariancePars= true;
					else if(thisFlagValue.compare("F")==0) fFixCovariancePars= false;
					else{
						ERROR_LOG("Invalid setting for fixCovariancePar, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("fixDeltaPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFixDeltaPars= true;
					else if(thisFlagValue.compare("F")==0) fFixDeltaPars= false;
					else{
						ERROR_LOG("Invalid setting for fixDeltaPars, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("fixNuPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFixNuPars= true;
					else if(thisFlagValue.compare("F")==0) fFixNuPars= false;
					else{
						ERROR_LOG("Invalid setting for fixNuPars, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("fixFractionPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFixFractionPars= true;
					else if(thisFlagValue.compare("F")==0) fFixFractionPars= false;
					else{
						ERROR_LOG("Invalid setting for fixFractionPar, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("forceDiagonalCovariance")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fForceDiagonalCovariance= true;
					else if(thisFlagValue.compare("F")==0) fForceDiagonalCovariance= false;
					else{
						ERROR_LOG("Invalid setting for forceDiagonalCovariance, use T or F!");
    				return -1;
					}
		  	}//close else if

				
				else if(descriptor.compare("useConstraints")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fConstraintAlphaScale >> fConstraintAlphaTolerance;			
					if(thisFlagValue.compare("T")==0) fUseConstraints= true;
					else if(thisFlagValue.compare("F")==0) fUseConstraints= false;
					else{
						ERROR_LOG("Invalid setting for useConstraint, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("useRandomRegenerationAfterStuck")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomRegenerationAfterStuck= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomRegenerationAfterStuck= false;
					else{
						ERROR_LOG("Invalid setting for useRandomRegenerationAfterStuck, use T or F!");
    				return -1;
					}
		  	}//close else if
				else if(descriptor.compare("useCovarianceEigenBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseCovarianceEigenBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseCovarianceEigenBoundConstraint= false;
					else{
						ERROR_LOG("Invalid setting for useCovarianceEigenBoundConstraint, use T or F!");
    				return -1;
					}
		  	}//close else if

				else if(descriptor.compare("covarianceEigenMinBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fCovarianceEigenMinBound.push_back(bound);	
				}//close else if
				else if(descriptor.compare("covarianceEigenMaxBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fCovarianceEigenMaxBound.push_back(bound);			
				}//close else if

				else if(descriptor.compare("useCovarianceBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseCovarianceBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseCovarianceBoundConstraint= false;
					else{
						ERROR_LOG("Invalid setting for useCovarianceBoundConstraint, use T or F!");
    				return -1;
					}
		  	}//close else if
				
				else if(descriptor.compare("covarianceMinBound")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					fCovarianceMinBound.push_back(covariancePar);	
				}//close else if
				else if(descriptor.compare("covarianceMaxBound")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					fCovarianceMaxBound.push_back(covariancePar);	
				}//close else if

				

				else if(descriptor.compare("useMeanBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseMeanBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseMeanBoundConstraint= false;
					else{
						ERROR_LOG("Invalid setting for useMeanBoundConstraint, use T or F!");
						return -1;
					}
		  	}//close else if

				else if(descriptor.compare("meanMinBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fMeanMinBound.push_back(bound);	
				}//close else if
				else if(descriptor.compare("meanMaxBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fMeanMaxBound.push_back(bound);			
				}//close else if


				else if(descriptor.compare("useDeltaBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseDeltaBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseDeltaBoundConstraint= false;
					else{
						ERROR_LOG("Invalid setting for useDeltaBoundConstraint, use T or F!");
						return -1;
					}
		  	}//close else if

				else if(descriptor.compare("deltaMinBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fDeltaMinBound.push_back(bound);	
				}//close else if
				else if(descriptor.compare("deltaMaxBound")==0){
					TMatrixD bound(1,fNDim);
					bound.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						bound(0,j)= thisEntry;	
					}
					fDeltaMaxBound.push_back(bound);			
				}//close else if


				else if(descriptor.compare("useNuBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseNuBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseNuBoundConstraint= false;
					else{
						ERROR_LOG("Invalid setting for useNuBoundConstraint, use T or F!");
						return -1;
					}
		  	}//close else if

				else if(descriptor.compare("nuMinBound")==0){
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						fNuMinBound.push_back(thisEntry);	
					}
				}//close else if
				else if(descriptor.compare("nuMaxBound")==0){
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						fNuMaxBound.push_back(thisEntry);	
					}
				}//close else if

				/*
				
				else if(descriptor.compare("useRandomStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomFromModelStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomFromModelStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomFromModelStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomFromModelStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomMeanStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomMeanStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomMeanStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomMeanStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomMeanDiffStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomMeanDiffStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomMeanDiffStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomMeanDiffStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomCovarianceStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomCovarianceStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomCovarianceStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomCovarianceStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if	
				else if(descriptor.compare("useRandomCovarianceEigenStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomCovarianceEigenStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomCovarianceEigenStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomCovarianceEigenStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if	
				else if(descriptor.compare("useRandomDeltaStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomDeltaStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomDeltaStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomDeltaStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomNuStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomNuStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomNuStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomNuStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useRandomFractionStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseRandomFractionStart= true;
					else if(thisFlagValue.compare("F")==0) fUseRandomFractionStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useRandomFractionStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useMeanImputationStart")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseMeanImputationStart= true;
					else if(thisFlagValue.compare("F")==0) fUseMeanImputationStart= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useMeanImputationStart, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				

				
				
				

				
				

				else if(descriptor.compare("useMeanConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseMeanConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseMeanConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useMeanConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				
				else if(descriptor.compare("useMeanDiffConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseMeanDiffConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseMeanDiffConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useMeanDiffConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if

				else if(descriptor.compare("useLocationConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseLocationConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseLocationConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useLocationConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useLocationDiffConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseLocationDiffConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseLocationDiffConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useMeanDiffConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if




				else if(descriptor.compare("useCovarianceConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseCovarianceConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseCovarianceConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useCovarianceConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				

				else if(descriptor.compare("useCovarianceEigenConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseCovarianceEigenConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseCovarianceEigenConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useCovarianceEigenConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				


		
				else if(descriptor.compare("useScaleMatrixConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseScaleMatrixConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseScaleMatrixConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useScaleMatrixConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useScaleMatrixBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseScaleMatrixBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseScaleMatrixBoundConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useScaleMatrixBoundConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if

				else if(descriptor.compare("useScaleMatrixEigenConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseScaleMatrixEigenConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseScaleMatrixEigenConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useScaleMatrixEigenConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useScaleMatrixEigenBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseScaleMatrixEigenBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseScaleMatrixEigenBoundConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useScaleMatrixEigenBoundConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if

				else if(descriptor.compare("useDeltaBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseDeltaBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseDeltaBoundConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useDeltaBoundConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("useNuBoundConstraint")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseNuBoundConstraint= true;
					else if(thisFlagValue.compare("F")==0) fUseNuBoundConstraint= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useNuBoundConstraint, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if

				else if(descriptor.compare("useAnnealing")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fAnnealingParStart >> fAnnealingParStep;			
					if(thisFlagValue.compare("T")==0) fUseAnnealing= true;
					else if(thisFlagValue.compare("F")==0) fUseAnnealing= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useAnnealing, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if


				else if(descriptor.compare("meanConstraintSign")==0){
					fMeanConstraintSign= new TMatrixD(fNDim,1);
					fMeanConstraintSign->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						std::string constrainSignSymbol;
						double constraintSign;
		    		line >> constrainSignSymbol;
						if(constrainSignSymbol.compare("<")==0) constraintSign= +1;
						else if(constrainSignSymbol.compare(">")==0) constraintSign= -1;
						else{
							string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for meanConstraintSign, use < or >!";
    					throw std::runtime_error(errMsg);
							exit(1);
						}
							
						(*fMeanConstraintSign)(j,0)= constraintSign;	
					}
				}//close else if

				else if(descriptor.compare("sigmaConstraintSign")==0){
					fSigmaConstraintSign= new TMatrixD(fNDim,1);
					fSigmaConstraintSign->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						std::string constrainSignSymbol;
						double constraintSign;
		    		line >> constrainSignSymbol;
						if(constrainSignSymbol.compare("<")==0) constraintSign= +1;
						else if(constrainSignSymbol.compare(">")==0) constraintSign= -1;
						else{
							string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for sigmaConstraintSign, use < or >!";
    					throw std::runtime_error(errMsg);
							exit(1);
						}
							
						(*fSigmaConstraintSign)(j,0)= constraintSign;	
					}
				}//close else if


				else if(descriptor.compare("meanToleranceMin")==0){
					fMeanParTolerance_min= new TMatrixD(fNDim,1);
					fMeanParTolerance_min->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fMeanParTolerance_min)(j,0)= thisEntry;	
					}
				}//close else if
				else if(descriptor.compare("meanToleranceMax")==0){
					fMeanParTolerance_max= new TMatrixD(fNDim,1);
					fMeanParTolerance_max->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fMeanParTolerance_max)(j,0)= thisEntry;	
					}
				}//close else if

				else if(descriptor.compare("meanDiffToleranceMin")==0){
					fMeanDiffParTolerance_min= new TMatrixD(fNDim,1);
					fMeanDiffParTolerance_min->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fMeanDiffParTolerance_min)(j,0)= thisEntry;	
					}
				}//close else if
				else if(descriptor.compare("meanDiffToleranceMax")==0){
					fMeanDiffParTolerance_max= new TMatrixD(fNDim,1);
					fMeanDiffParTolerance_max->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fMeanDiffParTolerance_max)(j,0)= thisEntry;	
					}
				}//close else if
	


				else if(descriptor.compare("covarianceToleranceMin")==0){
					fSigmaParTolerance_min= new TMatrixD(fNDim,fNDim);
					fSigmaParTolerance_min->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							(*fSigmaParTolerance_min)(j,l)= thisEntry;	
						}
					}
				}//close else if
				else if(descriptor.compare("covarianceToleranceMax")==0){
					fSigmaParTolerance_max= new TMatrixD(fNDim,fNDim);
					fSigmaParTolerance_max->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							(*fSigmaParTolerance_max)(j,l)= thisEntry;	
						}
					}
				}//close else if
		
				else if(descriptor.compare("deltaToleranceMin")==0){
					fDeltaParTolerance_min= new TMatrixD(fNDim,1);
					fDeltaParTolerance_min->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fDeltaParTolerance_min)(j,0)= thisEntry;	
					}
				}//close else if
				else if(descriptor.compare("deltaToleranceMax")==0){
					fDeltaParTolerance_max= new TMatrixD(fNDim,1);
					fDeltaParTolerance_max->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						(*fDeltaParTolerance_max)(j,0)= thisEntry;	
					}
				}//close else if


				else if(descriptor.compare("dataCutMin")==0){
					fDataCutMin= new TMatrixD(fNDim,1);
					fDataCutMin->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						std::string cutSymbol;
						double cutValue;
		    		line >> cutSymbol;
						if(cutSymbol.compare("+inf")==0) cutValue= +TMath::Infinity();
						else if(cutSymbol.compare("-inf")==0) cutValue= -TMath::Infinity();
						else{
							cutValue=	atof(cutSymbol.c_str());
						}
							
						(*fDataCutMin)(j,0)= cutValue;	
					}
				}//close else if
				else if(descriptor.compare("dataCutMax")==0){
					fDataCutMax= new TMatrixD(fNDim,1);
					fDataCutMax->Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						std::string cutSymbol;
						double cutValue;
		    		line >> cutSymbol;
						if(cutSymbol.compare("inf")==0) cutValue= +TMath::Infinity();
						else if(cutSymbol.compare("-inf")==0) cutValue= -TMath::Infinity();
						else{
							cutValue=	atof(cutSymbol.c_str());
						}
							
						(*fDataCutMax)(j,0)= cutValue;	
					}
				}//close else if

				else if(descriptor.compare("drawRange")==0){
					double minvalue, maxvalue, nbins;
					line >> descriptor >> minvalue >> maxvalue >> nbins;	
					
					fDrawMinRange.push_back(minvalue);
					fDrawMaxRange.push_back(maxvalue);
					fDrawNBins.push_back(nbins);

				}//close else if
				

				else if(descriptor.compare("meanStart")==0){
					TMatrixD meanPar(fNDim,1);
					meanPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		 line >> thisEntry;	
						meanPar(j,0)= thisEntry;	
					}
					fMeanStartPar.push_back(meanPar);	
				}//close else if
				else if(descriptor.compare("covarianceStart")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					//cout<<"*** C ***"<<endl;
					//covariancePar.Print();

					fCovarianceStartPar.push_back(covariancePar);	
				}//close else if
				else if(descriptor.compare("deltaStart")==0){
					TMatrixD deltaPar(fNDim,1);
					deltaPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						deltaPar(j,0)= thisEntry;	
					}
				
					fDeltaStartPar.push_back(deltaPar);	
				}//close else if

				else if(descriptor.compare("nuStart")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					
					fNuStartPar.push_back(thisEntry);

				}//close else if
				else if(descriptor.compare("fractionStart")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					
					fFractionStartPar.push_back(thisEntry);

				}//close else if
		
			

				else if(descriptor.compare("meanTrue")==0){
					TMatrixD meanPar(fNDim,1);
					meanPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						meanPar(j,0)= thisEntry;	
					}
					fMeanTruePar.push_back(meanPar);	
				}//close else if
				else if(descriptor.compare("covarianceTrue")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					fCovarianceTruePar.push_back(covariancePar);	
				}//close else if
				else if(descriptor.compare("deltaTrue")==0){
					TMatrixD deltaPar(fNDim,1);
					deltaPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						deltaPar(j,0)= thisEntry;	
					}
				
					fDeltaTruePar.push_back(deltaPar);	
				}//close else if

				else if(descriptor.compare("nuTrue")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					
					fNuTruePar.push_back(thisEntry);

				}//close else if
				else if(descriptor.compare("fractionTrue")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					
					fFractionTruePar.push_back(thisEntry);

				}//close else if


				else if(descriptor.compare("useMeanOffset")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue ;			
					if(thisFlagValue.compare("T")==0) fUseMeanOffset= true;
					else if(thisFlagValue.compare("F")==0) fUseMeanOffset= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useMeanOffset, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("meanOffset")==0){
					TMatrixD meanOffsetPar(fNDim,1);
					meanOffsetPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						meanOffsetPar(j,0)= thisEntry;	
					}
					fMeanOffsetPar.push_back(meanOffsetPar);	
				}//close else if


				else if(descriptor.compare("useDeltaOffset")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue ;			
					if(thisFlagValue.compare("T")==0) fUseDeltaOffset= true;
					else if(thisFlagValue.compare("F")==0) fUseDeltaOffset= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useDeltaOffset, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("deltaOffset")==0){
					TMatrixD deltaOffsetPar(fNDim,1);
					deltaOffsetPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						deltaOffsetPar(j,0)= thisEntry;	
					}
					fDeltaOffsetPar.push_back(deltaOffsetPar);	
				}//close else if



				else if(descriptor.compare("useSigmaOffset")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue ;			
					if(thisFlagValue.compare("T")==0) fUseSigmaOffset= true;
					else if(thisFlagValue.compare("F")==0) fUseSigmaOffset= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useSigmaOffset, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("sigmaOffset")==0){
					TMatrixD sigmaOffsetPar(fNDim,1);
					sigmaOffsetPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						sigmaOffsetPar(j,0)= thisEntry;	
					}
					fSigmaOffsetPar.push_back(sigmaOffsetPar);	
				}//close else if

	



				else if(descriptor.compare("meanModel")==0){
					TMatrixD meanPar(fNDim,1);
					meanPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						meanPar(j,0)= thisEntry;	
					}
					fMeanModelPar.push_back(meanPar);	
				}//close else if
				else if(descriptor.compare("covarianceModel")==0){
					TMatrixD covariancePar(fNDim,fNDim);
					covariancePar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						for(int l=0;l<fNDim;l++){
							double thisEntry;
		    		 	line >> thisEntry;	
							covariancePar(j,l)= thisEntry;	
						}
					}
					fCovarianceModelPar.push_back(covariancePar);	
				}//close else if
				else if(descriptor.compare("deltaModel")==0){
					TMatrixD deltaPar(fNDim,1);
					deltaPar.Zero();
			
					line >> descriptor;	
					for(int j=0;j<fNDim;j++){
						double thisEntry;
		    		line >> thisEntry;	
						deltaPar(j,0)= thisEntry;	
					}
				
					fDeltaModelPar.push_back(deltaPar);	
				}//close else if

				else if(descriptor.compare("nuModel")==0){
					double thisEntry;
					line >> descriptor >> thisEntry;
					
					fNuModelPar.push_back(thisEntry);

				}//close else if

				
				
				*/
				else{//config setting not defined
					ERROR_LOG("Descriptor " << descriptor.c_str()<<" not defined, bad config settings!");
					return -1;
				}

			}//close if descriptor
		}//close if buffer

		if (!in.good()) break;

	}//close while

	in.close();
	
	//## Print parsed info
	Print();

	return 0;

}//close ReadConfig()



void ConfigParser::Print()
{
	cout<<"################################"<<endl;
	cout<<"###     PARSED  SETTINGS   #####"<<endl;
	cout<<"################################"<<endl;
	cout<<"== MAIN OPTIONS =="<<endl;
	cout<<"InputFileName: "<<fInputFileName<<endl;
	cout<<"OutputFileName: "<<fOutputFileName<<endl;
	cout<<"NDim: "<<fNDim<<endl;
	cout<<"NComponents: "<<fNComponents<<endl;
	cout<<"NIterations: "<<fNIterations<<endl;
	cout<<"UseStoppingCriteria? "<<fUseStoppingCriteria<<"  Epsilon: "<<fEpsilon<<endl;
	cout<<"== START PAR OPTIONS =="<<endl;
	cout<<"ParInitMethod: "<<fParInitMethod<<endl;
	cout<<"RandomizeStartPars? "<<fRandomizeStartPars<<endl;
	for(size_t k=0;k<fFractionStartPars.size();k++) cout<<"fraction="<<fFractionStartPars[k]<<endl;
	cout<<"== CONSTRAINT PAR OPTIONS =="<<endl;
	cout<<"FixMeanPars? "<<fFixMeanPars<<endl;
	cout<<"FixCovariancePars? "<<fFixCovariancePars<<endl;
	cout<<"FixFractionPars? "<<fFixFractionPars<<endl;
	cout<<"ForceDiagonalCovariance? "<<fForceDiagonalCovariance<<endl;	
	cout<<"################################"<<endl;

}//close Print()


}//close namespace
