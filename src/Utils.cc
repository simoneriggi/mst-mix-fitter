/**
* @file Utils.cc
* @class Utils
* @brief Utility functions 
*
* Utility functions
* @author S. Riggi
* @date 17/09/2012
*/

#include <Utils.h>
#include <Logger.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
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
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>
#include <Math/RootFinder.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;


ClassImp(MSTMixFitter_ns::Utils)

namespace MSTMixFitter_ns {

RInside Utils::fR;

Utils::Utils(){

}

Utils::~Utils(){

}

int Utils::LoadRLibraries(std::vector<std::string> libraryNames)
{
	//Check libraries
	if(libraryNames.empty()){
		WARN_LOG("Empty library names, nothing to be done!");
		return -1;
	}

	//Load libraries
	for(size_t i=0;i<libraryNames.size();i++){
		std::string RCmd= Form("library('%s');",libraryNames[i].c_str());
		try{
			fR.parseEval(RCmd);	
		}
		catch(...){
			ERROR_LOG("Failed to load library "<<libraryNames[i]<<"!");
			return -1;
		}		
	}//end loop libraries

	return 0;

}//close LoadRLibraries()

int Utils::ClearRData()
{
	//## Clear R environment
	DEBUG_LOG("Clearing R environment...");
	try{
		fR.parseEvalQ("rm(list = ls(all = TRUE));");
	}
	catch(...){
		ERROR_LOG("Failed to clear R data!");
		return -1;
	}

	return 0;

}//close ClearRData()



TMatrixD* Utils::ConvertRVectToROOTMatrix(std::string RVect)
{
	//Check table name
	if(RVect==""){
		ERROR_LOG("Empty R vector name!");
		return nullptr;
	}

	//Store R table to NumericVector and then fill TMatrixD (double copy...not efficient!!!)
	TMatrixD* dataMatrix_ROOT= 0;
	try{
		Rcpp::NumericVector dataVect= fR.parseEval(RVect.c_str());
		long int N= fR.parseEval(Form("length(%s)",RVect.c_str()));
		//dataMatrix_ROOT= new TMatrixD(N,1);
		dataMatrix_ROOT= new TMatrixD(1,N);

		for(long int i=0;i<N;i++){
			//(*dataMatrix_ROOT)(i,0)= dataVect(i);
			(*dataMatrix_ROOT)(0,i)= dataVect(i);
		}
	}//close try block
	catch(...){
		ERROR_LOG("Failed to retrieve data table and relative size with imputed values in R!");
		return nullptr;
	}

	return dataMatrix_ROOT;

}//close ConvertRVectToROOTMatrix()


TMatrixD* Utils::ConvertRTableToROOTMatrix(std::string RTable)
{
	//Check table name
	if(RTable==""){
		ERROR_LOG("Empty R table name!");
		return nullptr;
	}

	//Store R table to Numeric matrix and then fill TMatrixD (double copy...not efficient!!!)
	TMatrixD* dataMatrix_ROOT= 0;
	try{
		Rcpp::NumericMatrix dataMatrix= fR.parseEval(RTable.c_str());
		long int N= fR.parseEval(Form("nrow(%s)",RTable.c_str()));
		long int NDim= fR.parseEval(Form("ncol(%s)",RTable.c_str()));
		
		dataMatrix_ROOT= new TMatrixD(N,NDim);
		for(long int i=0;i<N;i++){
			for(long int j=0;j<NDim;j++){
				(*dataMatrix_ROOT)(i,j)= dataMatrix(i,j);
			}
		}
	}//close try block
	catch(...){
		ERROR_LOG("Failed to retrieve data table and relative size with imputed values in R!");
		return nullptr;
	}

	return dataMatrix_ROOT;

}//close ConvertRTableToROOTMatrix()

int Utils::ImportMatrixInR(TMatrixD* dataMatrix,std::string dataname)
{
	//## Comvert matrix to R
	Rcpp::NumericMatrix* matrix_r= ConvertROOTMatrixToRMatrix(dataMatrix);
	if(!matrix_r){
		ERROR_LOG("Failed to convert ROOT matrix to R!");
		return -1;
	}

	//## Import in R prompt
	try{
		fR[dataname.c_str()]= *matrix_r;
	}
	catch(...){
		ERROR_LOG("Failed to import RNumeric matrix in R prompt!");
		return -1;
	}

	//## Delete R matrix
	if(matrix_r){
		delete matrix_r;
		matrix_r= 0;
	}

	return 0;

}//close ImportMatrixInR()


int Utils::ImportMatrixInR(std::vector<TMatrixD>& dataMatrix,std::string dataname)
{
	//## Comvert matrix to R
	Rcpp::NumericMatrix* matrix_r= ConvertROOTMatrixToRMatrix(dataMatrix);
	if(!matrix_r){
		ERROR_LOG("Failed to convert ROOT matrix to R!");
		return -1;
	}

	//## Import in R prompt
	try{
		fR[dataname.c_str()]= *matrix_r;
	}
	catch(...){
		ERROR_LOG("Failed to import RNumeric matrix in R prompt!");
		return -1;
	}

	//## Delete R matrix
	if(matrix_r){
		delete matrix_r;
		matrix_r= 0;
	}

	return 0;

}//close ImportMatrixInR()

Rcpp::NumericMatrix* Utils::ConvertROOTMatrixToRMatrix(TMatrixD* dataMatrix)
{
	//## Check data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return nullptr;
	}

	//## Create NumericMatrix
	long int nDim= dataMatrix->GetNcols();
	long int N= dataMatrix->GetNrows();
	Rcpp::NumericMatrix* matrix_r= new Rcpp::NumericMatrix(N,nDim);
	for(long int i=0;i<N;i++){
		for(long int j=0;j<nDim;j++){
			(*matrix_r)(i,j)= (*dataMatrix)(i,j);
		}//end loop dim
	}//end loop events

	return matrix_r;

}//close ConvertROOTMatrixToRMatrix()


Rcpp::NumericMatrix* Utils::ConvertROOTMatrixToRMatrix(std::vector<TMatrixD>& dataMatrix)
{
	//## Check data
	if(dataMatrix.empty()){
		ERROR_LOG("Empty data matrix list given!");
		return nullptr;
	}

	//## Create NumericMatrix
	long int nDim= dataMatrix[0].GetNcols();
	long int N= static_cast<long int>(dataMatrix.size());
	Rcpp::NumericMatrix* matrix_r= new Rcpp::NumericMatrix(N,nDim);
	for(long int i=0;i<N;i++){
		for(long int j=0;j<nDim;j++){
			(*matrix_r)(i,j)= dataMatrix[i](0,j);
		}//end loop dim
	}//end loop events

	return matrix_r;

}//close ConvertROOTMatrixToRMatrix()


TTree* Utils::MatrixToTTree(TMatrixD* dataMatrix,std::string treeName,std::string treeTitle)
{
	//## Check data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return nullptr;
	}
	int nDim= dataMatrix->GetNcols();
	long int N= dataMatrix->GetNrows();

	double x[nDim];
  TTree* tree = new TTree(treeName.c_str(),treeTitle.c_str()); 
  tree->Branch("x",x,Form("x[%d]/D",nDim));

  for (int i = 0; i<N;++i) { 
  	for (int j = 0; j < nDim; ++j) {
    	x[j] = (*dataMatrix)(i,j);
    }
  	tree->Fill();
  }

	return tree;

}//close MatrixToTTree()


int Utils::SetRandomMissingData(TMatrixD* dataMatrix,double missingDataFraction)
{
	//## Check data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return -1;
	}

	//Clone input data
	//TMatrixD* dataMatrixWithMiss= (TMatrixD*)dataMatrix->Clone();

	//## Set missing data at random
	int nDim= dataMatrix->GetNcols();
	long int N= dataMatrix->GetNrows();
	long int NElements= N*nDim;
	double NMissingElements= std::floor(missingDataFraction*NElements);
	
	long int missingCounter= 0;
	std::vector< std::vector<long int> > takenIndex;
	std::vector<long int> alltakenFlag;

	for(long int i=0;i<N;i++){
		takenIndex.push_back( std::vector<long int>() );
		alltakenFlag.push_back(0.);
		for(long int j=0;j<nDim;j++){
			takenIndex[i].push_back(0);
		}//end loop dim
	}//end loop events
	
	while(missingCounter<NMissingElements){

		if(missingCounter%100==0) INFO_LOG("--> "<<missingCounter<<"/"<<NMissingElements<<" events generated ...");

		long int eventId= static_cast<long int>(gRandom->Uniform(0,N));
		long int variableId= static_cast<long int>(gRandom->Uniform(0,nDim));

		//## Check if this combination has already be taken
		if(takenIndex[eventId][variableId]==0){
			takenIndex[eventId][variableId]= 1;//set as taken

			//## Now check if all patterns have been chosen
			bool isAllTaken= true;
			for(long int j=0;j<nDim;j++){
				if(j!=variableId && takenIndex[eventId][j]== 0) {
					isAllTaken= false;
					break;
				}
			}//end loop dim

			if(!isAllTaken){
				(*dataMatrix)(eventId,variableId)= TMath::SignalingNaN();
				missingCounter++;
			}
			else alltakenFlag[eventId]= true; 

		}//close if
	}//end loop while
	
	return 0;

}//close SetRandomMissingData()


TMatrixD* Utils::MakeRandomMissingData(TMatrixD* dataMatrix,double missingDataFraction)
{
	//## Check data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return nullptr;
	}

	//Clone input data
	TMatrixD* dataMatrixWithMiss= (TMatrixD*)dataMatrix->Clone();

	//Set missing data
	if(SetRandomMissingData(dataMatrixWithMiss,missingDataFraction)<0){
		ERROR_LOG("Failed to set missing data in input data matrix!");
		if(dataMatrixWithMiss){
			delete dataMatrixWithMiss;
			dataMatrixWithMiss= 0;
		}
		return nullptr;
	}
	
	return dataMatrixWithMiss;

}//close MakeRandomMissingData()


int Utils::DumpMatrixToAsciiFile(TMatrixD* dataMatrix,std::string filename)
{
	//Check filename
	if(filename==""){
		WARN_LOG("Empty output file specified!");
		return -1;
	}

	//Create output file
	FILE* fout= fopen(filename.c_str(),"w");

	//Write to file
	long int N= dataMatrix->GetNrows();
	long int NDim= dataMatrix->GetNcols();
	
	for(long int i=0;i<N;i++){
		for(long int j=0;j<NDim;j++){
			fprintf(fout,"%f  ",(*dataMatrix)(i,j));
		}//end loop dim
		fprintf(fout,"\n");
	}//end loop events

	//Close file
	fclose(fout);

	return 0;

}//close DumpMatrixToAsciiFile()




}//close namespace
