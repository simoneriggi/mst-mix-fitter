/**
* @file MSTMixtureFitter.cc
* @class MSTMixtureFitter
* @brief MSTMixtureFitter
*
* Fit a mixture of multivariate skew-t
* @author S. Riggi
* @date 18/09/2012
*/

#include <MSTMixtureFitter.h>
#include <MathUtils.h>
#include <ConfigParser.h>
#include <Logger.h>
#include <Utils.h>
#include <DataReader.h>
#include <KMeansClustering.h>

#include <TMinuit.h>
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
#include <TColor.h>
#include <TMatrixDEigen.h>

//#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(MSTMixFitter_ns::MSTMixtureFitter)

namespace MSTMixFitter_ns {

//Static members
int MSTMixtureFitter::fNDim;
int MSTMixtureFitter::fNComponents;
std::vector<TMatrixD> MSTMixtureFitter::fData;

MSTMixtureFitter::MSTMixtureFitter()
{
	fRTableName= "dataMatrix";

}//close constructor

MSTMixtureFitter::~MSTMixtureFitter()
{
	//Clear data
	ClearData();

}//close destructor

void MSTMixtureFitter::ClearData()
{
	//Clear data
	//...

}//close ClearData()

int MSTMixtureFitter::Init()
{
	INFO_LOG("Init data ...");

	//## Clear R environment
	INFO_LOG("Clearing R environment...");
	Utils::ClearRData();

	//### Initialize all necessary libraries
	INFO_LOG("Loading needed R libraries...");
	//std::vector<std::string> RLibraries{"mvtnorm","tmvtnorm","fMultivar","Matrix","moments"};
	std::vector<std::string> RLibraries{"Matrix"};
	if(Utils::LoadRLibraries(RLibraries)<0){
		ERROR_LOG("Failed to load one or more of these R libraries {Matrix}, check if they are installed!");
		return -1;
	}

	//## Initialize ROOT random generator
	delete gRandom;
	gRandom= new TRandom3(0);

	//## Get steering config from parser
	//- Main options
	fInputFileName= ConfigParser::fInputFileName;	
	fOutputFileName= ConfigParser::fOutputFileName;	
	fDataReadDelimiter= ConfigParser::fDataReadDelimiter;
	fNComponents= ConfigParser::fNComponents;
	fNDim= ConfigParser::fNDim;//this is overridden after reading the data

	//- Main EM options
	fNIterations= ConfigParser::fNIterations;
	fRunMinuitFitterAtConvergence= ConfigParser::fRunMinuitFitterAtConvergence;
	fFitter= ConfigParser::fFitterChoice;

	//- EM start values
	fParInitMethod= ConfigParser::fParInitMethod;
	fP_start= ConfigParser::fFractionStartPars;
	fNu_start= ConfigParser::fNuStartPars;
	fMu_start= ConfigParser::fMeanStartPars;
	fDelta_start= ConfigParser::fDeltaStartPars;
	fSigma_start= ConfigParser::fCovarianceStartPars;

	//- EM constraint options
	fFixMeanPars= ConfigParser::fFixMeanPars;
	fFixCovariancePars= ConfigParser::fFixCovariancePars;
	fFixDeltaPars= ConfigParser::fFixDeltaPars;
	fFixNuPars= ConfigParser::fFixNuPars;
	fFixFractionPars= ConfigParser::fFixFractionPars;
	fUseRandomRegenerationAfterStuck= ConfigParser::fUseRandomRegenerationAfterStuck;
	fUseConstraints= ConfigParser::fUseConstraints;
	fConstraintAlphaScale= ConfigParser::fConstraintAlphaScale;
	fConstraintAlphaTolerance= ConfigParser::fConstraintAlphaTolerance;
	
	fUseMeanBoundConstraint= ConfigParser::fUseMeanBoundConstraint;
	fMu_min= ConfigParser::fMeanMinBound;
	fMu_max= ConfigParser::fMeanMaxBound;
	
	fUseCovarianceBoundConstraint= ConfigParser::fUseCovarianceBoundConstraint;
	fSigma_min= ConfigParser::fCovarianceMinBound;
	fSigma_max= ConfigParser::fCovarianceMaxBound;

	fUseDeltaBoundConstraint= ConfigParser::fUseDeltaBoundConstraint;
	fDelta_min= ConfigParser::fDeltaMinBound;
	fDelta_max= ConfigParser::fDeltaMaxBound;

	fUseNuBoundConstraint= ConfigParser::fUseNuBoundConstraint;
	fNu_min= ConfigParser::fNuMinBound;
	fNu_max= ConfigParser::fNuMaxBound;

	/*
	fConstraintDataFileName= ConfigParser::fConstraintDataFileName;
	fUseMeanConstraint= ConfigParser::fUseMeanConstraint;
	fUseMeanDiffConstraint= ConfigParser::fUseMeanDiffConstraint;
	fUseCovarianceConstraint= ConfigParser::fUseCovarianceConstraint;
	fUseCovarianceBoundConstraint= ConfigParser::fUseCovarianceBoundConstraint;
	fUseCovarianceEigenConstraint= ConfigParser::fUseCovarianceEigenConstraint;
	fUseCovarianceEigenBoundConstraint= ConfigParser::fUseCovarianceEigenBoundConstraint;
	fUseLocationConstraint= ConfigParser::fUseLocationConstraint;
	fUseLocationDiffConstraint= ConfigParser::fUseLocationDiffConstraint;
	fUseScaleMatrixConstraint= ConfigParser::fUseScaleMatrixConstraint;
	fUseScaleMatrixBoundConstraint= ConfigParser::fUseScaleMatrixBoundConstraint;
	fUseScaleMatrixEigenConstraint= ConfigParser::fUseScaleMatrixEigenConstraint;
	fUseScaleMatrixEigenBoundConstraint= ConfigParser::fUseScaleMatrixEigenBoundConstraint;
	fUseDeltaBoundConstraint= ConfigParser::fUseDeltaBoundConstraint;
	fUseNuBoundConstraint= ConfigParser::fUseNuBoundConstraint;

	fSaveFullInfo= ConfigParser::fSaveFullInfo;

	fFixMeanDiffPar= false;

	fUseAnnealing= ConfigParser::fUseAnnealing;
	fAnnealingParStart= ConfigParser::fAnnealingParStart;	
	fAnnealingParStep= ConfigParser::fAnnealingParStep;

	fIsInteractiveRun= ConfigParser::fIsInteractiveRun;
	
	
	fRandomizeStartPar= ConfigParser::fUseRandomStart;
	fRandomizeStartMeanPar= ConfigParser::fUseRandomMeanStart;
	fRandomizeStartMeanDiffPar= ConfigParser::fUseRandomMeanDiffStart;
	fRandomizeStartCovariancePar= ConfigParser::fUseRandomCovarianceStart;
	fRandomizeStartCovarianceEigenPar= ConfigParser::fUseRandomCovarianceEigenStart;
	fRandomizeStartDeltaPar= ConfigParser::fUseRandomDeltaStart;
	fRandomizeStartNuPar= ConfigParser::fUseRandomNuStart;
	fRandomizeStartFractionPar= ConfigParser::fUseRandomFractionStart;
	fUseStoppingCriteria= ConfigParser::fUseStoppingCriteria;
	fEpsilon= ConfigParser::fEpsilon;

	fKsi_startFile= ConfigParser::fMeanStartPar;
	fSigma_startFile= ConfigParser::fCovarianceStartPar;
	fDelta_startFile= ConfigParser::fDeltaStartPar;
	fNu_startFile= ConfigParser::fNuStartPar;
	fP_startFile= ConfigParser::fFractionStartPar;
	fV_startFile= ConfigParser::fCovarianceStartPar;
	fSigmaEigen_startFile= ConfigParser::fCovarianceStartPar;


	fKsi_true= ConfigParser::fMeanTruePar;
	fSigma_true= ConfigParser::fCovarianceTruePar;
	fDelta_true= ConfigParser::fDeltaTruePar;
	fNu_true= ConfigParser::fNuTruePar;
	fP_true= ConfigParser::fFractionTruePar;

	fUseMeanOffset= ConfigParser::fUseMeanOffset;
	fKsiOffset= ConfigParser::fMeanOffsetPar;
	fUseDeltaOffset= ConfigParser::fUseDeltaOffset;
	fDeltaOffset= ConfigParser::fDeltaOffsetPar;

	fUseSigmaOffset= ConfigParser::fUseSigmaOffset;
	fSigmaOffset= ConfigParser::fSigmaOffsetPar;
		
	fClusterAGroups= ConfigParser::fClusterAGroups;
	fNClassificationGroups= ConfigParser::fNClassificationGroups;
	fClusterLnAMin= ConfigParser::fClusterLnAMin;
	fClusterLnAMax= ConfigParser::fClusterLnAMax;


	fMeanParTolerance_min= ConfigParser::fMeanParTolerance_min;	
	fMeanParTolerance_max= ConfigParser::fMeanParTolerance_max;
	fMeanDiffParTolerance_min= ConfigParser::fMeanDiffParTolerance_min;	
	fMeanDiffParTolerance_max= ConfigParser::fMeanDiffParTolerance_max;
	fSigmaParTolerance_min= ConfigParser::fSigmaParTolerance_min;
	fSigmaParTolerance_max= ConfigParser::fSigmaParTolerance_max;
	fDeltaParTolerance_min= ConfigParser::fDeltaParTolerance_min;
	fDeltaParTolerance_max= ConfigParser::fDeltaParTolerance_max;

	fMeanConstraintSign= ConfigParser::fMeanConstraintSign;	
	fSigmaConstraintSign= ConfigParser::fSigmaConstraintSign;	
	
	fDrawMinRange= ConfigParser::fDrawMinRange;
	fDrawMaxRange= ConfigParser::fDrawMaxRange;
	fDrawNBins= ConfigParser::fDrawNBins;
	

	//## Compute Ndof
	fNdof= 0;
	int nModelPars= 0;	
	if(!fFixMeanPar) nModelPars+= fNDim;
	if(!fFixCovariancePar) nModelPars+= fNDim*(fNDim+1)/2.;
	if(!fFixDeltaPar) nModelPars+= fNDim;
	if(!fFixNuPar) nModelPars+= fNDim;
	
	int nMixturePars= fNComponents*nModelPars;
	if(!fFixFractionPar) nMixturePars+= fNComponents-1;
	fNdof= nMixturePars;
	cout<<"MSTMixtureFitter::Init(): Fit Ndof="<<fNdof<<endl;

	//## Define data structures
	fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
	fOutputFile->cd();

	fFitInfo= new TTree("FitInfo","FitInfo");
	fFitInfo->Branch("NDim",&fNDim,"NDim/I");
	fFitInfo->Branch("NComponents",&fNComponents,"NComponents/I");
	fFitInfo->Branch("NIterations",&fNIterations,"NIterations/I");
	fFitInfo->Branch("Fraction",fFraction,"Fraction[NComponents]/D");
	fFitInfo->Branch("FractionErr",fFractionErr,"FractionErr[NComponents]/D");
	fFitInfo->Branch("XmaxMean",fXmaxMean,"XmaxMean[NComponents]/D");
	fFitInfo->Branch("XmaxMeanOffset",fXmaxMeanOffset,"XmaxMeanOffset[NComponents]/D");
	fFitInfo->Branch("NmuMean",fNmuMean,"NmuMean[NComponents]/D");
	fFitInfo->Branch("NmuMeanOffset",fNmuMeanOffset,"NmuMeanOffset[NComponents]/D");
	fFitInfo->Branch("XmaxVar",fXmaxVar,"XmaxVar[NComponents]/D");
	fFitInfo->Branch("XmaxVarOffset",fXmaxVarOffset,"XmaxVarOffset[NComponents]/D");	
	fFitInfo->Branch("NmuVar",fNmuVar,"NmuVar[NComponents]/D");
	fFitInfo->Branch("NmuVarOffset",fNmuVarOffset,"NmuVarOffset[NComponents]/D");	
	fFitInfo->Branch("XmaxNmuCov",fXmaxNmuCov,"XmaxNmuCov[NComponents]/D");
	fFitInfo->Branch("XmaxDelta",fXmaxDelta,"XmaxDelta[NComponents]/D");
	fFitInfo->Branch("XmaxDeltaOffset",fXmaxDeltaOffset,"XmaxDeltaOffset[NComponents]/D");	
	fFitInfo->Branch("NmuDelta",fNmuDelta,"NmuDelta[NComponents]/D");
	fFitInfo->Branch("NmuDeltaOffset",fNmuDeltaOffset,"NmuDeltaOffset[NComponents]/D");
	fFitInfo->Branch("Ndf",fNdf,"Ndf[NComponents]/D");
	fFitInfo->Branch("StartXmaxMean",fXmaxMean_start,"StartXmaxMean[NComponents]/D");
	fFitInfo->Branch("StartNmuMean",fNmuMean_start,"StartNmuMean[NComponents]/D");
	fFitInfo->Branch("StartXmaxVar",fXmaxVar_start,"StartXmaxVar[NComponents]/D");
	fFitInfo->Branch("StartNmuVar",fNmuVar_start,"StartNmuVar[NComponents]/D");
	fFitInfo->Branch("StartXmaxNmuCov",fXmaxNmuCov_start,"StartXmaxNmuCov[NComponents]/D");
	fFitInfo->Branch("StartXmaxDelta",fXmaxDelta_start,"StartXmaxDelta[NComponents]/D");
	fFitInfo->Branch("StartNmuDelta",fNmuDelta_start,"StartNmuDelta[NComponents]/D");
	fFitInfo->Branch("StartNdf",fNdf_start,"StartNdf[NComponents]/D");

	fFitInfo->Branch("TrueFraction",fFraction_true,"TrueFraction[NComponents]/D");
	fFitInfo->Branch("TrueXmaxMean",fXmaxMean_true,"TrueXmaxMean[NComponents]/D");
	fFitInfo->Branch("TrueNmuMean",fNmuMean_true,"TrueNmuMean[NComponents]/D");
	fFitInfo->Branch("TrueXmaxVar",fXmaxVar_true,"TrueXmaxVar[NComponents]/D");
	fFitInfo->Branch("TrueNmuVar",fNmuVar_true,"TrueNmuVar[NComponents]/D");
	fFitInfo->Branch("TrueXmaxNmuCov",fXmaxNmuCov_true,"TrueXmaxNmuCov[NComponents]/D");
	fFitInfo->Branch("TrueXmaxDelta",fXmaxDelta_true,"TrueXmaxDelta[NComponents]/D");
	fFitInfo->Branch("TrueNmuDelta",fNmuDelta_true,"TrueNmuDelta[NComponents]/D");
	fFitInfo->Branch("TrueNdf",fNdf_true,"TrueNdf[NComponents]/D");
	fFitInfo->Branch("TotClassificationEfficiency",&fTotClassificationEfficiency,"TotClassificationEfficiency/D");
	fFitInfo->Branch("ClassificationEfficiency",fClassificationEfficiency,"ClassificationEfficiency[NComponents]/D");
	fFitInfo->Branch("NEventPerGroup",fNEventPerGroup,"NEventPerGroup[NComponents]/D");
	fFitInfo->Branch("Likelihood",&fLogLikelihood,"Likelihood/D");
	fFitInfo->Branch("Ndof",&fNdof,"Ndof/D");
	fFitInfo->Branch("FitStatus",&fFitStatus,"FitStatus/I");
	
	fFitIterInfo= new TTree("FitIterInfo","FitIterInfo");	
	fFitIterInfo->Branch("IterNo",&fIterNo,"IterNo/I");
	fFitIterInfo->Branch("Likelihood",&fLogLikelihood,"Likelihood/D");
	fFitIterInfo->Branch("NComponents",&fNComponents,"NComponents/I");
	fFitIterInfo->Branch("Fraction",fFraction,"Fraction[NComponents]/D");	
	fFitIterInfo->Branch("XmaxMean",fXmaxMean,"XmaxMean[NComponents]/D");
	fFitIterInfo->Branch("NmuMean",fNmuMean,"NmuMean[NComponents]/D");
	fFitIterInfo->Branch("XmaxVar",fXmaxVar,"XmaxVar[NComponents]/D");
	fFitIterInfo->Branch("NmuVar",fNmuVar,"NmuVar[NComponents]/D");
	fFitIterInfo->Branch("XmaxNmuCov",fXmaxNmuCov,"XmaxNmuCov[NComponents]/D");
	fFitIterInfo->Branch("XmaxDelta",fXmaxDelta,"XmaxDelta[NComponents]/D");
	fFitIterInfo->Branch("NmuDelta",fNmuDelta,"NmuDelta[NComponents]/D");
	fFitIterInfo->Branch("Ndf",fNdf,"Ndf[NComponents]/D");

	fClassificationInfo= new TTree("ClassificationInfo","ClassificationInfo");	
	fClassificationInfo->Branch("TrueType",&fTrueType,"TrueType/I");
	fClassificationInfo->Branch("RecType",&fRecType,"RecType/I");


	//## Get constrain graph from file
	TFile* fConstraintFile= new TFile(fConstraintDataFileName.c_str(),"READ");
	if( !fConstraintFile || fConstraintFile->IsZombie() ){
		cerr<<"MSTMixtureFitter::Init(): Cannot open file with constrain data ...exit!"<<endl;
		exit(1);
	}

	TGraph* KsiConstraintGraph_epos[fNDim];
	TGraph* KsiConstraintGraph_sibyll[fNDim];
	TGraph* KsiConstraintGraph_qgsjetII[fNDim];

	TGraph* KsiConstraintGraph_min[fNDim];
	TGraph* KsiConstraintGraph_max[fNDim];

	TGraph* SigmaVarConstraintGraph_min[fNDim];
	TGraph* SigmaVarConstraintGraph_max[fNDim];
	TGraph* SigmaCovConstraintGraph_min[fNDim];
	TGraph* SigmaCovConstraintGraph_max[fNDim];
	TGraph* DeltaConstraintGraph_min[fNDim];
	TGraph* DeltaConstraintGraph_max[fNDim];
	TGraph* NuConstraintGraph_min= (TGraph*)fConstraintFile->Get("NuConstraint_min");
	TGraph* NuConstraintGraph_max= (TGraph*)fConstraintFile->Get("NuConstraint_max");;
	TGraph* SigmaEigenConstraintGraph_min[fNDim];
	TGraph* SigmaEigenConstraintGraph_max[fNDim];

	for(int j=0;j<fNDim;j++){
		KsiConstraintGraph_min[j]= (TGraph*)fConstraintFile->Get( Form("KsiConstraint_min_%d",j+1) );
		KsiConstraintGraph_max[j]= (TGraph*)fConstraintFile->Get( Form("KsiConstraint_max_%d",j+1) );
		
		SigmaVarConstraintGraph_min[j]= (TGraph*)fConstraintFile->Get(Form("SigmaVarConstraint_min_%d",j+1) );
		SigmaVarConstraintGraph_max[j]= (TGraph*)fConstraintFile->Get(Form("SigmaVarConstraint_max_%d",j+1) );

		SigmaCovConstraintGraph_min[j]= (TGraph*)fConstraintFile->Get( Form("SigmaCovConstraint_min_%d",j+1) );
		SigmaCovConstraintGraph_max[j]= (TGraph*)fConstraintFile->Get( Form("SigmaCovConstraint_max_%d",j+1) );		

		DeltaConstraintGraph_min[j]= (TGraph*)fConstraintFile->Get(Form("DeltaConstraint_min_%d",j+1) );
		DeltaConstraintGraph_max[j]= (TGraph*)fConstraintFile->Get(Form("DeltaConstraint_max_%d",j+1) );
	
		KsiConstraintGraph_epos[j]= (TGraph*)fConstraintFile->Get( Form("Ksi_epos_%d",j+1) );
		KsiConstraintGraph_sibyll[j]= (TGraph*)fConstraintFile->Get( Form("Ksi_sibyll_%d",j+1) );
		KsiConstraintGraph_qgsjetII[j]= (TGraph*)fConstraintFile->Get( Form("Ksi_qgsjetII_%d",j+1) );
		
		SigmaEigenConstraintGraph_min[j]= (TGraph*)fConstraintFile->Get(Form("EigenConstraint_min_%d",j+1) );
		SigmaEigenConstraintGraph_max[j]= (TGraph*)fConstraintFile->Get(Form("EigenConstraint_max_%d",j+1) );

	}//end loop dim



	cout<<"*** CONSTRAIN DATA ***"<<endl;
	for(int k=0;k<fNComponents;k++){
		TMatrixD Sigma_min(fNDim,fNDim);
		TMatrixD Sigma_max(fNDim,fNDim);

		TMatrixD Ksi_min(fNDim,1);
		TMatrixD Ksi_max(fNDim,1);
		
		TMatrixD Delta_min(fNDim,1);
		TMatrixD Delta_max(fNDim,1);

		TMatrixD SigmaEigen_min(fNDim,1);
		TMatrixD SigmaEigen_max(fNDim,1);

		double Nu_min= NuConstraintGraph_min->Eval(log(fClusterAGroups[k]));
		double Nu_max= NuConstraintGraph_max->Eval(log(fClusterAGroups[k]));
		
		for(int j=0;j<fNDim;j++){
	
			for(int l=0;l<fNDim;l++){
				if(j==l){
					//Sigma_min(j,j)= SigmaVarConstraintGraph_min[j]->Eval( log(fClusterAGroups[k]) );
					//Sigma_max(j,j)= SigmaVarConstraintGraph_max[j]->Eval( log(fClusterAGroups[k]) );
					Sigma_min(j,j)= (1-(*fSigmaParTolerance_min)(j,j))*SigmaVarConstraintGraph_min[j]->Eval( log(fClusterAGroups[k]) );
					Sigma_max(j,j)= (1+(*fSigmaParTolerance_max)(j,j))*SigmaVarConstraintGraph_max[j]->Eval( log(fClusterAGroups[k]) );
				}
				else{
					//Sigma_min(j,l)= SigmaCovConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
					//Sigma_max(j,l)= SigmaCovConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));
					Sigma_min(j,l)= SigmaCovConstraintGraph_min[j]->Eval(log(fClusterAGroups[k])) - fabs(SigmaCovConstraintGraph_min[j]->Eval(log(fClusterAGroups[k])))*(*fSigmaParTolerance_min)(j,l);
					Sigma_max(j,l)= SigmaCovConstraintGraph_max[j]->Eval(log(fClusterAGroups[k])) + fabs(SigmaCovConstraintGraph_max[j]->Eval(log(fClusterAGroups[k])))*(*fSigmaParTolerance_max)(j,l);
				}
			}

			//Ksi_min(j,0)= KsiConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
			//Ksi_max(j,0)= KsiConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));
			Ksi_min(j,0)= (1-(*fMeanParTolerance_min)(j,0))*KsiConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
			Ksi_max(j,0)= (1+(*fMeanParTolerance_max)(j,0))*KsiConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));

			//Delta_min(j,0)= DeltaConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
			//Delta_max(j,0)= DeltaConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));
			Delta_min(j,0)= (1-(*fDeltaParTolerance_min)(j,0))*DeltaConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
			Delta_max(j,0)= (1+(*fDeltaParTolerance_max)(j,0))*DeltaConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));
			
			SigmaEigen_min(j,0)= SigmaEigenConstraintGraph_min[j]->Eval(log(fClusterAGroups[k]));
			SigmaEigen_max(j,0)= SigmaEigenConstraintGraph_max[j]->Eval(log(fClusterAGroups[k]));

		}//end loop dim	
		

		fSigma_min.push_back(Sigma_min);
		fSigma_max.push_back(Sigma_max);

		fSigmaEigen_min.push_back(SigmaEigen_min);
		fSigmaEigen_max.push_back(SigmaEigen_max);

		fKsi_min.push_back(Ksi_min);
		fKsi_max.push_back(Ksi_max);

		fDelta_min.push_back(Delta_min);
		fDelta_max.push_back(Delta_max);

		fNu_min.push_back(Nu_min);
		fNu_max.push_back(Nu_max);


		cout<<"== Component "<<k+1<<" =="<<endl;	
		cout<<"Nu min/max= "<<Nu_min<<"/"<<Nu_max<<endl;

		cout<<"Ksi min= (";
		for(int j=0;j<fNDim-1;j++) cout<<Ksi_min(j,0)<<",";
		cout<<Ksi_min(fNDim-1,0)<<")"<<endl; 

		cout<<"Ksi max= (";
		for(int j=0;j<fNDim-1;j++) cout<<Ksi_max(j,0)<<",";
		cout<<Ksi_max(fNDim-1,0)<<")"<<endl; 
			
		cout<<"delta min = (";
		for(int j=0;j<fNDim-1;j++) cout<<Delta_min(j,0)<<",";
		cout<<Delta_min(fNDim-1,0)<<")"<<endl; 
		cout<<"delta max = (";
		for(int j=0;j<fNDim-1;j++) cout<<Delta_max(j,0)<<",";
		cout<<Delta_max(fNDim-1,0)<<")"<<endl; 
			
		cout<<"Sigma min = (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<Sigma_min(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<"Sigma max = (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<Sigma_max(j,l)<<",";
		}	
		cout<<")"<<endl;

		
	}//end loop components


	double meanDiffTolerance= 0;
	//double meanDiffTolerance= 1;
	//double meanDiffTolerance= 0.;

	for(int k=0;k<fNComponents;k++){
		TMatrixD KsiDiff_min(fNDim,1);
		TMatrixD KsiDiff_max(fNDim,1);

		for(int j=0;j<fNDim;j++){
			double ksiDiff_epos= KsiConstraintGraph_epos[j]->Eval(log(fClusterAGroups[k]))-KsiConstraintGraph_epos[j]->Eval(log(fClusterAGroups[0]));
			double ksiDiff_sibyll= KsiConstraintGraph_sibyll[j]->Eval(log(fClusterAGroups[k]))-KsiConstraintGraph_sibyll[j]->Eval(log(fClusterAGroups[0]));
			double ksiDiff_qgsjetII= KsiConstraintGraph_qgsjetII[j]->Eval(log(fClusterAGroups[k]))-KsiConstraintGraph_qgsjetII[j]->Eval(log(fClusterAGroups[0]));
			
			double minDiff= min(min(ksiDiff_epos,ksiDiff_sibyll),ksiDiff_qgsjetII);
			double maxDiff= max(max(ksiDiff_epos,ksiDiff_sibyll),ksiDiff_qgsjetII);
			double toleranceRegion= fabs(maxDiff-minDiff)*meanDiffTolerance;
			//KsiDiff_min(j,0)= minDiff - toleranceRegion;
			//KsiDiff_max(j,0)= maxDiff + toleranceRegion;
			KsiDiff_min(j,0)= minDiff*(1-(*fMeanDiffParTolerance_min)(j,0));
			KsiDiff_max(j,0)= maxDiff*(1+(*fMeanDiffParTolerance_max)(j,0));

			cout<<"COMPONENT "<<k+1<<"  j="<<j<<"  ksiDiffModels("<<ksiDiff_epos<<","<<ksiDiff_qgsjetII<<","<<ksiDiff_sibyll<<")  min/max: "<<KsiDiff_min(j,0)<<"/"<<KsiDiff_max(j,0)<<"  toleranceRegion="<<toleranceRegion<<"  fabs(maxDiff-minDiff)="<<fabs(maxDiff-minDiff)<<endl;
		}//end loop dim

		fKsiDiff_min.push_back(KsiDiff_min);
		fKsiDiff_max.push_back(KsiDiff_max);
	
	}//end loop components
	*/

	return 0;

}//close MSTMixtureFitter::Init()


int MSTMixtureFitter::ReadData()
{
	//## Read ascii data to ROOT matrix
	INFO_LOG("Start reading data ...");
	//int status= DataReader::ReadAscii(fData,fInputFileName,fDataReadDelimiter);
	int status= DataReader::ReadAscii(fData,fInputFileName,"");
	if(status<0){
		ERROR_LOG("Failed to read data from file "<<fInputFileName<<"!");
		return -1;
	}
	if(fData.empty()){
		ERROR_LOG("Empty data read from file "<<fInputFileName<<"!");
		return -1;
	}
	fNDim= fData[0].GetNcols();
	fN= static_cast<long int>(fData.size());
	INFO_LOG("Read "<<fN<<" x "<<fNDim<<" data table...");

	//## Import matrix also in R for later use
	INFO_LOG("Reading data table in R...");
	if(DataReader::ReadAsciiInR(fInputFileName,fDataReadDelimiter,fRTableName)<0){
		ERROR_LOG("Failed to read ascii file "<<fInputFileName<<" and import it as an R table with name "<<fRTableName<<"!");
		return -1;
	}

	//## Find data ranges/mean/cov matrix
	fMinDataRange.clear();
	fMinDataRange.resize(0);
	fMaxDataRange.clear();
	fMaxDataRange.resize(0);

	for(int j=0;j<fNDim;j++){
		fMinDataRange.push_back(1.e+99);
		fMaxDataRange.push_back(-1.e+99);
		fMeanData.push_back(0.);
		fVarianceData.push_back(0.);
	}

	for(size_t i=0;i<fData.size();i++){
		for(int j=0;j<fNDim;j++){
			if( (fData[i])(0,j)<fMinDataRange[j] ) fMinDataRange[j]= (fData[i])(0,j);
			if( (fData[i])(0,j)>fMaxDataRange[j] ) fMaxDataRange[j]= (fData[i])(0,j);	
			fMeanData[j]+= (fData[i])(0,j);
		}//end loop dim	
	}//end loop data entries
	
	for(int j=0;j<fNDim;j++) fMeanData[j]/= (double)(fN);
		
	for(size_t i=0;i<fData.size();i++){
		for(int j=0;j<fNDim;j++){
			fVarianceData[j]+= pow( (fData[i])(0,j)-fMeanData[j],2);	
		}
	}
	
	for(int j=0;j<fNDim;j++) fVarianceData[j]/= (double)(fN-1);


	std::stringstream ss;
	ss<<"Data range: ";
	for(int j=0;j<fNDim;j++){
		ss<<"DataVar"<<j+1<<"{range=["<<fMinDataRange[j]<<","<<fMaxDataRange[j]<<"], mean="<<fMeanData[j]<<", variance="<<fVarianceData[j]<<"}, ";
	}
	INFO_LOG(ss.str());
	

	/*
	//## Init data histos
	double histoRangeX= fMaxDataRange[0]-fMinDataRange[0];
	double histoRangeTolX= 0.5;
	double histoMinX= fMinDataRange[0] - histoRangeTolX*fabs(histoRangeX);
	double histoMaxX= fMaxDataRange[0] + histoRangeTolX*fabs(histoRangeX);

	double histoRangeY= fMaxDataRange[0]-fMinDataRange[0];
	double histoRangeTolY= 0.5;
	double histoMinY= fMinDataRange[0] - histoRangeTolY*fabs(histoRangeY);
	double histoMaxY= fMaxDataRange[0] + histoRangeTolY*fabs(histoRangeY);
	if(fNDim>1){
		histoRangeY= fMaxDataRange[1]-fMinDataRange[1];
		histoMinY= fMinDataRange[1] - histoRangeTolY*fabs(histoRangeY);
		histoMaxY= fMaxDataRange[1] + histoRangeTolY*fabs(histoRangeY);
	}
	

	fDataGraph= new TGraph;

	//fDataHisto= new TH2D("DataHisto","DataHisto",30,600,1200,20,0,2);
	//fDataHisto= new TH2D("DataHisto","DataHisto",20,0,2,20,0,2);
	//fDataHisto= new TH2D("DataHisto","DataHisto",30,600,1200,40,0,200);
	//fDataHisto= new TH2D("DataHisto","DataHisto",25,histoMinX,histoMaxX,25,histoMinY,histoMaxY);
	fDataHisto= new TH2D("DataHisto","DataHisto",fDrawNBins[0],fDrawMinRange[0],fDrawMaxRange[0],fDrawNBins[1],fDrawMinRange[1],fDrawMaxRange[1]);
	fDataHisto->SetLineColor(kBlack);
	fDataHisto->SetMarkerStyle(8);
	fDataHisto->SetMarkerSize(1.1);
	fDataHisto->Sumw2();
	
	//fFitHisto= new TH2D("FitHisto","FitHisto",30,600,1200,20,0,2);
	//fFitHisto= new TH2D("FitHisto","FitHisto",20,0,2,20,0,2);
	//fFitHisto= new TH2D("FitHisto","FitHisto",30,600,1200,20,0,200);
	//fFitHisto= new TH2D("FitHisto","FitHisto",25,histoMinX,histoMaxX,25,histoMinY,histoMaxY);	
	fFitHisto= new TH2D("FitHisto","FitHisto",fDrawNBins[0],fDrawMinRange[0],fDrawMaxRange[0],fDrawNBins[1],fDrawMinRange[1],fDrawMaxRange[1]);	
	fFitHisto->SetLineColor(kRed);
	fFitHisto->SetMarkerStyle(8);
	fFitHisto->SetMarkerSize(1.1);
	fFitHisto->Sumw2();

	//fDataXHisto= new TH1D("DataXHisto","DataXHisto",30,600,1200);
	//fDataXHisto= new TH1D("DataXHisto","DataXHisto",25,histoMinX,histoMaxX);
	fDataXHisto= new TH1D("DataXHisto","DataXHisto",fDrawNBins[0],fDrawMinRange[0],fDrawMaxRange[0]);
	fDataXHisto->SetLineColor(kBlack);
	fDataXHisto->SetMarkerStyle(8);
	fDataXHisto->SetMarkerSize(1.1);
	fDataXHisto->Sumw2();
	
	//fDataYHisto= new TH1D("DataYHisto","DataYHisto",40,0,200);
	//fDataYHisto= new TH1D("DataYHisto","DataYHisto",25,histoMinY,histoMaxY);
	fDataYHisto= new TH1D("DataYHisto","DataYHisto",fDrawNBins[1],fDrawMinRange[1],fDrawMaxRange[1]);
	fDataYHisto->SetLineColor(kBlack);
	fDataYHisto->SetMarkerStyle(8);
	fDataYHisto->SetMarkerSize(1.1);
	fDataYHisto->Sumw2();

	//fFitXHisto= new TH1D("FitXHisto","FitXHisto",30,600,1200);
	//fFitXHisto= new TH1D("FitXHisto","FitXHisto",25,histoMinX,histoMaxX);
	fFitXHisto= new TH1D("FitXHisto","FitXHisto",fDrawNBins[0],fDrawMinRange[0],fDrawMaxRange[0]);
	fFitXHisto->SetLineColor(kBlack);
	fFitXHisto->SetMarkerStyle(8);
	fFitXHisto->SetMarkerSize(1.1);
	fFitXHisto->Sumw2();
	
	//fFitYHisto= new TH1D("FitYHisto","FitYHisto",40,0,200);
	//fFitYHisto= new TH1D("FitYHisto","FitYHisto",25,histoMinY,histoMaxY);
	fFitYHisto= new TH1D("FitYHisto","FitYHisto",fDrawNBins[1],fDrawMinRange[1],fDrawMaxRange[1]);
	fFitYHisto->SetLineColor(kBlack);
	fFitYHisto->SetMarkerStyle(8);
	fFitYHisto->SetMarkerSize(1.1);
	fFitYHisto->Sumw2();

	for(int k=0;k<fNComponents;k++){
		int color= kBlack;
		if(k==0) color=	kRed;
		else if(k==1) color= kGreen+1;
		else if(k==2) color= kBlue+1;

		TString histoName= Form("FitXHisto_Component%d",k+1);
		//fFitXComponentsHisto[k]= new TH1D(histoName,histoName,25,histoMinX,histoMaxX);
		fFitXComponentsHisto[k]= new TH1D(histoName,histoName,fDrawNBins[0],fDrawMinRange[0],fDrawMaxRange[0]);
		fFitXComponentsHisto[k]->SetLineColor(color);
		fFitXComponentsHisto[k]->SetMarkerStyle(8);
		fFitXComponentsHisto[k]->SetMarkerSize(1.1);
		fFitXComponentsHisto[k]->Sumw2();

		histoName= Form("FitYHisto_Component%d",k+1);
		//fFitYComponentsHisto[k]= new TH1D(histoName,histoName,25,histoMinY,histoMaxY);
		fFitYComponentsHisto[k]= new TH1D(histoName,histoName,fDrawNBins[1],fDrawMinRange[1],fDrawMaxRange[1]);		
		fFitYComponentsHisto[k]->SetLineColor(color);
		fFitYComponentsHisto[k]->SetMarkerStyle(8);
		fFitYComponentsHisto[k]->SetMarkerSize(1.1);
		fFitYComponentsHisto[k]->Sumw2();
	}//end loop components
	*/

	return 0;	

}//close ReadData()

int MSTMixtureFitter::Run()
{
	//## Initialize data
	if(Init()<0){
		ERROR_LOG("Failed to initialize fitter data!");
		return -1;
	}

	//## Reading data
	if(ReadData()<0){
		ERROR_LOG("Failed to read data!");
		return -1;
	}

	//## Running fitter
	if(fFitter==eEM){
		if(RunEMFitter()<0){
			ERROR_LOG("Failed to run EM fitter!");
			return -1;
		}
	}
	else if(fFitter==eMINUIT){
		if(RunFitter(true)<0){
			ERROR_LOG("Failed to run MINUIT fitter!");
			return -1;
		}
	}
	else{
		ERROR_LOG("Invalid fitter ("<<fFitter<<") specified!");
		return -1;
	}

	//Save data
	//...

	return 0;

}//close Run()


int MSTMixtureFitter::RunEMFitter()
{
	//Initialize EM data
	if(RunEM_Init()<0){
		ERROR_LOG("Failed to initialize EM algorithm data!");
		return -1;
	}

	//## Start EM iteration loop
	INFO_LOG("Starting EM iteration loop (#"<<fNIterations<<" max niters) ...");
	fLogLikelihood= 0.;
	fIterNo= 0;
	fFitStatus= eInit;
	int badIterCounter= 0;

	for(int iter=0;iter<fNIterations;iter++)
	{
		fIterNo= iter;
		
		//############################################################
		//## E step: Compute the expectation e0=tau, LL
		//############################################################			
		double LL= 0;
		if(RunEM_EStep(LL)<0){
			ERROR_LOG("EM EStep failed at iter no. "<<iter+1<<"!");
			return -1;
		}
		double DeltaLogL_iter= LL-fLogLikelihood;
		fLogLikelihood= LL;
		fIterLogLikelihood.push_back(fLogLikelihood);
		if(DeltaLogL_iter>0){
			badIterCounter= 0;
			/*
			for(int k=0;k<fNComponents;k++){
				fP_best[k]= fP[k]; 
				fKsi_best[k]= fKsi[k];
				fMu_best[k]= fMu[k];
				fSigma_best[k]= fSigma[k];
				fOmega_best[k]= fOmega[k];
				fDelta_best[k]= fDelta[k];
				fNu_best[k]= fNu[k];
				fLikelihood_best= LL;
			}//end loop components
			*/
		}//close if
		else{
			fFitStatus= eTroubles;
			badIterCounter++;
		}
	
		//Check stopping criteria
		if(iter>3){
			double ak= (fIterLogLikelihood[iter]-fIterLogLikelihood[iter-1])/(fIterLogLikelihood[iter-1]-fIterLogLikelihood[iter-2]);
			double LogL_inf= fIterLogLikelihood[iter-1] + (fIterLogLikelihood[iter]-fIterLogLikelihood[iter-1])/(1.-ak);
			double DeltaLogL= LogL_inf-fIterLogLikelihood[iter];
			INFO_LOG("Checking stopping conditions at iter no. "<<iter<<": ak="<<ak<<", LogLikelihood_inf="<<LogL_inf<<", DeltaLogL_iter="<<DeltaLogL_iter<<", DeltaLogL="<<DeltaLogL);

			if(fUseStoppingCriteria && DeltaLogL>=0 && fabs(DeltaLogL)<fEpsilon){
				INFO_LOG("Stopping criteria reached (LL="<<LL<<", DeltaLogL="<<DeltaLogL<<"< eps="<<fEpsilon<<"), stop iteration.");
				fFitStatus= eConverged;
				break;
			}
		}//close if
		
		//############################################################
		//## M step: Update the fit parameters
		//############################################################
		RunEM_MStep();
		
		//############################################################
		//## Constrain Step
		//############################################################
		if(fUseConstraints && RunEM_ConstrainStep()<0){
			ERROR_LOG("Failed to run EM constraint step!");
			return -1;
		}

		//############################################################
		//## DUMP INFO
		//############################################################			
		cout<<"*** ITER NO. "<<iter<<" ***"<<endl;		
		PrintPars();
		cout<<"==> LL= "<<LL<<"  DeltaL="<<DeltaLogL_iter<<endl;
		cout<<"***************************"<<endl;
		cout<<endl;

	}//end loop EM iterations


	//## Run numerical fitter at convergence?
	if(fFitStatus==eConverged && fRunMinuitFitterAtConvergence){
		INFO_LOG("Running MINUIT fitter after EM convergence...");
		RunFitter(false);
	}
	

	return 0;

}//close RunEMFitter()


int MSTMixtureFitter::InitializeFitPars()
{
	//Initialize start fit parameters
	for(int k=0;k<fNComponents;k++){
		TMatrixD tmpVect(1,fNDim);
		TMatrixD tmpMatrix(fNDim,fNDim);		
	
		fKsi_start.push_back(tmpVect);
		fOmega_start.push_back(tmpMatrix);
		
		fKsi.push_back(tmpVect);
		fMu.push_back(tmpVect);
		fDelta.push_back(tmpVect);	
		fSigma.push_back(tmpMatrix);
		fSigmaEigen.push_back(tmpVect);
		fSigmaInv.push_back(tmpMatrix);
		fSigmaDet.push_back(0.);
		fOmega.push_back(tmpMatrix);
		fOmegaEigen.push_back(tmpMatrix);
		fOmegaInv.push_back(tmpMatrix);
		fOmegaDet.push_back(0.);
		fP.push_back(1./fNComponents);
		fNu.push_back(0.);

		fKsi_safe.push_back(tmpVect);
		fMu_safe.push_back(tmpVect);
		fDelta_safe.push_back(tmpVect);	
		fOmega_safe.push_back(tmpMatrix);
		fOmegaEigen_safe.push_back(tmpMatrix);		
		fSigma_safe.push_back(tmpMatrix);
		fSigmaEigen_safe.push_back(tmpVect);		
		fP_safe.push_back(1./fNComponents);
		fNu_safe.push_back(0.);
	}

	//## Initialize mixture pars
	int status= 0;
	if(fParInitMethod==eKMEANS){
		status= InitParsToKMeans();
	}
	else if(fParInitMethod==eRANDOM){	
		ERROR_LOG("Randomized parameter initialization not yet implemented!");
		status= -1;		
	}
	else if(fParInitMethod==eUSER){
		status= InitParsToUser();
	}
	else{
		ERROR_LOG("Invalid par initialization method given ("<<fParInitMethod<<")!");
		return -1;
	}

	if(status<0){
		ERROR_LOG("EM parameter initialization failed!");
		return -1;
	}

	
	//## Compute Sigma & Omega matrix eigenvalues
	INFO_LOG("Computing Sigma & Omega matrix eigenvalues...");
	for(int k=0;k<fNComponents;k++){
		TMatrixD eigenVects(fNDim,fNDim);
		INFO_LOG("Printing Omega matrix for component no. "<<k+1<<"...");
		fOmega[k].Print();
		if(MathUtils::ComputeSymMatrixEigenvalues(fOmegaEigen[k],eigenVects,fOmega[k])<0){
			WARN_LOG("Failed to compute scale matrix eigenvalues!");
			return -1;
		}

		INFO_LOG("Printing Sigma matrix for component no. "<<k+1<<"...");
		fSigma[k].Print();
			
		if(MathUtils::ComputeSymMatrixEigenvalues(fSigmaEigen[k],eigenVects,fSigma[k])<0){
			WARN_LOG("Failed to compute covariance matrix eigenvalues!");
			return -1;
		}
	}//end loop components

	//## Set safe values
	for(int k=0;k<fNComponents;k++){
		fKsi_safe[k]= fKsi[k];
		fMu_safe[k]= fMu[k];
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];
		fOmega_safe[k]= fOmega[k];
		//fOmegaEigen_safe[k]= fOmegaEigen[k];
		fDelta_safe[k]= fDelta[k];
		fNu_safe[k]= fNu[k];
		fP_safe[k]= fP[k];
	}//end loop components


	//Print parss
	INFO_LOG("Print starting pars...");
	PrintPars();

	return 0;

}//close InitializeFitPars()

int MSTMixtureFitter::RunEM_Init()
{		
	//## Init EM algorithm data	
	for(size_t i=0;i<fData.size();i++){
		fTau.push_back( std::vector<double>() );
		fE1.push_back( std::vector<double>() );
		fE2.push_back( std::vector<double>() );
		fE3.push_back( std::vector<double>() );
		fE4.push_back( std::vector<double>() );
	
		for(int k=0;k<fNComponents;k++){
			fTau[i].push_back(0.);
			fE1[i].push_back(0.);
			fE2[i].push_back(0.);
			fE3[i].push_back(0.);
			fE4[i].push_back(0.);
		}
	}//end loop data
	
	//## Initialize starting pars
	if(InitializeFitPars()<0){
		ERROR_LOG("Failed to initialize fit starting parameters!");
		return -1;
	}

	
	return 0;

}//close RunEM_Init()

int MSTMixtureFitter::RunEM_EStep(double& LL)
{
	INFO_LOG("Running EM E-step...");

	double beta= 1;//annealing term (no annealing set)

	//############################################################
	//## E step: Compute the expectation e0=tau, e1, e2, e3, e4
	//############################################################
	//Compute the inverse of the covariance matrix for each mixture
	for(int k=0;k<fNComponents;k++){
		//Compute Sigma inverse
		fSigmaInv[k]= TMatrixD(TMatrixD::kInverted, fSigma[k]);
		fSigmaDet[k]= fSigma[k].Determinant();	
  	if (fSigmaDet[k]<=0) {
			WARN_LOG("Covariance matrix inversion failed!");
		}

		//Compute Scale matrix inverse
		fOmegaInv[k]= TMatrixD(TMatrixD::kInverted, fOmega[k]);
		fOmegaDet[k]= fOmega[k].Determinant();	
  	if (fOmegaDet[k]<=0) {
			WARN_LOG("Scale matrix inversion failed!");
		}
	}//end loop mixtures


	//## Compute posterior probability & LL
	LL= 0.;
	double likelihood= 0.;

	for(long int i=0;i<fN;i++){
		double tauSum= 0.;	
		if(i%1000==0) DEBUG_LOG("--> "<<i<<" events processed...");
		
		for(int k=0;k<fNComponents;k++)
		{
			//## Compute a, b, c, C, D, R
			double a= MathUtils::MST_a( fOmegaInv[k],fDelta[k]);
			double b= MathUtils::MST_b( fData[i],fOmegaInv[k],fDelta[k],fKsi[k]);
			double c= MathUtils::MST_c( fData[i],fOmegaInv[k],fKsi[k],fNu[k]);
			double C= MathUtils::MST_C( fNu[k],fNDim,fOmegaDet[k]);
			double D= MathUtils::MST_D( fData[i],fOmegaInv[k],fKsi[k],fDelta[k],fNu[k]);
			double R= MathUtils::MST_R( fNu[k],fNDim);
			
			//## Compute e1, e2, e3, e4, tau
			double e1= MathUtils::MST_e1(a,b,D,R);
			double e2= MathUtils::MST_e2(a,b,D,R);
			double e3= MathUtils::MST_e3(a,b,c,D,R);
			double e4= MathUtils::MST_e4(a,b,D,R);
			double tauComponent= MathUtils::MST_tau(a,b,C,D,R,fP[k],beta);
				
			fE1[i][k]= e1;
			fE2[i][k]= e2;
			fE3[i][k]= e3;
			fE4[i][k]= e4;	
			fTau[i][k]= tauComponent;
			tauSum+= tauComponent;	
	
		}//end loop mixture components
			 
		LL+= log(tauSum);

		for(int k=0;k<fNComponents;k++) {
			if(tauSum>0) fTau[i][k]/= tauSum;
			else{
				ERROR_LOG("Negative or null tau sum ("<<tauSum<<")...exit!");
				return -1;
			}
		}//end loop mixture components		
	}//end loop events
	
	return 0;

}//close RunEM_EStep()

int MSTMixtureFitter::RunEM_MStep()
{
	INFO_LOG("Running EM M-step...");

	//############################################################
	//## M step: Update the fit parameters
	//############################################################
	for(int k=0;k<fNComponents;k++)
	{		
		double tauSum= 0.;
		TMatrixD ksiSum(1,fNDim);
		TMatrixD omegaSum(fNDim,fNDim);
		TMatrixD deltaSum(1,fNDim);
			
		TMatrixD diff(1,fNDim);
		TMatrixD diff_t(fNDim,1);
		TMatrixD delta_t(fNDim,1);
			
		double ksiNorm= 0.;
		double omegaNorm= 0.;
		double deltaNorm= 0.;
		double nuTerm= 0.;

		for(int i=0;i<fN;i++){	
			diff= (fData[i]-fKsi[k]);			
			diff_t= TMatrixD(TMatrixD::kTransposed, diff);
			delta_t= TMatrixD(TMatrixD::kTransposed, fDelta[k]);

			tauSum+= fTau[i][k];
			ksiSum+= (fData[i]*fE1[i][k]-fDelta[k]*fE2[i][k])*fTau[i][k];
			omegaSum+= (diff_t*diff*fE1[i][k] - fE2[i][k]*delta_t*diff - diff_t*fDelta[k]*fE2[i][k] + fE3[i][k]*delta_t*fDelta[k])*fTau[i][k];	
			deltaSum+= fTau[i][k]*fE2[i][k]*diff;

			ksiNorm+= fTau[i][k]*fE1[i][k];
			omegaNorm+= fTau[i][k];
			deltaNorm+= fTau[i][k]*fE3[i][k];
			nuTerm+= fTau[i][k]*(fE4[i][k]-fE1[i][k]);
		}//end loop events
	

		//## Update fit parameters
		//- Update fractions
		if(!fFixFractionPars) fP[k]= tauSum/(double)(fN);

		for(int j=0;j<fNDim;j++){

			//Update means
			if(!fFixMeanPars) {	
				(fKsi[k])(0,j)= ksiSum(0,j)/ksiNorm;					
			}

			//Update delta
			if(!fFixDeltaPars) {
				(fDelta[k])(0,j)= deltaSum(0,j)/deltaNorm;
			}
				
			//Update covariance				
			if(!fFixCovariancePars){
				for(int l=0;l<fNDim;l++){
					(fOmega[k])(j,l)= omegaSum(j,l)/omegaNorm;
				}	//end loop dim
			}//close if 
		}//end loop dim
			
		if(!fFixNuPars) {
			double previousNu= fNu[k];
			double nuRoot= 0;	
			int solverStatus= MathUtils::NuSolver(nuRoot,tauSum,nuTerm,previousNu);
			if(solverStatus<0){
				WARN_LOG("Failed to solve nu equation, setting nu to previous value!");
				fNu[k]= previousNu;
			}
			else{
				fNu[k]= nuRoot;
			}
		}
			
		//## Compute mean and covariance
		MathUtils::KsiToMu(fMu[k],fKsi[k],fDelta[k],fNu[k]);
		MathUtils::ScaleToCovarianceMatrix(fSigma[k],fOmega[k],fDelta[k],fNu[k]);
		
	}//end loop mixture components
		
	//## Check covariance matrix integrity
	INFO_LOG("Check covariance & scale matrix integrity after update...");
	if(!fFixCovariancePars && CheckCovariance()<0){
		WARN_LOG("Failed to perform checks in covariance matrix!");
		return -1;
	}
		
	INFO_LOG("Printing pars after EM M-step ...");
	PrintPars();
	

	return 0;

}//close RunEM_MStep()


int MSTMixtureFitter::CheckCovariance()
{
	//## Loop over components and check covariance matrix:
	//##   - must be symmetric
	//##   - must be positive def	
	//##   - set diagonal (if enabled)

	for(int k=0;k<fNComponents;k++)
	{
		//## Force diagonal covariance?
		if(fForceDiagonalCovariance){
			MathUtils::MakeDiagonalMatrix(fOmega[k]);
		}

		//## Force covariance to be symmetric and pos def
		//## Approximate to the nearest covariance matrix
		if(MathUtils::MakeSymmPosDefCovarianceMatrix(fOmega[k])<0){
			WARN_LOG("Failed to correct scale matrix for component no. "<<k+1<<"!");
			return -1;	
		}

		//## Update inverse & determinant	
		INFO_LOG("Updating scale matrix inverse...");
		fOmegaInv[k]= TMatrixD(TMatrixD::kInverted,fOmega[k]);
		fOmegaDet[k]= fOmega[k].Determinant();	
		if(fOmegaDet[k]<=1.e-16){
			WARN_LOG("Scale matrix for component "<<k+1<<" is not positive def (det="<<fOmegaDet[k]<<") ... setting to safe matrix!");
			fOmega[k]= fOmega_safe[k];
		}

		//## Update covariance matrix
		MathUtils::ScaleToCovarianceMatrix(fSigma[k],fOmega[k],fDelta[k],fNu[k]);
		fSigmaInv[k]= TMatrixD(TMatrixD::kInverted,fSigma[k]);
		fSigmaDet[k]= fSigma[k].Determinant();	
		
		//## Update eigenvalues & eigenvectors
		INFO_LOG("Update eigenvalues & eigenvectors of scale matrix ...");	
		fOmega[k].Print();
		TMatrixD eigenVects(fNDim,fNDim); 
		if(MathUtils::ComputeSymMatrixEigenvalues(fOmegaEigen[k],eigenVects,fOmega[k])<0){
			ERROR_LOG("Failed to compute covariance matrix eigenvalues/eigenvectors!");
			return -1;
		}

		INFO_LOG("Update eigenvalues & eigenvectors of covariance matrix ...");	
		fSigma[k].Print();
		if(MathUtils::ComputeSymMatrixEigenvalues(fSigmaEigen[k],eigenVects,fSigma[k])<0){
			ERROR_LOG("Failed to compute covariance matrix eigenvalues/eigenvectors!");
			return -1;
		}

		//## Set the current covariance as "safe"
		INFO_LOG("Set the current covariance as safe...");	
		fOmega_safe[k]= fOmega[k];
		fOmegaEigen_safe[k]= fOmegaEigen[k];
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];
		
	}//end loop components

	return 0;

}//close CheckCovariance()


int MSTMixtureFitter::InitParsToKMeans()
{
	//####################################
	//##  METHOD TO INIT FROM KMEANS
	//####################################
	//- Find clusters using kmeans (multiple starting points)
	//- Compute cluster sampling means, covariance & skewness
	//- Generate initial skew-t values using the method from Lee, McLachlan "Finite mixtures of multivariate skew t-distributions: some recent and new results" (sec. 5.1.3)
	//      Omega= CovMatrix + (a-1) diag(CovMatrix)

	double a= gRandom->Uniform(0,1);
	double nu_start= 40;//not justified in the papers

	INFO_LOG("Initializing mixture parameters with k-means...");

	//## Perform Kmeans clustering on pre-completed data
	KMeansClustering kmeans;
	if(kmeans.RunKMedians(fRTableName,fNComponents)<0){
		ERROR_LOG("KMeans clustering run failed!");
		return -1;
	}

	//## Retrieve data
	TMatrixD* clusterSizes= kmeans.GetClusterSizes();
	TMatrixD* clusterCenters= kmeans.GetClusterCenters();
	std::vector<TMatrixD*> clusterCovMatrixes= kmeans.GetClusterCovMatrixes();
	TMatrixD* clusterSkewness= kmeans.GetClusterSkewness();
	DEBUG_LOG("clusterSizes: "<<clusterSizes->GetNrows()<<" x "<<clusterSizes->GetNcols());
	DEBUG_LOG("clusterCenters: "<<clusterCenters->GetNrows()<<" x "<<clusterCenters->GetNcols());
		
	//## Compute cluster cov diagonal matrix
	std::vector<TMatrixD> CovMatrixDiag;
	std::vector<TMatrixD> CovMatrixDiagSqrRoot;
	for(int k=0;k<fNComponents;k++){
		CovMatrixDiag.push_back(TMatrixD(fNDim,fNDim));
		CovMatrixDiagSqrRoot.push_back(TMatrixD(fNDim,fNDim));
		for(int j=0;j<fNDim;j++){
			double w= (*clusterCovMatrixes[k])(j,j);	
			CovMatrixDiag[k](j,j)= w;
			CovMatrixDiagSqrRoot[k](j,j)= sqrt(w);
		}
	}

	//## Assign cluster weights as fraction start
	INFO_LOG("Assign the start fraction parameters according to kmeans cluster weights...");
	for(int k=0;k<fNComponents;k++){
		double clusterSize= (*clusterSizes)(k,0);
		fP[k]= clusterSize/(double)(fN);
	}//end loop components

	//## Assign nu start
	for(int k=0;k<fNComponents;k++){
		fNu[k]= nu_start;
	}//end loop components

	//## Assign delta start
	for(int k=0;k<fNComponents;k++){	
		for(int j=0;j<fNDim;j++){
			double skewness= (*clusterSkewness)(k,j);
			fDelta[k](0,j)= MathUtils::sign(skewness)*sqrt( (TMath::Pi()*(1.-a))/(TMath::Pi()-2.) ) * CovMatrixDiagSqrRoot[k](j,j);
		}
	}//end loop components


	//## Assign start scale matrix and covariance matrix
	INFO_LOG("Assign the start covariance parameters according to the kmeans cluster variances...");
	for(int k=0;k<fNComponents;k++){
		//fOmega[k]= (*clusterCovMatrixes[k]) - (a-1)*CovMatrixDiag[k];//in a paper the sign is minus
		fOmega[k]= (*clusterCovMatrixes[k]) + (a-1)*CovMatrixDiag[k];//in another paper the sign is plus
		MathUtils::ScaleToCovarianceMatrix(fSigma[k],fOmega[k],fDelta[k], fNu[k]);	
	}//end loop components
	
	
	//## Assign cluster centers as mean start
	INFO_LOG("Assign the start mean parameters according to the kmeans cluster centers ...");
	for(int k=0;k<fNComponents;k++){
		for(int j=0;j<fNDim;j++){
			fKsi[k](0,j)= (*clusterCenters)(k,j) - sqrt(2./TMath::Pi())*fDelta_start[k](0,j);  	
		}//end loop dim
		MathUtils::KsiToMu(fMu[k],fKsi[k],fDelta[k],fNu[k]);
	}//end loop components

	
	//## Print pars
	INFO_LOG("Printing par start values after kmeans...");
	PrintPars();

	return 0;

}//close InitParsToKMeans()


int MSTMixtureFitter::InitParsToUser()
{
	INFO_LOG("Initializing mixture parameters with user defaults...");

	//Check fractions weights and set to starting values
	int nFractionPars= static_cast<int>(fP_start.size());
	if(nFractionPars!=fNComponents){
		ERROR_LOG("Number of components weights given as arg ("<<nFractionPars<<") is different from nComponents ("<<fNComponents<<")!");
		return -1;
	}
	for(int k=0;k<fNComponents;k++){
		fP[k]= fP_start[k];
	}

	//Check nu pars and set to starting values
	int nNuPars= static_cast<int>(fNu_start.size());
	if(nNuPars!=fNComponents){
		ERROR_LOG("Number of components nu given as arg ("<<nNuPars<<") is different from nComponents ("<<fNComponents<<")!");
		return -1;
	}
	for(int k=0;k<fNComponents;k++){
		fNu[k]= fNu_start[k];
	}
	
	//Check component means
	INFO_LOG("Print start means...");
	for(size_t k=0;k<fMu_start.size();k++){
		fMu_start[k].Print();
	}

	int nMeanComponents= static_cast<int>(fMu_start.size());
	if(nMeanComponents != fNComponents){
		ERROR_LOG("Number of user mean components ("<<nMeanComponents<<") is different from nComponents ("<<fNComponents<<")!");
		return -1;
	}
	for(size_t k=0;k<(fMu_start).size();k++){
		int meanVectDim= (fMu_start)[k].GetNcols();
		if(meanVectDim!=fNDim){
			ERROR_LOG("User mean vector size for component no. "<<k+1<<" is not equal to nDim="<<fNDim<<"!");
			return -1;
		}
		fMu[k]= (fMu_start)[k];
	}//end loop components

	//Check component delta
	INFO_LOG("Print start delta...");
	for(size_t k=0;k<fDelta_start.size();k++){
		fDelta_start[k].Print();
	}

	int nDeltaComponents= static_cast<int>(fDelta_start.size());
	if(nDeltaComponents != fNComponents){
		ERROR_LOG("Number of user delta par components ("<<nDeltaComponents<<") is different from nComponents ("<<fNComponents<<")!");
		return -1;
	}
	for(size_t k=0;k<(fDelta_start).size();k++){
		int deltaVectDim= (fDelta_start)[k].GetNcols();
		if(deltaVectDim!=fNDim){
			ERROR_LOG("User delta vector size for component no. "<<k+1<<" is not equal to nDim="<<fNDim<<"!");
			return -1;
		}
		fDelta[k]= (fDelta_start)[k];
	}//end loop components
		
	//Check sigmas
	int nSigmaComponents= static_cast<int>(fSigma_start.size());
	if(nSigmaComponents != fNComponents){
		ERROR_LOG("Number of user sigma components ("<<nSigmaComponents<<") is different from nComponents ("<<fNComponents<<")!");
		return -1;
	}
	for(size_t k=0;k<(fSigma_start).size();k++){
		int nRows= (fSigma_start)[k].GetNrows();
		int nCols= (fSigma_start)[k].GetNcols();
		if(nRows!=fNDim || nCols!=fNDim){
			ERROR_LOG("User sigma matrix for component no. "<<k+1<<" has size different from nDim x nDim (nDim="<<fNDim<<")!");
			return -1;
		}
		fSigma[k]= (fSigma_start)[k];
	}//end loop components

	//Computing ksi & scale matrix from mean and covariance
	INFO_LOG("Computing location vector and scale matrix from mean & covariance ...");
	for(int k=0;k<fNComponents;k++){
		MathUtils::MuToKsi(fKsi[k],fMu[k],fDelta[k],fNu[k]);
		MathUtils::CovarianceToScaleMatrix(fOmega[k],fSigma[k],fDelta[k],fNu[k]);
	}

	PrintPars();

	return 0;

}//close InitParsToUser()



/*
void MSTMixtureFitter::EMRun(double beta){

	//## Start EM iteration loop
	fFitStatus= eInit;
	int badIterCounter= 0;
	int restartCounter= 0;
	fFixMeanDiffPar= false;
	fConstraintFirstMeanPar= false;
	
	for(int iter=0;iter<fNIterations;iter++){
		cout<<"*** ITER NO. "<<iter<<" ***"<<endl;		
		
		//############################################################
		//## E step: Compute the expectation e0=tau, e1, e2, e3, e4
		//############################################################
		cout<<"--> E Step"<<endl;
		//Compute the inverse of the covariance matrix for each mixture
		for(int k=0;k<fNComponents;k++){
			double SigmaDet= 0.;
			TMatrixD inverse(fNDim,fNDim);
			TMatrixD m(fNDim,fNDim);
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++){
					m(j,l)= (fSigma[k])(j,l);
				}			
			}
			
			inverse = m.Invert(&SigmaDet);
			
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++){
					(fSigmaInv[k])(j,l)= inverse(j,l);
				}			
			}

			fSigmaDet[k]= SigmaDet;	
  		if (SigmaDet<=0) {
				cerr<<"MSTMixtureFitter::EMFitter(): WARNING: Covariance matrix inversion failed!"<<endl;
			}
			
			cout<<"== Component "<<k+1<<" =="<<endl;
			cout<<"p= "<<fP[k]<<endl;
			cout<<"Nu= "<<fNu[k]<<endl;
			cout<<"Ksi= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fKsi[k])(j,0)<<",";
			cout<<(fKsi[k])(fNDim-1,0)<<")"<<endl; 
			cout<<"KsiDiff= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fKsiDiff[k])(j,0)<<",";
			cout<<(fKsiDiff[k])(fNDim-1,0)<<")"<<endl; 
			
			cout<<"delta= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fDelta[k])(j,0)<<",";
			cout<<(fDelta[k])(fNDim-1,0)<<")"<<endl; 
			
			cout<<"Sigma= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
			}	
			cout<<")"<<endl;

			cout<<"V= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fV[k])(j,l)<<",";
			}	
			cout<<")"<<endl;

			cout<<"SigmaInv= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fSigmaInv[k])(j,l)<<",";
			}	
			cout<<")"<<endl;

			cout<<"SigmaDet= "<<fSigmaDet[k]<<endl;
			cout<<"SigmaEigen= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fSigmaEigen[k])(j,l)<<",";
			}	
			cout<<")"<<endl;
			
		}//end loop mixtures


		double likelihood= 0.;

		for(int i=0;i<fN;i++){
			
			double tauSum= 0.;
			
			if(i%1000==0) cout<<"--> "<<i<<" events processed..."<<endl;

			
			for(int k=0;k<fNComponents;k++){

				//## Compute a, b, c, C, D, R
				double a= MathUtilities::a( fSigmaInv[k],fDelta[k]);
				double b= MathUtilities::b( fData[i],fSigmaInv[k],fDelta[k],fKsi[k]);
				double c= MathUtilities::c( fData[i],fSigmaInv[k],fKsi[k],fNu[k]);
				double C= MathUtilities::C( fNu[k],fNDim,fSigmaDet[k]);
				double D= MathUtilities::D( fData[i],fSigmaInv[k],fKsi[k],fDelta[k],fNu[k]);
				double R= MathUtilities::R( fNu[k],fNDim);
			
				//cout<<"Component "<<k+1<<" => a="<<a<<"  b="<<b<<endl;	

				//cout<<"Component "<<k+1<<" => (a,b,c,C,D,R)= ("<<a<<","<<b<<","<<c<<","<<C<<","<<D<<","<<R<<")"<<endl;	

				
				//## Compute e1, e2, e3, e4, tau
				double e1= MathUtilities::e1(a,b,D,R);
				double e2= MathUtilities::e2(a,b,D,R);
				double e3= MathUtilities::e3(a,b,c,D,R);
				double e4= MathUtilities::e4(a,b,D,R);
				double tauComponent= MathUtilities::tau(a,b,C,D,R,fP[k],beta);
				
				fE1[i][k]= e1;
				fE2[i][k]= e2;
				fE3[i][k]= e3;
				fE4[i][k]= e4;	
				fTau[i][k]= tauComponent;
				tauSum+= tauComponent;	

				//cout<<"Component "<<k+1<<" => (e1,e2,e3,e4,tau)= ("<<e1<<","<<e2<<","<<e3<<","<<e4<<","<<tauComponent<<")"<<endl;	

			}//end loop mixture components
			 
			likelihood+= log(tauSum);

			for(int k=0;k<fNComponents;k++) {
				if(tauSum>0) fTau[i][k]/= tauSum;
				else{
					cerr<<"MSTMixtureFitter::EMFitter(): ERROR: Negative or null tau sum ("<<tauSum<<")...exit!"<<endl;	
					exit(1);
				}
			}//end loop mixture components
			
		}//end loop events
	
		//## If deltaLogL>0 set current parameters to best
		double DeltaLogLikelihood= likelihood-fLogLikelihood;
		if(DeltaLogLikelihood>0){
			for(int k=0;k<fNComponents;k++){
				fP_best[k]= fP[k]; 
				fKsi_best[k]= fKsi[k];
				fKsiDiff_best[k]= fKsiDiff[k];
				fMu_best[k]= fMu[k];
				fMuDiff_best[k]= fMuDiff[k];
				fSigma_best[k]= fSigma[k];
				fV_best[k]= fV[k];
				fDelta_best[k]= fDelta[k];
				fNu_best[k]= fNu[k];
				fLikelihood_best= likelihood;
				badIterCounter= 0;
			}//end loop components
		}
		else{
			fFitStatus= eTroubles;
			badIterCounter++;
		}

		cout<<"==> LogLikelihood= "<<likelihood<<"  DeltaL="<<DeltaLogLikelihood<<endl;
		fLogLikelihood= likelihood;
		fIterNo= iter;
		//fFitIterInfo->Fill();

		fIterLogLikelihood[iter]= fLogLikelihood;
		if(iter>3){
			double ak= (fIterLogLikelihood[iter]-fIterLogLikelihood[iter-1])/(fIterLogLikelihood[iter-1]-fIterLogLikelihood[iter-2]);
			double LogL_inf= fIterLogLikelihood[iter-1] + (fIterLogLikelihood[iter]-fIterLogLikelihood[iter-1])/(1.-ak);
			double DeltaLogL= LogL_inf-fIterLogLikelihood[iter];
			cout<<"==> ak="<<ak<<"  LogLikelihood_inf= "<<LogL_inf<<"  DeltaLogL="<<DeltaLogL<<endl;

			if(fUseStoppingCriteria && DeltaLogL>=0 && fabs(DeltaLogL)<fEpsilon){
				cout<<"MSTMixtureFitter::EMFitter(): INFO: Stop criteria matched (DeltaLogL<eps)...exit iteration!"<<endl;
				fFitStatus= eConverged;
				break;
			}
		}

		//############################################################
		//## M step: Update the fit parameters
		//############################################################
		for(int k=0;k<fNComponents;k++){
		
			double tauSum= 0.;
			TMatrixD ksiSum(fNDim,1);
			ksiSum.Zero();
			TMatrixD sigmaSum(fNDim,fNDim);
			sigmaSum.Zero();
			TMatrixD deltaSum(fNDim,1);
			deltaSum.Zero();

			TMatrixD diff(fNDim,1);
			diff.Zero();
			TMatrixD diff_t(1,fNDim);
			diff_t.Zero();
			TMatrixD delta_t(1,fNDim);
			delta_t.Zero();

		
			double ksiNorm= 0.;
			double sigmaNorm= 0.;
			double deltaNorm= 0.;
			double nuTerm= 0.;

			for(int i=0;i<fN;i++){
				
				diff= (fData[i]-fKsi[k]);			
				diff_t= TMatrixD(TMatrixD::kTransposed, diff);
				delta_t= TMatrixD(TMatrixD::kTransposed, fDelta[k]);

				tauSum+= fTau[i][k];
				ksiSum+= (fData[i]*fE1[i][k]-fDelta[k]*fE2[i][k])*fTau[i][k];
				sigmaSum+= (diff*diff_t*fE1[i][k] - fE2[i][k]*fDelta[k]*diff_t - diff*delta_t*fE2[i][k] + fE3[i][k]*fDelta[k]*delta_t)*fTau[i][k];	
				//sigmaSum+= (diff*diff_t*fE1[i][k] - fE2[i][k]*fDelta[k]*diff_t*diff*delta_t + fE3[i][k]*fDelta[k]*delta_t)*fTau[i][k];	
				deltaSum+= fTau[i][k]*fE2[i][k]*diff;

				ksiNorm+= fTau[i][k]*fE1[i][k];
				sigmaNorm+= fTau[i][k];
				deltaNorm+= fTau[i][k]*fE3[i][k];
				nuTerm+= fTau[i][k]*(fE4[i][k]-fE1[i][k]);

			}//end loop events
	

			//## Update fit parameters
			if(!fFixFractionPar) fP[k]= tauSum/(double)(fN);

			for(int j=0;j<fNDim;j++){

				if(!fFixMeanPar) {
					
					//(fKsi[k])(j,0)= ksiSum(j,0)/ksiNorm;
						
					if(k==0){
						fKsi[k](j,0)= ksiSum(j,0)/ksiNorm;
						fKsiDiff[k](j,0)= 0;
					}
					else{
						if(!fFixMeanDiffPar) fKsiDiff[k](j,0)= ksiSum(j,0)/ksiNorm-fKsi_safe[0](j,0);
						fKsi[k](j,0)= fKsi[0](j,0)+fKsiDiff[k](j,0);
	
						//if(!fFixMeanDiffPar) fKsiDiff[k](j,0)= ksiSum(j,0)/ksiNorm-fKsi[0](j,0);
						//fKsi[k](j,0)= fKsi_safe[0](j,0)+fKsiDiff[k](j,0);
					}
					
				}//close if !fFixMeanPar


				if(!fFixDeltaPar) (fDelta[k])(j,0)= deltaSum(j,0)/deltaNorm;
								
				if(!fFixCovariancePar){
					for(int l=0;l<fNDim;l++){
						(fSigma[k])(j,l)= sigmaSum(j,l)/sigmaNorm;
					}	//end loop dim
				}//close if 
			}//end loop dim
			
			if(!fFixNuPar) {
				//MathUtilities::SetNuPar1(tauSum);
				//MathUtilities::SetNuPar2(nuTerm);
				double previousNu= fNu[k];	
				fNu[k]= MathUtilities::NuSolver(tauSum,nuTerm,previousNu);
			}
			
			//## Compute mean and covariance
			fMu[k]= MathUtilities::KsiToMu(fKsi[k],fDelta[k],fNu[k]);
			fMuDiff[k]= fMu[k]-fMu[0];
			fV[k]= MathUtilities::ScaleMatrixToCovarianceMatrix(fSigma[k],fDelta[k],fNu[k]);
			fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);
			fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);
				
		}//end loop mixture components
		
		
	
		cout<<"MSTMixtureFitter::EMFitter(): INFO: After ordering the groups..."<<endl;
		for(int k=0;k<fNComponents;k++){	
			cout<<"== Component "<<k+1<<" =="<<endl;
			cout<<"p= "<<fP[k]<<endl;
			cout<<"Nu= "<<fNu[k]<<endl;
			cout<<"Ksi= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fKsi[k])(j,0)<<",";
			cout<<(fKsi[k])(fNDim-1,0)<<")"<<endl; 
			cout<<"Mu= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(j,0)<<",";
			cout<<(fMu[k])(fNDim-1,0)<<")"<<endl; 

			cout<<"KsiDiff= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fKsiDiff[k])(j,0)<<",";
			cout<<(fKsiDiff[k])(fNDim-1,0)<<")"<<endl; 

			cout<<"Mu= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(j,0)<<",";
			cout<<(fMu[k])(fNDim-1,0)<<")"<<endl; 
			cout<<"MuDiff= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fMuDiff[k])(j,0)<<",";
			cout<<(fMuDiff[k])(fNDim-1,0)<<")"<<endl; 
			
			cout<<"delta= (";
			for(int j=0;j<fNDim-1;j++) cout<<(fDelta[k])(j,0)<<",";
			cout<<(fDelta[k])(fNDim-1,0)<<")"<<endl; 
			
			cout<<"Sigma= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
			}	
			cout<<")"<<endl;
			cout<<"V= (";
			for(int j=0;j<fNDim;j++){
				for(int l=0;l<fNDim;l++) cout<<(fV[k])(j,l)<<",";
			}	
			cout<<")"<<endl;
		}//end loop components
	
		//############################################################
		//## Constrain Step
		//############################################################
		//if(fUseConstraint) EMConstrain();	
		if(fUseConstraint) EMConstrain2();
			
		//## Store final parameters to file
		for(int k=0;k<fNComponents;k++){
			fFraction[k]= fP_best[k];
			fXmaxMean[k]= fKsi_best[k](0,0);
			if(fNDim>=2) fNmuMean[k]= fKsi_best[k](1,0);
			fXmaxVar[k]= fSigma_best[k](0,0);
			if(fNDim>=2) {
				fNmuVar[k]= fSigma_best[k](1,1);
				fXmaxNmuCov[k]= fSigma_best[k](0,1);
			}

			fXmaxDelta[k]= fDelta_best[k](0,0);
			if(fNDim>=2) fNmuDelta[k]= fDelta_best[k](1,0);
			fNdf[k]= fNu_best[k];			
		}

		fFitIterInfo->Fill();

		//#####################
		//### CHECK BAD ITER
		//#####################
		if(badIterCounter>=20){
			fConstraintFirstMeanPar= true;
			restartCounter++;

			if(restartCounter>=5){//## constrain all mean parameters starting from best optimum
				for(int k=0;k<fNComponents;k++){	
					for(int j=0;j<fNDim;j++){
						fKsi[k](j,0)= fKsi_best[k](j,0);
						fKsiDiff[k](j,0)= fKsiDiff_best[k](j,0);
					}	
				}

			}//close if
			else{//## try again regenerating mean diff to escape from stuck
				for(int j=0;j<fNDim;j++){
					fKsi[0](j,0)= fKsi_best[0](j,0);
				
					for(int k=1;k<fNComponents;k++){	
						double minBoundary= fKsiDiff_min[k](j,0);
						double maxBoundary= fKsiDiff_max[k](j,0);
						fKsiDiff[k](j,0)= gRandom->Uniform(minBoundary,maxBoundary);
						fKsi[k](j,0)= fKsi[0](j,0)+fKsiDiff[k](j,0);
					}//end loop components
				}//end loop dim
			}//close else
		
		}//close if

		//## Recompute again mean and covariance
		for(int k=0;k<fNComponents;k++){
			fMu[k]= MathUtilities::KsiToMu(fKsi[k],fDelta[k],fNu[k]);
			fMuDiff[k]= fMu[k]-fMu[0];
			fV[k]= MathUtilities::ScaleMatrixToCovarianceMatrix(fSigma[k],fDelta[k],fNu[k]);
			fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);
			fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);
		}

	}//end loop iterations


	//## Run numerical fitter?
	if(fRunMinuitFitterAtConvergence){
		INFO_LOG("Running MINUIT fitter after EM convergence...");
		RunFitter(false);
	}

}//close MSTMixtureFitter::EMRun()
*/


int MSTMixtureFitter::RunEM_ConstrainStep()
{
	INFO_LOG("Running EM constraint step...");

	//## Check if current parameters fullfil the group constraints
	cout<<"======================"<<endl;
	cout<<"== EMConstrain step =="<<endl;
	cout<<"======================"<<endl;

	/*	
	//############################
	//##   DELTA CONSTRAINTS    ##
	//############################
	cout<<"==> DELTA CONSTRAINT"<<endl;
	bool isBadDelta_minmax[fNComponents];
	bool isDeltaBoundConstraintViolated= false;
	std::vector<TMatrixD> alphaDeltaList_minmax;
	alphaDeltaList_minmax.clear();
	alphaDeltaList_minmax.resize(0);

	if(!fFixDeltaPar && fUseDeltaBoundConstraint){

		//## Check delta constraint violations
		for(int k=0;k<fNComponents;k++){
			isBadDelta_minmax[k]= false;

			TMatrixD alphaDelta_minmax(fNDim,1);
			alphaDelta_minmax.Zero();
	
			for(int j=0;j<fNDim;j++){
				alphaDelta_minmax(j,0)= 1;
				
				double constraintSign= 1;
					
				double Delta= fDelta[k](j,0);
				double safeDelta= fDelta_safe[k](j,0);
				double Deltamin= fDelta_min[k](j,0);
				double Deltamax= fDelta_max[k](j,0);
				double alphaDelta_minBound= (Deltamin-safeDelta)/(Delta-safeDelta);
				double alphaDelta_maxBound= (Deltamax-safeDelta)/(Delta-safeDelta);
		
				if(constraintSign*Delta<constraintSign*Deltamin){//min constraint violated
					alphaDelta_minmax(j,0)= alphaDelta_minBound;
					isBadDelta_minmax[k]= true;
					cout<<"INFO: Delta min constraint violated for component no. "<<k+1<<":  Delta(k)="<<Delta<<"  Deltamin(k)="<<Deltamin<<"  safeDelta(k)="<<safeDelta<<"  alpha="<<alphaDelta_minBound<<endl;
				}
				if(constraintSign*Delta>constraintSign*Deltamax){//max constraint violated
					alphaDelta_minmax(j,0)= alphaDelta_maxBound;
					isBadDelta_minmax[k]= true;
					cout<<"INFO: Delta max constraint violated for component no. "<<k+1<<":  Delta(k)="<<Delta<<"  Deltamax(k)="<<Deltamax<<"  safeDelta(k)="<<safeDelta<<"  alpha="<<alphaDelta_maxBound<<endl;
				}	
			}//end loop dim

			alphaDeltaList_minmax.push_back(alphaDelta_minmax);

		}//end loop components
	
	
		//## Update Delta values
		cout<<"INFO: Updating delta parameter..."<<endl;			
		std::vector<bool> hasToBeGenerated;
		hasToBeGenerated.assign(fNComponents,false);
		bool hasToRegenerateDelta= false;

		for(int k=0;k<fNComponents;k++){
			if(!isBadDelta_minmax[k]) continue;

			double alphaDeltaMin_minmax= alphaDeltaList_minmax[k].Min();
			double alphaDeltaOpt= alphaDeltaMin_minmax/fConstraintAlphaScale;
		
			cout<<"INFO: Component no. "<<k+1<<": alphaDeltaMin="<<alphaDeltaMin_minmax<<"  alphaDeltaOpt="<<alphaDeltaOpt<<endl;

			if(alphaDeltaOpt<fConstraintAlphaTolerance){
				hasToRegenerateDelta= true;
				hasToBeGenerated[k]= true;	
				cout<<"INFO: Delta par is stuck in bound constraint ("<<alphaDeltaOpt<<"<"<<fConstraintAlphaTolerance<<")"<<endl;	
			}
			
			//## Update delta values
			for(int j=0;j<fNDim;j++){
				double Delta= fDelta[k](j,0);
				double safeDelta= fDelta_safe[k](j,0);
				fDelta[k](j,0)= (1.-alphaDeltaOpt)*safeDelta + alphaDeltaOpt*Delta;
			}//end loop dim
			
		}//end loop components

		if(hasToRegenerateDelta && fUseRandomRegenerationAfterStuck){
			RandomizeDeltaPar(hasToBeGenerated);	
		}

	}//close if	(!fFixDeltaPar && fUseDeltaBoundConstraint)
	//########################################################



	//##########################
	//##   NU CONSTRAINTS    ##
	//##########################
	cout<<endl;
	cout<<"==> NU CONSTRAINT"<<endl;
	bool isBadNu_minmax[fNComponents];
	bool isNuBoundConstraintViolated= false;
	std::vector<double> alphaNuList_minmax;
	alphaNuList_minmax.clear();
	alphaNuList_minmax.resize(0);

	if(!fFixNuPar && fUseNuBoundConstraint){
		double constraintSign= 1;

		for(int k=0;k<fNComponents;k++){
			isBadNu_minmax[k]= false;
			alphaNuList_minmax.push_back(1);

			double Nu= fNu[k];
			double safeNu= fNu_safe[k];
			double Numin= fNu_min[k];
			double Numax= fNu_max[k];
			double alphaNu_minBound= (Numin-safeNu)/(Nu-safeNu);
			double alphaNu_maxBound= (Numax-safeNu)/(Nu-safeNu);
		
			if(constraintSign*Nu<constraintSign*Numin){//min constraint violated
				alphaNuList_minmax[k]= alphaNu_minBound;
				isBadNu_minmax[k]= true;
				isNuBoundConstraintViolated= true;
				cout<<"INFO: Nu min constraint violated for component no. "<<k+1<<":  Nu(k)="<<Nu<<"  Numin(k)="<<Numin<<"  safeNu(k)="<<safeNu<<"  alpha="<<alphaNu_minBound<<endl;
			}
			if(constraintSign*Nu>constraintSign*Numax){//max constraint violated
				alphaNuList_minmax[k]= alphaNu_maxBound;
				isBadNu_minmax[k]= true;
				cout<<"INFO: Nu max constraint violated for component no. "<<k+1<<":  Nu(k)="<<Nu<<"  Numax(k)="<<Numax<<"  safeNu(k)="<<safeNu<<"  alpha="<<alphaNu_maxBound<<endl;
			}
		}//end loop components	

	
		//## Update Nu values
		cout<<"INFO: Updating nu parameter..."<<endl;	
		std::vector<bool> hasToBeGenerated;
		hasToBeGenerated.assign(fNComponents,false);
		bool hasToRegenerateNu= false;

		for(int k=0;k<fNComponents;k++){
			if(!isBadNu_minmax[k]) continue;

			double alphaNuOpt= alphaNuList_minmax[k]/fConstraintAlphaScale;
			cout<<"INFO: Component no. "<<k+1<<": alphaNuMin="<<alphaNuList_minmax[k]<<"  alphaNuOpt="<<alphaNuOpt<<endl;

			if(alphaNuOpt<fConstraintAlphaTolerance){
				hasToRegenerateNu= true;
				hasToBeGenerated[k]= true;		
				cout<<"INFO: Nu is stuck in bound constraint ("<<alphaNuOpt<<"<"<<fConstraintAlphaTolerance<<")"<<endl;		
			}
			
			//## Update nu values
			double Nu= fNu[k];
			double safeNu= fNu_safe[k];
			fNu[k]= (1.-alphaNuOpt)*safeNu + alphaNuOpt*Nu;
			

		}//end loop components

		if(hasToRegenerateNu && fUseRandomRegenerationAfterStuck){
			RandomizeNuPar(hasToBeGenerated);	
		}

	}//close if (!fFixNuPar && fUseNuBoundConstraint)
	//########################################################



	//##########################
	//##   KSI CONSTRAINTS    ##
	//##########################
	cout<<endl;
	cout<<"==> MEAN CONSTRAINT"<<endl;
	//## Recompute mean value after delta & nu update
	cout<<"INFO: Recompute mean parameter after delta & nu update..."<<endl;	
	for(int k=0;k<fNComponents;k++){
		fMu[k]= MathUtilities::KsiToMu(fKsi[k],fDelta[k],fNu[k]);
		fMuDiff[k]= fMu[k]-fMu[0];
	}
	
	//## MU DIFF CONSTRAINT
	if(!fFixMeanPar && fUseMeanDiffConstraint){

		while(true){
			bool isBadDiffMu= false;
			bool isBadDiffMu_minmax[fNComponents];
			bool isDiffMuBoundConstraintViolated= false;
			std::vector<TMatrixD> alphaDiffMuList_minmax;
			alphaDiffMuList_minmax.clear();
			alphaDiffMuList_minmax.resize(0);

			for(int k=0;k<fNComponents;k++){
				isBadDiffMu_minmax[k]= false;

				TMatrixD alphaDiffMu_minmax(fNDim,1);
				alphaDiffMu_minmax.Zero();
	
				for(int j=0;j<fNDim;j++){
					alphaDiffMu_minmax(j,0)= 1;
				
					double constraintSign= 1;
					
					double a_k= fKsiDiff_min[k](j,0);
					double b_k= fKsiDiff_max[k](j,0);
					
					double diffmu= fMuDiff[k](j,0);
					double safediffmu= fMuDiff_safe[k](j,0);
				
					double alphaMuDiff_minBound= (a_k-safediffmu)/(diffmu-safediffmu);
					double alphaMuDiff_maxBound= (b_k-safediffmu)/(diffmu-safediffmu);

					if(constraintSign*diffmu<constraintSign*a_k){//min constraint violated
						alphaDiffMu_minmax(j,0)= alphaMuDiff_minBound;
						isBadDiffMu_minmax[k]= true;
						isDiffMuBoundConstraintViolated= true;
						cout<<"INFO: Mean min diff constraint violated for component no. "<<k+1<<":  diffMu="<<diffmu<<"  diffMuMin="<<a_k<<"  safediffMu(k)="<<safediffmu<<"  alpha="<<alphaMuDiff_minBound<<endl;
					}
					
					if(constraintSign*diffmu>b_k){//max constraint violated
						alphaDiffMu_minmax(j,0)= alphaMuDiff_maxBound;
						isBadDiffMu_minmax[k]= true;
						isDiffMuBoundConstraintViolated= true;
						cout<<"INFO: Mean max diff constraint constraint violated for component no. "<<k+1<<":  diffMu="<<diffmu<<"  diffMuMax="<<b_k<<"  safediffMu(k)="<<safediffmu<<"  alpha="<<alphaMuDiff_maxBound<<endl;
					}	
					
				}//end loop dim

				alphaDiffMuList_minmax.push_back(alphaDiffMu_minmax);

			}//end loop components


			
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateMu= false;

			cout<<"INFO: Updating mu parameter..."<<endl;	
			for(int k=1;k<fNComponents;k++){
			
				double alphaDiffMuOpt= 1;
				double alphaDiffMuMin_minmax= 1;
	
				if(isBadDiffMu_minmax[k]){
					alphaDiffMuMin_minmax= alphaDiffMuList_minmax[k].Min();
					alphaDiffMuOpt= alphaDiffMuMin_minmax/fConstraintAlphaScale;	
					cout<<"INFO: Component No. "<<k+1<<":  alphaDiffMuMin="<<alphaDiffMuMin_minmax<<"  alphaDiffMuOpt="<<alphaDiffMuOpt<<endl;
				}

				
				if(alphaDiffMuOpt<fConstraintAlphaTolerance){
					hasToRegenerateMu= true;
					hasToBeGenerated[k]= true;	
					cout<<"INFO: Mu diff for component "<<k+1<<" is stuck in bound constraint ("<<alphaDiffMuOpt<<"<"<<fConstraintAlphaTolerance<<")"<<endl;	
				}

				// Update mu values				
				for(int j=0;j<fNDim;j++){
					double mu_k= fMu[k](j,0);
					double safemu_k= fMu_safe[k](j,0);	
					double diffmu_k= fMuDiff[k](j,0);
					double safediffmu_k= fMuDiff_safe[k](j,0);		

					if(fConstraintFirstMeanPar) fMu[0](j,0)= (1.-alphaDiffMuOpt)*fMu_safe[0](j,0) + alphaDiffMuOpt*fMu[0](j,0);	
					if(!fFixMeanDiffPar) fMuDiff[k](j,0)= (1.-alphaDiffMuOpt)*safediffmu_k + alphaDiffMuOpt*diffmu_k;	
					fMu[k](j,0)= fMu[0](j,0)+fMuDiff[k](j,0);
				}//end loop dim

				//Update ksi value
				fKsi[k]= MathUtilities::MuToKsi(fMu[k],fDelta[k],fNu[k]);
				fKsiDiff[k]= fKsi[k]-fKsi[0];

			}//end loop components


			if(hasToRegenerateMu && fUseRandomRegenerationAfterStuck){
				cout<<"INFO: Regenerating mu pars..."<<endl;
				RandomizeMeanPar(hasToBeGenerated);	
			}

			//## Check exit condition
			if(!isDiffMuBoundConstraintViolated) break;

		}//end loop while

	}//close if !fFixMeanPar && fUseMeanDiffConstraint		


	//## LOCATION DIFF PAR CONSTRAINT
	if(!fFixMeanPar && fUseLocationDiffConstraint){
		
		while(true){
			bool isBadDiffKsi= false;
			bool isBadDiffKsi_minmax[fNComponents];
			bool isDiffKsiBoundConstraintViolated= false;
			std::vector<TMatrixD> alphaDiffKsiList_minmax;
			alphaDiffKsiList_minmax.clear();
			alphaDiffKsiList_minmax.resize(0);

			for(int k=0;k<fNComponents;k++){
				isBadDiffKsi_minmax[k]= false;

				TMatrixD alphaDiffKsi_minmax(fNDim,1);
				alphaDiffKsi_minmax.Zero();
	
				for(int j=0;j<fNDim;j++){
					alphaDiffKsi_minmax(j,0)= 1;
				
					double constraintSign= 1;
					
					double a_k= fKsiDiff_min[k](j,0);
					double b_k= fKsiDiff_max[k](j,0);
					
					double diffksi= fKsiDiff[k](j,0);
					double safediffksi= fKsiDiff_safe[k](j,0);
				
					double alphaKsiDiff_minBound= (a_k-safediffksi)/(diffksi-safediffksi);
					double alphaKsiDiff_maxBound= (b_k-safediffksi)/(diffksi-safediffksi);

					if(constraintSign*diffksi<constraintSign*a_k){//min constraint violated
						alphaDiffKsi_minmax(j,0)= alphaKsiDiff_minBound;
						isBadDiffKsi_minmax[k]= true;
						isDiffKsiBoundConstraintViolated= true;
						cout<<"INFO: Location min diff constraint violated for component no. "<<k+1<<":  diffKsi="<<diffksi<<"  diffKsiMin="<<a_k<<"  safediffKsi(k)="<<safediffksi<<"  alpha="<<alphaKsiDiff_minBound<<endl;
					}
					
					if(constraintSign*diffksi>b_k){//max constraint violated
						alphaDiffKsi_minmax(j,0)= alphaKsiDiff_maxBound;
						isBadDiffKsi_minmax[k]= true;
						isDiffKsiBoundConstraintViolated= true;
						cout<<"INFO: Location max diff constraint violated for component no. "<<k+1<<":  diffKsi="<<diffksi<<"  diffKsiMax="<<b_k<<"  safediffKsi(k)="<<safediffksi<<"  alpha="<<alphaKsiDiff_maxBound<<endl;
					}	
					
				}//end loop dim

				alphaDiffKsiList_minmax.push_back(alphaDiffKsi_minmax);

			}//end loop components

			cout<<"INFO: Updating location parameters..."<<endl;
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateKsi= false;

			for(int k=1;k<fNComponents;k++){
			
				double alphaDiffKsiOpt= 1;
				double alphaDiffKsiMin_minmax= 1;
	
				if(isBadDiffKsi_minmax[k]){
					alphaDiffKsiMin_minmax= alphaDiffKsiList_minmax[k].Min();
					alphaDiffKsiOpt= alphaDiffKsiMin_minmax/fConstraintAlphaScale;
				}

				cout<<"INFO: Component no. "<<k+1<<": alphaDiffKsiMin="<<alphaDiffKsiMin_minmax<<"  alphaDiffKsiOpt="<<alphaDiffKsiOpt<<endl;

				if(alphaDiffKsiOpt<fConstraintAlphaTolerance){
					hasToRegenerateKsi= true;
					hasToBeGenerated[k]= true;	
					cout<<"INFO: Ksi diff for component "<<k+1<<" is stuck in bound constraint ("<<alphaDiffKsiOpt<<"<"<<fConstraintAlphaTolerance<<")"<<endl;	
				}

				// Update ksi values			
				for(int j=0;j<fNDim;j++){
					double ksi_k= fKsi[k](j,0);
					double safeksi_k= fKsi_safe[k](j,0);	
					double diffksi_k= fKsiDiff[k](j,0);
					double safediffksi_k= fKsiDiff_safe[k](j,0);		

					if(fConstraintFirstMeanPar) fKsi[0](j,0)= (1.-alphaDiffKsiOpt)*fKsi_safe[0](j,0) + alphaDiffKsiOpt*fKsi[0](j,0);	
					if(!fFixMeanDiffPar) fKsiDiff[k](j,0)= (1.-alphaDiffKsiOpt)*safediffksi_k + alphaDiffKsiOpt*diffksi_k;	
					fKsi[k](j,0)= fKsi[0](j,0)+fKsiDiff[k](j,0);
				}//end loop dim
	
				fMu[k]= MathUtilities::KsiToMu(fKsi[k],fDelta[k],fNu[k]);
				fMuDiff[k]= fMu[0]-fMu[k];
	
			}//end loop components


			if(hasToRegenerateKsi && fUseRandomRegenerationAfterStuck){
				RandomizeMeanPar(hasToBeGenerated);	
			}

			if(!isDiffKsiBoundConstraintViolated) break;

		}//end while loop

	}//close if !fFixMeanPar && fUseLocationDiffConstraint	

	
	//##################################
	//##   SIGMA CONSTRAINTS          ##
	//##################################
	cout<<endl;
	cout<<"==> COVARIANCE CONSTRAINT"<<endl;

	//## Update covariance matrix after delta and nu update
	cout<<"INFO: Recompute covariance matrix after delta and nu update..."<<endl;
	for(int k=0;k<fNComponents;k++) {
		fV[k]= MathUtilities::ScaleMatrixToCovarianceMatrix(fSigma[k], fDelta[k], fNu[k]);	
		fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);

		cout<<"== Component "<<k+1<<" =="<<endl;
		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<"SigmaEigen= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigmaEigen[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<"V= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fV[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<"VEigen= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fVEigen[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
	}//end loop components


	//## Apply constraints?
	if(!fFixCovariancePar){

		//## Calculate Sigma & V eigenvalues-eigenvectors
		cout<<"INFO: Compute covariance matrix eigenvalues and eigenvectors ..."<<endl;
		std::vector<TMatrixD> fSigmaEigenvectors;
		fSigmaEigenvectors.clear();
		fSigmaEigenvectors.resize(0);

		std::vector<TMatrixD> fSigmaEigenvectorsInv;
		fSigmaEigenvectorsInv.clear();
		fSigmaEigenvectorsInv.resize(0);		
	
		
		std::vector<TMatrixD> fVEigenvectors;
		fVEigenvectors.clear();
		fVEigenvectors.resize(0);

		std::vector<TMatrixD> fVEigenvectorsInv;
		fVEigenvectorsInv.clear();
		fVEigenvectorsInv.resize(0);

		for(int k=0;k<fNComponents;k++){
			//scale matrix eigen
			TMatrixDEigen SigmaDecomposition(fSigma[k]);
			TMatrixD SigmaEigenvalues(fNDim,fNDim);
			SigmaEigenvalues= SigmaDecomposition.GetEigenValues();
			fSigmaEigen[k]= SigmaEigenvalues;

			TMatrixD SigmaEigenvectors(fNDim,fNDim);
			SigmaEigenvectors= SigmaDecomposition.GetEigenVectors();
			fSigmaEigenvectors.push_back(SigmaEigenvectors);
			TMatrixD SigmaEigenvectorsInv(fNDim,fNDim);	
			SigmaEigenvectorsInv= SigmaEigenvectors;
			double det= 0;
			SigmaEigenvectorsInv= SigmaEigenvectorsInv.Invert(&det);	
			if (det<=0) {
				cerr<<"WARNING: Scale matrix eigenvectors inversion failed (det="<<det<<")!"<<endl;
			}
			fSigmaEigenvectorsInv.push_back(SigmaEigenvectorsInv);	
		
			//covariance eigen
			TMatrixDEigen VDecomposition(fV[k]);
			TMatrixD VEigenvalues(fNDim,fNDim);
			VEigenvalues= VDecomposition.GetEigenValues();
			fVEigen[k]= SigmaEigenvalues;

			TMatrixD VEigenvectors(fNDim,fNDim);
			VEigenvectors= VDecomposition.GetEigenVectors();
			fVEigenvectors.push_back(VEigenvectors);

			TMatrixD VEigenvectorsInv(fNDim,fNDim);	
			VEigenvectorsInv= VEigenvectors;
			det= 0;			
			VEigenvectorsInv= VEigenvectorsInv.Invert(&det);
			if (det<=0) {
				cerr<<"WARNING: Covariance matrix eigenvectors inversion failed (V("<<fV[k](0,0)<<","<<fV[k](0,1)<<","<<fV[k](1,1)<<") VEigen("<<fVEigen[k](0,0)<<","<<fVEigen[k](0,1)<<","<<fVEigen[k](1,1)<<")  VEigenvect("<<VEigenvectors(0,0)<<","<<VEigenvectors(0,1)<<","<<VEigenvectors(1,0)<<","<<VEigenvectors(1,1)<<")  det="<<det<<")!"<<endl;
			}

			fVEigenvectorsInv.push_back(VEigenvectorsInv);
		}//end loop components





		//## COVARIANCE CONSTRAINT
		//group constraint
		std::vector<double> alphaVList;
		alphaVList.clear();
		alphaVList.resize(0);

		bool isBadV= false;
		bool isVGroupConstraintViolated= false;

		if(fUseCovarianceConstraint){

			for(int j=0;j<fNDim;j++){
				double constraintSign= 1;
			
				for(int k=0;k<fNComponents-1;k++){
					double V_k= fV[k](j,j);
					double V_k_1= fV[k+1](j,j);
					double safeV_k= fV_safe[k](j,j);
					double safeV_k_1= fV_safe[k+1](j,j);
					double denom= 1.-(V_k_1-V_k)/(safeV_k_1-safeV_k);
					double alphaV= 1./denom;
		
					if(constraintSign*V_k>constraintSign*V_k_1) continue;//constraint satisfied...skip
					isBadV= true;
					isVGroupConstraintViolated= true;
					alphaVList.push_back(alphaV);

					cout<<"INFO: Covariance constraint violated for component no. "<<k+1<<":  V(k)="<<V_k<<"  V(k+1)="<<V_k_1<<"  safeV(k)="<<safeV_k<<"  safeV(k+1)="<<safeV_k_1<<"  alpha="<<alphaV<<endl;
		
				}//end loop components	
			}//end loop dim

		}//close if fUseCovarianceConstraint
	

		//min-max constraint
		bool isBadV_minmax[fNComponents];
		bool isVBoundConstraintViolated= false;
		std::vector<TMatrixD> alphaVList_minmax;
		alphaVList_minmax.clear();
		alphaVList_minmax.resize(0);

		if(fUseCovarianceBoundConstraint){
			for(int k=0;k<fNComponents;k++){
				isBadV_minmax[k]= false;

				TMatrixD alphaV_minmax(fNDim,fNDim);
				alphaV_minmax.Zero();

				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double constraintSign= 1;
						alphaV_minmax[j][l]= 1;

						//## Do not apply constraint on covariance terms
						if(j!=l) continue;

						double V= fV[k](j,l);
						double safeV= fV_safe[k](j,l);
						double Vmin= fSigma_min[k](j,l);
						double Vmax= fSigma_max[k](j,l);
						double alphaV_minBound= (Vmin-safeV)/(V-safeV);
						double alphaV_maxBound= (Vmax-safeV)/(V-safeV);
		
						if(constraintSign*V<constraintSign*Vmin){//min constraint violated
							alphaV_minmax[j][l]= alphaV_minBound;
							isBadV_minmax[k]= true;
							isVBoundConstraintViolated= true;
							if(j==l) cout<<"INFO: Variance min constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmin(k)="<<Vmin<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_minBound<<endl;	
							else cout<<"INFO: Covariance min constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmin(k)="<<Vmin<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_minBound<<endl;
						}
						if(constraintSign*V>constraintSign*Vmax){//max constraint violated
							alphaV_minmax[j][l]= alphaV_maxBound;
							isBadV_minmax[k]= true;
							isVBoundConstraintViolated= true;
							if(j==l) cout<<"INFO: Variance max constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmax(k)="<<Vmax<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_maxBound<<endl;
							else cout<<"INFO: Covariance max constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmax(k)="<<Vmax<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_maxBound<<endl;
						}
					}//end loop dim
				}//end loop dim

				alphaVList_minmax.push_back(alphaV_minmax);
			}//end loop components

		}//close if fUseCovarianceBoundConstraint

	
		//Update Sigma values
		if( fUseCovarianceConstraint || fUseCovarianceBoundConstraint ){
			cout<<"INFO: Updating covariance matrix..."<<endl;

			double alphaVMin= 1;
			cout<<"INFO: alphaVList(";
			if(fUseCovarianceConstraint && isVGroupConstraintViolated){
				for(unsigned int i=0;i<alphaVList.size();i++){
					cout<<alphaVList[i]<<",";
					if(alphaVList[i]<alphaVMin) alphaVMin= alphaVList[i];
				}			
				cout<<")"<<endl;
			}	

			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateSigma= false;

			for(int k=0;k<fNComponents;k++){
				if(!isBadV && !isBadV_minmax[k]) continue;
		
				double alphaVMin_minmax= alphaVList_minmax[k].Min();
				double alphaVOpt= min(alphaVMin_minmax,alphaVMin)/fConstraintAlphaScale;
				cout<<"INFO: Component no. "<<k+1<<":  alphaVMin="<<alphaVMin<<"  alphaVMin_minmax="<<alphaVMin_minmax<<"  alphaVOpt="<<alphaVOpt<<endl;

	
				if(alphaVOpt<fConstraintAlphaTolerance) {
					cout<<"INFO: Sigma for component "<<k+1<<" is stuck in constraint ("<<alphaVOpt<<"<"<<fConstraintAlphaTolerance<<")...regenerate!"<<endl;
					hasToRegenerateSigma= true;
					hasToBeGenerated[k]= true;
				}

				//## Update sigma values
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double V= fV[k](j,l);
						double safeV= fV_safe[k](j,l);
						fV[k](j,l)= (1.-alphaVOpt)*safeV + alphaVOpt*V;
					}//end loop dim
				}//end loop dim
					
			}//end loop components

			if(hasToRegenerateSigma && fUseRandomRegenerationAfterStuck){
				RandomizeSigmaPar(hasToBeGenerated);	
			}
		
			//## Re-Calculate covariance and sigma starting from the updated eigenvalues
			for(int k=0;k<fNComponents;k++){
				fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);
				fSigma[k]= MathUtilities::CovarianceMatrixToScaleMatrix(fV[k],fDelta[k],fNu[k]);
				fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);				
			}//end loop components

		}//close if (fUseCovarianceConstraint || fUseCovarianceBoundConstraint)
		




		//## SCALE MATRIX CONSTRAINT
		//group constraint
		std::vector<double> alphaSigmaList;
		alphaSigmaList.clear();
		alphaSigmaList.resize(0);

		bool isBadSigma= false;
		bool isSigmaGroupConstraintViolated= false;

		if(fUseScaleMatrixConstraint){

			for(int j=0;j<fNDim;j++){
				double constraintSign= 1;
			
				for(int k=0;k<fNComponents-1;k++){
					double V_k= fSigma[k](j,j);
					double V_k_1= fSigma[k+1](j,j);
					double safeV_k= fSigma_safe[k](j,j);
					double safeV_k_1= fSigma_safe[k+1](j,j);
					double denom= 1.-(V_k_1-V_k)/(safeV_k_1-safeV_k);
					double alphaV= 1./denom;
		
					if(constraintSign*V_k>constraintSign*V_k_1) continue;//constraint satisfied...skip
					isBadSigma= true;
					isSigmaGroupConstraintViolated= true;
					alphaSigmaList.push_back(alphaV);

					cout<<"INFO: Covariance constraint violated for component no. "<<k+1<<":  V(k)="<<V_k<<"  V(k+1)="<<V_k_1<<"  safeV(k)="<<safeV_k<<"  safeV(k+1)="<<safeV_k_1<<"  alpha="<<alphaV<<endl;
		
				}//end loop components	
			}//end loop dim

		}//close if fUseCovarianceConstraint
	

		//min-max constraint
		bool isBadSigma_minmax[fNComponents];
		bool isSigmaBoundConstraintViolated= false;
		std::vector<TMatrixD> alphaSigmaList_minmax;
		alphaSigmaList_minmax.clear();
		alphaSigmaList_minmax.resize(0);

		if(fUseScaleMatrixBoundConstraint){
			for(int k=0;k<fNComponents;k++){
				isBadSigma_minmax[k]= false;

				TMatrixD alphaV_minmax(fNDim,fNDim);
				alphaV_minmax.Zero();

				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double constraintSign= 1;
						alphaV_minmax[j][l]= 1;

						//## Do not apply constraint on covariance terms
						if(j!=l) continue;

						double V= fSigma[k](j,l);
						double safeV= fSigma_safe[k](j,l);
						double Vmin= fSigma_min[k](j,l);
						double Vmax= fSigma_max[k](j,l);
						double alphaV_minBound= (Vmin-safeV)/(V-safeV);
						double alphaV_maxBound= (Vmax-safeV)/(V-safeV);
		
						if(constraintSign*V<constraintSign*Vmin){//min constraint violated
							alphaV_minmax[j][l]= alphaV_minBound;
							isBadSigma_minmax[k]= true;
							isSigmaBoundConstraintViolated= true;
							if(j==l) cout<<"INFO: Variance min constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmin(k)="<<Vmin<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_minBound<<endl;	
							else cout<<"INFO: Covariance min constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmin(k)="<<Vmin<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_minBound<<endl;
						}
						if(constraintSign*V>constraintSign*Vmax){//max constraint violated
							alphaV_minmax[j][l]= alphaV_maxBound;
							isBadSigma_minmax[k]= true;
							isSigmaBoundConstraintViolated= true;
							if(j==l) cout<<"INFO: Variance max constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmax(k)="<<Vmax<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_maxBound<<endl;
							else cout<<"INFO: Covariance max constraint violated for component no. "<<k+1<<":  V(k)="<<V<<"  Vmax(k)="<<Vmax<<"  safeV(k)="<<safeV<<"  alpha="<<alphaV_maxBound<<endl;
						}
					}//end loop dim
				}//end loop dim

				alphaSigmaList_minmax.push_back(alphaV_minmax);
			}//end loop components

		}//close if fUseScaleMatrixBoundConstraint

	
		//Update Sigma values
		if( fUseScaleMatrixConstraint || fUseScaleMatrixBoundConstraint ){
			cout<<"INFO: Updating covariance matrix..."<<endl;

			double alphaVMin= 1;
			cout<<"INFO: alphaSigmaList(";
			if(fUseCovarianceConstraint && isSigmaGroupConstraintViolated){
				for(unsigned int i=0;i<alphaSigmaList.size();i++){
					cout<<alphaSigmaList[i]<<",";
					if(alphaSigmaList[i]<alphaVMin) alphaVMin= alphaSigmaList[i];
				}			
				cout<<")"<<endl;
			}	

			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateSigma= false;

			for(int k=0;k<fNComponents;k++){
				if(!isBadSigma && !isBadSigma_minmax[k]) continue;
		
				double alphaVMin_minmax= alphaSigmaList_minmax[k].Min();
				double alphaVOpt= min(alphaVMin_minmax,alphaVMin)/fConstraintAlphaScale;
				cout<<"INFO: Component no. "<<k+1<<":  alphaVMin="<<alphaVMin<<"  alphaVMin_minmax="<<alphaVMin_minmax<<"  alphaVOpt="<<alphaVOpt<<endl;

	
				if(alphaVOpt<fConstraintAlphaTolerance) {
					cout<<"INFO: Sigma for component "<<k+1<<" is stuck in constraint ("<<alphaVOpt<<"<"<<fConstraintAlphaTolerance<<")...regenerate!"<<endl;
					hasToRegenerateSigma= true;
					hasToBeGenerated[k]= true;
				}

				//## Update sigma values
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double V= fSigma[k](j,l);
						double safeV= fSigma_safe[k](j,l);
						fSigma[k](j,l)= (1.-alphaVOpt)*safeV + alphaVOpt*V;
					}//end loop dim
				}//end loop dim
					
			}//end loop components

			if(hasToRegenerateSigma && fUseRandomRegenerationAfterStuck){
				RandomizeSigmaPar(hasToBeGenerated);	
			}
		
			//## Re-Calculate covariance and sigma starting from the updated eigenvalues
			for(int k=0;k<fNComponents;k++){
				fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);
				fV[k]= MathUtilities::ScaleMatrixToCovarianceMatrix(fSigma[k],fDelta[k],fNu[k]);
				fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);				
			}//end loop components

		}//close if (fUseScaleMatrixConstraint || fUseScaleMatrixBoundConstraint)
		

		//## COVARIANCE EIGEN CONSTRAINT
		std::vector<double> alphaVEigenList;
		alphaVEigenList.clear();
		alphaVEigenList.resize(0);
		bool isBadVEigen= false;
		bool isVEigenGroupConstraintViolated= false;

		if(fUseCovarianceEigenConstraint){
	
			for(int j=0;j<fNDim;j++){
				double constraintSign= 1;
			
				for(int k=0;k<fNComponents-1;k++){
					double Lambda_k= fVEigen[k](j,j);
					double Lambda_k_1= fVEigen[k+1](j,j);
					double safeLambda_k= fVEigen_safe[k](j,j);
					double safeLambda_k_1= fVEigen_safe[k+1](j,j);
					double denom= 1.-(Lambda_k_1-Lambda_k)/(safeLambda_k_1-safeLambda_k);
					double alphaLambda= 1./denom;
		
					if(constraintSign*Lambda_k>constraintSign*Lambda_k_1) continue;//constraint satisfied...skip
					isBadVEigen= true;
					isVEigenGroupConstraintViolated= true;

					alphaVEigenList.push_back(alphaLambda);

					cout<<"INFO: Covariance eigen constraint violated for component no. "<<k+1<<": l(k)="<<Lambda_k<<"  l(k+1)="<<Lambda_k_1<<"  safel(k)="<<safeLambda_k<<"  safel(k+1)="<<safeLambda_k_1<<"  alpha="<<alphaLambda<<endl;
		
				}//end loop components	
			}//end loop dim

		}//close if fUseCovarianceEigenConstraint
	

		//## COVARIANCE EIGEN BOUND CONSTRAINT
		bool isBadVEigen_minmax[fNComponents];
		bool isVEigenBoundConstraintViolated= false;
		std::vector<TMatrixD> alphaVEigenList_minmax;
		alphaVEigenList_minmax.clear();
		alphaVEigenList_minmax.resize(0);

		if(fUseCovarianceEigenBoundConstraint){

			cout<<"INFO: Check for covariance bound constraint..."<<endl;

			for(int k=0;k<fNComponents;k++){
				isBadVEigen_minmax[k]= false;

				TMatrixD alphaLambda_minmax(fNDim,fNDim);
				alphaLambda_minmax.Zero();

				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double constraintSign= 1;
						alphaLambda_minmax(j,l)= 1;

						//## Do not apply constraint on covariance terms
						if(j!=l) continue;

						double Lambda= fVEigen[k](j,l);
						double safeLambda= fVEigen_safe[k](j,l);
						double Lambdamin= fSigmaEigen_min[k](j,0);
						double Lambdamax= fSigmaEigen_max[k](j,0);
						double alphaLambda_minBound= (Lambdamin-safeLambda)/(Lambda-safeLambda);
						double alphaLambda_maxBound= (Lambdamax-safeLambda)/(Lambda-safeLambda);
		
						if(constraintSign*Lambda<constraintSign*Lambdamin){//min constraint violated
							alphaLambda_minmax(j,l)= alphaLambda_minBound;
							isBadVEigen_minmax[k]= true;
							isVEigenBoundConstraintViolated= true;
							cout<<"INFO: Covariance min eigen bound constraint violated for component no. "<<k+1<<":  l(k)="<<Lambda<<"  Lambdamin(k)="<<Lambdamin<<"  safel(k)="<<safeLambda<<"  alpha="<<alphaLambda_minBound<<endl;	
						}
						if(constraintSign*Lambda>constraintSign*Lambdamax){//max constraint violated
							alphaLambda_minmax(j,l)= alphaLambda_maxBound;
							isBadVEigen_minmax[k]= true;
							isVEigenBoundConstraintViolated= true;
							cout<<"INFO: Covariance max eigen bound constraint violated for component no. "<<k+1<<":  l(k)="<<Lambda<<"  Lambdamax(k)="<<Lambdamax<<"  safel(k)="<<safeLambda<<"  alpha="<<alphaLambda_maxBound<<endl;		
						}
					}//end loop dim
				}//end loop dim
				alphaVEigenList_minmax.push_back(alphaLambda_minmax);
			
			}//end loop components
		}//close if fUseCovarianceBoundConstraint

	
		//## Update covariance eigen values		
		if(fUseCovarianceEigenConstraint || fUseCovarianceEigenBoundConstraint){
			cout<<"INFO: Updating covariance matrix..."<<endl;

			double alphaEigenMin= 1;
			if(fUseCovarianceEigenConstraint && isVEigenGroupConstraintViolated){
				cout<<"INFO: alphaVEigenList (";
				for(unsigned int i=0;i<alphaVEigenList.size();i++){
					cout<<alphaVEigenList[i]<<",";
					if(alphaVEigenList[i]<alphaEigenMin) alphaEigenMin= alphaVEigenList[i];
				}				
				cout<<")"<<endl;
			}	
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateSigma= false;

			for(int k=0;k<fNComponents;k++){
				if(!isBadVEigen && !isBadVEigen_minmax[k]) continue;
				double alphaEigenMin_minmax= alphaVEigenList_minmax[k].Min();
				double alphaEigenOpt= min(alphaEigenMin_minmax,alphaEigenMin)/fConstraintAlphaScale;
				cout<<"INFO: Component no. "<<k+1<<": alphaEigenMin="<<alphaEigenMin<<"  alphaEigenMin_minmax="<<alphaEigenMin_minmax<<"  alphaEigenOpt="<<alphaEigenOpt<<endl;

				if(alphaEigenOpt<fConstraintAlphaTolerance) {
					cout<<"INFO: Covariance eigen for component "<<k+1<<" is stuck in constraint ("<<alphaEigenOpt<<"<"<<fConstraintAlphaTolerance<<")...regenerate!"<<endl;
					hasToRegenerateSigma= true;
					hasToBeGenerated[k]= true;
				}
			
				//## Update sigma eigen values
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double Lambda= fVEigen[k](j,l);
						double safeLambda= fVEigen_safe[k](j,l);
						fVEigen[k](j,l)= (1.-alphaEigenOpt)*safeLambda + alphaEigenOpt*Lambda;
					}//end loop dim
				}//end loop dim
					
			}//end loop components

			if(hasToRegenerateSigma && fUseRandomRegenerationAfterStuck){
				cout<<"INFO: Randomize sigma pars..."<<endl;
				RandomizeSigmaPar(hasToBeGenerated);	
			}

			//## Re-Calculate covariance and sigma starting from the updated eigenvalues
			for(int k=0;k<fNComponents;k++){
				fV[k]= fVEigenvectors[k]*fVEigen[k]*fVEigenvectorsInv[k];
				fSigma[k]= MathUtilities::CovarianceMatrixToScaleMatrix(fV[k],fDelta[k],fNu[k]);
				fSigmaEigen[k]= MathUtilities::GetEigenDecomposition(fSigma[k]);
			}//end loop components
			
		}//close if (fUseCovarianceEigenConstraint || fUseCovarianceEigenBoundConstraint)

		

		//## SCALE MATRIX EIGEN CONSTRAINT
		std::vector<double> alphaEigenList;
		alphaEigenList.clear();
		alphaEigenList.resize(0);
		bool isBadSigmaEigen= false;
		bool isSigmaEigenGroupConstraintViolated= false;
	
		if(fUseScaleMatrixEigenConstraint){

			for(int j=0;j<fNDim;j++){
				double constraintSign= 1;
			
				for(int k=0;k<fNComponents-1;k++){
					double Lambda_k= fSigmaEigen[k](j,j);
					double Lambda_k_1= fSigmaEigen[k+1](j,j);
					double safeLambda_k= fSigmaEigen_safe[k](j,j);
					double safeLambda_k_1= fSigmaEigen_safe[k+1](j,j);
					double denom= 1.-(Lambda_k_1-Lambda_k)/(safeLambda_k_1-safeLambda_k);
					double alphaLambda= 1./denom;
		
					if(constraintSign*Lambda_k>constraintSign*Lambda_k_1) continue;//constraint satisfied...skip
					isBadSigmaEigen= true;
					isSigmaEigenGroupConstraintViolated= true;

					alphaEigenList.push_back(alphaLambda);

					cout<<"INFO: Scale matrix eigen constraint violated for component no. "<<k+1<<":  l(k)="<<Lambda_k<<"  l(k+1)="<<Lambda_k_1<<"  safel(k)="<<safeLambda_k<<"  safel(k+1)="<<safeLambda_k_1<<"  alpha="<<alphaLambda<<endl;
		
				}//end loop components	
			}//end loop dim
		}//close if fUseScaleMatrixEigenConstraint
	

		//## SCALE MATRIX EIGEN BOUND CONSTRAINT
		bool isBadSigmaEigen_minmax[fNComponents];
		bool isSigmaEigenBoundConstraintViolated= false;
		std::vector<TMatrixD> alphaEigenList_minmax;
		alphaEigenList_minmax.clear();
		alphaEigenList_minmax.resize(0);

		if(fUseScaleMatrixEigenBoundConstraint){
			for(int k=0;k<fNComponents;k++){
				isBadSigmaEigen_minmax[k]= false;

				TMatrixD alphaLambda_minmax(fNDim,fNDim);
				alphaLambda_minmax.Zero();

				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double constraintSign= 1;
						alphaLambda_minmax(j,l)= 1;

						//## Do not apply constraint on covariance terms
						if(j!=l) continue;

						double Lambda= fSigmaEigen[k](j,l);
						double safeLambda= fSigmaEigen_safe[k](j,l);
						double Lambdamin= fSigmaEigen_min[k](j,0);
						double Lambdamax= fSigmaEigen_max[k](j,0);
						double alphaLambda_minBound= (Lambdamin-safeLambda)/(Lambda-safeLambda);
						double alphaLambda_maxBound= (Lambdamax-safeLambda)/(Lambda-safeLambda);
		
						if(constraintSign*Lambda<constraintSign*Lambdamin){//min constraint violated
							alphaLambda_minmax(j,l)= alphaLambda_minBound;
							isBadSigmaEigen_minmax[k]= true;
							isSigmaEigenBoundConstraintViolated= true;
							cout<<"INFO: Scale matrix eigen min bound constraint violated for component no. "<<k+1<<":  l(k)="<<Lambda<<"  Lambdamin(k)="<<Lambdamin<<"  safel(k)="<<safeLambda<<"  alpha="<<alphaLambda_minBound<<endl;	
						}

						if(constraintSign*Lambda>constraintSign*Lambdamax){//max constraint violated
							alphaLambda_minmax(j,l)= alphaLambda_maxBound;
							isBadSigmaEigen_minmax[k]= true;
							isSigmaEigenBoundConstraintViolated= true;
							cout<<"INFO: Scale matrix eigen max bound constraint violated for component no. "<<k+1<<":  l(k)="<<Lambda<<"  Lambdamax(k)="<<Lambdamax<<"  safel(k)="<<safeLambda<<"  alpha="<<alphaLambda_maxBound<<endl;	
						}
					}//end loop dim
				}//end loop dim

				alphaEigenList_minmax.push_back(alphaLambda_minmax);
			
			}//end loop components
		}//close if fUseScaleMatrixEigenBoundConstraint

	
		//## Update Sigma eigen values
		if(fUseScaleMatrixEigenConstraint || fUseScaleMatrixEigenBoundConstraint){
			cout<<"INFO: Updating scale matrix..."<<endl;
		
			double alphaEigenMin= 1;
			cout<<"INFO: alphaEigenList(";
			if(fUseScaleMatrixEigenConstraint && isSigmaEigenGroupConstraintViolated){
				for(unsigned int i=0;i<alphaEigenList.size();i++){
					cout<<alphaEigenList[i]<<",";
					if(alphaEigenList[i]<alphaEigenMin) alphaEigenMin= alphaEigenList[i];
				}	
				cout<<")"<<endl;		
			}	
			std::vector<bool> hasToBeGenerated;
			hasToBeGenerated.assign(fNComponents,false);
			bool hasToRegenerateSigma= false;

			for(int k=0;k<fNComponents;k++){
				if(!isBadSigmaEigen && !isBadSigmaEigen_minmax[k]) continue;
				double alphaEigenMin_minmax= alphaEigenList_minmax[k].Min();
				double alphaEigenOpt= min(alphaEigenMin_minmax,alphaEigenMin)/fConstraintAlphaScale;
				cout<<"INFO: Component no. "<<k+1<<": alphaEigenMin="<<alphaEigenMin<<"  alphaEigenMin_minmax="<<alphaEigenMin_minmax<<"  alphaEigenOpt="<<alphaEigenOpt<<endl;

				if(alphaEigenOpt<fConstraintAlphaTolerance) {
					cout<<"INFO: Scale matrix eigen for component "<<k+1<<" is stuck in constraint ("<<alphaEigenOpt<<"<"<<fConstraintAlphaTolerance<<")...regenerate!"<<endl;
					hasToRegenerateSigma= true;
					hasToBeGenerated[k]= true;
				}
			
				//## Update sigma eigen values
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						double Lambda= fSigmaEigen[k](j,l);
						double safeLambda= fSigmaEigen_safe[k](j,l);
						fSigmaEigen[k](j,l)= (1.-alphaEigenOpt)*safeLambda + alphaEigenOpt*Lambda;
					}//end loop dim
				}//end loop dim		
			}//end loop components

			if(hasToRegenerateSigma && fUseRandomRegenerationAfterStuck){
				RandomizeSigmaPar(hasToBeGenerated);	
			}

			//## Re-Calculate Sigma and V starting from the updated eigenvalues
			for(int k=0;k<fNComponents;k++){
				fSigma[k]= fSigmaEigenvectors[k]*fSigmaEigen[k]*fSigmaEigenvectorsInv[k];
				fV[k]= MathUtilities::ScaleMatrixToCovarianceMatrix(fSigma[k],fDelta[k],fNu[k]);
				fVEigen[k]= MathUtilities::GetEigenDecomposition(fV[k]);
			}//end loop components
		
		}//close if (fUseScaleMatrixEigenConstraint || fUseScaleMatrixEigenBoundConstraint)

	}//close if !fFixCovariancePar



	//## Set the current constrained set as "safe"
	cout<<"==> PAR UPDATE AFTER CONSTRAIN STEP"<<endl;
	for(int k=0;k<fNComponents;k++){
		fKsi_safe[k]= fKsi[k];	
		fKsiDiff_safe[k]= fKsiDiff[k];
		fMu_safe[k]= fMu[k];	
		fMuDiff_safe[k]= fMuDiff[k];
		fSigma_safe[k]= fSigma[k];
		fSigmaEigen_safe[k]= fSigmaEigen[k];
		fV_safe[k]= fV[k];
		fVEigen_safe[k]= fVEigen[k];
		fDelta_safe[k]= fDelta[k];
		fNu_safe[k]= fNu[k];
		fP_safe[k]= fP[k];
		
		cout<<"== Component no. "<<k+1<<" =="<<endl;
		cout<<"p= "<<fP[k]<<endl;
		cout<<"Nu= "<<fNu[k]<<endl;
		cout<<"Ksi= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fKsi[k])(j,0)<<",";
		cout<<(fKsi[k])(fNDim-1,0)<<")"<<endl; 
		cout<<"Mu= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(j,0)<<",";
		cout<<(fKsi[k])(fNDim-1,0)<<")"<<endl; 
		cout<<"KsiDiff= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fKsiDiff[k])(j,0)<<",";
		cout<<(fKsiDiff[k])(fNDim-1,0)<<")"<<endl; 
		cout<<"MuDiff= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fMuDiff[k])(j,0)<<",";
		cout<<(fMuDiff[k])(fNDim-1,0)<<")"<<endl; 
	
		cout<<"delta= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fDelta[k])(j,0)<<",";
		cout<<(fDelta[k])(fNDim-1,0)<<")"<<endl; 
			
		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;

		cout<<"SigmaEigen= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigmaEigen[k])(j,l)<<",";
		}	
		cout<<")"<<endl;

		cout<<"V= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fV[k])(j,l)<<",";
		}	
		cout<<")"<<endl;

		cout<<"VEigen= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fVEigen[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
		cout<<"================"<<endl;

	}//end loop components
	*/

	return 0;

}//close RunEM_ConstrainStep()





int MSTMixtureFitter::RunFitter(bool initializePars)
{
	
	//================================================
	//====      INITIALIZE FIT START PARS
	//================================================
	//## Initialize starting pars
	if(initializePars){
		INFO_LOG("Initializing fit starting pars ...");
		if(InitializeFitPars()<0){
			ERROR_LOG("Failed to initialize fit starting parameters!");
			return -1;
		}
	}
	else{
		INFO_LOG("Using present fit parameters (e.g. computed with EM algorithm)...");
	}

	//================================================
	//====         INIT FITTER
	//================================================
	//Compute number of fit parameters
	int nPar= 0;	
	nPar+= fNComponents-1;//fraction
	nPar+= fNComponents*fNDim;//mean pars
	nPar+= fNComponents*fNDim*(fNDim+1)/2.;//sigma pars
	nPar+= fNComponents*fNDim;//delta pars
	nPar+= fNComponents;//nu par

	//Initialize minuit
	INFO_LOG("Initializing MINUIT fitter...");
	int par_counter= 0;
	int ierflag;
  double arglist[2];  

	TMinuit theMinuit(nPar); 
  theMinuit.SetPrintLevel(0);
  theMinuit.SetMaxIterations(fNIterations);
  theMinuit.SetFCN(MathUtils::MST_LikelihoodFcn);
  theMinuit.mnexcm("SET NOW", arglist ,0,ierflag);
  
	//## Set fraction
	if(fNComponents>1){
		INFO_LOG("Setting fraction fitting pars...");
		double parStep= 0.1;
		double parMin= 0.;
		double parMax= 1.;
  	double parStart[fNComponents-1];
  
		//Constraint (fractions must sum unity)
		parStart[0]= fP[0];
		double ConstrFactor= 1.-parStart[0];
		for(int k=1;k<fNComponents-1;k++){					
			parStart[k]= fP[k]/ConstrFactor;//equal values
			ConstrFactor*= (1.-parStart[k]);
		}

		for(int k=0;k<fNComponents-1;k++){
			TString parName= Form("f%d",k+1);
			theMinuit.mnparm(par_counter,parName,parStart[k],parStep,parMin,parMax,ierflag);
			if(fFixFractionPars) theMinuit.FixParameter(par_counter);
  		par_counter++;
		}
	}//close if

	//## Set ksi
	INFO_LOG("Setting Ksi start pars...");
	for(int k=0;k<fNComponents;k++){
		TMatrixD ksi_min;
		TMatrixD ksi_max;
		TMatrixD ksi;
		MathUtils::MuToKsi(ksi,fMu[k],fDelta[k],fNu[k]);
		if(fUseMeanBoundConstraint){
			MathUtils::MuToKsi(ksi_min,fMu_min[k],fDelta[k],fNu[k]);
			MathUtils::MuToKsi(ksi_max,fMu_max[k],fDelta[k],fNu[k]);
		}
	
		for(int j=0;j<fNDim;j++){
			TString parName= Form("Ksi%d_%d",k+1,j+1);
			double parStart= ksi(0,j);
			double parStep= 0.1;
			double parMin= 0;
			double parMax= 0;
			if(fUseMeanBoundConstraint){
				parMin= std::min(ksi_min(0,j),ksi_max(0,j));
				parMax= std::max(ksi_min(0,j),ksi_max(0,j));
			}
			theMinuit.mnparm(par_counter,parName,parStart,parStep,parMin,parMax,ierflag);
			if(fFixMeanPars) theMinuit.FixParameter(par_counter);
			par_counter++;
		}//end loop dims
	}//end loop components
	

	
	//## Set Sigma
	INFO_LOG("Setting Sigma start pars...");
	for(int k=0;k<fNComponents;k++){
		TMatrixD Omega;
		TMatrixD Omega_min;
		TMatrixD Omega_max;
		MathUtils::CovarianceToScaleMatrix(Omega,fSigma[k],fDelta[k],fNu[k]);
		if(fUseCovarianceBoundConstraint){
			MathUtils::CovarianceToScaleMatrix(Omega_min,fSigma_min[k],fDelta[k],fNu[k]);
			MathUtils::CovarianceToScaleMatrix(Omega_max,fSigma_max[k],fDelta[k],fNu[k]);		
		}

		for(int j=0;j<fNDim;j++){
			for(int l=j;l<fNDim;l++){
				double parStart= fOmega[k](j,l);
				double parStep= 0.1;
				double parMin= 0;
				double parMax= 0;
				if(fUseCovarianceBoundConstraint){
					parMin= std::min(Omega_min(j,l),Omega_max(j,l));
					parMin= std::max(Omega_min(j,l),Omega_max(j,l));
				}
				TString parName= Form("omega%d_%d%d",k+1,j+1,l+1);
				theMinuit.mnparm(par_counter,parName,parStart,parStep,parMin,parMax,ierflag);
				if(fFixCovariancePars) theMinuit.FixParameter(par_counter);
				par_counter++;	
			}
		}
	}//end loop components
  


	//## Set delta
	INFO_LOG("Setting delta start pars...");
	for(int k=0;k<fNComponents;k++){
		for(int j=0;j<fNDim;j++){
			double parStart= fDelta[k](0,j);
			double parStep= 0.1;
			double parMin= 0;
			double parMax= 0;
			if(fUseDeltaBoundConstraint){
				parMin= std::min(fDelta_min[k](0,j),fDelta_max[k](0,j));
				parMax= std::max(fDelta_min[k](0,j),fDelta_max[k](0,j));
			}
			TString parName= Form("delta%d_%d",k+1,j+1);
			theMinuit.mnparm(par_counter,parName,parStart,parStep,parMin,parMax,ierflag);
			if(fFixDeltaPars) theMinuit.FixParameter(par_counter);
			par_counter++;
		}
	}//end loop components

	//## Set nu
	INFO_LOG("Setting nu start pars ...");
	for(int k=0;k<fNComponents;k++){
		double parStart= fNu[k];
		double parStep= 1.;
		double parMin= 0;
		double parMax= 0;
		if(fUseNuBoundConstraint){
			parMin= std::min(fNu_min[k],fNu_max[k]);
			parMax= std::max(fNu_min[k],fNu_max[k]);
		}
		TString parName= Form("nu%d",k+1);
		theMinuit.mnparm(par_counter,parName,parStart,parStep,parMin,parMax,ierflag);
		if(fFixNuPars) theMinuit.FixParameter(par_counter);
  	par_counter++;
	}//end loop components

	//================================================
	//====         FIT DATA
	//================================================
	INFO_LOG("Fitting data...");

	//SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
  arglist[0]= 1;
	theMinuit.mnexcm("SET STR",arglist,1,ierflag);

	//##SET MINUIT ERROR
  arglist[0]= 0.5;//0.5 likelihood, 1 ChiSquare
  theMinuit.mnexcm("SET ERR",arglist,1,ierflag);

	//## SET MINUIT COMMAND
	arglist[0]= 10000;
	arglist[1]= 0.1;
	
  //theMinuit.mnexcm("CALL FCN", arglist ,0,ierflag);//just call FCN
  theMinuit.mnexcm("MINIMIZE", arglist, 2, ierflag);//use MINIMIZE for scan
	//theMinuit.mnexcm("MIGRAD", arglist, 2, ierflag);//use MIGRAD for scan
	//theMinuit.mnexcm("MINOS", arglist, 2, ierflag);//use MINIMIZE for scan

	//================================================
	//====        RETRIEVE FIT RESULTS
	//================================================	
	int fitStatus= theMinuit.GetStatus();
	fitStatus= ierflag;

  double LikelihoodMin, edm, errdef;
  int nvpar, nparx, icstat;
  theMinuit.mnstat(LikelihoodMin, edm, errdef, nvpar, nparx, icstat);
	fLogLikelihood= LikelihoodMin;

	//## Get fitted pars covariance matrix
	double covMatrix[nPar][nPar];		
	theMinuit.mnemat(&covMatrix[0][0],nPar);

	//TMatrixDSym mat(nPar); 
	//theMinuit.mnemat( mat.GetMatrixArray(), nPar);
  	
  cout<<"== COV MATRIX =="<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
	//TMatrixD FractionCovarianceMatrix(nPar,nPar);
			
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i

	/*
	for(int i=0;i<fNComponents-1;i++){
		for(int j=0;j<fNComponents-1;j++){
		 	FractionCovarianceMatrix(i,j)= covMatrix[i][j];	
		  //if(j==nPar-1) cout<<FractionCovarianceMatrix[i][j]<<endl;
		  //else cout<<FractionCovarianceMatrix[i][j]<<"  ";
		}//end loop j
	}//end loop i
	*/

	//## Get fitted parameters  
	double p,ep;
	par_counter= 0;
	std::vector<double> fittedPars;
	std::vector<double> fittedParErrors;
	for(int i=0;i<nPar;i++){
		theMinuit.GetParameter(i,p,ep);	
		fittedPars.push_back(p);
		fittedParErrors.push_back(ep);
	}

	//- Fraction
	if(fNComponents>1){	
		fP[0]= fittedPars[0];
		double ConstrFactor= (1.-fittedPars[0]);
		for(int k=1;k<fNComponents;k++){
			fP[k]= fittedPars[k]*ConstrFactor;
			ConstrFactor*= (1.-fittedPars[k]);
		}
		fP[fNComponents-1]= ConstrFactor;
		par_counter+= fNComponents-1;
	}//close if
	else{
		fP[0]= 1;
	}

	//- Ksi pars
	for(int k=0;k<fNComponents;k++) {
		for(int j=0;j<fNDim;j++){
			fKsi[k](0,j)= fittedPars[par_counter];
			par_counter++;
		}
	}

	//- Omega pars
	for(int k=0;k<fNComponents;k++) {
		for(int j=0;j<fNDim;j++){
			for(int l=j;l<fNDim;l++){
				double p= fittedPars[par_counter];
				fOmega[k](j,l)= p;
				fOmega[k](l,j)= p;
				par_counter++;
			}
		}
	}

	//- Delta pars
	for(int k=0;k<fNComponents;k++) {
		for(int j=0;j<fNDim;j++){
			fDelta[k](0,j)= fittedPars[par_counter];
			par_counter++;
		}
	}
	
	//- Nu pars
	for(int k=0;k<fNComponents;k++) {
		fNu[k]= fittedPars[par_counter]; 
		par_counter++;
	}

	//Compute mean & covariance
	for(int k=0;k<fNComponents;k++){
		MathUtils::KsiToMu(fMu[k],fKsi[k],fDelta[k],fNu[k]);
		MathUtils::ScaleToCovarianceMatrix(fSigma[k],fOmega[k],fDelta[k],fNu[k]);
	}

	//Print parameters
	INFO_LOG("Print fitted parameters");
	PrintPars();

	return 0;	

}//close RunFitter()


/*
std::vector<double> MSTMixtureFitter::GetParamUncertainty(std::vector<double> parameters,TMatrixD CovarianceMatrix,int Npar){

	unsigned int Npar_fit= parameters.size();
	for(int j=0;j<Npar_fit;j++) cout<<parameters[j]<<endl;
		
	int Ndim= CovarianceMatrix.GetNrows();

	//check dimensions
	if(Npar_fit!=Ndim){
		cerr<<"Fitter::GetParamUncertainty(): Npar_fit!=Ndim"<<endl;
		exit(1);
	}//close if

	
	std::vector<double> parErrors;
	parErrors.resize(Npar);
	for(unsigned int i=0;i<parErrors.size();i++) parErrors[i]=0.;

	//calculate derivative matrix
	TMatrixD DerMatrix(Npar,Npar_fit);
	for(unsigned int i=0;i<Npar;i++){
		for(unsigned int j=0;j<Npar_fit;j++){
			DerMatrix(i,j)= GetDerivativeMatrixElement(i,j,parameters);
		}//close for j
	}//close for i

	//print der matrix
	cout<<"Printing Derivative Matrix"<<endl;
	DerMatrix.Print();
	
	//cout<<"Printing Covariance Matrix"<<endl;
	//CovarianceMatrix.Print();
	
	TMatrixD DerMatrixTransp(Npar_fit,Npar);
	DerMatrixTransp.Transpose(DerMatrix);//transpose of DerivativeMatrix

	TMatrixD tmp(Npar,Npar_fit);
	tmp.Mult(DerMatrix,CovarianceMatrix);

	TMatrixD FinalErrMatrix(Npar,Npar);
	FinalErrMatrix.Mult(tmp,DerMatrixTransp); //store Final Correct Covariance Matrix

  cout<<"Printing Error Matrix"<<endl;
	FinalErrMatrix.Print();
	
	//get diagonal
	for(int i=0;i<Npar;i++) parErrors[i]= sqrt(FinalErrMatrix(i,i));

	return parErrors;

}//close MSTMixtureFitter::GetParamUncertainty()
*/

/*
double MSTMixtureFitter::GetDerivativeMatrixElement(int i,int j,std::vector<double> parameters){

	//index i is relative to real mass fractions
	//index j is relative to fit mass fraction params
	int si= (int)(i/fNComponents); 
	int ii= i-si*fNComponents;
	int sj= (int)(j/(fNComponents-1));
	int jj= j-sj*(fNComponents-1);

	//calculate matrix elements according to the general indexes
	double matrixValue;

	//mass fit params
	std::vector<double> MassParams;
	MassParams.clear();
	MassParams.resize(fNComponents-1);
	std::vector<double> MassFract;
	MassFract.clear();
	MassFract.resize(fNComponents);
	double ConstrFactor;

	//copy from full vector
	for(unsigned int k=0;k<MassParams.size();k++){		
		MassParams[k]= parameters[k];
	}

	MassFract[0]= MassParams[0];
  ConstrFactor= (1.-MassParams[0]);

	for(int k=1;k<fNComponents-1;k++){   
    MassFract[k]= MassParams[k]*ConstrFactor;
    ConstrFactor*= (1.-MassParams[k]);   
  }//close for i
	MassFract[fNComponents-1]= ConstrFactor;
	

	if(i<fNComponents && j<fNComponents-1){
		//fract vs fract params
		if(j==i && j==0) matrixValue= 1.;
		else if(j==i && j!=0) matrixValue= MassFract[i]/MassParams[j];
		else if(j>i) matrixValue= 0.;	
		else if(j<i) matrixValue= -MassFract[i]/(1.-MassParams[j]); 
	}//close if
	else if(i<fNComponents && j>=fNComponents-1){
		//fract vs other params
		matrixValue= 0.;	
	}
	else if(i>=fNComponents && j<fNComponents-1){
		//other params vs fract params
		matrixValue= 0.;	
	}
	else if(i>=fNComponents && j>=fNComponents-1){
		//other params vs other params
		if(j==i-1) matrixValue= 1.;
		else matrixValue= 0.;	
	}
	else{
		cerr<<"Invalid index for derivative calculation!..please check!"<<endl;
		exit(1);
	}


	return matrixValue;

}//close MSTMixtureFitter::GetDerivativeMatrixElement()
*/


/*
void MSTMixtureFitter::RandomizePar(std::vector<bool> hasToBeGenerated){

	//## Randomize fit parameters according to the constraints
	cout<<"MSTMixtureFitter::RandomizePar(): Set random starting values..."<<endl;
		
	//## Generate p
	if(fRandomizeStartFractionPar){
		for(int k=0;k<fNComponents;k++){
			//fP_start[k]= 1./(double)(fNComponents);
			fP[k]= 1./(double)(fNComponents);
		}//end loop components
	}

	//## Generate delta
	if(fRandomizeStartDeltaPar) RandomizeDeltaPar(hasToBeGenerated); 
			
	//## Generate nu
	if(fRandomizeStartNuPar) RandomizeNuPar(hasToBeGenerated); 

	//## Generate random ksi	
	if(fRandomizeStartMeanPar) RandomizeMeanPar(hasToBeGenerated);
	
	//## Generate Sigma
	if(fRandomizeStartCovariancePar) RandomizeSigmaPar(hasToBeGenerated); 
	

}//close MSTMixtureFitter::RandomizePar()
*/

/*
void MSTMixtureFitter::RandomizeMeanPar(std::vector<bool> hasToBeGenerated){

	cout<<"MSTMixtureFitter::RandomizeMeanPar(): INFO: Generating random mean...";

	//## Generate only the components present in componentIdList
	if(hasToBeGenerated.size()<=0){
		cerr<<"MSTMixtureFitter::RandomizeMeanPar(): ERROR: vector with hasToBeGenerated flag is empty!"<<endl;
		exit(1);
	}
	
	for(int j=0;j<fNDim;j++){
		
		double constraintSign= 1;
		if(j==1) constraintSign= -1;

		//double constraintSign= fMeanConstraintSign(j,0);

		double meanSigmaGeneration= 3.*sqrt(fVarianceData[j]);	
		double minBoundary= fMeanData[j]-meanSigmaGeneration;
		double maxBoundary= fMeanData[j]+meanSigmaGeneration;
				
		while(true){

			//## Generate ksi 0
			if(hasToBeGenerated[0]){
				
				bool isGoodRandomStart= false;
				while(!isGoodRandomStart){
	
					double rand= gRandom->Uniform(minBoundary,maxBoundary);

					if(fUseMeanDiffConstraint || fUseMeanConstraint){//rand value is the mean
						fMu[0](j,0)= rand;
						fMuDiff[0](j,0)= 0.;
						isGoodRandomStart= true;		
						if(fUseMeanConstraint && (fMu[0](j,0)<fKsi_min[0](j,0) || fMu[0](j,0)>fKsi_max[0](j,0))) isGoodRandomStart= false;
					}	
					
					if(fUseLocationDiffConstraint || fUseLocationConstraint){
						fKsi[0](j,0)= rand;
						fKsiDiff[0](j,0)= 0.;	
						isGoodRandomStart= true;							
						if(fUseLocationConstraint && (fKsi[0](j,0)<fKsi_min[0](j,0) || fKsi[0](j,0)>fKsi_max[0](j,0))) isGoodRandomStart= false;
					}

				}//end while

			}//close if

			//## Generate ksi diff
			for(int k=1;k<fNComponents;k++){	

				if(!hasToBeGenerated[k]) continue;

				double minBoundary= fKsiDiff_min[k](j,0);
				double maxBoundary= fKsiDiff_max[k](j,0);

				if(fUseMeanConstraint || fUseMeanDiffConstraint){//rand value is the mean
					fMuDiff[k](j,0)= gRandom->Uniform(minBoundary,maxBoundary);
					fMu[k](j,0)= fMu[0](j,0)+fMuDiff[k](j,0);
				}
				if(fUseLocationConstraint || fUseLocationDiffConstraint){
					fKsiDiff[k](j,0)= gRandom->Uniform(minBoundary,maxBoundary);
					fKsi[k](j,0)= fKsi[0](j,0)+fKsiDiff[k](j,0);
				}

			}//end loop components

			//## Check ksi order
			bool isGoodKsiSet= true;
			for(int k=0;k<fNComponents-1;k++){

				double ksi_k= 0;
				double ksi_k_1= 0;
				if(fUseMeanConstraint || fUseMeanDiffConstraint){
					ksi_k= fMu[k](j,0);
					ksi_k_1= fMu[k+1](j,0);
				}
				if(fUseLocationConstraint || fUseLocationDiffConstraint){
					ksi_k= fKsi[k](j,0);
					ksi_k_1= fKsi[k+1](j,0);
				}

				if(constraintSign*ksi_k<constraintSign*ksi_k_1){
					isGoodKsiSet= false;
					break;
				} 
			}//end loop components

			if(isGoodKsiSet) break;		

		}//end while
	}//end loop dim


	cout<<"done!"<<endl;
	cout<<"== RANDOM MEAN PARS =="<<endl;
	for(int k=0;k<fNComponents;k++){	
		cout<<"--> Component no. "<<k+1<<endl;
		cout<<"Ksi= (";
		for(int j=0;j<fNDim-1;j++) cout<<fKsi[k](j,0)<<",";
		cout<<fKsi[k](fNDim-1,0)<<")"<<endl; 

		cout<<"Mu= (";
		for(int j=0;j<fNDim-1;j++) cout<<fMu[k](j,0)<<",";
		cout<<fMu[k](fNDim-1,0)<<")"<<endl; 		
	}
	cout<<"======================"<<endl;
	cout<<"end!"<<endl;
	
		
}//close MSTMixtureFitter::RandomizeMeanPar()
*/


/*
void MSTMixtureFitter::RandomizeSigmaPar(std::vector<bool> hasToBeGenerated){


	cout<<"MSTMixtureFitter::RandomizeSigmaPar(): INFO: Generating random covariance...";

	
	//## Generate only the components present in componentIdList
	if(hasToBeGenerated.size()<=0){
		cerr<<"MSTMixtureFitter::RandomizeSigmaPar(): ERROR: vector with hasToBeGenerated flag is empty!"<<endl;
		exit(1);
	}

	while(true){
	
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++){
		
				double parValue[fNComponents];

				for(int k=0;k<fNComponents;k++){

					if(!hasToBeGenerated[k]) continue;

					fSigma[k](j,l)= gRandom->Uniform(0,fVarianceData[j]);

					if(fUseCovarianceConstraint || fUseCovarianceBoundConstraint){					
						double minBoundary= fSigma_min[k](j,l);
						double maxBoundary= fSigma_max[k](j,l);					
						fV[k](j,l)= 0.;

						if(k==0 || j!=l) fV[k](j,l)= gRandom->Uniform(minBoundary,maxBoundary);
						else if(k>0 && j==l) fV[k](j,l)= gRandom->Uniform(minBoundary,min(fV[k-1](j,l),maxBoundary));
						else{
							cerr<<"ERROR!"<<endl;
							exit(1);
						}
					}//close if
					

					if(fUseScaleMatrixConstraint || fUseScaleMatrixBoundConstraint){					
						double minBoundary= fSigma_min[k](j,l);
						double maxBoundary= fSigma_max[k](j,l);					
						fSigma[k](j,l)= 0.;

						if(k==0 || j!=l) fSigma[k](j,l)= gRandom->Uniform(minBoundary,maxBoundary);
						else if(k>0 && j==l) fSigma[k](j,l)= gRandom->Uniform(minBoundary,min(fSigma[k-1](j,l),maxBoundary));
						else{
							cerr<<"ERROR!"<<endl;
							exit(1);
						}
					}//close if
									

				}//end loop components

			}//end loop dim
	
		}//end loop dim

		//## Set equal covariances and find eigenvalues
		bool isSigmaPosDef= true;	

		if(fUseScaleMatrixConstraint || fUseScaleMatrixBoundConstraint){					
					
			for(int k=0;k<fNComponents;k++){
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						if(j==l) continue;
						fSigma[k](j,l)= fSigma[k](l,j);
					}//end loop dim
				}//end loop dim

				double Det= fSigma[k].Determinant();
				if(Det<0) {
					isSigmaPosDef= false;	
					break;
				}	
			}//end loop components

		}//close if		
	
		bool isVPosDef= true;	

		if(fUseCovarianceConstraint || fUseCovarianceBoundConstraint){					
					
			for(int k=0;k<fNComponents;k++){
				for(int j=0;j<fNDim;j++){
					for(int l=0;l<fNDim;l++){
						if(j==l) continue;
						fV[k](j,l)= fV[k](l,j);
					}//end loop dim
				}//end loop dim

				double Det= fV[k].Determinant();
				if(Det<0) {
					isVPosDef= false;	
					break;
				}	
			}//end loop components

		}//close if		

	
		//## Check eigen constraints
		bool isGoodSigmaEigen= true;

		if(fUseScaleMatrixEigenConstraint || fUseScaleMatrixEigenBoundConstraint){
			
			for(int k=0;k<fNComponents;k++){

				//## Spectral decomposition
				TMatrixDEigen SigmaDecomposition(fSigma[k]);
				TMatrixD SigmaEigenvalues(fNDim,fNDim);
				SigmaEigenvalues= SigmaDecomposition.GetEigenValues();
				fSigmaEigen[k]= SigmaEigenvalues;
				
				for(int j=0;j<fNDim;j++){
					double eigen= fSigmaEigen[k](j,j);
					double minBoundary= fSigmaEigen_min[k](j,0);
					double maxBoundary= fSigmaEigen_max[k](j,0);
					if(k==0){
						if(eigen<minBoundary || eigen>maxBoundary){
							isGoodSigmaEigen= false;
							break;
						}
					}
					else{
						if(eigen<minBoundary || eigen>maxBoundary || eigen>fSigmaEigen[k-1](j,j)){
							isGoodSigmaEigen= false;
							break;
						}
					}
				}//end loop dim
			}//end loop components
		}//close if


		bool isGoodVEigen= true;

		if(fUseCovarianceEigenConstraint || fUseCovarianceEigenBoundConstraint){
			
			for(int k=0;k<fNComponents;k++){

				//## Spectral decomposition
				TMatrixDEigen SigmaDecomposition(fV[k]);
				TMatrixD SigmaEigenvalues(fNDim,fNDim);
				SigmaEigenvalues= SigmaDecomposition.GetEigenValues();
				fVEigen[k]= SigmaEigenvalues;

				//cout<<"== Component no. "<<k+1<<" ==="<<endl;
				//cout<<"VEigen= (";
				//for(int j=0;j<fNDim;j++){
				//	for(int l=0;l<fNDim;l++) cout<<(fVEigen[k])(j,l)<<",";
				//}	
				//cout<<")"<<endl;
	

				for(int j=0;j<fNDim;j++){
					double eigen= fVEigen[k](j,j);
					double minBoundary= fSigmaEigen_min[k](j,0);
					double maxBoundary= fSigmaEigen_max[k](j,0);
					if(k==0){
						if(eigen<minBoundary || eigen>maxBoundary){
							isGoodVEigen= false;
							break;
						}
					}
					else{
						if(eigen<minBoundary || eigen>maxBoundary || eigen>fVEigen[k-1](j,j)){
							isGoodVEigen= false;
							break;
						}
					}
				}//end loop dim
			}//end loop components
		}//close if
	
		if(fUseScaleMatrixConstraint || fUseScaleMatrixBoundConstraint){
			if(isSigmaPosDef && isGoodSigmaEigen) break;
		}
	
		if(fUseCovarianceConstraint || fUseCovarianceBoundConstraint){
			if(isVPosDef && isGoodVEigen) break;
		}
		
		//cout<<"isVPosDef? "<<isVPosDef<<"  isGoodVEigen? "<<isGoodVEigen<<endl;
	
	}//close while

	
	cout<<"done!"<<endl;

}//close MSTMixtureFitter::RandomizeSigmaPar()
*/



/*
void MSTMixtureFitter::RandomizeDeltaPar(std::vector<bool> hasToBeGenerated){
	
	cout<<"MSTMixtureFitter::RandomizeDeltaPar(): INFO: Generating random delta...";

	//## Generate only the components present in componentIdList
	if(hasToBeGenerated.size()<=0){
		cerr<<"MSTMixtureFitter::RandomizeSigmaPar(): ERROR: vector with hasToBeGenerated flag is empty!"<<endl;
		exit(1);
	}

	for(int j=0;j<fNDim;j++){
		double parValue[fNComponents];

		for(int k=0;k<fNComponents;k++){
			if(!hasToBeGenerated[k]) continue;

			double minBoundary= fDelta_min[k](j,0);
			double maxBoundary= fDelta_max[k](j,0);
			//parValue[k]= gRandom->Uniform(minBoundary,maxBoundary);
			fDelta[k](j,0)= gRandom->Uniform(minBoundary,maxBoundary);
			
			//fDelta_start[k](j,0)= parValue[k];
			//fDelta[k](j,0)= parValue[k];
		}//end loop components

	}//end loop dim

		cout<<"done!"<<endl;

}//close MSTMixtureFitter::RandomizeDeltaPar()
*/

/*
void MSTMixtureFitter::RandomizeNuPar(std::vector<bool> hasToBeGenerated){

	//## Generate only the components present in componentIdList
	if(hasToBeGenerated.size()<=0){
		cerr<<"MSTMixtureFitter::RandomizeSigmaPar(): ERROR: vector with hasToBeGenerated flag is empty!"<<endl;
		exit(1);
	}

	double nuValue[fNComponents];

	for(int k=0;k<fNComponents;k++){
		if(!hasToBeGenerated[k]) continue;

		double minBoundary= fNu_min[k];
		double maxBoundary= fNu_max[k];
		//nuValue[k]= gRandom->Uniform(minBoundary,maxBoundary);
		fNu[k]= gRandom->Uniform(minBoundary,maxBoundary);	
		//fNu_start[k]= nuValue[k];
		//fNu[k]= nuValue[k];
	}//end loop components


}//close MSTMixtureFitter::RandomizeNuPar()
*/



void MSTMixtureFitter::PrintPars()
{
	for(int k=0;k<fNComponents;k++){
		cout<<"== Component "<<k+1<<" =="<<endl;
		cout<<"p= "<<fP[k]<<endl;
		cout<<"Nu= "<<fNu[k]<<endl;
		cout<<"Ksi= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fKsi[k])(0,j)<<",";
		cout<<(fKsi[k])(0,fNDim-1)<<")"<<endl; 
		cout<<"Mu= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fMu[k])(0,j)<<",";
		cout<<(fMu[k])(0,fNDim-1)<<")"<<endl; 
		cout<<"delta= (";
		for(int j=0;j<fNDim-1;j++) cout<<(fDelta[k])(0,j)<<",";
		cout<<(fDelta[k])(0,fNDim-1)<<")"<<endl; 
				
		cout<<"Omega= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fOmega[k])(j,l)<<",";
		}	
		cout<<")"<<endl;

		cout<<"Sigma= (";
		for(int j=0;j<fNDim;j++){
			for(int l=0;l<fNDim;l++) cout<<(fSigma[k])(j,l)<<",";
		}	
		cout<<")"<<endl;
	}//end loop components

}//close PrintPars()

}//close namespace
