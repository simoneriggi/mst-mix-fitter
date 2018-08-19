/**
* @file MSTMixtureFitter.h
* @class MSTMixtureFitter
* @brief MSTMixtureFitter
*
* Fit a mixture of multivariate skew-t
* @author S. Riggi
* @date 18/09/2012
*/



#ifndef _MST_MIXTURE_FITTER_h
#define _MST_MIXTURE_FITTER_h 1

#include <TObject.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>


#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

namespace MSTMixFitter_ns {

class MSTMixtureFitter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MSTMixtureFitter();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~MSTMixtureFitter();

		
		enum FitterModel {eEM=1, eMINUIT=2};
		enum FitStatus {eInit=0, eConverged=1, eTroubles=-1};
		enum ParInitMethod {
			eKMEANS= 1,
			eUSER= 2,
			eRANDOM= 3
		};

	public:
		/**
		* \brief Class destructor: free allocated memory
		*/
		int Run();

	private:
		/**
		* \brief Read input data
		*/
		int ReadData();
		/**
		* \brief Initialize data
		*/
		int Init();	
		/**
		* \brief Clear allocated data
		*/
		void ClearData();
		/**
		* \brief Run EM fitter
		*/
		int RunEMFitter();	

		/**
		* \brief Run Minuit fitter
		*/
		int RunFitter(bool initializePars=true);	

		/**
		* \brief Initialize fitter pars
		*/
		int InitializeFitPars();	

		/**
		* \brief Run EM init step
		*/
		int RunEM_Init();	

		/**
		* \brief Run EM E-step
		*/
		int RunEM_EStep(double& LL);	
		/**
		* \brief Run EM M-step
		*/
		int RunEM_MStep();
		/**
		* \brief Run EM constraint step
		*/
		int RunEM_ConstrainStep();	

		/**
		* \brief Initialize starting pars to KMEANS solution
		*/	
		int InitParsToKMeans();
		
		/**
		* \brief Init mixture parameters to user provided values
		*/
		int InitParsToUser();
		
		/**
		* \brief Print parameters
		*/	
		void PrintPars();
		/**
		* \brief Check scale matrix over EM fitting
		*/
		int CheckCovariance();
	
		/*			
		void RandomizePar(std::vector<bool> hasToBeGenerated);
		void RandomizeMeanPar(std::vector<bool> hasToBeGenerated);
		void RandomizeSigmaPar(std::vector<bool> hasToBeGenerated);
		void RandomizeDeltaPar(std::vector<bool> hasToBeGenerated);
		void RandomizeNuPar(std::vector<bool> hasToBeGenerated);
		*/

	private:

		static void LikelihoodFcn(int& nPar,double* const grad, double& value,double* const par,const int iFlag);

		/*
		double GetDerivativeMatrixElement(int i,int j,std::vector<double> parameters);
		std::vector<double> GetParamUncertainty(std::vector<double> parameters,TMatrixD CovarianceMatrix,int Npar);
		*/

	protected:
		static int fNComponents;
		static int fNDim;
		static std::vector<TMatrixD> fData;//data (each event of size 1 x fNDim)
		
	private:
		//- Main parameters
		long int fN;
		bool fIsInteractiveRun;
		int fFitter;
		
		//- Data read vars
		std::string fInputFileName;
		std::string fOutputFileName;
		std::string fDataReadDelimiter;
		std::string fRTableName;

		//- EM main options
		bool fUseStoppingCriteria;
		double fEpsilon;
		int fNIterations;	
		bool fRunMinuitFitterAtConvergence;		

		//- EM Constraint options
		//std::string fConstraintDataFileName;
		bool fUseConstraints;	
		double fConstraintAlphaScale;
		double fConstraintAlphaTolerance;
		bool fUseRandomRegenerationAfterStuck;
		bool fForceDiagonalCovariance;
		bool fFixMeanPars;
		bool fFixCovariancePars;
		bool fFixDeltaPars;
		bool fFixNuPars;
		bool fFixFractionPars;
		bool fUseMeanBoundConstraint;
		std::vector<TMatrixD> fMu_min;
		std::vector<TMatrixD> fMu_max;
		std::vector<TMatrixD> fKsi_min;
		std::vector<TMatrixD> fKsi_max;
		bool fUseCovarianceBoundConstraint;
		std::vector<TMatrixD> fSigma_min;
		std::vector<TMatrixD> fSigma_max;
		std::vector<TMatrixD> fOmega_min;
		std::vector<TMatrixD> fOmega_max;
		bool fUseCovarianceEigenBoundConstraint;			
		std::vector<TMatrixD> fSigmaEigen_min;
		std::vector<TMatrixD> fSigmaEigen_max;
		bool fUseDeltaBoundConstraint;
		std::vector<TMatrixD> fDelta_min;
		std::vector<TMatrixD> fDelta_max;
		bool fUseNuBoundConstraint;
		std::vector<double> fNu_min;
		std::vector<double> fNu_max;

		/*
		bool fUseMeanConstraint;
		bool fUseMeanDiffConstraint;
		bool fUseLocationConstraint;
		bool fUseLocationDiffConstraint;
		bool fUseCovarianceConstraint;
		bool fUseCovarianceBoundConstraint;
		bool fUseCovarianceEigenConstraint;
		bool fUseCovarianceEigenBoundConstraint;
		bool fUseScaleMatrixConstraint;
		bool fUseScaleMatrixBoundConstraint;
		bool fUseScaleMatrixEigenConstraint;
		bool fUseScaleMatrixEigenBoundConstraint;

		bool fUseDeltaBoundConstraint;
		bool fUseNuBoundConstraint;
		bool fRandomizeStartPar;
		bool fRandomizeStartMeanPar;
		bool fRandomizeStartMeanDiffPar;
		bool fRandomizeStartCovariancePar;
		bool fRandomizeStartCovarianceEigenPar;
		bool fRandomizeStartDeltaPar;
		bool fRandomizeStartNuPar;
		bool fRandomizeStartFractionPar;
		bool fUseAnnealing;
		double fAnnealingParStart;	
		double fAnnealingParStep;
		bool fConstraintFirstMeanPar;
		*/
		
	
		//- Data matrix
		std::vector<double> fMinDataRange;
		std::vector<double> fMaxDataRange;
		std::vector<double> fMeanData;
		std::vector<double> fVarianceData;

		//- EM fitted parameters
		int fFitStatus;
		long int fIterNo;
		double fLogLikelihood;
		std::vector<double> fIterLogLikelihood;
		std::vector< std::vector<double> > fTau;
		std::vector< std::vector<double> > fE1;
		std::vector< std::vector<double> > fE2;
		std::vector< std::vector<double> > fE3;
		std::vector< std::vector<double> > fE4;
		std::vector<double> fNu;//degrees of freedom (fNMixtures x 1)
		std::vector<double> fP;//mixture weights
		std::vector<TMatrixD> fKsi;//location for each mixture of size (1 x fNDim)
		
		std::vector<TMatrixD> fSigma;//covariance matrix for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fSigmaInv;
		std::vector<double> fSigmaDet;
		std::vector<TMatrixD> fSigmaEigen;
		std::vector<TMatrixD> fOmega;//scale matrix for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fOmegaInv;
		std::vector<double> fOmegaDet;
		std::vector<TMatrixD> fOmegaEigen;
		std::vector<TMatrixD> fDelta;//skewness for each mixture of size (1 x fNDim)
		std::vector<TMatrixD> fMu;//mean for each mixture of size (1 x fNDim)
		
		//- EM starting parameters
		int fParInitMethod;
		std::vector<double> fNu_start;//start degrees of freedom (fNMixtures x 1)
		std::vector<double> fP_start;//start mixture weights
		std::vector<TMatrixD> fKsi_start;//location for each mixture component of size (1 x fNDim)
		//std::vector<TMatrixD> fKsiDiff_start;//location for each mixture component of size (1 x fNDim)
		std::vector<TMatrixD> fSigma_start;//covariance for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fDelta_start;//scale for each mixture of size (1 x fNDim)
		std::vector<TMatrixD> fMu_start;//location for each mixture component of size (1 x fNDim)
		//std::vector<TMatrixD> fV_start;//covariance for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fOmega_start;//covariance for each mixture of size (fNDim x fNDim)
	
		
		std::vector<double> fNu_safe;//start degrees of freedom (fNMixtures x 1)
		std::vector<double> fP_safe;
		std::vector<TMatrixD> fKsi_safe;//location for each mixture component of size (1 x fNDim)
		std::vector<TMatrixD> fMu_safe;//mean for each mixture component of size (1 x fNDim)
		std::vector<TMatrixD> fSigma_safe;//covariance for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fSigmaEigen_safe;
		std::vector<TMatrixD> fOmega_safe;//scale matrix for each mixture of size (fNDim x fNDim)
		std::vector<TMatrixD> fOmegaEigen_safe;
		std::vector<TMatrixD> fDelta_safe;//delta par for each mixture of size (1 x fNDim)
		

	friend class MathUtils;

	
	ClassDef(MSTMixtureFitter,1)

};//close class

}//close namespace  

#endif
