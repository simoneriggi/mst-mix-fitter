/**
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/
#ifndef _CONFIGPARSER_H_
#define _CONFIGPARSER_H_


#include <TMatrixD.h>
#include <vector>
#include <string>

namespace MSTMixFitter_ns {

class ConfigParser {
  
	public:
  
		/** 
		\brief Class constructor
 		*/
  	ConfigParser();

		/** 
		\brief Class destructor
 		*/
  	virtual ~ConfigParser();
	
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		int ReadConfig(std::string filename);


	private:

		/** 
		\brief Print parsed information
 		*/
		void Print();

	public:
	
		//## File info	
		static std::string fConfigFileName;
		static std::string fInputFileName;
		static std::string fOutputFileName;
		static std::string fDataReadDelimiter;
		
		//## Algorithm info
		static int fNDim;
		static int fNComponents;
		static int fNIterations;
		static bool fRunMinuitFitterAtConvergence;
		static int fFitterChoice;	
		static bool fUseStoppingCriteria;	
		static double fEpsilon;
		
		//## Start user parameters
		static int fParInitMethod;
		static bool fRandomizeStartPars;
		static bool fRandomizeStartCovariancePars;
		static bool fRandomizeStartMeanPars;
		static std::vector<TMatrixD> fMeanStartPars;
		static std::vector<TMatrixD> fDeltaStartPars;		
		static std::vector<TMatrixD> fCovarianceStartPars;
		static std::vector<double> fFractionStartPars;	
		static std::vector<double> fNuStartPars;

		//## Constraint parameters
		static bool fFixMeanPars;
		static bool fFixDeltaPars;
		static bool fFixNuPars;
		static bool fFixCovariancePars;
		static bool fFixFractionPars;
		static bool fForceDiagonalCovariance;
		static bool fUseConstraints;
		static double fConstraintAlphaScale;
		static double fConstraintAlphaTolerance;
		static bool fUseRandomRegenerationAfterStuck;
		static bool fUseCovarianceEigenBoundConstraint;
		static std::vector<TMatrixD> fCovarianceEigenMinBound;
		static std::vector<TMatrixD> fCovarianceEigenMaxBound;
		static bool fUseCovarianceBoundConstraint;
		static std::vector<TMatrixD> fCovarianceMinBound;
		static std::vector<TMatrixD> fCovarianceMaxBound;
		static bool fUseMeanBoundConstraint;
		static std::vector<TMatrixD> fMeanMinBound;
		static std::vector<TMatrixD> fMeanMaxBound;
		static bool fUseDeltaBoundConstraint;
		static std::vector<TMatrixD> fDeltaMinBound;
		static std::vector<TMatrixD> fDeltaMaxBound;
		static bool fUseNuBoundConstraint;
		static std::vector<double> fNuMinBound;
		static std::vector<double> fNuMaxBound;

		/*
		static bool fUseDeltaLStoppingCriteria;	
		static double fDeltaLEpsilon;
	
	
		static bool fUseRandomStart;
		static bool fUseRandomFromModelStart;	
		static bool fUseRandomMeanStart;
		static bool fUseRandomMeanDiffStart;
		static bool fUseRandomCovarianceStart;
		static bool fUseRandomCovarianceEigenStart;
		static bool fUseRandomDeltaStart;
		static bool fUseRandomNuStart;
		static bool fUseRandomFractionStart;
		static bool fUseMeanImputationStart;

			
		static bool fUseMeanConstraint;
		
		static bool fUseMeanDiffConstraint;
		static bool fUseCovarianceConstraint;
		static bool fUseCovarianceEigenConstraint;
		
		static bool fUseScaleMatrixConstraint;
		static bool fUseScaleMatrixBoundConstraint;
		static bool fUseScaleMatrixEigenConstraint;
		static bool fUseScaleMatrixEigenBoundConstraint;
		static bool fUseLocationConstraint;
		static bool fUseLocationDiffConstraint;
		
		
		static TMatrixD* fMeanConstraintSign;
		static TMatrixD* fSigmaConstraintSign;
	
	
		*/

	friend class MNMixtureClustering;
		
};//close class

}//close namespace 

#endif
 
