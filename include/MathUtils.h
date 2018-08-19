/**
* @file MathUtils.h
* @class MathUtils
* @brief Mathematical utility functions
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/

#ifndef _MATH_UTILS_h
#define _MATH_UTILS_h 1

#include <Logger.h>
#include <Utils.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVector3.h>

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"



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

class MathUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MathUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MathUtils();

		

	public:
		
		/**
		* \brief Force input matrix to be symmetric
		*/
		static int MakeSymmetricMatrix(TMatrixD& C);

		
		/**
		* \brief Compute symmetric matrix eigenvectors and eigenvalues
		*/
		static int ComputeSymMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,const TMatrixD& inputMatrix);

		/**
		* \brief Compute matrix eigenvectors and eigenvalues
		*/
		static int ComputeMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& inputMatrix,bool forceSymmetric=false);

		/**
		* \brief Compute matrix eigenvectors and eigenvalues
		*/
		static int ComputeMatrixEigenvaluesInR(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& inputMatrix,bool forceSymmetric=false);

		/**
		* \brief Force input matrix to be diagonal
		*/
		static void MakeDiagonalMatrix(TMatrixD& C);

		/**
		* \brief Correct covariance matrix
		*/
		static int MakeSymmPosDefCovarianceMatrix(TMatrixD& covMatrix);

		/**
		* \brief Get diagonal matrix
		*/
		static TMatrixD* GetDiagonalMatrix(TMatrixD* dataMatrix);


		/**
		* \brief Compute covariance matrix of a stored R table and return it as a ROOT TMatrix
		*/
		static TMatrixD* ComputeCovarianceMatrixFromRTable(std::string RTable,std::string covMatrixRName="covMatrix");

		/**
		* \brief Compute R table column means and return it as a ROOT TMatrix
		*/
		static TMatrixD* ComputeRTableColMeans(std::string RTable,std::string colMeansRName="colMeansVect");


		/**
		* \brief Convert skew-t ksi par to mean 
		*/
		static int KsiToMu(TMatrixD& mu,const TMatrixD& ksi, const TMatrixD& delta, double nu);

		/**
		* \brief Convert skew-t mean par to ksi 
		*/
		static int MuToKsi(TMatrixD& ksi,const TMatrixD& mu, const TMatrixD& delta, double nu);
		/**
		* \brief Convert skew-t scale matrix par to covariance matrix 
		*/
		static int ScaleToCovarianceMatrix(TMatrixD& Sigma,const TMatrixD& Omega, const TMatrixD& delta, double nu);
		/**
		* \brief Convert covariance matrix to skew-t scale matix par
		*/
		static int CovarianceToScaleMatrix(TMatrixD& Omega,const TMatrixD& Sigma, const TMatrixD& delta, double nu);

		/**
		* \brief Skew-t computation par: a (see ref paper)
		*/
		static double MST_a(const TMatrixD& SigmaInv,const TMatrixD& delta);
		/**
		* \brief Skew-t computation par: b (see ref paper)
		*/
		static double MST_b(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& delta,const TMatrixD& ksi);
		/**
		* \brief Skew-t computation par: c (see ref paper)
		*/
		static double MST_c(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& ksi,double nu);
		/**
		* \brief Skew-t computation par: C (see ref paper)
		*/
		static double MST_C(double nu,double p,double SigmaDet);
		/**
		* \brief Skew-t computation par: D (see ref paper)
		*/
		static double MST_D(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& ksi,const TMatrixD& delta,double nu);
		/**
		* \brief Skew-t computation par: R (see ref paper)
		*/
		static double MST_R(double nu,double p){return (nu+p+1)/2;}

		/**
		* \brief Compute integral of (1+x^2)^-E in [alpha,+inf] analytically
		*/
		double I1(double E,double alpha);

		/**
		* \brief Compute integral of (1+x^2)^-E in [alpha,+inf] numerically
		*/
		static double I1Num(double E,double alpha);
		/**
		* \brief Compute integral of x(1+x^2)^-E in [alpha,+inf] analytically
		*/
		static double I2(double E,double alpha);
	
		/**
		* \brief Compute integral of x(1+x^2)^-E in [alpha,+inf] numerically
		*/
		static double I2Num(double E,double alpha);

		/**
		* \brief Compute integral of log(1+x^2)(1+x^2)^-E in [alpha,+inf] numerically
		*/
		static double I3(double E,double alpha);
		/**
		* \brief Skew-t computation par: e1 (see ref paper)
		*/
		static double MST_e1(double a,double b,double D,double R);
		/**
		* \brief Skew-t computation par: e2 (see ref paper)
		*/
		static double MST_e2(double a,double b,double D,double R);
		/**
		* \brief Skew-t computation par: e3 (see ref paper)
		*/
		static double MST_e3(double a,double b,double c,double D,double R);
		/**
		* \brief Skew-t computation par: e4 (see ref paper)
		*/
		static double MST_e4(double a,double b,double D,double R);

		/**
		* \brief Skew-t density function
		*/
		static double MSTDensityFcn(double a,double b,double C,double D,double R);

		/**
		* \brief Skew-t computation par: tau (see ref paper)
		*/
		static double MST_tau(double a,double b,double C,double D,double R,double pi,double beta=1);

		/**
		* \brief Solve nu equation for roots
		*/
		static int NuSolver(double& rootValue,double p1,double p2,double start,int maxIter=1000,double relTolerance=1.e-12,double absTolerance= 1.e-12);

		/**
		* \brief Define nu equation to be solved
		*/
		//static double NuEquationFcn(double x);
		/**
		* \brief Define nu equation derivative to be solved
		*/
		//static double NuEquationDerivativeFcn(double x);

		/**
		* \brief MST maximum likelihood function for MINUIT fitting
		*/
		static void MST_LikelihoodFcn(int& nPar,double* const grad, double& value,double* const par,const int iFlag);

		/**
		* \brief Sign operator
		*/
		template <typename T> 
		static int sign(T val) {
    	return (T(0) < val) - (val < T(0));
		}

	private:
	
		//static double fNuPar1;
		//static double fNuPar2;


		friend class MSTMixtureFitter;

	ClassDef(MathUtils,1)

};//close class



class IntegrandI1: public ROOT::Math::IParametricFunctionOneDim
{
	private:
  	const double *pars;
 
	public:
  	double DoEvalPar(double x,const double* p) const
 		{
    	return pow(1+x*x,-p[0]);
   	}
 
   	ROOT::Math::IBaseFunctionOneDim* Clone() const
   	{
    	return new IntegrandI1();
   	}
 
   	const double* Parameters() const
   	{
    	return pars;
   	}
 
   	void SetParameters(const double* p)
   	{
    	pars = p;
   	}
 
   	unsigned int NPar() const
   	{
    	return 1;
   	}

};//close class IntegrandI1()

class IntegrandI2: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
      return x*pow(1+x*x,-p[0]);
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new IntegrandI2();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 1;
   }
};

class IntegrandI3: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
      return log(1+x*x)*pow(1+x*x,-p[0]);
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new IntegrandI3();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 1;
   }
};


class NuEquation: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
			//fcn= (log(x/2)-psi(x/2)+1)*tauSum + k
			double tauSum= p[0];
			double k= p[1];
			double digamma= 0.;
			try{
				Utils::fR["nu_half"]= x/2.;
				Rcpp::NumericVector digammaVect = Utils::fR.parseEval(std::string("digamma(nu_half)"));
				digamma= digammaVect[0];
			}
			catch(...){
				ERROR_LOG("Failed to compute digamma function in R, returning 0!");
				return 0;
			}

			double fcnValue= (log(x/2.) - digamma + 1)* tauSum + k;
			DEBUG_LOG("nu="<<x<<"  tauSum="<<tauSum<<"  k="<<k<<"  log(x/2.)="<<log(x/2.)<<"  digamma="<<digamma<<"  fac="<<(log(x/2.) - digamma + 1)<<"  fcn="<<fcnValue);
      return fcnValue;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new NuEquation();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 2;
   }
};


}//close namespace 

#endif

