/**
* @file MathUtils.cc
* @class MathUtils
* @brief Mathematical utility functions
*
* Useful math functions
* @author S. Riggi
* @date 17/09/2012
*/


#include <MathUtils.h>
#include <Utils.h>
#include <MSTMixtureFitter.h>

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
#include <TMatrix.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <Math/WrappedTF1.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>
#include <Math/RootFinder.h>
#include <Math/DistSampler.h>
#include <Math/DistSamplerOptions.h>
#include <Math/MinimizerOptions.h>
#include <Math/Factory.h>

#include <RInside.h>                    // for the embedded R via RInside


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(MSTMixFitter_ns::MathUtils)

namespace MSTMixFitter_ns {

//double MathUtils::fNuPar1;
//double MathUtils::fNuPar2;


MathUtils::MathUtils(){

}

MathUtils::~MathUtils(){

}



int MathUtils::MakeSymmetricMatrix(TMatrixD& C)
{
	//Check if already symmetric
	if(C.IsSymmetric()){
		INFO_LOG("Matrix is already symmetric, nothing to be done...");
		return 0;
	}

	//Force matrix symmetric
	TMatrixD C_t(TMatrixD::kTransposed,C);
	TMatrixD C_sym = 0.5*(C + C_t);
	C= C_sym;

	return 0;

}//close MakeSymmetricMatrix()


int MathUtils::ComputeSymMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,const TMatrixD& C)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();

	//Resize eigen matrix
	eigenVals.ResizeTo(1,nCols);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		ERROR_LOG("Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)");
		return -1;
	}
	
	//Check symmetric
	if(!C.IsSymmetric()){
		ERROR_LOG("Input matrix is not symmetric!");
		return -1;
	}

	//Compute eigenvalues & eigenvect
	TMatrixDEigen matrixDecomposition(C);
	TVectorD eigenVals_re= matrixDecomposition.GetEigenValuesRe();
	eigenVects= matrixDecomposition.GetEigenVectors();
	for(int i=0;i<eigenVals_re.GetNrows();i++) eigenVals(0,i)= eigenVals_re(i);

	return 0;

}//close ComputeSymMatrixEigenvalues()

int MathUtils::ComputeMatrixEigenvalues(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& C,bool forceSymmetric)
{
	//## NB: If matrix is not simmetric then the eigenvalue matrix D is block
	//       diagonal with the real eigenvalues in 1-by-1 blocks and any complex
	//       eigenvalues, a + i*b, in 2-by-2 blocks, [a, b; -b, a]
	//       If matrix is simmetric the eigenvalue matrix is diagonal with real eigenvalue on the diagonal

	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();

	//Resize eigen matrix
	eigenVals.ResizeTo(nRows,nCols);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		ERROR_LOG("Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)");
		return -1;
	}
	
	//Make symmetric
	if(forceSymmetric && MakeSymmetricMatrix(C)){
		ERROR_LOG("Failed to force simmetric matrix!");
		return -1;
	}

	//Compute eigenvalues & eigenvect
	TMatrixDEigen matrixDecomposition(C);
	eigenVals= matrixDecomposition.GetEigenValues();
	eigenVects= matrixDecomposition.GetEigenVectors();

	return 0;

}//close ComputeMatrixEigenvalues()

int MathUtils::ComputeMatrixEigenvaluesInR(TMatrixD& eigenVals,TMatrixD& eigenVects,TMatrixD& C,bool forceSymmetric)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();
	
	//Resize eigen matrix
	eigenVals.ResizeTo(nCols,1);
	eigenVects.ResizeTo(nRows,nCols);

	//Check for square matrix
	if(nCols!=nRows){
		ERROR_LOG("Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)");
		return -1;
	}

	//Import matrix in R
	std::string matrixRName= "sigma";
	if(Utils::ImportMatrixInR(&C,matrixRName)<0){
		ERROR_LOG("Failed to import matrix in R!");
		return -1;
	}

	//Force symmetric?
	if(forceSymmetric){
		try{
			Utils::fR.parseEvalQ(Form("forceSymmetric(%s);",matrixRName.c_str()));
		}
		catch(...){
			ERROR_LOG("Failed to force symmetric matrix in R!");
			return -1;
		}
	}//close if

	//Compute eigenvalues & vectors
	std::stringstream ss;
	ss<<"eigenRes <- eigen("<<matrixRName<<");";
	std::string RCmd= ss.str();
	try{
		Utils::fR.parseEvalQ(RCmd);
		Utils::fR.parseEvalQ(Form("rm(%s);",matrixRName.c_str()));//remove tmp sigma
	}
	catch(...){
		ERROR_LOG("Failed to compute eigenvalues & eigenvectors of given matrix!");
		return -1;
	}

	//Retrieve results
	try{
		Rcpp::NumericVector sigmaEigenValues= Utils::fR.parseEval("eigenRes$values;");
		Rcpp::NumericMatrix sigmaEigenVectors= Utils::fR.parseEval("eigenRes$vectors;");

		for (int l=0; l<nCols; l++) {
			eigenVals(0,l)= sigmaEigenValues(l);
			for (int j=0; j<nCols; j++) {
				double w= sigmaEigenVectors(l,j);
				eigenVects(l,j)= w;
			}
		}
		
		Utils::fR.parseEvalQ("rm(eigenRes);");

	}//close try
	catch(...){
		ERROR_LOG("Failed to retrieve eigenvalues/eigenvector matrix!");
		return -1;
	}

	return 0;
	
}//close ComputeMatrixEigenvalues()

void MathUtils::MakeDiagonalMatrix(TMatrixD& C)
{
	long int nCols= C.GetNcols();
	long int nRows= C.GetNrows();
	
	for(int j=0;j<nRows;j++){
		for(int l=0;l<nCols;l++){
			if(j==l) continue;
			C(j,l)= 0.;
		}//end loop dim
	}//end loop dim

}//close MakeDiagonalMatrix()

int MathUtils::MakeSymmPosDefCovarianceMatrix(TMatrixD& covMatrix)
{
	long int nCols= covMatrix.GetNcols();
	long int nRows= covMatrix.GetNrows();
	
	//Check for square matrix
	if(nCols!=nRows){
		ERROR_LOG("Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a square matrix!)");
		return -1;
	}

	//Import matrix in R
	std::string covMatrixRName= "sigma";
	if(Utils::ImportMatrixInR(&covMatrix,covMatrixRName)<0){
		ERROR_LOG("Failed to import covariance matrix in R!");
		return -1;
	}

	//## Force covariance to be symmetric and pos def
	//## Approximate to the nearest covariance matrix
	std::stringstream ss;
	ss<<"res <- nearPD("<<covMatrixRName<<", corr=FALSE, do2eigen=FALSE, ensureSymmetry= TRUE);";
	std::string RCmd= ss.str();

	try{
		Utils::fR.parseEvalQ(RCmd);

		//Remove tmp sigma
		Utils::fR.parseEvalQ(Form("rm(%s);",covMatrixRName.c_str()));
	}
	catch(...){
		ERROR_LOG("Failed to approximate covariance to nearest symm & pos-def matrix!");
		return -1;
	}
	
	//## Get corrected matrix and re-assign to Sigma
	try{
		Rcpp::NumericMatrix covMatrix_corr= Utils::fR.parseEval("as.matrix(res$mat);");
		for (int l=0; l<nCols; l++) {
			for (int j=0; j<nCols; j++) {
				double w= covMatrix_corr(l,j);
				covMatrix(l,j)= w;
			}
		}
	}//close try
	catch(...){
		ERROR_LOG("Failed to retrieve corrected cov matrix and update given matrix!");
		return -1;
	}

	return 0;

}//close MakeSymmPosDefCovarianceMatrix()


TMatrixD* MathUtils::GetDiagonalMatrix(TMatrixD* dataMatrix)
{
	//Check data matrix
	if(!dataMatrix){
		ERROR_LOG("Null ptr to matrix given!");
		return nullptr;
	}
	long int nCols= dataMatrix->GetNcols();
	long int nRows= dataMatrix->GetNrows();
	
	//Check for square matrix
	if(nCols!=nRows){
		ERROR_LOG("Number of matrix cols ("<<nCols<<") different from number of rows ("<<nRows<<") (hint: you must pass a simmetric matrix!)");
		return nullptr;
	}

	//Create matrix
	TMatrixD* diagMatrix= new TMatrixD(nRows,nRows);
	diagMatrix->Zero();

	//Fill diagonal values
	for(long int i=0;i<nRows;i++){
		(*diagMatrix)(i,i)= (*dataMatrix)(i,i);
	}

	return diagMatrix;

}//close GetDiagonalMatrix()

TMatrixD* MathUtils::ComputeRTableColMeans(std::string RTable,std::string colMeansRName)
{
	//Check table name
	if(RTable==""){
		ERROR_LOG("Empty R table name!");
		return nullptr;
	}

	//Compute col means matrix	
	std::stringstream ss;
	ss<<colMeansRName<<" <- colMeans("<<RTable<<")";
	std::string RCmd= ss.str(); 
	try{
		Utils::fR.parseEval(RCmd.c_str());
	}
	catch(...){
		ERROR_LOG("Failed to compute data column means in R!");
		return nullptr;
	}

	//Convert data to TMatrixD
	TMatrixD* colMeans= Utils::ConvertRVectToROOTMatrix(colMeansRName);
	if(!colMeans){
		ERROR_LOG("Failed to convert R vector to ROOT!");
		return nullptr;
	}

	return colMeans;

}//close ComputeRTableColMeans()


TMatrixD* MathUtils::ComputeCovarianceMatrixFromRTable(std::string RTable,std::string covMatrixRName)
{
	//Check table name
	if(RTable==""){
		ERROR_LOG("Empty R table name!");
		return nullptr;
	}

	//Compute covariance matrix	
	std::stringstream ss;
	ss<<covMatrixRName<<" <- cov("<<RTable<<")";
	std::string RCmd= ss.str(); 

	try{
		Utils::fR.parseEval(RCmd.c_str());
	}
	catch(...){
		ERROR_LOG("Failed to compute data covariance matrix in R!");
		return nullptr;
	}

	//Convert data to TMatrixD
	TMatrixD* covMatrix= Utils::ConvertRTableToROOTMatrix(covMatrixRName);
	if(!covMatrix){
		ERROR_LOG("Failed to convert cov matrix from R to ROOT!");
		return nullptr;
	}

	return covMatrix;

}//close ComputeCovarianceMatrixFromRTable()



int MathUtils::KsiToMu(TMatrixD& mu,const TMatrixD& ksi, const TMatrixD& delta, double nu)
{
	//## Ref: Kim et al, "Moments of random vectors with skew t distribution and their quadratic forms"
	//NB: ksi is the skew-t location par
	//    mu is the first moment
	int nRows= ksi.GetNrows();
	int nCols= ksi.GetNcols();
	mu.ResizeTo(nRows,nCols);
	
	double A= sqrt(nu/TMath::Pi())*TMath::Gamma((nu-1)/2)/TMath::Gamma(nu/2);
	mu= ksi + A*delta;
	DEBUG_LOG("mu ("<<mu.GetNrows()<<"x"<<mu.GetNcols()<<"), delta("<<delta.GetNrows()<<"x"<<delta.GetNcols()<<"), ksi("<<ksi.GetNrows()<<"x"<<ksi.GetNcols()<<")");
	
	return 0;

}//close KsiToMu()

int MathUtils::MuToKsi(TMatrixD& ksi,const TMatrixD& mu,const TMatrixD& delta, double nu)
{
	int nRows= mu.GetNrows();
	int nCols= mu.GetNcols();
	ksi.ResizeTo(nRows,nCols);
	
	double A= sqrt(nu/TMath::Pi())*TMath::Gamma((nu-1)/2)/TMath::Gamma(nu/2);
	
	ksi= mu - A*delta;
	DEBUG_LOG("mu ("<<mu.GetNrows()<<"x"<<mu.GetNcols()<<"), delta("<<delta.GetNrows()<<"x"<<delta.GetNcols()<<"), ksi("<<ksi.GetNrows()<<"x"<<ksi.GetNcols()<<")");
	
	return 0;

}//close MuToKsi()

int MathUtils::ScaleToCovarianceMatrix(TMatrixD& Sigma,const TMatrixD& Omega, const TMatrixD& delta, double nu)
{
	//## Ref: Kim et al, "Moments of random vectors with skew t distribution and their quadratic forms"
	//NB: Omega is the skew-t scale par (scale matrix)
	//    Sigma is the second moment (covariance matrix)
	int nRows= Omega.GetNrows();
	int nCols= Omega.GetNcols();
	Sigma.ResizeTo(nRows,nCols);
	TMatrixD deltaT= TMatrixD(TMatrixD::kTransposed, delta);
	
	double A= sqrt(nu/TMath::Pi())*TMath::Gamma((nu-1)/2)/TMath::Gamma(nu/2);
	Sigma= nu/(nu-2)*(Omega+deltaT*delta) - A*A*deltaT*delta;

	DEBUG_LOG("Omega ("<<Omega.GetNrows()<<"x"<<Omega.GetNcols()<<"), delta("<<delta.GetNrows()<<"x"<<delta.GetNcols()<<"), Sigma("<<Sigma.GetNrows()<<"x"<<Sigma.GetNcols()<<")");
	
	return 0;

}//close ScaleToCovarianceMatrix()


int MathUtils::CovarianceToScaleMatrix(TMatrixD& Omega,const TMatrixD& Sigma, const TMatrixD& delta, double nu)
{
	int nRows= Sigma.GetNrows();
	int nCols= Sigma.GetNcols();
	Omega.ResizeTo(nRows,nCols);

	TMatrixD deltaT= TMatrixD(TMatrixD::kTransposed, delta);
	double A= sqrt(nu/TMath::Pi())*TMath::Gamma((nu-1)/2)/TMath::Gamma(nu/2);
	Omega= (nu-2)/nu*(Sigma + A*A*deltaT*delta) - deltaT*delta;

	DEBUG_LOG("Omega ("<<Omega.GetNrows()<<"x"<<Omega.GetNcols()<<"), delta("<<delta.GetNrows()<<"x"<<delta.GetNcols()<<"), Sigma("<<Sigma.GetNrows()<<"x"<<Sigma.GetNcols()<<")");
	
	return 0;

}//close CovarianceToScaleMatrix()


double MathUtils::MST_a(const TMatrixD& SigmaInv,const TMatrixD& delta)
{
	//## SigmaInv: size (p x p)
	//## delta: size (1 x p)
	TMatrixD delta_t= TMatrixD(TMatrixD::kTransposed, delta);
	TMatrixD prod= delta*SigmaInv*delta_t;
	double value= (1.+prod(0,0))/2.;
	
	return value;

}//close MST_a()

double MathUtils::MST_b(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& delta,const TMatrixD& ksi)
{
	//## SigmaInv: size (p x p)
	//## delta: size (1 x p)
	//## y: size (1 x p)
	//## ksi: size (1 x p)
	TMatrixD prod= (y-ksi)*SigmaInv*TMatrixD(TMatrixD::kTransposed, delta);
	double value= -prod(0,0)/2.;

	return value;

}//close MST_b()

double MathUtils::MST_c(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& ksi,double nu)
{	
	//## SigmaInv: size (p x p)
	//## y: size (1 x p)
	//## ksi: size (1 x p)
 	TMatrixD prod= (y-ksi)*SigmaInv*TMatrixD(TMatrixD::kTransposed, y-ksi);
	double value= prod(0,0)/2. + nu/2;
	
	return value;

}//close MST_c()

double MathUtils::MST_C(double nu,double p,double SigmaDet)
{	
	double value= pow(nu/2,nu/2)/( pow(2*TMath::Pi(),(p+1)/2.) * pow(SigmaDet,1./2.) * TMath::Gamma(nu/2.) );
	return value;

}//close MST_C()


double MathUtils::MST_D(const TMatrixD& y,const TMatrixD& SigmaInv,const TMatrixD& ksi,const TMatrixD& delta,double nu)
{
	//## D= c − b^2/a
	double value= MST_c(y,SigmaInv,ksi,nu) - pow(MST_b(y,SigmaInv,delta,ksi),2)/MST_a(SigmaInv,delta);
	
	return value;

}//close MST_D()

double MathUtils::I1(double E,double alpha)
{
	double hyperGeom= ROOT::Math::hyperg(1./2.,E,3./2.,-alpha*alpha);
	double integralValue= (TMath::Pi()*TMath::Gamma(2*E-1))/(pow(TMath::Gamma(E),2)*pow(2,2*E-1))-alpha*hyperGeom; 

	return integralValue;

}//close I1()

double MathUtils::I1Num(double E,double alpha)
{	
	//##################################################################
	//## Perform numerical integration of (1+x^2)^-E in [alpha,+inf]
	//##################################################################
	//## Create the parametric function
	double par[1]= {E};
	
	IntegrandI1 fcn;
	fcn.SetParameters(par);

  //## Create the Integrator
  ROOT::Math::Integrator integrator(ROOT::Math::IntegrationOneDim::kADAPTIVE);

  //## Set parameters of the integration
  integrator.SetFunction(fcn, false);
	integrator.SetAbsTolerance(1.e-12);
	integrator.SetRelTolerance(1.e-12);
 	
	double integralValue= integrator.IntegralUp(alpha);
  
	return integralValue;
	
}//close I1Num()

double MathUtils::I2Num(double E,double alpha)
{	
	//##################################################################
	//## Perform numerical integration of x(1+x^2)^-E in [alpha,+inf]
	//##################################################################

	//## Create the parametric function
	double par[1]= {E};
	
	IntegrandI2 fcn;
	fcn.SetParameters(par);

  //## Create the Integrator
  ROOT::Math::Integrator integrator(ROOT::Math::IntegrationOneDim::kADAPTIVE);

  //## Set parameters of the integration
  integrator.SetFunction(fcn, false);
	integrator.SetAbsTolerance(1.e-12);
	integrator.SetRelTolerance(1.e-12);
 	
	double integralValue= integrator.IntegralUp(alpha);
  
	return integralValue;
	

}//close I2Num()

double MathUtils::I2(double E,double alpha){

	double integralValue= -pow(1+alpha*alpha,1-E)/(2*(1-E));
	
	return integralValue;

}//close I2()


double MathUtils::I3(double E,double alpha)
{
	//#########################################################################
	//## Perform numerical integration of log(1+x^2)(1+x^2)^-E in [alpha,+inf]
	//#########################################################################

	//## Create the parametric function
	double par[1]= {E};
	
	IntegrandI3 fcn;
	fcn.SetParameters(par);

  //## Create the Integrator
  ROOT::Math::Integrator integrator(ROOT::Math::IntegrationOneDim::kADAPTIVE);

  //## Set parameters of the integration
  integrator.SetFunction(fcn, false);
	integrator.SetAbsTolerance(1.e-12);
	integrator.SetRelTolerance(1.e-12);
 	
	double integralValue= integrator.IntegralUp(alpha);
  
	return integralValue;

}//close I3



double MathUtils::MST_e1(double a,double b,double D,double R)
{
	double arg= b/sqrt(D*a);
	//double num= R*I1(R+1,arg);
	//double denom= D*I1(R,arg);
	double num= R*I1Num(R+1,arg);
	double denom= D*I1Num(R,arg);
	double value= num/denom;
	
	return value;

}//close MST_e1()

double MathUtils::MST_e2(double a,double b,double D,double R)
{
	double arg= b/sqrt(D*a);
	double num= R*I2(R+1,arg);
	//double denom= sqrt(D*a)*I1(R,arg);
	double denom= sqrt(D*a)*I1Num(R,arg);
	double value= num/denom - b/a*MST_e1(a,b,D,R);
	
	return value;

}//close MST_e2()


double MathUtils::MST_e3(double a,double b,double c,double D,double R)
{
	double value= (R-2*b*MST_e2(a,b,D,R)-c*MST_e1(a,b,D,R))/a;
	
	return value;

}//close MST_e3()

double MathUtils::MST_e4(double a,double b,double D,double R)
{
	double arg= b/sqrt(D*a);

	//## Compute digamma function digamma(R)
	double value= 0;
	try{
		Utils::fR["x"]= R;
		Rcpp::NumericVector digammaVect = Utils::fR.parseEval(std::string("digamma(x)"));
		double digamma= digammaVect[0];

		//value= digamma - log(D) + I3(R,arg)/I1(R,arg);
		value= digamma - log(D) + I3(R,arg)/I1Num(R,arg);
	}
	catch(...){
		ERROR_LOG("Failed to compute digamma function, returning zero!");
		return 0;
	}

	return value;

}//close MST_e4()

double MathUtils::MSTDensityFcn(double a,double b,double C,double D,double R)
{
	double arg= b/sqrt(D*a);
	//double value= 2*C*TMath::Gamma(R)*pow(D,1./2.-R)/sqrt(a)*I1(R,arg);
	double value= 2*C*TMath::Gamma(R)*pow(D,1./2.-R)/sqrt(a)*I1Num(R,arg);

	return value;

}//close MSTDensityFcn()


double MathUtils::MST_tau(double a,double b,double C,double D,double R,double pi,double beta)
{
	double f= MSTDensityFcn(a,b,C,D,R);
	double value= pow(pi*f,beta);

	return value;

}//close MST_tau()

/*
double MathUtils::NuEquationFcn(double x) 
{ 
	double nu= x;
  double tauSum= fNuPar1;
	double k= fNuPar2;
	double digamma= 0.;
	try{
		Utils::fR["nu_half"]= nu/2.;
		Rcpp::NumericVector digammaVect = Utils::fR.parseEval(std::string("digamma(nu_half)"));
		digamma= digammaVect[0];
	}
	catch(...){
		ERROR_LOG("Failed to compute digamma function in R, returning 0!");
		return 0;
	}

	double fcnValue= (log(nu/2.) - digamma + 1)* tauSum + k;
	DEBUG_LOG("nu="<<nu<<"  tauSum="<<tauSum<<"  k="<<k<<"  log(nu/2.)="<<log(nu/2.)<<"  digamma="<<digamma<<"  fcn="<<fcnValue);

  return fcnValue;	
		
}//close NuEquationFcn()
 
double MathUtils::NuEquationDerivativeFcn(double x) 
{ 
	double nu= x;
  double tauSum= fNuPar1;
	//double k= fNuPar2;
	double trigamma= 0.;
	try{
		Utils::fR["nu_half"]= nu/2.;
		Rcpp::NumericVector trigammaVect = Utils::fR.parseEval(std::string("trigamma(nu_half)"));
		trigamma= trigammaVect[0];
	}
	catch(...){
		ERROR_LOG("Failed to compute trigamma function in R, returning 0!");
		return 0;
	}
	double fcnValue= (1./nu-1./2.*trigamma)*tauSum/2.;
			
	return fcnValue;	

}//close NuEquationDerivativeFcn()
*/

int MathUtils::NuSolver(double& rootValue,double p1,double p2,double start,int maxIter,double relTolerance,double absTolerance)
{
	//## Create the parametric function
	double par[2]= {p1,p2};
	//fNuPar1= p1;
	//fNuPar2= p2;
	
	NuEquation fcn;
	fcn.SetParameters(par);

	//## Search proper range such that fcn has opposite signs
	bool hasFoundRootRange= false;
	double minNu= 1.e-6;
	double maxNu= 100000;
	double stepNu= 1;
	int nsearch= (int)((maxNu-minNu)/stepNu);
	
	double startFcnValue= fcn.DoEvalPar(minNu,par);	
	double endNu= minNu;
	double currentFcnVal= startFcnValue;
	double currentNu= minNu;
	
	for(int i=0;i<nsearch;i++)
	{
		currentFcnVal= fcn.DoEvalPar(currentNu,par);
		DEBUG_LOG("currentNu="<<currentNu<<"  currentFcnVal="<<currentFcnVal<<"  startFcnValue="<<startFcnValue);
		if(currentFcnVal*startFcnValue==-1){
			hasFoundRootRange= true;
			endNu= currentNu;
			break;
		}
		currentNu+= stepNu;
		
	}//end search

	if(!hasFoundRootRange){
		ERROR_LOG("Cannot determine proper range for root finder!");
		return -1;
	}
	DEBUG_LOG("Root finder range ["<<minNu<<","<<endNu<<"]  fcn=["<<startFcnValue<<","<<currentFcnVal<<"] ...");
	

	//## Create the RootFinder
  //ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kGSL_NEWTON);//do not need range as based on derivative
	//ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kGSL_BISECTION);//need range
	ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);//need range
	rootFinder.SetFunction(fcn,minNu,endNu);
	//rootFinder.SetFunction(fcn,start);

	//## Solve and get root	
	rootFinder.Solve(maxIter,absTolerance,relTolerance);
	rootValue= rootFinder.Root();
  
	return 0;

}//close NuSolver()


void MathUtils::MST_LikelihoodFcn(int& nPar,double* const grad, double& value,double* const par,const int iFlag) 
{  
  value= 0;
	int par_counter=0;
	int nComponents= MSTMixtureFitter::fNComponents;
	int nDim= MSTMixtureFitter::fNDim;
	long int N= static_cast<long int>((MSTMixtureFitter::fData).size());

	//Set p
	DEBUG_LOG("Get fraction pars...");
	double fractionParList[nComponents];
	for(int k=0;k<nComponents;k++) fractionParList[k]= 1;

	if(nComponents>1){
		//***************************************************
  	//****  METHOD TO COSTRAIN THE FRACTIONS IN THE FIT
  	//***************************************************
  	// Unconstrained fractions are only constrained in [0,1] in Minuit set parameters
  	// Real fractions must sum to unity. Use method in "Statistical Methods in Data Analysis - W. J. Metzger"
 		// fract0= fract0_uncostr
  	// fract1= fract1_uncostr * (1-fract0_uncostr)
  	// fract2= fract2_uncostr * (1-fract0_uncostr) * (1-fract1_uncostr)
  	// ...
  	// fractN= (1-fract0_uncostr) * (1-fract1_uncostr) * ... * (1-fractN-1_uncostr)
  	double UnconstrainedFraction[nComponents-1];
		for(int i=0; i< nComponents-1; i++) UnconstrainedFraction[i]=0.;

  	UnconstrainedFraction[0]= *(par+0);
  	fractionParList[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
  	for(int i=1; i< nComponents-1; i++){
  		UnconstrainedFraction[i]= *(par+i);
    	fractionParList[i]= UnconstrainedFraction[i]*ConstrFactor;
    	ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}
  	fractionParList[nComponents-1]= ConstrFactor;		

		par_counter+= nComponents-1;
	}//close if

	//set ksi 
	DEBUG_LOG("Get Ksi pars...");
	std::vector<TMatrixD> ksiParList;
	for(int k=0;k<nComponents;k++){
		ksiParList.push_back(TMatrixD(1,nDim));		
		for(int j=0;j<nDim;j++){
			ksiParList[k](0,j)= *(par+par_counter);	
			par_counter++;
		}//end loop dim
	}//end loop components

	//set Omega
	DEBUG_LOG("Get Omega pars...");
	std::vector<TMatrixD> OmegaParList;
	std::vector<TMatrixD> OmegaInvParList;
	double OmegaDet[nComponents];

	for(int k=0;k<nComponents;k++){		
		OmegaParList.push_back(TMatrixD(nDim,nDim));
		OmegaInvParList.push_back(TMatrixD(nDim,nDim));
		for(int j=0;j<nDim;j++){
			for(int l=j;l<nDim;l++){
				double w= *(par+par_counter);
				OmegaParList[k](j,l)= w;
				OmegaParList[k](l,j)= w;
				par_counter++;
			}
		}//end loop dims
		OmegaInvParList[k]= TMatrixD(TMatrixD::kInverted,OmegaParList[k]);
		OmegaDet[k]= OmegaParList[k].Determinant();
		if(OmegaDet[k]<=0){
			WARN_LOG("Covariance inversion failed (det="<<OmegaDet[k]<<")!");
		}	
		
	}//end loop components
	
	//Set delta
	DEBUG_LOG("Get delta pars...");
	std::vector<TMatrixD> deltaParList;
	for(int k=0;k<nComponents;k++){
		deltaParList.push_back(TMatrixD(1,nDim));	
		for(int j=0;j<nDim;j++){
			deltaParList[k](0,j)= *(par+par_counter);	
			par_counter++;
		}//end loop dim
	}//end loop components
	
	//set nu
	DEBUG_LOG("Get nu pars...");	
	double nuParList[nComponents];
	for(int k=0;k<nComponents;k++){
		nuParList[k]= *(par+par_counter);
		par_counter++;
	}	
	
	//## Compute LL
	DEBUG_LOG("Compute LL...");
	
	double likelihood= 0.;

	for(long int i=0;i<N;i++){
		double tauSum= 0.;			
		for(int k=0;k<nComponents;k++){

			//## Compute a, b, c, C, D, R
			double a= MathUtils::MST_a( OmegaInvParList[k],deltaParList[k]);
			double b= MathUtils::MST_b( MSTMixtureFitter::fData[i],OmegaInvParList[k],deltaParList[k],ksiParList[k]);
			double c= MathUtils::MST_c( MSTMixtureFitter::fData[i],OmegaInvParList[k],ksiParList[k],nuParList[k]);
			double C= MathUtils::MST_C( nuParList[k],nDim,OmegaDet[k]);
			double D= MathUtils::MST_D( MSTMixtureFitter::fData[i],OmegaInvParList[k],ksiParList[k],deltaParList[k],nuParList[k]);
			double R= MathUtils::MST_R( nuParList[k],nDim);
				
			//## Compute e1, e2, e3, e4, tau
			double tauComponent= MathUtils::MST_tau(a,b,C,D,R,fractionParList[k],1);
			tauSum+= tauComponent;	

		}//end loop mixture components
			 
		likelihood+= log(tauSum);

	}//end loop events

	value= -likelihood;
	//value= likelihood;

	cout<<"*************************"<<endl;
  cout<<"****  CURRENT FIT     ***"<<endl;
  cout<<"**************************"<<endl;
	cout<<"L="<<likelihood<<endl;
	cout<<"Fractions (";
	for(int k=0;k<nComponents;k++){	
		cout<<fractionParList[k]<<",";
	}//end loop components
	cout<<")"<<endl;

	for(int k=0;k<nComponents;k++){	
		cout<<"Ksi"<<k<<" (";
		for(int j=0;j<nDim;j++){
			cout<<ksiParList[k](0,j)<<",";
		}
		cout<<")"<<endl;
	}

	for(int k=0;k<nComponents;k++){	
		cout<<"Omega"<<k<<"= (";
		for(int j=0;j<nDim;j++){
			for(int l=0;l<nDim;l++){
				cout<<OmegaParList[k](j,l)<<",";
			}
		}
		cout<<")"<<endl;
	}

	for(int k=0;k<nComponents;k++){		
		cout<<"Delta"<<k<<"= (";
		for(int j=0;j<nDim;j++){
			cout<<deltaParList[k](0,j)<<",";
		}
		cout<<")"<<endl;
	}
	cout<<"Nu (";
	for(int k=0;k<nComponents;k++){	
		cout<<nuParList[k]<<",";
	}//end loop components
	cout<<")"<<endl;
	cout<<"*************************"<<endl;
	
}//close MST_LikelihoodFcn()

}//close namespace 
