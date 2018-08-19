/**
* @file KMeansClustering.cc
* @class KMeansClustering
* @brief KMeansClustering
*
* Perform kmeans clustering
* @author S. Riggi
* @date 30/07/2013
*/

#include <KMeansClustering.h>
#include <DataReader.h>
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
#include <TColor.h>

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

ClassImp(MSTMixFitter_ns::KMeansClustering)

namespace MSTMixFitter_ns {


KMeansClustering::KMeansClustering()
{
	fClusterIds= 0;
	fClusterCenters= 0;
	fClusterSkewness= 0;
	fClusterSizes= 0;
	fClusterCovMatrixes.clear();
}

KMeansClustering::~KMeansClustering()
{
	//Delete allocated data
	ClearData();

}//close destructor

void KMeansClustering::ClearData()
{
	if(fClusterSizes){
		delete fClusterSizes;
		fClusterSizes= 0;
	}
	if(fClusterIds){
		delete fClusterIds;
		fClusterIds= 0;
	}

	if(fClusterCenters){
		delete fClusterCenters;
		fClusterCenters= 0;
	}

	if(fClusterSkewness){
		delete fClusterSkewness;
		fClusterSkewness= 0;
	}

	for(size_t i=0;i<fClusterCovMatrixes.size();i++){
		if(fClusterCovMatrixes[i]){
			delete fClusterCovMatrixes[i];		
			fClusterCovMatrixes[i]= 0;
		}
	}
	fClusterCovMatrixes.clear();

}//close Clear()

int KMeansClustering::RunKMeans(std::string RDataName,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(nComponents<=0){
		ERROR_LOG("Invalid numbr of components arg given (must be >0)");
		return -1;
	}

	//Clear existing data
	ClearData();

	//Load moments package
	std::vector<std::string> RLibraries{"moments"};
	if(Utils::LoadRLibraries(RLibraries)<0){
		ERROR_LOG("Failed to load one or more of these R libraries {moments}, check if they are installed!");
		return -1;
	}
	
	
	//Define commands to be run in R
	std::stringstream ss;
	ss<<"clusterResults <- kmeans("<<RDataName<<", centers="<<nComponents<<", iter.max="<<maxIter<<", nstart="<<nRandomInit<<");";
	std::string RCmd= ss.str();

	//Run kmeans clustering
	INFO_LOG("Running kmeans command in R: "<<RCmd);
	try{
		Utils::fR.parseEval(RCmd);
	}
	catch(...){
		ERROR_LOG("Kmeans run failed in R!");
		return -1;
	}

	//## Retrieve parameters	
	//- Cluster membership: vector of integers (from 1:k) indicating the cluster to which each point is allocated
	INFO_LOG("Retrieve kmeans clusters ...");
	fClusterIds= Utils::ConvertRVectToROOTMatrix("clusterResults$cluster");
	if(!fClusterIds){
		ERROR_LOG("Failed to retrieve cluster membership vector and convert it to ROOT!");
		return -1;
	}

	//- Centers: a matrix of cluster centers
	INFO_LOG("Retrieve kmeans cluster centers ...");
	fClusterCenters= Utils::ConvertRTableToROOTMatrix("clusterResults$centers");
	if(!fClusterCenters){
		cerr<<"ERROR: Failed to retrieve cluster center matrix and convert it to ROOT!"<<endl;
		return -1;
	}

	//- Number of data in each cluster
	INFO_LOG("Retrieve kmeans cluster centers ...");
	fClusterSizes= Utils::ConvertRVectToROOTMatrix("clusterResults$size");
	if(!fClusterSizes){
		ERROR_LOG("Failed to retrieve cluster size vector and convert it to ROOT!");
		return -1;
	}

	//## Compute sample covariance matrix for each clusters
	INFO_LOG("Computing sample covariance of each clusters...");
	for(int i=0;i<nComponents;i++){
		//Define cmd
		//cov(data[which(clusterResults$cluster==clusterId),])
		std::string covMatrixName= Form("covMatrix_cluster%d",i+1);
		RCmd= Form("%s <- cov(%s[which(clusterResults$cluster==%d),])",covMatrixName.c_str(),RDataName.c_str(),i+1);
		DEBUG_LOG("Running cmd: "<<RCmd);

		//Exec cmd
		try{
			Utils::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			ERROR_LOG("Failed to compute covariance matrix for data belonging to cluster no. "<<i+1<<"!");
			return -1;				
		}

		//Fill TMatrixD
		TMatrixD* clusterCovMatrix= Utils::ConvertRTableToROOTMatrix(covMatrixName);
		if(!clusterCovMatrix){
			ERROR_LOG("Failed to retrieve cluster cov matrix and convert it to ROOT!");
			return -1;
		}		
		fClusterCovMatrixes.push_back(clusterCovMatrix);
	}//end loop components

	//## Compute sample skewness for each clusters
	INFO_LOG("Computing sample skewness of each clusters...");
	for(int i=0;i<nComponents;i++){
		//Define cmd
		std::string skewnessVectorName= Form("skewnessVec_cluster%d",i+1);
		RCmd= Form("%s <- skewness(%s[which(clusterResults$cluster==%d),])",skewnessVectorName.c_str(),RDataName.c_str(),i+1);
		DEBUG_LOG("Running cmd: "<<RCmd);

		//Exec cmd
		try{
			Utils::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			ERROR_LOG("Failed to compute skewness for data belonging to cluster no. "<<i+1<<"!");
			return -1;				
		}

		//Retrieve skewness vector
		fClusterSkewness= Utils::ConvertRVectToROOTMatrix(skewnessVectorName);
		if(!fClusterSkewness){
			ERROR_LOG("Failed to retrieve cluster skewness vector and convert it to ROOT!");
			return -1;
		}
	}//end loop components


	return 0;

}//close RunKMeans()


int KMeansClustering::RunKMeans(TMatrixD* dataMatrix,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return -1;
	}
	if(nComponents<=0){
		ERROR_LOG("Invalid numbr of components arg given (must be >0)");
		return -1;
	}

	//## Import data matrix in R
	if(Utils::ImportMatrixInR(dataMatrix,"dataMatrix")<0){
		ERROR_LOG("Failed to import data matrix in R!");
		return -1;
	}

	//## Run kmeans on imported data
	return RunKMeans("dataMatrix",nComponents,maxIter,nRandomInit);

}//close RunKMeans()


int KMeansClustering::RunKMedians(std::string RDataName,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(nComponents<=0){
		cerr<<"ERROR: Invalid numbr of components arg given (must be >0)"<<endl;
		return -1;
	}

	//Clear existing data
	ClearData();

	//Load library flexclust
	std::vector<std::string> RLibNames= {"flexclust","moments"};
	if(Utils::LoadRLibraries(RLibNames)<0){
		ERROR_LOG("Failed to load R libraries {flexclust,moments}, check if they are installed!");
		return -1;
	}
	
	//Define commands to be run in R
	std::stringstream ss;
	ss<<"clusterResults <- kcca("<<RDataName<<",k="<<nComponents<<", family=kccaFamily(\"kmedians\"), save.data=TRUE)";
	std::string RCmd= ss.str();

	//Run kmeans clustering
	INFO_LOG("Running kmedians command in R: "<<RCmd);
	try{
		Utils::fR.parseEval(RCmd);
	}
	catch(...){
		ERROR_LOG("Kmedians run failed in R!");
		return -1;
	}

	//## Retrieve parameters	
	//- Cluster membership: vector of integers (from 1:k) indicating the cluster to which each point is allocated
	INFO_LOG("Retrieve kmedians clusters ...");
	fClusterIds= Utils::ConvertRVectToROOTMatrix("attributes(clusterResults)$cluster");
	if(!fClusterIds){
		ERROR_LOG("Failed to retrieve cluster membership vector and convert it to ROOT!");
		return -1;
	}

	//- Centers: a matrix of cluster centers
	INFO_LOG("Retrieve kmeans cluster centers ...");
	fClusterCenters= Utils::ConvertRTableToROOTMatrix("attributes(clusterResults)$centers");
	if(!fClusterCenters){
		ERROR_LOG("Failed to retrieve cluster center matrix and convert it to ROOT!");
		return -1;
	}

	//- Number of data in each cluster
	INFO_LOG("Retrieve kmeans cluster centers ...");
	fClusterSizes= Utils::ConvertRVectToROOTMatrix("attributes(clusterResults)$clusinfo$size");
	if(!fClusterSizes){
		ERROR_LOG("Failed to retrieve cluster size vector and convert it to ROOT!");
		return -1;
	}

	//## Compute sample covariance matrix for each clusters
	INFO_LOG("Computing sample covariance of each clusters...");
	for(int i=0;i<nComponents;i++){
		//Define cmd
		//cov(data[which(clusterResults$cluster==clusterId),])
		std::string covMatrixName= Form("covMatrix_cluster%d",i+1);
		RCmd= Form("%s <- cov(%s[which(attributes(clusterResults)$cluster==%d),])",covMatrixName.c_str(),RDataName.c_str(),i+1);
		DEBUG_LOG("unning cmd: "<<RCmd);

		//Exec cmd
		try{
			Utils::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			ERROR_LOG("Failed to compute covariance matrix for data belonging to cluster no. "<<i+1<<"!");
			return -1;				
		}

		//Fill TMatrixD
		TMatrixD* clusterCovMatrix= Utils::ConvertRTableToROOTMatrix(covMatrixName);
		if(!clusterCovMatrix){
			ERROR_LOG("Failed to retrieve cluster cov matrix and convert it to ROOT!");
			return -1;
		}		
		fClusterCovMatrixes.push_back(clusterCovMatrix);
	}//end loop components

	//## Compute sample skewness for each clusters
	INFO_LOG("Computing sample skewness of each clusters...");
	for(int i=0;i<nComponents;i++){
		//Define cmd
		std::string skewnessVectorName= Form("skewnessVec_cluster%d",i+1);
		RCmd= Form("%s <- skewness(%s[which(clusterResults$cluster==%d),])",skewnessVectorName.c_str(),RDataName.c_str(),i+1);
		DEBUG_LOG("Running cmd: "<<RCmd);

		//Exec cmd
		try{
			Utils::fR.parseEval(RCmd.c_str());
		}
		catch(...){
			ERROR_LOG("Failed to compute skewness for data belonging to cluster no. "<<i+1<<"!");
			return -1;				
		}

		//Retrieve skewness vector
		fClusterSkewness= Utils::ConvertRVectToROOTMatrix(skewnessVectorName);
		if(!fClusterSkewness){
			ERROR_LOG("Failed to retrieve cluster skewness vector and convert it to ROOT!");
			return -1;
		}
	}//end loop components

	return 0;

}//close RunKMedians()


int KMeansClustering::RunKMedians(TMatrixD* dataMatrix,int nComponents,int maxIter,int nRandomInit)
{
	//Check input data
	if(!dataMatrix){
		ERROR_LOG("Null ptr to data matrix given!");
		return -1;
	}
	if(nComponents<=0){
		ERROR_LOG("Invalid numbr of components arg given (must be >0)");
		return -1;
	}

	//## Import data matrix in R
	if(Utils::ImportMatrixInR(dataMatrix,"dataMatrix")<0){
		ERROR_LOG("Failed to import data matrix in R!");
		return -1;
	}

	//## Run kmeans on imported data
	return RunKMedians("dataMatrix",nComponents,maxIter,nRandomInit);

}//close RunKMedians()



}//close namespace

