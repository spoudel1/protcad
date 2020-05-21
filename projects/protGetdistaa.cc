//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protGetdistaa    *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include "protein.h"
#include <time.h>
#include<map>
using namespace std;
string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Saf","Hem"};
//--Program setup-------------------------------------------------------------
//Below is the maximum accessible surface area (ASA) of all 20 types of amino acid residues (likely doesn't include the water radius - need to check this though) calculated using tabulateSurfaceArea
map<string,double>maxsolv={{"A",94.10},{"C",112.74},{"D",133.59},{"E",157.47},{"F",185.48},{"G",94.73},{"H",176.27},{"I",137.03},{"K",178.33},{"L",151.53},{"M",150.08},{"N",137.53},{"P",117.59},{"Q",150.54},{"R",207.13},{"S",106.33},{"T",119.83},{"V",124.97},{"W",212.63},{"Y",216.55}};
//standard charges of all amino acids
map<string,double>chargeaa={{"A",0.00},{"C",0.00},{"D",-1.00},{"E",-1.00},{"F",0.00},{"G",0.00},{"H",0.00},{"I",0.00},{"K",1.00},{"L",0.00},{"M",0.00},{"N",0.00},{"P",0.00},{"Q",0.00},{"R",1.00},{"S",0.00},{"T",0.00},{"V",0.00},{"W",0.00},{"Y",0.00}};
//get median dielectric of all the residues
double getmedianDielectric(protein* &_prot)
{
	_prot->updateDielectrics();
	double median, resD;
	vector <double> resDielectric;
	for (UInt i = 0; i < _prot->getNumChains(); i++)
	{
		for (UInt j = 0; j < _prot->getNumResidues(i); j++)
		{
			resD = _prot->getDielectric(i,j);
			resDielectric.push_back(resD);
		}
	}
	
	size_t size = resDielectric.size();
	sort(resDielectric.begin(), resDielectric.end());
	if (size % 2 == 0)
	{
		median = (resDielectric[size / 2 - 1] + resDielectric[size / 2]) / 2;
	}
	else
	{
		median = resDielectric[size / 2];
	}
	return median;
}
//get standard dielectric of all the residues
double getstandarddev(protein* &_prot)
{
	_prot->updateDielectrics();
	double stddev=0.0, resD, dielectric_sum=0.0;
	UInt count=0;
	vector <double> resDielectric;
	for (UInt i = 0; i < _prot->getNumChains(); i++)
	{
		for (UInt j = 0; j < _prot->getNumResidues(i); j++)
		{
			resD = _prot->getDielectric(i,j);
			dielectric_sum+=resD;
			count++;
			resDielectric.push_back(resD);
		}
	}
	double meanDielectric = resD/count;
	for (int k=0; k< resDielectric.size();k++){
		stddev+=pow(resDielectric[k]-meanDielectric,2);
	}
	return (meanDielectric*meanDielectric) - sqrt(stddev/resDielectric.size());
}
//get average dielectric of all the residues
double getavdielectric(protein* &_prot){
	double aveDielectric, sumDielectric=0.0;
	UInt counter = 0;
	_prot->updateDielectrics();
	for (UInt i = 0; i < _prot->getNumChains(); i++)
	{
		for (UInt j = 0; j < _prot->getNumResidues(i); j++)
		{
			sumDielectric += _prot->getDielectric(i,j);
			counter++;
		}
	}
	aveDielectric = sumDielectric/counter;
	return aveDielectric;
}
//this function will check for the 30% accessible surface area (ASA) threshold compared to the maximum ASA
int checkthreshold(string aa, double sufarea){
	map<string,double>::iterator itr;
	int present=0;
    for (itr=maxsolv.begin();itr!=maxsolv.end();itr++){
		if (aa==itr->first && sufarea>(0.05*itr->second)){
			present=1;
			break;
		}
    }
	return present;
}
int main (int argc, char* argv[])
{	

	if (argc !=5)
	{
		cout << "protGetdistaa <inFile.pdb> centralatomtype maxdistance mindistance" << endl;
		exit(1);
	}
	string infile = argv[1];
//	UInt centralatom=atoi(argv[2]);
	string centralatom=argv[2];
	UInt maxdistance = atoi(argv[3]);
	UInt mindistance = atoi(argv[4]);
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	atomIterator theIter1(bundle);
	atom* pAtom;

	bundle->updateDielectrics();
	dblVec central_coords, CA_coords;
	UInt chainNum = bundle->getNumChains();
	double averageDielectric=getavdielectric(bundle);
	double medianDielectric=getmedianDielectric(bundle);
	double stddevDielectric=getstandarddev(bundle);

	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{
			if (bundle->getTypeStringFromResNum(i,j)==centralatom){
			//if (bundle->getResNum(i,j)==centralatom){
			//	central_coords = bundle->getCoords(i, j, "CA");
				central_coords = bundle->getCoords(i, j, 0);
				break;
			}
		}
	}
	vector<string>selectaa; //this vector will store all the residues within specific angstrom from the central atom
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{
	//		bundle->initializeSpherePoints(i,j);
	//		bundle->removeSpherePoints(i,j);
	//		cout<<aminoAcidString[bundle->getTypeFromResNum(i,j)]<<" "<<bundle->tabulateSurfaceArea(i,j)<<" exposed surface"<<endl;
	//		checkthreshold(aminoAcidString[bundle->getTypeFromResNum(i,j)],bundle->tabulateSurfaceArea(i,j));
			if (!bundle->isCofactor(i,j)){
				UInt atomNum = bundle->getNumAtoms(i,j);
				double fdist=CMath::distance(central_coords, bundle->getCoords(i, j, 0)); //distance from the central atom
				UInt final_i=i;
				UInt final_j=j;
				for (UInt k = 1; k < atomNum; k ++)
				{
					CA_coords = bundle->getCoords(i, j, k);
					double dist = CMath::distance(central_coords, CA_coords);
					if (dist<fdist){
						fdist=dist;
						final_i=i;
						final_j=j;
					}
				}
				if (fdist<=maxdistance && fdist >=mindistance){
					bundle->initializeSpherePoints(final_i,final_j); //initialize the sphere
					bundle->removeSpherePoints(final_i,final_j); //remove any overlapping points in the sphere
					if (checkthreshold(aminoAcidString[bundle->getTypeFromResNum(final_i,final_j)],bundle->tabulateSurfaceArea(final_i,final_j))==1){
						selectaa.push_back(aminoAcidString[bundle->getTypeFromResNum(i,j)]); //uncomment it if using print option 1 (i.e., single run)
						cout<<aminoAcidString[bundle->getTypeFromResNum(i,j)]<<" "<<bundle->getResNum(i,j)<<" "<<fdist<<endl;
						//	string tpair=(aminoAcidString[bundle->getTypeFromResNum(i,j)])+"#"+to_string(dist);
						//	cout<<tpair<<endl;
	//					selectaa.push_back((aminoAcidString[bundle->getTypeFromResNum(final_i,final_j)])+"#"+to_string(fdist));//only use this if using the second option to print (i.e, multiple angs...)
	//					cout<<(aminoAcidString[bundle->getTypeFromResNum(final_i,final_j)])+"#"+to_string(fdist)<<" "<<bundle->getResNum(final_i,final_j)<<endl;
					}
				}
			}
		}
	}
	
	//Print options:
	//1. for single run of specific distance
	
	int netcharge=0;
	if (selectaa.size()>=5){
		for (int d=0;d<selectaa.size();d++){
			map<string,double>::iterator itr1;
    		for (itr1=chargeaa.begin();itr1!=chargeaa.end();itr1++){
				if (selectaa[d]==itr1->first){
					netcharge+=itr1->second;
				}
  	  		}
		}
		cout<<argv[1]<<" "<<mindistance<<"-"<<maxdistance<<" Total number of residues: "<<selectaa.size()<<" net charge: "<<netcharge<<" and average charge: "<<(netcharge/double(selectaa.size()))<<endl;
	}else{
		cout<<argv[1]<<" number of residues less than the specified threshold (i.e., 5 residues minimum)"<<endl;
	}
	/*
	//2. for multiple angstrom use this
	for (int e=mindistance;e<maxdistance;e++){
		for (int f=1;f<=(maxdistance-mindistance);f++){
			int aacount=0;
			int netcharge=0;
			if ((e+f)<=maxdistance){
				for (int d=0;d<selectaa.size();d++){
					string tdist=selectaa[d].substr(selectaa[d].find('#')+1);
					string tempaa=selectaa[d].substr(0,selectaa[d].find('#'));
					if (stod(tdist)>=e && stod(tdist)<=(e+f)){
						map<string,double>::iterator itr1;
    					for (itr1=chargeaa.begin();itr1!=chargeaa.end();itr1++){
							if (tempaa==itr1->first){
								netcharge+=itr1->second;
								aacount+=1;
							}
  	 	 				}
					}
				}
				if (aacount>=5){
					cout<<argv[1]<<" "<<e<<"-"<<e+f<<" Total number of residues: "<<aacount<<" net charge: "<<netcharge<<" and average charge: "<<(netcharge/double(aacount))<<endl;
				}else{
					cout<<argv[1]<<" "<<e<<"-"<<e+f<<" Total number of residues: "<<aacount<<" net charge: "<<"NA"<<" and average charge: "<<"NA"<<endl;
				}
			}
		}
	}*/


}
