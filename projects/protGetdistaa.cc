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
string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Saf","Hem"};
//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{	

	if (argc !=4)
	{
		cout << "protGetdistaa <inFile.pdb> maxdistance mindistance" << endl;
		exit(1);
	}
	string infile = argv[1];
	UInt maxdistance = atoi(argv[2]);
	UInt mindistance = atoi(argv[3]);
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	dblVec MG_coords, CA_coords;
	UInt chainNum = bundle->getNumChains();
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{
			if (bundle->getTypeStringFromResNum(i,j)=="MG2"){
				MG_coords = bundle->getCoords(i, j, "MG");
			}
		}
	}
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{
			if (!bundle->isNotAminoAcid(i,j)){
				UInt atomNum = bundle->getNumAtoms(i,j);
				for (UInt k = 0; k < atomNum; k ++)
				{
					CA_coords = bundle->getCoords(i, j, k);
					double dist = CMath::distance(MG_coords, CA_coords);
					if (dist<maxdistance && dist >mindistance){
						cout<<aminoAcidString[bundle->getTypeFromResNum(i,j)]<<" "<<bundle->getResNum(i,j)<<endl;
						break;
					}
				}
			}
		}
	}
}
//	dblVec Ocoords = frame->getCoords(0, 90, "HH");
  //      dblVec Ncoords = frame->getCoords(0, 14, "NE2");
    //    dist = CMath::distance(Ocoords, Ncoords);
     //   cout << a << " " << (a*.02) << " " << dist << endl;
