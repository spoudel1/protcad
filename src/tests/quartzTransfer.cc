#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>

UInt countHbonds(protein* _prot);
double getEnergy(protein* _prot);
UInt countCAOclashes(protein* _prot);
string  intToString(int _num);
void saveState(protein* _prot, string _prefix, UInt _fileNum, UInt _step, double _beta, double _score);
bool isFixed(UIntVec _fixedPositions, UInt _res);
bool flipChirality(protein* _prot, UInt _res);

double getProb(double _oldScore, double _score, double _beta);
vector <dblVec> getLigandCoords(string _fileName);

int main (int argc, char* argv[])
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};


	//read command line
	if (argc < 2)
	{
		cout << "primalSoup inputfile" << endl;
		exit(1);
	}

	//read inputfile

	string inputFile = argv[1];
	ifstream inFile;
	inFile.open(inputFile.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}
	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	string inputFileName = parsedStrings[0];

	// read in prot structure
	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);


	residue::setCutoffDistance(4.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.9);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	ran ranNumber;
	UInt randomSeed;

	if (argc == 3)
	{
		getline(inFile, currentLine, '\n'); // read in but skip
		sscanf(argv[2], "%u", &randomSeed); // use second command line arg instead
		ranNumber.setSeed(randomSeed);
	}
	else
	{
		getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[0].c_str(), "%u", &randomSeed);
		ranNumber.setSeed(randomSeed);
	}
	UInt iterations, history;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &iterations); sscanf(parsedStrings[1].c_str(), "%u", &history);

	double startBeta, endBeta;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &startBeta); sscanf(parsedStrings[1].c_str(), "%lf", &endBeta);

	double probLAla; // enantiomeric excess of L-Ala
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &probLAla);





	UInt size = prot->getNumResidues(0);
	for (UInt i = 0; i < size; i ++)
	{
		prot->activateForRepacking(0,i);
		prot->setPhi(0,i,180);
		prot->setPsi(0,i,180);
		if (ranNumber.getNext() > 0.5)  // randomize initial sequence
			prot->mutate(0, i, 20);
		else
			prot->mutate(0, i, 0);
	}

	UIntVec fixedPositions;
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings=Parse::parse(currentLine);
		UInt tempRes, tempID;
		sscanf(parsedStrings[0].c_str(), "%u", &tempRes);
		sscanf(parsedStrings[1].c_str(), "%u", &tempID);
		fixedPositions.push_back(tempRes);
		prot->mutate(0, tempRes, tempID);
	}
	
	// center dimer on origin
	dblVec center = prot->getBackBoneCentroid(0);
	prot->translate(0, -1.0* center);
	center = prot->getBackBoneCentroid(1);
	prot->translate(1, -1.0* center);
	// move mobile monomer away

	double xtrans, ytrans, ztrans, irot1, irot2, irot3;
	xtrans = 3 + 6 * ranNumber.getNext();
	ytrans = 3 + 6 * ranNumber.getNext();
	ztrans = 3 + 6 * ranNumber.getNext();
	irot1 = -180 + 360*ranNumber.getNext();
	irot2 = -180 + 360*ranNumber.getNext();
	irot3 = -180 + 360*ranNumber.getNext();

	prot->eulerRotate(0, irot1, irot2, irot3);
	prot->translate(0, xtrans, ytrans, ztrans);
	//doMove(prot, 0, radius, rot1, rot2, rot3);

	pdbWriter(prot, "start.pdb");


	prot->silenceMessages();
	double lowestScore = 1000;
	double oldScore = lowestScore;
	UInt lowCounter = 1;

	UInt counter = 0;
	for (UInt i = 0; i < iterations; i ++)
	{
		counter ++;   // counter for storing history files
		double beta = startBeta + i * (endBeta - startBeta) / iterations;  // linear cooling

		if ( 0.80 > ranNumber.getNext()) // then deal with intrachain optimization
		{
			UInt res = UInt(ranNumber.getNext()*(size+1));
			if (res >= size) res = size - 1;

			UInt choice = 0;  // 0 is change dihedral  (default)  , 1 is change id, 2 is change both
			if (isFixed(fixedPositions,res))
			{
				choice = 0;
			}
			else
			{
				double tempProb = ranNumber.getNext();
				if (tempProb < probLAla && prot->getTypeFromResNum(0,res) == 20) // currently D, then flip to L
				{
					if (ranNumber.getNext() < 0.5) choice = 1;
					else choice = 2;
				}
				else if (tempProb >= probLAla && prot->getTypeFromResNum(0,res) == 0) // currently L, then flip to D
				{
					if(ranNumber.getNext() < 0.5) choice = 1;
					else choice = 2;
				}
			}

			// CHOICE 0 - change dihedral only

			if (choice == 0)
			{
				double oldPhi = prot->getPhi(0,res);
				double oldPsi = prot->getPsi(0,res);
				double newPhi = 180 - 360*ranNumber.getNext();
				double newPsi = 180 - 360*ranNumber.getNext();
				prot->setPhi(0,res,newPhi);
				prot->setPsi(0,res,newPsi);

				double score = getEnergy(prot);
				if (score < oldScore)
				{
					oldScore = score;
					if (score < lowestScore)
					{
						lowestScore = score;
						lowCounter ++;
						saveState(prot, "low", lowCounter, i, beta, score);
						pdbWriter(prot, "final.pdb");
					}
				}
				else
				{
					double metropolis = getProb(oldScore, score, beta);
					if ( metropolis > ranNumber.getNext() )
					{
						oldScore = score;
					}
					else
					{
						prot->setPhi(0, res, oldPhi);
						prot->setPsi(0, res, oldPsi);
					}
				}
			}

			// CHOICE 1 - flip chirality

			if (choice == 1)
			{
				flipChirality(prot,res);
				double score = getEnergy(prot);
				if (score < oldScore)
				{
					oldScore = score;
					if (score < lowestScore)
					{
						lowestScore = score;
						lowCounter ++;
						saveState(prot, "low", lowCounter, i, beta, score);
						pdbWriter(prot, "final.pdb");
					}
				}
				else
				{
					double metropolis = getProb(oldScore, score, beta);
					if ( metropolis > ranNumber.getNext() )
					{
						oldScore = score;
					}
					else
					{
						flipChirality(prot, res); // if rejected, flip back
					}
				}
			}

			// CHOICE 2 - change both simultaneously

			if (choice == 2)
			{
				flipChirality(prot, res);
				double oldPhi = prot->getPhi(0,res);
				double oldPsi = prot->getPsi(0,res);
				double newPhi = 180 - 360*ranNumber.getNext();
				double newPsi = 180 - 360*ranNumber.getNext();
				prot->setPhi(0,res,newPhi);
				prot->setPsi(0,res,newPsi);

				double score = getEnergy(prot);
				if (score < oldScore)
				{
					oldScore = score;
					if (score < lowestScore)
					{
						lowestScore = score;
						lowCounter ++;
						saveState(prot, "low", lowCounter, i, beta, score);
						pdbWriter(prot, "final.pdb");
					}
				}
				else
				{
					double metropolis = getProb(oldScore, score, beta);
					if ( metropolis > ranNumber.getNext() )
					{
						oldScore = score;
					}
					else
					{
						flipChirality(prot, res);
						prot->setPhi(0, res, oldPhi);
						prot->setPsi(0, res, oldPsi);
					}
				}
			}
		}
		else  // translate/rotate
		{
			dblVec centroid = prot->getBackBoneCentroid(0);
			dblVec newPos(3);
			newPos[0] = centroid[0] - 1 + 2.0 * ranNumber.getNext();
			newPos[1] = centroid[1] - 1 + 2.0 * ranNumber.getNext();
			newPos[2] = centroid[2] - 1 + 2.0 * ranNumber.getNext();

			double radius = sqrt(CMath::dotProduct(newPos,newPos));
			if (radius > 12.0)
			{
				newPos = newPos * 12.0/radius;
			}
			if (radius < 3.0)
			{
				newPos = newPos * 3.0/radius;
			}

			double rot1 = - 15 + 30.0 * ranNumber.getNext();
			double rot2 = - 15 + 30.0 * ranNumber.getNext();
			double rot3 = - 15 + 30.0 * ranNumber.getNext();

			prot->translate(0, -1.0*centroid);
			prot->eulerRotate(0, rot1, rot2, rot3);
			prot->translate(0, newPos);
			double score = getEnergy(prot);
			if (score < oldScore)
			{
				oldScore = score;
				if (score < lowestScore)
				{
					lowestScore = score;
					lowCounter ++;
					saveState(prot, "low", lowCounter, i, beta, score);
					pdbWriter(prot, "final.pdb");
				}
			}
			else
			{
				double metropolis = getProb(oldScore, score, beta);
				if ( metropolis > ranNumber.getNext() )
				{
					oldScore = score;
				}
				else
				{
					prot->translate(0, -1.0* newPos);
					prot->undoEulerRotate(0, rot1, rot2, rot3);
					prot->translate(0, centroid);
				}
			}
		}
		cout << "Score " << i << " " << oldScore << " temp " << beta << " lowest " << lowestScore << endl;
		if (counter == history)
		{
			counter = 0;
			saveState(prot, "out", i, i, beta, oldScore);
		}
	}

	pdbWriter(prot,"end.pdb");
	return 0;
}

double getProb(double _oldScore, double _score, double _beta)
{
	double dScore = _oldScore - _score;
	double prob = pow(2.718,  dScore * _beta);
	return prob;
}

double getEnergy(protein* _prot)
{

	UInt nhbs = countHbonds(_prot);
	UInt nclashes = countCAOclashes(_prot);
	double vdwE = _prot->intraEnergy();

	dblVec centroid = _prot->getBackBoneCentroid(0);
	double radius = sqrt(CMath::dotProduct(centroid, centroid));

	double force = -1 * (3 - radius);

	double score = force + vdwE + 5.0*(float)nclashes - 5.0*(float)nhbs;
	return score;
}

UInt countCAOclashes(protein* _prot)
{
	UInt size = _prot->getNumResidues(0);
	UInt clashes = 0;
	vector < dblVec > CBcoords;
	vector < dblVec > Ocoords;

	for (UInt i = 0; i < size-1; i ++) // do not count clashes on last residue as we cannot change its PSI anyways
	{
		dblVec tempcoord;
		UInt type = _prot->getTypeFromResNum(0,i);
		if (type == 0 || type == 20)
		{
			tempcoord =_prot->getCoords(0, i, 4);
			CBcoords.push_back(tempcoord);
		}
		tempcoord =_prot->getCoords(0, i, 3);
		Ocoords.push_back(tempcoord);
	}

	for (UInt i = 0; i < CBcoords.size(); i ++)
	{
		for (UInt j = 0; j < Ocoords.size(); j ++)
		{
			dblVec temp = CBcoords[i] - Ocoords[j];
			double distance = sqrt(CMath::dotProduct(temp,temp));
			int resSpace = (i - j);
			if (distance < 3.0 && abs(resSpace) < 2 )
			{
				clashes ++;
			}
		}
	}
	return clashes;
}



UInt countHbonds(protein* _prot)
{
	UInt size = _prot->getNumResidues(0);
	UInt numHbonds = 0;
	vector < dblVec > Ncoords;
	vector < dblVec > Ocoords;
	vector < dblVec > Ccoords;
	vector < dblVec > CAcoords;

	for (UInt i = 0; i < size; i ++)
	{
		dblVec tempcoord =_prot->getCoords(0, i, 0);
		Ncoords.push_back(tempcoord);
		tempcoord =_prot->getCoords(0, i, 3);
		Ocoords.push_back(tempcoord);
		tempcoord =_prot->getCoords(0, i, 1);
		CAcoords.push_back(tempcoord);
		tempcoord = _prot->getCoords(0, i, 2);
		Ccoords.push_back(tempcoord);
	}

	// for first residue = no angular constraint
	for (UInt j = 0; j < Ocoords.size(); j ++)
	{
		dblVec temp = Ncoords[0] - Ocoords[j];
		double distance = sqrt(CMath::dotProduct(temp,temp));
		int resSpace = (0 - j);
		if (distance < 3.2 && abs(resSpace) > 2 )
		{
			numHbonds ++;
		}
	}


	for (UInt i = 1; i < Ncoords.size(); i ++)
	{
		for (UInt j = 0; j < Ocoords.size(); j ++)
		{
			int resSpace = (i - j);
			if ( abs(resSpace) > 2)
			{
				dblVec NO = Ncoords[i] - Ocoords[j];
				double distance = sqrt(CMath::dotProduct(NO,NO));

				dblVec pseudoAtom = (Ccoords[i-1] + CAcoords[i])/2.0;
				dblVec NH = Ncoords[i] - pseudoAtom;

				double magNH = sqrt(CMath::dotProduct(NH,NH));
				double angle = acos( CMath::dotProduct(NO,NH) / (magNH * distance) );
				angle = angle * 180 / 3.14159;
				if (distance < 3.2 && angle > 90.0)
				{
					numHbonds ++;
				}
			}
		}
	}
	return numHbonds;
}

string intToString(int _num)
{
	ostringstream myStream;
	myStream << _num << flush;
	return (myStream.str());
}

void saveState(protein* _prot, string _prefix, UInt _fileNum, UInt _step, double _beta, double _score)
{

	string pad = "";
	if     (_fileNum < 10) pad = "000000";
	else if(_fileNum < 100) pad = "00000";
	else if(_fileNum < 1000) pad = "0000";
	else if(_fileNum < 10000) pad = "000";
	else if(_fileNum < 100000) pad = "00";
	else if(_fileNum < 1000000) pad = "0";


	string filename = _prefix + pad + intToString(_fileNum) + ".pdb";
	pdbWriter(_prot, filename);
	string logfile = _prefix + ".log";
	ofstream fout (logfile.c_str() ,ios::app);
	fout << "step " << _step << " beta " << _beta << " score " << _score << " nhbs " <<
		countHbonds(_prot)<< " clashes " <<  countCAOclashes(_prot) << " sequence ";
	UInt size = _prot->getNumResidues(0);
	for (UInt n = 0; n < size; n ++)
	{
		if (_prot->getTypeFromResNum(0,n) == 0)
		{
			fout << "L ";
		}
		else if (_prot->getTypeFromResNum(0, n) == 20)
		{
			fout << "D ";
		}
		else
		{
			fout << _prot->getTypeStringFromResNum(0,n) << " ";
		}
	}
	for (UInt n = 0; n < size; n ++)
	{
		fout << _prot->getPhi(0,n) << " ";
		fout << _prot->getPsi(0,n) << " ; ";
	}
	fout << endl;

	fout.close();
	return;
}

bool flipChirality(protein* _prot, UInt _res)
{
	if (_prot->getTypeFromResNum(0,_res) == 0)
	{
		_prot->mutate(0, _res, 20);
		return true;
	}
	else if (_prot->getTypeFromResNum(0, _res) == 20)
	{
		_prot->mutate(0, _res, 0);
		return true;
	}
	return false;
}

bool isFixed(UIntVec _fixedPositions, UInt _res)
{
	bool flag = false;
	for (UInt i = 0; i < _fixedPositions.size(); i ++)
	{
		if (_res == _fixedPositions[i])
		{
			flag = true;
		}
	}

	return flag;
}

vector <dblVec> getLigandCoords(string _inputFileName)
{
    ifstream inFile;
    inFile.open(_inputFileName.c_str());
    if (!inFile)
    {
        cout << "Unable to find or open file" << endl;
        exit(1);
    }

    string currentLine;
    vector <string> parsedStrings;
    parsedStrings.resize(0);
    vector <dblVec> ligCoords;
    while(getline(inFile, currentLine, '\n'))
    {
        parsedStrings = parse(currentLine);
        if (parsedStrings[0] == "ATOM")
        {
            double x,y,z;
            sscanf(parsedStrings[5].c_str(), "%lf", &x);
            sscanf(parsedStrings[6].c_str(), "%lf", &y);
            sscanf(parsedStrings[7].c_str(), "%lf", &z);
            dblVec thisAtom(3);
            thisAtom[0] = x; thisAtom[1] = y ; thisAtom[2] = z;
            cout << "xyz " << x << " " << y << " " << z << " " << endl;
            ligCoords.push_back(thisAtom);
        }
    }
    inFile.close();
    return ligCoords;
}

double getProteinLigandEnergy(protein* _prot, vector <dblVec> _ligCoords, double _r0, double _eps, vector <UIntVec> _exclude)
{
    double energy = 0.0;
    for (UInt i = 0; i < _prot->getNumChains(); i ++)
    {
        for (UInt j = 0; j < _prot->getNumResidues(i); j ++)
        {
            bool excludeRes = false;
            for (UInt n = 0; n < _exclude.size(); n ++)
            {
                if (_exclude[n][0] == i && _exclude[n][1] == j)
                {
                    excludeRes = true;
                }
            }
            if (!excludeRes)
            {
                for (UInt k = 0; k < _prot->getNumAtoms(i,j); k ++)
                {
                    dblVec atomCoords = _prot->getCoords(i,j,k);
                    for (UInt l = 0; l < _ligCoords.size(); l ++)
                    {
                        dblVec diff = atomCoords - _ligCoords[l];
                        double distance = sqrt(CMath::dotProduct(diff,diff));
                        double distratio = _r0/distance;
                        energy += _eps*(pow(distratio,12) - 2* pow(distratio,6));
                    }
                }
            }
        }
    }

    return energy;
}

