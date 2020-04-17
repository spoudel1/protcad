/**********************************************************************
  *********************************************************************
 	filename: residue_basic.cpp
 	contents: class residue implemention
 	static variables defined: howMany
 	functions defined:
				constructors
				destructors
				generateAtoms(..)
				buildConnectivity(..)
				addAtom(..)
				queryChildren(..)
				setHydrogensOn(..)
  *********************************************************************
  *********************************************************************/

#include "residue.h"
typedef vector<atom*>::iterator iterATOM;
vector<residueTemplate> residue::dataBase;
bool residue::dataBaseBuilt = false;
bool residue::polarizableElec = true;
double residue::temperature = 300.0;
double residue::HsolvationFactor = 1.0;
double residue::EsolvationFactor = 1.0;
double residue::cutoffDistance = 7.0;
double residue::cutoffDistanceSquared = residue::cutoffDistance*residue::cutoffDistance;
double residue::cutoffCubeVolume = pow((residue::cutoffDistance*2),3);
double residue::dielectricWidth = 2.8; // diameter of water
double residue::KT = KB*residue::getTemperature();

void residue::setupDataBase()
{	if (!dataBaseBuilt)
	{	residue* dummyRes = new residue(1);
		delete dummyRes;
	}
	return;
}

void residue::setupDataBase(const bool _Hflag)
{	if (!dataBaseBuilt)
	{	residue* dummyRes = new residue(1,_Hflag);
		delete dummyRes;
	}
	return;
}

void residue::setupDataBase(const bool _Hflag, const bool _HPflag)
{	if (!dataBaseBuilt)
	{	residue* dummyRes = new residue(1,_Hflag,_HPflag);
		delete dummyRes;
	}
	return;
}

UInt residue::howMany = 0;


// Constructors and Utilities


residue::residue()
{
	hydrogensOn = false;
	polarHydrogensOn= false;
	if(!dataBaseBuilt)
	{
		//cout << "default build: ";
		buildDataBase();
	}
	itsAtoms.resize(0);
	itsSidechainDihedralAngles.resize(0);
	howMany ++;
}

residue::residue(const UInt _itsType)
{
	hydrogensOn = false;
	polarHydrogensOn = false;
	if(!dataBaseBuilt)
	{
		//cout << "itsType build:  ";
		buildDataBase();
	}
	itsType = _itsType;
	// hydrogensOn defaults to false if not explicity set
	initializeAtomsAndConnectivity();
	initializeSidechainDihedralAngles();
}

residue::residue(const UInt _itsType, const bool _hFlag)
{
	hydrogensOn = _hFlag;
	polarHydrogensOn = false;
	if(!dataBaseBuilt)
	{
		//cout << "itstype and hflag build:  ";
		buildDataBase();
	}
	itsType = _itsType;
	initializeAtomsAndConnectivity();
	initializeSidechainDihedralAngles();
	if(hydrogensOn)
	{
		initializePolarHDihedralAngle();
	}
}

residue::residue(const UInt _itsType, const bool _hFlag, const bool _hPFlag)
{
	hydrogensOn = _hFlag;
	polarHydrogensOn = _hPFlag;

	if(hydrogensOn && polarHydrogensOn)
	{ cout << "error!! only one hFlag can be true.\n";
	  cout << "Setting only polarHydrogensOn as true.\n";
	  polarHydrogensOn = true;
	  hydrogensOn = false;
	}

	if(!dataBaseBuilt)
	{
		if(polarHydrogensOn)
		{
			cout << "itstype and hPflag build:  \n";
		}
		if(hydrogensOn)
		{
			cout << "itstype and Hflag build:  \n";
		}
		buildDataBase();
	}
	itsType = _itsType;
	initializeAtomsAndConnectivity();
	initializeSidechainDihedralAngles();
	if(hydrogensOn || polarHydrogensOn)
	{
		initializePolarHDihedralAngle();
	}
}

residue::residue(const string& _aaType)
{
	hydrogensOn = true;
	polarHydrogensOn = false;
	if(!dataBaseBuilt)
	{
		cout << "string aaType build:  ";
		buildDataBase();
	}
	// hydrogensOn defaults to true if not explicity set
	defineType(_aaType);
	residue(itsType, hydrogensOn);
}

residue::residue(const string& _aaType, const bool _hFlag)
{
	hydrogensOn = _hFlag;
	polarHydrogensOn = false;
	if(!dataBaseBuilt)
	{
		cout << "string aaType and hflag build:  ";
		buildDataBase();
	}
	defineType(_aaType);
	residue(itsType,_hFlag);
}

residue::residue(const string& _aaType, const bool _hFlag, const bool _hPFlag)
{
	hydrogensOn = _hFlag;
	polarHydrogensOn = _hPFlag;

	if(hydrogensOn && polarHydrogensOn)
	{
		cout << "Error!  Only hydrogensOn or polarHydrogensOn can be true, not both.\n";
		cout << "Setting polarHydrogensOn = true.\n";
		polarHydrogensOn = true;
		hydrogensOn = false;
	}

	if(!dataBaseBuilt)
	{
		if(hydrogensOn)
		{ 	cout << "string aaType and hflag build:  \n"; }
		if(polarHydrogensOn)
		{	cout << "string aaType and HPflag build:  \n"; }
		buildDataBase();
	}
	defineType(_aaType);
	residue(itsType,_hFlag,_hPFlag);
}

residue::residue(const residue& _rhs)
{
	hydrogensOn = _rhs.hydrogensOn;
	polarHydrogensOn = _rhs.polarHydrogensOn;
	if(!dataBaseBuilt)
	{
		//cout << "deepcopy build:  ";
		buildDataBase();
	}
	itsType = _rhs.itsType;
	isArtificiallyBuilt = _rhs.isArtificiallyBuilt;
	for(UInt i=0; i<_rhs.itsAtoms.size(); i++)
	{	itsAtoms.push_back(new atom(*(_rhs.itsAtoms[i])));
	}
	itsResNum = _rhs.itsResNum;
	//buildConnectivity();
	buildConnectivityNew();
	itsSidechainDihedralAngles = _rhs.itsSidechainDihedralAngles;
	if(hydrogensOn || polarHydrogensOn)
	{ itsPolarHDihedralAngle = _rhs.itsPolarHDihedralAngle; }
	howMany++;
}

void residue::initializeAtomsAndConnectivity()
{
	isArtificiallyBuilt = false;
	generateAtoms();
	//buildConnectivity();
	buildConnectivityNew();
#ifdef __RES_DEBUG
	for(UInt i=0; i<itsAtoms.size(); i++)
	{	cout << itsAtoms[i]->getName() << " connect to: ";
		itsAtoms[i]->queryChildren();
		cout << endl;
	}
#endif
	howMany ++;
}

// cannot be called before residueTemplate is initialized
void residue::initializeSidechainDihedralAngles()
{
	//cout << getType() << " ";
	DouVec chis;
	UInt numbranchpoints = dataBase[itsType].chiDefinitions.size();
	//cout << "numbranchpoints = " << numbranchpoints << endl;
	for(UInt bpt=0; bpt<numbranchpoints; bpt++)
	{
	//	cout << "branchpoint " << bpt << " ";
		chis.resize(0);
	//	cout << "numDihdralAngles = " << getNumDihedralAngles(itsType,bpt) << endl;
		for(UInt angle=0; angle<getNumDihedralAngles(itsType,bpt); angle++)
		{
	//		cout << "angle " << angle << endl;
			chis.push_back(getChi(bpt, angle));
		}
		itsSidechainDihedralAngles.push_back(chis);
	}
}

void residue::initializePolarHDihedralAngle()
{
	if(dataBase[itsType].getHasPolarHRotamers())
	{	itsPolarHDihedralAngle = getPolarHChi(); }
	else {itsPolarHDihedralAngle = 999.0; }
}

// ************************************************************************
// ************************************************************************
// 	Destructor
// ************************************************************************
// ************************************************************************

residue::~residue()
{
#ifdef __RES_DEBUG
//	cout << "Residue destructor called " << endl;
#endif
	for (UInt i=0; i<itsAtoms.size(); i++)
	{	delete itsAtoms[i];
	}
	howMany--;
}

// ************************************************************************
// ************************************************************************
// 	Generate the Atoms
// ************************************************************************
// ************************************************************************

void residue::generateAtoms()
{
	//cout << dataBase[itsType].atomList.size() << endl;
	for(UInt i=0; i < dataBase[itsType].atomList.size(); i++)
	{
		itsAtoms.push_back(new atom(dataBase[itsType].atomList[i]));
	}
}

// ************************************************************************
// ************************************************************************
// 	Build the Connvectivity
// ************************************************************************
// ************************************************************************

// NOTE !!! that the atom index in the pdb file and connectivity file
// starts from 1
// HOWEVER, in the vector or array, the indexing starts from 0 !!!
// This is why there is a "-1" being carried around in the code below.
/*
void residue::buildConnectivity()
{	UInt tempInt;
	vector<UInt> intVector;
	intVector.resize(0);
	cout << "itsType = " << itsType << endl;
	for (UInt i=0; i<dataBase[itsType].connectivity.size(); i++)
	{	cout << dataBase[itsType].connectivity[i] << " ";
	}
	cout << endl;
	for(UInt i=0; i< dataBase[itsType].connectivity.size(); i++)
	{	tempInt = dataBase[itsType].connectivity[i];
		if (tempInt == 0)
		{	// we're starting a new entry
			if (intVector.size() !=0)
			{	 // add the first child
				if (intVector.size() > 1)
				{	if (intVector[0]-1 < itsAtoms.size())
					{	if(intVector[1]-1 < itsAtoms.size())
						{
							if (!itsAtoms[intVector[1]-1]->getParent() &&
							    !itsAtoms[intVector[1]-1]->getPreviousSib() &&
							    !itsAtoms[intVector[1]-1]->isHeadNode())
							{
							 itsAtoms[intVector[0]-1]->setChild
							(itsAtoms[intVector[1]-1]);
							}
						}
						else
						{	intVector.resize(0);
							continue;
						}
					}
				}
				// recursively add the next child as the sibling to
				// the previous child
				for(UInt j=2; j<intVector.size(); j++)
				{	if(intVector[j-1]-1 < itsAtoms.size())
					{	if(intVector[j]-1 < itsAtoms.size())
						{
							if (!itsAtoms[intVector[j]-1]->getParent() &&
							    !itsAtoms[intVector[j]-1]->getPreviousSib() &&
							    !itsAtoms[intVector[j]-1]->isHeadNode())
							{
							 itsAtoms[intVector[j-1]-1]->setNextSib
							(itsAtoms[intVector[j]-1]);
							}
						}
						else
						{	intVector.resize(0);
							continue;
						}
					}
				}
				intVector.resize(0);
			}
		}
		else
		{	intVector.push_back(tempInt);
		}
	}
}
*/

void residue::buildConnectivityNew()
{	vector<vector<UInt> > connectivityVector;
	vector<UInt> perAtomConnectivityVector;
	// transform linear connectivity vector into a more usable format
	for (UInt i=0; i<dataBase[itsType].connectivity.size(); i++)
	{
		UInt atomNumber = dataBase[itsType].connectivity[i];
		if (atomNumber == 0)
		{
			if (i != 0 )
			{	connectivityVector.push_back(perAtomConnectivityVector);
			}
			perAtomConnectivityVector.resize(0);
			continue;
		}
		perAtomConnectivityVector.push_back(atomNumber);
	}
	connectivityVector.push_back(perAtomConnectivityVector);


#ifdef CONNECTIVITY_DEBUG
	cout << "BUILD CONNECTIVITY NEW" << endl;
	for (UInt i=0; i<connectivityVector.size(); i++)
	{	for (UInt j=0; j<connectivityVector[i].size(); j++)
		{	cout << connectivityVector[i][j] << " ";
		}
		cout << endl;
	}
	cout << "END BUILD CONNECTIVITY NEW" << endl;
#endif

	UInt firstChildPlaceholder;
	for (UInt i=0; i<connectivityVector.size(); i++)
	{
		// We only need to do something if the perAtomConnectivityVector is larger than 1 field
		UInt perAtomConnectivityVectorSize = connectivityVector[i].size();
		if ( perAtomConnectivityVectorSize > 1 )
		{
			firstChildPlaceholder = 0;
			// Parent atom occurs at j==0
			// Step through the other atoms to find the first "valid" first child
			for (UInt j=1; j<connectivityVector[i].size(); j++)
			{
				if (connectivityVector[i][j]-1 < itsAtoms.size())
				{
					if (!itsAtoms[connectivityVector[i][j]-1]->getParent() &&
					    !itsAtoms[connectivityVector[i][j]-1]->getPreviousSib() &&
					    !itsAtoms[connectivityVector[i][j]-1]->isHeadNode())
					{
						firstChildPlaceholder = j;
						break;
					}
				}
			}
			// If, at this point, the firstChildPlaceholder is still set to zero,
			// there are no valid atoms, so continue to the top of the FOR loop
			// and deal with the next atom
			if (firstChildPlaceholder == 0)
				continue;
			// Now set the first child
			itsAtoms[connectivityVector[i][0]-1]->setChild(itsAtoms[connectivityVector[i][firstChildPlaceholder]-1]);
			// If there are still more atoms, make them siblings of the first child
			if (firstChildPlaceholder < perAtomConnectivityVectorSize)
			{
				for (UInt j=firstChildPlaceholder+1; j<perAtomConnectivityVectorSize; j++)
				{
					if ( connectivityVector[i][j]-1 < itsAtoms.size())
					{
						if (!itsAtoms[connectivityVector[i][j]-1]->getParent() &&
						    !itsAtoms[connectivityVector[i][j]-1]->getPreviousSib() &&
						    !itsAtoms[connectivityVector[i][j]-1]->isHeadNode())
						{
							itsAtoms[connectivityVector[i][j-1]-1]->setNextSib(itsAtoms[connectivityVector[i][j]-1]);
						}
					}
				}
			}
		}
	}
}


// ************************************************************************
// ************************************************************************
// Atoms Related
// ************************************************************************
// ************************************************************************

void residue::addAtom(pdbAtom& _theData)
{	for(UInt i=0; i < getAtomNameBaseSize(itsType); i++)
	{	if(_theData.getItem(pdbAtomName) == getAtomNameBaseItem(itsType,i))
		{	itsAtoms[i]->setCoords(_theData.getAtomCoord());
			itsAtoms[i]->setSerialNumber(_theData.getSerial());
			itsAtoms[i]->setChainID(_theData.getChainID());
			itsAtoms[i]->setOccupancy(_theData.getOccupancy());
			itsAtoms[i]->setTempFactor(_theData.getTempFactor());
			itsAtoms[i]->setCharge(_theData.getCharge());
			itsAtoms[i]->setIsFullySpecified();
			return;
		}
	}
#ifdef _RESIDUE_DEBUG
	cout << "Couldn't add atom " << _theData.getItem(pdbAtomName);
	cout << "  " << _theData.getItem(sSerial) << endl;
#endif
}

void residue::addAtom(PDBAtomRecord& _theRecord)
{	for(UInt i=0; i < getAtomNameBaseSize(itsType); i++)
	{	if(_theRecord.getAtomName() == getAtomNameBaseItem(itsType,i))
		{
			itsAtoms[i]->setCoords(_theRecord.getAtomCoord());
			itsAtoms[i]->setSerialNumber(_theRecord.getSerial());
			// rigamarole to convert string into a char
			const char* pTheChainID = (_theRecord.getChainID()).c_str();
			char theChainID = *pTheChainID;
			itsAtoms[i]->setChainID(theChainID);
			itsAtoms[i]->setOccupancy(_theRecord.getOccupancy());
			itsAtoms[i]->setTempFactor(_theRecord.getTempFactor());
			itsAtoms[i]->setCharge(0.0);
			itsAtoms[i]->setIsFullySpecified();
			return;
		}
	}
#ifdef _RESIDUE_DEBUG
	cout << "Couldn't add atom " << _theRecord.getAtomName();
	cout << "  " << _theRecord.getSerial() << endl;
#endif
}

void residue::deleteAtom(const UInt _atomIndex)
{
	if (_atomIndex < itsAtoms.size())
	{
		// First, ensure that we're not deleting an
		// atom that's crucial for connectivity
		if (itsAtoms[_atomIndex]->getNumChildren() !=0 )
		{
			cout << "Can't delete atom ";
			cout << getType() << " ";
			cout << itsAtoms[_atomIndex]->getName();
			cout << " because it is not an end atom" << endl;
			return;
		}

		// If I'm going to delete it, i first need to patch up the
		// parent atom to which it is connected

		treeNode* theParent = itsAtoms[_atomIndex]-> getParent();
		treeNode* thePreviousSib = itsAtoms[_atomIndex]->getPreviousSib();
		//cout << getType() << " ";
		//cout << itsAtoms[_atomIndex]->getName();
		//cout << "  pTheParent = " << theParent << " thePreviousSib = " << thePreviousSib << endl;

		// For now, let's only allow the deletion of certain
		// atoms, such as OXT
		if (itsAtoms[_atomIndex]->getName() == "OXT")
		{	// this should be allowed
			if (theParent)
			{	theParent->setChild(0);
			}
			if (thePreviousSib)
			{	thePreviousSib->setChild(0);
			}
			delete itsAtoms[_atomIndex];
			// now we need to "collapse" itsAtoms so that
			// we're not seeing a null pointer in the place
			// of an atom (and so that the functions which
			// ask the residue how many atoms it has are not
			// fooled)
			iterATOM firstAtom;
			firstAtom = itsAtoms.begin();
			itsAtoms.erase(firstAtom + _atomIndex);
			// now the vector no longer contains the null pointer
			// and itsAtoms has been modified such that it looks
			// as though the atom never existed in the first place
		}
		// We're going to allow ASP and GLU to be unprotonated
        if ( (itsAtoms[_atomIndex]->getName() == "HD2" && getType() == "ASP") ||
                     (itsAtoms[_atomIndex]->getName() == "HE2" && getType() == "GLU") )
		{	// this should be allowed
			if (theParent)
			{	theParent->setChild(0);
			}
			if (thePreviousSib)
			{	thePreviousSib->setChild(0);
			}
			delete itsAtoms[_atomIndex];
			iterATOM firstAtom;
			firstAtom = itsAtoms.begin();
			itsAtoms.erase(firstAtom + _atomIndex);
		}

		// We're going to allow HIS to be unprotonated
		if ( (itsAtoms[_atomIndex]->getName() == "HD1" && getType() == "HIS") ||
                     (itsAtoms[_atomIndex]->getName() == "HE2" && getType() == "HIS") )
		{	// this should be allowed
			if (theParent)
			{	theParent->setChild(0);
			}
			if (thePreviousSib)
			{	thePreviousSib->setChild(0);
			}
			delete itsAtoms[_atomIndex];
            iterATOM firstAtom;
			firstAtom = itsAtoms.begin();
			itsAtoms.erase(firstAtom + _atomIndex);
		}

		// We're going to allow ARG to be unprotonated
		if ( (itsAtoms[_atomIndex]->getName() == "HE" && getType() == "ARG") )
		{	// this should be allowed
			if (theParent)
			{	theParent->setChild(0);
			}
			if (thePreviousSib)
			{	thePreviousSib->setChild(0);
			}
			delete itsAtoms[_atomIndex];
			iterATOM firstAtom;
			firstAtom = itsAtoms.begin();
			itsAtoms.erase(firstAtom + _atomIndex);
		}
	}
	else
	{
		cout << "Index out of range in residue::deleteAtom" << endl;
	}

}

void residue::queryChildren(UInt start)
{	itsAtoms[start]->queryChildren();
}


// ************************************************************************
// ************************************************************************
// 	Flags Related
// ************************************************************************
// ************************************************************************

void residue::setHydrogensOn(const bool _hydrogensOn )
{	if(_hydrogensOn != hydrogensOn)
	{	if(hydrogensOn)
		{	hydrogensOn = _hydrogensOn;
			// add code to build hydrogens here
		}
		else
		{	hydrogensOn = _hydrogensOn;
			// add code to build hydrogens here
		}
	}

}

void residue::setPolarHydrogensOn(const bool _polarHydrogensOn)
{	if(_polarHydrogensOn!= polarHydrogensOn)
	{	if(polarHydrogensOn)
		{	polarHydrogensOn = _polarHydrogensOn;
			// build polar H atoms here
		}
		else
		{	polarHydrogensOn = _polarHydrogensOn;
			// build here
		}
	}
}

// ************************************************************************
// ************************************************************************
//	Build the residue dataBase
// ************************************************************************
// ************************************************************************

void residue::buildDataBase()
{
	dataBase.resize(0);
	// first, build the residueTemplates
	buildResidueDataBaseAminoAcids();
	buildDihedralDataBase();
	buildDataBaseFromPrep();
	interpretBondingPattern();
	buildRotamerLib();
	dataBaseBuilt = true;

	// now that database has been constructed, build electrostatics forcefields
	residueTemplate::itsAmberElec.buildElectrostatics();
	residueTemplate::itsHelixPropensity.buildDatabase();
	//cout << "Done!!!" << endl;
}

// *****************************************************************************
// *****************************************************************************
//	this function is ONLY for buildDataBase(..)
// *****************************************************************************
// *****************************************************************************

void residue::buildResidueDataBaseAminoAcids()
{	// initialize the dataBase
	//cout << "Initializing the residue type database...." << endl;

	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	
	string aaLib = "/data/mol.lib";
	string iFile;
	ifstream inFile;

	// first build the amino acids
	iFile = path + aaLib;
	inFile.open(iFile.c_str());

	if(!inFile)
	{	cout << "Error: unable to open input file: "
			 << iFile << endl;
		exit (1);
	}
	
	string linebuffer;

	// Used to assign the type index
	int resTypeIndex = 0;
	// used to assign the atom name within each residue
	int atomNameIndex = 0;
	
	residueTemplate tempTemplate;
	// a strange string
	const string strange = "*^{POLHGYHG><MDEZF!";
	string currentResName = strange;
	
	while (getline(inFile,linebuffer,'\n'))
	{
		PDBAtomRecord currentRecord(linebuffer);
		if (currentResName != currentRecord.getResName())
	 	{	if (currentResName != strange)
			{	dataBase.push_back(tempTemplate);
				resTypeIndex++;
				tempTemplate.reset();
			}
			currentResName = currentRecord.getResName();
			tempTemplate.typeIndex = resTypeIndex;
			tempTemplate.typeString = currentResName;
			atomNameIndex = 0;
		//	cout << currentResName << " " << resTypeIndex << endl;
		}
		atom tempAtom(currentRecord);
		tempAtom.itsResType = resTypeIndex;
		tempTemplate.atomNameList.push_back(currentRecord.getAtomName());
		// set the corresponding index atom name
		tempAtom.itsName = atomNameIndex;
		atomNameIndex++;
		tempTemplate.atomList.push_back(tempAtom);
		// initialize the isMainChain database
		tempTemplate.isMainChain.push_back(false);
	}
}

// *****************************************************************************
// *****************************************************************************
//  buildDataBaseCC is ONLY for buildDataBase
// *****************************************************************************
// *****************************************************************************

void residue::buildDihedralDataBase()
{	
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	
	string aaDat = "/data/mol/";

	path += aaDat;
	string iFile;
	ifstream inFile;

	// Pointer is used to make the process efficient
	// But any modification should be conducted very carefully!!!
	residueTemplate* pCurrentResTemplate = 0;

	for(UInt i=0; i<dataBase.size(); i++)
	{	pCurrentResTemplate = &dataBase[i];
		//cout << pCurrentResTemplate->typeString << endl;
		iFile = path + pCurrentResTemplate->typeString + ".dat";
		inFile.open(iFile.c_str());
		//if(!inFile)
		//{	cout << "Error: unable to open input file: "
		//	     << iFile << endl;
		//	exit (1);
		//}
		
		// the available headers without spaces in between
		string mcn = "mainchainatoms";
		string bpt = "sidechainbranchpoints";
		string chi = "defineddihedralangles";
		
		char ch;
		// used to convert a char to a string
		string strCh;
		strCh.resize(1);
		// string buffer
		string strBuf;
		strBuf.resize(0);
		// header buffer
		string header;
		header.resize(0);
		// atom index buffer
		int tempInt;
		// strings buffer
		StrVec strVect;
		strVect.resize(0);

		// '!' designates a header
		// contents follow the header
		// '#' designates a comment -- ignored
		// comments are not allowed to follow the header in the same line

		while(inFile.get(ch))
		{	if(ch == '!')
			{	// build a new header
				header.resize(0);
				// headers should be within one line
				while(inFile.get(ch) && ch != '\n')
				{	// ignore all the spaces and the tabs if there is any
					if(ch != '\t' && ch != ' ')
					{	strCh[0] = ch;
						header.append(strCh);
					}
				}
			}
			else if(ch == '#')
			{	// ignore the comments
				for(; inFile.get(ch) && ch !='\n'; );
			}
			else if(header == mcn)
			{	if(ch != ' ' && ch != '\t' && ch != '\n')
				{	strCh[0] = ch;
					strBuf.append(strCh);
				}
				else
				{	if(strBuf.size() != 0)
					{	tempInt = pCurrentResTemplate->getAtomIndexOf(strBuf);
						if(tempInt != -1)
						{	pCurrentResTemplate->mainChain.push_back(tempInt);
						}
						strBuf.resize(0);
					}
				}
			}
			else if(header == bpt)
			{	if(ch !=' ' && ch != '\t' && ch != '\n')
				{	strCh[0] = ch;
					strBuf.append(strCh);
				}
				else
				{	if(strBuf.size() != 0)
					{	tempInt = pCurrentResTemplate->getAtomIndexOf(strBuf);
						{	if(tempInt != -1)
							{	pCurrentResTemplate->branchPoints.push_back(tempInt);
							}
						}
						strBuf.resize(0);
					}
				}
			}
			else if(header == chi)
			{	if(ch != ' ' && ch !='\t' && ch !='\n')
				{	strCh[0] = ch;
					strBuf.append(strCh);
				}
				else
				{	if(strBuf.size() != 0)
					{	strVect.push_back(strBuf);
						if(ch == '\n')
						{	// process the strVect
							pCurrentResTemplate->addChiDefinitions(strVect);
							// reset the strVect after the processing
							strVect.resize(0);
						}
						strBuf.resize(0);
					}
				}
			}
		}
		inFile.close();
		inFile.clear();
	}
}

void residue::buildDataBaseFromPrep()
{	
	string itsFileName = "amber.prep";
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	path += "/data/";
	string iFile = path + itsFileName;
	ifstream inFile;
	string currentLine;
	StrVec parsedStrings;
	vector< vector< string > > atomNames;
	vector< vector< string > > loopedAtoms;
	vector< vector< int > > bondIndices;
	vector< vector< int > > bondedIndices;
	vector< vector< string > > vdwTypes;
	vector< string > resNames;
	int bondIndex, bondedIndex;
	parsedStrings.resize(0);

	inFile.open(iFile.c_str());
	if (!inFile)
	{
		cout << "Error: unable to open input file: " ;
		cout << iFile << endl;
		exit (1);
	}
	bool loop = false;
	while (getline (inFile, currentLine))
	{
		parsedStrings=Parse::parse(currentLine);
		if (loop && parsedStrings.size() == 0){loop = false;}
		if (parsedStrings.size() > 1){
			if (loop){
				loopedAtoms[loopedAtoms.size()-1].push_back(parsedStrings[0]);
				loopedAtoms[loopedAtoms.size()-1].push_back(parsedStrings[1]);
			}
			if (parsedStrings[1].size() == 3){
				if (parsedStrings[1].compare("INT") == 0){
					string tmpStr;
					string tmpChar;
					tmpChar.resize(1);
					tmpStr.resize(0);
					for (UInt i = 0; i < parsedStrings[0].size(); i++) // loop starts after the # mark
					{
						tmpChar = parsedStrings[0][i];
						tmpStr.append(tmpChar);
					}
					resNames.push_back(tmpStr);
					StrVec tempVec;
					loopedAtoms.push_back(tempVec);
				}
			}
		}
		else if (parsedStrings.size() == 1){
			if (parsedStrings[0].compare("LOOP") == 0){loop = true;}
		}
		if (parsedStrings.size() == 11)
		{
			if (parsedStrings[2][0] != 'D'){
			
				// Get bonding index
				bondIndex = 0;
				sscanf(parsedStrings[0].c_str(), "%d", &bondIndex);
				
				// Get bonded index
				bondedIndex = 0;
				sscanf(parsedStrings[4].c_str(), "%d", &bondedIndex);
				
				if (resNames.size() != atomNames.size()) // if there are fewer atom vectors than their are residues
				{
					vector <string> atomList;
					atomList.push_back(parsedStrings[1]);
					atomNames.push_back(atomList);
					vector <string> vdwList;
					vdwList.push_back(parsedStrings[2]);
					vdwTypes.push_back(vdwList);
					vector <int> bondList;
					bondList.push_back(bondIndex);
					bondIndices.push_back(bondList);
					vector <int> bondedList;
					bondedList.push_back(bondedIndex);
					bondedIndices.push_back(bondedList);
				}
				else
				{
					UInt size = atomNames.size();
					atomNames[size - 1].push_back(parsedStrings[1]);
					vdwTypes[size - 1].push_back(parsedStrings[2]);
					bondIndices[size - 1].push_back(bondIndex);
					bondedIndices[size - 1].push_back(bondedIndex);
				}
			}
		}
		parsedStrings.resize(0);
	}
	inFile.close();
	inFile.clear();

	// Pointer is used to make the process efficient
	// But any modification should be conducted very carefully!!!
	residueTemplate* pCurrentResTemplate = 0;

	for(UInt i=0; i<dataBase.size(); i++)
	{	pCurrentResTemplate = &dataBase[i];
		string res = pCurrentResTemplate->typeString;
		int resIndex;
		for (UInt j = 0; j < resNames.size(); j++)
		{
			if (res.compare(resNames[j]) == 0){resIndex = j; break;}
		}
		
		// string buffer
		string strBuf;
		strBuf.resize(0);
		// atom index buffer
		int tempInt;
		// strings buffer
		StrVec typeVect;
		typeVect.resize(0);

		// build mainchain atoms as backbone atoms or first four atoms in mol.lib
		if (pCurrentResTemplate->mainChain.size() == 0){
			tempInt = pCurrentResTemplate->getAtomIndexOf("N");
			if(tempInt != -1){pCurrentResTemplate->mainChain.push_back(tempInt);}
			else{pCurrentResTemplate->mainChain.push_back(0);}
			tempInt = pCurrentResTemplate->getAtomIndexOf("CA");
			if(tempInt != -1){pCurrentResTemplate->mainChain.push_back(tempInt);}
			else{pCurrentResTemplate->mainChain.push_back(1);}
			tempInt = pCurrentResTemplate->getAtomIndexOf("C");
			if(tempInt != -1){pCurrentResTemplate->mainChain.push_back(tempInt);}
			else{pCurrentResTemplate->mainChain.push_back(2);}
			tempInt = pCurrentResTemplate->getAtomIndexOf("O");
			if(tempInt != -1){pCurrentResTemplate->mainChain.push_back(tempInt);}
			else{pCurrentResTemplate->mainChain.push_back(3);}
		}
		
		// build branchpoint atom as alpha carbon or second atom in mol.lib
		if (pCurrentResTemplate->branchPoints.size() == 0){
			tempInt = pCurrentResTemplate->getAtomIndexOf("CA");
			if(tempInt != -1){pCurrentResTemplate->branchPoints.push_back(tempInt);}
			else{pCurrentResTemplate->branchPoints.push_back(1);}
		}
		
		// build vdw type, connectivity from amber prep file
		for (UInt j = 0; j < atomNames[resIndex].size(); j++)
		{
			int atomIndex;
			string atomName = pCurrentResTemplate->getAtomNameOf(j);
			for (UInt k = 0; k < atomNames[resIndex].size(); k++)
			{
				if (atomName.compare(atomNames[resIndex][k]) == 0){atomIndex = k; break;}
			}
			
			// Load vdw type into atom
			typeVect.push_back(atomName);
			typeVect.push_back(vdwTypes[resIndex][atomIndex]);
			typeVect.push_back(vdwTypes[resIndex][atomIndex]);
			pCurrentResTemplate->addAtomTypeDefinitions(typeVect);
			typeVect.resize(0);
			
			// add to connectivity vector for each atom
			pCurrentResTemplate->connectivity.push_back(0);
			tempInt = pCurrentResTemplate->getAtomIndexOf(atomName);
			if(tempInt != -1){
				pCurrentResTemplate->connectivity.push_back(tempInt+1);
			}
			for (UInt k = 0; k < bondedIndices[resIndex].size(); k++)
			{
				int otherAtomIndex;
				string otherAtomName = pCurrentResTemplate->getAtomNameOf(k);
				for (UInt l = 0; l < atomNames[resIndex].size(); l++)
				{
					if (otherAtomName.compare(atomNames[resIndex][l]) == 0){otherAtomIndex = l; break;}
				}
				if (bondIndices[resIndex][atomIndex] == bondedIndices[resIndex][otherAtomIndex])
				{
					tempInt = pCurrentResTemplate->getAtomIndexOf(atomNames[resIndex][otherAtomIndex]);
					if(tempInt != -1){
						pCurrentResTemplate->connectivity.push_back(tempInt+1);
					}
				}
			}
			if (loopedAtoms[resIndex].size() > 1){
				for (UInt k = 1; k < loopedAtoms[resIndex].size(); k++)
				{
					if (atomName.compare(loopedAtoms[resIndex][k]) == 0)
					{
						tempInt = pCurrentResTemplate->getAtomIndexOf(loopedAtoms[resIndex][k-1]);
						if(tempInt != -1){
							pCurrentResTemplate->connectivity.push_back(tempInt+1);
						}
					}
					k++;
				}
			}
		}
	}
}

/*void residue::buildAutoRotamerLib()
{
	
	rotamerLib* temp;
	UIntVec angles;
	for(UInt i=0; i<dataBase.size(); i++)
	{	// if chis are defined
		if( dataBase[i].chiDefinitionsNonempty() )
		{
			temp = new rotamerLib(dataBase[i].branchPoints.size());
			// since it is a pointer, the changes below done to
			// temp is the same as the ones done to rotamerlib
			ASSERT(temp != 0);
			dataBase[i].itsRotamerLibs.push_back(temp);
			for (UInt j = 0; j < dataBase[i].branchPoints.size(); j++)
			{
				for () 
						temp->addRotamer(j, angles);
						angles.resize(0);
		}
	}			
}*/

void residue::buildRotamerLib()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	
	string aaLib = "/data/rotamerLib/";

	path += aaLib;
	rotamerLib* temp;
	string filename;
	string iFile;
	ifstream inFile;
	char ch;
	string strCh;
	strCh.resize(1);
	string strBuf;
	strBuf.resize(0);
	StrVec strVec;
	strVec.resize(0);
	const StrVec& jason = strVec;

	for(UInt i=0; i<dataBase.size(); i++)
	{	// if chis are defined
		if( dataBase[i].chiDefinitionsNonempty() )
		{	filename = dataBase[i].typeString + ".rot";
			iFile = path + filename;
			inFile.open(iFile.c_str());
			if(!inFile)
			{	cout << "Warning: rotamer library for ";
				cout << dataBase[i].typeString;
				cout << " not found " << endl;
				// then skip
				continue;
			}
		
			while(inFile.get(ch))
			{
				if(ch == '!')
				{
					for(; inFile.get(ch) && ch != '\n'; );
					temp = new rotamerLib(dataBase[i].branchPoints.size());
					// since it is a pointer, the changes below done to
					// temp is the same as the ones done to rotamerlib
					ASSERT(temp != 0);
					dataBase[i].itsRotamerLibs.push_back(temp);
				}
				else if(ch != '\n')
				{	if(ch != ' ' && ch != '\t')
					{	strCh[0]=ch;
						strBuf.append(strCh);
					}
					else
					{	if(strBuf.size() != 0)
						{	strVec.push_back(strBuf);
							strBuf.resize(0);
						}
					}
				}
				else
				{
/*	
					cout << "Im in the else phrase" << endl;
					cout << "strBuf.size() = " << strBuf.size() << endl;
					cout << "strVec.size() = " << strVec.size() << endl;
*/
					if( strBuf.size() != 0)
					{	strVec.push_back(strBuf);
						strBuf.resize(0);
					}
					if(strVec.size() != 0)
					{	
				//		cout << strVec.size() << "  ";
				//		cout << temp << endl;
						temp->addRotamer(jason);
					}
					strVec.resize(0);
				}
			}
			inFile.close();
			inFile.clear();
		}
	}
}

		
// *****************************************************************************
// *****************************************************************************
//	dataBase information accessors
// *****************************************************************************
// *****************************************************************************

UInt residue::getDataBaseSize()
{	return dataBase.size();
}

string residue::getDataBaseItem(const UInt _itemIndex)
{	if(_itemIndex < dataBase.size())
	{	return dataBase[_itemIndex].typeString;
	}
	else
	{	cout << " _itemIndex is incompatible with residue::dataBase " << endl;
		return "UNK";
	}
}


// *****************************************************************
// *****************************************************************
//	ATOM NAME BASE ACCESS
// *****************************************************************
// *****************************************************************

UInt residue::getAtomNameBaseSize(const UInt _resType)
{	return dataBase[_resType].atomNameList.size();
}

string residue::getAtomNameBaseItem(const UInt _resType, const UInt _index)
{	if(_resType < dataBase.size())
	{	if(_index < dataBase[_resType].atomNameList.size())
		{	return dataBase[_resType].atomNameList[_index];
		}
		else
		{	cout << "index is incompatible with atomNameList " << endl;
			return "UNK";
		}
	}
	else
	{	cout << "index is larger than size of atomNameList" << endl;
		return "UNK";
	}
}

// ************************************************************************
// ************************************************************************
// 	Define or Query Type of a Residue
// ************************************************************************
// ************************************************************************

void residue::defineType( const string& _aaType)
{	if(!dataBaseBuilt)
	{	buildDataBase();
	}
	setType( _aaType);
}


void residue::setType( const string& _aaType)
{	bool flag = false;

	for(UInt i=0; i<dataBase.size(); i++)
	{	if(_aaType == dataBase[i].typeString)
		{	itsType = i;
			flag = true;
		}
	}

	if(!flag)
	{	cout << "Error: _aaType: [" << _aaType << "]"
			 << " is incompatible with dataBase " << endl;
		cout << "       residue type is not set or changed " << endl;
		cout << "Error reported by: residue::setType(const string&) "
			 << endl;
	}
}

string residue::getType() const
{	if(itsType < dataBase.size())
	{	return dataBase[itsType].typeString;
	}
	else
	{	cout << "Error: itsType incompatible with dataBase " << endl;
		cout << "Error reported by residue::getType() " << endl;
		return "UNK";
	}
	// dummy return to keep compiler shut-up
	return "UNK";
}

string residue::getType(UInt resType)
{	if(resType < dataBase.size())
	{	return dataBase[resType].typeString;
	}
	else
	{	cout << "Error: itsType incompatible with dataBase " << endl;
		cout << "Error reported by residue::getType() " << endl;
		return "UNK";
	}
	// dummy return to keep compiler shut-up
	return "UNK";
}

void residue::interpretBondingPattern()
{	residueTemplate* pCurrentResTemplate;
	vector <UInt> theConnectivityVector;
	vector< vector < UInt> > theBondingPatternVector;
	for (UInt i=0; i<dataBase.size(); i++)
	{	pCurrentResTemplate = &dataBase[i];
		theConnectivityVector.resize(0);
		theConnectivityVector = pCurrentResTemplate->connectivity;
		theBondingPatternVector.resize(0);
		// start the interpretation
		// First, check how many atoms we're dealing with...
		UInt numatoms = pCurrentResTemplate->atomList.size();
		//cout << i << " " << numatoms << endl;
		//for (UInt j=0; j<theConnectivityVector.size(); j++)
		//{	cout << theConnectivityVector[j] << " ";
		//}
		//cout << endl;
		for (UInt j=0; j<numatoms; j++)
		{	vector<UInt> tempBPvector;
			theBondingPatternVector.push_back(tempBPvector);
		}
		vector<UInt> tempVec;
		UInt counter = 0;
		for (UInt j=0; j<theConnectivityVector.size();j++)
		{	
			if (theConnectivityVector[j] == 0)
			{
				if (tempVec.size())
				{	
#ifdef _CONNECTIVITY_DEBUG
					cout << counter << " : ";
#endif
					for (UInt k=1;k<tempVec.size();k++)
					{	
#ifdef _CONNECTIVITY_DEBUG
						cout << tempVec[k] << " ";
#endif
						theBondingPatternVector[counter].push_back(tempVec[k]);
						theBondingPatternVector[tempVec[k]].push_back(counter);
					}
#ifdef _CONNECTIVITY_DEBUG
					cout << endl;
#endif
					tempVec.resize(0);
					counter++;
				}
			}
			else
			{	tempVec.push_back(theConnectivityVector[j] - 1);
			}
		}
#ifdef _CONNECTIVITY_DEBUG
		cout << "---------------------" << endl;
		for (UInt j=0; j<theBondingPatternVector.size(); j++)
		{	
			cout << j << " : ";
			for (UInt k=0; k<theBondingPatternVector[j].size(); k++)
			{		cout << theBondingPatternVector[j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
#endif
		pCurrentResTemplate->itsBondingPattern = theBondingPatternVector;
	}
}

//#define _MUTATE_DEBUG
//#define _CHI_DEBUG

// ***********************************************************************
// ************************************************************************
// Mutate function
// ************************************************************************
// ************************************************************************
residue* residue::mutate(const UInt _newTypeIndex)
{
	residue* newAA = new residue( _newTypeIndex, hydrogensOn );
	bool betapivot = false;
	if (itsType == newAA->itsType){betapivot = true;}
	UInt numbpt = getNumBpt(_newTypeIndex);
	for (UInt i=0; i<numbpt; i++)
	{
		dblVec Coord_A_target = getMainChain(i)->getCoords();
		dblVec Coord_C_target = getMainChain(i+2)->getCoords();
		dblVec Coord_B_target = getMainChain(i+1)->getCoords();

#ifdef _MUTATE_DEBUG
		cout << "Coord_A_target : " << Coord_A_target << endl;
		cout << "Coord_B_target : " << Coord_B_target << endl;
		cout << "Coord_C_target : " << Coord_C_target << endl;
#endif

		dblVec Coord_A_new(3);
		dblVec Coord_B_new(3);
		dblVec Coord_C_new(3);
		Coord_B_new = newAA->getMainChain(i+1)->getCoords();

		// This equivalences the coordinates of the branchpoint atoms
		// in the new amino acid to that of the original amino acid
		newAA->translate(Coord_B_target - Coord_B_new);

		Coord_A_new = newAA->getMainChain(i)->getCoords();
		Coord_C_new = newAA->getMainChain(i+2)->getCoords();
		Coord_B_new = newAA->getMainChain(i+1)->getCoords();

#ifdef _MUTATE_DEBUG
		cout << "After translation : " << endl;
		cout << "Coord_A_new : " << Coord_A_new << endl;
		cout << "Coord_B_new : " << Coord_B_new << endl;
		cout << "Coord_C_new : " << Coord_C_new<< endl;
#endif
		//Using a backbone atom plane determine necessary rotation matrix
		//to align mutant residue to old residue's backbone position
		vector<dblVec> eTarget;
		dblVec vec_B_A_target = Coord_A_target - Coord_B_target;
		double norm = sqrt(CMath::dotProduct(vec_B_A_target,vec_B_A_target));
		eTarget.push_back( vec_B_A_target / norm );
		dblVec vec_B_C_target = Coord_C_target - Coord_B_target;
		dblVec jason = vec_B_C_target - CMath::dotProduct(vec_B_C_target,eTarget[0]) * eTarget[0];
		norm = sqrt(CMath::dotProduct(jason,jason));
		eTarget.push_back(jason/norm);
		eTarget.push_back(CMath::cross(eTarget[0],eTarget[1]));

		vector<dblVec> eNew;
		dblVec tempVec = Coord_A_new- Coord_B_new;
		norm = sqrt(CMath::dotProduct(tempVec,tempVec));
		eNew.push_back(tempVec/norm);
		tempVec = Coord_C_new - Coord_B_new;
		tempVec = tempVec - CMath::dotProduct(tempVec,eNew[0])*eNew[0];
		norm = sqrt(CMath::dotProduct(tempVec,tempVec));
		eNew.push_back(tempVec/norm);
		eNew.push_back(CMath::cross(eNew[0],eNew[1]));

		dblMat RotMat(3,3,0.0);
		for(UInt index1=0; index1<3; index1++)
		{	for(UInt index2=0; index2<3; index2++)
			{	RotMat[index2][index1] = 0;
				for(UInt index3=0; index3<3; index3++)
				{	RotMat[index2][index1] += eNew[index3][index1]*eTarget[index3][index2];
				}
			}
		}

#ifdef _MUTATE_DEBUG
		cout << " e1 : " << eTarget[0] << endl;
		cout << " e2 : " << eTarget[1] << endl;
		cout << " e3 : " << eTarget[2] << endl;
		cout << " eNew : " << endl;
		cout << " e1 : " << eNew[0] << endl;
		cout << " e2 : " << eNew[1] << endl;
		cout << " e3 : " << eNew[2] << endl;
		cout << "rotation matrix : " << endl;
		cout << RotMat << endl;
#endif

		// newAA->getMainChain(i+1) is the pivot-point atom.
		// we want to rotate all the sidechain atoms which branch
		// from that atom only.
		atom* pivotAtom = newAA->getMainChain(i+1);

		if (i==0)
		{
			newAA->rotate_new(pivotAtom,newAA->getMainChain(0), RotMat);
		}
		else
		{
			newAA->rotate_new(pivotAtom, RotMat);
		}
	}

	//position new residue with proper position context
	if (pItsNextRes)
	{
		newAA->setNextRes(pItsNextRes);
		pItsNextRes->setPrevRes(newAA);
	}
	else
	{
		newAA->setNextRes(pItsNextRes);
	}
	if (pItsPrevRes)
	{
		newAA->setPrevRes(pItsPrevRes);
		pItsPrevRes->setNextRes(newAA);
	}
	else
	{
		newAA->setPrevRes(pItsPrevRes);
	}
	newAA->setResNum(itsResNum);

	// ensure that we have not modified the main chain
	// coordinates in any way
	for (UInt i=0;i<dataBase[itsType].mainChain.size(); i++)
	{
		newAA->getMainChain(i)->setCoords( getMainChain(i)->getCoords());
	}

	// When mutating the same residue in place of old
	// make sure same amino acid has near identical Calpha-Cbeta angle
	// and Cbeta position relative to backbone
	if (betapivot)
	{
		if (itsAtoms.size() > 4){
			if (itsAtoms[4]->getName() == "CB"){
				//newAA->getAtom(4)->setCoords(itsAtoms[4]->getCoords());
				newAA->setBetaChi(getBetaChi());
			}
		}
	}
	
	
	// now, since we're changing the backbone coordinates, make
	// sure that the amide hydrogen H is in the right place...
	// all other hydrogens should have been taken care of by the
	// code above which sets the coordinates of the branchpoint atoms.
	if (hydrogensOn)
	{
		newAA->alignAmideProtonToBackbone();
	}
	newAA->setMoved();
	return newAA;
}

residue* residue::mutateNew(const UInt _newTypeIndex) // is generalized to support non-backbone mutation using fortan alignment algorithm and supports terminal residues
{
	residue* newAA; int maxsize = 3; bool sameAA = false;
	if (itsType == _newTypeIndex && (isD(itsType) || isL(itsType))){maxsize = 4; sameAA = true;}
	if (!pItsPrevRes && (isL(_newTypeIndex) || isD(_newTypeIndex) || isG(_newTypeIndex))){newAA = new residue( _newTypeIndex+Nterm, true);} //Nterminal mutation
	if (!pItsNextRes && (isL(_newTypeIndex) || isD(_newTypeIndex) || isG(_newTypeIndex))){newAA = new residue( _newTypeIndex+Cterm, true);} //Cterminal mutation
	else{newAA = new residue( _newTypeIndex, true); }
	
    // Load atoms for alignment onto position into double array for fit
	double coord1[maxsize*3], coord2[maxsize*3]; int list1[maxsize]; int list2[maxsize]; dblVec oldCoords(3), newCoords(3);
	for (int i=0; i<maxsize; i++)
	{	
		if (i == 3 && sameAA){
			oldCoords = itsAtoms[4]->getCoords();
			newCoords = newAA->itsAtoms[4]->getCoords();
		}
		else{
			oldCoords = getMainChain(i)->getCoords();
			newCoords = newAA->getMainChain(i)->getCoords();
		}
		for (int j=0; j<3; j++)
		{
			coord1[(i*3) + j] = oldCoords[j];
			coord2[(i*3) + j] = newCoords[j];
		}
		list1[i] = i+1;
		list2[i] = i+1;
	}
	
	// Calculate best fit of atoms to get rotation matrix and centroid from fortran algorithm based on Machlachlan
	double rotmat[9]; double centroid1[3]; double centroid2[3]; double rmsd = 0; double coord3[maxsize*3]; int ierr = 0;
	bestfit_(coord1, &maxsize, coord2, &maxsize, &maxsize, coord3, list1, list2, &rmsd, &ierr, rotmat, centroid1, centroid2);
	
	// Load rotation vector into rotation matrix and double array into double Vector
	dblMat rotMat(3,3,3);
    for (UInt i=0; i<3; i++)
    {	for (UInt j=0; j<3; j++)
		{
			rotMat[i][j] = rotmat[(j*3) + i];
		}
    }
	dblVec centroidOne(3); centroidOne[0] = centroid1[0]; centroidOne[1] = centroid1[1]; centroidOne[2] = centroid1[2];
	dblVec centroidTwo(3); centroidTwo[0] = centroid2[0]*-1; centroidTwo[1] = centroid2[1]*-1; centroidTwo[2] = centroid2[2]*-1;
    
	//translate to origin, rotate using rotation matrix from fit and move to position of starting residue
	newAA->translate(centroidTwo); newAA->transform(rotMat); newAA->translate(centroidOne);

	//position new residue with proper position context
	if (pItsNextRes){ newAA->setNextRes(pItsNextRes); pItsNextRes->setPrevRes(newAA);}
	else{newAA->setNextRes(pItsNextRes);}
	if (pItsPrevRes){newAA->setPrevRes(pItsPrevRes); pItsPrevRes->setNextRes(newAA);}
	else{newAA->setPrevRes(pItsPrevRes);}
	newAA->setResNum(itsResNum);

	// move carbonyl oxygen, amide proton into position and activate for energy calculation
	if(!isCofactor(_newTypeIndex)){
		if (pItsNextRes){newAA->getMainChain(3)->setCoords(getMainChain(3)->getCoords());}
		if (pItsPrevRes){newAA->alignAmideProtonToBackbone();}
	}
	newAA->setMoved();

	return newAA;
}

residue* residue::superimposeGLY()
{
	residue* newAA = new residue(7,true);

	dblVec Coord_A_target = getMainChain(0)->getCoords();
	dblVec Coord_B_target = getMainChain(1)->getCoords();
	dblVec Coord_C_target = getMainChain(2)->getCoords();

	dblVec Coord_A_new(3);
	dblVec Coord_B_new(3);
	dblVec Coord_C_new(3);

	Coord_B_new = newAA->getMainChain(1)->getCoords();

	// This equivalences the coordinates of the branchpoint atoms
	// in the new amino acid to that of the original amino acid
	newAA->translate(Coord_B_target - Coord_B_new);

	Coord_A_new = newAA->getMainChain(0)->getCoords();
	Coord_B_new = newAA->getMainChain(1)->getCoords();
	Coord_C_new = newAA->getMainChain(2)->getCoords();

	vector<dblVec> eTarget;
	dblVec vec_B_A_target = Coord_A_target - Coord_B_target;
	// norm = length of vec_B_A_target
	double norm = sqrt(CMath::dotProduct(vec_B_A_target,vec_B_A_target));
	eTarget.push_back( vec_B_A_target / norm );
	dblVec vec_B_C_target = Coord_C_target - Coord_B_target;
	dblVec jason = vec_B_C_target - CMath::dotProduct(vec_B_C_target,eTarget[0]) * eTarget[0];
	norm = sqrt(CMath::dotProduct(jason,jason));
	eTarget.push_back(jason/norm);
	eTarget.push_back(CMath::cross(eTarget[0],eTarget[1]));

	vector<dblVec> eNew;
	dblVec tempVec = Coord_A_new- Coord_B_new;
	norm = sqrt(CMath::dotProduct(tempVec,tempVec));
	eNew.push_back(tempVec/norm);
	tempVec = Coord_C_new - Coord_B_new;
	tempVec = tempVec - CMath::dotProduct(tempVec,eNew[0])*eNew[0];
	norm = sqrt(CMath::dotProduct(tempVec,tempVec));
	eNew.push_back(tempVec/norm);
	eNew.push_back(CMath::cross(eNew[0],eNew[1]));

	dblMat RotMat(3,3,0.0);
	for(UInt index1=0; index1<3; index1++)
	{	for(UInt index2=0; index2<3; index2++)
		{	RotMat[index2][index1] = 0;
			for(UInt index3=0; index3<3; index3++)
			{	RotMat[index2][index1] += eNew[index3][index1]*eTarget[index3][index2];
			}
		}
	}
	// we want to rotate all the sidechain atoms which branch
	// from that atom only.
	atom* pivotAtom = newAA->getMainChain(1);
	newAA->rotate_new(pivotAtom,newAA->getMainChain(0), RotMat);


	// ensure that we have not modified the main chain
	// coordinates in any way

	for (UInt i=0;i<dataBase[itsType].mainChain.size(); i++)
	{
		newAA->getMainChain(i)->setCoords( getMainChain(i)->getCoords());
	}
	return newAA;
}


residue* residue::fixBrokenResidue()
{
	// First, make sure that backbone is fully defined
	// Assumption:: if the coordinates of the main chain have
	// not been changed from their defaults (in the residue
	// type base) they are missing.  In this case, we can't really
	// fix the residue simply using the mutagenesis method.

	bool OK = true;
	residue* fixedres = 0;
	for (UInt i=0;i<dataBase[itsType].mainChain.size(); i++)
	{	
		if ( getMainChain(i)->getCoords() == dataBase[itsType].atomList[dataBase[itsType].mainChain[i]].getCoords() )
		{
			OK = false;
			break;
		}
	}
	if (OK)		
	{
		fixedres = mutate(itsType);
		cout << "Residue fixed successfully" << endl;
	}
	else
	{
		cout << "Backbone Atoms Missing: cannot fix this residue" << endl;
		fixedres = this;
	}
	return fixedres;
}

void residue::alignAmideProtonToBackbone()
{
	if (pItsPrevRes)
	{	
		// get the index of "H" if it exists
		string amideProton = "H";
		UInt HIndex = 99;
		for (UInt i = 0; i < itsAtoms.size(); i++)
		{
			string itsName = itsAtoms[i]->getName();
			if (itsName == amideProton)
			{
				HIndex = i;
			}
		}
		if (HIndex != 99)
		{
			atom* pHNew = getAtom(HIndex);
			dblVec newCoords(3);
			dblVec prevCarbon = pItsPrevRes->getCoords("C");
			dblVec prevOxygen = pItsPrevRes->getCoords("O");
			dblVec itsNitrogen = getCoords("N");
			newCoords[0] = ((prevCarbon[0]-prevOxygen[0])*0.8)+itsNitrogen[0];
			newCoords[1] = ((prevCarbon[1]-prevOxygen[1])*0.8)+itsNitrogen[1];
			newCoords[2] = ((prevCarbon[2]-prevOxygen[2])*0.8)+itsNitrogen[2];
			pHNew->setCoords(newCoords); // translate Hydrogen to same plane as C=O of prev residue and axis of Nitrogen
		}
	}
}
// ************************************************************************
// ************************************************************************
// Dihedral Related
// ************************************************************************
// ************************************************************************
void residue::setRotamer(const UInt _lib, const UInt _bpt, const UInt _rotamer)
{	//cout << " set rotamer at residue level " << endl;
	itsSidechainDihedralAngles[_bpt] = dataBase[itsType].itsRotamerLibs[_lib]->getAngles(_bpt,_rotamer);
	for (UInt i=0;i<itsSidechainDihedralAngles[_bpt].size();i++)
	{	setChi(_bpt,i,itsSidechainDihedralAngles[_bpt][i]);
	}
	//calculateSidechainDihedralAngles();
	setMoved();
}

void residue::setRotamer(const UInt _bpt, const DouVec _rotamer)
{	//cout << " set rotamer at residue level " << endl;
	itsSidechainDihedralAngles[_bpt] = _rotamer;
	for (UInt i=0;i<itsSidechainDihedralAngles[_bpt].size();i++)
	{	setChi(_bpt,i,itsSidechainDihedralAngles[_bpt][i]);
	}
	//calculateSidechainDihedralAngles();
	setMoved();
}

void residue::setRotamerWithCheck(const UInt _lib, const UInt _bpt, const UInt _rotamer)
{	//cout << " set rotamer (with checking) at residue level " << endl;
	if ( dataBase[itsType].chiDefinitions.size() > _bpt &&
		dataBase[itsType].chiDefinitions[_bpt].size() >= 4)
	{	if (dataBase[itsType].itsRotamerLibs[_lib]->rotamersExist(_bpt,_rotamer))
		{
			itsSidechainDihedralAngles[_bpt] = dataBase[itsType].itsRotamerLibs[_lib]->getAngles(_bpt,_rotamer);
			for (UInt i=0;i<itsSidechainDihedralAngles[_bpt].size();i++)
			{
			//	cout << "setRotamerWithCheck invoked successfully." << i << ", " << itsSidechainDihedralAngles[_bpt][i] << endl;
			//	for (UInt j = 0; j < itsSidechainDihedralAngles[_bpt].size(); j++) cout << " " << itsSidechainDihedralAngles[_bpt][j];
			//	cout <<endl;
			setChi(_bpt,i,itsSidechainDihedralAngles[_bpt][i]);
              //  for (UInt j = 0; j < itsSidechainDihedralAngles[_bpt].size(); j++) cout << " " << itsSidechainDihedralAngles[_bpt][j];
              //    cout <<endl;
			}
		}
		//else cout << "ERROR in setRotamerWithCheck...\n rotamer " << _bpt << ", " << _rotamer << " does not exist." << endl;
	}
	//else cout << "ERROR in setRotamerWithCheck...\n bpt " << _bpt << " does not exist." << endl;
	//calculateSidechainDihedralAngles();
	setMoved();
}

void residue::setPolarHRotamer(UInt _rotamerIndex)
{
	setPolarHChi(_rotamerIndex);
	calculatePolarHDihedralAngle();
}

void residue::setPolarHRotamerWithCheck(UInt _rotamerIndex)
{
	setPolarHChi(_rotamerIndex);
	calculatePolarHDihedralAngle();
}

DouVec residue::setRotamerWithCheckTest(const UInt _lib, const UInt _bpt, const UInt _rotamer)
{	//cout << " set rotamer (with checking) at residue level " << endl;
	DouVec theAngles;
	theAngles.resize(0);
	if ( dataBase[itsType].chiDefinitions.size() > _bpt &&
		dataBase[itsType].chiDefinitions[_bpt].size() >= 4)
	{	if (dataBase[itsType].itsRotamerLibs[_lib]->rotamersExist(_bpt,_rotamer))
		{	theAngles = itsSidechainDihedralAngles[_bpt];
			itsSidechainDihedralAngles[_bpt] = dataBase[itsType].itsRotamerLibs[_lib]->getAngles(_bpt,_rotamer);
			UInt chiSize = itsSidechainDihedralAngles[_bpt].size();
			for (UInt i=0;i<chiSize;i++)
			{
				//cout << "setRotamerWithCheckTest invoked successfully." << i << ", " << itsSidechainDihedralAngles[_bpt][i] << endl;
				setChi(_bpt,i,itsSidechainDihedralAngles[_bpt][i]);
			}
		}
		else cout << "ERROR in setRotamerWithCheckTest...\n rotamer " << _bpt << ", " << _rotamer << " does not exist." << endl;
    }
	else cout << "ERROR in setRotamerWithCheckTest...\n bpt " << _bpt << " does not exist." << endl;
	calculateSidechainDihedralAngles();
	return theAngles;
	setMoved();
}

void residue::setBetaChi(const double _angle)
{
	if (_angle != 1000.0)
	{
		double currentBetaChi = getBetaChi();
		ASSERT(currentBetaChi < 1e5 && currentBetaChi > -1e5);
		double diff = _angle - currentBetaChi;
		if (itsType == 19){rotate(0,4, diff);}//proline
		else{rotate(0,1, diff);}
		setMoved();
	}
}

void residue::setChi(const UInt _bpt, const UInt _index, const double _angle)
{	double currentChi = getChi(_bpt,_index);
	//cout << "cc:" << currentChi << " ";
	ASSERT(currentChi < 1e5 && currentChi > -1e5);
	double diff = _angle - currentChi;
	setChiByDelta(_bpt, _index, diff);
	setMoved();
}

void residue::setChiByDelta(const UInt _bpt, const UInt _index, const double _angleDelta)
{	if( _bpt < dataBase[itsType].chiDefinitions.size() &&
		_index <= dataBase[itsType].chiDefinitions[_bpt].size()-3)
	{	rotate( dataBase[itsType].chiDefinitions[_bpt][_index+1],
				dataBase[itsType].chiDefinitions[_bpt][_index+2],
				_angleDelta);
	}
}

void residue::setChi(const UInt _index, const double _angle)
{	UInt bpt = 0;
	setChi(bpt, _index, _angle);
}

void residue::setPolarHChi(const UInt _rotamerIndex)
{
	vector<UInt> atoms = dataBase[itsType].getAtomsOfPolarHChi();
	double originalChiAngle = calculateDihedral(atoms);
	double newChiAngle;
	if(atoms[0])
	{
		switch(_rotamerIndex)
		{	case 1:	newChiAngle = 60.0;
				break;
			case 2: newChiAngle = 180.0;
				break;
			case 3: newChiAngle = -60.0;
				break;
			default: {cout << "hit default in residue::setPolarHChi.\n"; newChiAngle=0;}
		}
		if(newChiAngle!=originalChiAngle)
		{	double diff = newChiAngle - originalChiAngle;
			setPolarHChiByDelta(atoms[1],atoms[2],diff);
		}
	}
}

void residue::setPolarHChiByDelta(const UInt _atom1, const UInt _atom2, const double _angle)
{
	rotate(_atom1,_atom2,_angle);
}

double residue::getPolarHChi() const
{
	vector<UInt> theAtomIndices = dataBase[itsType].getAtomsOfPolarHChi();
	return calculateDihedral(theAtomIndices);
	return 0;
}


double residue::getChi(const UInt _bpt, const UInt _index) const
{       
	vector<UInt> theAtomIndices = dataBase[itsType].getAtomsOfChi(_bpt, _index);
	double theAngle = calculateDihedral(theAtomIndices);
	return theAngle;
}

double residue::getChi(const UInt _index) const
{
	UInt bpt = 0;
	return getChi(bpt, _index);
}

double residue::getBetaChi()
{	
	if(itsAtoms.size() > 4)
	{
		if(pItsPrevRes != 0 && (itsAtoms[4]->getName() == "CB" || itsAtoms[4]->getName() == "CD"))
		{
			vector< dblVec > quadVect(4);
			quadVect[0] = pItsPrevRes->getMainChain(2)->getCoords();
			quadVect[1] = itsAtoms[0]->getCoords();
			quadVect[2] = itsAtoms[1]->getCoords();
			if (itsType == 19) {quadVect[3] = itsAtoms[6]->getCoords();} //proline beta carbon in reverse bond order
			else{quadVect[3] = itsAtoms[4]->getCoords();}
			return CMath::dihedral(quadVect[0], quadVect[1], quadVect[2], quadVect[3]);
		}	
	}
	return 1000.0;
}

double residue::getBetaChiR()
{	
	if(itsAtoms.size() > 4)
	{
		if(pItsNextRes != 0 && (itsAtoms[4]->getName() == "CB" || itsAtoms[4]->getName() == "CD"))
		{
			vector< dblVec > quadVect(4);
			quadVect[0] = pItsNextRes->getMainChain(0)->getCoords();
			quadVect[1] = itsAtoms[2]->getCoords();
			quadVect[2] = itsAtoms[1]->getCoords();
			if (itsType == 19) {quadVect[3] = itsAtoms[6]->getCoords();} //proline beta carbon in reverse bond order
			else{quadVect[3] = itsAtoms[4]->getCoords();}
			return CMath::dihedral(quadVect[0], quadVect[1], quadVect[2], quadVect[3]);
		}	
	}
	return 1000.0;
}

double residue::calculateDihedral(const vector<UInt>& _quad) const
{	if( _quad.size() != 4)
	{	return 1000.0;
	}
	vector< dblVec > quadVect;
	for(UInt i=0; i < 4; i++)
	{
		dblVec coords = itsAtoms[_quad[i]]->getCoords();
		quadVect.push_back(coords);
	}
#ifdef _RESIDUE_GEOMETRY_DEBUG
	cout << endl << "index        coords" << endl;
	for (UInt i=0; i<4; i++)
	{	cout << _quad[i] << " " << quadVect[i][0] << " ";
		cout << quadVect[i][1] << " " << quadVect[i][2] << endl;
	}
#endif

	return CMath::dihedral(quadVect[0], quadVect[1], quadVect[2], quadVect[3]);
}

double residue::calculateDihedral(vector<atom*>& _quad) const
{	if( _quad.size() != 4)
	{	return 1000.0;
	}
	vector<dblVec> quadVect;
	quadVect.resize(4);
	for(UInt i=0; i<quadVect.size(); i++)
	{
		quadVect[i].newsize(3);
		quadVect[i] = (_quad[i])->getCoords();
	}
	//cout << endl << "residue::calculateDihedral(vector(atom*)) called" << endl;
	return CMath::dihedral(quadVect[0], quadVect[1], quadVect[2], quadVect[3]);
}

vector <double> residue::getBackboneAngles()
{
	vector <double> angles;
	angles.push_back(getPhi());
	angles.push_back(getPsi());
	return angles;
}
		
double residue::getPhi()
{
	double tempdouble;
	if (pItsPrevRes != 0)
	{
		vector<atom*> fourAtomPointers;
		UInt i = dataBase[pItsPrevRes->getTypeIndex()].mainChain.size()-1;
		fourAtomPointers.push_back(pItsPrevRes->getMainChain(i-1));
		fourAtomPointers.push_back(getMainChain(0));
		fourAtomPointers.push_back(getMainChain(1));
		fourAtomPointers.push_back(getMainChain(2));
		tempdouble = calculateDihedral(fourAtomPointers);
	}
	else
	{
		//cout << "Cannot calcualate PHI for this amino acid" << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

double residue::getPsi()
{
	double tempdouble;
	if (pItsNextRes != 0)
	{
		vector<atom*> fourAtomPointers;
		UInt i = dataBase[itsType].mainChain.size()-1;
		fourAtomPointers.push_back(getMainChain(i-3));
		fourAtomPointers.push_back(getMainChain(i-2));
		fourAtomPointers.push_back(getMainChain(i-1));
		fourAtomPointers.push_back(pItsNextRes->getMainChain(0));
		tempdouble = calculateDihedral(fourAtomPointers);
	}
	else
	{
		//cout << "Cannot calcualate PSI for this amino acid" << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

double residue::getAngle(UInt angleType)
{
    double tempdouble = 0.0;
	if (angleType == 0) //phi
	{	
		if (pItsPrevRes != 0)
		{
			vector<atom*> fourAtomPointers;
			UInt i = dataBase[pItsPrevRes->getTypeIndex()].mainChain.size()-1;
			fourAtomPointers.push_back(pItsPrevRes->getMainChain(i-1));
			fourAtomPointers.push_back(getMainChain(0));
			fourAtomPointers.push_back(getMainChain(1));
			fourAtomPointers.push_back(getMainChain(2));
			tempdouble = calculateDihedral(fourAtomPointers);
		}
		else
		{
			//cout << "Cannot calcualate PHI for this amino acid" << endl;
			tempdouble = 1000.0;
		}
	}
	if (angleType == 1) //psi
	{
		if (pItsNextRes != 0)
		{
			vector<atom*> fourAtomPointers;
			UInt i = dataBase[itsType].mainChain.size()-1;
			fourAtomPointers.push_back(getMainChain(i-3));
			fourAtomPointers.push_back(getMainChain(i-2));
			fourAtomPointers.push_back(getMainChain(i-1));
			fourAtomPointers.push_back(pItsNextRes->getMainChain(0));
			tempdouble = calculateDihedral(fourAtomPointers);
		}
		else
		{
			//cout << "Cannot calcualate PSI for this amino acid" << endl;
			tempdouble = 1000.0;
		}
	}
	return tempdouble;
}

int residue::setPhi(double _phi)
{
	if (pItsPrevRes != 0)
	{
		double currentPhi = getPhi();
		ASSERT(currentPhi < 1e5 && currentPhi > -1e5);
		double angle = _phi - currentPhi;
		rotate(getMainChain(0), getMainChain(1), angle, true);
		//cout << "Phi " << _phi << " old " << currentPhi << " delta " << angle << " current " << getPhi() << endl;
	}
	else
	{
		//cout << "Cannot set PHI for the first amino acid in a chain." << endl;
		return -1;
	}
	return 0;
}

int residue::setPsi(double _psi)
{
	if (pItsNextRes != 0)
	{
		double currentPsi = getPsi();
		ASSERT (currentPsi < 1e5 && currentPsi > -1e5);
		double angle = _psi - currentPsi;
		UInt i = dataBase[itsType].mainChain.size()-1;
		rotate(getMainChain(i-2), getMainChain(i-1), angle, true);
		//cout << "Psi " << _psi << " old " << currentPsi << " delta " << angle << " current " << getPsi() << endl;
	}
	else
	{
		/*double currentPsi = getPsi();
		ASSERT (currentPsi < 1e5 && currentPsi > -1e5);
		double angle = _psi - currentPsi;
		UInt i = dataBase[itsType].mainChain.size()-1;
		rotate(getMainChain(i-2), getMainChain(i-1), 180-angle, false);*/
		return -1;
	}
	return 0;
}

int residue::setDihedral(double _dihedral, UInt _angleType, UInt _direction)
{
	if (_angleType == 0) //phi
	{
		if (pItsPrevRes != 0)
		{
			double currentPhi = getPhi();
			ASSERT(currentPhi < 1e5 && currentPhi > -1e5);
			double deltaTheta = _dihedral - currentPhi;
			if (_direction == 0)
			{
				rotateDihedral(getMainChain(0), getMainChain(1), deltaTheta, _angleType, _direction); //NtoC
			}
			if (_direction == 1)
			{
				rotateDihedral(getMainChain(1), getMainChain(0), deltaTheta, _angleType, _direction); //CtoN
			}
		}
		else
		{
			//cout << "Cannot set PHI for the first amino acid in a chain." << endl;
			return -1;
		}
	}
	if (_angleType == 1) //psi
	{
		if (pItsNextRes != 0)
		{
			double currentPsi = getPsi();
			ASSERT (currentPsi < 1e5 && currentPsi > -1e5);
			double deltaTheta = _dihedral - currentPsi;
			UInt i = dataBase[itsType].mainChain.size()-1;
			if (_direction == 0)
			{
				rotateDihedral(getMainChain(i-2), getMainChain(i-1), deltaTheta, _angleType, _direction); //NtoC
			}
			if (_direction == 1)				
			{
				rotateDihedral(getMainChain(i-1), getMainChain(i-2), deltaTheta, _angleType, _direction); //CtoN
			}
		}
		else
		{
			//cout << "Cannot set PSI for the last amino acid in a chain." << endl;
			return -1;
		}
	}
	return 0;
}

double residue::getOmega()
{
	double tempdouble;
	if (pItsNextRes != 0)
	{
		UInt i = dataBase[itsType].mainChain.size()-1;
		vector<atom*> fourAtomPointers;
		fourAtomPointers.push_back(getMainChain(i-4));
		fourAtomPointers.push_back(getMainChain(i-3));
		fourAtomPointers.push_back(getMainChain(i-2));
		fourAtomPointers.push_back(getMainChain(i-1));
		tempdouble = calculateDihedral(fourAtomPointers);
	}
	else
	{
		//cout << "Cannot calcualate OMEGA for this amino acid" << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

void residue::setOmega(double _omega)
{
	if (pItsNextRes != 0)
	{
		UInt i = dataBase[itsType].mainChain.size()-1;
		double currentOmega = getOmega();
		ASSERT(currentOmega < 1e5 && currentOmega > -1e5);
		double angle = _omega - currentOmega;
		rotate(getMainChain(i-3), getMainChain(i-2), angle, false);
	}
	return;
}

double residue::getAmide()
{
	double tempdouble;
	if (pItsNextRes != 0)
	{
		vector<atom*> fourAtomPointers;
		UInt  i = dataBase[itsType].mainChain.size()-1;
		fourAtomPointers.push_back(getMainChain(i-2));
		fourAtomPointers.push_back(getMainChain(i-1));
		fourAtomPointers.push_back(pItsNextRes->getMainChain(0));
		fourAtomPointers.push_back(pItsNextRes->getMainChain(1));
		tempdouble = calculateDihedral(fourAtomPointers);
	}
	else
	{
		cout << "Cannot calcualate AMIDE for this amino acid" << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

atom* residue::getMainChain(UInt _index)
{
	return itsAtoms[dataBase[itsType].mainChain[_index]];
}

UInt residue::getNumBpt(const UInt _resTypeIndex)
{
	UInt theNumber = 0;
	if (_resTypeIndex < getDataBaseSize() )
	{
		theNumber =  dataBase[_resTypeIndex].branchPoints.size();
	}
	return theNumber;
}

UInt residue::getNumDihedralAngles(const UInt _resTypeIndex, const UInt _bpt)
{
	int chisDefined = dataBase[_resTypeIndex].getNumberOfChis(_bpt);
	if (chisDefined >=0)
	{
		return chisDefined;
	}
	return 0;
}

void residue::calculateSidechainDihedralAngles()
{
	UInt branchpoints = getNumBpt(itsType);
	for (UInt i=0; i<branchpoints; i++)
	{	UInt dihedrals = getNumDihedralAngles(itsType,i);
		double angle;
		for (UInt j=0; j<dihedrals; j++)
		{	angle = getChi(i,j);
			itsSidechainDihedralAngles[i][j] = angle;
		}
	}
}

vector< vector< double > > residue::randContinuousSidechainConformation()
{	
	vector < vector <double> > sideChainDihedralAngles;
	vector <double> chis;
	double angle;
	UInt branchpoints = getNumBpt(itsType);
	for (UInt i=0; i<branchpoints; i++)
	{	
		chis.clear();
		UInt dihedrals = getNumDihedralAngles(itsType,i);
		for (UInt j=0; j<dihedrals; j++)
		{
			angle = (rand() % 360)-179;
			chis.push_back(angle);
		}
		sideChainDihedralAngles.push_back(chis);
	}
	return sideChainDihedralAngles;
}

void residue::calculatePolarHDihedralAngle()
{
	if(dataBase[itsType].getHasPolarHRotamers())
	{  itsPolarHDihedralAngle = getPolarHChi();  }
}

vector< vector< double > > residue::getSidechainDihedralAngles()
{
	calculateSidechainDihedralAngles();
	 return itsSidechainDihedralAngles;
}

// ************************************************************************
// ************************************************************************
// Rotation and Translation and more .....
// ************************************************************************
// ************************************************************************

dblVec residue::getCoords(const string _atomName)
{
    dblVec coords(3);
    for (UInt i = 0; i < itsAtoms.size(); i++)
    {
        //cout << itsAtoms[i]->getName() << " ";
        string itsName = itsAtoms[i]->getName();
        if (itsName == _atomName)
        {
            coords[0] = itsAtoms[i]->getX();
            coords[1] = itsAtoms[i]->getY();
            coords[2] = itsAtoms[i]->getZ();
        }
    }
    if (coords.size() != 3) cout << "ATOM " << _atomName << " not found ..." << endl;
    return coords;
}

void residue::printCoords() const
{	for (UInt i=0; i<itsAtoms.size(); i++)
	{	cout << itsAtoms[i]->getName() << ' ';
		cout << itsAtoms[i]->getX() << ' ';
		cout << itsAtoms[i]->getY() << ' ';
		cout << itsAtoms[i]->getZ() << endl;
	}
}

void residue::rotate(UInt _first, UInt _second, double _theta)
{	bool backboneRotation = false;
	atom* pAtom1 = itsAtoms[_first];
	atom* pAtom2 = itsAtoms[_second];
	if (dataBase[itsType].isMainChain[_first] &&
	    dataBase[itsType].isMainChain[_second])
	{	backboneRotation = true;
	}
	rotate(pAtom1, pAtom2, _theta, backboneRotation);
	setMoved();
}

void residue::rotateDihedral(atom* _pAtom1, atom* _pAtom2, double _deltaTheta,  UInt _angleType, UInt _direction)
{
	setMoved();
	dblVec toOrigin = _pAtom1->getCoords() * (-1.0);
	dblVec backHome = _pAtom1->getCoords();

	if (_direction == 0)
	{	
		_pAtom1->translate(toOrigin);
		_pAtom2->translate(toOrigin);
		_pAtom2->translateChildren(toOrigin);
		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslateWithDirection(toOrigin, _direction);
		}
	}
	if (_direction == 1)
	{
		recursiveTranslateWithDirection(toOrigin, _direction);
	}
	
	dblVec atomCoords = _pAtom2->getCoords();
	dblMat R(3,3,0.0);
	R = CMath::rotationMatrix(atomCoords, _deltaTheta);
	
	if (_direction == 0)
	{	_pAtom2->transformChildren(R);
		if (pItsNextRes)
		{	
			pItsNextRes->recursiveTransform(R);
		}
	}
	if (_direction == 1)
	{
		if (_angleType == 0){
			if (pItsPrevRes)
			{	
				pItsPrevRes->recursiveTransformR(R);
			}
		}
		else{
			recursiveTransformR(R);
			R = CMath::rotationMatrix(atomCoords, _deltaTheta*-1);
			_pAtom1->transformChildren(R);
		}
	}

	if (_direction == 0)
	{	
		_pAtom1->translate(backHome);
		_pAtom2->translate(backHome);
		_pAtom2->translateChildren(backHome);
		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslateWithDirection(backHome, _direction);
		}
	}
	if (_direction == 1)
	{
		recursiveTranslateWithDirection(backHome, _direction);
	}
	if(_angleType == 0 && _direction == 1)
	{alignAmideProtonToBackbone();}
	setMoved();
}

void residue::rotate(atom* _pAtom1, atom* _pAtom2, double _theta,
					 bool backboneRotation)
{	// the coords of _atom1 give us the translation vector to get
	// _atom1 to the origin
#ifdef __RES_DEBUG
	cout << "original coords" << endl;
	cout << _pAtom1->getName() << " " << _pAtom1->getCoords() << endl;
	cout << _pAtom2->getName() << " " << _pAtom2->getCoords() << endl;
	_pAtom2->queryChildrensCoords();
#endif

	dblVec toOrigin = _pAtom1->getCoords() * (-1.0);
	dblVec backHome = _pAtom1->getCoords();

	// Now use the toOrigin translation to translate
	// _atom1, _atom2, and all children of _atom2

	_pAtom1->translate(toOrigin);
	_pAtom2->translate(toOrigin);
	_pAtom2->translateChildren(toOrigin);
	if (backboneRotation)
	{	//cout << "Backbone rotation!" << endl;
		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslate(toOrigin);
		}
	}

#ifdef __RES_DEBUG
	cout << "after translation to origin" << endl;
	cout << _pAtom1->getName() << " " << _pAtom1->getCoords() << endl;
	cout << _pAtom2->getName() << " " << _pAtom2->getCoords() << endl;
	_pAtom2->queryChildrensCoords();
#endif
	//calcualate rotation matrix based on vector defined by
	//point q1 (second point)

	dblMat R(3,3,0.0);
	R = CMath::rotationMatrix(_pAtom2->getCoords(), _theta);

#ifdef __RES_DEBUG
	cout << "rotation matrix" << endl;
	cout << R;
#endif
	_pAtom2->transformChildren(R);

	if (backboneRotation)
	{	//some code to perform rotation on other residues
		if (pItsNextRes)
		{	pItsNextRes->recursiveTransform(R);
		}
	}

#ifdef __RES_DEBUG
	cout << "after rotation has been applied" << endl;
	_pAtom2->queryChildrensCoords();
#endif

	//Finally, translate all the atoms back via the backHome vector
	_pAtom1->translate(backHome);
	_pAtom2->translate(backHome);
	_pAtom2->translateChildren(backHome);

	if (backboneRotation)
	{	//some code to perform rotation on other residues
		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslate(backHome);
		}
	}

#ifdef __RES_DEBUG
	cout << "after translation back Home" << endl;
	cout << _pAtom1->getName() << " " << _pAtom1->getCoords() << endl;
	cout << _pAtom2->getName() << " " << _pAtom2->getCoords() << endl;
	_pAtom2->queryChildrensCoords();
#endif
	setMoved();
}

void residue::rotate(const point& _point, const dblMat& _RMatrix )
{	atom* pAtom = itsAtoms[0];
	// the coords of _pPoint give us the translation vector to get
	// _pPoint to the origin

	dblVec toOrigin = _point.getCoords() * (-1.0);
	dblVec backHome = _point.getCoords();

	// Now use the toOrigin translation to translate
	// _atom1, _atom2, and all children of _atom2

	pAtom->translate(toOrigin);
	pAtom->translateChildren(toOrigin);

	pAtom->transform(_RMatrix);
	pAtom->transformChildren(_RMatrix);

#ifdef __RES_DEBUG
	cout << "after rotation has been applied" << endl;
	pAtom->queryChildrensCoords();
#endif

//Finally, translate all the atoms back via the backHome vector
	pAtom->translate(backHome);
	pAtom->translateChildren(backHome);

#ifdef __RES_DEBUG
	cout << "after translation back Home" << endl;    
	cout << _pAtom->getName() << " " << _pAtom->getCoords() << endl;
	pAtom->queryChildrensCoords();
#endif
	setMoved();
}

void residue::rotate_new(atom* _pivotAtom, const dblMat& _RMatrix)
{
	dblVec toOrigin = _pivotAtom->getCoords() * (-1.0);
	dblVec backHome =  toOrigin * (-1.0);
	
	_pivotAtom->translate(toOrigin);
	_pivotAtom->translateChildren(toOrigin);
	_pivotAtom->transform(_RMatrix);
	_pivotAtom->transformChildren(_RMatrix);

	_pivotAtom->translate(backHome);
	_pivotAtom->translateChildren(backHome);
	setMoved();
}

void residue::rotate_new(atom* _pivotAtom, atom* _firstAtom, const dblMat& _RMatrix)
{
	dblVec toOrigin = _pivotAtom->getCoords() * (-1.0);
	dblVec backHome =  toOrigin * (-1.0);
	
	_firstAtom->translate(toOrigin);
	_firstAtom->translateChildren(toOrigin);
	_firstAtom->transform(_RMatrix);
	_firstAtom->transformChildren(_RMatrix);

	_firstAtom->translate(backHome);
	_firstAtom->translateChildren(backHome);
	setMoved();
}


void residue::rotate(const point& _point, const dblVec& _R_axis,const double _theta)
{	atom* pAtom = itsAtoms[0];
	// the coords of _pPoint give us the translation vector to get
	// _pPoint to the origin

	dblVec toOrigin = _point.getCoords() * (-1.0);
	dblVec backHome = _point.getCoords();

	// Now use the toOrigin translation to translate
	// _atom1, _atom2, and all children of _atom2

	pAtom->translate(toOrigin);
	pAtom->translateChildren(toOrigin);

	//calculate rotation matrix based on _R_axis vector

	dblMat R(3,3,0.0);
	R = CMath::rotationMatrix(_R_axis, _theta);
#ifdef __RES_DEBUG
	cout << "rotation matrix" << endl;                            
	cout << R;                                                    
#endif
                                                                      
	pAtom->transform(R);
	pAtom->transformChildren(R);

#ifdef __RES_DEBUG                                                    
	cout << "after rotation has been applied" << endl;            
	pAtom->queryChildrensCoords();                               
#endif

	//Finally, translate all the atoms back via the backHome vector
	pAtom->translate(backHome);                                  
	pAtom->translateChildren(backHome);                          

#ifdef __RES_DEBUG                                                    
	cout << "after translation back Home" << endl;                
	cout << _pAtom->getName() << " " << _pAtom->getCoords() << endl;
	pAtom->queryChildrensCoords();                               
#endif                                                                
	setMoved();
}                                       

void residue::rotate(atom* _pAtom, const dblVec& _R_axis, const double _theta)
{	// the coords of _atom1 give us the translation vector to get
	// _atom1 to the origin

#ifdef __RES_DEBUG
	cout << "original coords" << endl;
	cout << _pAtom->getName() << " " << _pAtom->getCoords() << endl;
	_pAtom->queryChildrensCoords();
#endif

	dblVec toOrigin = _pAtom->getCoords() * (-1.0);
	dblVec backHome = _pAtom->getCoords();

	// Now use the toOrigin translation to translate
	// _atom1, _atom2, and all children of _atom2

	_pAtom->translate(toOrigin);
	_pAtom->translateChildren(toOrigin);

	//calculate rotation matrix based on vector defined by
	//point q1 (second point)

	dblMat R(3,3,0.0);
	R = CMath::rotationMatrix(_R_axis, _theta);

#ifdef __RES_DEBUG
	cout << "rotation matrix" << endl;                            
	cout << R;                                                    
#endif
	_pAtom->transformChildren(R);                                 

#ifdef __RES_DEBUG                                                    
	cout << "after rotation has been applied" << endl;            
	_pAtom->queryChildrensCoords();                               
#endif                                                                

	//Finally, translate all the atoms back via the backHome vector
	_pAtom->translate(backHome);                                  
	_pAtom->translateChildren(backHome);                          

#ifdef __RES_DEBUG                                                    
	cout << "after translation back Home" << endl;                
	cout << _pAtom->getName() << " " << _pAtom->getCoords() << endl;
	_pAtom->queryChildrensCoords();
#endif                                                                
	setMoved();
} 
                                                    
void residue::translate(const dblVec& _dblVec)
{
    for (UInt i=0; i < itsAtoms.size(); i++)
	{	itsAtoms[i]->translate(_dblVec);
	}
	setMoved();
}


void residue::recursiveTranslateWithDirection(dblVec& _dblVec, UInt _direction)
{	
    setMoved();
	translate(_dblVec);
	if (_direction == 0)
	{
		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslateWithDirection(_dblVec, _direction);
		}
	}
	if (_direction == 1)
	{
		if (pItsPrevRes)
		{	pItsPrevRes->recursiveTranslateWithDirection(_dblVec, _direction);
		}
	}
}

void residue::recursiveTranslate(dblVec& _dblVec)
{
    setMoved();
	translate(_dblVec);
	if (pItsNextRes)
	{	pItsNextRes->recursiveTranslate(_dblVec);
	}
}

void residue::recursiveTransform(dblMat& _dblMat)
{
    setMoved();
	transform(_dblMat);
	if (pItsNextRes)
	{	pItsNextRes->recursiveTransform(_dblMat);
	}
}

void residue::recursiveTransformR(dblMat& _dblMat)
{
    setMoved();
	transform(_dblMat);
	if (pItsPrevRes)
	{	pItsPrevRes->recursiveTransformR(_dblMat);
	}
	setMoved();
}

void residue::transform(const dblMat& _dblMat)
{	for (UInt i=0; i < itsAtoms.size(); i++)
	{	itsAtoms[i]->transform(_dblMat);
	}
    setMoved();
}



void residue::printMainChain() const
{	if(itsType < dataBase.size())
	{ 	cout << " main chain atoms: ";
		for(UInt i=0; i< dataBase[itsType].mainChain.size(); i++)
		{	cout << itsAtoms[ dataBase[itsType].mainChain[i] ]->getName();
			cout << " ";
		}
		cout << endl;
	}
	else
	{	cout << "Error: itsType incompatible with dataBase " << endl;
		cout << "Error reported by residue::printMainChain() " << endl;
	}
}

void residue::printBranchPoints() const
{	if(itsType < dataBase.size())
	{	cout << " main chain branch points: ";
		for(UInt i=0; i< dataBase[itsType].branchPoints.size(); i++)
		{	cout << itsAtoms[dataBase[itsType].branchPoints[i]]->getName();
			cout << ' ';
		}
		cout << endl;
	}
	else
	{	cout << "Error: itsType incompatible with dataBase " << endl;
		cout << "Error reported by residue::printBranchPoints() " << endl;
	}
}


void residue::makeAtomSilent(const UInt _atomIndex)
{
	if (_atomIndex >= 0 && _atomIndex < itsAtoms.size())
	{
		itsAtoms[_atomIndex]->makeAtomSilent();
	}
	else
		cout << "ERROR in residue::makeAtomSilent( ... )\n\t" << _atomIndex << " is out of range as an atom number..." << endl;
	return;
}

void residue::makeResidueSilent()
{
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		itsAtoms[i]->makeAtomSilent();
	}
	return;
}

double residue::getIntraEnergy(const UInt _atom1, residue* _other, const UInt _atom2)
{
	if (_atom1 >=0 && _atom1 < itsAtoms.size() && _atom2 >=0 && _atom2 < _other->itsAtoms.size())
	{
		double vdwEnergy = 0.0;
		double elecEnergy = 0.0;
		double intraEnergy = 0.0;
		bool isBonded = false;
		if (this == _other)
			isBonded = isSeparatedByFewBonds(_atom1, _atom2);
		else
			isBonded = isSeparatedByFewBonds(this, _atom1, _other, _atom2);
		if (!isBonded)
		{
			double distanceSquared = itsAtoms[_atom1]->distanceSquared(_other->itsAtoms[_atom2]);
			cout << "distance " <<sqrt(distanceSquared) << endl;
			// calculate VDW Energy
			int index1, index2;
			if (hydrogensOn)
			{
				index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[_atom1][0];
				index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[_atom2][0];
			}
			else
			{
				index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[_atom1][1];
				index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[_atom2][1];
			}
			cout << "indices " << index1 << " " << index2 << endl;
			vdwEnergy = residueTemplate::getVDWEnergySQ(index1,index2,distanceSquared);
			// calculate Elec Energy
			UInt resType1 = itsType;
			UInt resType2 = _other->itsType;
			elecEnergy = residueTemplate::getAmberElecEnergySQ(resType1, _atom1, resType2, _atom2, distanceSquared);
		}
		cout << "ELEC energy:  " << elecEnergy << endl;
		cout << "VDW energy:  " << vdwEnergy << endl;
		intraEnergy = elecEnergy + vdwEnergy;
		cout << "INTRA energy:  " << intraEnergy << endl;
		return intraEnergy;
	}
	else
	{
		cout << "ERROR i}n residue::getIntraEnergy(...)" << endl;
		cout << "atom index or indices out of range." << endl;
		exit(1);
	}
}

double residue::intraEnergy()
{	//double distance;
	double distanceSquared;
	int index1;
	int index2;
	double intraEnergy = 0.0;
	double vdwEnergy = 0.0;
	double pmfEnergy = 0.0;
	double amberElecEnergy = 0.0;
    bool twoBonds;

	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			for(UInt j=i+1; j<itsAtoms.size(); j++)
			{
				if (!itsAtoms[j]->getSilentStatus())
				{
					//distance = itsAtoms[i]->distance(itsAtoms[j]);
					distanceSquared = itsAtoms[i]->distanceSquared(itsAtoms[j]);
					twoBonds = isSeparatedByOneOrTwoBonds(i,j);
					// ** intra AMBER vdW
					if (residueTemplate::itsAmberVDW.getScaleFactor() != 0.0 )
					{
						if (hydrogensOn)
						{
								index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
								index2 = dataBase[itsType].itsAtomEnergyTypeDefinitions[j][0];
						}
						else
						{
								index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][1];
								index2 = dataBase[itsType].itsAtomEnergyTypeDefinitions[j][1];
						}
						if (!twoBonds)
						{
							double tempvdwEnergy = residueTemplate::getVDWEnergySQ(index1,index2,distanceSquared);
							vdwEnergy += tempvdwEnergy;
						}
					}

					// ** intra AMBER Electrostatics
					if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0)
					{
						if (!twoBonds)
						{
							UInt resType1 = itsType;
							UInt atomType1 = i;
							UInt resType2 = itsType;
							UInt atomType2 = j;
							double tempAmberElecEnergy = residueTemplate::getAmberElecEnergySQ(resType1, atomType1, resType2, atomType2, distanceSquared);
							amberElecEnergy += tempAmberElecEnergy;
						}
					}
				}
			}
		}
	}
//	if (vdwEnergy > 1.0e4) cout << "intra: " << itsType << vdwEnergy << endl;
	intraEnergy = pmfEnergy + vdwEnergy + amberElecEnergy;
	if (intraEnergy < -1e5)
	{
		cout << "Extremely Low intraEnergy value of " << intraEnergy << endl;
	}
	if (aaBaseline::getScaleFactor() != 0.0)
	{
		intraEnergy += residueTemplate::getAABaselineEnergy(getType());
	}
	return intraEnergy;
}

// protEnergy functions (new (2019) complete solvation and forcefeild)-----------------------------------------------

double residue::intraSoluteEnergy()
{	
	double intraEnergy = 0.0;
	bool threeBonds;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			// ** get solvationEnergy
			if (EsolvationFactor != 0.0 || HsolvationFactor != 0.0){
				double tempSolvEnergy = calculateSolvationEnergy(i);
				intraEnergy += tempSolvEnergy;
			}
			for(UInt j=i+1; j<itsAtoms.size(); j++)
			{
				if (!itsAtoms[j]->getSilentStatus())
				{
					threeBonds = isSeparatedByFewBonds(i,j);
					if (!threeBonds)
					{
						// ** get distance
						double distanceSquared = itsAtoms[i]->distanceSquared(itsAtoms[j]);

						// ** intra AMBER vdW
						if (residueTemplate::itsAmberVDW.getScaleFactor() != 0.0)
						{
							int index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
							int index2 = dataBase[itsType].itsAtomEnergyTypeDefinitions[j][0];
							double tempvdwEnergy = residueTemplate::getVDWEnergySQ(index1,index2,distanceSquared);
							intraEnergy += tempvdwEnergy;
						}

						// ** intra Electrostatics
						if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0)
						{
							// ** get dielectric average of atoms
							double dielectric = (itsAtoms[i]->getDielectric() + itsAtoms[j]->getDielectric()) * 0.5;
							
							//recalculate the dielectric using the Maxwell Garnett mixing formula to include the polarizability of the pairwise dipole inclusion
							if (polarizableElec){
								dielectric = maxwellGarnettApproximation(i, j, dielectric, distanceSquared);
							}
							
							// calculate coulombic energy with effective dielectric
							double tempAmberElecEnergy = residueTemplate::getAmberElecSoluteEnergySQ(itsType, i, itsType, j, distanceSquared, dielectric);
							intraEnergy += tempAmberElecEnergy;
						}
					}
				}
			}
		}
	}
	return intraEnergy;
}

double residue::interSoluteEnergy(residue* _other)
{
	double interEnergy = 0.0;
	bool threeBonds;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			for(UInt j=0; j<_other->itsAtoms.size(); j++)
			{
				if (!_other->itsAtoms[j]->getSilentStatus())
				{
					threeBonds = isSeparatedByThreeBackboneBonds(i,_other,j);
					if (!threeBonds)
					{
						double distanceSquared = itsAtoms[i]->inCubeWithDistSQ(_other->itsAtoms[j], cutoffDistance);
						if (distanceSquared <= cutoffDistanceSquared)
						{
							// ** inter AMBER Electrostatics
							if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0)
							{
								// ** get dielectric average
								double dielectric = (itsAtoms[i]->getDielectric() + _other->itsAtoms[j]->getDielectric()) * 0.5;
								
								//recalculate the dielectric using the Maxwell Garnett mixing formula to include the polarizability of the pairwise dipole-dipole inclusion of hbonds
								if (polarizableElec){
									dielectric = maxwellGarnettApproximation(i, _other, j, dielectric, distanceSquared);
								}
								// calculate coulombic energy with effective dielectric
								double tempAmberElecEnergy = residueTemplate::getAmberElecSoluteEnergySQ(itsType, i, _other->itsType, j, distanceSquared, dielectric);
								interEnergy += tempAmberElecEnergy;
							}
							// ** inter AMBER vdW
							if (residueTemplate::itsAmberVDW.getScaleFactor() != 0.0)
							{
								int index1, index2;
								index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
								index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][0];
								double tempvdwEnergy = residueTemplate::getVDWEnergySQ(index1, index2, distanceSquared);
								interEnergy += tempvdwEnergy;
							}
						}
					}
				}
			}
		}
	}
	return interEnergy;
}

double residue::calculateSolvationEnergy(UInt _atomIndex)
{	// note: Requires update of dielectrics at protein level to be accurate for water count and local dielctric. Meant to be part of protEnergy().
	double soluteSolventEnthalpy = 0.0;
	double soluteSolventEntropy = 0.0;

	// First estimate water occupancy around solute atom in solvent volume shells of total proximal solute atom excluded volume
	int atomVDWtype = dataBase[itsType].itsAtomEnergyTypeDefinitions[_atomIndex][0];
	double solvatedRadius = residueTemplate::getVDWRadius(atomVDWtype)+1.4; //atom radius + water radius
	int shellWaters = itsAtoms[_atomIndex]->getNumberofWaters();;
	if (shellWaters > 0){
		// Polar solvation
		if (EsolvationFactor != 0.0)
		{	// Electrostatic enthalpy estimate of solute atom and solvent
			// with variable dielectric and water occupancy estimate
			// Born Electrostatic solvation Still WC, et al J Am Chem Soc 1990
			double atomDielectric = itsAtoms[_atomIndex]->getDielectric();
			double charge = residueTemplate::itsAmberElec.getItsCharge(itsType, _atomIndex);
			soluteSolventEnthalpy += (-(KC/2)*(charge*charge)/(solvatedRadius*atomDielectric))*shellWaters*EsolvationFactor;
		}
	
		// Non-Polar solvation
		if (HsolvationFactor != 0.0)
		{	// Lennard Jones dipole packing ethalpy estimate assuming ideal interaction of solute atom and solvent
			// TIP3P VDW water interaction R. W. Impey, and M. L. Klein, J. Chem. Phys. 79 (1983) 926-935
			double tempvdwEnergy = residueTemplate::getVDWWaterEnergy(atomVDWtype);
			soluteSolventEnthalpy += tempvdwEnergy*shellWaters;
	
			// Solvent Entropy loss estimate due to lack of ideal water lattice hydrogen bond network formation (hydrophobic effect)
			// Gill Hydrophobic solvation  S.J.Gill, S.F.Dec. J Phys. Chem. 1985
			soluteSolventEntropy = (-temperature*KB*log(pow(0.5,shellWaters)))*HsolvationFactor;
		}
	}
	//Total atom solvation Energy
	double solvationEnergy = soluteSolventEntropy+soluteSolventEnthalpy;
	itsAtoms[_atomIndex]->setSolvationEnergy(solvationEnergy);
	return solvationEnergy;
}

double residue::getSolvationEnergy()
{
	double solvationEnergy = 0.0;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		solvationEnergy += itsAtoms[i]->getSolvationEnergy();
	}
	return solvationEnergy;
}

double residue::getDielectric()
{
	double dielectricTotal = 0.0;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		dielectricTotal += itsAtoms[i]->getDielectric();
	}
	return dielectricTotal/itsAtoms.size();
}

void residue::polarizability()
{	
	bool inCube;
	int vdwIndexI, vdwIndexJ;
	double solvatedRadius;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			//--inlude self volume
			vdwIndexI = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
			solvatedRadius = residueTemplate::getVDWRadius(vdwIndexI)+dielectricWidth;
			itsAtoms[i]->sumEnvVol(residueTemplate::getVolume(vdwIndexI)/2);
			for(UInt j=i+1; j<itsAtoms.size(); j++)
			{
				if (!itsAtoms[j]->getSilentStatus())
				{
					inCube = itsAtoms[i]->inCube(itsAtoms[j], solvatedRadius);
					if (inCube)
					{
						//i sum environment j
						vdwIndexJ = dataBase[itsType].itsAtomEnergyTypeDefinitions[j][0];
						itsAtoms[i]->sumEnvVol(residueTemplate::getVolume(vdwIndexJ)/2);
						
						//j sum environment i
						itsAtoms[j]->sumEnvVol(residueTemplate::getVolume(vdwIndexI)/2);
					}
				}
			}
		}
	}
}

void residue::polarizability(residue* _other)
{	
	bool inCube, resI = getMoved(0), resJ = _other->getMoved(0);
	int vdwIndexI, vdwIndexJ;
	double solvatedRadius;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			vdwIndexI = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
			solvatedRadius = residueTemplate::getVDWRadius(vdwIndexI)+dielectricWidth;
			for(UInt j=0; j<_other->itsAtoms.size(); j++)
			{
				if (!_other->itsAtoms[j]->getSilentStatus())
				{
					inCube = itsAtoms[i]->inCube(_other->itsAtoms[j], solvatedRadius);
					if (inCube)
					{
						//i sum environment j
						if (resI){
							vdwIndexJ = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][0];
							itsAtoms[i]->sumEnvVol(residueTemplate::getVolume(vdwIndexJ)/2);
						}
						//j sum environment i
						if (resJ){
							_other->itsAtoms[j]->sumEnvVol(residueTemplate::getVolume(vdwIndexI)/2);
						}
					}
				}
			}
		}
	}
}

void residue::calculateDielectrics()
{
	double envVol, totalWaterVol, dielectric, waters;
	double waterPol = residueTemplate::getPolarizability(52);
	double waterVol = residueTemplate::getVolume(52);
	double pol, vol, solvatedRadius;
	UInt vdwIndexI;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			// calculate local dielectric for atom
			vdwIndexI = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
			solvatedRadius = residueTemplate::getVDWRadius(vdwIndexI)+dielectricWidth;
			vol = pow((solvatedRadius*2),3); pol=0.0; waters=0.0;
			envVol = itsAtoms[i]->getEnvVol();
			totalWaterVol = vol-envVol;
			if (totalWaterVol > 0){
				waters = (totalWaterVol/waterVol); pol = waters*waterPol;
			}

			// Solve for the effective dielectric with the Lorentz local field correction
			dielectric =1+(8*PI/3)*(pol)/1-(4*PI/3)*(pol);
			itsAtoms[i]->setDielectric(dielectric);
			itsAtoms[i]->setNumberofWaters(waters);
		}
	}
}

double residue::maxwellGarnettApproximation(UInt _atomIndex1, UInt _atomIndex2, double _dielectric, double _distanceSquared)
{	//Polarizable electrostatics model via a dipole-dipole polarization effect on the medium
	//Vadim A. Markel 1244 Vol. 33, No. 7 / July 2016 / J Opt Soc Amer
	
	// Check for hbonds and metal ligation where polarization is significant and worth calculating
	if( ((itsAtoms[_atomIndex1]->getType() == "H" || itsAtoms[_atomIndex1]->getType() == "FE" || itsAtoms[_atomIndex1]->getType() == "NI") &&
	     (itsAtoms[_atomIndex2]->getType() == "O" || itsAtoms[_atomIndex2]->getType() == "S"  || itsAtoms[_atomIndex2]->getType() == "N")) ||
		((itsAtoms[_atomIndex2]->getType() == "H" || itsAtoms[_atomIndex2]->getType() == "FE" || itsAtoms[_atomIndex2]->getType() == "NI") &&
	     (itsAtoms[_atomIndex1]->getType() == "O" || itsAtoms[_atomIndex1]->getType() == "S"  || itsAtoms[_atomIndex1]->getType() == "N"))
	  )
	{
		//get dipole-dipole polarization
		double pol = approximateDipoleDipolePolarization(_atomIndex1, _atomIndex2);
		double vol = 4/3*PI*pow((sqrt(_distanceSquared)/2),3);
		
		//recalculate the dielectric using the Maxwell Garnett mixing formula to include the polarizability of the dipole inclusion over the volume of inclusion
		double dielectric = _dielectric+4*PI*(pol/vol)/1-(4*PI/3*_dielectric)*(pol/vol);
		if (dielectric < 1){dielectric = 1;}
		return dielectric;
	}
	return _dielectric;
}

double residue::maxwellGarnettApproximation(UInt _atomIndex1, residue* _other, UInt _atomIndex2,  double _dielectric, double _distanceSquared)
{	//Polarizable electrostatics model via a dipole-dipole polarization effect on the medium
	//Vadim A. Markel 1244 Vol. 33, No. 7 / July 2016 / J Opt Soc Amer
	
	// Check for hbonds and metal ligation where polarization is significant and worth calculating
	if( ((itsAtoms[_atomIndex1]->getType() == "H" || itsAtoms[_atomIndex1]->getType() == "FE" || itsAtoms[_atomIndex1]->getType() == "NI") &&
	     (_other->itsAtoms[_atomIndex2]->getType() == "O" || _other->itsAtoms[_atomIndex2]->getType() == "S"  || _other->itsAtoms[_atomIndex2]->getType() == "N")) ||
		((_other->itsAtoms[_atomIndex2]->getType() == "H" || _other->itsAtoms[_atomIndex2]->getType() == "FE" || _other->itsAtoms[_atomIndex2]->getType() == "NI") &&
	     (itsAtoms[_atomIndex1]->getType() == "O" || itsAtoms[_atomIndex1]->getType() == "S"  || itsAtoms[_atomIndex1]->getType() == "N"))
	  )
	{	//get dipole-dipole polarization
		double pol = approximateDipoleDipolePolarization(_atomIndex1, _other, _atomIndex2);
		double vol = 4/3*PI*pow((sqrt(_distanceSquared)/2),3);
		
		//recalculate the dielectric using the Maxwell Garnett mixing formula to include the polarizability of the dipole inclusion over the volume of inclusion
		double dielectric = _dielectric+4*PI*(pol/vol)/1-(4*PI/3*_dielectric)*(pol/vol);
		if (dielectric < 1){dielectric = 1;}
		return dielectric;
	}
	return _dielectric;
}

double residue::approximateDipoleDipolePolarization(UInt _atomIndex1, UInt _atomIndex2)
{	//Approximate the polarizability of inclusion in medium due to the sum of the dipole angle product of the static polarizabilities
	
	//Identify Atoms in Dipoles
	dblVec a1Coords = getCoords(_atomIndex1); dblVec p1Coords;
	bool dipole1 = false;
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		if (isBonded(_atomIndex1, i)){
			p1Coords = getCoords(i);
			dipole1 = true;
			break;
		}
	}
	dblVec a2Coords = getCoords(_atomIndex2); dblVec p2Coords;
	bool dipole2 = false;
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		if (isBonded(_atomIndex2, i)){
			p2Coords = getCoords(i);
			dipole2 = true;
			break;
		}
	}
	double statpol1 = residueTemplate::getPolarizability(dataBase[itsType].itsAtomEnergyTypeDefinitions[_atomIndex1][0]); 
	double statpol2 = residueTemplate::getPolarizability(dataBase[itsType].itsAtomEnergyTypeDefinitions[_atomIndex2][0]);
	double pol1, pol2;
	if(dipole1){pol1 = fabs(cos(CMath::angle(p1Coords,a1Coords,a2Coords)*PI/180.0))*statpol1;} else{pol1 = statpol1;}
	if(dipole2){pol2 = fabs(cos(CMath::angle(p2Coords,a2Coords,a1Coords)*PI/180.0))*statpol2;} else{pol2 = statpol2;}
	return pol1+pol2;
}

double residue::approximateDipoleDipolePolarization(UInt _atomIndex1, residue* _other, UInt _atomIndex2)
{	//Approximate the polarizability of inclusion in medium due to the sum of the dipole angle product of the static polarizabilities
	
	//Identify Atoms in Dipoles
	dblVec a1Coords = getCoords(_atomIndex1); dblVec p1Coords;
	bool dipole1 = false;
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		if (isBonded(_atomIndex1, i)){
			p1Coords = getCoords(i);
			dipole1 = true;
			break;
		}
	}
	dblVec a2Coords = _other->getCoords(_atomIndex2); dblVec p2Coords;
	bool dipole2 = false;
	for (UInt i = 0; i < _other->itsAtoms.size(); i++)
	{
		if (_other->isBonded(_atomIndex2, i)){
			p2Coords = _other->getCoords(i);
			dipole2 = true;
			break;
		}
	}
	double statpol1 = residueTemplate::getPolarizability(dataBase[itsType].itsAtomEnergyTypeDefinitions[_atomIndex1][0]); 
	double statpol2 = residueTemplate::getPolarizability(dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[_atomIndex2][0]);
	double pol1, pol2;
	if(dipole1){pol1 = fabs(cos(CMath::angle(p1Coords,a1Coords,a2Coords)*PI/180.0))*statpol1;} else{pol1 = statpol1;}
	if(dipole2){pol2 = fabs(cos(CMath::angle(p2Coords,a2Coords,a1Coords)*PI/180.0))*statpol2;} else{pol2 = statpol2;}
	return pol1+pol2;
}

void residue::updateMovedDependence(residue* _other, UInt _EorC)
{	
	bool inCube;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			for(UInt j=0; j<_other->itsAtoms.size(); j++)
			{
				if (!_other->itsAtoms[j]->getSilentStatus())
				{
					inCube = itsAtoms[i]->inCube(_other->itsAtoms[j], cutoffDistance);
					if (inCube)
					{
						_other->setMoved(true, _EorC);
						_other->setCheckMovedDependence(false, _EorC);
						return;
					}
				}
			}
		}
	}
}

// end protEnergy functions ---------------------------------------------------------------------

void residue::listConnectivity()
{
    for(UInt i=0; i<itsAtoms.size(); i++)
    {
        cout << itsAtoms[i]->getName();
        for(UInt j=i+1; j<itsAtoms.size(); j++)
        {
            // ** get distance
            double distanceSquared = itsAtoms[i]->distanceSquared(itsAtoms[j]);
            if (distanceSquared < 4)
            {
                cout << " " << itsAtoms[j]->getName();
            }
        }
        cout << endl;
    }
}


double residue::interEnergy(residue* _other)
{
	double distanceSquared;
	int index1;
	int index2;
	double interEnergy = 0.0;
	double vdwEnergy = 0.0;
	double pmfEnergy = 0.0;
	double amberElecEnergy = 0.0;
	bool threeBonds;
	bool withinCutoff;
	for(UInt i=0; i<itsAtoms.size(); i++)
	{
		if (!itsAtoms[i]->getSilentStatus())
		{
			for(UInt j=0; j<_other->itsAtoms.size(); j++)
			{
				if (!_other->itsAtoms[j]->getSilentStatus())
				{
					withinCutoff = itsAtoms[i]->inCutoffSQ(_other->itsAtoms[j], cutoffDistance, cutoffDistanceSquared);
					if (withinCutoff)
					{
						distanceSquared = itsAtoms[i]->distanceSquared(_other->itsAtoms[j]);
						threeBonds = isSeparatedByFewBonds(this,i,_other,j);
						
						// ** inter AMBER Electrostatics
						if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0)
						{
							if (!threeBonds)
							{
								UInt resType1 = itsType;
								UInt resType2 = _other->itsType;
								UInt index1 = i;
								UInt index2 = j;
								double tempAmberElecEnergy = residueTemplate::getAmberElecEnergySQ(resType1, index1, resType2, index2, distanceSquared);
								amberElecEnergy += tempAmberElecEnergy;
							}
						}

						// ** inter AMBER vdW
						if (residueTemplate::itsAmberVDW.getScaleFactor() != 0.0)
						{	if (hydrogensOn)
							{	index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
								index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][0];
							}
							else
							{	index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][1];
								index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][1];
							}
							if (!threeBonds)
							{
								double tempvdwEnergy = residueTemplate::getVDWEnergySQ(index1, index2, distanceSquared);
								vdwEnergy += tempvdwEnergy;
								
				//				cout << i << " " << j << " " << tempvdwEnergy << endl;
							
							}
						}
					}
				}
			}
		}
	}
//	if (vdwEnergy > 1.0e4) cout << "inter: " << itsType << " "<< _other->itsType << " " << vdwEnergy << endl;
	interEnergy = pmfEnergy + vdwEnergy + amberElecEnergy;
	if (interEnergy < -1e5)
	{	cout << "Extremely Low interEnergy value of " << interEnergy << endl;
	}
#ifdef _INTER_ENERGY_DEBUG
	cout << "inter : " << interEnergy << " " << endl ;
#endif
	return interEnergy;
}

double residue::calculateHCA_O_hBondEnergy(residue* _other)
{
	atom* hca1 = getAtom(5);
	atom* hca2 = getAtom(6);
	atom* ca = getAtom(1);

	atom* c = _other->getAtom(2);
	atom* o = _other->getAtom(3);

	dblVec HCA1 = hca1->getCoords();
	dblVec HCA2 = hca2->getCoords();
	dblVec CA = ca->getCoords();

	dblVec C = c->getCoords();
	dblVec O = o->getCoords();

	double dist = CMath::distance(CA,O);
	double theta1 = CMath::angle(CA,HCA1,O);
	double theta2 = CMath::angle(CA,HCA2,O);
	double phi1 = CMath::angle(HCA1,O,C);
	double phi2 = CMath::angle(HCA2,O,C);

	double energy = 0.0;
	energy += 2.5*(5*pow(3.43/dist,12) - 6*pow(3.43/dist,10))*(cos(theta1)*cos(phi1) + cos(theta2)*cos(phi2));

	return energy;
}

UInt residue::getNumHardClashes()
{
	UInt numClashes = 0;
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = i + 1; j < itsAtoms.size(); j ++)
		{
			if (isClash(i,j)) numClashes++;
		}
	}
	return numClashes;
}
UInt residue::getNumHardClashes(residue* _other)
{
	UInt numClashes = 0;
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = 0; j < _other->getNumAtoms(); j ++)
		{
			if (isClash(i, _other, j)) numClashes++;
		} 
	}
	return numClashes;
}

UInt residue::getNumHardBackboneClashes(residue* _other)
{
	UInt numClashes = 0;
	if (itsAtoms.size() > 3 && _other->itsAtoms.size()> 3){
		UInt atomsI =4, atomsJ =4;
		if (itsAtoms[4]->getName() == "CB"){atomsI = 5;}

		if (_other->itsAtoms[4]->getName() == "CB"){atomsJ = 5;}
		for (UInt i = 0; i < atomsI; i ++)
		{
			for (UInt j = 0; j < atomsJ; j ++)
			{
				if (isClash(i, _other, j)) numClashes++;
			} 
		}
	}
	return numClashes;
}

bool residue::isClash(UInt _index1, UInt _index2)
{
	if (isSeparatedByFewBonds(_index1, _index2)) {return false;}
	double minDist = getRadius(_index1)+getRadius(_index2);
	double cubeLength = minDist/1.414213562; //vdw contact distance / sqrt(2) (in a square within circle for fast hard clash)
	if (itsAtoms[_index1]->inCube(itsAtoms[_index2], cubeLength)) {return true;}
	return false;
}

bool residue::isClash(UInt _index1, residue* _other, UInt _index2)
{
	if (isSeparatedByThreeBackboneBonds(_index1, _other, _index2)) {return false;}
	double minDist = getRadius(_index1)+_other->getRadius(_index2);
	double cubeLength = minDist/1.414213562; //vdw contact distance / sqrt(2) (in a square within circle for fast hard clash)
	if (itsAtoms[_index1]->inCube(_other->itsAtoms[_index2], cubeLength)) {return true;}
	return false;
}

bool residue::inCube(const residue* _other, double _cutoff)
{
	if (fabs(itsAtoms[1]->itsCoords[0] - _other->itsAtoms[1]->getX()) > _cutoff) return false;
	if (fabs(itsAtoms[1]->itsCoords[1] - _other->itsAtoms[1]->getY()) > _cutoff) return false;
	if (fabs(itsAtoms[1]->itsCoords[2] - _other->itsAtoms[1]->getZ()) > _cutoff) return false;
	return true;
}


bool residue::isSeparatedByOneOrTwoBonds(UInt _index1, UInt _index2)
{
	// Note: all interaction up to and including 1-3 interactions are skipped

	vector< UInt > bondedList1;
	vector< UInt > bondedList2;
	UInt sizeOfList1;
	UInt sizeOfList2;

	bondedList1 = dataBase[itsType].getBondingPattern(_index1);
	if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == _index2)
			{	return true;
			}
			bondedList2 = dataBase[itsType].getBondingPattern(bondedList1[i]);
			if ( (sizeOfList2 = bondedList2.size()) )
			{
				for (UInt j=0; j< sizeOfList2; j++)
				{
					if (bondedList2[j] == _index2)
					{	return true;
					}
				}
			}
		}
	}
	return false;
}

bool residue::isSeparatedByFewBonds(UInt _index1, UInt _index2)
{
	// Note: all interaction up to and including 1-4 interactions are skipped

	vector< UInt > bondedList1;
	vector< UInt > bondedList2;
	vector< UInt > bondedList3;
	UInt sizeOfList1;
	UInt sizeOfList2;
	UInt sizeOfList3;

	bondedList1 = dataBase[itsType].getBondingPattern(_index1);

	if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == _index2)
			{	return true;
			}
			bondedList2 = dataBase[itsType].getBondingPattern(bondedList1[i]);
			if ( (sizeOfList2 = bondedList2.size()) )
			{
				for (UInt j=0; j< sizeOfList2; j++)
				{
					if (bondedList2[j] == _index2)
					{	return true;
					}
					bondedList3 = dataBase[itsType].getBondingPattern(bondedList2[j]);
					if ( (sizeOfList3 = bondedList3.size()) )
					{
						for (UInt k=0; k< sizeOfList3; k++)
						{
							if (bondedList3[k] == _index2)
							{	return true;
							}
						}
					}
				}
			}
		}
	}
	return false;
}

UInt residue::getBondSeparation(UInt _index1, UInt _index2)
{
    // Note: all interaction up to and including 1-4 interactions are skipped

    vector< UInt > bondedList1;
    vector< UInt > bondedList2;
    vector< UInt > bondedList3;
    UInt sizeOfList1;
    UInt sizeOfList2;
    UInt sizeOfList3;

    bondedList1 = dataBase[itsType].getBondingPattern(_index1);

    if ( (sizeOfList1 = bondedList1.size()) )
    {
        for (UInt i=0; i < sizeOfList1; i++)
        {
            if (bondedList1[i] == _index2)
            {	return 0;
            }
            bondedList2 = dataBase[itsType].getBondingPattern(bondedList1[i]);
            if ( (sizeOfList2 = bondedList2.size()) )
            {
                for (UInt j=0; j< sizeOfList2; j++)
                {
                    if (bondedList2[j] == _index2)
                    {	return 1;
                    }
                    bondedList3 = dataBase[itsType].getBondingPattern(bondedList2[j]);
                    if ( (sizeOfList3 = bondedList3.size()) )
                    {
                        for (UInt k=0; k< sizeOfList3; k++)
                        {
                            if (bondedList3[k] == _index2)
                            {	return 2;
                            }
                        }
                    }
                }
            }
        }
    }
    return 99;
}

bool residue::isSeparatedByThreeBackboneBonds(UInt _index1, residue* _pRes2, UInt _index2)
{
	// first, check if they are sequential residues
	int theOrder = 0;
	// is _pRes1 the residue N-terminal to _pRes2?
	if ( this == _pRes2->getPrevRes())
		theOrder = 1;
	// if _pRes1 the residue C-terminal to _pRes2?
	if ( this == _pRes2->getNextRes())
		theOrder = -1;
	if (theOrder == 0)
	{
		return false;
	}
	string name1, name2;
	UInt atom1, atom2;
	if (theOrder == -1){
		name2 = itsAtoms[_index1]->getName(); // c terminal
		name1 = _pRes2->itsAtoms[_index2]->getName(); //n terminal
		atom2 = _index1;
		atom1 = _index2;
	}
	else{
		name1 = itsAtoms[_index1]->getName(); // n terminal
		name2 = _pRes2->itsAtoms[_index2]->getName(); // c terminal
		atom1 = _index1;
		atom2 = _index2;
	}
	
	if (atom1 == 0 && atom2 == 0){
		return true;
	}
	if (atom1 == 1 && (atom2 == 0 || atom2 == 1 || name2 == "H")){
		return true;
	}
	if (atom1 == 2 && (atom2 == 0 || atom2 == 1 || atom2 == 2 || name2 == "H" || name2 == "CB" || name2 == "HA" || name2 == "HA2" || name2 == "HA3")){
		return true;
	}
	if (atom1 == 3 && (atom2 == 0 || atom2 == 1 || name2 == "H")){
		return true;
	}
	return false;
}

bool residue::isSeparatedByOneOrTwoBackboneBonds(UInt _index1, residue* _pRes2, UInt _index2)
{
	// first, check if they are sequential residues
	int theOrder = 0;
	// is _pRes1 the residue N-terminal to _pRes2?
	if ( this == _pRes2->getPrevRes())
		theOrder = 1;
	// if _pRes1 the residue C-terminal to _pRes2?
	if ( this == _pRes2->getNextRes())
		theOrder = -1;
	if (theOrder == 0)
	{
		return false;
	}
	string name1, name2;
	UInt atom1, atom2;
	if (theOrder == -1){
		name2 = itsAtoms[_index1]->getName(); // c terminal
		name1 = _pRes2->itsAtoms[_index2]->getName(); //n terminal
		atom2 = _index1;
		atom1 = _index2;
	}
	else{
		name1 = itsAtoms[_index1]->getName(); // n terminal
		name2 = _pRes2->itsAtoms[_index2]->getName(); // c terminal
		atom1 = _index1;
		atom2 = _index2;
	}
	
	if (atom2 == 0 || atom2 == 1 || name2 == "H")
	{
		if (atom1 > 0 && atom1 < 4){
			if ((name2 == "H" || atom2 == 1) && atom1 == 2){return true;}
			if (atom2 == 0) {return true;}
			return false;
		}
		else{return false;}
	}
	else{return false;}
}


bool residue::isSeparatedByOneOrTwoBonds(UInt _index1, residue* _pRes2, UInt _index2)
{
	// first, check if they are sequential residues
	int theOrder = 0;
	// is _pRes1 the residue N-terminal to _pRes2?
	if ( this == _pRes2->getPrevRes())
		theOrder = 1;
	// if _pRes1 the residue C-terminal to _pRes2?
	if ( this == _pRes2->getNextRes())
		theOrder = -1;
	if (theOrder == 0)
	{
		return false;
	}
	
	// ok, now we know we've got two sequential amino acids.
	// find out how far the atom in question in the N-term
	// amino acid is from the carboxyl carbon (the end).
	UInt Cindex;
	vector< UInt > bondedList1;
	vector< UInt > bondedList2;
	UInt sizeOfList1;
	UInt sizeOfList2;
	UInt typeOfRes1;
	UInt typeOfRes2;
	UInt atomIndex1;
	UInt atomIndex2;

	if (theOrder == 1)
	{
		typeOfRes1 = getTypeIndex();
		atomIndex1 = _index1;
		typeOfRes2 = _pRes2->getTypeIndex();
		atomIndex2 = _index2;
	}
	else
	{
		typeOfRes1 = _pRes2->getTypeIndex();
		atomIndex1 = _index2;
		typeOfRes2 = getTypeIndex();
		atomIndex2 = _index1;
	}

	// find the index of "C" should be at mainchain end -2
	UInt mcsize = dataBase[typeOfRes1].mainChain.size();
	if (mcsize > 2)
	{	Cindex = mcsize -2;
	}
	else
	{	cout << "Error! mainChain size too small" << endl;
		return false;
	}

	bondedList1 = dataBase[typeOfRes1].getBondingPattern(atomIndex1);
	UInt firstIterationLevel = 10;
	if ( atomIndex1 == Cindex)
	{
		firstIterationLevel = 0;
	}
	else if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == Cindex)
			{	firstIterationLevel = 1;
			}
			else{
				bondedList2 = dataBase[typeOfRes1].getBondingPattern(bondedList1[i]);
				if ( (sizeOfList2 = bondedList2.size()) )
				{
					for (UInt j=0; j< sizeOfList2; j++)
					{
						if (bondedList2[j] == Cindex)
						{	firstIterationLevel = 2;
						}
					}
				}
			}
		}
	}

	bondedList1 = dataBase[typeOfRes2].getBondingPattern(0);
	UInt secondIterationLevel = 10;

	if (atomIndex2 == 0)
	{	secondIterationLevel = 0;
	}
	else if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == atomIndex2)
			{	secondIterationLevel = 1;
			}
			else{
				bondedList2 = dataBase[typeOfRes2].getBondingPattern(bondedList1[i]);
				if ( (sizeOfList2 = bondedList2.size()) )
				{
					for (UInt j=0; j< sizeOfList2; j++)
					{
						if (bondedList2[j] == atomIndex2)
						{	secondIterationLevel = 2;
						}
					}
				}
			}
		}
	}
	UInt numBonds = firstIterationLevel + secondIterationLevel;
	if (numBonds < 2 )
	{
		return true;
	}
	return false;
}

bool residue::isSeparatedByFewBonds(residue* _pRes1, UInt _index1, residue*
 _pRes2, UInt _index2)
{
	// Note: all interaction up to and including 1-4 interactions are skipped

	// first, check if they are sequential residues
	int theOrder = 0;
	// is _pRes1 the residue N-terminal to _pRes2?
	if ( _pRes1 == _pRes2->getPrevRes())
		theOrder = 1;
	// if _pRes1 the residue C-terminal to _pRes2?
	if ( _pRes1 == _pRes2->getNextRes())
		theOrder = -1;

	// if they're not sequential, there's no way that their
	// atoms can be within 3 bonds.... except if they are
	// Cysteine sulfurs involved in a disulfide bond....

	// Check for disulfide
	if ((_pRes1->getType() == "CYX" || _pRes1->getType() == "CXD") && _index1 == 5 )
	{	if ((_pRes2->getType() == "CYX" ||_pRes1->getType() == "CXD") && _index2 == 5)
		{   atom* pAtom1=_pRes1->getAtom(_index1);
			atom* pAtom2=_pRes2->getAtom(_index2);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}

	// Check for ASP/GLU - OEC
	if ((_pRes1->getType() == "ASP" && (_index1 == 6 || _index1 == 7)) || (_pRes1->getType() == "GLU" && (_index1 == 8 || _index1 == 7)))
	{	if (_pRes2->getType() == "OEC" && (_index2 == 0 || _index2 == 2 || _index2 == 4 || _index2 == 6 || _index2 == 8))
		{ 	atom* pAtom1=_pRes1->getAtom(_index1);
			atom* pAtom2=_pRes2->getAtom(_index2);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}
	if ((_pRes2->getType() == "ASP" && (_index2 == 6 || _index2 == 7)) || (_pRes2->getType() == "GLU" && (_index2 == 8 || _index2 == 7)))
	{	if (_pRes1->getType() == "OEC" && (_index1 == 0 || _index1 == 2 || _index1 == 4 || _index1 == 6 || _index1 == 8))
		{ 	atom* pAtom1=_pRes1->getAtom(_index1);
			atom* pAtom2=_pRes2->getAtom(_index2);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}

	// Check for SF4 Cysteine bond
	if ((_pRes1->getType() == "CYF" || _pRes1->getType() == "CFD") && _index1 == 5 )
	{	if (_pRes2->getType() == "CSF" && (_index2 == 6 || _index2 == 8 || _index2 == 10 || _index2 == 12))
		{ 	atom* pAtom1=_pRes1->getAtom(_index1);
			atom* pAtom2=_pRes2->getAtom(_index2);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}
	if ((_pRes2->getType() == "CYF" || _pRes2->getType() == "CFD") && _index2 == 5)
	{	if (_pRes1->getType() == "CSF" && (_index1 == 6 || _index1 == 8 || _index1 == 10 || _index1 == 12))
		{ 	atom* pAtom1=_pRes2->getAtom(_index2);
			atom* pAtom2=_pRes1->getAtom(_index1);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}

	if ((_pRes1->getType() == "CYF" || _pRes1->getType() == "CFD") && _index1 == 5 )
	{	if (_pRes2->getType() == "SF4")
		{ 	atom* pAtom1=_pRes1->getAtom(_index1);
			atom* pAtom2=_pRes2->getAtom(_index2);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}
	if ((_pRes2->getType() == "CYF" || _pRes2->getType() == "CFD") && _index2 == 5)
	{	if (_pRes1->getType() == "SF4")
		{ 	atom* pAtom1=_pRes2->getAtom(_index2);
			atom* pAtom2=_pRes1->getAtom(_index1);
			if (pAtom1->distance(pAtom2) < 3.0)
			{
				return true;
			}
		}
	}

    // Check for HIS-Heme bond
    if ((_pRes1->getType() == "HID" || _pRes1->getType() == "HDD") && _index1 > 5 )
    {	if (_pRes2->getType() == "HEM")
        { 	atom* pAtom1=_pRes1->getAtom(_index1);
            atom* pAtom2=_pRes2->getAtom(_index2);
            if (pAtom1->distance(pAtom2) < 3.0)
            {
                return true;
            }
        }
    }
    if ((_pRes2->getType() == "HID" || _pRes2->getType() == "HDD") && _index2 > 5)
    {	if (_pRes1->getType() == "HEM")
        { 	atom* pAtom1=_pRes2->getAtom(_index2);
            atom* pAtom2=_pRes1->getAtom(_index1);
            if (pAtom1->distance(pAtom2) < 3.0)
            {
                return true;
            }
        }
    }

	// Check for peptide bond
	if ((_index1 == 0 && _index2 == 2) || (_index2 == 0 && _index1 == 2))
	{
		atom* pAtom1=_pRes1->getAtom(_index1);
		atom* pAtom2=_pRes2->getAtom(_index2);
		if (pAtom1->distance(pAtom2) < 2.0)
		{
			return true;
		}
	}

	if (theOrder == 0)
	{
#ifdef _SKIP_DEBUG
		cout << "Not sequential amino acids" << endl;
#endif
		return false;
	}

#ifdef _SKIP_DEBUG
	cout << "Order = " << theOrder << "  ";
#endif

	// ok, now we know we've got two sequential amino acids.
	// find out how far the atom in question in the N-term
	// amino acid is from the carboxyl carbon (the end).
	UInt Cindex;
	vector< UInt > bondedList1;
	vector< UInt > bondedList2;
	UInt sizeOfList1;
	UInt sizeOfList2;
	UInt typeOfRes1;
	UInt typeOfRes2;
	UInt atomIndex1;
	UInt atomIndex2;

	if (theOrder == 1)
	{
		typeOfRes1 = _pRes1->getTypeIndex();
		atomIndex1 = _index1;
		typeOfRes2 = _pRes2->getTypeIndex();
		atomIndex2 = _index2;
	}
	else
	{
		typeOfRes1 = _pRes2->getTypeIndex();
		atomIndex1 = _index2;
		typeOfRes2 = _pRes1->getTypeIndex();
		atomIndex2 = _index1;
	}

	// find the index of "C" should be at mainchain end -2
	UInt mcsize = dataBase[typeOfRes1].mainChain.size();
	if (mcsize > 2)
	{	Cindex = mcsize -2;
#ifdef _SKIP_DEBUG
		cout << "Cindex=" << Cindex << " ";
#endif
	}
	else
	{	cout << "Error! mainChain size too small" << endl;
		return false;
	}

	bondedList1 = dataBase[typeOfRes1].getBondingPattern(atomIndex1);
	UInt firstIterationLevel = 10;
	if ( atomIndex1 == Cindex)
	{
		firstIterationLevel = 0;
#ifdef _SKIP_DEBUG
		cout << "| FirstAtom is C | ";
#endif
	}
	else if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == Cindex)
			{	firstIterationLevel = 1;
			}
			bondedList2 = dataBase[typeOfRes1].getBondingPattern(bondedList1[i]);
			if ( (sizeOfList2 = bondedList2.size()) )
			{
				for (UInt j=0; j< sizeOfList2; j++)
				{
					if (bondedList2[j] == Cindex)
					{	firstIterationLevel = 2;
					}
				}
			}
		}
#ifdef _SKIP_DEBUG
		cout << " Steps To C From atom1 =" << firstIterationLevel << "  ";
#endif
	}

	bondedList1 = dataBase[typeOfRes2].getBondingPattern(0);
	UInt secondIterationLevel = 10;

	if (atomIndex2 == 0)
	{	secondIterationLevel = 0;
#ifdef _SKIP_DEBUG
		cout << " |atom2 is N| ";
#endif
	}
	else if ( (sizeOfList1 = bondedList1.size()) )
	{
		for (UInt i=0; i < sizeOfList1; i++)
		{
			if (bondedList1[i] == atomIndex2)
			{	secondIterationLevel = 1;
			}
			bondedList2 = dataBase[typeOfRes2].getBondingPattern(bondedList1[i]);
			if ( (sizeOfList2 = bondedList2.size()) )
			{
				for (UInt j=0; j< sizeOfList2; j++)
				{
					if (bondedList2[j] == atomIndex2)
					{	secondIterationLevel = 2;
					}
				}
			}
		}
#ifdef _SKIP_DEBUG
		cout << " Steps To N from atom2 =" << secondIterationLevel << " ";
#endif
	}

	UInt numBonds = firstIterationLevel + secondIterationLevel;
#ifdef _SKIP_DEBUG
	cout << " numbonds=" << numBonds << " ";
#endif
	if (numBonds < 3 )
	{
#ifdef _SKIP_DEBUG
		cout << " skipping " << endl;
#endif
		return true;
	}

	return false;
}

UInt residue::getBondSeparation(residue* _pRes1, UInt _index1, residue*
 _pRes2, UInt _index2)
{
    // Note: all interaction up to and including 1-4 interactions are skipped

    // first, check if they are sequential residues
    int theOrder = 0;
    // is _pRes1 the residue N-terminal to _pRes2?
    if ( _pRes1 == _pRes2->getPrevRes())
        theOrder = 1;
    // if _pRes1 the residue C-terminal to _pRes2?
    if ( _pRes1 == _pRes2->getNextRes())
        theOrder = -1;

    // if they're not sequential, there's no way that their
    // atoms can be within 3 bonds.... except if they are
    // Cysteine sulfurs involved in a disulfide bond....

    // Check for disulfide
    if (_pRes1->getType() == "CYX" || _pRes1->getType() == "CXD" )
    {	if (_pRes2->getType() == "CYX" ||_pRes1->getType() == "CXD" )
        {   atom* pAtom1=_pRes1->getAtom(_index1);
            atom* pAtom2=_pRes2->getAtom(_index2);
            if (pAtom1->distance(pAtom2) < 3.0)
            {
                return 1;
            }
        }
    }

    // Check for SF4 Cysteine bond
    if (_pRes1->getType() == "CYF" || _pRes1->getType() == "CFD" )
    {	if (_pRes2->getType() == "CSF")
        { 	atom* pAtom1=_pRes1->getAtom(_index1);
            atom* pAtom2=_pRes2->getAtom(_index2);
            if (pAtom1->distance(pAtom2) < 4.0)
            {
                return 1;
            }
        }
    }
    if (_pRes2->getType() == "CYF" || _pRes2->getType() == "CFD" )
    {	if (_pRes1->getType() == "CSF")
        { 	atom* pAtom1=_pRes2->getAtom(_index2);
            atom* pAtom2=_pRes1->getAtom(_index1);
            if (pAtom1->distance(pAtom2) < 4.0)
            {
                return 1;
            }
        }
    }

    // Check for peptide bond
    if ((_index1 == 0 && _index2 == 2) || (_index2 == 0 && _index1 == 2))
    {
        atom* pAtom1=_pRes1->getAtom(_index1);
        atom* pAtom2=_pRes2->getAtom(_index2);
        if (pAtom1->distance(pAtom2) < 1.6)
        {
            return 1;
        }
    }

    if (theOrder == 0)
    {
#ifdef _SKIP_DEBUG
        cout << "Not sequential amino acids" << endl;
#endif
        return 99;
    }

#ifdef _SKIP_DEBUG
    cout << "Order = " << theOrder << "  ";
#endif

    // ok, now we know we've got two sequential amino acids.
    // find out how far the atom in question in the N-term
    // amino acid is from the carboxyl carbon (the end).
    UInt Cindex;
    vector< UInt > bondedList1;
    vector< UInt > bondedList2;
    UInt sizeOfList1;
    UInt sizeOfList2;
    UInt typeOfRes1;
    UInt typeOfRes2;
    UInt atomIndex1;
    UInt atomIndex2;

    if (theOrder == 1)
    {
        typeOfRes1 = _pRes1->getTypeIndex();
        atomIndex1 = _index1;
        typeOfRes2 = _pRes2->getTypeIndex();
        atomIndex2 = _index2;
    }
    else
    {
        typeOfRes1 = _pRes2->getTypeIndex();
        atomIndex1 = _index2;
        typeOfRes2 = _pRes1->getTypeIndex();
        atomIndex2 = _index1;
    }

    // find the index of "C" should be at mainchain end -2
    UInt mcsize = dataBase[typeOfRes1].mainChain.size();
    if (mcsize > 2)
    {	Cindex = mcsize -2;
#ifdef _SKIP_DEBUG
        cout << "Cindex=" << Cindex << " ";
#endif
    }
    else
    {	cout << "Error! mainChain size too small" << endl;
        return 99;
    }

    bondedList1 = dataBase[typeOfRes1].getBondingPattern(atomIndex1);
    UInt firstIterationLevel = 10;
    if ( atomIndex1 == Cindex)
    {
        firstIterationLevel = 0;
#ifdef _SKIP_DEBUG
        cout << "| FirstAtom is C | ";
#endif
    }
    else if ( (sizeOfList1 = bondedList1.size()) )
    {
        for (UInt i=0; i < sizeOfList1; i++)
        {
            if (bondedList1[i] == Cindex)
            {	firstIterationLevel = 1;
            }
            bondedList2 = dataBase[typeOfRes1].getBondingPattern(bondedList1[i]);
            if ( (sizeOfList2 = bondedList2.size()) )
            {
                for (UInt j=0; j< sizeOfList2; j++)
                {
                    if (bondedList2[j] == Cindex)
                    {	firstIterationLevel = 2;
                    }
                }
            }
        }
#ifdef _SKIP_DEBUG
        cout << " Steps To C From atom1 =" << firstIterationLevel << "  ";
#endif
    }

    bondedList1 = dataBase[typeOfRes2].getBondingPattern(0);
    UInt secondIterationLevel = 10;

    if (atomIndex2 == 0)
    {	secondIterationLevel = 0;
#ifdef _SKIP_DEBUG
        cout << " |atom2 is N| ";
#endif
    }
    else if ( (sizeOfList1 = bondedList1.size()) )
    {
        for (UInt i=0; i < sizeOfList1; i++)
        {
            if (bondedList1[i] == atomIndex2)
            {	secondIterationLevel = 1;
            }
            bondedList2 = dataBase[typeOfRes2].getBondingPattern(bondedList1[i]);
            if ( (sizeOfList2 = bondedList2.size()) )
            {
                for (UInt j=0; j< sizeOfList2; j++)
                {
                    if (bondedList2[j] == atomIndex2)
                    {	secondIterationLevel = 2;
                    }
                }
            }
        }
#ifdef _SKIP_DEBUG
        cout << " Steps To N from atom2 =" << secondIterationLevel << " ";
#endif
    }

    UInt numBonds = firstIterationLevel + secondIterationLevel;
#ifdef _SKIP_DEBUG
    cout << " numbonds=" << numBonds << " ";
#endif
    return numBonds;
}

bool residue::notHydrogen(UInt _atomIndex)
{
	string atomType = getTypeStringFromAtomNum(_atomIndex);
	if (atomType != "H")
	{
		return true;
	}
	return false;
}


bool residue::isBonded(UInt _index1, UInt _index2)
{
	vector< UInt > bondedList;
	bondedList = dataBase[itsType].getBondingPattern(_index1);
	UInt sizeOfList;
	// intended assignment in if statement
	if ( (sizeOfList = bondedList.size()) )
	{	
		for (UInt i=0; i< sizeOfList; i++)
		{	
			if (bondedList[i] == _index2)
			{	
				return true;
			}
		}
	}
	else
	{	
		return false;
	}
	return false;
}


bool residue::isBonded(atom* _pAtom1, atom* _pAtom2)
{
	bool bonded = false;
	atom* pAtomTemp1;
	atom* pAtomTemp2;
	pAtomTemp1 = static_cast<atom*>(_pAtom1->getChild());
#ifdef _BOND_DEBUG
	cout << _pAtom1->getName() << " ";
	cout << _pAtom2->getName();
#endif
	if (pAtomTemp1 == _pAtom2)
	{	bonded = true;
#ifdef _BOND_DEBUG
		cout << " bonded atom!" << endl;
#endif
		return bonded;
	}
	pAtomTemp2 = static_cast<atom*>(_pAtom1->getParent());
	if (pAtomTemp2 == _pAtom2)
	{	bonded = true;
#ifdef _BOND_DEBUG
		cout << " bonded atom!" << endl;
#endif
		return bonded;
	}

	if (!pAtomTemp1)
	{
#ifdef _BOND_DEBUG
		cout << endl;
#endif
		return bonded;
	}
	atom* pAtomTemp3;
	// check through sibs of temp1
	for (;(pAtomTemp3 = static_cast<atom*>(pAtomTemp1->getNextSib()));)
	{
		if (pAtomTemp3 == _pAtom2)
		{	bonded = true;
#ifdef _BOND_DEBUG
			cout << " bonded atom!" << endl;
#endif
			return bonded;
		}
		pAtomTemp1 = pAtomTemp3;
	}

	if (!pAtomTemp2)
	{
#ifdef _BOND_DEBUG
		cout << endl;
#endif
		return bonded;
	}

	atom* pAtomTemp4;
	// check through sibs of temp1
	for (;(pAtomTemp4 = static_cast<atom*>(pAtomTemp2->getNextSib()));)
	{
		if (pAtomTemp4 == _pAtom2)
		{	bonded = true;
#ifdef _BOND_DEBUG
			cout << " bonded atom!" << endl;
#endif
			return bonded;
		}
		pAtomTemp2 = pAtomTemp4;
	}

#ifdef _BOND_DEBUG
	cout << endl;
#endif
	return bonded;
}

double residue::getSelfEnergy(residue* _other)
{
	double selfEnergy = 0.0;
	UInt index1, index2;
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		bool isSideChain = true;  // for this energy term, CA is considered to be SIDECHAIN
		for (UInt counter = dataBase[itsType].mainChain[0]; counter < dataBase[itsType].mainChain.size(); counter++)
		{
			if (i == dataBase[itsType].mainChain[counter] && itsAtoms[i]->getName() != "CA") isSideChain = false;
		}
		if (isSideChain)
		{
			for (UInt j = dataBase[_other->itsType].mainChain[0]; j < dataBase[_other->itsType].mainChain.size(); j++)
			{
				bool withinCutoff = itsAtoms[i]->inCutoffSQ(_other->itsAtoms[j], cutoffDistance, cutoffDistanceSquared);
				bool isBonded;
		if (0 == 0)//_other->itsAtoms[j]->getName() != "CA")
		{
					if (this != _other)
					{
						isBonded = isSeparatedByFewBonds(this, i, _other, j);
					}
					else
						isBonded = isSeparatedByFewBonds(i,j);
					if (!isBonded && !itsAtoms[j]->getSilentStatus() && withinCutoff && residueTemplate::itsAmberVDW.getScaleFactor() != 0.0)
					{
						double distanceSquared = itsAtoms[i]->distanceSquared(_other->itsAtoms[j]);
						if (hydrogensOn)
						{
							index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][0];
							index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][0];
						}
						else
						{
							index1 = dataBase[itsType].itsAtomEnergyTypeDefinitions[i][1];
							index2 = dataBase[_other->itsType].itsAtomEnergyTypeDefinitions[j][1];
						}
				UInt resType1 = itsType;
				UInt resType2 = _other->itsType;
				double tempElecEnergy = residueTemplate::getAmberElecEnergySQ(resType1, i, resType2, j, distanceSquared);
						double tempvdwEnergy = residueTemplate::getVDWEnergySQ(index1, index2, distanceSquared);
						selfEnergy += tempvdwEnergy + tempElecEnergy;
						// atom* temp = _other->itsAtoms[j];
						//cout <<  itsAtoms[i]->getName() << " to " << temp->getName()  << " " << selfEnergy <<endl;
					}
		}
			}
		}
	}
	return selfEnergy;
}

bool residue::getMoved(UInt EorC)
{
	if (EorC == 0){return movedE;}
	if (EorC == 1){return movedC;}
	else {return movedB;}
}

void residue::setMoved()
{
		movedE = true;
		movedC = true;
		movedB = true;
		setEnergy(0.0);
		clearEnvironment();
		setClashes(0);
		setBackboneClashes(0);
		setCheckMovedDependence(true, 0);
		setCheckMovedDependence(true, 1);
		setCheckMovedDependence(true, 2);
}

void residue::setMoved(bool _moved, UInt _EorC)
{
	if(_EorC == 0){
		movedE = _moved;
		if (_moved){
			setEnergy(0.0);
			clearEnvironment();
			setCheckMovedDependence(true, _EorC);
		}
	}
	if(_EorC == 1){
		movedC = _moved;
		if (_moved){
			setClashes(0);
			setCheckMovedDependence(true, _EorC);
		}
	}
	if(_EorC == 2){
		movedB = _moved;
		if (_moved){
			setBackboneClashes(0);
			setCheckMovedDependence(true, _EorC);
		}
	}
	if (!_moved){
		setCheckMovedDependence(false, _EorC);
	}
}

bool residue::getCheckMovedDependence(UInt _EorC)
{
	if(_EorC == 0){return dependentMoveE;}
	if(_EorC == 1){return dependentMoveC;}
	else{return dependentMoveB;}
}

void residue::setCheckMovedDependence(bool _check, UInt _EorC)
{
	if(_EorC == 0){
		dependentMoveE = _check;
	}
	if(_EorC == 1){
		dependentMoveC = _check;
	}
	if(_EorC == 2){
		dependentMoveB = _check;
	}
}

void residue::clearEnvironment()
{
	for (UInt i=0; i < itsAtoms.size(); i++)
	{
		itsAtoms[i]->setEnvVol(0.0);
	}
}

void residue::setClashes(UInt _clashes)
{
	clashes = _clashes;
}

void residue::sumClashes(UInt _clashes)
{
	clashes += _clashes;
}

void residue::setBackboneClashes(UInt _clashes)
{
	clashesB = _clashes;
}

void residue::sumBackboneClashes(UInt _clashes)
{
	clashesB += _clashes;
}

void residue::setEnergy(double _Energy)
{
	Energy = _Energy;
}

void residue::sumEnergy(double _Energy)
{
	Energy += _Energy;
}

void residue::setResiduesPerTurnType(UInt _RPT)
{
	RPT = _RPT;
	for (UInt i = 0; i < itsAtoms.size(); i++)
	{
		itsAtoms[i]->setRPTType(_RPT);
	}
}

double residue::getVolume(UInt _method)
{
	double itsVolume = 0.0;

	if (_method == 0)
		itsVolume = wodakVolume();
	else
	{
		cout << "that volume method does not exist\n";
		return -1;
	}
	return itsVolume;
}

double residue::netCharge()
{
	double nCharge = 0.0;
	string residueType = getDataBaseItem(itsType);
	if (residueType == "ARG") nCharge =  1.0;
	if (residueType == "ASP") nCharge = -1.0;
	if (residueType == "GLU") nCharge = -1.0;
	if (residueType == "HIN") nCharge = -1.0;
	if (residueType == "HIP") nCharge =  1.0;
	if (residueType == "LYS") nCharge =  1.0;
	if (residueType == "ARD") nCharge =  1.0;
	if (residueType == "APD") nCharge = -1.0;
	if (residueType == "GUD") nCharge = -1.0;
	if (residueType == "HND") nCharge = -1.0;
	if (residueType == "HPD") nCharge =  1.0;
	if (residueType == "LYD") nCharge =  1.0;
	return nCharge;
}


// wodakVolume() returns mean residue volume in cubic Angstroms based on parameterization
// in Pontius, Richelle, Wodak, JMB v 264 p 121 (1996).  Volumes refer specifically to
// Table 2, column 1
double residue::wodakVolume()
{
	double volume = 0;
	string residueType = getDataBaseItem(itsType);
	if (residueType == "ALA") volume = 91.5;
	if (residueType == "ARG") volume = 196.1;
	if (residueType == "ASN") volume = 138.3;
	if (residueType == "ASP") volume = 135.2;
	if (residueType == "ASH") volume = 135.2;
	if (residueType == "CYS") volume = 114.4;
	if (residueType == "CYX") volume = 114.4;
	if (residueType == "GLN") volume = 156.4;
	if (residueType == "GLU") volume = 154.6;
	if (residueType == "GLH") volume = 154.6;
	if (residueType == "GLY") volume = 67.5;
	if (residueType == "HID") volume = 163.2;
	if (residueType == "HIE") volume = 163.2;
	if (residueType == "HIN") volume = 163.2;
	if (residueType == "HIP") volume = 163.2;
	if (residueType == "ILE") volume = 162.6;
	if (residueType == "LEU") volume = 163.4;
	if (residueType == "LYS") volume = 162.5;
	if (residueType == "MET") volume = 165.9;
	if (residueType == "PHE") volume = 198.8;
	if (residueType == "PRO") volume = 123.4;
	if (residueType == "SER") volume = 102.0;
	if (residueType == "THR") volume = 126.0;
	if (residueType == "TRP") volume = 237.2;
	if (residueType == "TYR") volume = 209.8;
	if (residueType == "VAL") volume = 138.4;
	if (residueType == "ALD") volume = 91.5;
	if (residueType == "ARD") volume = 196.1;
	if (residueType == "AND") volume = 138.3;
	if (residueType == "APD") volume = 135.2;
	if (residueType == "AHD") volume = 135.2;
	if (residueType == "CYD") volume = 114.4;
	if (residueType == "CXD") volume = 114.4;
	if (residueType == "GND") volume = 156.4;
	if (residueType == "GUD") volume = 154.6;
	if (residueType == "GHD") volume = 154.6;
	if (residueType == "GLY") volume = 67.5;
	if (residueType == "HDD") volume = 163.2;
	if (residueType == "HED") volume = 163.2;
	if (residueType == "HND") volume = 163.2;
	if (residueType == "HPD") volume = 163.2;
	if (residueType == "ILD") volume = 162.6;
	if (residueType == "LED") volume = 163.4;
	if (residueType == "LYD") volume = 162.5;
	if (residueType == "MED") volume = 165.9;
	if (residueType == "PHD") volume = 198.8;
	if (residueType == "PRD") volume = 123.4;
	if (residueType == "SED") volume = 102.0;
	if (residueType == "THD") volume = 126.0;
	if (residueType == "TRD") volume = 237.2;
	if (residueType == "TYD") volume = 209.8;
	if (residueType == "VAD") volume = 138.4;
	return volume;
}

void residue::coilcoil(const double _pitch)
{
    for (UInt i = 0; i < itsAtoms.size(); i ++)
    {
        double x = itsAtoms[i]->getX();
        double y = itsAtoms[i]->getY();
        double z = itsAtoms[i]->getZ();

        double theta = -6.28318 * z / _pitch;

        double xcoil = x * cos(theta) - y * sin(theta);
        double ycoil = x * sin(theta) + y * cos(theta);
        double zcoil = z;

        itsAtoms[i]->setCoords(xcoil, ycoil, zcoil);
        }
    return;
}

void residue::initializeSpherePoints()
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		itsAtoms[i]->initializeSpherePoints();
	}
	return;
}

void residue::removeIntraResidueSpherePoints()
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = 0; j < itsAtoms.size(); j ++)
		{
			itsAtoms[i]->removeSpherePoints(itsAtoms[j]);
		}
	}

	return;
}

void residue::removeInterResidueSpherePoints(residue* _other)
{
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		for (UInt j = 0; j < _other->itsAtoms.size(); j ++)
		{
			itsAtoms[i]->removeSpherePoints(_other->itsAtoms[j]);
		}
	}
	return;
}

double residue::tabulateSurfaceArea()
{
	double surfaceArea = 0.0;
	for (UInt i = 0; i < itsAtoms.size(); i ++)
	{
		surfaceArea += itsAtoms[i]->calculateExposedSASA();
	}
	return surfaceArea;
}

double residue::tabulateSurfaceArea(UInt _atomIndex)
{
	double surfaceArea;
	surfaceArea = itsAtoms[_atomIndex]->calculateExposedSASA();
	return surfaceArea;
}

dblVec residue::getBackBoneCentroid()
{
	dblVec centroid(3);
	atom* temp = getAtom(1);  // CA
	centroid = temp->getCoords();
	return centroid;
}
