// Extra Class for TsScorePDB4DNA
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:											      *
// * Copyright (c) <2015>, <Harald Paganetti>                         *
// * All rights reserved.                                             *
// * Redistribution    and   use in   source and   binary    forms,   *
// * with or without modification, are permitted provided that the    *
// * following conditions are met:                                    *
// *                                                                  *
// *                                                                  *
// * 1. Redistributions of source code must retain the above          *
// * copyright notice, this                                           *
// * list of conditions and the following disclaimer.                 *
// * 2. Redistributions in binary form must reproduce the above       *
// * copyright notice, this list of conditions and the following      *
// * disclaimer in the documentation and/or other materials provided  *
// * with the distribution.                                           *
// *                                                                  *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
// * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
// * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
// * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
// * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR             *
// * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT *
// * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF *
// * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  *
// * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      *
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
// * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF   *
// * THE POSSIBILITY OF SUCH DAMAGE.                                  *
// *                                                                  *
// * The views and conclusions contained in the software and          *
// * documentation are those of the authors and should not be         *
// * interpreted as representing official policies, either expressed  *
// * or implied, of the FreeBSD Project.                              *
// *                                                                  *
// * Contacts: Jan Schuemann, jschuemann@mgh.harvard.edu              *
// *           Harald Paganetti, hpaganetti@mgh.harvard.edu           *
// *                                                                  *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
// 
/// \file PDBlib.cc
/// \brief Implementation file for PDBlib class

//define if the program is running with Geant4
#define GEANT4

#include "PDBlib.hh"
#ifdef GEANT4
//Specific to Geant4, globals.hh is used for G4cout
#include "globals.hh"
#endif
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <string>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PDBlib::PDBlib():fNbNucleotidsPerStrand(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * \brief    Load a PDB file into memory
 * \details  molecule (polymer,?), 
 *                  read line, key words
 * \param    filename    string for filename
 * \param    verbosity
 * \return   List of molecules
 */

Molecule * PDBlib::Load( const string &filename,unsigned short int &isDNA,
    unsigned short int verbose=0)
{
  string sLine = "";
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
#ifdef GEANT4
    G4cout<<"PDBlib::load >> file "<<filename<<" not found !!!!"<<G4endl;
#else
    cout << "PDBlib::load >> file "<<filename<<" not found !!!!"<<endl;
#endif
  }
  else
  {
    int nbAtomTot=0;        // total number of atoms
    int nbAtom=0;           // total number of atoms in a residue
    int numAtomInRes=0;     // number of an atom in a residue
    int nbResidue=0;        // total number of residues
    int nbMolecule=0;       // total number of molecule

    Atom * AtomicOld = NULL;
    Atom * AtomicNext = NULL;

    int serial;         //Atom serial number
    string atomName;    //Atom name
    string element;     //Element Symbol
    string resName;     //Residue name for this atom
    double x,y,z;        //Orthogonal coordinates in Angstroms
    double occupancy;    //Occupancy

    Residue * residueOld = NULL;
    Residue * residueFirst = NULL;
    Residue * residueNext = NULL;

    Molecule * moleculeFirst = NULL;
    Molecule * moleculeOld = NULL;
    Molecule * moleculeNext = NULL;

    /////////////////////////////////
    //NEW variable to draw a fitting cylinder if z oriented
    //=> fitting box
    double minGlobZ,maxGlobZ;
    double minGlobX,maxGlobX;
    double minGlobY,maxGlobY;
    double minX,maxX,minY,maxY,minZ,maxZ; //Sort of 'mother volume' box

#ifdef GEANT4
    minGlobZ=-2000000000;
    minGlobX=-2000000000;
    minGlobY=-2000000000;
    maxGlobZ=2000000000;
    maxGlobX=2000000000;
    maxGlobY=2000000000;
    minX=-2000000000;
    minY=-2000000000;
    minZ=-2000000000;
    maxX=2000000000;
    maxY=2000000000;
    maxZ=2000000000;
#else
    minGlobZ=numeric_limits<double>::min();
    minGlobX=numeric_limits<double>::min();
    minGlobY=numeric_limits<double>::min();
    maxGlobZ=numeric_limits<double>::max();
    maxGlobX=numeric_limits<double>::max();
    maxGlobY=numeric_limits<double>::max();
    minX=numeric_limits<double>::min();
    minY=numeric_limits<double>::min();
    minZ=numeric_limits<double>::min();
    maxX=numeric_limits<double>::max();
    maxY=numeric_limits<double>::max();
    maxZ=numeric_limits<double>::max();
#endif

    int lastResSeq=-1;
    int resSeq=0;

    string nameMol;

    unsigned short int numModel=0;
    unsigned short int model=0;

    //Number of TER
    int ter=0;
    //Number of TER (chain) max for this file
#ifdef GEANT4
    int terMax=2000000000;
#else
    int terMax=numeric_limits<int>::max();
#endif
    if (!infile.eof())
    {
      getline(infile, sLine);
      std::size_t found = sLine.find("DNA");
      if (found!=std::string::npos)
      {
        terMax=2;
        isDNA=1;
      }
      else
        isDNA=0;
      //If PDB file have not a header line
      found = sLine.find("HEADER");
      if (found==std::string::npos)
      {
        infile.close();
        infile.open(filename.c_str());
        //Specific to Geant4, use std::cout instead
#ifdef GEANT4
        G4cout<<"PDBlib::load >> No header found !!!!"<<G4endl;
#else
        cout<<"PDBlib::load >> No header found !!!!"<<endl;
#endif
      }

    }

    while (!infile.eof())
    {
      getline(infile, sLine);
      //In the case of several models
      if ((sLine.substr(0,6)).compare("NUMMDL") == 0)
      {
      istringstream ((sLine.substr(10,4))) >> numModel;
      }
      if ((numModel > 0) && ( (sLine.substr(0,6)).compare("MODEL ") == 0))
      {
        istringstream ((sLine.substr(10,4))) >> model;
        //////////////////////////////////////////
        if (model > 1) break;
        //////////////////////////////////////////
      }
      //For verification of residue sequence
      if ((sLine.substr(0,6)).compare("SEQRES") == 0)
      {
        //Create list of molecule here
#ifdef GEANT4
        if (verbose > 1) G4cout << sLine << G4endl;
#else
        if (verbose > 1) cout << sLine << endl;
#endif
      }

      //Coordinate section
      if ((numModel > 0) && ( (sLine.substr(0,6)).compare("ENDMDL") == 0))
      {;
      }
      else if ((sLine.substr(0,6)).compare("TER   ") == 0) //3 spaces
      {
        //////////////////////////////////////////
        //Currently retrieve only the two first chains(TER) => two DNA strands
        ter++;
        if (ter > terMax) break;
        //////////////////////////////////////////

        //if (verbose>1) G4cout << sLine << G4endl;
        /************ Begin TER ******************/
        lastResSeq=-1;
        resSeq=0;

        AtomicOld->SetNext(NULL);
        residueOld->SetNext(NULL);
        residueOld->fNbAtom=nbAtom;

        //Molecule creation:
        if (moleculeOld == NULL)
        {
          nameMol=filename;//+numModel
          moleculeOld =new Molecule(nameMol,nbMolecule);
          moleculeOld->SetFirst(residueFirst);
          moleculeOld->fNbResidue=nbResidue;
          moleculeFirst = moleculeOld;
        }
        else
        {
          moleculeNext = new Molecule(nameMol,nbMolecule);
          moleculeOld->SetNext(moleculeNext);
          moleculeNext ->SetFirst(residueFirst);
          moleculeNext ->fNbResidue=nbResidue;
          moleculeOld = moleculeNext;
        }

        nbMolecule++;
        moleculeOld->SetNext(NULL);
        moleculeOld->fCenterX = (int)((minX + maxX) /2.);
        moleculeOld->fCenterY = (int)((minY + maxY) /2.);
        moleculeOld->fCenterZ = (int)((minZ + maxZ) /2.);
        moleculeOld->fMaxGlobZ = maxGlobZ;
        moleculeOld->fMinGlobZ = minGlobZ;
        moleculeOld->fMaxGlobX = maxGlobX;
        moleculeOld->fMinGlobX = minGlobX;
        moleculeOld->fMaxGlobY = maxGlobY;
        moleculeOld->fMinGlobY = minGlobY;

#ifdef GEANT4
        minGlobZ=-2000000000;
        minGlobX=-2000000000;
        minGlobY=-2000000000;
        maxGlobZ=2000000000;
        maxGlobX=2000000000;
        maxGlobY=2000000000;
        minX=-2000000000;
        minY=-2000000000;
        minZ=-2000000000;
        maxX=2000000000;
        maxY=2000000000;
        maxZ=2000000000;
#else
        minGlobZ=numeric_limits<double>::min();
        minGlobX=numeric_limits<double>::min();
        minGlobY=numeric_limits<double>::min();
        maxGlobZ=numeric_limits<double>::max();
        maxGlobX=numeric_limits<double>::max();
        maxGlobY=numeric_limits<double>::max();
        minX=numeric_limits<double>::min();
        minY=numeric_limits<double>::min();
        minZ=numeric_limits<double>::min();
        maxX=numeric_limits<double>::max();
        maxY=numeric_limits<double>::max();
        maxZ=numeric_limits<double>::max();
#endif

        nbAtom=0;
        numAtomInRes=0;
        nbResidue=0;
        AtomicOld = NULL;
        AtomicNext = NULL;
        residueOld = NULL;
        residueFirst = NULL;
        residueNext = NULL;

        ///////////// End TER ///////////////////
      }
      else if ((sLine.substr(0,6)).compare("ATOM  ") == 0)
      {

        /************ Begin ATOM ******************/
        //serial
        istringstream ((sLine.substr(6,5))) >> serial;

        //To be improved
        //atomName :
        atomName=sLine.substr(12,4);
        if (atomName.substr(0,1).compare(" ") == 0 )
          element=sLine.substr(13,1);
        else
          element=sLine.substr(12,1);

        // set Van der Waals radius expressed in Angstrom
        double vdwRadius;
        if(element=="H")
        {
          vdwRadius=1.2;
        }
        else if (element=="C")
        {
          vdwRadius=1.7;
        }
        else if (element=="N")
        {
          vdwRadius=1.55;
        }
        else if (element=="O")
        {
          vdwRadius=1.52;
        }
        else if (element=="P")
        {
          vdwRadius=1.8;
        }
        else if (element=="S")
        {
          vdwRadius=1.8;
        }
        else
        {
#ifdef GEANT4
          G4cout << "Element not recognized : " << element << G4endl;
          G4cout << "Stop now" << G4endl;
#else
          cout << "Element not recognized : " << element << endl;
          cout << "Stop now" << endl;
#endif
          exit(1);
        }

        {
          nbAtomTot++;

          //resName :
          resName=sLine.substr(17,3);
          //resSeq :
          istringstream ((sLine.substr(22,4))) >> resSeq;
          //x,y,z :
          stringstream ((sLine.substr(30,8))) >> x;
          stringstream ((sLine.substr(38,8))) >> y;
          stringstream ((sLine.substr(46,8))) >> z;
          //occupancy :
          occupancy=1.;

          if (minGlobZ  < z) minGlobZ=z;
          if (maxGlobZ > z) maxGlobZ=z;
          if (minGlobX  < x) minGlobX=x;
          if (maxGlobX > x) maxGlobX=x;
          if (minGlobY  < y) minGlobY=y;
          if (maxGlobY > y) maxGlobY=y;
          if (minX  > x) minX=x;
          if (maxX < x) maxX=x;
          if (minY  > y) minY=y;
          if (maxY < y) maxY=y;
          if (minZ  > z) minZ=z;
          if (maxZ < z) maxZ=z;

          //treatment for Atoms:
          if (AtomicOld == NULL)
          {
            AtomicOld = new Atom(serial,
                                 atomName,
                                 "",
                                 0,
                                 0,
                                 x,y,z,
                                 vdwRadius,
                                 occupancy,
                                 0,
                                 element);
            AtomicOld->SetNext(NULL);//If only one Atom inside residue
          }
          else
          {
            AtomicNext = new Atom(serial,
                                  atomName,
                                  "",
                                  0,
                                  0,
                                  x,y,z,
                                  vdwRadius,
                                  occupancy,
                                  0,
                                  element);
            AtomicOld->SetNext(AtomicNext);
            AtomicOld = AtomicNext;
          }
          nbAtom++;
        }//END if (element!="H")
        /****************************Begin Residue************************/
        //treatment for residues:
        if (residueOld == NULL)
        {
#ifdef GEANT4
          if (verbose>2) G4cout << "residueOld == NULL"<<G4endl;
#else
          if (verbose>2) cout << "residueOld == NULL"<<endl;
#endif
          AtomicOld->fNumInRes=0;
          residueOld = new Residue(resName,resSeq);
          residueOld->SetFirst(AtomicOld);
          residueOld->SetNext(NULL);
          residueFirst = residueOld;
          lastResSeq=resSeq;
          nbResidue++;
        }
        else
        {
          if (lastResSeq==resSeq)
          {
            numAtomInRes++;
            AtomicOld->fNumInRes=numAtomInRes;
          }
          else
          {
            numAtomInRes=0;
            AtomicOld->fNumInRes=numAtomInRes;

            residueNext = new Residue(resName,resSeq);
            residueNext->SetFirst(AtomicOld);
            residueOld->SetNext(residueNext);
            residueOld->fNbAtom=nbAtom-1;

            nbAtom=1;
            residueOld = residueNext;
            lastResSeq=resSeq;
            nbResidue++;
          }
        }
        /////////////////////////End Residue////////////
        ///////////// End ATOM ///////////////////
      }//END if Atom
    }//END while (!infile.eof())

    infile.close();
    return moleculeFirst;
  }//END else if (!infile)
  return NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * \brief    Compute barycenters
 * \details  Compute barycenters and its coordinate 
 *                  for nucleotides
 * \param    Molecule *    moleculeList
 * \return   Barycenter *
 */
Barycenter * PDBlib::ComputeNucleotideBarycenters(Molecule * moleculeListTemp)
{
  ///////////////////////////////////////////////////////
  //Placement and physical volume construction from memory
  Barycenter * BarycenterFirst = NULL;
  Barycenter * BarycenterOld = NULL;
  Barycenter * BarycenterNext = NULL;

  //Residue (Base, Phosphate,sugar) list
  Residue *residueListTemp;
  //Atom list inside a residu
  Atom *AtomTemp;

  int k=0;
  int j_old=0;

  while (moleculeListTemp)
  {
    residueListTemp = moleculeListTemp->GetFirst();

    k++;
    int j=0;

    //Check numerotation style (1->n per strand or 1->2n for two strand)
    int correctNumerotationNumber=0;
    if (k==2 && residueListTemp->fResSeq > 1)
    {
      correctNumerotationNumber=residueListTemp->fResSeq;
    }

    while (residueListTemp)
    {
      AtomTemp=residueListTemp->GetFirst();
      j++;

      //Correction consequently to numerotation check
      if (correctNumerotationNumber!=0)
      {
        residueListTemp->fResSeq=residueListTemp->fResSeq-
            correctNumerotationNumber+1;
      }

      //Barycenter computation
      double baryX=0.,baryY=0.,baryZ=0.;
      double baryBaseX=0.,baryBaseY=0.,baryBaseZ=0.;
      double barySugX=0.,barySugY=0.,barySugZ=0.;
      double baryPhosX=0.,baryPhosY=0.,baryPhosZ=0.;
      unsigned short int nbAtomInBase=0;

      for (int i=0 ; i < residueListTemp->fNbAtom ; i++)
      {
        //Compute barycenter of the nucletotide
        baryX+=AtomTemp->fX;
        baryY+=AtomTemp->fY;
        baryZ+=AtomTemp->fZ;
        //Compute barycenters for Base Sugar Phosphat
        if (residueListTemp->fResSeq == 1)
        {
          if (i==0)
          {
            baryPhosX+=AtomTemp->fX;
            baryPhosY+=AtomTemp->fY;
            baryPhosZ+=AtomTemp->fZ;
          }
          else if (i<8)
          {
            barySugX+=AtomTemp->fX;
            barySugY+=AtomTemp->fY;
            barySugZ+=AtomTemp->fZ;
          }
          else
          {
            //hydrogen are placed at the end of the residue in a PDB file
            //We don't want them for this calculation
            if (AtomTemp->fElement!="H"){
              baryBaseX+=AtomTemp->fX;
              baryBaseY+=AtomTemp->fY;
              baryBaseZ+=AtomTemp->fZ;
              nbAtomInBase++;}
          }
        }
        else
        {
          if (i<4)
          {
            baryPhosX+=AtomTemp->fX;
            baryPhosY+=AtomTemp->fY;
            baryPhosZ+=AtomTemp->fZ;
          }
          else if (i<11)
          {
            barySugX+=AtomTemp->fX;
            barySugY+=AtomTemp->fY;
            barySugZ+=AtomTemp->fZ;
          }
          else
          {
            //hydrogen are placed at the end of the residue in a PDB file
            //We don't want them for this calculation
            if (AtomTemp->fElement!="H"){// break;
              baryBaseX+=AtomTemp->fX;
              baryBaseY+=AtomTemp->fY;
              baryBaseZ+=AtomTemp->fZ;
              nbAtomInBase++;}
          }
        }
        AtomTemp=AtomTemp->GetNext();
      }//end of for (  i=0 ; i < residueListTemp->nbAtom ; i++)

      baryX = baryX / (double)residueListTemp->fNbAtom;
      baryY = baryY / (double)residueListTemp->fNbAtom;
      baryZ = baryZ / (double)residueListTemp->fNbAtom;

      if (residueListTemp->fResSeq != 1) //Special case first Phosphate
      {
        baryPhosX = baryPhosX / 4.;
        baryPhosY = baryPhosY / 4.;
        baryPhosZ = baryPhosZ / 4.;
      }
      barySugX = barySugX / 7.;
      barySugY = barySugY / 7.;
      barySugZ = barySugZ / 7.;
      baryBaseX = baryBaseX / (double)nbAtomInBase;
      baryBaseY = baryBaseY / (double)nbAtomInBase;
      baryBaseZ = baryBaseZ / (double)nbAtomInBase;

      //Barycenter creation:
      if (BarycenterOld == NULL)
      {
        BarycenterOld =new Barycenter(j+j_old,baryX,baryY,baryZ, //j [1..n]
                                      baryBaseX,baryBaseY,baryBaseZ,
                                      barySugX,barySugY,barySugZ,
                                      baryPhosX,baryPhosY,baryPhosZ);
        BarycenterFirst = BarycenterOld;
      }
      else
      {
        BarycenterNext = new Barycenter(j+j_old,baryX,baryY,baryZ,
                                        baryBaseX,baryBaseY,baryBaseZ,
                                        barySugX,barySugY,barySugZ,
                                        baryPhosX,baryPhosY,baryPhosZ);
        BarycenterOld -> SetNext(BarycenterNext);
        BarycenterOld = BarycenterNext;
      }

      /////////////////////////////////////////////////
      //distance computation between all atoms inside
      //a residue and the barycenter
      AtomTemp=residueListTemp->GetFirst();
      double dT3Dp;
      double max=0.;
      for (int ii=0 ; ii < residueListTemp->fNbAtom ; ii++)
      {
        dT3Dp = DistanceTwo3Dpoints(AtomTemp->fX,BarycenterOld->fCenterX,
                                    AtomTemp->fY,BarycenterOld->fCenterY,
                                    AtomTemp->fZ,BarycenterOld->fCenterZ);
        BarycenterOld->SetDistance(ii,dT3Dp);
        if (dT3Dp>max) max=dT3Dp;
        AtomTemp=AtomTemp->GetNext();
      }//end of for (  i=0 ; i < residueListTemp->nbAtom ; i++)

      BarycenterOld->SetRadius(max+1.8);
      residueListTemp=residueListTemp->GetNext();

    }//end of while sur residueListTemp

    j_old+=j;

    ///molecs->push_back(*moleculeListTemp);
    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while sur moleculeListTemp

  if(BarycenterNext!=NULL)
  {
    BarycenterNext -> SetNext(NULL);
  }

  return BarycenterFirst;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * \brief    the corresponding bounding volume parameters
 * \details  the corresponding bounding volume parameters
 *            to build a box from atoms coordinates
 */
void PDBlib::ComputeBoundingVolumeParams(Molecule *moleculeListTemp,
    double &dX,double &dY,double &dZ,  //Dimensions for bounding volume
    double &tX,double &tY,double &tZ)  //Translation for bounding volume
{
  double minminX,minminY,minminZ;    //minimum minimorum
  double maxmaxX,maxmaxY,maxmaxZ;    //maximum maximorum

#ifdef GEANT4
  minminX=2000000000;
  minminY=2000000000;
  minminZ=2000000000;
  maxmaxX=-2000000000;
  maxmaxY=-2000000000;
  maxmaxZ=-2000000000;
#else
  minminX=numeric_limits<double>::max();
  minminY=numeric_limits<double>::max();
  minminZ=numeric_limits<double>::max();
  maxmaxX=numeric_limits<double>::min();
  maxmaxY=numeric_limits<double>::min();
  maxmaxZ=numeric_limits<double>::min();
#endif

  while (moleculeListTemp)
  {
    if (minminX > moleculeListTemp->fMaxGlobX)
    {
      minminX = moleculeListTemp->fMaxGlobX;
    }
    if (minminY > moleculeListTemp->fMaxGlobY)
    {
      minminY = moleculeListTemp->fMaxGlobY;
    }
    if (minminZ > moleculeListTemp->fMaxGlobZ)
    {
      minminZ = moleculeListTemp->fMaxGlobZ;
    }
    if (maxmaxX < moleculeListTemp->fMinGlobX)
    {
      maxmaxX = moleculeListTemp->fMinGlobX;
    }
    if (maxmaxY < moleculeListTemp->fMinGlobY)
    {
      maxmaxY = moleculeListTemp->fMinGlobY;
    }
    if (maxmaxZ < moleculeListTemp->fMinGlobZ)
    {
      maxmaxZ = moleculeListTemp->fMinGlobZ;
    }

    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while sur moleculeListTemp

  dX=(maxmaxX-minminX)/2.+1.8;  //1.8 => size of biggest radius for atoms
  dY=(maxmaxY-minminY)/2.+1.8;
  dZ=(maxmaxZ-minminZ)/2.+1.8;

  tX=minminX+(maxmaxX-minminX)/2.;
  tY=minminY+(maxmaxY-minminY)/2.;
  tZ=minminZ+(maxmaxZ-minminZ)/2.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * \brief    Compute number of nucleotide per strand
 * \details  Compute number of nucleotide per strand 
 *                  for DNA
 */
void PDBlib::ComputeNbNucleotidsPerStrand(Molecule * moleculeListTemp)
{
  Residue *residueListTemp;

  int k=0;
  int j_old=0;

  while (moleculeListTemp)
  {
    residueListTemp = moleculeListTemp->GetFirst();

    k++;
    int j=0;

    while (residueListTemp)
    {
      j++;
      residueListTemp=residueListTemp->GetNext();
    }//end of while sur residueListTemp

    j_old+=j;
    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while sur moleculeListTemp

  fNbNucleotidsPerStrand=j_old/2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * \brief    Compute barycenters
 * \details  Compute barycenters and its coordinate 
 *                  for nucleotides
 * \param    Molecule *    moleculeList
 * \return   Barycenter *
 */
unsigned short int PDBlib::ComputeMatchEdepDNA(Barycenter *BarycenterList,
    Molecule *moleculeListTemp,
    double x, double y,double z,
    int &numStrand, int &numNucleotid, int &codeResidue)
{
  unsigned short int matchFound=0;

  Molecule *mLTsavedPointer = moleculeListTemp;
  Barycenter *BLsavedPointer = BarycenterList;

  short int strandNum=0;//Strand number
  int residueNum=1;//Residue (nucleotide) number
  string baseName;//Base name [A,C,T,G]
  unsigned short int BSP=2;//Base (default value), Sugar, Phosphat

  double smallestDist;//smallest dist Atom <-> edep coordinates
  double distEdepDNA;
  double distEdepAtom;

  //Residue (Base, Phosphate,suggar) list
  Residue *residueListTemp;
  //Atom list inside a residue
  Atom *AtomTemp;

  int k = 0;//Molecule number
  moleculeListTemp = mLTsavedPointer;
  BarycenterList = BLsavedPointer;

  smallestDist=33.0;//Sufficiently large value
  while (moleculeListTemp)
  {
    k++;
    residueListTemp = moleculeListTemp->GetFirst();

    int j = 0;//Residue number

#ifdef GEANT4
    int j_save=2000000000;//Saved res. number if match found
#else
    int j_save = numeric_limits<int>::max();//Saved res. number if match found
#endif

    while (residueListTemp)
    {
      j++;

      if (j - j_save > 2 ) break;

      distEdepDNA=DistanceTwo3Dpoints(x,BarycenterList->fCenterX,
                                      y,BarycenterList->fCenterY,
                                      z,BarycenterList->fCenterZ);
      if (distEdepDNA < BarycenterList->GetRadius())
      {
        //Find the closest atom
        //Compute distance between energy deposited and atoms for a residue
        //if distance <1.8 then match OK but search inside 2 next residues
        AtomTemp=residueListTemp->GetFirst();
        for (int iii=0 ; iii < residueListTemp->fNbAtom ; iii++)
        {

          distEdepAtom=DistanceTwo3Dpoints(x,AtomTemp->GetX(),
                                           y,AtomTemp->GetY(),
                                           z,AtomTemp->GetZ());

          if ((distEdepAtom < AtomTemp->GetVanDerWaalsRadius())
              && (smallestDist > distEdepAtom))
          {
            strandNum=k;

            if (k==2){
              residueNum = fNbNucleotidsPerStrand + 1 -
                  residueListTemp->fResSeq;
            }
            else
            {
              residueNum = residueListTemp->fResSeq;
            }

            baseName = (residueListTemp->fResName)[2];
            if (residueListTemp->fResSeq == 1)
            {
              if (iii == 0) BSP = 0;//"Phosphate"
              else if (iii < 8) BSP = 1;//"Sugar"
              else BSP = 2;//"Base"
            }
            else
            {
              if (iii < 4) BSP = 0;//"Phosphate"
              else if (iii < 11) BSP = 1;//"Sugar"
              else BSP = 2;//"Base"
            }

            smallestDist=distEdepAtom;

#ifdef GEANT4
            int j_max_value=2000000000;
#else
            int j_max_value = numeric_limits<int>::max();
#endif
            if (j_save == j_max_value) j_save = j;
            matchFound = 1;
          }
          AtomTemp=AtomTemp->GetNext();
        }//end of for (  iii=0 ; iii < residueListTemp->nbAtom ; iii++)
      }//end for if (distEdepDNA < BarycenterList->GetRadius())
      BarycenterList=BarycenterList->GetNext();
      residueListTemp=residueListTemp->GetNext();
    }//end of while sur residueListTemp
    moleculeListTemp=moleculeListTemp->GetNext();
  }//end of while sur moleculeListTemp

  numStrand = strandNum;
  numNucleotid = residueNum;
  codeResidue = BSP;

  return matchFound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PDBlib::DistanceTwo3Dpoints(double xA,double xB,double yA,double yB,
    double zA,double zB)
{
  return std::sqrt ( (xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB) );
}
