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
// --------------------------------------------------------------
// Authors: E. Delage
// november 2013
// --------------------------------------------------------------
//
// $Id$
//
/// \file PDBmolecule.hh
/// \brief Definition of the PDBmolecule class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MOLECULE_H
#define MOLECULE_H

#include "PDBresidue.hh"

#include <iostream>

using namespace std;

class Residue;

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Molecule
{
public:
  //! First constructor
  Molecule();
  //! Second constructor
  Molecule(string resName,int mNum);
  //! Destructor
  ~Molecule() {};

  //! information about molecule (not implemented)
  //void PrintInfo();
  //! Get the next molecule
  Molecule *GetNext();
  //! Get the first Residue
  Residue *GetFirst();
  //! Get number Molecule
  int GetID();
  //! Set the next Molecule
  void SetNext(Molecule *);
  //! Set the first Residue
  void SetFirst(Residue *);

  string fMolName;   //!< Molecule name
  int fMolNum;       //!< Molecule number

  double fMinGlobZ;   //Cylinder length => min Z
  double fMaxGlobZ;
  double fMinGlobX;   //Radius => min X
  double fMaxGlobX;
  double fMinGlobY;   //=> min Y
  double fMaxGlobY;

  int fCenterX;      //!< "X center" of this Molecule (for rotation...)
  int fCenterY;      //!< "Y center" of this Molecule (for rotation...)
  int fCenterZ;//!< "Z center" of this Molecule (for rotation...)
  int fDistCenterMax;//!< dist from center to most away most of the molecule
  int fNbResidue;        //!< Number of residue inside the molecule

private:
  Molecule *fpNext;//!< Header of the next Molecule (usage before vector)
  Residue *fpFirst;//!< Header of the first Residue (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
