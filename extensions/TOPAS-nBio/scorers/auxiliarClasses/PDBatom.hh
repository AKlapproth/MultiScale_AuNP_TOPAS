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
/// \file PDBatom.hh
/// \brief Definition of the Atom class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ATOM_H
#define ATOM_H

#include <iostream>
using namespace std;

//! Atom Class
/*!
 * This Class define Atom model ... 
 */
class Atom
{
public:
  //! constructor with initialization
  Atom(int serial,string name,string resName,int numInRes,int resSeq,
      double xInit,double yInit,double zInit,
      double radius,
      double occupancy, double tempFactor, string element);
  //! Empty destructor
  ~Atom()
  {
  };

  //! Returns the next Atom
  Atom *GetNext();
  //! Return the X position for the Atom
  double GetX();
  //! Return the Y position for the Atom
  double GetY();
  //! Return the Z position for the Atom
  double GetZ();
  //! Return the Atom's ID
  int GetID();
  //! Return name of the atom
  string GetName();
  //! Return name of the element
  string GetElementName();
  //! Return name of the atom
  double GetVanDerWaalsRadius();
  //! Set the next atom
  void SetNext(Atom *);

  int fSerial;       //!< its serial number
  int fNumInRes;     //!< its number in residue sequence
  string fName;      //!< Atom name
  string fResName;   //!< Residue name
  int fResSeq;       //!< Residue sequence number
  double fX;          //!< X orthogonal coordinates in Angstroms
  double fY;          //!< Y orthogonal coordinates in Angstroms
  double fZ;          //!< Z orthogonal coordinates in Angstroms
  double fVdwRadius;  // Vand der Waals Radius in Angstrom
  double fOccupancy;  //!< Occupancy for the Atom
  string fElement;   //!< Element symbol extracted from 'atom name'
  double fTempFactor; //!< Temperature factor for the Atom

private:
  Atom * fpNext;       //!< Pointer to the next Atom
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
