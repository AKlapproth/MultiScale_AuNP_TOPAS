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
/// \file PDBlib.hh
/// \brief Definition of the PDBlib class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PDBlib_h
#define PDBlib_h 1

#include "PDBbarycenter.hh"
#include "PDBmolecule.hh"
#include <vector>

//! PDBlib Class
/*!
 * This Class define Molecule model ... 
 */
class PDBlib
{
public:
  //! First constructor
  PDBlib();
  //! Destructor
  ~PDBlib() {};

  //! Load PDB file into memory
  Molecule* Load( const string &filename,
      unsigned short int &isDNA,
      unsigned short int verbose);

  // All declarations below are 'DNA specific'
  // Just comment those lines if you need to use this code elsewhere.

  //! Compute nucleotide barycenter from memory
  Barycenter* ComputeNucleotideBarycenters(Molecule *moleculeListTemp);

  //! Compute the corresponding bounding volume parameters
  void ComputeBoundingVolumeParams(Molecule *moleculeListTemp,
      double &dX,double &dY,double &dZ,       //Dimensions for bounding volume
      double &tX,double &tY,double &tZ);      //Translation for bounding volume

  //! Compute number of nucleotide per strand
  void ComputeNbNucleotidsPerStrand(Molecule * moleculeListTemp);

  //! Compute if energy is deposited in per atom
  unsigned short int ComputeMatchEdepDNA(Barycenter *,Molecule *,
      double x, double y,double z,
      int &numStrand, int &numNucleotid, int &codeResidue);

private:
  //! return distance between two 3D points
  double DistanceTwo3Dpoints(double xA,double xB,
      double yA,double yB,
      double zA,double zB);

  //! Number of nucleotid per strand
  int fNbNucleotidsPerStrand;
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
