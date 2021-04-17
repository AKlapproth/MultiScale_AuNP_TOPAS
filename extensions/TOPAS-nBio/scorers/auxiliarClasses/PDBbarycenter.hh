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
/// \file PDBbarycenter.hh
/// \brief Definition of the Barycenter class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BARY_H
#define BARY_H

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Barycenter
{
public:
  //! First constructor
  Barycenter();
  //! Second constructor
  Barycenter(int bNum,double x,double y, double z,//Nucleotide bar. coordinates
      double Bx,double By, double Bz, //Base bar. coordinates
      double Sx,double Sy, double Sz, //Sugar bar. coordinates
      double Px,double Py, double Pz);//Phosphate bar. coordinates
  //! Destructor
  ~Barycenter() {};

  //! Get the next Barycenter
  Barycenter *GetNext();
  //! Get the first
  //Residue *GetFirst();
  //! Get number Barycenter
  int GetID();
  //! Set the next Barycenter
  void SetNext(Barycenter *);
  //! Set the distance between atom i and nucleotide barycenter
  void SetDistance(int i, double);
  //! Get the distance between atom i and nucleotide barycenter
  double GetDistance(int i);
  //! Set the distance between the farther atom and nucleotide barycenter
  void SetRadius(double );
  //! Get the distance between the farther atom and nucleotide barycenter
  double GetRadius();

  int fBaryNum;//!< Barycenter number
  double fDistanceTab[33];//!< distance table [0..32] (11 hydrogens!)
  double fRadius;

  double fCenterX;            //!< "X coordinate" of this nucelotide Barycenter
  double fCenterY;            //!< "Y coordinate" of this nucelotide Barycenter
  double fCenterZ;            //!< "Z coordinate" of this nucelotide Barycenter

  double fCenterBaseX;        //!< "X coordinate" of this Base Barycenter
  double fCenterBaseY;        //!< "Y coordinate" of this Base Barycenter
  double fCenterBaseZ;        //!< "Z coordinate" of this Base Barycenter

  double fCenterSugarX;       //!< "X coordinate" of this Sugar Barycenter
  double fCenterSugarY;       //!< "Y coordinate" of this Sugar Barycenter
  double fCenterSugarZ;       //!< "Z coordinate" of this Sugar Barycenter

  double fCenterPhosphateX;   //!< "X coordinate" of this Phosphate Barycenter
  double fCenterPhosphateY;   //!< "Y coordinate" of this Phosphate Barycenter
  double fCenterPhosphateZ;   //!< "Z coordinate" of this Phosphate Barycenter

private:
  Barycenter *fpNext;    //!< Header of the next Molecule (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

