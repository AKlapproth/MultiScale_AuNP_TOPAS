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
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file ChromosomeParameterisation.hh
/// \brief Definition of the ChromosomeParameterisation class


#ifndef ChromosomeParameterisation_H
#define ChromosomeParameterisation_H 1

#include <G4VPVParameterisation.hh>
#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4Torus.hh>
#include <G4Trap.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>
#include <G4Ellipsoid.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include <vector>

class G4VPhysicalVolume;
class G4Box;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;
class G4Ellipsoid;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChromosomeParameterisation : public G4VPVParameterisation
{
public:
  ChromosomeParameterisation(const char* filename);

  virtual ~ChromosomeParameterisation();
  inline int GetNumRosettes()
  {
    return fPositions.size();
  }

  virtual void
  ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

  virtual void
  ComputeDimensions(G4Tubs& rosette,
      const G4int copyNo,
      const G4VPhysicalVolume* physVol) const;

private:
  // Dummy declarations to get rid of warnings ...
  virtual void
  ComputeDimensions(G4Box&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Trd&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Trap&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Sphere&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Orb&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Torus&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Para&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Hype&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Polycone&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Polyhedra&, const G4int, const G4VPhysicalVolume*) const
  {
  }
  virtual void
  ComputeDimensions(G4Ellipsoid&, const G4int, const G4VPhysicalVolume*) const
  {
  }

  std::vector<G4ThreeVector*> fPositions;
  std::vector<G4RotationMatrix*> fRotations;
};

#endif
