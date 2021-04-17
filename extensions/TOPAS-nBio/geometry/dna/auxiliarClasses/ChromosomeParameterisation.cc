// Extra Class for TsDNA
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
/// \file ChromosomeParameterisation.cc
/// \brief Implementation of the ChromosomeParameterisation class

#include <ChromosomeParameterisation.hh>

// STD
#include <fstream>

// G4
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4VPhysicalVolume.hh>

using namespace std;
using CLHEP::nanometer;
using CLHEP::degree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeParameterisation::ChromosomeParameterisation(const char* filename):
    G4VPVParameterisation()
{
  ifstream f(filename, ios::in);
  if (!f)
    return;

  fPositions.reserve(100);
  fRotations.reserve(100);

  while (!f.eof())
  {
    double x, y, z, rot;
    f >> x >> y >> z >> rot;
    fPositions.push_back(new G4ThreeVector(x * nanometer,
                                           y * nanometer,
                                           z * nanometer));
    fRotations.push_back(new G4RotationMatrix(0, 0, rot * degree));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeParameterisation::~ChromosomeParameterisation()
{
  unsigned int i;
  for (i = 0; i < fPositions.size(); i++)
    delete fPositions[i];
  for (i = 0; i < fRotations.size(); i++)
    delete fRotations[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeParameterisation::ComputeDimensions(
    G4Tubs& /*rosette*/,
    const G4int /*copyNo*/,
    const G4VPhysicalVolume* /*physVol*/) const
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeParameterisation::ComputeTransformation(
    const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // see passive vs active method to specify
  // this transformation
  physVol->SetTranslation(*fPositions[copyNo]);
  physVol->SetRotation(fRotations[copyNo]);
}
