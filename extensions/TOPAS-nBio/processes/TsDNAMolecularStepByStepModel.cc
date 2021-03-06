// Extra Class for TsEmDNAChemistry
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TsDNAMolecularStepByStepModel.cc 94218 2015-11-09 08:24:48Z gcosmo $
//

#include "TsDNAMolecularStepByStepModel.hh"
#include "TsDNASmoluchowskiReactionModel.hh"

#include "globals.hh"
#include "G4DNAMolecularReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Molecule.hh"
#include "G4ReferenceCast.hh"
#include "G4VDNAReactionModel.hh"

TsDNAMolecularStepByStepModel::TsDNAMolecularStepByStepModel(const G4String& name) :
    G4VITStepModel(name),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
  fpTimeStepper = new G4DNAMoleculeEncounterStepper();
  fpReactionProcess = new G4DNAMolecularReaction();

  fType1 = G4Molecule::ITType();
  fType2 = G4Molecule::ITType();
  fReactionModel = 0;
}

TsDNAMolecularStepByStepModel::~TsDNAMolecularStepByStepModel()
{
  if(fReactionModel) delete fReactionModel;
}

TsDNAMolecularStepByStepModel& TsDNAMolecularStepByStepModel::operator=(const TsDNAMolecularStepByStepModel& right)
{
  G4ExceptionDescription exceptionDescription("Use copy constructor rather than assignement operator.");
  G4Exception("TsDNAMolecularStepByStepModel::operator=(const TsDNAMolecularStepByStepModel&)",
              "TsDNAMolecularStepByStepModel001",
              FatalErrorInArgument,
              exceptionDescription);

  if(&right == this) return *this;
  return *this; // avoid warnings
}

TsDNAMolecularStepByStepModel::TsDNAMolecularStepByStepModel(const TsDNAMolecularStepByStepModel& right) :
    G4VITStepModel(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
  fpReactionTable = right.fpReactionTable;
  if(right.fReactionModel)
  {
    fReactionModel = right.fReactionModel->Clone();
    ((G4DNAMolecularReaction*) fpReactionProcess)->SetReactionModel(fReactionModel);
    ((G4DNAMoleculeEncounterStepper*) fpTimeStepper)->SetReactionModel(fReactionModel);
  }
  else fReactionModel = 0;
}

void TsDNAMolecularStepByStepModel::Initialize()
{
  if(fpReactionTable == 0)
  {
    SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
  }

  if(fReactionModel == 0)
  {
    fReactionModel = new TsDNASmoluchowskiReactionModel();
  }

  fReactionModel->SetReactionTable((const G4DNAMolecularReactionTable*) fpReactionTable);

  ((G4DNAMolecularReaction*) fpReactionProcess)->SetReactionModel(fReactionModel);
  ((G4DNAMoleculeEncounterStepper*) fpTimeStepper)->SetReactionModel(fReactionModel);

  G4VITStepModel::Initialize();
}

void TsDNAMolecularStepByStepModel::PrintInfo()
{
#ifdef G4VERBOSE
  G4cout << "DNAMolecularStepByStepModel will be used" << G4endl;
#endif
}
