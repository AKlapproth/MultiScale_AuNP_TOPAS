//
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
#ifndef TsEmDNAChemistry_hh
#define TsEmDNAChemistry_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include "globals.hh"
#include <map>

class TsParameterManager;

class G4DNAMolecularReactionTable;

class TsEmDNAChemistry
: public G4VUserChemistryList,
  public G4VPhysicsConstructor
{
    
public:
    
    explicit TsEmDNAChemistry();
    TsEmDNAChemistry(TsParameterManager* pM);
    
    virtual ~TsEmDNAChemistry();
    
    virtual void ConstructParticle()
    {
        ConstructMolecule();
    }
    
    virtual void ConstructMolecule();
    virtual void ConstructProcess();
    
    virtual void ConstructDissociationChannels();
    virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable);
    virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable);
    
    G4bool IsWaterMolecularConfigurationActivated();
    
private:
    void DefineParameters();
    G4bool MoleculeExists(G4String);
    G4String GetFullParmName(G4String);
    void Quit(G4String, G4String);
    
private:
    TsParameterManager* fPm;
    std::map<G4String, G4String> fExistingMolecules;
    G4String fName;
    
    std::map<G4String, G4double> fDiffusionCoefficients;
    std::vector< std::vector<G4String> > fReactionSpecies;
    std::vector< std::vector<G4String> > fReactionProducts;
    std::vector<G4double> fReactionRates;
    
    G4bool fSetWaterConfiguration;
    
    G4double fIonizationStates;
    G4double fA1B1DissociativeDecay;
    G4double fA1B1Relaxation;
    G4double fB1A1AutoIonization;
    G4double fB1A1DissociativeDecay;
    G4double fB1A1Relaxation;
    G4double fRydDiffAutoIonization;
    G4double fRydDiffRelaxation;
    G4double fDissociativeAttachment;

    std::vector<G4String> fScavengedMolecules;
    std::vector<G4double> fScavengedCapacities;
    std::vector<G4bool> fScavengerHasProducts;
	std::vector<G4String> fRemoveInMaterialTheseMolecules;
	G4String fRemoveInMaterial;
};

#endif
