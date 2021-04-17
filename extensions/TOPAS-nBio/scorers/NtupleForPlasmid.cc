// Scorer for NtupleForPlasmid
//
/// ********************************************************************
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

//An example ntuple scorer for the DNA circular plasmid geometry.
//File may be edited to score different quantities.

#include "NtupleForPlasmid.hh"

#include "G4PSDirectionFlag.hh"

#include "TsVGeometryComponent.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

NtupleForPlasmid::NtupleForPlasmid(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //Define ntuple columns

    fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnI(&fVolume, "Volume");
    fNtuple->RegisterColumnI(&fRepNum, "RepNum");
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
    fNtuple->RegisterColumnI(&fProcess, "Process Type");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");

    //fNtuple->RegisterColumnI(&fSectionID, "Section ID");
    //fNtuple->RegisterColumnI(&fStrandID, "Strand ID");
    //fNtuple->RegisterColumnI(&fSegmentID, "Segment ID");
    //fNtuple->RegisterColumnS(&fVolume, "Volume");
}


NtupleForPlasmid::~NtupleForPlasmid() {;}


G4bool NtupleForPlasmid::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    ResolveSolid(aStep);

    fEnergyDep      = aStep->GetTotalEnergyDeposit();

    G4StepPoint* theStepPoint=aStep->GetTrack()->GetStep()->GetPreStepPoint();
    G4TouchableHandle theTouchable = theStepPoint->GetTouchableHandle();
    G4String volumeName = theTouchable->GetVolume()->GetName();

    G4double volumeNumb = theTouchable->GetReplicaNumber();
    G4int VolNum;

    if (volumeName == "BP") VolNum = 1;
    if (volumeName == "Sugar1") VolNum = 2;
    if (volumeName == "Sugar2") VolNum = 3;

    if ((fEnergyDep > 0) && ((volumeName == "Sugar1") || (volumeName == "Sugar2") || (volumeName == "BP"))) {

        //Get position
        G4ThreeVector pos = theStepPoint->GetPosition();
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();

        //Get kinetic energy
        fEnergy = theStepPoint->GetKineticEnergy();

        //Get energy deposition
        fEnergyDep = aStep->GetTrack()->GetStep()->GetTotalEnergyDeposit();

        //Get particle type
        fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

        //Get Volume name and information
        fVolume = VolNum;
        fRepNum = volumeNumb;

        G4String processName = aStep->GetTrack()->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

        //Define process type:
        if (processName=="msc") flagProcess =1;
        if (processName=="e-_G4DNAElastic")		flagProcess =2;
        if (processName=="e-_G4DNAExcitation")		flagProcess =3;
        if (processName=="e-_G4DNAIonisation")		flagProcess =4;
        if (processName=="e-_G4DNAAttachment")		flagProcess =5;
        if (processName=="e-_G4DNAVibExcitation")	flagProcess =6;
        if (processName=="eCapture")			flagProcess =7;

        if (processName=="proton_G4DNAExcitation")	flagProcess =8;
        if (processName=="proton_G4DNAIonisation")	flagProcess =9;
        if (processName=="proton_G4DNAChargeDecrease")	flagProcess =10;

        if (processName=="hydrogen_G4DNAExcitation")	 flagProcess =11;
        if (processName=="hydrogen_G4DNAIonisation")	 flagProcess =12;
        if (processName=="hydrogen_G4DNAChargeIncrease")flagProcess =13;

        if (processName=="alpha_G4DNAExcitation")	flagProcess =14;
        if (processName=="alpha_G4DNAIonisation")	flagProcess =15;
        if (processName=="alpha_G4DNAChargeDecrease")	flagProcess =16;

        if (processName=="alpha+_G4DNAExcitation")	flagProcess =17;
        if (processName=="alpha+_G4DNAIonisation")	flagProcess =18;
        if (processName=="alpha+_G4DNAChargeDecrease")	flagProcess =19;
        if (processName=="alpha+_G4DNAChargeIncrease")	flagProcess =20;

        if (processName=="helium_G4DNAExcitation")	flagProcess =21;
        if (processName=="helium_G4DNAIonisation")	flagProcess =22;
        if (processName=="helium_G4DNAChargeIncrease")	flagProcess =23;

        if (processName=="hIoni")	flagProcess =24;
        if (processName=="eIoni")	flagProcess =25;

        fProcess = flagProcess;

        // Get IDs
        fTrackID = aStep->GetTrack()->GetTrackID();
        fRunID   = GetRunID();
        fEventID = GetEventID();

        fNtuple->Fill();
        return true;
    }
    return false;
}
