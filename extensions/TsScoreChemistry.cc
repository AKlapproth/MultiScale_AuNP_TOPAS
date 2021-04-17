// Scorer for TsScoreChemistry
// Author: Alexander Klapproth
// Last edit by Alexander Klapproth: 04/17/2021
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                 *
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
//

// This scorer is only viable for spherical components!
// Scored particles are duplicated and the duplicates are randomly rotated around the center of the component.

#include "TsScoreChemistry.hh"

#include "TsVGeometryComponent.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsScoreChemistry::TsScoreChemistry(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
									 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
, fPosX(0), fPosY(0), fPosZ(0), fPType(0)
, fTrackID(0)
, fKillAfterPhaseSpace(false), fInNucleus(false)
{

	if (fPm->ParameterExists(GetFullParmName("KillAfterPhaseSpace")))
		fKillAfterPhaseSpace = fPm->GetBooleanParameter(GetFullParmName("KillAfterPhaseSpace"));


    fIncludeNucleusFlag = false;
    if ( fPm->ParameterExists(GetFullParmName("IncludeNucleusFlag")) ) {
        fIncludeNucleusFlag = fPm->GetBooleanParameter(GetFullParmName("IncludeNucleusFlag"));
    }

    // offset of the Nucleus center from the origin of World
	fTransX = 0*mm;
	fTransY = 0*mm;
	fTransZ = 0*mm;
	if (fPm->ParameterExists(GetFullParmName("OffsetX")))
		fTransX = fPm->GetDoubleParameter(GetFullParmName("OffsetX"),"Length");
	if (fPm->ParameterExists(GetFullParmName("OffsetY")))
		fTransY = fPm->GetDoubleParameter(GetFullParmName("OffsetY"),"Length");
	if (fPm->ParameterExists(GetFullParmName("OffsetZ")))
		fTransZ = fPm->GetDoubleParameter(GetFullParmName("OffsetZ"),"Length");

	fNucleusRad = 6.9*um;
	if (fPm->ParameterExists(GetFullParmName("NucleusRadius")))
		fNucleusRad = fPm->GetDoubleParameter(GetFullParmName("NucleusRadius"),"Length");

	fNtuple->RegisterColumnI(&fPType, "Particle Type Flag");

	if (fIncludeNucleusFlag)
		fNtuple->RegisterColumnB(&fInNucleus, "Flag to tell if the history originated inside the Nucleus");
}


TsScoreChemistry::~TsScoreChemistry()
{;}


G4bool TsScoreChemistry::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	//only score chemical histories
	fTrackID = aStep->GetTrack()->GetTrackID();
	if (fTrackID >= 0)
		return false;

	//only score new tracks
	if (aStep->GetTrack()->GetCurrentStepNumber() != 1) {
		return false;
	}

	if (fIncludeNucleusFlag){

		G4StepPoint* theStepPoint=0;
		theStepPoint = aStep->GetPreStepPoint();

		G4ThreeVector pos       = theStepPoint->GetPosition();

		fPosX           = pos.x() - fTransX;
		fPosY           = pos.y() - fTransY;
		fPosZ           = pos.z() - fTransZ;

		G4double NucTest = sqrt(pow(fPosX,2) + pow(fPosY,2) + pow(fPosZ,2));
		fInNucleus = (NucTest < fNucleusRad);
		
	}

	fParticleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

    if (fParticleName == "H2O")
		fPType = 0;
	else if (fParticleName == "H")
		fPType = 1;
	else if (fParticleName == "OH")
		fPType = 2;
	else if (fParticleName == "H_2")
		fPType = 3;
	else if (fParticleName == "H2O2")
		fPType = 4;
	else if (fParticleName == "H3O")
		fPType = 5;
	else if (fParticleName == "e_aq")
		fPType = 6;
	else {
		G4cout << "Unknown chemical particle type: " << fParticleName << G4endl;
		return false;
	}

	fNtuple->Fill();

	if ( fKillAfterPhaseSpace ) aStep->GetTrack()->SetTrackStatus(fStopAndKill);

	return true;
}