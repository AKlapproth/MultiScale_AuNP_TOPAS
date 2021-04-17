// Scorer for TsScorePhaseSpaceVoxelTumor
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

#include "TsScorePhaseSpaceVoxelTumor.hh"

#include "TsVGeometryComponent.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"
#include "G4PhysicalConstants.hh"

TsScorePhaseSpaceVoxelTumor::TsScorePhaseSpaceVoxelTumor(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
									 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
, fPrePosX(0), fPrePosY(0), fPrePosZ(0), fCosX(0), fCosY(0), fCosZIsNegative(0), fEnergy(0), fWeight(0), fPType(0)
, fIsNewHistory(0), fRunID(0), fEventID(0), fTrackID(0), fParentID(0), fCharge(0), fCreatorProcess(""), fVertexKE(0)
, fVertexPosX(0), fVertexPosY(0), fVertexPosZ(0), fVertexCosX(0), fVertexCosY(0), fVertexCosZ(0)
, fTOPASTime(0), fTimeOfFlight(0), fSeedPart1(0), fSeedPart2(0), fSeedPart3(0), fSeedPart4(0)
, fSignedEnergy(0), fSignedPType(0)
, fPosX(0), fPosY(0), fPosZ(0)
, fKillAfterPhaseSpace(false), fIncludeCreatorProcess(false), fIncludeSeed(false), fOutputToLimited(false)
, fPrevRunID(-1), fPrevEventID(-1), fNumberOfHistoriesThatMadeItToPhaseSpace(0)
{
	SetUnit("");

	if (fPm->ParameterExists(GetFullParmName("KillAfterPhaseSpace")))
		fKillAfterPhaseSpace = fPm->GetBooleanParameter(GetFullParmName("KillAfterPhaseSpace"));

	fNbOfRepeats = 1;
	if (fPm->ParameterExists(GetFullParmName("NbOfRepeats")))
        fNbOfRepeats  = fPm->GetIntegerParameter(GetFullParmName("NbOfRepeats"));

	if (fOutFileType == "limited")
		fOutputToLimited = true;

	if (!fOutputToLimited) {
		// Mandatory data output
		fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
		fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
		fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
		fNtuple->RegisterColumnF(&fCosX, "Direction Cosine X", "");
		fNtuple->RegisterColumnF(&fCosY, "Direction Cosine Y", "");
		fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
		fNtuple->RegisterColumnF(&fWeight, "Weight", "");
		fNtuple->RegisterColumnI(&fPType, "Particle Type (in PDG Format)");
		fNtuple->RegisterColumnB(&fCosZIsNegative, "Flag to tell if Third Direction Cosine is Negative (1 means true)");
		fNtuple->RegisterColumnB(&fIsNewHistory, "Flag to tell if this is the First Scored Particle from this History (1 means true)");

		// Optional data output
		if (fPm->ParameterExists(GetFullParmName("IncludeTOPASTime")) && fPm->GetBooleanParameter(GetFullParmName("IncludeTOPASTime")))
			fNtuple->RegisterColumnF(&fTOPASTime, "TOPAS Time", "s");
		if (fPm->ParameterExists(GetFullParmName("IncludeTimeOfFlight")) && fPm->GetBooleanParameter(GetFullParmName("IncludeTimeOfFlight")))
			fNtuple->RegisterColumnF(&fTimeOfFlight, "Time of Flight", "ns");
		if (fPm->ParameterExists(GetFullParmName("IncludeRunID")) && fPm->GetBooleanParameter(GetFullParmName("IncludeRunID")))
			fNtuple->RegisterColumnI(&fRunID, "Run ID");
		if (fPm->ParameterExists(GetFullParmName("IncludeEventID")) && fPm->GetBooleanParameter(GetFullParmName("IncludeEventID")))
			fNtuple->RegisterColumnI(&fEventID, "Event ID");
		if (fPm->ParameterExists(GetFullParmName("IncludeTrackID")) && fPm->GetBooleanParameter(GetFullParmName("IncludeTrackID")))
			fNtuple->RegisterColumnI(&fTrackID, "Track ID");
		if (fPm->ParameterExists(GetFullParmName("IncludeParentID")) && fPm->GetBooleanParameter(GetFullParmName("IncludeParentID")))
			fNtuple->RegisterColumnI(&fParentID, "Parent ID");
		if (fPm->ParameterExists(GetFullParmName("IncludeCharge")) && fPm->GetBooleanParameter(GetFullParmName("IncludeCharge")))
			fNtuple->RegisterColumnF(&fCharge, "Charge", "e+");
		if (fPm->ParameterExists(GetFullParmName("IncludeCreatorProcess")) && fPm->GetBooleanParameter(GetFullParmName("IncludeCreatorProcess"))) {
			fIncludeCreatorProcess = true;
			fNtuple->RegisterColumnS(&fCreatorProcess, "Creator Process Name");
		}
		if (fPm->ParameterExists(GetFullParmName("IncludeVertexInfo")) && fPm->GetBooleanParameter(GetFullParmName("IncludeVertexInfo"))) {
			fNtuple->RegisterColumnF(&fVertexKE, "Initial Kinetic Energy", "MeV");
			fNtuple->RegisterColumnF(&fVertexPosX, "Vertex Position X", "cm");
			fNtuple->RegisterColumnF(&fVertexPosY, "Vertex Position Y", "cm");
			fNtuple->RegisterColumnF(&fVertexPosZ, "Vertex Position Z", "cm");
			fNtuple->RegisterColumnF(&fVertexCosX, "Initial Direction Cosine X", "");
			fNtuple->RegisterColumnF(&fVertexCosY, "Initial Direction Cosine Y", "");
			fNtuple->RegisterColumnF(&fVertexCosZ, "Initial Direction Cosine Z", "");
		}
		if (fPm->ParameterExists(GetFullParmName("IncludeSeed")) && fPm->GetBooleanParameter(GetFullParmName("IncludeSeed"))) {
			fIncludeSeed = true;
			fNtuple->RegisterColumnI(&fSeedPart1, "Seed Part 1");
			fNtuple->RegisterColumnI(&fSeedPart2, "Seed Part 2");
			fNtuple->RegisterColumnI(&fSeedPart3, "Seed Part 3");
			fNtuple->RegisterColumnI(&fSeedPart4, "Seed Part 4");
		}
	} else { // limited format
		fNtuple->RegisterColumnI8(&fSignedPType, "Particle Type (sign from z direction)");
		fNtuple->RegisterColumnF(&fSignedEnergy, "Energy (-ve if new history)", "MeV");
		fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
		fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
		fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
		fNtuple->RegisterColumnF(&fCosX, "Direction Cosine X", "");
		fNtuple->RegisterColumnF(&fCosY, "Direction Cosine Y", "");
		fNtuple->RegisterColumnF(&fWeight, "Weight", "");
	}


	ReadTumorParameters();
}


TsScorePhaseSpaceVoxelTumor::~TsScorePhaseSpaceVoxelTumor()
{;}


G4bool TsScorePhaseSpaceVoxelTumor::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	//ResolveSolid(aStep);

	if (aStep->GetTrack()->GetVolume() == 0 || aStep->GetTrack()->GetNextVolume() == 0){
		return false;
	}

	G4StepPoint* theStepPoint = aStep->GetPreStepPoint();
	G4ThreeVector pos = theStepPoint->GetPosition();
	fPrePosX = pos.x();
	fPrePosY = pos.y();
	fPrePosZ = pos.z();
	G4ThreeVector postpos = aStep->GetPostStepPoint()->GetPosition();
	fPostPosX = postpos.x();
	fPostPosY = postpos.y();
	fPostPosZ = postpos.z();
	
	G4TouchableHistory* touchablePost = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
	fMaterialName = touchablePost->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();

	if (fMaterialName == "G4_WATER"){

		G4ThreeVector pos_pre = G4ThreeVector(fPrePosX,fPrePosY,fPrePosZ)/cm;
		G4ThreeVector shift_pre = pos_pre - TumorCenter;

	    G4double elltest_pre = pow(shift_pre[0]/TumorRadX,2) +
	    	pow(shift_pre[1]/TumorRadY,2) +
	    	pow(shift_pre[2]/TumorRadZ,2);

		G4ThreeVector pos_post = G4ThreeVector(fPostPosX,fPostPosY,fPostPosZ)/cm;
		G4ThreeVector shift_post = pos_post - TumorCenter;

	    G4double elltest_post = pow(shift_post[0]/TumorRadX,2) +
	    	pow(shift_post[1]/TumorRadY,2) +
	    	pow(shift_post[2]/TumorRadZ,2);

	    if (!(elltest_post < 1 && elltest_pre >= 1)) return false;
	}
	else return false;

	G4ThreeVector mom       = theStepPoint->GetMomentumDirection();
	G4ThreeVector vertexPos = aStep->GetTrack()->GetVertexPosition();
	G4ThreeVector vertexMom = aStep->GetTrack()->GetVertexMomentumDirection();

	fPType          = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
	fEnergy	        = theStepPoint->GetKineticEnergy();
	fWeight	        = theStepPoint->GetWeight();
	fRunID	        = GetRunID();
	fEventID        = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	fTrackID        = aStep->GetTrack()->GetTrackID();
	fParentID       = aStep->GetTrack()->GetParentID();
	fCharge         = aStep->GetTrack()->GetDynamicParticle()->GetCharge();
	
	if (fIncludeCreatorProcess) {
		const G4VProcess* creatorProcess = aStep->GetTrack()->GetCreatorProcess();
		if (creatorProcess)
			fCreatorProcess = creatorProcess->GetProcessName();
		else
			fCreatorProcess = "Primary";
	}
	
	fVertexKE       = aStep->GetTrack()->GetVertexKineticEnergy();
	fVertexPosX     = vertexPos.x();
	fVertexPosY     = vertexPos.y();
	fVertexPosZ     = vertexPos.z();
	fVertexCosX     = vertexMom.x();
	fVertexCosY     = vertexMom.y();
	fVertexCosZ     = vertexMom.z();
	fTimeOfFlight   = aStep->GetTrack()->GetGlobalTime();
	fTOPASTime      = GetTime();

	if (fIncludeSeed) {
		G4Tokenizer next(GetRandomNumberStatusForThisEvent());
		next();
		next();
		G4String token = next();
		fSeedPart1 = G4UIcommand::ConvertToInt(token);
		token = next();
		fSeedPart2 = G4UIcommand::ConvertToInt(token);
		token = next();
		fSeedPart3 = G4UIcommand::ConvertToInt(token);
		token = next();
		fSeedPart4 = G4UIcommand::ConvertToInt(token);
	}

	// Check if this is a new history
	if (fEventID != fPrevEventID || fRunID != fPrevRunID) {
		fIsNewHistory = true;
		fNumberOfHistoriesThatMadeItToPhaseSpace += fNbOfRepeats;
		fPrevEventID = fEventID;
		fPrevRunID = fRunID;
	} else {
		fIsNewHistory = false;
	}

	if (fOutputToLimited) {
		switch(fPType)
		{
			case 22:
				fPType = 1;  // gamma
				break;
			case 11:
				fPType = 2;  // electron
				break;
			case -11:
				fPType = 3;  // positron
				break;
			case 2112:
				fPType = 4;  // neutron
				break;
			case 2212:
				fPType = 5;  // proton
				break;
			default:
				return false;
		}
		fSignedEnergy = fIsNewHistory ? -fEnergy : fEnergy;
		fSignedPType = fCosZIsNegative ? -fPType : fPType;
	}

	// Record some additional statistics
	G4int particleEncoding = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

	G4ThreeVector randomPosition;
	G4ThreeVector randomMomentum;
	if (fNbOfRepeats > 1){
		// works only for sources shooting towards (0,1,0)
		// randomly distrubute points within elliptical source range and straighten momentum
		for (G4int i = 0; i < fNbOfRepeats; i++){ 
			G4bool foundPoint = false;
			G4double randX, randY, randZ, elltest_XZ;
			while(!foundPoint){
				randX = 2*TumorRadX*G4UniformRand() - TumorRadX;
				randZ = 2*TumorRadZ*G4UniformRand() - TumorRadZ;
				elltest_XZ = pow(randX/TumorRadX,2) + pow(randZ/TumorRadZ,2);
				if(elltest_XZ < 1)
					foundPoint = true;
			}
			randY = (-1)*TumorRadY*sqrt(1 - elltest_XZ);

			fPosX           = randX*cm;
			fPosY           = randY*cm;
			fPosZ           = randZ*cm;
			fCosX           = 0;
			fCosY           = 1;
			fCosZIsNegative = false;

			fNtuple->Fill();
			fNumberOfParticles[particleEncoding]++;
		}
	} else {
		fPosX           = pos.x();
		fPosY           = pos.y();
		fPosZ           = pos.z();
		fCosX           = mom.x();
		fCosY           = mom.y();
		fCosZIsNegative = mom.z() < 0.;
		fNtuple->Fill();
		fNumberOfParticles[particleEncoding]++;
	}

	if (fEnergy > fMaximumKE[particleEncoding]) fMaximumKE[particleEncoding] = fEnergy;
	if (fMinimumKE[particleEncoding]==0.0) fMinimumKE[particleEncoding] = fEnergy;
	else if (fEnergy < fMinimumKE[particleEncoding]) fMinimumKE[particleEncoding] = fEnergy;

	if ( fKillAfterPhaseSpace ) aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	return true;
}


// called for each worker at the end of a run
void TsScorePhaseSpaceVoxelTumor::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

	// Absorb additional statistics
	TsScorePhaseSpaceVoxelTumor* workerPhaseSpaceScorer = dynamic_cast<TsScorePhaseSpaceVoxelTumor*>(workerScorer);

	fNumberOfHistoriesThatMadeItToPhaseSpace += workerPhaseSpaceScorer->fNumberOfHistoriesThatMadeItToPhaseSpace;

	std::map<G4int,G4long>::iterator wIter1;
	std::map<G4int,G4long>::iterator mIter1;
	for ( wIter1 = workerPhaseSpaceScorer->fNumberOfParticles.begin(); wIter1 != workerPhaseSpaceScorer->fNumberOfParticles.end(); wIter1++) {
		mIter1 = fNumberOfParticles.find( wIter1->first);
		if (mIter1 == fNumberOfParticles.end())
			fNumberOfParticles.insert(std::pair<G4int, G4long>( wIter1->first, wIter1->second));
		else
			mIter1->second += wIter1->second;
	}

	std::map<G4int,G4double>::iterator wIter2;
	std::map<G4int,G4double>::iterator mIter2;
	for ( wIter2 = workerPhaseSpaceScorer->fMinimumKE.begin(); wIter2 != workerPhaseSpaceScorer->fMinimumKE.end(); wIter2++) {
		mIter2 = fMinimumKE.find( wIter2->first);
		if (mIter2 == fMinimumKE.end())
			fMinimumKE.insert(std::pair<G4int, G4double>( wIter2->first, wIter2->second));
		else if (wIter2->second < mIter2->second)
			mIter2->second = wIter2->second;
	}

	for ( wIter2 = workerPhaseSpaceScorer->fMaximumKE.begin(); wIter2 != workerPhaseSpaceScorer->fMaximumKE.end(); wIter2++) {
		mIter2 = fMaximumKE.find( wIter2->first);
		if (mIter2 == fMaximumKE.end())
			fMaximumKE.insert(std::pair<G4int, G4double>( wIter2->first, wIter2->second));
		else if (wIter2->second > mIter2->second)
			mIter2->second = wIter2->second;
	}

	// Clear additional statistics of worker scorer
	workerPhaseSpaceScorer->fNumberOfHistoriesThatMadeItToPhaseSpace = 0;
	workerPhaseSpaceScorer->fNumberOfParticles.clear();
	workerPhaseSpaceScorer->fMinimumKE.clear();
	workerPhaseSpaceScorer->fMaximumKE.clear();
}


// called for master only at the end of a run
void TsScorePhaseSpaceVoxelTumor::Output()
{
	G4long totalNumberOfParticles = 0;
	std::map<G4int,G4long>::iterator itr;
	for (itr=fNumberOfParticles.begin();itr!=fNumberOfParticles.end();++itr) {
		totalNumberOfParticles += itr->second;
	}

	if (fOutputToLimited) {
		std::ostringstream header;

		header << "$TITLE:" << G4endl;
		header << "TOPAS Phase Space in \"limited\" format. " <<
		"Should only be used when it is necessary to read or write from restrictive older codes." << G4endl;

		header << "$RECORD_CONTENTS:" << G4endl;
		header << "    1     // X is stored ?" << G4endl;
		header << "    1     // Y is stored ?" << G4endl;
		header << "    1     // Z is stored ?" << G4endl;
		header << "    1     // U is stored ?" << G4endl;
		header << "    1     // V is stored ?" << G4endl;
		header << "    1     // W is stored ?" << G4endl;
		header << "    1     // Weight is stored ?" << G4endl;
		header << "    0     // Extra floats stored ?" << G4endl;
		header << "    0     // Extra longs stored ?" << G4endl;

		header << "$RECORD_LENGTH:" << G4endl;
		G4int recordLength = 7*sizeof(G4float) + 1;
		header << recordLength << G4endl;

		header << "$ORIG_HISTORIES:" << G4endl;
		header << GetScoredHistories() << G4endl;

		header << "$PARTICLES:" << G4endl;
		header << totalNumberOfParticles << G4endl;

		header << "$EXTRA_FLOATS:" << G4endl;
		header << "0" << G4endl;

		header << "$EXTRA_INTS:" << G4endl;
		header << "0" << G4endl;

		fNtuple->fHeaderPrefix = header.str();
		fNtuple->SuppressColumnDescription(true);
		fNtuple->Write();
	}


	// Collect prefix statistics
	std::ostringstream prefix;
	prefix << "Number of Original Histories: " << GetScoredHistories() << G4endl;
	prefix << "Number of Original Histories that Reached Phase Space: " << fNumberOfHistoriesThatMadeItToPhaseSpace << G4endl;
	prefix << "Number of Scored Particles: " << totalNumberOfParticles << G4endl;

	// Collect suffix statistics
	std::ostringstream suffix;
	for (itr=fNumberOfParticles.begin();itr!=fNumberOfParticles.end();++itr) {
		if (itr->first)
			suffix << "Number of " << G4ParticleTable::GetParticleTable()->FindParticle(itr->first)->GetParticleName() << ": " << itr->second << G4endl;
		else
			suffix << "Number of particles with PDG code zero: " << itr->second << G4endl;
	}
	suffix << std::endl;
	std::map<G4int,G4double>::iterator itr_d;
	for (itr_d=fMinimumKE.begin();itr_d!=fMinimumKE.end();++itr_d) {
		if (itr_d->first)
			suffix << "Minimum Kinetic Energy of " << G4ParticleTable::GetParticleTable()->FindParticle(itr_d->first)->GetParticleName() << ": " << itr_d->second/MeV << " MeV" << G4endl;
		else
			suffix << "Minimum Kinetic Energy of particles with PDG code zero: " << itr_d->second/MeV << " MeV" << G4endl;
	}
	suffix << std::endl;
	for (itr_d=fMaximumKE.begin();itr_d!=fMaximumKE.end();++itr_d) {
		if (itr_d->first)
			suffix << "Maximum Kinetic Energy of " << G4ParticleTable::GetParticleTable()->FindParticle(itr_d->first)->GetParticleName() << ": " << itr_d->second/MeV << " MeV" << G4endl;
		else
			suffix << "Maximum Kinetic Energy of particles with PDG code zero: " << itr_d->second/MeV << " MeV" << G4endl;
	}

	if (!fOutputToLimited) {
		std::ostringstream title;
		if (fOutFileType == "ascii")
			title << "TOPAS ASCII Phase Space" << G4endl << G4endl;
		else if (fOutFileType == "binary")
			title << "TOPAS Binary Phase Space" << G4endl << G4endl;


		fNtuple->fHeaderPrefix = title.str() + prefix.str();
		fNtuple->fHeaderSuffix = suffix.str();
		fNtuple->Write();
	}

	// report additional statistics
	G4cout << G4endl;
	G4cout << "Scorer: " << GetNameWithSplitId() << G4endl;
	if (fNtuple->HasHeaderFile())
		G4cout << "Header   has been written to file: " << fNtuple->GetHeaderFileName() << G4endl;
	G4cout << "Contents has been written to file: " << fNtuple->GetDataFileName() << G4endl;
	G4cout << "Scored on surface: " << fComponent->GetName() << "/" << GetSurfaceName() << G4endl;
	G4cout << G4endl;
	G4cout << prefix.str() << G4endl;
	G4cout << suffix.str() << G4endl;


	if (fHitsWithNoIncidentParticle > 0) {
		G4cout << "Warning, one of your filters was based on the energy or momentum of the parent particle that was incident on the scoring component," << G4endl;
		G4cout << "but at least one hit resulted from a primary that was already inside the scoring component when it was generated." << G4endl;
		G4cout << "Such hits are always left out by this incident particle filter. Total number of such hits:" << fHitsWithNoIncidentParticle << G4endl;
	}

	UpdateFileNameForUpcomingRun();
}


void TsScorePhaseSpaceVoxelTumor::Clear()
{
	fScoredHistories = 0;
	fNumberOfHistoriesThatMadeItToPhaseSpace = 0;
	fNumberOfParticles.clear();
	fMinimumKE.clear();
	fMaximumKE.clear();
}


void TsScorePhaseSpaceVoxelTumor::ReadTumorParameters()
{
  TumorRadX = 0*cm;
  TumorRadY = 0*cm;
  TumorRadZ = 0*cm;

  TumorCtr = fPm->GetDoubleVector(GetFullParmName("TumorCenterPos"),"Length");
  TumorCenter = G4ThreeVector(TumorCtr[0], TumorCtr[1], TumorCtr[2])/cm;

  TumorRadX = (fPm->GetDoubleParameter(GetFullParmName("TumorRadiusX"),"Length"))/cm;
  TumorRadY = (fPm->GetDoubleParameter(GetFullParmName("TumorRadiusY"),"Length"))/cm;
  TumorRadZ = (fPm->GetDoubleParameter(GetFullParmName("TumorRadiusZ"),"Length"))/cm;
}
