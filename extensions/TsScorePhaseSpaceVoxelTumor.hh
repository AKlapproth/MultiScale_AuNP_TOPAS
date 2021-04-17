//
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

#ifndef TsScorePhaseSpaceVoxelTumor_hh
#define TsScorePhaseSpaceVoxelTumor_hh

#include "TsVNtupleScorer.hh"

#include <stdint.h>

class TsScorePhaseSpaceVoxelTumor : public TsVNtupleScorer
{
public:
	TsScorePhaseSpaceVoxelTumor(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScorePhaseSpaceVoxelTumor();

	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);

    void ReadTumorParameters();

protected:
	void Output();
	void Clear();

	// Output variables
	G4int fNbOfRepeats;
	G4float fPosX;
	G4float fPosY;
	G4float fPosZ;
	G4float fCosX;
	G4float fCosY;
	G4bool fCosZIsNegative;
	G4float fEnergy;
	G4float fWeight;
	G4int fPType;
	G4bool fIsNewHistory;
	G4int fRunID;
	G4int fEventID;
	G4int fTrackID;
	G4int fParentID;
	G4float fCharge;
	G4String fCreatorProcess;
	G4float fVertexKE;
	G4float fVertexPosX;
	G4float fVertexPosY;
	G4float fVertexPosZ;
	G4float fVertexCosX;
	G4float fVertexCosY;
	G4float fVertexCosZ;
	G4float fTOPASTime;
	G4float fTimeOfFlight;
	G4int fSeedPart1;
	G4int fSeedPart2;
	G4int fSeedPart3;
	G4int fSeedPart4;
	G4float fSignedEnergy;
	int8_t fSignedPType;


	// Other properites
	G4bool fKillAfterPhaseSpace;
	G4bool fIncludeCreatorProcess;
	G4bool fIncludeSeed;
	G4bool fOutputToLimited;

	G4int fPrevRunID;
	G4int fPrevEventID;

	G4long fNumberOfHistoriesThatMadeItToPhaseSpace;
	std::map<G4int, G4long> fNumberOfParticles; // # of particles of each type
	std::map<G4int, G4double> fMinimumKE;  // min KE of each particle type
	std::map<G4int, G4double> fMaximumKE;  // max KE of each particle type

	G4float  fPrePosX,fPrePosY,fPrePosZ;
	G4float  fShiftPosX,fShiftPosY,fShiftPosZ;
	G4float  fPostPosX,fPostPosY,fPostPosZ;

	G4String fMaterialName;

	G4double* TumorCtr;
    G4ThreeVector TumorCenter;
    G4double TumorRadX,TumorRadY,TumorRadZ;
};

#endif
