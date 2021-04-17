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

#ifndef TsScoreChemistry_hh
#define TsScoreChemistry_hh

#include "TsVNtupleScorer.hh"

#include <stdint.h>

class TsScoreChemistry : public TsVNtupleScorer
{
public:
	TsScoreChemistry(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreChemistry();

	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

protected:

	// Output variables
	G4float fPosX;
	G4float fPosY;
	G4float fPosZ;
	G4float fTransX;
	G4float fTransY;
	G4float fTransZ;
	G4int fPType;
	G4int fTrackID;
	G4bool fInNucleus;
	G4double fNucleusRad;
	G4String fParticleName;
	
	// Other properites
	G4bool fKillAfterPhaseSpace;
	G4bool fIncludeNucleusFlag;

	G4int fPrevTrackID;

	G4long fNumberOfHistoriesThatMadeItToPhaseSpace;
	std::map<G4int, G4long> fNumberOfParticles; // # of particles of each type
	std::map<G4int, G4double> fMinimumKE;  // min KE of each particle type
	std::map<G4int, G4double> fMaximumKE;  // max KE of each particle type
};

#endif
