//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreMoleculeTuple_hh
#define TsScoreMoleculeTuple_hh

#include "TsVNtupleScorer.hh"

class TsScoreMoleculeTuple : public TsVNtupleScorer
{
public:
	TsScoreMoleculeTuple(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	
	virtual ~TsScoreMoleculeTuple();
	
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	
private:
	G4int    fEvt;
	G4String fParticleName;
	G4String fProcessName;
	G4String fVolumeName;
	G4int    fVolumeCopyNumber;
	G4float  fPosX;
	G4float  fPosY;
	G4float  fPosZ;
	G4float  fKineticEnergy;
	G4float  fEnergyDeposited;
	G4float  fTime;
	G4int    fTrackID;
	G4int    fParentAID;
	G4int    fParentBID;
	G4int    fMoleculeID;
	G4int    fStepNumber;
	G4float  fVertexPositionX;
	G4float  fVertexPositionY;
	G4float  fVertexPositionZ;
	
	G4double fTimeCut;
	
private:
	G4bool fIncludeChemistry;
	G4bool fIncludePhysics;
	G4bool fIncludeKineticEnergy;
	G4bool fIncludeEventID;
	G4bool fIncludeTrackID;
	G4bool fIncludeParentID;
	G4bool fIncludeStepNumber;
	G4bool fIncludeParticleName;
	G4bool fIncludeProcessName;
	G4bool fIncludeVolumeName;
	G4bool fIncludeVolumeCopyNumber;
	G4bool fIncludeGlobalTime;
	G4bool fIncludeEnergyDeposited;
	G4bool fIncludeVertex;

};
#endif

