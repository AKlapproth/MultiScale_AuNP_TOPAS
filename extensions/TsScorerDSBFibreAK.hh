// Scorer for DSBFibreAK
// Author: Hongyu Zhu & Alexander Klapproth
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

#ifndef TsScorerDSBFibreAK_hh
#define TsScorerDSBFibreAK_hh

#include "TsVNtupleScorer.hh"
#include "TsHitsRecord.hh"
#include "TsDefineDamage.hh"


#include <stdint.h>


class TsScorerDSBFibreAK : public TsVNtupleScorer
{
public:
    TsScorerDSBFibreAK(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                      G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsScorerDSBFibreAK();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
    //void UserHookForEndOfEvent();
    void UserHookForEndOfRun();
    void GetGeoinfo();
    void Analyze( vector<TsHitsRecord*> Hits, G4int EventID);
    G4int CalculateVoxelID (G4ThreeVector position);
    void FindPlaceInChromosome(G4int VoxelID, G4int FibreID,G4int bpIDinFiber, G4int &VoxelNumInSphere,G4int &ChromosomeID, G4int &bpIDinChromosome );

protected:
    void AccumulateEvent();

 
    
private:
    TsParameterManager* fPm;

    G4double fScoringRadius;
    G4double fScoringTransX;
    G4double fScoringTransY;
    G4double fScoringTransZ;

    G4bool fHistoneAsScavenger;
    G4bool fOHDamageWholeDNA;
    
    G4double fPosX;
    G4double fPosY;
    G4double fPosZ;

    G4int hits_counter;

    G4int    fEventID;
    G4int    fVoxelID;
    G4int    fCopyNumber;
    G4int    fChromosomeID;
    G4int    fBasePairID;
    G4double fEnergyDeposit;
    G4double fParticleEnergy;
    G4double fParticleVelocity;
    G4String fVolumeName;
    G4String fParticleName;
    // G4String fComponentName;

    G4double fProbabilityOfOHInteraction;
    G4double fProbabilityOfOHDamage;

    G4double fDamageThreshold;
    G4bool   fUseLinearProbabilitythreshold;
    G4double fLinearProbability_lower_limit;
    G4double fLinearProbability_upper_limit;
    G4int    fDSBSeparation;
    G4bool   fDefineComplexity;
    G4int    fComplexitySeparation;
    G4bool   fExcludeShortFragment;
    G4int    fLowerFragmentDetectionThreshold;
    G4int    fUpperFragmentDetectionThreshold;
    G4bool   fOnlyIncludeDSBinSDD;
    G4double fDosePerExposoure;
    G4int    NumberOfHistoriesInRun;
    

    G4double fDoseInThisExposoure;
    G4bool   fCalZetaBeta;
    G4int    fExposoureID;
    G4double fEdep;
    G4double fTravelDistance;
    G4double fTrackAveragedLET;
    G4double fZetaBeta_sq;
    std::vector<TsHitsRecord*> Hits;

    TsDefineDamage * aDefineDamage;
    std::vector< vector<TsHitsRecord*>>HitsofEvents;
    std::vector< G4double> EdepofEvents;
    std::vector< G4double> TravelDistanceofEvents;

    // Geometry infomation
    G4double fFiberDNAContent;  //bp
    G4double fVoxelDNAContent;  //bp
    G4double fTotalDNAContent;  //bp
    G4double fNucleusMass;      //g
    G4int    fNumofVoxel;
    G4int    fVoxel3Drepeat;
    G4double fVoxelSize;
    G4double fVoxelContainerminXYZ;

    G4int chromosomeNum;
    std::vector<G4int> fChromosomeDNAContent;
    std::vector<G4double> fChromosomeDNAContentSum;

    G4int TotalVoxelNum;
    std::vector<G4int> fVoxelNumInBox;
    std::vector<G4int> fVoxelNumInSphere; 

    std::vector<G4int> fCH_ID;
    std::vector<G4int> fVoxel_ID;


    G4bool WithinDNAVolumme;
    G4bool WithinDirReactionVolumme;
    G4bool WithinIndirReactionVolumme;

    std::vector< G4String> ParticleOcurrances;
    std::vector< G4String> GeometryOcurrances;
    std::vector< G4int> TrackIDOcurrances;
    
};




#endif