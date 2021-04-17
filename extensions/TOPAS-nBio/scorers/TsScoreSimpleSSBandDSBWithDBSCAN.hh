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

#ifndef TsScoreSimpleSSBandDSBWithDBSCAN_hh
#define TsScoreSimpleSSBandDSBWithDBSCAN_hh

#include "TsVNtupleScorer.hh"

#include <map>

class G4Material;
class ClusteringAlgo;

class TsScoreSimpleSSBandDSBWithDBSCAN : public TsVNtupleScorer
{
public:
    TsScoreSimpleSSBandDSBWithDBSCAN(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsScoreSimpleSSBandDSBWithDBSCAN();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    TsMaterialManager* fmM;
    std::vector<ClusteringAlgo*> fVClustering;
    G4int fBasePairDepth;
    G4Material* fStrand1Material;
    G4Material* fStrand2Material;
    
    G4int fNbOfAlgo;
    
    G4int fEventID;
    G4int fSSB;
    G4int fDSB;
    G4int fCSB;
    G4int fSizeX;
    G4int fSizeY;
    
};
#endif
