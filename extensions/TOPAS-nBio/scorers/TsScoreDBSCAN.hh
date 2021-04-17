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

#ifndef TsScoreDBSCAN_hh
#define TsScoreDBSCAN_hh

#include "TsVNtupleScorer.hh"

#include <map>

class ClusteringAlgo;

class TsScoreDBSCAN : public TsVNtupleScorer
{
public:
    TsScoreDBSCAN(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsScoreDBSCAN();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    std::vector<ClusteringAlgo*> fVClustering;
    G4int fNbOfAlgo;
    
    G4int fEventID;
    G4int fSSB;
    G4int fDSB;
    G4int fCSB;
    G4int fSizeX;
    G4int fSizeY;
    
};
#endif
