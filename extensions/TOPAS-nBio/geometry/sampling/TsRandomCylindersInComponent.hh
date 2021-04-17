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

#ifndef TsRandomCylindersInComponent_hh
#define TsRandomCylindersInComponent_hh

#include "TsVGeometryComponent.hh"

class G4VSolid;

class TsRandomCylindersInComponent : public TsVGeometryComponent
{
public:
	TsRandomCylindersInComponent(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
								 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsRandomCylindersInComponent();
	
	G4VPhysicalVolume* Construct();
	
private:
	G4VSolid* fEnvelope;
	G4bool fCreate;
	G4double fTRMax;
	G4double fTHL;
	G4double fRMax;
	G4double fHL;
	G4double fHLX;
	G4double fHLY;
	G4double fHLZ;
	
	G4String fEnvelopeType;
	G4String fOutputFile;
	G4int    fNumberOfCylinders;
	
};

#endif

