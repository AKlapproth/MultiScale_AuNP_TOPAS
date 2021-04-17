//
#ifndef TsEmDNAPhysics_h
#define TsEmDNAPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class TsParameterManager;

class TsEmDNAPhysics : public G4VPhysicsConstructor
{
public:
    
    explicit TsEmDNAPhysics(G4int ver=1, const G4String& name="");
    TsEmDNAPhysics(TsParameterManager* pM);
    
    virtual ~TsEmDNAPhysics();
    
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
private:
    TsParameterManager* fPm;
    G4int  verbose;
};

#endif
