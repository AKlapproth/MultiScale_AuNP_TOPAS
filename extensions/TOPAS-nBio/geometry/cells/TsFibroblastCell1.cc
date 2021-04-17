// Component for TsFibroblastCell1
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:											      *
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


#include "TsFibroblastCell1.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsFibroblastCell1::TsFibroblastCell1(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsFibroblastCell1::~TsFibroblastCell1()
{
}


G4VPhysicalVolume* TsFibroblastCell1::Construct()
{
	BeginConstruction();

    //***********************************************************************
    //              Envelope Geometry : Cell
    //***********************************************************************
    
    //Base of the cell:
    std::vector<G4TwoVector> cellpoly(9);
    cellpoly[0] = G4TwoVector(24*micrometer, 24*micrometer);
    cellpoly[1] = G4TwoVector(15*micrometer, 15*micrometer);
    cellpoly[2] = G4TwoVector(-24*micrometer, 24*micrometer);
    cellpoly[3] = G4TwoVector(-15*micrometer, -10*micrometer);
    cellpoly[4] = G4TwoVector(-24*micrometer,-24*micrometer);
    cellpoly[5] = G4TwoVector(-20*micrometer, -27*micrometer);
    cellpoly[6] = G4TwoVector(3*micrometer, -18*micrometer);
    cellpoly[7] = G4TwoVector(26*micrometer, -26*micrometer);
    cellpoly[8] = G4TwoVector(15*micrometer, 5*micrometer);
    
    // Height of cell (z):
    G4double hz = 10*micrometer;  //half length along z

    // Used for defining the area that mitochondria are distributed in
    G4double CellRadius = 10*micrometer;
    
    G4ExtrudedSolid* gFibroCell = new G4ExtrudedSolid(fName,
                                                      cellpoly,
                                                      hz,
                                                      G4TwoVector(0,0), 0.5, G4TwoVector(0,0), 1.0); // shape along z
    
    G4LogicalVolume* fEnvelopeLog = CreateLogicalVolume(gFibroCell);
    G4VPhysicalVolume* fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    // Nucleus
    G4double NuclRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NuclRadius");
    if (fPm->ParameterExists(name)) {
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
        
        //***************************
        // Subcomponent: Nucleus
        //***************************
        
        G4String subComponentName1 = "Nucleus";
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
        
        G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);
    }
    
    
    //Mitochondria
    name = GetFullParmName("Mitochondria/NbOfMito");
    if (fPm->ParameterExists(name)) {
        
        //number of mitochondria
        const G4int NbOfMito  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NbOfMito") );
        
        //Semi-axis lengths of the ellpsoid
        G4double EllA = 0.5*micrometer;
        G4double EllB = 0.3*micrometer;
        G4double EllC = 0.9*micrometer;
        
        name=GetFullParmName("Mitochondria/a");
        if (fPm->ParameterExists(name)){EllA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/a"), "Length" );}
        
        name=GetFullParmName("Mitochondria/b");
        if (fPm->ParameterExists(name)){EllB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/b"), "Length" );}
        
        name=GetFullParmName("Mitochondria/c");
        if (fPm->ParameterExists(name)){EllC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/c"), "Length" );}
        
        //*******************************
        // Subcomponent: Mitochondria
        //*******************************
        G4String subComponentName2 = "Mitochondria";
        G4Ellipsoid* gMito = new G4Ellipsoid("gMito", EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
    
        //Randomly distribute mitochondria throughout cell volume outside nucleus (default)
        for (int j = 0; j < NbOfMito; j++){
    
            G4bool Overlap = true;
            while (Overlap == true){
        
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(CellRadius - NuclRadius);
                G4double phi = G4UniformRand()*2*pi;
                G4double psi = G4UniformRand()*2*pi;
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
    
                x = (NuclRadius + dr)* std::cos(u) * std::sin(v);
                y = (NuclRadius + dr)* std::sin(u) * std::sin(v);
                z = (NuclRadius + dr)* std::cos(v);

                G4ThreeVector* position = new G4ThreeVector(x,y,z);
    
                G4RotationMatrix* rotm = new G4RotationMatrix();
    
                rotm->rotateX(psi);
                rotm->rotateY(phi);
    
                G4VPhysicalVolume* pMito = CreatePhysicalVolume(subComponentName2, j, true, lMito, rotm, position, fEnvelopePhys);
            
                G4bool OverlapCheck = pMito->CheckOverlaps();
            
                if (OverlapCheck == false){break;}
                if (OverlapCheck == true){delete pMito;}
            }
        }
    
    }
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
