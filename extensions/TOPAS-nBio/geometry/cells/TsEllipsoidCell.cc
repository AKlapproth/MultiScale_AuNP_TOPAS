// Component for TsEllipsoidCell
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
//

#include "TsEllipsoidCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsEllipsoidCell::TsEllipsoidCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsEllipsoidCell::~TsEllipsoidCell()
{
}


G4VPhysicalVolume* TsEllipsoidCell::Construct()
{
	BeginConstruction();
    
    //User defined parameters:
    //User can specify the cell semi-axis lengths or use the default values.
    
    //Set default values for cell dimensions
    G4double xSA = 20 * um;
    G4double ySA = 10 * um;
    G4double zSA = 15 * um;
    
    G4String name = GetFullParmName("xSemiAxis");
    if (fPm->ParameterExists(name)) {
        
        //User specified cell size values
        xSA = fPm->GetDoubleParameter(GetFullParmName("xSemiAxis"), "Length");
        ySA = fPm->GetDoubleParameter(GetFullParmName("ySemiAxis"), "Length");
        zSA = fPm->GetDoubleParameter(GetFullParmName("zSemiAxis"), "Length");
        
    }
	
    //***********************************************************************
    //              Envelope Geometry : ellipsoid cell
    //***********************************************************************
    
	G4Ellipsoid* envelopeSolid = new G4Ellipsoid(fName, xSA, ySA, zSA);
	G4LogicalVolume* envelopeLog = CreateLogicalVolume(envelopeSolid);
	G4VPhysicalVolume* fEnvelopePhys = CreatePhysicalVolume(envelopeLog);
    
    //**************************************************************************
    // Optional Components : include a nucleus and/or mitochondria in the cell
    //**************************************************************************
    
    
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    G4double NuclRadius = 0.0*um;
    name = GetFullParmName("Nucleus/NuclRadius");
    if (fPm->ParameterExists(name)) {
        
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
        G4String subComponentName1 = "Nucleus";
        
        G4RotationMatrix* rotNuc = new G4RotationMatrix();
        
        rotNuc->rotateX(0);
        rotNuc->rotateY(0);
        
        
        G4double transNucX = 0 * um;
        G4double transNucY = 0 * um;
        G4double transNucZ = 0 * um;
        
        G4String name1 = GetFullParmName("Nucleus/transNucX");
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
            if (transNucX > xSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/transNucY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
            if (transNucY > ySA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/transNucZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
            if (transNucZ > zSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        G4ThreeVector* NucPos = new G4ThreeVector(transNucX,transNucY,transNucZ);
        
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
        //G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);
        G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, rotNuc, NucPos, fEnvelopePhys);
    }
    
    

    //*******************************
    // Subcomponent: Mitochondria
    //*******************************
    
    name = GetFullParmName("Mitochondria/NbOfMito");
    if (fPm->ParameterExists(name)) {
        
        //number of mitochondria
        const G4int NbOfMito  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NbOfMito") );
        
        //Semi-axis lengths of the ellpsoid/mitochondria (default values if none are specified)
        G4double EllA = 0.5*micrometer;
        G4double EllB = 0.3*micrometer;
        G4double EllC = 0.9*micrometer;
        
        name=GetFullParmName("Mitochondria/a");
        if (fPm->ParameterExists(name)){EllA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/a"), "Length" );}
        
        name=GetFullParmName("Mitochondria/b");
        if (fPm->ParameterExists(name)){EllB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/b"), "Length" );}
        
        name=GetFullParmName("Mitochondria/c");
        if (fPm->ParameterExists(name)){EllC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/c"), "Length" );}
        
        
        G4String subComponentName2 = "Mitochondria";
        G4Ellipsoid* gMito = new G4Ellipsoid("gMito", EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
        
        //Find biggest semi-axis length
        G4double CellRadius;
        if ((xSA >= ySA) && (xSA >= zSA)){CellRadius = xSA;}
        if ((ySA >= xSA) && (ySA >= zSA)){CellRadius = ySA;}
        if ((zSA >= xSA) && (zSA >= ySA)){CellRadius = zSA;}
        
        
        //Randomly distribute mitochondria throughout cell volume
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
