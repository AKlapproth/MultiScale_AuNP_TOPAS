// Component for TsPlasmid
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                  *
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


// Geometry for a circular DNA plasmid, user specifies the number of base pairs in the ring.
// The plasmid consists of a ring (containing all geometric components) and
// boxes arranged within the ring (each containing the base pair component).
// Each DNA segment consists of a spherical base pair and the surrounding sugar backbone (2 quarter spheres)


#include "TsPlasmid.hh"
#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcommand.hh"


TsPlasmid::TsPlasmid(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                     TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsPlasmid::~TsPlasmid()
{
}


void TsPlasmid::ResolveParameters() {
	fNumberOfBasePairs = fPm->GetIntegerParameter(GetFullParmName("NumberOfBasePairs"));
	fRMin = (0.34 * nm) * (fNumberOfBasePairs) / twopi;
	fRMax = fRMin + 2.4 * nm;
}


G4VPhysicalVolume* TsPlasmid::Construct()
{
    BeginConstruction();
    
    //****************************************************************************
    //                             Ring for plasmid (envelope)
    //****************************************************************************
	
	G4Tubs* envelope = new G4Tubs(fName, fRMin, fRMax, 1.2*nm, 0.0*deg, 360*deg);
    fEnvelopeLog = CreateLogicalVolume(envelope);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
    //****************************************************************************
    //                             Boxes for base pairs
    //****************************************************************************

    G4Box* gBox = new G4Box("BasePair", 1.185*nm, 1.185*nm, 0.17*nm);
	G4LogicalVolume* lBox = CreateLogicalVolume(gBox);
    
    //**************************************************************************
    //                 Subcomponent 1: Base pair
    //**************************************************************************
    //Base pair - a cylinder of radius 0.5 nm and length 0.34 nm.
    G4String subComponent1 = "Base";
    G4Tubs* gBp1 = new G4Tubs(subComponent1, 0, 0.5*nm, 0.17*nm, 0.0*deg, 360.0*deg);
    G4LogicalVolume* lBp1 = CreateLogicalVolume(subComponent1, gBp1);
	
    //**************************************************************************
    //                 Subcomponent 1: Sugar phosphate
    //**************************************************************************
    //Phosphodiester group - two sugars each consisting of quarter cylinders
    //The sugars are wrapped around the base pair
    G4String subComponent2 = "Sugar1";
    G4Tubs* gSugarPhosphate1 = new G4Tubs(subComponent2, 0.5*nm, 1.185*nm, 0.17*nm, 0*deg, 90*deg);
    G4LogicalVolume* lSugarPhosphate1 = CreateLogicalVolume(subComponent2, gSugarPhosphate1);
	
    G4String subComponent3 = "Sugar2";
    G4Tubs* gSugarPhosphate2 = new G4Tubs(subComponent3, 0.5*nm, 1.185*nm, 0.17*nm, 180*deg, 90*deg);
    G4LogicalVolume* lSugarPhosphate2 = CreateLogicalVolume(subComponent3, gSugarPhosphate2);
	
    //*************************************************************************
    //              Place components in ring
    //*************************************************************************
	
    for (G4int j = 0; j < fNumberOfBasePairs ; j++) {
        G4double phi;
        phi = 360.*deg/fNumberOfBasePairs*j;
		
        G4RotationMatrix rotm  = G4RotationMatrix();
        rotm.rotateX(90*deg);
        rotm.rotateY(36*deg*j);
        rotm.rotateZ(phi);
        G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
        G4ThreeVector position = (fRMin+1.2*nm)*uz;
        G4Transform3D transform = G4Transform3D(rotm,position);
        
        //place boxes within the ring
        G4String aName = fName + "_" + G4UIcommand::ConvertToString(j);
		G4VPhysicalVolume* pBasePair = CreatePhysicalVolume("BasePair",
															j, true, lBox,
															transform,
															fEnvelopePhys);
		if ( j == 0 ) {
		  CreatePhysicalVolume("Base", lBp1, pBasePair);
		  CreatePhysicalVolume("Sugar1", lSugarPhosphate1, pBasePair);
		  CreatePhysicalVolume("Sugar2", lSugarPhosphate2, pBasePair);
		}
    }
    
    InstantiateChildren(fEnvelopePhys);
    
    return fEnvelopePhys;
}



