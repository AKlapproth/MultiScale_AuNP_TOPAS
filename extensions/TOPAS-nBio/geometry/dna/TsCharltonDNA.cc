// Component for TsCharltonDNA
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

// DNA segment consisting of a central cylinder (base pair) and the sugarphosphate backbone
// The DNA backbone consists of two half cylinders wrapped around the base pair
// The sugar phosphate cylinders are rotated by 30 degrees to imitate the double helix structure of DNA
// Model is based on that in Charlton, Nikjoo & Humm (1989) Int J Radiat Biol 56(1), 1-19

#include "TsCharltonDNA.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"

#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

#include "G4Box.hh"

TsCharltonDNA::TsCharltonDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsCharltonDNA::~TsCharltonDNA()
{
}


G4VPhysicalVolume* TsCharltonDNA::Construct()
{
	BeginConstruction();
    
    G4double BoxHLX = fPm->GetDoubleParameter(GetFullParmName("BoxHLX"), "Length");
    G4double BoxHLY = fPm->GetDoubleParameter(GetFullParmName("BoxHLY"), "Length");
    G4double BoxHLZ = fPm->GetDoubleParameter(GetFullParmName("BoxHLZ"), "Length");
    
    //****************************************************************************
    //              Envelope geometry: Box containing entire strand
    //****************************************************************************
	
    G4Box* envelopeSolid = new G4Box("box", BoxHLX, BoxHLY, BoxHLZ);
    fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    //****************************************************************************
    //              Subcomponent 1: DNA base-pair
    //****************************************************************************
    
    //Dimensions are taken from Charlton et al. (1989)
    
    G4String subComponentName1 = "BasePair";
    G4Tubs* gBasePair = new G4Tubs("basepair", 0, 0.5*nm, 1.02*nm, 0.0*deg, 360.0*deg);
    G4LogicalVolume* lBasePair = CreateLogicalVolume(gBasePair);
    G4VPhysicalVolume* pBasePair = CreatePhysicalVolume(subComponentName1, lBasePair, fEnvelopePhys);
    
    //****************************************************************************
    //              Subcomponent 2: DNA sugar-phosphate backbone strands
    //****************************************************************************
    
    //Strand 1
    
    G4String subComponentName2 = "SugarPhosphate";
    G4Tubs* gSugarPhosphate = new G4Tubs("sugarphosphate", 0.5*nm, 1.15*nm, 0.17*nm, 0*deg, 180*deg);
    G4LogicalVolume* lSugarPhosphate = CreateLogicalVolume(gSugarPhosphate);
    
    
    //Strand 2
    
    G4String subComponentName3 = "SugarPhosphateStrand2_";
    G4Tubs* gSugarPhosphate2 = new G4Tubs("sugarphosphatestrand2_", 0.5*nm, 1.15*nm, 0.17*nm, 180*deg, 180*deg);
    G4LogicalVolume* lSugarPhosphate2 = CreateLogicalVolume(gSugarPhosphate2);
    
    
    //****************************************************************************
    //              Rotation of strands around the base pair
    //****************************************************************************
    
    G4double x = 0.0;
    G4double y = 0.0;
    G4double z0 = -1.02*nm + 0.17*nm;
    
    for (int j = 0; j < 6; j++){
        
        G4double theta = 36*deg*j;
        G4double z = z0 + j*0.34*nm;
        
        G4ThreeVector* position = new G4ThreeVector(x, y, z);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot ->rotateZ(theta);
        
        G4VPhysicalVolume* pSugarPhosphate = CreatePhysicalVolume(subComponentName2, j, true, lSugarPhosphate, rot, position, fEnvelopePhys);
        CreatePhysicalVolume(subComponentName3, j, true, lSugarPhosphate2, rot, position, fEnvelopePhys);
    }
    
    //Visualization attribute for each component:
    
    G4VisAttributes* SugarVis = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    SugarVis->SetForceSolid(true);
    lSugarPhosphate->SetVisAttributes(SugarVis);
    
    G4VisAttributes* BasePairVis = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    BasePairVis->SetForceSolid(true);
    lBasePair->SetVisAttributes(BasePairVis);
    
    G4VisAttributes* BasePairVis1 = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    BasePairVis1->SetForceSolid(true);
    lSugarPhosphate2->SetVisAttributes(BasePairVis1);
    
    G4VisAttributes* BoxVis = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    BoxVis->SetVisibility(false);
    fEnvelopeLog->SetVisAttributes(BoxVis);
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
