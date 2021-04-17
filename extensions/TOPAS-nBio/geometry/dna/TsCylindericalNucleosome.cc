// Component for TsCylindericalNucleosome
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

#include "TsCylindericalNucleosome.hh"

#include "TsParameterManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Tubs.hh"

TsCylindericalNucleosome::TsCylindericalNucleosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsCylindericalNucleosome::~TsCylindericalNucleosome()
{
}


G4VPhysicalVolume* TsCylindericalNucleosome::Construct()
{
	BeginConstruction();

    G4String name = GetFullParmName("NucleosomeRadius");
    G4double NucleosomeRadius;
    if (!fPm->ParameterExists(name)) {
        NucleosomeRadius = 5*nm;
    }
    else{
       NucleosomeRadius = fPm->GetDoubleParameter(name, "Length");
    }
    
    name = GetFullParmName("NucleosomeLength");
    G4double NucleosomeLength;
    if (!fPm->ParameterExists(name)) {
        NucleosomeLength = 2.5*nm;
    }
    else{
        NucleosomeLength = fPm->GetDoubleParameter(name, "Length");
    }
    
	G4Tubs* envelopeSolid = new G4Tubs(fName, 0, NucleosomeRadius, NucleosomeLength, 0.0*deg, 360.0*deg);
	G4LogicalVolume* envelopeLog = CreateLogicalVolume(envelopeSolid);
	G4VPhysicalVolume* envelopePhys = CreatePhysicalVolume(envelopeLog);
    
	InstantiateChildren(envelopePhys);
	
	return envelopePhys;
}
