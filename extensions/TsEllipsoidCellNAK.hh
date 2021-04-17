
// Author: Alexander Klapproth
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
//

#ifndef TsEllipsoidCellNAK_hh
#define TsEllipsoidCellNAK_hh

#include "TsVGeometryComponent.hh"
#include "G4PVParameterised.hh"
#include "VoxelParameterisation.hh"
#include "G4NistManager.hh"

using namespace CLHEP;
using namespace std;


class TsEllipsoidCellNAK : public TsVGeometryComponent
{    
public:
	TsEllipsoidCellNAK(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsEllipsoidCellNAK();
	
	G4VPhysicalVolume* Construct();

    void SetBasicInfo();

    void int2str(const int &int_temp,string &string_temp);

private:

    G4bool fNucleusExists;
    G4double NuclRadius;
        
    G4double transNucX;
    G4double transNucY;
    G4double transNucZ;
    
    G4LogicalVolume * FiberEnvelopeLogic;
    G4LogicalVolume * FiberLogic;
    std::vector<G4VPhysicalVolume*> pvLoop;
    std::vector<G4LogicalVolume*> lvLoop;
    std::vector<G4double> fx;
    std::vector<G4double> fy;
    std::vector<G4double> fz;
    
    // control parameters
    G4double FibreRadius;
    G4double FibreLength;
    G4double VoxelLength;
    G4double ParaContainerHalfSize;
    G4int HilbertCurveLayer;
    G4int HilbertCurve3DRepeat;   
    G4bool CheckOverlap;
    G4bool testmode;
    G4bool AddHydrationShell;
    G4double HydrationShellThickness;
    G4double ToleranceSpace;
    
     
    G4double fNucleusMass;
    G4int CountFibers;

    //ParameterisationInfo
    VoxelParameterisation* param;
    G4int fnVoxelX, fnVoxelY, fnVoxelZ, fnVoxels;
    G4double fvoxelDimX, fvoxelDimY, fvoxelDimZ;
    G4double fminX, fmaxX, fminY, fmaxY, fminZ, fmaxZ;
    G4double voxelHalfDimX,  voxelHalfDimY, voxelHalfDimZ;
    std::vector<G4Material*> fMaterials; 

// ---------- AuNPs ----------

    std::vector<G4ThreeVector> FibreMidpoints;
    std::vector<G4String> FibreDirections;
    std::vector<G4ThreeVector> NPpositions;
    std::vector<G4VPhysicalVolume*> NPLoop;
    G4LogicalVolume* lAuNPN;
    G4LogicalVolume* lCore;

// ---------- Fibre ----------

    G4LogicalVolume* Fibre_Construct();

    void Fibre_BuildHistones(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                       G4double ChromatinRadius,
                       G4double ChromatinLength);

    void Fibre_BuildDNA(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails);


    G4LogicalVolume *lBase1, *lBase2, *lBack1, *lBack2;
    G4LogicalVolume *lWater1, *lWater2;
    G4LogicalVolume* aEnvelopeLog ;
    G4VPhysicalVolume* aEnvelopePhys;
    G4Material * water;
    G4Material * mat_DNA_base;
    G4Material * mat_DNA_backbone;
    G4Material * mat_DNA_hydrationShell;
    G4Material * mat_Histone;
    G4bool fAddHydrationShell; 
    G4double fHydrationShellThickness;
    G4double fChromatinRadius; 
    G4double fChromatinLength;

    G4double BuildSphere;



    void Fibre_SetDNAVolumes(G4bool BuildHalfCyl,
                       G4bool BuildQuartCyl,
                       G4bool BuildSphere);
    void Fibre_GenerateDNAPath(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                         vector<G4ThreeVector> &path);
    void Fibre_SegmentDNAPath(vector<G4ThreeVector> &path);
    void Fibre_PlaceDNASphere(vector<G4ThreeVector> &path);
    void Fibre_PlaceDNA(vector<G4ThreeVector> &path);
    void Fibre_ApplyRotation(G4ThreeVector &Rotated,
                       G4ThreeVector &Position,
                       G4RotationMatrix *Rot);
    void Fibre_Bezier(G4ThreeVector &start,
                G4ThreeVector &MidPoint1,
                G4ThreeVector & MidPoint2,
                G4ThreeVector &end,
                vector<G4ThreeVector> &path,
                G4int nSteps);

};

#endif
