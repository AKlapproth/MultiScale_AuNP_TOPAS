
// Author: Hongyu Zhu & Alexander Klapproth
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

#ifndef TsVoxelParameterisation_HH
#define TsVoxelParameterisation_HH

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"

#include "TsVoxelParameterisation.hh"

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdint.h>

class TsVoxelParameterisation : public G4VPVParameterisation
{
public:
    TsVoxelParameterisation();
    virtual ~TsVoxelParameterisation();

    virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;
    virtual void ComputeDimensions(G4Box& trackerLayer, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
    virtual G4VSolid* ComputeSolid(const G4int copyNo, G4VPhysicalVolume* physVol) const;
    virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable* parentTouch);


    G4ThreeVector GetTranslation(const G4int copyNo) const;
    // G4String GetCompartmentName(G4int copyNo);
    
    void ComputeVoxelIndices(const G4int copyNo, G4int& nx, G4int& ny, G4int& nz) const;// Convert the copyNo to voxel numbers in x, y and z.
    
    void SetVoxelDimensions(G4double halfx, G4double halfy, G4double halfz);
    void SetContainerDimensions(G4double halfx, G4double halfy, G4double halfz);
    void SetNoVoxel(G4int fnx, G4int fny, G4int fnz);
    void CheckCopyNo(const G4int copyNo) const;
    void SetMaterials(std::vector<G4Material*> Mat){
        fMaterials = Mat; }
    void SetMatrixWithFlags(std::vector <G4int> v_mat){
        fAllVoxelFlags = v_mat; }
    void SetTumorParameters(G4int tccn, G4double trX, G4double trY, G4double trZ);
    void SetAirExists(G4bool ae){
      AirExists = ae; }
    void SetMaterialPositions(G4int MatPos[22]);
    void SetMaterialFlags(G4int MatFla[22]);
    void SetMaterialColors(G4double MatCol[22][3]);
    void SetSkinIsVisible(G4bool sv){
      SkinIsVisible = sv; }
    void SetTissIsVisible(G4bool tv){
      OthTissIsVisible = tv; }

    void SetAreaZ(G4int z1, G4int z2);



protected:
    G4double fVoxelHalfX,fVoxelHalfY,fVoxelHalfZ;
      // Half dimension of voxels (assume they are boxes).
    G4double fContainerWallX, fContainerWallY, fContainerWallZ;
      // Save position of container wall for speed-up.
    G4int fNoVoxelX,fNoVoxelY,fNoVoxelZ;
      // Number of voxel in x, y and z dimensions.
    G4int fNoVoxelXY;
      // Number of voxels in x times number of voxels in y (for speed-up).
    G4int fNoVoxel;
      // Total number of voxels (for speed-up).
    std::vector <G4Material*> fMaterials;
      // Vector containing all used materials
    std::vector <G4int> fAllVoxelFlags;
      // Vector containing the flags for all voxels
    G4bool AirExists;
    G4int fMaterialFlags[22];
      // Flags for each material
    G4int fMaterialPositions[22];
      // Position of each material within fMaterials
    G4double fMaterialColors[22][3];
      // Vector containing 13 colors

    G4int fTumorCenterCopyNo;
    G4double fTumorRadiusX, fTumorRadiusY, fTumorRadiusZ;

    G4bool SkinIsVisible = "true";
    G4bool OthTissIsVisible = "true";

    G4String fCompartmentName;

    G4int fZ1, fZ2;
    G4bool TumorExists = false;
    G4bool AreaZSet = false;

    G4bool TumorOverlap = false;
};

#endif // TsVoxelParameterisation_HH
