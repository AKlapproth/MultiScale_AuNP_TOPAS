// Author: Hongyu Zhu
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

#ifndef VOXELPARAMETERISATION_HH
#define VOXELPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "TsVGeometryComponent.hh"

#include "VoxelParameterisation.hh"

class VoxelParameterisation : public G4VPVParameterisation
{
public:
    VoxelParameterisation();
    virtual ~VoxelParameterisation();
    virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

    
    void ComputeVoxelIndices(const G4int copyNo, size_t& nx,size_t& ny, size_t& nz ) const;// Convert the copyNo to voxel numbers in x, y and z.
    void SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetContainerDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetNoVoxel( size_t nx, size_t ny, size_t nz );
    void CheckCopyNo( const G4int copyNo ) const;
    G4ThreeVector GetTranslation(const G4int copyNo ) const;


protected:

    G4double fVoxelHalfX,fVoxelHalfY,fVoxelHalfZ;
    G4double fContainerWallX, fContainerWallY, fContainerWallZ;
    size_t fNoVoxelX,fNoVoxelY,fNoVoxelZ;
    size_t fNoVoxelXY;
    size_t fNoVoxel;

};

#endif // VOXELPARAMETERISATION_HH
