// Extra Class for TsFibreDNA
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

#include "VoxelParameterisation.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "TsParameterManager.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"


VoxelParameterisation::VoxelParameterisation()
    : G4VPVParameterisation()
{}

VoxelParameterisation::~VoxelParameterisation()
{}

void VoxelParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector trans = GetTranslation( copyNo );
  physVol->SetTranslation( trans );
}


G4ThreeVector VoxelParameterisation::GetTranslation(const G4int copyNo ) const
{
  CheckCopyNo( copyNo );

  size_t nx;
  size_t ny;
  size_t nz;

  ComputeVoxelIndices( copyNo, nx, ny, nz );

  G4ThreeVector trans( (2*nx+1)*fVoxelHalfX - fContainerWallX,
                       (2*ny+1)*fVoxelHalfY - fContainerWallY,
                       (2*nz+1)*fVoxelHalfZ - fContainerWallZ);

  return trans;
}

void VoxelParameterisation::ComputeVoxelIndices(const G4int copyNo, size_t& nx,
                            size_t& ny, size_t& nz ) const
{
  CheckCopyNo( copyNo );
  nx = size_t(copyNo%fNoVoxelX);
  ny = size_t( (copyNo/fNoVoxelX)%fNoVoxelY );
  nz = size_t(copyNo/fNoVoxelXY);
}

void VoxelParameterisation::SetNoVoxel( size_t nx, size_t ny, size_t nz )
{
  fNoVoxelX = nx; 
  fNoVoxelY = ny; 
  fNoVoxelZ = nz; 
  fNoVoxelXY = nx*ny; 
  fNoVoxel = nx*ny*nz;
}

void VoxelParameterisation::SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz )
{
  fVoxelHalfX = halfx; 
  fVoxelHalfY = halfy; 
  fVoxelHalfZ = halfz; 
}

void VoxelParameterisation::SetContainerDimensions( G4double halfx, G4double halfy, G4double halfz )
{
    fContainerWallX = halfx; 
    fContainerWallY = halfy; 
    fContainerWallZ = halfz; 
}

void VoxelParameterisation::CheckCopyNo( const G4int copyNo ) const
{ 
  if( copyNo < 0 || copyNo >= G4int(fNoVoxel) )
  {
    std::ostringstream message;
    message << "Copy number is negative or too big!" << G4endl
            << "        Copy number: " << copyNo << G4endl
            << "        Total number of voxels: " << fNoVoxel;
    G4Exception("VoxPhantomParameterisation::CheckCopyNo()",
                "GeomNav0002", FatalErrorInArgument, message);
  }
}

