// Extra Class for TsVoxelParameterisation
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

#include "TsVoxelParameterisation.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdint.h>


TsVoxelParameterisation::TsVoxelParameterisation()
    : G4VPVParameterisation()
{}

TsVoxelParameterisation::~TsVoxelParameterisation()
{}

void TsVoxelParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

    G4ThreeVector trans = GetTranslation( copyNo );
    physVol->SetTranslation( trans );
        
}

void TsVoxelParameterisation::ComputeDimensions
(G4Box& trackerLayer, const G4int copyNo, const G4VPhysicalVolume* physVol) const
{
  trackerLayer.SetXHalfLength(fVoxelHalfX);
  trackerLayer.SetYHalfLength(fVoxelHalfY);
  trackerLayer.SetZHalfLength(fVoxelHalfZ);
}

G4VSolid* TsVoxelParameterisation::ComputeSolid
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4VSolid* boxo = new G4Box("Box3000",fVoxelHalfX,fVoxelHalfY,fVoxelHalfZ);
  return boxo;
}

G4Material* TsVoxelParameterisation::ComputeMaterial
(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable *parentTouch=0)
{
  CheckCopyNo(copyNo);
  G4Material* mate;
  G4int flag = fAllVoxelFlags[copyNo];
  G4bool found = false;
  G4int posi;

  G4VisAttributes * VoxelVisAtt = new G4VisAttributes();
  VoxelVisAtt->SetForceWireframe(true);
  VoxelVisAtt->SetForceAuxEdgeVisible(true);

  G4bool inTumor = false;
  G4bool aroundTumor = false;

  // insert tumor
  if (TumorExists){
    G4int Nx,Ny,Nz;
    ComputeVoxelIndices(copyNo,Nx,Ny,Nz);
    G4ThreeVector trans_here = GetTranslation(copyNo)/mm;
    G4ThreeVector trans_tum = GetTranslation(fTumorCenterCopyNo)/mm;
    G4ThreeVector shift_tum = trans_here - trans_tum;
    G4double elltest = (shift_tum[0]/fTumorRadiusX)*(shift_tum[0]/fTumorRadiusX) +
      (shift_tum[1]/fTumorRadiusY)*(shift_tum[1]/fTumorRadiusY) +
      (shift_tum[2]/fTumorRadiusZ)*(shift_tum[2]/fTumorRadiusZ);

    //G4cout << trans_tum << G4endl;

    G4double VoxelDiameter = sqrt((fVoxelHalfX/mm)*(fVoxelHalfX/mm) + 
      (fVoxelHalfY/mm)*(fVoxelHalfY/mm) + (fVoxelHalfZ/mm)*(fVoxelHalfZ/mm));
    if (VoxelDiameter < .3)
      VoxelDiameter = .3; // skin thickness is max(voxelDiameter,0.3mm)
    G4double TumRadExtendedX = fTumorRadiusX + VoxelDiameter;
    G4double TumRadExtendedY = fTumorRadiusY + VoxelDiameter;
    G4double TumRadExtendedZ = fTumorRadiusZ + VoxelDiameter;
    G4double elltest2 = (shift_tum[0]/TumRadExtendedX)*(shift_tum[0]/TumRadExtendedX) +
      (shift_tum[1]/TumRadExtendedY)*(shift_tum[1]/TumRadExtendedY) +
      (shift_tum[2]/TumRadExtendedZ)*(shift_tum[2]/TumRadExtendedZ);

    if (elltest <= 1)
      inTumor = true;
    else if (elltest2 <= 1)
      aroundTumor = true;
  } 


  for (G4int i = 0; i < 22; i ++){
    if (flag == fMaterialFlags[i]){

      if (inTumor){
        if (i == 0 || i == 1 || i == 12) // only other tissue, skin and air can be replaced by tumor tissue
          i = 21;
        else if(!TumorOverlap){
          G4cout << G4endl << " -###- Warning: Tumor overlaps at least with compartment " <<
            i << ". -###- " << G4endl << G4endl;
          TumorOverlap = true;
        }
      }
      else if (aroundTumor && i == 0)
        i = 12; // air around tumor is replaced by skin

      posi = fMaterialPositions[i];
      //G4cout << "posi: " << posi << G4endl;
      mate = fMaterials[posi];
      if (found == true)
        G4cout << "Flag " << flag <<
        " assigned at least twice. Please do not use negative flags." << G4endl;
      else
        found = true;
      //G4cout << "Voxel No " << copyNo << "; Flag " << flag << G4endl;
      if (AreaZSet){
        G4int nx,ny,nz;
        ComputeVoxelIndices(copyNo,nx,ny,nz);
        if (nz < fZ1 || nz > fZ2)
          VoxelVisAtt->SetVisibility(false); // only show z-slices within AreaZ
        else
          VoxelVisAtt->SetVisibility(true);
      }
      if (i == 0 && AirExists)
        VoxelVisAtt->SetVisibility(false); // make air invisible
      else if (i == 12 && !SkinIsVisible)
        VoxelVisAtt->SetVisibility(false); // make skin invisible
      else if (i == 1 && !OthTissIsVisible)
        VoxelVisAtt->SetVisibility(false); // make other tissue invisible
      else if (!AreaZSet)
        VoxelVisAtt->SetVisibility(true);

      const G4Colour& Col =
        G4Colour(fMaterialColors[i][0],fMaterialColors[i][1],fMaterialColors[i][2]);
      VoxelVisAtt->SetColour(Col);
    }
  }
  if (!found) G4cerr << "Flag could not be found; Flag = " << flag
                     << "; copyNo = " << copyNo
                     << "; Voxel count = " << fAllVoxelFlags.size() << G4endl;
  
  physVol->GetLogicalVolume()->SetVisAttributes(VoxelVisAtt);
  return mate;
}

G4ThreeVector TsVoxelParameterisation::GetTranslation(const G4int copyNo) const
{
  CheckCopyNo(copyNo);

  G4int nx,ny,nz;
  ComputeVoxelIndices(copyNo,nx,ny,nz);

  G4ThreeVector trans( (2*nx+1)*fVoxelHalfX - fContainerWallX,
                       (2*ny+1)*fVoxelHalfY - fContainerWallY,
                       (2*nz+1)*fVoxelHalfZ - fContainerWallZ);

  return trans;
}

void TsVoxelParameterisation::ComputeVoxelIndices
(const G4int copyNo, G4int& nx, G4int& ny, G4int& nz) const
{
  CheckCopyNo(copyNo);
  nx = copyNo%fNoVoxelX;
  ny = (copyNo/fNoVoxelX)%fNoVoxelY;
  nz = copyNo/fNoVoxelXY;
}

void TsVoxelParameterisation::SetNoVoxel(G4int fnx, G4int fny, G4int fnz)
{
  fNoVoxelX = fnx; 
  fNoVoxelY = fny; 
  fNoVoxelZ = fnz; 
  fNoVoxelXY = fnx*fny;
  fNoVoxel = fnx*fny*fnz;
}

void TsVoxelParameterisation::SetVoxelDimensions
(G4double halfx, G4double halfy, G4double halfz)
{
  fVoxelHalfX = halfx; 
  fVoxelHalfY = halfy; 
  fVoxelHalfZ = halfz; 
}

void TsVoxelParameterisation::SetContainerDimensions
(G4double halfx, G4double halfy, G4double halfz)
{
  fContainerWallX = halfx; 
  fContainerWallY = halfy; 
  fContainerWallZ = halfz; 
}

void TsVoxelParameterisation::SetMaterialPositions(G4int MatPos[22])
{
  for (G4int i=0;i<22;i++){
    fMaterialPositions[i] = MatPos[i];
  }
}

void TsVoxelParameterisation::SetMaterialFlags(G4int MatFla[22])
{
  for (G4int i=0;i<22;i++){
    fMaterialFlags[i] = MatFla[i];
  }
}

void TsVoxelParameterisation::SetTumorParameters
  (G4int tccn, G4double trX, G4double trY, G4double trZ)
{
  fTumorCenterCopyNo = tccn;

  fTumorRadiusX = trX;
  fTumorRadiusY = trY;
  fTumorRadiusZ = trZ;

  TumorExists = true;
}

void TsVoxelParameterisation::SetMaterialColors(G4double MatCol[22][3])
{
  for (G4int i=0;i<22;i++){
    for (G4int j=0;j<3;j++){
        fMaterialColors[i][j] = MatCol[i][j];
      }
  }
}

void TsVoxelParameterisation::CheckCopyNo( const G4int copyNo ) const
{ 
  if(copyNo < 0 || copyNo >= fNoVoxel)
  {
    std::ostringstream message;
    message << "Copy number is negative or too big!" << G4endl
            << "        Copy number: " << copyNo << G4endl
            << "        Total number of voxels: " << fNoVoxel;
    G4Exception("VoxPhantomParameterisation::CheckCopyNo()",
                "GeomNav0002", FatalErrorInArgument, message);
  }
}

void TsVoxelParameterisation::SetAreaZ(G4int z1, G4int z2){
    fZ1 = z1;
    fZ2 = z2;
    AreaZSet = true;
}

