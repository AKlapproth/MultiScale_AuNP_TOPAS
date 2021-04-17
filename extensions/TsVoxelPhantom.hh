//
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

#ifndef TsVoxelPhantom_hh
#define TsVoxelPhantom_hh

#include "TsVGeometryComponent.hh"
#include "G4PVParameterised.hh"
#include "TsVoxelParameterisation.hh"


class TsVoxelPhantom : public TsVGeometryComponent
{    
public:
	TsVoxelPhantom(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsVoxelPhantom();
	
	G4VPhysicalVolume* Construct();
    
    void DefineMaterials();

    void ReadTumorParameters();
    
private:

    std::vector <G4int> v_mat;
  
    G4bool CheckOverlap;
    G4bool testmode;
    G4bool skinVis;
    G4bool otisVis;

    G4bool AirExists;
    G4bool OtherTissueExists;
    G4bool SkeletonExists;
    G4bool HeartExists;
    G4bool LungExists;
    G4bool LiverExists;
    G4bool StomachExists;
    G4bool KidneyExists;
    G4bool SmIntestineExists;
    G4bool SpleenExists;
    G4bool BladderExists;
    G4bool TestesExists;
    G4bool SkinExists;
    G4bool GallBladderExists;
    G4bool BrainExists;
    G4bool ThyroidExists;
    G4bool PancreasExists;
    G4bool VasDefExists;
    G4bool LaIntestineExists;
    G4bool AirwayExists;
    G4bool BloodExists;

    G4bool TumorExists;

    G4bool AreaZExists;

    G4int AirAt = -1;
    G4int OtherTissueAt = -1;
    G4int SkeletonAt = -1;
    G4int HeartAt = -1;
    G4int LungAt = -1;
    G4int LiverAt = -1;
    G4int StomachAt = -1;
    G4int KidneyAt = -1;
    G4int SmIntestineAt = -1;
    G4int SpleenAt = -1;
    G4int BladderAt = -1;
    G4int TestesAt = -1;
    G4int SkinAt = -1;
    G4int GallBladderAt = -1;
    G4int BrainAt = -1;
    G4int ThyroidAt = -1;
    G4int PancreasAt = -1;
    G4int VasDefAt = -1;
    G4int LaIntestineAt = -1;
    G4int AirwayAt = -1;
    G4int BloodAt = -1;
    G4int TumorAt = -1;

    G4int TumorCtrX,TumorCtrY,TumorCtrZ;
    G4double TumorRadX,TumorRadY,TumorRadZ;

    G4int* AreaZ;

// MaterialPositions conatins the position of each defined
// material flag in the MaterialFlags array. This regards the
// case that not all flags are assigned in the parameter file.
    G4int MaterialPositions[22];
    G4int MaterialFlags[22];

    G4double ContainerHalfSizeX;
    G4double ContainerHalfSizeY;
    G4double ContainerHalfSizeZ;

    G4double VoxelXLength;
	G4double VoxelYLength;
	G4double VoxelZLength;

    G4double VoxelXHalfLength;
	G4double VoxelYHalfLength;
	G4double VoxelZHalfLength;

    G4int VoxelXCount;
	G4int VoxelYCount;
	G4int VoxelZCount;

    G4int VoxelCount;

    G4int Flag_Air;
    G4int Flag_OtherTissue;
    G4int Flag_Skeleton;
    G4int Flag_Heart;
    G4int Flag_Lung;
    G4int Flag_Liver;
    G4int Flag_Stomach;
    G4int Flag_Kidney;
    G4int Flag_SmIntestine;
    G4int Flag_Spleen;
    G4int Flag_Bladder;
    G4int Flag_Testes;
    G4int Flag_Skin;
    G4int Flag_GallBladder;
    G4int Flag_Brain;
    G4int Flag_Thyroid;
    G4int Flag_Pancreas;
    G4int Flag_VasDef;
    G4int Flag_LaIntestine;
    G4int Flag_Airway;
    G4int Flag_Blood;
    G4int Flag_Tumor;

    G4int count;

    static TsVoxelParameterisation* param;

    std::vector<G4Material*> MaterialVector; 

    G4int fMatNb;

    G4double fColors[22][3];

};

#endif
