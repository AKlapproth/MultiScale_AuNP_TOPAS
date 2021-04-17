// Component for TsVoxelPhantom
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
//******************************************************************


#include "TsVoxelPhantom.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"
#include "Randomize.hh"


#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdint.h>

using namespace std;

TsVoxelParameterisation* TsVoxelPhantom::param;


TsVoxelPhantom::TsVoxelPhantom(TsParameterManager* pM, TsExtensionManager* eM,
	TsMaterialManager* mM, TsGeometryManager* gM,
	TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{}


TsVoxelPhantom::~TsVoxelPhantom()
{}


G4VPhysicalVolume* TsVoxelPhantom::Construct()
{
	BeginConstruction();

  // User defined parameters.
  // Geometry parameters
  VoxelXLength = .71*mm;
  if (fPm->ParameterExists(GetFullParmName("VoxelXLength"))){
    VoxelXLength = fPm->GetDoubleParameter(GetFullParmName("VoxelXLength"),"Length");
  }

  VoxelYLength = .71*mm;
  if (fPm->ParameterExists(GetFullParmName("VoxelYLength"))){
    VoxelYLength = fPm->GetDoubleParameter(GetFullParmName("VoxelYLength"),"Length");
  }

  VoxelZLength = .71*mm;
  if (fPm->ParameterExists(GetFullParmName("VoxelZLength"))){
    VoxelZLength = fPm->GetDoubleParameter(GetFullParmName("VoxelZLength"),"Length");
  }

  VoxelXCount = 1;
  if (fPm->ParameterExists(GetFullParmName("VoxelXCount"))){
    VoxelXCount = fPm->GetIntegerParameter(GetFullParmName("VoxelXCount"));
  }

  VoxelYCount = 1;
  if (fPm->ParameterExists(GetFullParmName("VoxelYCount"))){
    VoxelYCount = fPm->GetIntegerParameter(GetFullParmName("VoxelYCount"));
  }

  VoxelZCount = 1;
  if (fPm->ParameterExists(GetFullParmName("VoxelZCount"))){
    VoxelZCount = fPm->GetIntegerParameter(GetFullParmName("VoxelZCount"));
  }

  // Flags for Materials
  AirExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Air"));
  OtherTissueExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_OtherTissue"));
  SkeletonExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Skeleton"));
  HeartExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Heart"));
  LungExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Lung"));
  LiverExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Liver"));
  StomachExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Stomach"));
  KidneyExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Kidney"));
  SmIntestineExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Small_Intestine"));
  SpleenExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Spleen"));
  BladderExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Bladder"));
  TestesExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Testes"));
  SkinExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Skin"));
  GallBladderExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Gall_Bladder"));
  BrainExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Brain"));
  ThyroidExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Thyroid"));
  PancreasExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Pancreas"));
  VasDefExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_VasDef"));
  LaIntestineExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Large_Intestine"));
  AirwayExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Airway"));
  BloodExists = fPm->ParameterExists(GetFullParmName("FlagForMaterial_Blood"));

  TumorExists = false;
  if (fPm->ParameterExists(GetFullParmName("Tumor")))
    TumorExists = fPm->GetBooleanParameter(GetFullParmName("Tumor"));

  AreaZExists = fPm->ParameterExists(GetFullParmName("ShowAreaZ"));

  Flag_Air = -1;
  if (AirExists){
    Flag_Air = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Air"));
  }

  Flag_OtherTissue = -2;
  if (OtherTissueExists){
    Flag_OtherTissue = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_OtherTissue"));
  }

  Flag_Skeleton = -3;
  if (SkeletonExists){
    Flag_Skeleton = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Skeleton"));
  }

  Flag_Heart = -4;
  if (HeartExists){
    Flag_Heart = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Heart"));
  }

  Flag_Lung = -5;
  if (LungExists){
    Flag_Lung = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Lung"));
  }

  Flag_Liver = -6;
  if (LiverExists){
    Flag_Liver = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Liver"));
  }

  Flag_Stomach = -7;
  if (StomachExists){
    Flag_Stomach = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Stomach"));
  }

  Flag_Kidney = -8;
  if (KidneyExists){
    Flag_Kidney = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Kidney"));
  }

  Flag_SmIntestine = -9;
  if (SmIntestineExists){
    Flag_SmIntestine = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Small_Intestine"));
  }

  Flag_Spleen = -10;
  if (SpleenExists){
    Flag_Spleen = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Spleen"));
  }

  Flag_Bladder = -11;
  if (BladderExists){
    Flag_Bladder = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Bladder"));
  }

  Flag_Testes = -12;
  if (TestesExists){
    Flag_Testes = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Testes"));
  }

  Flag_Skin = -13;
  if (SkinExists){
    Flag_Skin = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Skin"));
  }

  Flag_GallBladder = -14;
  if (GallBladderExists){
    Flag_GallBladder = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Gall_Bladder"));
  }

  Flag_Brain = -15;
  if (BrainExists){
    Flag_Brain = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Brain"));
  }

  Flag_Thyroid = -16;
  if (ThyroidExists){
    Flag_Thyroid = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Thyroid"));
  }

  Flag_Pancreas = -17;
  if (PancreasExists){
    Flag_Pancreas = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Pancreas"));
  }

  Flag_VasDef = -18;
  if (VasDefExists){
    Flag_VasDef = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_VasDef"));
  }

  Flag_LaIntestine = -19;
  if (LaIntestineExists){
    Flag_LaIntestine = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Large_Intestine"));
  }

  Flag_Airway = -20;
  if (AirwayExists){
    Flag_Airway = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Airway"));
  }

  Flag_Blood = -21;
  if (BloodExists){
    Flag_Blood = fPm->GetIntegerParameter(GetFullParmName("FlagForMaterial_Blood"));
  }

  Flag_Tumor = -100;

  G4int MaterialFlags[22] = {Flag_Air,Flag_OtherTissue,Flag_Skeleton,Flag_Heart,Flag_Lung,
    Flag_Liver,Flag_Stomach,Flag_Kidney,Flag_SmIntestine,Flag_Spleen,Flag_Bladder,Flag_Testes,
    Flag_Skin,Flag_GallBladder,Flag_Brain,Flag_Thyroid,Flag_Pancreas,Flag_VasDef,Flag_LaIntestine,
    Flag_Airway,Flag_Blood,Flag_Tumor};


  // Read Binary File
  G4String fname = GetFullParmName("FileName");
  if (!fPm->ParameterExists(fname)) {
      G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
      G4cerr << "Parameter " << fname << " has to be specified for this geometry." << G4endl;
      exit(1);
  }
  G4String FileName = fPm->GetStringParameter(fname);

  CheckOverlap = false;
  if (fPm->ParameterExists(GetFullParmName("CheckOverlap"))){
    CheckOverlap = fPm->GetBooleanParameter(GetFullParmName("CheckOverlap"));
  } 

  testmode = false;
  if (fPm->ParameterExists(GetFullParmName("testmode"))){
    testmode = fPm->GetBooleanParameter(GetFullParmName("testmode"));
  }

  skinVis = true;
  if (fPm->ParameterExists(GetFullParmName("SkinIsVisible"))){
    skinVis = fPm->GetBooleanParameter(GetFullParmName("SkinIsVisible"));
  }

  otisVis = true;
  if (fPm->ParameterExists(GetFullParmName("OtherTissueIsVisible"))){
    otisVis = fPm->GetBooleanParameter(GetFullParmName("OtherTissueIsVisible"));
  }

  if (TumorExists)
    ReadTumorParameters();

  G4int TumorCenterCopyNo = VoxelXCount*VoxelYCount*TumorCtrZ + VoxelXCount*TumorCtrY + TumorCtrX;

  if (AreaZExists){
    AreaZ = fPm->GetIntegerVector(GetFullParmName("ShowAreaZ"));
  }

  if (testmode){
    VoxelXCount = 1;
    VoxelYCount = 1;
    VoxelZCount = 1;
  }

  VoxelXHalfLength = VoxelXLength/2;
  VoxelYHalfLength = VoxelYLength/2;
  VoxelZHalfLength = VoxelZLength/2;

  VoxelCount = VoxelXCount*VoxelYCount*VoxelZCount;

  //****************************************************************************
  //                        Read  Voxel  data  file                           //
  //****************************************************************************
  const char* filename = FileName;
  v_mat.clear();
  if (testmode){
    // v_mat.assign(4,Flag_Air);
    // v_mat.push_back(Flag_Heart);
    // v_mat.push_back(Flag_Liver);
    // v_mat.push_back(Flag_Liver);
    v_mat.push_back(Flag_Heart);
  }
  G4int here;
  count = 0;

  if (!testmode){
    G4cout << G4endl << " # Reading " << VoxelCount <<
    	" voxels from " << filename << " # " << G4endl << G4endl;
    
    ifstream f(filename, ios::in | ios::binary);

    if (f.is_open()) {
    	while (count < VoxelCount){
	    	f.read((char *) &here, sizeof here);
	        v_mat.push_back(here);
	        count++;
	    }
    	if (count != VoxelCount){
    		G4cerr << "ERROR: Amount of binary file elements is smaller than voxel count." << G4endl <<
    		"File elements : " << count << "; Voxel Count: " << VoxelCount << G4endl;
    		exit(1);
    	}

    }
    else {
        G4cerr << "ERROR: Unable to open file " << FileName << G4endl;
        exit(1);
    }

    f.close();

    G4cout << G4endl << " # Binary file was read successfully. # " << G4endl << G4endl;
  }

  //****************************************************************************
  //                               Set geometry size                          //
  //****************************************************************************
  
  ContainerHalfSizeX = VoxelXHalfLength*VoxelXCount;
	ContainerHalfSizeY = VoxelYHalfLength*VoxelYCount;
	ContainerHalfSizeZ = VoxelZHalfLength*VoxelZCount;

  G4cout << "*********************************************************************************"<<G4endl;
  G4cout << "Voxel X length: " << VoxelXLength/mm << " mm " << G4endl;
  G4cout << "Voxel Y length: " << VoxelYLength/mm << " mm " << G4endl;
  G4cout << "Voxel Z length: " << VoxelZLength/mm << " mm " << G4endl;
  G4cout << "Voxel container X length: " << ContainerHalfSizeX*2/mm << " mm " << G4endl;
  G4cout << "Voxel container Y length: " << ContainerHalfSizeY*2/mm << " mm " << G4endl;
  G4cout << "Voxel container Z length: " << ContainerHalfSizeZ*2/mm << " mm " << G4endl;
  //G4cout<<"*********************************************************************************"<<G4endl;

  //****************************************************************************
  //                                Parameterise                              //
  //****************************************************************************

  DefineMaterials();
  G4cout << fMatNb << " Materials defined." << G4endl;

  //DefineColors();

  G4double fColors[22][3] = {
  {1,1,1},    // Air
  {.8,.8,.8}, // Other Tissue
  {.9,.9,.7}, // Skeleton
  {.8,0,0},    // Heart
  {.7,0,0},   // Lung
  {.5,.3,.3}, // Liver
  {0,.7,0},   // Stomach
  {0,0,.7},   // Kidney
  {0,.8,.8},  // Small Intestine
  {1,1,0},    // Spleen
  {1,0,1},    // Bladder
  {.9,.3,0},  // Testes
  {.3,.3,.3}, // Skin
  {0,1,0},    // Gall Bladder
  {1,.5,1},   // Brain
  {0,.2,1},   // Thyroid
  {.4,.4,0},  // Pancreas
  {.7,0,.7},    // Vas Deferens
  {0,.5,.5},    // Large Intestine
  {1,1,1},    // Airway
  {1,0,0},    // Blood
  {1,.4,0}  // Tumor
  };

////----- Create parameterisation and set
  param = new TsVoxelParameterisation();
  param->SetVoxelDimensions(VoxelXHalfLength, VoxelYHalfLength, VoxelZHalfLength); 
  param->SetNoVoxel(VoxelXCount, VoxelYCount, VoxelZCount); 
  param->SetContainerDimensions(ContainerHalfSizeX, ContainerHalfSizeY, ContainerHalfSizeZ);
  param->SetMaterialPositions(MaterialPositions);
  param->SetAirExists(AirExists);
  param->SetMaterialFlags(MaterialFlags);
  param->SetMaterials(MaterialVector);
  param->SetMatrixWithFlags(v_mat);
  param->SetMaterialColors(fColors);
  if (TumorExists)
    param->SetTumorParameters(TumorCenterCopyNo, TumorRadX, TumorRadY, TumorRadZ);
  param->SetSkinIsVisible(skinVis);
  param->SetTissIsVisible(otisVis);
  if (AreaZExists)
    param->SetAreaZ(AreaZ[0],AreaZ[1]);

  if (!skinVis)
    G4cout << G4endl << " # Skin is not going to be displayed. #" << G4endl;

  if (!otisVis)
    G4cout << " # Undefined Tissue is not going to be displayed. #" << G4endl << G4endl;

  //****************************************************************************
  //                             Build basic geometry                         //
  //****************************************************************************

  // Container
  G4Box* box_solid = new G4Box("ContainerBox", ContainerHalfSizeX, ContainerHalfSizeY, ContainerHalfSizeZ);
  fEnvelopeLog = CreateLogicalVolume(box_solid);
  fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
  G4Colour blue (0.0, 0.0, 1.0);
  G4VisAttributes* boxVis = new G4VisAttributes(blue);
  boxVis->SetVisibility(false);
  fEnvelopeLog->SetVisAttributes(boxVis);

  G4Box           * voxel_solid = new G4Box( "Voxel", VoxelXHalfLength, VoxelYHalfLength, VoxelZHalfLength);
  G4LogicalVolume * voxel_logic = new G4LogicalVolume(voxel_solid,MaterialVector[0],"VoxelLogical",0,0,0);


  G4VPhysicalVolume*  TheVoxelPhantom = CreatePhysicalVolume("VoxelPhantom",voxel_logic, fEnvelopePhys, kXAxis, VoxelCount, param);

  InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}



void  TsVoxelPhantom::DefineMaterials()
{
  MaterialVector.clear();

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* mat_air   = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* mat_tissue   = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
  G4Material* mat_skeleton   = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4Material* mat_heart   = nist->BuildMaterialWithNewDensity("Heart_Material", "G4_WATER", 1.06 *g/cm/cm/cm);
  G4Material* mat_lung   = nist->FindOrBuildMaterial("G4_LUNG_ICRP");
  G4Material* mat_liver   = nist->BuildMaterialWithNewDensity("Liver_Material", "G4_WATER", 1.05 *g/cm/cm/cm);
  G4Material* mat_stomach   = nist->BuildMaterialWithNewDensity("Stomach_Material", "G4_WATER", 1.04 *g/cm/cm/cm);
  G4Material* mat_kidney   = nist->BuildMaterialWithNewDensity("Kidney_Material", "G4_WATER", 1.05 *g/cm/cm/cm);
  G4Material* mat_Sintestine   = nist->BuildMaterialWithNewDensity("Sm_Intest_Material", "G4_WATER", 1.04 *g/cm/cm/cm);
  G4Material* mat_spleen   = nist->BuildMaterialWithNewDensity("Spleen_Material", "G4_WATER", 1.06 *g/cm/cm/cm);
  G4Material* mat_bladder   = nist->BuildMaterialWithNewDensity("Bladder_Material", "G4_WATER", 1.04 *g/cm/cm/cm);
  G4Material* mat_testes   = nist->FindOrBuildMaterial("G4_TESTIS_ICRP");
  G4Material* mat_skin   = nist->FindOrBuildMaterial("G4_SKIN_ICRP");
  G4Material* mat_gallbladder   = nist->BuildMaterialWithNewDensity("Gall_Bladder_Material", "G4_WATER", 1.03 *g/cm/cm/cm);
  G4Material* mat_brain   = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");
  G4Material* mat_thyroid   = nist->BuildMaterialWithNewDensity("Thyroid_Material", "G4_WATER", 1.05 *g/cm/cm/cm);
  G4Material* mat_pancreas   = nist->BuildMaterialWithNewDensity("Pancreas_Material", "G4_WATER", 1.05 *g/cm/cm/cm);
  G4Material* mat_vasdef   = nist->BuildMaterialWithNewDensity("VasDef_Material", "G4_WATER", 1.04 *g/cm/cm/cm);
  G4Material* mat_Lintestine   = nist->BuildMaterialWithNewDensity("La_Intest_Material", "G4_WATER", 1.04 *g/cm/cm/cm);
  G4Material* mat_airway   = nist->BuildMaterialWithNewDensity("Airway_Material", "G4_AIR", 0.00120479 *g/cm/cm/cm);
  G4Material* mat_blood   = nist->BuildMaterialWithNewDensity("Blood_Material", "G4_WATER", 1 *g/cm/cm/cm);
  G4Material* mat_tumor   = nist->FindOrBuildMaterial("G4_WATER");

  G4Material* water   = nist->FindOrBuildMaterial("G4_WATER");

  G4int Mcounter = 0;

  G4int AirAt,OtherTissueAt,SkeletonAt,HeartAt,LungAt,LiverAt,StomachAt,KidneyAt,
    SmIntestineAt,SpleenAt,BladderAt,TestesAt,SkinAt,GallBladderAt,BrainAt,
    ThyroidAt,PancreasAt,VasDefAt,LaIntestineAt,AirwayAt,BloodAt,TumorAt;

  if (AirExists){
    MaterialVector.push_back(mat_air);
    AirAt = Mcounter;
    Mcounter++;
  }
  if (OtherTissueExists){
    MaterialVector.push_back(mat_tissue);
    OtherTissueAt = Mcounter;
    Mcounter++;
  }
  if (SkeletonExists){
    MaterialVector.push_back(mat_skeleton);
    SkeletonAt = Mcounter;
    Mcounter++;
  }
  if (HeartExists){
    MaterialVector.push_back(mat_heart);
    HeartAt = Mcounter;
    Mcounter++;
  }
  if (LungExists){
    MaterialVector.push_back(mat_lung);
    LungAt = Mcounter;
    Mcounter++;
  }
  if (LiverExists){
    MaterialVector.push_back(mat_liver);
    LiverAt = Mcounter;
    Mcounter++;
  }
  if (StomachExists){
    MaterialVector.push_back(mat_stomach);
    StomachAt = Mcounter;
    Mcounter++;
  }
  if (KidneyExists){
    MaterialVector.push_back(mat_kidney);
    KidneyAt = Mcounter;
    Mcounter++;
  }
  if (SmIntestineExists){
    MaterialVector.push_back(mat_Sintestine);
    SmIntestineAt = Mcounter;
    Mcounter++;
  }
  if (SpleenExists){
    MaterialVector.push_back(mat_spleen);
    SpleenAt = Mcounter;
    Mcounter++;
  }
  if (BladderExists){
    MaterialVector.push_back(mat_bladder);
    BladderAt = Mcounter;
    Mcounter++;
  }
  if (TestesExists){
    MaterialVector.push_back(mat_testes);
    TestesAt = Mcounter;
    Mcounter++;
  }
  if (SkinExists){
    MaterialVector.push_back(mat_skin);
    SkinAt = Mcounter;
    Mcounter++;
  }
  if (GallBladderExists){
    MaterialVector.push_back(mat_gallbladder);
    GallBladderAt = Mcounter;
    Mcounter++;
  }
  if (BrainExists){
    MaterialVector.push_back(mat_brain);
    BrainAt = Mcounter;
    Mcounter++;
  }
  if (ThyroidExists){
    MaterialVector.push_back(mat_thyroid);
    ThyroidAt = Mcounter;
    Mcounter++;
  }
  if (PancreasExists){
    MaterialVector.push_back(mat_pancreas);
    PancreasAt = Mcounter;
    Mcounter++;
  }
  if (VasDefExists){
    MaterialVector.push_back(mat_vasdef);
    VasDefAt = Mcounter;
    Mcounter++;
  }
  if (LaIntestineExists){
    MaterialVector.push_back(mat_Lintestine);
    LaIntestineAt = Mcounter;
    Mcounter++;
  }
  if (AirwayExists){
    MaterialVector.push_back(mat_airway);
    AirwayAt = Mcounter;
    Mcounter++;
  }
  if (BloodExists){
    MaterialVector.push_back(mat_blood);
    BloodAt = Mcounter;
    Mcounter++;
  }
  if (TumorExists){
    MaterialVector.push_back(mat_tumor);
    TumorAt = Mcounter;
    Mcounter++;
  }

  fMatNb = MaterialVector.size();
  MaterialPositions[0] = AirAt;
  MaterialPositions[1] = OtherTissueAt;
  MaterialPositions[2] = SkeletonAt;
  MaterialPositions[3] = HeartAt;
  MaterialPositions[4] = LungAt;
  MaterialPositions[5] = LiverAt;
  MaterialPositions[6] = StomachAt;
  MaterialPositions[7] = KidneyAt;
  MaterialPositions[8] = SmIntestineAt;
  MaterialPositions[9] = SpleenAt;
  MaterialPositions[10] = BladderAt;
  MaterialPositions[11] = TestesAt;
  MaterialPositions[12] = SkinAt;
  MaterialPositions[13] = GallBladderAt;
  MaterialPositions[14] = BrainAt;
  MaterialPositions[15] = ThyroidAt;
  MaterialPositions[16] = PancreasAt;
  MaterialPositions[17] = VasDefAt;
  MaterialPositions[18] = LaIntestineAt;
  MaterialPositions[19] = AirwayAt;
  MaterialPositions[20] = BloodAt;
  MaterialPositions[21] = TumorAt;
}

void TsVoxelPhantom::ReadTumorParameters()
{
  TumorCtrX = 0;
  TumorCtrY = 0;
  TumorCtrZ = 0;

  TumorRadX = 0*mm;
  TumorRadY = 0*mm;
  TumorRadZ = 0*mm;

  if (fPm->ParameterExists(GetFullParmName("TumorCenterX"))){
    TumorCtrX = fPm->GetIntegerParameter(GetFullParmName("TumorCenterX"));
  }
  if (fPm->ParameterExists(GetFullParmName("TumorCenterY"))){
    TumorCtrY = fPm->GetIntegerParameter(GetFullParmName("TumorCenterY"));
  }
  if (fPm->ParameterExists(GetFullParmName("TumorCenterZ"))){
    TumorCtrZ = fPm->GetIntegerParameter(GetFullParmName("TumorCenterZ"));
  }

  if (fPm->ParameterExists(GetFullParmName("TumorRadiusX"))){
    TumorRadX = fPm->GetDoubleParameter(GetFullParmName("TumorRadiusX"),"Length");
  }
  if (fPm->ParameterExists(GetFullParmName("TumorRadiusY"))){
    TumorRadY = fPm->GetDoubleParameter(GetFullParmName("TumorRadiusY"),"Length");
  }
  if (fPm->ParameterExists(GetFullParmName("TumorRadiusZ"))){
    TumorRadZ = fPm->GetDoubleParameter(GetFullParmName("TumorRadiusZ"),"Length");
  }
}

