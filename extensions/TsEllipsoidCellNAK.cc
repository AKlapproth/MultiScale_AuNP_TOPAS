// Component for TsEllipsoidCellNAK
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

#include "TsEllipsoidCellNAK.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "G4TwoVector.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
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

void TsEllipsoidCellNAK::int2str(const int &int_temp,string &string_temp) 
{  
        stringstream stream;  
        stream<<int_temp;  
        string_temp=stream.str();
}  



TsEllipsoidCellNAK::TsEllipsoidCellNAK(
    TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
	TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsEllipsoidCellNAK::~TsEllipsoidCellNAK()
{
}


G4VPhysicalVolume* TsEllipsoidCellNAK::Construct()
{
	BeginConstruction();
    
    //User defined parameters:
    //User can specify the cell semi-axis lengths or use the default values.
    
    //Set default values for cell dimensions
    G4double xSA = 20 * um;
    G4double ySA = 10 * um;
    G4double zSA = 15 * um;
    
    G4String name = GetFullParmName("xSemiAxis");
    if (fPm->ParameterExists(name)) {
        
        //User specified cell size values
        xSA = fPm->GetDoubleParameter(GetFullParmName("xSemiAxis"), "Length");
        ySA = fPm->GetDoubleParameter(GetFullParmName("ySemiAxis"), "Length");
        zSA = fPm->GetDoubleParameter(GetFullParmName("zSemiAxis"), "Length");
        
    }

    // Nucleus parameters

    FibreRadius = (37.0879/2.0) *nm;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/FibreRadius"))){
      FibreRadius = fPm->GetDoubleParameter(GetFullParmName("Nucleus/FibreRadius"),"Length");
    }

    FibreLength = 120 *nm; 
    if (fPm->ParameterExists(GetFullParmName("Nucleus/FibreLength"))){
      FibreLength = fPm->GetDoubleParameter(GetFullParmName("Nucleus/FibreLength"),"Length");
    }

    NuclRadius = 0*um;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/NucleusRadius"))){
      NuclRadius = fPm->GetDoubleParameter(GetFullParmName("Nucleus/NucleusRadius"),"Length");
    }  

    HilbertCurveLayer = 1;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/HilbertCurveLayer"))){
      HilbertCurveLayer = fPm->GetIntegerParameter(GetFullParmName("Nucleus/HilbertCurveLayer"));
    } 

   HilbertCurve3DRepeat = 1;  // repeat hilbert curve n times on X/Y/Z directions
    if (fPm->ParameterExists(GetFullParmName("Nucleus/HilbertCurve3DRepeat"))){
      HilbertCurve3DRepeat = fPm->GetIntegerParameter(GetFullParmName("Nucleus/HilbertCurve3DRepeat"));
    } 

    G4String FileName;
    fNucleusExists = true;
    if (!fPm->ParameterExists(GetFullParmName("Nucleus/FileName"))) {
        fNucleusExists = false;
    }
    else
        FileName = fPm->GetStringParameter(GetFullParmName("Nucleus/FileName"));

    CheckOverlap = false;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/CheckOverlap"))){
      CheckOverlap = fPm->GetBooleanParameter(GetFullParmName("Nucleus/CheckOverlap"));
    } 

    testmode = false;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/testmode"))){
      testmode = fPm->GetBooleanParameter(GetFullParmName("Nucleus/testmode"));
    } 

    AddHydrationShell = true;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/AddHydrationShell"))){
      AddHydrationShell = fPm->GetBooleanParameter(GetFullParmName("Nucleus/AddHydrationShell"));
    } 

    HydrationShellThickness = 0.16*nm;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/HydrationShellThickness"))){
      HydrationShellThickness = fPm->GetDoubleParameter(GetFullParmName("Nucleus/HydrationShellThickness"),"Length");
    }

    ToleranceSpace = .001*um;
    if (fPm->ParameterExists(GetFullParmName("Nucleus/ToleranceSpace"))){
      ToleranceSpace = fPm->GetDoubleParameter(GetFullParmName("Nucleus/ToleranceSpace"),"Length");
    }

    transNucX = 0*um;
    transNucY = 0*um;
    transNucZ = 0*um;
    name = GetFullParmName("Nucleus/NuclTransX");
    if (fPm->ParameterExists(name)) {
        
        G4String name1 = GetFullParmName("Nucleus/NuclTransX");
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
            if (transNucX > xSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/NuclTransY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
            if (transNucY > ySA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/NuclTransZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
            if (transNucZ > zSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
    }
    
    G4ThreeVector* NucleusTranslation = new G4ThreeVector(transNucX, transNucY, transNucZ);
    G4RotationMatrix* NucleusRotation = new G4RotationMatrix(0,0,0);
	
    //***********************************************************************
    //              Envelope Geometry : ellipsoid cell
    //***********************************************************************
    
	G4Ellipsoid* envelopeSolid = new G4Ellipsoid(fName, xSA, ySA, zSA);
	G4LogicalVolume* envelopeLog = CreateLogicalVolume(envelopeSolid);
	G4VPhysicalVolume* fEnvelopePhys = CreatePhysicalVolume(envelopeLog);


        
    //Find biggest semi-axis length of the cell
    G4double CellRadius;
    if ((xSA >= ySA) && (xSA >= zSA)){CellRadius = xSA;}
    else if ((ySA >= xSA) && (ySA >= zSA)){CellRadius = ySA;}
    else {CellRadius = zSA;}


    //*******************************
    // Subcomponent: (hybrid) AuNPs
    //*******************************

    G4int NbOfAuNPsInCytoplasm = 0;
    name = GetFullParmName("AuNPs/NbOfAuNPsInCytoplasm");
    if (fPm->ParameterExists(name)) {

        // Nanoparticle placement requires the nucleus center to be at (0,0,0)!

        if ((pow(transNucX,2) + pow(transNucY,2) + pow(transNucZ,2)) > 0){
            G4cout << "Nanoparticle placement requires the nucleus center to be at (0,0,0)!" << G4endl;
            exit(1);
        }
        
        //Number of AuNPs in the cytoplasm
        NbOfAuNPsInCytoplasm  =
            fPm->GetIntegerParameter(GetFullParmName("AuNPs/NbOfAuNPsInCytoplasm"));

        G4bool MaxDistExists = false;
        G4double MaxDistFromNucleus;
        if (fPm->ParameterExists(GetFullParmName("AuNPs/MaxDistFromNucleus"))){
            MaxDistExists = true;
            MaxDistFromNucleus = fPm->GetDoubleParameter(GetFullParmName("AuNPs/MaxDistFromNucleus"),"Length");
        }
        
        //Radius of the AuNPs (default values if none are specified)
        G4double AuNPRadius = 2*nanometer;
        
        name=GetFullParmName("AuNPs/Rmax");
        if (fPm->ParameterExists(name))
            {AuNPRadius = fPm->GetDoubleParameter(GetFullParmName("AuNPs/Rmax"), "Length" );}
        
        
        G4String subComponentName3 = "AuNPC";
        G4Orb* gAuNPC = new G4Orb("gAuNPC", AuNPRadius);
        G4LogicalVolume* lAuNPC = CreateLogicalVolume(subComponentName3, gAuNPC);

        G4bool ThereIsACore = false;
        G4String subsubComponentName;
        G4LogicalVolume* lCore;
        name=GetFullParmName("AuNPs/Rcore");
        if (fPm->ParameterExists(name)) {
            ThereIsACore = true;

            G4String symbol;
            G4int n_atoms,n_comps;
            G4double z,w,d_Fe3O4;
            G4Element* iron =
            new G4Element("Iron", symbol="Fe", z=26, w=55.845*g/mole);
            G4Element* oxygen =
            new G4Element("Oxygen", symbol="O", z=8, w=15.9994*g/mole);
            G4Material * Fe3O4 =
            new G4Material("Ironoxide", d_Fe3O4=5.17*g/cm3,n_comps=2);

            Fe3O4->AddElement(iron, n_atoms=3);
            Fe3O4->AddElement(oxygen, n_atoms=4);

            subsubComponentName = "AuNPCore";
            G4double Rcore = fPm->GetDoubleParameter(GetFullParmName("AuNPs/Rcore"), "Length" );
            G4Orb* gCore = new G4Orb("gCore",Rcore);
            lCore = new G4LogicalVolume(gCore,Fe3O4,"CoreLogical",0,0,0);
        }
    
        //Randomly distribute AuNPs throughout cytoplasm
        
        G4double xC,yC,zC;
        G4double dist_sq_C;
        G4bool intersectC;
        G4int IterationCounter;

        static const int statC = floor(NbOfAuNPsInCytoplasm/100);
        // static const int statC = NbOfAuNPsInCytoplasm;
        auto NPpositions0 = new G4double[statC][3];
        auto NPpositionsA = new G4double[statC][3];
        auto NPpositionsB = new G4double[statC][3];

        G4cout << " +++ Placing Nanoparticles ... +++ " << G4endl << G4endl;

        // to improve performance, NPs are divided into 50 sectors
        // place NPs in 2 of 50 sectors per iteration
        for (G4int sector = 0; sector < 50; sector++){

            for (G4int jA = 0; jA < statC; jA++){

                intersectC = true;
                IterationCounter = 0;

                while (intersectC){
          
                    intersectC = false;
                    IterationCounter++;

                    if (IterationCounter > 10000000){
                        G4cout << "ERROR: Apparently too little space for AuNPs" << FileName << G4endl;
                        exit(1);
                    }

                    G4double uC = 2*pi*2*sector/100 + G4UniformRand()*2*pi/100;
                    G4double vC = std::acos(2*G4UniformRand()-1);
                    G4double drC;
                    if (MaxDistExists)
                        drC = NuclRadius + AuNPRadius + G4UniformRand()*(MaxDistFromNucleus - 2*AuNPRadius);
                    else
                        drC = NuclRadius + AuNPRadius + G4UniformRand()*(CellRadius - NuclRadius - 2*AuNPRadius);
                    
                    xC = (AuNPRadius + drC)* std::cos(uC) * std::sin(vC);
                    yC = (AuNPRadius + drC)* std::sin(uC) * std::sin(vC);
                    zC = (AuNPRadius + drC)* std::cos(vC);
        
                    if (!intersectC){
                        for (G4int countA=0;countA<jA;countA++){
                            dist_sq_C = (xC - NPpositionsA[countA][0])*(xC - NPpositionsA[countA][0]) + 
                                (yC - NPpositionsA[countA][1])*(yC - NPpositionsA[countA][1]) + 
                                (zC - NPpositionsA[countA][2])*(zC - NPpositionsA[countA][2]);

                            if (dist_sq_C < 4*AuNPRadius*AuNPRadius){
                                intersectC = true;
                                break;
                            } 
                        }
                    }

                    if (!intersectC && sector > 0){
                        for (G4int countB=0;countB<statC;countB++){
                            dist_sq_C = (xC - NPpositionsB[countB][0])*(xC - NPpositionsB[countB][0]) + 
                                (yC - NPpositionsB[countB][1])*(yC - NPpositionsB[countB][1]) + 
                                (zC - NPpositionsB[countB][2])*(zC - NPpositionsB[countB][2]);

                            if (dist_sq_C < 4*AuNPRadius*AuNPRadius){
                                intersectC = true;
                                break;
                            } 
                        }
                    }

                    if (!intersectC){
                        G4ThreeVector* positionC = new G4ThreeVector(xC,yC,zC);

                        G4RotationMatrix* rotationC = new G4RotationMatrix(0,0,0);

                        G4VPhysicalVolume* pAuNPC =
                            CreatePhysicalVolume(
                                subComponentName3, jA, true, lAuNPC, rotationC , positionC, fEnvelopePhys);

                        NPpositionsA[jA][0] = xC;
                        NPpositionsA[jA][1] = yC;
                        NPpositionsA[jA][2] = zC;
                        if (sector == 0){
                            NPpositions0[jA][0] = xC;
                            NPpositions0[jA][1] = yC;
                            NPpositions0[jA][2] = zC;
                        
                        }

                        if (jA%100 == 99)
                            G4cout << (2*statC*sector + jA + 1) << "NPs placed." << G4endl;

                        if (ThereIsACore){
                            G4VPhysicalVolume* pAuNPCore =
                            CreatePhysicalVolume(subsubComponentName, lCore, pAuNPC);
                        }
                    }
                }
            }

            for (G4int jB = 0; jB < statC; jB++){

                intersectC = true;

                while (intersectC){
          
                    intersectC = false;

                    G4double uC = 2*pi*2*sector/100 + 2*pi/100 + G4UniformRand()*2*pi/100;
                    G4double vC = std::acos(2*G4UniformRand()-1);
                    G4double drC;
                    if (MaxDistExists)
                        drC = NuclRadius + AuNPRadius + G4UniformRand()*(MaxDistFromNucleus - 2*AuNPRadius);
                    else
                        drC = NuclRadius + AuNPRadius + G4UniformRand()*(CellRadius - NuclRadius - 2*AuNPRadius);
                    
                    xC = (AuNPRadius + drC)* std::cos(uC) * std::sin(vC);
                    yC = (AuNPRadius + drC)* std::sin(uC) * std::sin(vC);
                    zC = (AuNPRadius + drC)* std::cos(vC);
        
                    if (!intersectC){
                        for (G4int countB=0;countB<jB;countB++){
                            dist_sq_C = (xC - NPpositionsB[countB][0])*(xC - NPpositionsB[countB][0]) + 
                                (yC - NPpositionsB[countB][1])*(yC - NPpositionsB[countB][1]) + 
                                (zC - NPpositionsB[countB][2])*(zC - NPpositionsB[countB][2]);

                            if (dist_sq_C < 4*AuNPRadius*AuNPRadius){
                                intersectC = true;
                                break;
                            } 
                        }
                    }
        
                    if (sector == 49 && !intersectC){
                        for (G4int count0=0;count0<statC;count0++){
                            dist_sq_C = (xC - NPpositions0[count0][0])*(xC - NPpositions0[count0][0]) + 
                                (yC - NPpositions0[count0][1])*(yC - NPpositions0[count0][1]) + 
                                (zC - NPpositions0[count0][2])*(zC - NPpositions0[count0][2]);

                            if (dist_sq_C < 4*AuNPRadius*AuNPRadius){
                                intersectC = true;
                                break;
                            } 
                        }
                    }

                    if (!intersectC){
                        for (G4int countA=0;countA<statC;countA++){
                            dist_sq_C = (xC - NPpositionsA[countA][0])*(xC - NPpositionsA[countA][0]) + 
                                (yC - NPpositionsA[countA][1])*(yC - NPpositionsA[countA][1]) + 
                                (zC - NPpositionsA[countA][2])*(zC - NPpositionsA[countA][2]);

                            if (dist_sq_C < 4*AuNPRadius*AuNPRadius){
                                intersectC = true;
                                break;
                            } 
                        }
                    }

                    if (!intersectC){
                        G4ThreeVector* positionC = new G4ThreeVector(xC,yC,zC);

                        G4RotationMatrix* rotationC = new G4RotationMatrix(0,0,0);

                        G4VPhysicalVolume* pAuNPC =
                            CreatePhysicalVolume(
                                subComponentName3, jB, true, lAuNPC, rotationC , positionC, fEnvelopePhys);
                        NPpositionsB[jB][0] = xC;
                        NPpositionsB[jB][1] = yC;
                        NPpositionsB[jB][2] = zC;

                        if (jB%100 == 99)
                            G4cout << (2*statC*sector + statC + jB + 1) << "NPs placed." << G4endl;

                        if (ThereIsACore){
                            G4VPhysicalVolume* pAuNPCore =
                            CreatePhysicalVolume(subsubComponentName, lCore, pAuNPC);
                        }
                    }
                }
            }
        }
    }



	
    //****************************************************************************
    //****************************************************************************
    //                ===              Nucleus             ===                  //
    //****************************************************************************
    //****************************************************************************

    if (fNucleusExists){

    //****************************************************************************
    //                     Read  Hilbert space filling  data                    //
    //****************************************************************************
    FileName = "./supportFiles/"+FileName;
    const char* filename = FileName;
    std::string line = "";
    ifstream f(filename, ios::in);
    G4double x, y, z;
 
    if (f.is_open()) {
        while (f >> x >> y >> z){
            fx.push_back(x);
            fy.push_back(y);
            fz.push_back(z);
        }
    }
    else {
        G4cout << "ERROR: Unable to open file " << FileName << G4endl;
        exit(1);
    }

    G4int TotalPoints = fx.size();
    G4cout << "Number of points " << TotalPoints << G4endl;

    //****************************************************************************
    //                               Set geometry size                          //
    //****************************************************************************
    
    G4double FiberEnvelopeRadius = 18.55*nm;    //Fiber radius
    G4double FiberEnvelopeLength = 120 *nm; 
    G4double HiberterPointDistance = 40*nm+FiberEnvelopeLength;

    G4int HilbertFold = 1;
    if(TotalPoints==8) HilbertFold = 1;
    if(TotalPoints==64) HilbertFold = 3;
    if(TotalPoints==512) HilbertFold = 7;
    if(TotalPoints==4096) HilbertFold = 15;
    VoxelLength = HiberterPointDistance*HilbertFold + HilbertCurveLayer*FiberEnvelopeRadius*2 + ToleranceSpace; //add tolerance
    ParaContainerHalfSize = VoxelLength*HilbertCurve3DRepeat/2;

    G4cout<<"*********************************************************************************"<<G4endl;
    G4cout<<"The Hilbert curve data file name: "<<FileName<<G4endl;
    G4cout<<"The nuclues was devided into "<<HilbertCurve3DRepeat<<" subdomains(voxels) on X/Y/Z diection respectively."<<G4endl;
    G4cout<<"Voxel size :"<< VoxelLength/um<<" um "<<G4endl;
    G4cout<<"Voxel container size :"<< ParaContainerHalfSize*2/um<<" um "<<G4endl;
    G4cout<<"Reuse the Hilbert curve "<<HilbertCurveLayer<<" times to fill each subdomain(voxels)."<<G4endl;
    G4cout<<"HiberterPointDistance = "<<HiberterPointDistance/um<<" um"<<G4endl;
    G4cout<<"*********************************************************************************"<<G4endl;

    //****************************************************************************
    //                                Parameterise                              //
    //****************************************************************************
    SetBasicInfo();

    ////----- Create parameterisation and set
    param = new VoxelParameterisation();
    param->SetVoxelDimensions( VoxelLength/2, VoxelLength/2, VoxelLength/2 ); 
    param->SetNoVoxel( HilbertCurve3DRepeat, HilbertCurve3DRepeat, HilbertCurve3DRepeat ); 
    param->SetContainerDimensions(ParaContainerHalfSize, ParaContainerHalfSize, ParaContainerHalfSize);

    //****************************************************************************
    //                             Build basic geometry                         //
    //****************************************************************************
    // Nucleus box 
    //G4Box* Nucleus_solid = new G4Box( "nucl",ParaContainerHalfSize, ParaContainerHalfSize, ParaContainerHalfSize);
    G4Orb * Nucleus_solid = new G4Orb("nucl", NuclRadius);
    G4LogicalVolume* fNuclEnvelopeLog= CreateLogicalVolume(Nucleus_solid);
    G4Colour  blue    (0.0, 0.0, 1.0) ;
    G4VisAttributes* Vis = new G4VisAttributes( blue );
    Vis->SetVisibility(false);
    fNuclEnvelopeLog->SetVisAttributes(Vis);
    

    // Fiber Envelope
    G4double length = std::sqrt( pow(fx[2]-fx[1],2) + pow(fy[2]-fy[1],2) + pow(fz[2]-fz[1],2));
    G4double scaleFactor = HiberterPointDistance/length; // scaleFactor used to get fiber lengths
    length = scaleFactor*length;
    G4Tubs* FiberEnvelopeSolid = new G4Tubs("Chromoloop", 0, FiberEnvelopeRadius, FiberEnvelopeLength/2, 0*deg, 360*deg);
    FiberEnvelopeLogic = CreateLogicalVolume(FiberEnvelopeSolid);
    G4VisAttributes* FiberEnvelopeVis = new G4VisAttributes( G4Colour(0.6,0.0,0.4) ); //
    FiberEnvelopeVis->SetVisibility(true);
    FiberEnvelopeVis->SetForceWireframe(true);
    FiberEnvelopeVis->SetForceAuxEdgeVisible(true);
    FiberEnvelopeLogic->SetVisAttributes(FiberEnvelopeVis);

    // Logic volume of fiber
    if(!testmode){
        fAddHydrationShell = AddHydrationShell;
        fHydrationShellThickness = HydrationShellThickness;
        fChromatinRadius = FibreRadius;
        fChromatinLength = FibreLength/2;
        FiberLogic = Fibre_Construct();
    }

 
    //****************************************************************************
    //                             Build loop geometry                          //
    //****************************************************************************
    G4int maxSize = HilbertCurveLayer*TotalPoints*(G4int)pow(HilbertCurve3DRepeat,3);
    pvLoop.resize(maxSize);
    lvLoop.resize(maxSize);
    FibreMidpoints.resize(maxSize);
    FibreDirections.resize(maxSize);

    //----- Define voxel logical volume
    G4Box           * voxel_solid = new G4Box( "Voxel", voxelHalfDimX, voxelHalfDimY, voxelHalfDimZ);
    G4LogicalVolume * voxel_logic = CreateLogicalVolume(voxel_solid); 
    G4Colour  yellow  (1.0, 1.0, 0.0) ;
    G4VisAttributes* voxelVis = new G4VisAttributes( G4Colour(0.0,0.8,0.7) ); //
    voxelVis->SetVisibility(true);
    voxelVis->SetForceWireframe(true);
    voxelVis->SetForceAuxEdgeVisible(true);
    voxel_logic->SetVisAttributes(voxelVis);


    CountFibers = 0;
    G4double xshift = -FiberEnvelopeRadius*(HilbertCurveLayer-1);
    G4double yshift = -FiberEnvelopeRadius*(HilbertCurveLayer-1);
    G4double zshift = -FiberEnvelopeRadius*(HilbertCurveLayer-1);

    G4ThreeVector* midpoint  = new G4ThreeVector();

    G4String current_dir;

    for (G4int jj = 0; jj < HilbertCurveLayer; jj++)
    {
        //G4cout << "jj = " << jj << G4endl;
        G4int max_loopthelayer = fx.size();
        for (G4int loopthelayer = 1; loopthelayer < max_loopthelayer; loopthelayer++) 
        {
            //G4cout << "loopthelayer = " << loopthelayer << G4endl;
            G4double layershift = FiberEnvelopeRadius*2.01;
            G4int i=0;
            if(jj%2==0) i= loopthelayer;
            if(jj%2==1) i= max_loopthelayer-loopthelayer; 
    
            G4double midpoint_x = jj*layershift + xshift + (fx[i]+fx[i-1])/2*scaleFactor;
            G4double midpoint_y = jj*layershift + yshift + (fy[i]+fy[i-1])/2*scaleFactor;
            G4double midpoint_z = jj*layershift + zshift + (fz[i]+fz[i-1])/2*scaleFactor;
            midpoint->setX(midpoint_x);
            midpoint->setY(midpoint_y);
            midpoint->setZ(midpoint_z);

            G4ThreeVector MidP = G4ThreeVector(midpoint_x,midpoint_y,midpoint_z);
            FibreMidpoints[CountFibers] = MidP;

            G4RotationMatrix* rotLoop = new G4RotationMatrix();
            G4ThreeVector direction = G4ThreeVector(fx[i]-fx[i-1],fy[i]-fy[i-1],fz[i]-fz[i-1]);

            if ((direction.x() == 0) && (direction.y() == 0)){
                rotLoop->rotateZ(0*deg);
                current_dir = "z";
            }
            else if ((direction.y() == 0) && (direction.z() == 0)){
                rotLoop->rotateY(90*deg); //rotate along x
                current_dir = "x";
            }
            else if ((direction.x() == 0) && (direction.z() == 0)){
                rotLoop->rotateX(90*deg);
                current_dir = "y";
            }
            else G4cout << "Two consecutive Hilbert paths have the same direction." << G4endl;

            FibreDirections[CountFibers] = current_dir;


            if(testmode){
                G4String volumeName = "Chromoloop";
                // int2str(CountFibers, volumeName);
                // volumeName = "Chromoloop"+volumeName;
                
                // lvLoop[CountFibers] = CreateLogicalVolume(FiberEnvelopeSolid);
                pvLoop[CountFibers] = CreatePhysicalVolume(
                                                        volumeName,        //name
                                                        CountFibers,       //copy number
                                                        true,              //many
                                                        FiberEnvelopeLogic,//logical volume
                                                        rotLoop,      // rotation
                                                        midpoint,      //translation 
                                                        voxel_logic); // mother
            }
            else{
                G4String volumeName ;
                int2str(CountFibers, volumeName);
                volumeName = "Fiber"+volumeName;
                
                // lvLoop[CountFibers] = CreateLogicalVolume(FiberEnvelopeSolid);
                pvLoop[CountFibers] = CreatePhysicalVolume(
                                                        volumeName,        //name
                                                        CountFibers,       //copy number
                                                        true,              //many
                                                        FiberLogic,//logical volume
                                                        rotLoop,      // rotation
                                                        midpoint,      //translation 
                                                        voxel_logic); // mother

            }
            
            if(CheckOverlap)
            {
                //pvLoop[CountFibers]->CheckOverlaps();
                G4cout<<"checking CountFibers="<<CountFibers<<G4endl;
                if( pvLoop[CountFibers]->CheckOverlaps(1000, 0, false)) 
                {
                    ofstream outfile;
                    outfile.open("overlap.txt",std::ios::out|std::ios::app);
                    outfile<<"In "<<" th subdomain "<<jj+1<<" th layer "<<CountFibers<<" th volume detected overlap\n";
                    outfile.close(); 
                }

            }

            CountFibers++;
        }
    }

    //*******************************
    // Subcomponent: (hybrid) AuNPs inside nucleus
    //*******************************

    G4int NbOfAuNPsInNucleus = 0;
    if (fPm->ParameterExists(GetFullParmName("AuNPs/NbOfAuNPsInNucleus")))
        NbOfAuNPsInNucleus = fPm->GetIntegerParameter(GetFullParmName("AuNPs/NbOfAuNPsInNucleus"));

    if (NbOfAuNPsInNucleus > 0){

        //Number of AuNPs in the nucleus
        NbOfAuNPsInNucleus = fPm->GetIntegerParameter(GetFullParmName("AuNPs/NbOfAuNPsInNucleus"));
        
        //Radius of the AuNPs (default values if none are specified)
        G4double AuNPRadius = 2*nanometer;
        
        if (fPm->ParameterExists(GetFullParmName("AuNPs/Rmax")))
            {AuNPRadius = fPm->GetDoubleParameter(GetFullParmName("AuNPs/Rmax"), "Length" );}

        G4Orb* gAuNP = new G4Orb("gAuNP", AuNPRadius);
        lAuNPN = CreateLogicalVolume("AuNPs", gAuNP);

        G4bool ThereIsACore = false;
        if (fPm->ParameterExists(GetFullParmName("AuNPs/Rcore"))) {
            ThereIsACore = true;

            G4String symbol;
            G4int n_atoms,n_comps;
            G4double z,w,d_Fe3O4;
            G4Element* iron =
            new G4Element("Iron", symbol="Fe", z=26, w=55.845*g/mole);
            G4Element* oxygen =
            new G4Element("Oxygen", symbol="O", z=8, w=15.9994*g/mole);
            G4Material * Fe3O4 =
            new G4Material("Ironoxide", d_Fe3O4=5.17*g/cm3,n_comps=2);

            Fe3O4->AddElement(iron, n_atoms=3);
            Fe3O4->AddElement(oxygen, n_atoms=4);

            G4double Rcore = fPm->GetDoubleParameter(GetFullParmName("AuNPs/Rcore"), "Length" );
            G4Orb* gCore = new G4Orb("gCore",Rcore);
            lCore = new G4LogicalVolume(gCore,Fe3O4,"CoreLogical",0,0,0);
        }

        //Randomly distribute hybrid AuNPs throughout nucleus
        
        G4double xN,yN,zN;
        G4double dist_sq_N;
        G4ThreeVector current_midpt;
        G4double xM,yM,zM;
        G4double circtest;
        G4double disttest;

        NPpositions.clear();
        G4double EnvelopeNucleusRatio = pow(ParaContainerHalfSize*2,3)/(pow(NuclRadius,3)*3.14159*4/3);
        G4int AuNPsPerVoxel = ceil(EnvelopeNucleusRatio*NbOfAuNPsInNucleus/pow(HilbertCurve3DRepeat,3));

        G4cout << "Inserting " << AuNPsPerVoxel << " nanoparticles per Voxel." << G4endl;
        G4cout<<"*********************************************************************************"<<G4endl;

        G4ThreeVector current_pos;
        NPLoop.resize(AuNPsPerVoxel);
        G4double reducedLength = voxelHalfDimX - AuNPRadius;

        for (G4int jN = 0; jN < AuNPsPerVoxel; jN++){

            G4bool intersectN = true;

            while (intersectN){
      
                intersectN = false;

                xN = (G4UniformRand()*reducedLength*2 - reducedLength);
                yN = (G4UniformRand()*reducedLength*2 - reducedLength);
                zN = (G4UniformRand()*reducedLength*2 - reducedLength);

                // check for overlapse with each fibre
                for (G4int lF=0;lF<CountFibers;lF++){
                    current_midpt = FibreMidpoints[lF];
                    xM = current_midpt[0];
                    yM = current_midpt[1];
                    zM = current_midpt[2];

                    current_dir = FibreDirections[lF];

                    if (current_dir == "x"){
                        circtest = sqrt(pow((yN - yM),2) + pow((zN - zM),2));
                        if (circtest < (AuNPRadius + FiberEnvelopeRadius)){
                            disttest = sqrt(pow((xN - xM),2));
                            if (disttest < (AuNPRadius + FiberEnvelopeLength)){
                                intersectN = true;
                                break;
                            }

                        }
                    }

                    if (current_dir == "y"){
                        circtest = sqrt(pow((xN - xM),2) + pow((zN - zM),2));
                        if (circtest < (AuNPRadius + FiberEnvelopeRadius)){
                            disttest = sqrt(pow((yN - yM),2));
                            if (disttest < (AuNPRadius + FiberEnvelopeLength)){
                                intersectN = true;
                                break;
                            }

                        }
                    }

                    if (current_dir == "z"){
                        circtest = sqrt(pow((xN - xM),2) + pow((yN - yM),2));
                        if (circtest < (AuNPRadius + FiberEnvelopeRadius)){
                            disttest = sqrt(pow((zN - zM),2));
                            if (disttest < (AuNPRadius + FiberEnvelopeLength)){
                                intersectN = true;
                                break;
                            }

                        }
                    }

                }

                // check for overlapse with other AuNPs
                for (G4int lN=0;lN<jN;lN++){
                    current_pos = NPpositions[lN];
                    dist_sq_N = (xN - current_pos[0])*(xN - current_pos[0]) + 
                        (yN - current_pos[1])*(yN - current_pos[1]) + 
                        (zN - current_pos[2])*(zN - current_pos[2]);

                    if (dist_sq_N < 4*AuNPRadius*AuNPRadius){
                        intersectN = true;
                        break;
                    }
                }

                if (!intersectN){
                    G4ThreeVector* positionN = new G4ThreeVector(xN,yN,zN);
                    G4RotationMatrix* rotationN = new G4RotationMatrix(0,0,0);
                    G4String aunp_name = "AuNPN";

                    NPLoop[jN] = CreatePhysicalVolume(aunp_name,    // String identifier for this placement
                                                    jN,             // Integer which identifies this placement
                                                    false,          // For future use. Can be set to false
                                                    lAuNPN,         // The associated Logical Volume
                                                    rotationN,      // Rotation with respect to its mother volume
                                                    positionN,      // Translation with respect to its mother volume
                                                    voxel_logic     // The associated mother volume
                                                    );

                    NPpositions.push_back(G4ThreeVector(xN,yN,zN));

                    if (ThereIsACore){
                        G4VPhysicalVolume* pAuNPCore =
                        CreatePhysicalVolume("AuNPCore", lCore, NPLoop[jN]);
                    }
                }
            }
        }
    }


    G4VPhysicalVolume* fNuclEnvelopePhys = CreatePhysicalVolume("Nucleus",         // name
                                            0,                  // copy number
                                            false,              // many
                                            fNuclEnvelopeLog,   // logical volume
                                            NucleusRotation,    // rotation
                                            NucleusTranslation, // translation 
                                            fEnvelopePhys);     // mother

    G4VPhysicalVolume * DNAContent = CreatePhysicalVolume("DNAContent", voxel_logic, fNuclEnvelopePhys, kXAxis, fnVoxels, param);  
    
    // //----- Define empty volume. for test 
    // G4Box              * emptyvoxel_solid = new G4Box( "Voxel", voxelHalfDimX, voxelHalfDimY, voxelHalfDimZ);
    // G4LogicalVolume    * emptyvoxel_logic = new G4LogicalVolume(emptyvoxel_solid,fMaterials[0],"emptyVoxelLogical",0,0,0); 
    // G4VPhysicalVolume  * DNAContent = new G4PVParameterised("DNAContent",emptyvoxel_logic,fNuclEnvelopeLog,kXAxis, fnVoxels, param); 
    // emptyvoxel_logic->SetVisAttributes(voxelVis); 

    }

    InstantiateChildren(fEnvelopePhys);
    return fEnvelopePhys;
}



void  TsEllipsoidCellNAK::SetBasicInfo()
{

    fnVoxelX=HilbertCurve3DRepeat; 
    fnVoxelY=HilbertCurve3DRepeat;
    fnVoxelZ=HilbertCurve3DRepeat; 
    fnVoxels = fnVoxelX*fnVoxelY*fnVoxelZ;
    fvoxelDimX = VoxelLength; 
    fvoxelDimY = VoxelLength; 
    fvoxelDimZ = VoxelLength;

    voxelHalfDimX = fvoxelDimX/2;
    voxelHalfDimY = fvoxelDimY/2;
    voxelHalfDimZ = fvoxelDimY/2;

    fminX = -fnVoxelX*fvoxelDimX/2; fmaxX = fnVoxelX*fvoxelDimX/2;
    fminY = -fnVoxelY*fvoxelDimY/2, fmaxY = fnVoxelY*fvoxelDimY/2;
    fminZ = -fnVoxelZ*fvoxelDimZ/2; fmaxZ = fnVoxelZ*fvoxelDimZ/2;

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* water  = nist->FindOrBuildMaterial("G4_WATER");
    G4Material* H2O_mod= nist -> BuildMaterialWithNewDensity("G4_WATER_Denser", "G4_WATER", 1.407 *g/cm/cm/cm);

    fMaterials.push_back(water); 
    fMaterials.push_back(H2O_mod); 

}

G4LogicalVolume* TsEllipsoidCellNAK::Fibre_Construct()
{
    BuildSphere = false;
    G4NistManager* nist = G4NistManager::Instance();
    water  = nist->FindOrBuildMaterial("G4_WATER");
    mat_Histone            = nist -> BuildMaterialWithNewDensity("mat_Histone", "G4_WATER", 1. *g/cm/cm/cm);
    mat_DNA_base           = nist -> BuildMaterialWithNewDensity("mat_DNA_base", "G4_WATER", 1. *g/cm/cm/cm);
    mat_DNA_backbone       = nist -> BuildMaterialWithNewDensity("mat_DNA_backbone", "G4_WATER", 1.407 *g/cm/cm/cm);
    mat_DNA_hydrationShell = nist -> BuildMaterialWithNewDensity("mat_DNA_hydrationShell", "G4_WATER", 1. *g/cm/cm/cm);

    G4double ChromatinEnvelopeRadius = (37.1/2.0)*nm;
    G4double ChromatinEnvelopeLength = (120/2.0)*nm ;

    //****************************************************************************
    //                   Cylinder containing Histone & DNA (envelope)           //
    //****************************************************************************

    G4Tubs* gCylinder = new G4Tubs("Chromatin",
                                   0,
                                   ChromatinEnvelopeRadius,
                                   ChromatinEnvelopeLength,
                                   180.0 * deg,
                                   360.0 * deg);

    aEnvelopeLog = new G4LogicalVolume( gCylinder,water,"Chromatin",0,0,0);   


    //vis
    G4VisAttributes* EnvelopeVis = new G4VisAttributes( G4Colour(0.6,0.0,0.4) ); //
    EnvelopeVis->SetVisibility(true);
    //EnvelopeVis->SetForceWireframe(true);
    EnvelopeVis->SetForceAuxEdgeVisible(true);
    aEnvelopeLog->SetVisAttributes(EnvelopeVis);


    //****************************************************************************
    //                              Histones                                    //
    //****************************************************************************
    //Record the position and rotation of each Histone placed
    vector<pair<G4ThreeVector, G4RotationMatrix*>> HistoneDetails;
    Fibre_BuildHistones(HistoneDetails, FibreRadius, FibreLength/2);


    //****************************************************************************
    //                                DNA                                       //
    //****************************************************************************
    //Builds DNA around and between the Histones
    Fibre_BuildDNA(HistoneDetails);

    return aEnvelopeLog;
}

//build the Histones
void TsEllipsoidCellNAK::Fibre_BuildHistones(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                            G4double ChromatinRadius,
                            G4double ChromatinLength)
{
    //****************************************************************************
    //                              Histones
    //****************************************************************************

    //place the histones in a left handed solenoid conformation

    G4int nb_HistPerTurn = 6;
    G4double HistoneRadius = 3.3*nm;
    G4double HistoneLength = (5.7 / 2.0)*nm;

    G4String SubComp = "Histone";

    G4Tubs* gHistone = new G4Tubs("Histone",
                                  0,
                                  HistoneRadius,
                                  HistoneLength,
                                  180*deg,
                                  360*deg);


    G4LogicalVolume* lHistone = new G4LogicalVolume( gHistone,mat_Histone,SubComp,0,0,0); 

    HistoneRadius+=3.0*nm;  //extra spacing (3nm) for the double helix


    //Generate path for nucleosomes -- solenoid
    G4double x = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * sin(0.0));
    G4double y = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * cos(0.0));
    G4double z = (-ChromatinLength + HistoneRadius);

    G4int nNuc = (6.0*2.0*(-z)) / (2.0*HistoneRadius); //how many histones fit along chromatin z
    G4double zStep = (2.0*HistoneRadius)/nb_HistPerTurn;
    G4ThreeVector position (x,y,z);

    G4double theta = 0.0;
    G4double thetaStep = (2.0*pi)/(G4double)nb_HistPerTurn;

    G4int built=0;
    for (G4int i=0;i<nNuc;i++){
        theta+=thetaStep;
        position[0]=((ChromatinRadius-HistoneRadius) * sin(theta));
        position[1]=((ChromatinRadius-HistoneRadius) * cos(theta));
        position[2]=z;
        z+=zStep;

        if (position[2]+HistoneRadius>=ChromatinLength){
            nNuc=built;
            break;
        }

        G4RotationMatrix *HistoneRotation = new G4RotationMatrix();
        HistoneRotation->rotateZ((-120.0+((G4double)(i)*360.0/(G4double)nb_HistPerTurn))*deg);
        HistoneRotation->rotateY(90.0*deg);

        //record the position and rotation of the histone
        pair<G4ThreeVector,G4RotationMatrix*> Details;
        Details.first=position;
        Details.second=HistoneRotation;
        HistoneDetails.push_back(Details);

        G4ThreeVector* pos_hist = new G4ThreeVector();
        pos_hist->setX(position[0]);
        pos_hist->setY(position[1]);
        pos_hist->setZ(position[2]);

        G4String histone_name = "Histone";

        CreatePhysicalVolume(
                            histone_name,          // String identifier for this placement
                            i,                  // Integer which identifies this placement
                            false,               // For future use. Can be set to false
                            lHistone,           // The associated Logical Volume
                            HistoneRotation,    // Rotation with respect to its mother volume
                            pos_hist,           // Translation with respect to its mother volume
                            aEnvelopeLog        // The associated mother volume
                            );

        built++;
    }

    G4VisAttributes* HistoneVis = new G4VisAttributes(G4Colour(0.,0.6,0.));
    HistoneVis->SetVisibility(true);
    HistoneVis->SetForceSolid(true);
    lHistone->SetVisAttributes(HistoneVis);
}

//build the DNA
void TsEllipsoidCellNAK::Fibre_BuildDNA(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails)
{
    //****************************************************************************
    //                                DNA
    //****************************************************************************

    G4bool BuildHalfCyl=false;
    G4bool BuildQuartCyl=false;
    G4bool BuildSphere=false;

    G4String DNAModel = "HalfCyl";

    if (DNAModel=="HalfCyl"){BuildHalfCyl=true;}
    else if (DNAModel=="Sphere"){BuildSphere=true;}
    else if (DNAModel=="QuartCyl"){BuildQuartCyl=true;}

    Fibre_SetDNAVolumes(BuildHalfCyl, BuildQuartCyl, BuildSphere);

    vector<G4ThreeVector> DNAPath;
    Fibre_GenerateDNAPath(HistoneDetails, DNAPath);

    Fibre_SegmentDNAPath(DNAPath);

    if (BuildSphere){
        Fibre_PlaceDNASphere(DNAPath);

    } else {
        Fibre_PlaceDNA(DNAPath);
    }

}

//set up DNA volumes
void TsEllipsoidCellNAK::Fibre_SetDNAVolumes(G4bool BuildHalfCyl,
                            G4bool BuildQuartCyl,
                            G4bool BuildSphere)
{
    // A choice to build the DNA volumes as Half Cylinders, Quarter Cylinders, or Spheres
    // Quarter cylinders by default
    // A denser water for DNA backbone

    //SubCompartment names
    G4String back1 = "back1";
    G4String back2 = "back2";
    G4String base1 = "base1";
    G4String base2 = "base2";

    //sphere DNA
    if (BuildSphere){
        G4Sphere* gDNA_base = new G4Sphere("DNA_base",
                                           0*nm,
                                           0.208*nm,
                                           0*deg,
                                           360*deg,
                                           0*deg,
                                           180*deg);

        lBase1  = new G4LogicalVolume( gDNA_base,mat_DNA_base,base1,0,0,0); 
        lBase2  = new G4LogicalVolume( gDNA_base,mat_DNA_base,base2,0,0,0); 


        G4Sphere* gDNA_backbone = new G4Sphere("DNA_backbone",
                                               0*nm,
                                               0.24*nm,
                                               0*deg,
                                               360*deg,
                                               0*deg,
                                               180*deg);

        lBack1  = new G4LogicalVolume( gDNA_backbone,mat_DNA_backbone,back1,0,0,0); 
        lBack2  = new G4LogicalVolume( gDNA_backbone,mat_DNA_backbone,back2,0,0,0); 

        // ************************** build water layer **************************
        G4Sphere* gWater = new G4Sphere("DNA_WaterLayer",
                                               0*nm,
                                               0.24*nm+fHydrationShellThickness,
                                               0*deg,
                                               360*deg,
                                               0*deg,
                                               180*deg);

        lWater1 = new G4LogicalVolume( gWater,mat_DNA_hydrationShell,"water1",0,0,0); 
        lWater2 = new G4LogicalVolume( gWater,mat_DNA_hydrationShell,"water2",0,0,0); 
    }
    //half cylinder
    else if (BuildHalfCyl){
        G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        180.*deg,
                                        180.*deg);
        G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        0.*deg,
                                        180.*deg);

        lBase1  = new G4LogicalVolume( gDNA_base1,mat_DNA_base,base1,0,0,0); 
        lBase2  = new G4LogicalVolume( gDNA_base2,mat_DNA_base,base2,0,0,0); 

        G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            180.*deg);
        G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            180.*deg);

        lBack1  = new G4LogicalVolume( gDNA_backbone1,mat_DNA_backbone,back1,0,0,0); 
        lBack2  = new G4LogicalVolume( gDNA_backbone2,mat_DNA_backbone,back2,0,0,0); 

         // ************************** build water layer **************************
         G4Tubs* gWater1 = new G4Tubs("WaterLayer1",
                                            1.15*nm,
                                            1.15*nm+fHydrationShellThickness,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            180.*deg);
        G4Tubs* gWater2 = new G4Tubs("WaterLayer2",
                                            1.15*nm,
                                            1.15*nm+fHydrationShellThickness,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            180*deg);

        lWater1 = new G4LogicalVolume( gWater1,mat_DNA_hydrationShell,"water1",0,0,0); 
        lWater2 = new G4LogicalVolume( gWater2,mat_DNA_hydrationShell,"water2",0,0,0); 

    }
    //quarter cylinder
    else if (BuildQuartCyl){
        G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        135.*deg,
                                        180.*deg);
        G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        -45.*deg,
                                        180.*deg);

        lBase1  = new G4LogicalVolume( gDNA_base1,mat_DNA_base,base1,0,0,0); 
        lBase2  = new G4LogicalVolume( gDNA_base2,mat_DNA_base,base2,0,0,0); 

        G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            90.*deg);
        G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            90*deg);

        lBack1  = new G4LogicalVolume( gDNA_backbone1,mat_DNA_backbone,back1,0,0,0); 
        lBack2  = new G4LogicalVolume( gDNA_backbone2,mat_DNA_backbone,back2,0,0,0); 

        // ************************** build water layer **************************
         G4Tubs* gWater1 = new G4Tubs("WaterLayer1",
                                            1.15*nm,
                                            1.15*nm+fHydrationShellThickness,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            90.*deg);
        G4Tubs* gWater2 = new G4Tubs("WaterLayer2",
                                            1.15*nm,
                                            1.15*nm+fHydrationShellThickness,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            90*deg);

        lWater1 = new G4LogicalVolume( gWater1,mat_DNA_hydrationShell,"water1",0,0,0); 
        lWater2 = new G4LogicalVolume( gWater2,mat_DNA_hydrationShell,"water1",0,0,0); 

    }


    //vis
    G4VisAttributes* Back1Vis = new G4VisAttributes(G4Colour(1.0, 0.1, 0.02));
    Back1Vis->SetVisibility(true);
    Back1Vis->SetForceSolid(true);
    lBack1->SetVisAttributes(Back1Vis);

    G4VisAttributes * Back2Vis = new G4VisAttributes(G4Colour(0.0, 0.2, 1.0));
    Back2Vis->SetVisibility(true);
    Back2Vis->SetForceSolid(true);
    lBack2 -> SetVisAttributes(Back2Vis);

    G4VisAttributes * Base1Vis = new G4VisAttributes(G4Colour(0.02, 0.92, 1.0));
    Base1Vis->SetVisibility(true);
    Base1Vis->SetForceSolid(true);
    lBase1 -> SetVisAttributes(Base1Vis);

    G4VisAttributes * Base2Vis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.8));
    Base2Vis->SetVisibility(true);
    Base2Vis->SetForceSolid(true);
    lBase2 -> SetVisAttributes(Base2Vis);

    if(fAddHydrationShell)    {

    G4VisAttributes * Water1Vis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    Water1Vis->SetVisibility(true);
    Water1Vis->SetForceSolid(true);
    lWater1 -> SetVisAttributes(Water1Vis);

    G4VisAttributes * Water2Vis = new G4VisAttributes(G4Colour(0.0, 0.4, 1.0));
    Water2Vis->SetVisibility(true);
    Water2Vis->SetForceSolid(true);
    lWater2 -> SetVisAttributes(Water2Vis);
    }

}

//generate a path for DNA around the histones
void TsEllipsoidCellNAK::Fibre_GenerateDNAPath(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                              vector<G4ThreeVector> &path)
{
    //generate a path of 3vector points that defines the centre of the dbl helix

    G4double TraceBP=15.0; //the number of BP to trace back for bezier curve (more = bendier or smoother)
    G4double nCoils = 1.65; //#coils around histone
    G4int nHistones=HistoneDetails.size();
    G4double offset = (pi/2.0) - (2.0*pi*nCoils);
    G4double zHist = (5.7 / 2.0)*nm;
    G4double Radius = 3.3*nm + 1.15*nm + 0.1*nm; //central radius of double helix around histone

    G4int nSteps=200; //go around the histone in 200 steps (segmented later)

    for (G4int i=0;i<nHistones;i++){
        G4RotationMatrix *rot=HistoneDetails[i].second;

        //histone DNA
        for (G4int j=0;j<nSteps;j++){
            G4double angle = ((G4double)j * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
            //take a pt around the histone
            G4ThreeVector pt( (Radius * sin(angle)),
                             (Radius * cos(angle)),
                             (-zHist + ((G4double)j/(G4double)nSteps)*2.0*(zHist)));

            //rotate the step to the histone
            G4ThreeVector RotPt;
            RotPt[0] = pt.x()*rot->xx() + pt.y()*rot->yx() + pt.z()*rot->zx();
            RotPt[1] = pt.x()*rot->xy() + pt.y()*rot->yy() + pt.z()*rot->zy();
            RotPt[2] = pt.x()*rot->xz() + pt.y()*rot->yz() + pt.z()*rot->zz();

            //translate to the histone position
            RotPt += HistoneDetails[i].first;

            path.push_back(RotPt);
        }


        //link to nextPt -- off histone (start from end of histone i -> start of histone i+1)
        if (i!=nHistones-1){
            G4ThreeVector start = path[path.size()-1]; //last pt on histone

            //calc 1st point on next histone
            G4double angle = (pi/2.0)+offset;
            G4ThreeVector pt( (Radius * sin(angle)),
                             (Radius * cos(angle)),
                             -zHist );
            G4RotationMatrix *nextRot = HistoneDetails[i+1].second;
            G4ThreeVector RotPt;
            RotPt[0] = pt.x()*nextRot->xx() + pt.y()*nextRot->yx() + pt.z()*nextRot->zx();
            RotPt[1] = pt.x()*nextRot->xy() + pt.y()*nextRot->yy() + pt.z()*nextRot->zy();
            RotPt[2] = pt.x()*nextRot->xz() + pt.y()*nextRot->yz() + pt.z()*nextRot->zz();
            G4ThreeVector end = RotPt + HistoneDetails[i+1].first; //1st pt on next histone

            //calc 2nd point on next histone
            angle= (1.0 * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
            G4ThreeVector secpt( (Radius * sin(angle)),
                                (Radius * cos(angle)),
                                (-zHist + (1.0/(G4double)nSteps)*2.0*(zHist)));
            RotPt[0] = secpt.x()*nextRot->xx() + secpt.y()*nextRot->yx() + secpt.z()*nextRot->zx();
            RotPt[1] = secpt.x()*nextRot->xy() + secpt.y()*nextRot->yy() + secpt.z()*nextRot->zy();
            RotPt[2] = secpt.x()*nextRot->xz() + secpt.y()*nextRot->yz() + secpt.z()*nextRot->zz();
            secpt = RotPt + HistoneDetails[i+1].first; //2nd pt on next histone

            G4ThreeVector StartDir = (start - path[path.size()-2]).unit(); //direction coming off histone
            G4ThreeVector MidPoint1 = ((TraceBP*0.3)*StartDir); //trace it forward x bp
            MidPoint1[0] *= nm;
            MidPoint1[1] *= nm;
            MidPoint1[2] *= nm;
            MidPoint1 += start;

            G4ThreeVector EndDir = (secpt-end).unit();
            G4ThreeVector MidPoint2 = ((TraceBP*0.3)*EndDir); //trace it forward x bp
            MidPoint2[0] *= nm;
            MidPoint2[1] *= nm;
            MidPoint2[2] *= nm;
            MidPoint2 = end - MidPoint2;

            //link the 4 points with Bezier in 100 steps
            Fibre_Bezier(start, MidPoint1, MidPoint2, end, path, 100);
        }
    }

}

//Bezier curve to smoothly join 2 points (bending through 2 other points)
void TsEllipsoidCellNAK::Fibre_Bezier(G4ThreeVector &start,
                     G4ThreeVector &MidPoint1,
                     G4ThreeVector &MidPoint2,
                     G4ThreeVector &end,
                     vector<G4ThreeVector> &path,
                     G4int nSteps)
{
    G4double j=0.0;
    for (G4int k=0;k<nSteps-1;k++){
        j+=1.0/(G4double)nSteps;
        G4ThreeVector BezPt = pow(1.0-j,3.0)*start + 3.0*pow(1.0-j,2.0)*j*MidPoint1 +
        3.0*(1.0-j)*j*j*MidPoint2 + j*j*j*end;
        path.push_back(BezPt);
    }
}

//split the DNA path into 0.34nm steps
void TsEllipsoidCellNAK::Fibre_SegmentDNAPath(vector<G4ThreeVector> &path)
{
    vector<G4ThreeVector> newPath;
    double rise=0.34*nm;
    int nPts=path.size();
    int counter=0;
    newPath.push_back(path[0]);
    for (int i=0;i<nPts-1;i++){
        G4ThreeVector vector = (path[i+1]-newPath[counter]);
        double length = vector.mag();
        int nDiv = length/rise;
        G4ThreeVector unit = vector.unit();
        unit *= rise;
        for (int j=0;j<nDiv;j++){
            G4ThreeVector nuwe = newPath[counter] + unit;
            newPath.push_back(nuwe);
            counter++;
        }
    }

    newPath=path;
    //path.clear();

}

//use DNA path to place sphere DNA volumes
void TsEllipsoidCellNAK::Fibre_PlaceDNASphere(vector<G4ThreeVector> &newPath)
{
    G4double helixRadius = 1.0*nm;
    G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    G4int nBP=newPath.size();
    G4double rBack=helixRadius - 0.24*nm;
    G4double rBase=rBack - 0.24*nm - 0.208*nm;

    for (int bp=0; bp<nBP-1; bp++){
        G4double angle1 = (G4double)bp * rotPair;
        G4double angle2 = angle1+pi;// + (120.0*pi/180.0); //offset for strand2 (major and minor groove)

        //temporary positions
        G4ThreeVector back1temp = G4ThreeVector((rBack*cos(angle1)), (rBack*sin(angle1)), 0.0);
        G4ThreeVector back2temp = G4ThreeVector((rBack*cos(angle2)), (rBack*sin(angle2)), 0.0);
        G4ThreeVector base1temp = G4ThreeVector((rBase*cos(angle1)), (rBase*sin(angle1)), 0.0);
        G4ThreeVector base2temp = G4ThreeVector((rBase*cos(angle2)), (rBase*sin(angle2)), 0.0);

        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,-1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm

        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);

        //Apply rotation
        G4RotationMatrix *rot = new G4RotationMatrix;
        rot->rotate(AngBetween, cross);

        Fibre_ApplyRotation(back1, back1temp, rot);
        Fibre_ApplyRotation(back2, back2temp, rot);
        Fibre_ApplyRotation(base1, base1temp, rot);
        Fibre_ApplyRotation(base2, base2temp, rot);

        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];

        G4int bpID = bp + 1;

        //SubCompartment names
        G4String back1_name = "back1";
        G4String back2_name = "back2";
        G4String base1_name = "base1";
        G4String base2_name = "base2";
        G4String water1_name = "water1";
        G4String water2_name = "water2";

        G4ThreeVector *posBack1 = &back1;
        CreatePhysicalVolume(back1_name, bpID, false, lBack1, 0, posBack1, aEnvelopeLog);
        G4ThreeVector *posBase1 = &base1;
        CreatePhysicalVolume(base1_name, bpID, false, lBase1, 0, posBase1, aEnvelopeLog);
        G4ThreeVector *posBack2 = &back2;
        CreatePhysicalVolume(back2_name, bpID, false, lBack2, 0, posBack2, aEnvelopeLog);
        G4ThreeVector *posBase2 = &base2;
        CreatePhysicalVolume(base2_name, bpID, false, lBase2, 0, posBase2, aEnvelopeLog);  

        if(fAddHydrationShell)
        {
        CreatePhysicalVolume(water1_name, bpID, false, lWater1, 0, posBack1, aEnvelopeLog);
        CreatePhysicalVolume(water2_name, bpID, false, lWater2, 0, posBack2, aEnvelopeLog);
         }

    }

    G4cout<<G4endl<<"DNA SPHERE BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
}

//rotate a 3vector by a rot matrix and return new coordinates
void TsEllipsoidCellNAK::Fibre_ApplyRotation(G4ThreeVector &rotated, G4ThreeVector &vector, G4RotationMatrix *rot)
{
    rotated[0] = vector[0]*rot->xx() + vector[1]*rot->yx() + vector[2]*rot->zx();
    rotated[1] = vector[0]*rot->xy() + vector[1]*rot->yy() + vector[2]*rot->zy();
    rotated[2] = vector[0]*rot->xz() + vector[1]*rot->yz() + vector[2]*rot->zz();

}

//Use DNA path to place halfCyl or quartCyl DNA volumes
void TsEllipsoidCellNAK::Fibre_PlaceDNA(vector<G4ThreeVector> &newPath)
{
    G4cout << "PlaceDNA start" << G4endl;
    G4double helixRadius = 1.0*nm;
    G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    G4int nBP=newPath.size();

    for (G4int bp=0; bp<nBP-1; bp++){

        //Position of base + back in xy
        //Definitely gives right handed coil (checked) -- left handed in -ve z?
        G4double angle1 = -(G4double)bp * rotPair;

        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm

        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);

        //Apply rotation
        G4RotationMatrix *rot1 = new G4RotationMatrix();
        rot1->rotate(AngBetween, cross);
        rot1->rotateZ(angle1);

        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];

        G4int bpID = bp + 1;

        //SubCompartment names
        G4String back1_name = "back1";
        G4String back2_name = "back2";
        G4String base1_name = "base1";
        G4String base2_name = "base2";
        G4String water1_name = "water1";
        G4String water2_name = "water2";

        G4ThreeVector *posBack1 = &back1;
        CreatePhysicalVolume(back1_name, bpID, false, lBack1, rot1, posBack1, aEnvelopeLog);
        G4ThreeVector *posBase1 = &base1;
        CreatePhysicalVolume(base1_name, bpID, false, lBase1, rot1, posBase1, aEnvelopeLog);
        G4ThreeVector *posBack2 = &back2;
        CreatePhysicalVolume(back2_name, bpID, false, lBack2, rot1, posBack2, aEnvelopeLog);
        G4ThreeVector *posBase2 = &base2;
        CreatePhysicalVolume(base2_name, bpID, false, lBase2, rot1, posBase2, aEnvelopeLog);


        
        if(fAddHydrationShell)
        {
        CreatePhysicalVolume(water1_name, bpID, false, lWater1, rot1, posBack1, aEnvelopeLog);
        CreatePhysicalVolume(water2_name, bpID, false, lWater2, rot1, posBack2, aEnvelopeLog);
        }
    }

    G4cout<<G4endl<<"DNA BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
}




