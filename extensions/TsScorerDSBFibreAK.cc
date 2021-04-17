// Scorer for DSBFibreAK
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


#include "TsScorerDSBFibreAK.hh"
#include "TsChromosome.hh"
#include "TsHitsRecord.hh"
#include "TsDefineDamage.hh"

#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Molecule.hh"
#include "Randomize.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <stdint.h>
using namespace std;


TsScorerDSBFibreAK::TsScorerDSBFibreAK(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM)
{
    // Initialize
    fTrackAveragedLET = 0 ;
    fTravelDistance =0 ;
    fEdep =0 ;
    fDoseInThisExposoure = 0;
    fExposoureID = 0;
    fCalZetaBeta=false;
    fZetaBeta_sq=0;

    //Define ntuple columns
    SetUnit("");
    // fNtuple->RegisterColumnI(&fEventID, "Event ID");
    // fNtuple->RegisterColumnD(&fPosX, "Position X", "um");
    // fNtuple->RegisterColumnD(&fPosY, "Position Y", "um");
    // fNtuple->RegisterColumnD(&fPosZ, "Position Z", "um");
    // fNtuple->RegisterColumnI(&fVoxelID, "fVoxelID ID");
    // fNtuple->RegisterColumnI(&fChromosomeID, "Chromosome ID");
    // fNtuple->RegisterColumnI(&fBasePairID, "BasePair number");
    // fNtuple->RegisterColumnS(&fVolumeName, "Volume name");
    fNtuple->RegisterColumnS(&fParticleName, "Particle name ");
    // fNtuple->RegisterColumnD(&fParticleEnergy, "Particle energy ", "eV");
    fNtuple->RegisterColumnD(&fEnergyDeposit, "Energy deposit", "eV");
    

    NumberOfHistoriesInRun =  fPm->GetIntegerParameter(GetFullParmName("NumberOfHistoriesInRun"));

    fProbabilityOfOHDamage = 0.4;
    if ( fPm->ParameterExists(GetFullParmName("ProbabilityForOHToCauseDamage")) ) {
        fProbabilityOfOHDamage = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityForOHToCauseDamage"));
    }

    fOHDamageWholeDNA = false;
    if ( fPm->ParameterExists(GetFullParmName("OHDamageWholeDNA")) ) {
        fOHDamageWholeDNA = fPm->GetBooleanParameter(GetFullParmName("OHDamageWholeDNA"));
    }

    fDamageThreshold = 17.5*eV;
    if ( fPm->ParameterExists(GetFullParmName("DamageThreshold")) ) {
        fDamageThreshold = fPm->GetDoubleParameter(GetFullParmName("DamageThreshold"),"Energy");
    }
    
    fUseLinearProbabilitythreshold = false;
    if ( fPm->ParameterExists(GetFullParmName("UseLinearProbabilitythreshold")) ) {
        fUseLinearProbabilitythreshold = fPm->GetBooleanParameter(GetFullParmName("UseLinearProbabilitythreshold"));
    }

    fLinearProbability_lower_limit= 5*eV;
        if ( fPm->ParameterExists(GetFullParmName("LinearProbability_lower_limit")) ) {
        fLinearProbability_lower_limit = fPm->GetDoubleParameter(GetFullParmName("LinearProbability_lower_limit"),"Energy");
    }
    fLinearProbability_upper_limit= 37.5*eV;
        if ( fPm->ParameterExists(GetFullParmName("LinearProbability_upper_limit")) ) {
        fLinearProbability_upper_limit = fPm->GetDoubleParameter(GetFullParmName("LinearProbability_upper_limit"),"Energy");
    }
    
    fDSBSeparation = 10; //10 BP
    if ( fPm->ParameterExists(GetFullParmName("DSBSeparation")) ) {
        fDSBSeparation = fPm->GetIntegerParameter(GetFullParmName("DSBSeparation"));
    }

    fDefineComplexity = true;
    if ( fPm->ParameterExists(GetFullParmName("DefineComlexity")) ) {
        fDefineComplexity = fPm->GetBooleanParameter(GetFullParmName("DefineComlexity"));
    }

    fComplexitySeparation = 10; //10 BP
    if ( fPm->ParameterExists(GetFullParmName("ComplexitySeparation")) ) {
        fComplexitySeparation = fPm->GetIntegerParameter(GetFullParmName("ComplexitySeparation"));
    }

    fExcludeShortFragment = true;
    if ( fPm->ParameterExists(GetFullParmName("ExcludeShortFragment")) ) {
        fExcludeShortFragment = fPm->GetBooleanParameter(GetFullParmName("ExcludeShortFragment"));
    }

    fLowerFragmentDetectionThreshold = 0; //In unit of basepair
    if ( fPm->ParameterExists(GetFullParmName("LowerFragmentDetectionThreshold")) ) {
         fLowerFragmentDetectionThreshold = fPm->GetIntegerParameter(GetFullParmName("LowerFragmentDetectionThreshold"));
    }

    fUpperFragmentDetectionThreshold = 1*3E8; //In unit of basepair, maximum chromosome size
    if ( fPm->ParameterExists(GetFullParmName("UpperFragmentDetectionThreshold")) ) {
         fUpperFragmentDetectionThreshold = fPm->GetIntegerParameter(GetFullParmName("UpperFragmentDetectionThreshold"));
    }

    fOnlyIncludeDSBinSDD = true;
    if ( fPm->ParameterExists(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD")) ) {
        fOnlyIncludeDSBinSDD = fPm->GetBooleanParameter(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD"));
    }

    fDosePerExposoure = 1 ;//In unit of Gy
    if ( fPm->ParameterExists(GetFullParmName("DosePerExposoure")) ) {
        fDosePerExposoure = fPm->GetUnitlessParameter(GetFullParmName("DosePerExposoure"));
    }

    fScoringRadius = 4.65*um ;
    if ( fPm->ParameterExists(GetFullParmName("ScoringRadius")) ) {
        fScoringRadius = fPm->GetDoubleParameter(GetFullParmName("ScoringRadius"),"Length");
    }

    fScoringTransX = 0*um;
    fScoringTransY = 0*um;
    fScoringTransZ = 0*um;

    if (fPm->ParameterExists(GetFullParmName("ScoringTransX"))) {
        
        G4String name1 = GetFullParmName("ScoringTransX");
        if (fPm -> ParameterExists(name1)){
            fScoringTransX = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("ScoringTransY");
        if (fPm -> ParameterExists(name1)){
            fScoringTransY = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("ScoringTransZ");
        if (fPm -> ParameterExists(name1)){
            fScoringTransZ = fPm->GetDoubleParameter(name1, "Length");
        }
    }

    fHistoneAsScavenger = true;
    if ( fPm->ParameterExists(GetFullParmName("HistoneAsScavenger")) ) {
        fHistoneAsScavenger = fPm->GetBooleanParameter(GetFullParmName("HistoneAsScavenger"));
    }


    GetGeoinfo();
    aDefineDamage = new TsDefineDamage();
    aDefineDamage->SetDamageThreshold(fDamageThreshold);
    aDefineDamage->SetUseLinearProbabilitythreshold(fUseLinearProbabilitythreshold);
    aDefineDamage->SetLinearProbability_lower_limit(fLinearProbability_lower_limit);
    aDefineDamage->SetLinearProbability_upper_limit(fLinearProbability_upper_limit);
    aDefineDamage->SetDSBSeparation(fDSBSeparation);
    aDefineDamage->SetDefineComplexity(fDefineComplexity);
    aDefineDamage->SetComplexitySeparation(fComplexitySeparation);
    aDefineDamage->SetExcludeShortFragment(fExcludeShortFragment);
    aDefineDamage->SetLowerFragmentDetectionThreshold( fLowerFragmentDetectionThreshold);
    aDefineDamage->SetUpperFragmentDetectionThreshold( fUpperFragmentDetectionThreshold);
    aDefineDamage->SetNucleusMass(fNucleusMass);
    aDefineDamage->SetTotalDNAContent(fTotalDNAContent);
    aDefineDamage->SetChromosomeDNAContent(fChromosomeDNAContent);
    aDefineDamage->SetOnlyIncludeDSBinSDD(fOnlyIncludeDSBinSDD);
    aDefineDamage->OutputSDDHeader("DNAdamage");

    // ********************************************************************************
    // Print parameters
    // ********************************************************************************
    G4cout<<"*********************************************************************************"<<G4endl;
    G4cout << "fProbabilityOfOHDamage = "<<fProbabilityOfOHDamage<<G4endl;
    G4cout << "fOHDamageWholeDNA = "<<fOHDamageWholeDNA<<G4endl;
    if (!fUseLinearProbabilitythreshold)
    G4cout << "fDamageThreshold = "<<fDamageThreshold/eV<<" eV"<<G4endl;
    else
    G4cout << "fDamageThreshold = "<<fLinearProbability_lower_limit/eV<<" - "<<fLinearProbability_upper_limit/eV<<" eV"<<G4endl;
    G4cout << "fDSBSeparation = "<<fDSBSeparation<<" bp"<<G4endl;
    G4cout << "fLowerFragmentDetectionThreshold = "<< fLowerFragmentDetectionThreshold<<G4endl;
    G4cout << "fUpperFragmentDetectionThreshold = "<< fUpperFragmentDetectionThreshold<<G4endl;
    G4cout<<"*********************************************************************************"<<G4endl;

    // ParticleOcurrances.push_back("gamma");

    hits_counter = 0;
}


TsScorerDSBFibreAK::~TsScorerDSBFibreAK()
{;}


G4bool TsScorerDSBFibreAK::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
  
    G4Track *           aTrack       = aStep->GetTrack();
    G4StepPoint *       theStepPoint = aStep->GetPreStepPoint();
    G4ThreeVector pos = theStepPoint->GetPosition();

    // check if within nucleus  
    G4bool WithInNucleus =
        ((pow((pos.x() - fScoringTransX),2)
        + pow((pos.y() - fScoringTransY),2)
        + pow((pos.z() - fScoringTransZ),2)) < pow(fScoringRadius,2));

    if(!WithInNucleus)
        return false;

    fEdep += aStep->GetTotalEnergyDeposit();
    fParticleName     = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();

    G4int fTrackID = aTrack->GetTrackID();
    G4int fParentID = aTrack->GetParentID();

    if(fParticleName == "proton")
        fTravelDistance += aStep->GetStepLength();

    // calculate (Zeff/beta)^2
    if (!fCalZetaBeta)
    {
        fParticleVelocity = aTrack->GetVelocity();
        G4double Lightspeed = 3*1e8*m/s;
        G4double a_beta= fParticleVelocity/Lightspeed;
        G4double z_eff = 1-exp(-125*a_beta);
        fZetaBeta_sq = pow(z_eff/a_beta, 2);
        G4cout<<"ParticleVelocity= "<<fParticleVelocity*s/m<<" m/sec;  "
              <<"beta = "<<a_beta<<"  "<<"(Z_eff/β)^2 ="<<fZetaBeta_sq<<G4endl;
        fCalZetaBeta = true;
    }

    G4String  VolumeName   = aTrack->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
    WithinDNAVolumme       = (VolumeName =="back1" || VolumeName =="back2" ||
                              VolumeName =="base1" || VolumeName =="base2" || 
                              VolumeName =="water1" || VolumeName =="water2"||
                              VolumeName =="Histone");

    G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

    WithinDirReactionVolumme = (VolumeName =="back1" || VolumeName =="back2" || VolumeName =="water1" || VolumeName =="water2");
    if(fOHDamageWholeDNA)
    WithinIndirReactionVolumme = (VolumeName =="back1" || VolumeName =="back2" || VolumeName =="base1" || VolumeName =="base2");
    else
    WithinIndirReactionVolumme = (VolumeName =="back1" || VolumeName =="back2");
 
    G4int VoxelID=0,  FibreID=0, bpIDinFiber=0; 
    G4int VoxelNumInSphere=0;
    G4int ChromosomeID=0;
    G4int bpIDinChromosome=0;

    if( WithinDNAVolumme) 
    { 
        G4VPhysicalVolume * MotherPhyVol = theStepPoint->GetTouchableHandle()->GetVolume(1);
        G4VPhysicalVolume * GrandMaPhyVol = theStepPoint->GetTouchableHandle()->GetVolume(2);  
        VoxelID = CalculateVoxelID ( pos);
        FibreID = MotherPhyVol->GetCopyNo();
        bpIDinFiber = aTrack->GetVolume()->GetCopyNo();
        FindPlaceInChromosome( VoxelID,  FibreID, bpIDinFiber, VoxelNumInSphere, ChromosomeID,  bpIDinChromosome );
        if(VoxelNumInSphere == -1 || VoxelID>=pow(fVoxel3Drepeat,3))
            return false;            

        fEventID = GetEventID()+ NumberOfHistoriesInRun*GetRunID();
        fEnergyDeposit = aStep->GetTotalEnergyDeposit();
        fParticleEnergy = aTrack->GetKineticEnergy();
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();
        fVolumeName = VolumeName; 
        fVoxelID    = VoxelNumInSphere;
        fBasePairID = bpIDinChromosome;
        fChromosomeID = ChromosomeID;
        if (fChromosomeID < 0) {
            return 1;
        }
    
        TsHitsRecord * hit = new TsHitsRecord();        
        hit->SetEventID(fEventID);
		hit->SetPosX(fPosX);
		hit->SetPosY(fPosY);
		hit->SetPosZ(fPosZ);
        hit->SetPosition( G4ThreeVector(fPosX, fPosY, fPosZ) );
        hit->SetBasePairID(fBasePairID);
        hit->SetChromosomeID(fChromosomeID);
        hit->SetVolumeName(fVolumeName); 
        hit->SetEdep(fEnergyDeposit);
        hit->SetParticleEnergy(fParticleEnergy);
        hit->SetParticleName(fParticleName);
        
       // **************************************** Direct damage ****************************************
        if ( 0 <= aTrack->GetTrackID() && 0 < aStep->GetTotalEnergyDeposit() && WithinDirReactionVolumme) 
        { 
            fNtuple->Fill();
            Hits.push_back(hit);
            return true;    
        } 

        // **************************************** Indirect damage ***************************************
        if (0 > aTrack->GetTrackID())
        {
            G4String SpeciesName = GetMolecule(aTrack)->GetName();
            G4bool IsSpeciesToKill = (SpeciesName=="OH^0" || SpeciesName=="e_aq^-1" || SpeciesName=="H^0");
            G4bool IsHydroxyl      = (SpeciesName=="OH^0");
            G4bool JustEnterVolume = (aStep->GetPreStepPoint()->GetStepStatus() != fGeomBoundary);

            // ************ Kill ALL spercies generated inside DNA volume ************
            if(aStep->IsFirstStepInVolume() ) 
            {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		        delete hit; 
                return false;
            }
            // ******************* Make the hydroxyl damage the DNA *******************
            else if (JustEnterVolume && WithinIndirReactionVolumme && IsHydroxyl)
            {
                hit->SetEdep(-0.001 * eV);
                hit->SetIsDirectDamage (false);

                G4bool reacted = true;
                if ( reacted && G4UniformRand() < fProbabilityOfOHDamage ) 
                {
                    Hits.push_back(hit);      
                }
                else delete hit;
                aStep->GetTrack()->SetTrackStatus(fStopAndKill); 
                return true;
            }
            // ************  Kill {OH•, e_aq- and H•} diffuse into DNA volume ************
            else if (!JustEnterVolume && IsSpeciesToKill)
            {
		        delete hit;
                if(!fHistoneAsScavenger&& VolumeName =="Histone")
                    return false;
                else{
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill); 
                    return false; 
                }
            }   
        }
        delete hit;
    } 
 
    return false;
}

void TsScorerDSBFibreAK::AccumulateEvent()
{
    G4int EventIDinAllRun = GetEventID()+NumberOfHistoriesInRun*GetRunID();
    HitsofEvents.push_back(Hits);
    EdepofEvents.push_back(fEdep);
    TravelDistanceofEvents.push_back(fTravelDistance);
    
    Hits.clear();
    fEdep=0;
    fTravelDistance=0;

}

void TsScorerDSBFibreAK::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) 
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsScorerDSBFibreAK* workerMTScorer = dynamic_cast<TsScorerDSBFibreAK*>(workerScorer);

    fZetaBeta_sq= workerMTScorer->fZetaBeta_sq;

    for(G4int i = 0; i < workerMTScorer->HitsofEvents.size(); i++)
        HitsofEvents.push_back(workerMTScorer->HitsofEvents[i]);
    workerMTScorer->HitsofEvents.clear();
    
    for (G4int i = 0; i < workerMTScorer->EdepofEvents.size(); i++)
        EdepofEvents.push_back(workerMTScorer->EdepofEvents[i]);
    workerMTScorer->EdepofEvents.clear();    

    for (G4int i = 0; i < workerMTScorer->TravelDistanceofEvents.size(); i++)
        TravelDistanceofEvents.push_back(workerMTScorer->TravelDistanceofEvents[i]);
    workerMTScorer->TravelDistanceofEvents.clear();
}

void TsScorerDSBFibreAK::Analyze( vector<TsHitsRecord*> Hits, G4int EventIDinAllRun)
{
    G4double accumulatedEdep = 0;
    G4double accumulatedTravelDistance =0;
    for(G4int i = 0; i <= EventIDinAllRun; i++){
        accumulatedEdep += EdepofEvents[i];
        accumulatedTravelDistance += TravelDistanceofEvents[i];
    }

    fDoseInThisExposoure = (1.6E-13*accumulatedEdep/MeV)/(fNucleusMass/1000);  //Gy
    fTrackAveragedLET = (accumulatedEdep/keV)/(accumulatedTravelDistance/um); //keV/um
    aDefineDamage->SetEdep(accumulatedEdep/MeV);
    aDefineDamage->SetLET (fTrackAveragedLET);
    aDefineDamage->SetZetaBeta_sq(fZetaBeta_sq);
    aDefineDamage->SetEventID(EventIDinAllRun);

    if( fDoseInThisExposoure >= fExposoureID*fDosePerExposoure)
    {
        fExposoureID ++;
        G4cout <<"Start new exposure("<< fDosePerExposoure<<" Gy per exposure)!"
               <<"ExposoureID = "<<fExposoureID<< " Event ID = "<<EventIDinAllRun<<G4endl;
    }


        

    // ***************************  seperate hits on different chromosome ***************************//
    vector<int> ChromosomeCopyNum;
    for (G4int i=0; i<Hits.size();i++) {
        ChromosomeCopyNum.push_back(Hits[i]->GetChromosomeID());
    }

    sort(ChromosomeCopyNum.begin(), ChromosomeCopyNum.end());
    ChromosomeCopyNum.erase( unique( ChromosomeCopyNum.begin(), ChromosomeCopyNum.end() ), ChromosomeCopyNum.end() );


    for(G4int i=0;i<ChromosomeCopyNum.size(); i++)
    {
        vector<TsHitsRecord*>  HitsOnThisChromosome;
        for(G4int j= 0; j<Hits.size(); j++)
        {
            if( Hits[j]->GetChromosomeID() == ChromosomeCopyNum[i])
                HitsOnThisChromosome.push_back(Hits[j]);
        }
        
        vector<TsHitsRecord*>   SBonStrand1,  SBonStrand2;
        vector<TsHitsRecord*>  SSBonStrand1, SSBonStrand2;
        vector <pair<TsHitsRecord*, TsHitsRecord*>> DSBpairs;
        aDefineDamage->SeparateHitsOnDifferentDNAStrands(HitsOnThisChromosome, SBonStrand1, true);
        aDefineDamage->SeparateHitsOnDifferentDNAStrands(HitsOnThisChromosome, SBonStrand2, false);
        aDefineDamage->DefineDSBorSSB(SBonStrand1,SBonStrand2, DSBpairs, SSBonStrand1, SSBonStrand2); 
        aDefineDamage->OutputDNAdamageTuple(SBonStrand1,SBonStrand2, "DNAdamage");
        aDefineDamage->OutputSDDFile(fOnlyIncludeDSBinSDD, DSBpairs, SSBonStrand1, SSBonStrand2, "DNAdamage", EventIDinAllRun, fExposoureID, ChromosomeCopyNum[i]);

        if(fExcludeShortFragment)
            aDefineDamage->ExcludeShortDNAFragments(DSBpairs, SSBonStrand1, SSBonStrand2, ChromosomeCopyNum[i]);
        

        HitsOnThisChromosome.clear();
        SBonStrand1.clear();
        SBonStrand2.clear();
        DSBpairs.clear();
        SSBonStrand1.clear();
        SSBonStrand2.clear();

    }
    aDefineDamage->OutputDNAdamageTupleHeader("DNAdamage");
    aDefineDamage->OutputDNAdamageSummary("DNAdamage_sum");
    ChromosomeCopyNum.clear();


}

void TsScorerDSBFibreAK::UserHookForEndOfRun()
{
    G4cout<<"\n\n-------------------------------------"<<G4endl;
    G4cout<<"HitsofEvents size:"<< HitsofEvents.size()<<G4endl;
    for(G4int id = 0; id < HitsofEvents.size(); id++)
        Analyze( HitsofEvents[id], id);  
 
}



void TsScorerDSBFibreAK::GetGeoinfo()
{
    string folder ="./supportFiles/";

    // ---------------------------------------- //
    //       Read geometry information file
    // ---------------------------------------- //
    string filename = folder+"Geoinfo.txt";
    ifstream infile;
    infile.open(filename);
    G4String input;

    if(!infile.is_open())
    {
        G4cout<<"ERROR: Unable to open file  "<<filename <<G4endl;
        exit(0);
    }
    else
        G4cout<<"Reading "<<filename<<G4endl;

	while(infile>>input)						
    {
		if(input == "DNAContentperFiber(bp):")  infile >> fFiberDNAContent;  //bp
        if(input == "DNAContentperVoxel(bp):")  infile >> fVoxelDNAContent;  //bp
		if(input == "DNAContentinTotal(bp):")  infile >> fTotalDNAContent;  //bp
		if(input == "NucluesMass(g):")  infile >> fNucleusMass;     //g
        if(input == "NumberofVoxels(subdomains):")  infile >> fNumofVoxel;
        if(input == "Voxel3DrepeatTimes:")  infile >> fVoxel3Drepeat;
        if(input == "VoxelSize(um):")  infile >> fVoxelSize;
        fVoxelSize = fVoxelSize*um;
        fVoxelContainerminXYZ = -fVoxelSize*fVoxel3Drepeat/2 ; // HalfSize
        G4cout << "Voxel size: " << fVoxelSize/um << G4endl;
       
	 }
    infile.close();


    // ---------------------------------------- //
    //           Get DNA content 
    // ---------------------------------------- //
    TsChromosome* HumanFibroblastNucleus = new TsChromosome();
    fChromosomeDNAContent = HumanFibroblastNucleus->GetDNAcontentofEachChromosome_BP(); //in unit of bp
    delete HumanFibroblastNucleus;


    // ---------------------------------------- //
    //           Read copy number file
    // ---------------------------------------- //
    filename = folder+"CopyNoTable.txt";
    infile.open(filename);

    if(!infile.is_open())    {
        G4cout<<"ERROR: Unable to open file  "<<filename <<G4endl;
        exit(0);
    }
    else
        G4cout<<"Reading "<<filename<<G4endl;

    G4int cpbox=0;
    G4int cpsphere=0;
    string line;
    while (getline(infile, line))
    {    
        stringstream stream(line.data());
        stream >>cpbox >>cpsphere ;  	
        fVoxelNumInBox.push_back(cpbox);   
        fVoxelNumInSphere.push_back(cpsphere);     
    }
    infile.close();

    
    // ---------------------------------------- //
    //   Read chromosomeID and voxel ID 
    // ---------------------------------------- //
    filename = folder+"signedCHVoxel.txt";
    infile.open(filename);

    if(!infile.is_open())    {
        G4cout<<"ERROR: Unable to open file  "<<filename <<G4endl;
        exit(0);
    }
    else
        G4cout<<"Reading "<<filename<<G4endl;

    G4int ID_X=0, ID_Y=0, ID_Z=0, CH_ID=0, Voxel_ID=0;
    while (getline(infile, line))
    {    
        stringstream stream(line.data());
        stream >>cpsphere >>ID_X >>ID_Y >>ID_Z >>CH_ID >>Voxel_ID;  	
        fCH_ID.push_back(CH_ID);   
        fVoxel_ID.push_back(Voxel_ID);     
    }
    infile.close();

}

void TsScorerDSBFibreAK::FindPlaceInChromosome(G4int VoxelID, G4int FibreID,G4int bpIDinFiber,G4int &VoxelNumInSphere, G4int &ChromosomeID, G4int &bpIDinChromosome )
{
    G4int VoxelIDinNucleus=fVoxelNumInSphere[VoxelID];
    G4int VoxelIDonCH= fVoxel_ID[VoxelIDinNucleus];

    VoxelNumInSphere = VoxelIDinNucleus;
    ChromosomeID     = fCH_ID[VoxelIDinNucleus];
    bpIDinChromosome = VoxelIDonCH*fVoxelDNAContent+ FibreID*fFiberDNAContent + bpIDinFiber; // notice voxelID start from 0
}

G4int TsScorerDSBFibreAK::CalculateVoxelID (G4ThreeVector position)
{
    G4int CalculatedVoxelID = 0;
    G4int nx =0, ny=0, nz=0;
    nx = floor ((position.x()-fVoxelContainerminXYZ)/fVoxelSize);
    ny = floor ((position.y()-fVoxelContainerminXYZ)/fVoxelSize);
    nz = floor ((position.z()-fVoxelContainerminXYZ)/fVoxelSize);

    CalculatedVoxelID = nz*fVoxel3Drepeat*fVoxel3Drepeat+ny*fVoxel3Drepeat+nx;

    return CalculatedVoxelID;

}
