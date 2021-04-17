// Extra Class for DSBFibre
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

#include <stdint.h>
#include "TsDefineDamage.hh"
#include "TsHitsRecord.hh"
#include "TsChromosome.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>


using namespace CLHEP;
using namespace std;



TsDefineDamage::TsDefineDamage()
{
    fDamageSiteSize =10; // 10 basepairs
    TsChromosome* HumanFibroblastNucleus = new TsChromosome();
    fChromosomeDNAContent = HumanFibroblastNucleus->GetDNAcontentofEachChromosome_BP(); //in unit of bp
    fTotalDNAContent      = HumanFibroblastNucleus->GetTotalDNAContent_MBP()*1E6; //in unit of bp
    delete HumanFibroblastNucleus;


    fEdep =0;
    fLET  =0;
    fZetaBeta_sq=0;
    fEventID =0;
    fNucleusMass =0;
    fZetaBeta_sq =0;


    numSB =0 ;
    numSB_dir =0 ;
    numSB_indir =0 ;
    numSSB =0 ;
    numSSB_dir =0 ;
    numSSB_indir =0 ;
    numDSB =0 ;
    numDSB_dir =0 ;
    numDSB_indir =0 ;
    numDSB_hybrid=0;
    numSSB_P =0 ;
    numDSB_P =0 ;
    numDSB_PP =0 ;

    Excluded_numSB =0 ;
    Excluded_numSB_dir =0 ;
    Excluded_numSB_indir =0 ;
    Excluded_numSSB =0 ;
    Excluded_numSSB_dir =0 ;
    Excluded_numSSB_indir =0 ;
    Excluded_numDSB =0 ;
    Excluded_numDSB_dir =0 ;
    Excluded_numDSB_indir =0 ;
    Excluded_numDSB_hybrid=0;
    Excluded_numSSB_P =0 ;
    Excluded_numDSB_P =0 ;
    Excluded_numDSB_PP =0 ;

    lastOutputEventID = -1;
    lastOutputExposureID = -1;
    lastDSBonlyOutputEventID = -1;
    lastDSBonlyOutputExposureID = -1;
}
    
TsDefineDamage::~TsDefineDamage()
{}


void TsDefineDamage::SeparateHitsOnDifferentDNAStrands(vector<TsHitsRecord*> Hits, vector<TsHitsRecord*> &HitsBack, G4double IsStrand1)
{
    // Merge different hits on the same backbone (and neighboring hydration shell)
    // Accumulate energy deposition and only record one effective hit for each damaged backbone

    if (Hits.size()<=0)
		return;

    G4String backName = "back1";
    G4String hydrationName = "water1";
	if (!IsStrand1){
        backName = "back2";
        hydrationName = "water2";
    }

    vector<TsHitsRecord*> HitsOnThisStrand;
    vector<int> DamagedBasepairID;
    for (unsigned index=0; index<Hits.size();index++)  
    {
        if (Hits[index]->GetVolumeName() == backName || Hits[index]->GetVolumeName() == hydrationName)
        {
            HitsOnThisStrand.push_back(Hits[index]);
            DamagedBasepairID.push_back(Hits[index]->GetBasePairID());
        }	
    }
    sort(DamagedBasepairID.begin(), DamagedBasepairID.end());
    DamagedBasepairID.erase( unique( DamagedBasepairID.begin(), DamagedBasepairID.end() ), DamagedBasepairID.end() );
    
    double TotalEdepInOneCopy = 0;
    int index =0;
    for(unsigned iter=0; iter<DamagedBasepairID.size(); iter++)
    {
        for(unsigned i =0; i<HitsOnThisStrand.size(); i++)
        {
            if(HitsOnThisStrand[i]->GetBasePairID()==DamagedBasepairID[iter]){
                TotalEdepInOneCopy +=HitsOnThisStrand[i]->GetEdep();
                index = i;
            }
        }

        if(CasueDirectDamge(TotalEdepInOneCopy) || TotalEdepInOneCopy<0)
        {
            TsHitsRecord* EffectiveDamage = HitsOnThisStrand[index];
            EffectiveDamage->SetEdep(TotalEdepInOneCopy);
            HitsBack.push_back(EffectiveDamage);
        }                                       
        TotalEdepInOneCopy = 0;
    }

    HitsOnThisStrand.clear();
    DamagedBasepairID.clear();
	
}

void TsDefineDamage::DefineDSBType (TsHitsRecord * &hit1, TsHitsRecord * &hit2)
{
    G4bool IsDir_1 = hit1->GetIsDirectDamage();
    G4bool IsDir_2 = hit2->GetIsDirectDamage();
    if( IsDir_1 && IsDir_2) 
    {
        hit1->SetDSBType(0);
        hit2->SetDSBType(0);
    }
        
    else if (!IsDir_1 && !IsDir_2) 
    {
        hit1->SetDSBType(1);
        hit2->SetDSBType(1);
    }
        
    else 
    {
        hit1->SetDSBType(2);
        hit2->SetDSBType(2);
    }
        
}


G4bool TsDefineDamage::CasueDirectDamge(G4double EnergyDeposition)
{
    if(!fUseLinearProbabilitythreshold )
    {
        if(EnergyDeposition>fDamageThreshold)
            return true;
    }
    else
    {
        if (EnergyDeposition>=fLinearProbability_upper_limit)
            return true;
        else if (EnergyDeposition<=fLinearProbability_lower_limit)
            return false;
        else 
        {
            G4double samplePro = G4UniformRand();
            G4double acceptPro = (EnergyDeposition-fLinearProbability_lower_limit)/(fLinearProbability_upper_limit-fLinearProbability_lower_limit);
            if (samplePro < acceptPro)
                return true;
        }
    }
    return false;
}


void TsDefineDamage::GetSSBonStrand(vector<TsHitsRecord*> Hits, vector<TsHitsRecord*> & SSBonStrand)
{
    for (int iter=0; iter<Hits.size(); iter++) // find SSB
    {
        G4bool   isDSB = Hits[iter]->GetFlagDSB();
        if( !isDSB )
        {
            Hits[iter]->SetFlagSSB(true);
            SSBonStrand.push_back(Hits[iter]);
        }
    }
}

void AddCountofDifferentTypeSSB(vector<TsHitsRecord*> SSBonStrand, G4int &SSB_num, G4int &dir_SSB_num, G4int &indir_SSB_num)
{
    SSB_num += SSBonStrand.size();
    for(G4int i=0; i<SSBonStrand.size(); i++)
    {
        if(SSBonStrand[i]->GetIsDirectDamage())
            dir_SSB_num ++;
        else
            indir_SSB_num++;  
    }
}

void AddCountofDifferentTypeDSB (vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, int &NumTotal, int &NumDir, int &NumIndir, int &NumHybird)
{
    NumTotal += DSB_pairs.size();
    for(G4int i=0; i<DSB_pairs.size(); i++ )
    {
        G4bool IsDir_1 = DSB_pairs[i].first->GetIsDirectDamage();
        G4bool IsDir_2 = DSB_pairs[i].second->GetIsDirectDamage();
        if( IsDir_1 && IsDir_2) 
            NumDir++;    
        else if (!IsDir_1 && !IsDir_2) 
            NumIndir++;    
        else 
            NumHybird++;     
    }
}


void TsDefineDamage::DefineDSBorSSB(vector<TsHitsRecord*> HitsBack1, vector<TsHitsRecord*> HitsBack2,
                                    vector <pair<TsHitsRecord*, TsHitsRecord*>> &DSB_pairs,
                                    vector<TsHitsRecord*> & SSBonStrand1,vector<TsHitsRecord*> & SSBonStrand2)
{
    
	if (HitsBack1.size()<=0 && HitsBack2.size()<=0)
		return;
	else
	{
		int cp1=0;
		int cp2=0;
		bool isDSB=false;
		bool isSSB=false;
		bool isBreakBack1 = false;
        bool isBreakBack2 = false;

		//---------------------------- Find DSB ----------------------------
         
		for (int iter=0; iter<HitsBack1.size(); iter++) // find DSB
		{
			cp1=HitsBack1[iter]->GetBasePairID();
            for(int loop=0; loop<HitsBack2.size(); loop++)
            {
                cp2=HitsBack2[loop]->GetBasePairID();
                if( abs(cp1-cp2)<=fDSBSeparation)
                {
                    bool DSBflag1 = HitsBack1[iter]->GetFlagDSB();
                    bool DSBflag2 = HitsBack2[loop]->GetFlagDSB();
                    if( DSBflag1==false && DSBflag2==false) //haven't been classified as DSB with former hits
                    {
                        HitsBack1[iter]->SetFlagDSB(true);
                        HitsBack2[loop]->SetFlagDSB(true);

                        DefineDSBType (HitsBack1[iter], HitsBack2[loop]);
                        pair<TsHitsRecord*, TsHitsRecord*> OnePairDSB;
                        OnePairDSB.first=HitsBack1[iter];
                        OnePairDSB.second=HitsBack2[loop];
                        DSB_pairs.push_back(OnePairDSB);                         
                    }
                }
            }
		}

		// ---------------------------- Find SSB ----------------------------
        GetSSBonStrand(HitsBack1, SSBonStrand1);
        GetSSBonStrand(HitsBack2, SSBonStrand2);


        // ---------------------------- calculate damage number ----------------------------
        AddCountofDifferentTypeDSB(DSB_pairs,    numDSB, numDSB_dir, numDSB_indir, numDSB_hybrid );
        AddCountofDifferentTypeSSB(SSBonStrand1, numSSB, numSSB_dir, numSSB_indir);
        AddCountofDifferentTypeSSB(SSBonStrand2, numSSB, numSSB_dir, numSSB_indir);
        numSB       = numSSB       + numDSB*2;
        numSB_dir   = numSSB_dir   + numDSB_dir*2   + numDSB_hybrid;
        numSB_indir = numSSB_indir + numDSB_indir*2 + numDSB_hybrid;

	}
}



void TsDefineDamage::ExcludeShortDNAFragments(vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, vector<TsHitsRecord*> SSBonStrand1,vector<TsHitsRecord*> SSBonStrand2, G4int ChromosomeID)
{
    // split chromosomes according to DSB pair position
    TsChromosome * processChromosome = new TsChromosome();
    vector<vector<G4int> > SplitedChromosome = processChromosome->SplitChromosome(DSB_pairs);
    delete processChromosome;

    G4int InvisableFragmentNum =0;

    vector<TsHitsRecord*> SSBDamages;
    SSBDamages.insert(SSBDamages.end(),SSBonStrand1.begin(),SSBonStrand1.end());
    SSBDamages.insert(SSBDamages.end(),SSBonStrand2.begin(),SSBonStrand2.end());
    
    vector<G4int> AllFreeEndofThisChromosome = SplitedChromosome[ChromosomeID-1];
    for (G4int i=1; i<AllFreeEndofThisChromosome.size(); i++ )
    {
        G4int  FrangmentSize = AllFreeEndofThisChromosome[i] - AllFreeEndofThisChromosome[i-1];
        G4bool IsInvisableFragment = (FrangmentSize<fLowerFragmentDetectionThreshold || FrangmentSize>fUpperFragmentDetectionThreshold );
        if (IsInvisableFragment)
            InvisableFragmentNum++;

        // -- exclude SSBs on the short fragments
    vector<TsHitsRecord*>::iterator it;
        for(it=SSBDamages.begin();it!=SSBDamages.end();)
        {
            G4int BasepairID   = (*it)->GetBasePairID();
            G4bool SSBOnThisFragment = (BasepairID>=AllFreeEndofThisChromosome[i-1] &&  BasepairID<= AllFreeEndofThisChromosome[i]);
            if(IsInvisableFragment && SSBOnThisFragment )
                it=SSBDamages.erase(it);
            else
                ++it;
        }

        // -- exclude DSBs on the short fragments
        vector <pair<TsHitsRecord*, TsHitsRecord*>> ::iterator iter;
        for(iter=DSB_pairs.begin();iter!=DSB_pairs.end();)
        {
            G4int BasepairID   = (*iter).first->GetBasePairID();
            G4bool DSBOnThisFragment = (BasepairID == AllFreeEndofThisChromosome[i]);
            if(IsInvisableFragment && DSBOnThisFragment )
                iter=DSB_pairs.erase(iter);
            else
                ++iter;
        }
    }

    AddCountofDifferentTypeDSB(DSB_pairs,  Excluded_numDSB, Excluded_numDSB_dir, Excluded_numDSB_indir, Excluded_numDSB_hybrid );
    AddCountofDifferentTypeSSB(SSBDamages, Excluded_numSSB, Excluded_numSSB_dir, Excluded_numSSB_indir);
    Excluded_numSB       = Excluded_numSSB       + Excluded_numDSB*2;
    Excluded_numSB_dir   = Excluded_numSSB_dir   + Excluded_numDSB_dir*2   + Excluded_numDSB_hybrid;
    Excluded_numSB_indir = Excluded_numSSB_indir + Excluded_numDSB_indir*2 + Excluded_numDSB_hybrid;

    SplitedChromosome.clear();
    AllFreeEndofThisChromosome.clear();
    SSBDamages.clear();

}



void TsDefineDamage::OutputDNAdamageTupleHeader(string filename)
{

	ofstream outfile;
    outfile.open(filename+".header",std::ios::out);

    outfile << "Columns of data are as follows:\n"
            << "1: Event ID" <<"\n"
            << "2: Chromosom ID " <<"\n"
            << "3: bp ID on chromosome "<<"\n"
            << "4: Volume name" <<"\n"
            << "5: Particel name" <<"\n"
            << "6: Particel Energy [eV]" <<"\n"
            << "7: Energy deposit [eV]" <<"\n"
            << "8: FlagSSB    (0-Not SSB, 1-Is SSB)" <<"\n"
            << "9: FlagDSB    (0-Not DSB, 1-Is DSB)" <<"\n"
            << "10: DSBType   (-1-Not DSB, 0-dir DSB, 1-indir DSB, 2-hybrid DSB)" <<"\n";
	outfile.close();      

}

void TsDefineDamage::OutputDNAdamageSummary(string filename)
{
    G4double NucleusDose = (fEdep*1.6E-13)/(fNucleusMass/1000); // MeV/g —> Gy
	ofstream outfile;
    outfile.open(filename+".header",std::ios::out);

    outfile <<"Recoreded Energy deposition(MeV): "<<fEdep<<"\n"
            <<"(Zeff/β)^2 :                      "<<fZetaBeta_sq<<"\n"
            <<"LET(keV/um):                      "<<fLET<<"\n"
            <<"Dose in nucleusmass(Gy):          "<<NucleusDose<<"\n"
            <<"DNA content in nucleus (Gbp):     "<<fTotalDNAContent/1E9<<"\n"
            <<"Simulated events:                 "<<fEventID+1<<"\n\n";
    outfile <<"************************************************"<<"\n"
            <<"************** Yied of DNA damage **************" <<"\n"
            <<"Total SB:           "<<  numSB       <<"\n"
            <<"Total SB_direct:    "<<  numSB_dir   <<"\n"
            <<"Total SB_indirect:  "<<  numSB_indir <<"\n"
            <<"Total SSB:          "<<  numSSB      <<"\n"
            <<"Total SSB_direct:   "<<  numSSB_dir  <<"\n"
            <<"Total SSB_indirect: "<<  numSSB_indir<<"\n"
            <<"Total DSB:          "<<  numDSB    <<"\n"
            <<"Total DSB_direct:   "<<  numDSB_dir    <<"\n"
            <<"Total DSB_indirect: "<<  numDSB_indir  <<"\n"
            <<"Total DSB_hybrid:   "<<  numDSB_hybrid <<"\n"
            <<"Total SSB+:  "<<  numSSB_P  <<"\n"
            <<"Total DSB+:  "<<  numDSB_P  <<"\n"
            <<"Total DSB++: "<<  numDSB_PP <<"\n\n"
            <<"SSB/Gy/Gbp : "<<   numSSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"DSB/Gy/Gbp : "<<   numDSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"************************************************" <<"\n\n";
    if(fExcludeShortFragment) 
    outfile <<"************************************************"<<"\n"
            <<"**  Yied of DNA damage (ExcludeShortFragment) **" <<"\n"
            <<"Total SB:           "<<  Excluded_numSB       <<"\n"
            <<"Total SB_direct:    "<<  Excluded_numSB_dir   <<"\n"
            <<"Total SB_indirect:  "<<  Excluded_numSB_indir <<"\n"
            <<"Total SSB:          "<<  Excluded_numSSB  <<"\n"
            <<"Total SSB_direct:   "<<  Excluded_numSSB_dir  <<"\n"
            <<"Total SSB_indirect: "<<  Excluded_numSSB_indir  <<"\n"
            <<"Total DSB:          "<<  Excluded_numDSB    <<"\n"
            <<"Total DSB_direct:   "<<  Excluded_numDSB_dir    <<"\n"
            <<"Total DSB_indirect: "<<  Excluded_numDSB_indir  <<"\n"
            <<"Total DSB_hybrid:   "<<  Excluded_numDSB_hybrid <<"\n"
            <<"Total SSB+:  "<<   Excluded_numSSB_P  <<"\n"
            <<"Total DSB+:  "<<   Excluded_numDSB_P  <<"\n"
            <<"Total DSB++: "<<   Excluded_numDSB_PP <<"\n"
            <<"SSB/Gy/Gbp : "<<   Excluded_numSSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"DSB/Gy/Gbp : "<<   Excluded_numDSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"************************************************" <<"\n";
	outfile.close();

}


void TsDefineDamage::OutputDNAdamageTuple(vector<TsHitsRecord*> HitsBack1, vector<TsHitsRecord*> HitsBack2, string filename)
{
    vector<TsHitsRecord*> Hit;
	ofstream outfile;
    outfile.open(filename+".phsp",std::ios::out|std::ios::app);
    

	for(int backind=1; backind<=2; backind++)
	{
		if (backind==1) 
			Hit=HitsBack1;
		else
			Hit=HitsBack2;

		for(unsigned i=0; i<Hit.size(); i++) 
		{   
			if(Hit[i]->GetFlagSSB()==true || Hit[i]->GetFlagDSB()==true)
            {
                outfile << setiosflags(ios::fixed)<<setprecision(0)<<std::setw(3) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetEventID() <<"    "
                        << setiosflags(ios::fixed)<<setprecision(0)<<std::setw(3) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetChromosomeID()<<"    "
                        << setiosflags(ios::fixed)<<setprecision(0)<<std::setw(12)<<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetBasePairID() <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(6) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetVolumeName() <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(12)<<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetParticleName() <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(18)<<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetParticleEnergy()/eV  <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(12)<<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetEdep()/eV <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagSSB() <<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagDSB()<<"    "
                        << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetDSBType()<<"    ";
                
                outfile << "\n";
            }
		}
	}
	outfile.close(); 
    Hit.clear();     

}


int GetStranID(string name)
{
    if (name == "back1" || name == "water1")
        return 1;
    return 4;
}


void TsDefineDamage::SeparateDamageSitesOnOneStrand(vector <TsHitsRecord*> Hits, vector <G4ThreeVector> &XYZPosition, vector <double  > &ChromosomePosition, 
                                              vector <double  > &Cause, vector <int > &DamageSiteDSBflag, vector< vector <vector <int > > >  &DamageSites)
{
    if (Hits.size()==0)
        return;

    int StartIDofDamageSite = Hits[0]->GetBasePairID();
    int EndIDofDamageSite   = StartIDofDamageSite+fDamageSiteSize-1;
    int DamageType     = 0;
    int localBasePairID = 0;
    int BasePairID     = 0;
    bool DirectDamageFlag   = false;
    bool IndirectDamageFlag = false;
    G4ThreeVector aXYZPos ;

    int iter=0;
    while( iter<Hits.size()  )
    {
        
        vector <vector <int > > aDamageSite;            
        for(int i=iter; i<Hits.size(); i++)
        {  
            BasePairID = Hits[i]->GetBasePairID();
            if( BasePairID >= StartIDofDamageSite && BasePairID<= EndIDofDamageSite )
            {
                if(Hits[i]->GetIsDirectDamage() ) {DamageType = 1; DirectDamageFlag =true;}
                if(!Hits[i]->GetIsDirectDamage()) {DamageType = 2; IndirectDamageFlag =true;}
                localBasePairID =BasePairID-StartIDofDamageSite +1;
                vector <int > aDamage = {GetStranID(Hits[i]->GetVolumeName()), localBasePairID,  DamageType };
                aDamageSite.push_back(aDamage);
                aXYZPos    = Hits[i]->GetPosition();
                iter++;  
            }
            else 
                break;          
        } 

        int ChromoID = Hits[0]->GetChromosomeID();
        int ChromoLength = fChromosomeDNAContent[ChromoID-1];
        double ChromoPos = (double)BasePairID/ChromoLength;

        bool aCause =0;
             if (DirectDamageFlag == true && IndirectDamageFlag == false)  aCause =0;
        else if (DirectDamageFlag == false && IndirectDamageFlag == true)  aCause =1;
        else if (DirectDamageFlag == true  && IndirectDamageFlag == true)  aCause =2;


        XYZPosition.push_back(aXYZPos);            // Field 2
        ChromosomePosition.push_back(ChromoPos);   // Field 4
        Cause.push_back(aCause);                   // Field 5
        DamageSiteDSBflag.push_back(0);            // Field 6, DSB flag
        DamageSites.push_back(aDamageSite);        // Field 7

        if (iter == Hits.size())  
            break;     

        StartIDofDamageSite = Hits[iter]->GetBasePairID();
        EndIDofDamageSite   = StartIDofDamageSite + fDamageSiteSize-1;
        G4cout<<"iter="<<iter<<", DamageSites.size()="<<DamageSites.size()<<G4endl;
    }
}


void TsDefineDamage::OutputSDDFile(bool OnlyIncludeDSBinSDD, vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, 
                            vector<TsHitsRecord*> SSBonStrand1,vector<TsHitsRecord*> SSBonStrand2,
                            string filename, G4int EventID, G4int ExposureID, G4int ChromosomeID)
{
    // Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al. A New Standard DNA Damage (SDD) Data Format[J]. Radiation research, 2018.
    // Data was generated with 7 SDD filed.
    // Filed 1: Calssification (2 Int) : Is damage associated with new primary particle or new exposure? And event ID.
    // Filed 2: X,Y,Z (3x3 Floats): Spatial X, Y, Z coordinates and extent (unit: um)
    // Filed 3: Chromosome IDs (4 Int) : ID of chromosome/chromatid where damage occurred and on which arm (long/short) or specification of non-nuclear DNA type.
    // Filed 4: Chromosome position (Float): Location of damage within chromosome
    // Filed 5: Cause (3 Int): Cause of damage - direct or indirect and number
    // Filed 6: Damage type (3 int): Types of damage at site (Base damage, SSB, DSB)
    // Filed 7: Full break spec (.../3 int/...) Full description of strand break structure

    vector <G4ThreeVector> XYZPosition;             // For field 2
    vector <double  > ChromosomePosition;           // For filed 4
    vector <double  > Cause;                        // For filed 5
    vector <int > DamageSiteDSBflag;                // For field 6
    vector< vector <vector <int > > > DamageSites;  // For field 7 
    
    vector<TsHitsRecord*> SSBDamages;
    SSBDamages.insert(SSBDamages.end(),SSBonStrand1.begin(),SSBonStrand1.end());
    SSBDamages.insert(SSBDamages.end(),SSBonStrand2.begin(),SSBonStrand2.end());

    for( int i=0; i<DSB_pairs.size(); i++)
    {
        vector <vector <int > > aDamageSite;
        pair<TsHitsRecord*, TsHitsRecord*> onePairDSB = DSB_pairs[i]; 
        int DamageSiteStartID = min(onePairDSB.first->GetBasePairID(), onePairDSB.second->GetBasePairID());
        int DamageSiteEndID   = DamageSiteStartID + fDamageSiteSize-1;
        int BasePairID = 0;
        int DamageType = 0;
        bool DirectDamageFlag   = false;
        bool IndirectDamageFlag = false;

        // ------ Get all SSB/DSB damages on this site 
        vector<TsHitsRecord*> HitsOnThisSite ;
        HitsOnThisSite.push_back(onePairDSB.first);
        HitsOnThisSite.push_back(onePairDSB.second);

        vector<TsHitsRecord*>::iterator it;
        for(it=SSBDamages.begin();it!=SSBDamages.end();)
        {
            BasePairID = (*it)->GetBasePairID();
            if(BasePairID>= DamageSiteStartID && BasePairID<=DamageSiteEndID )
            {
                HitsOnThisSite.push_back(*it); 
                it=SSBDamages.erase(it);
            }
            else
                ++it;
        }

        // ------ Get damage spectrum in this site 
        for (int i=0; i<HitsOnThisSite.size(); i++)
        {
            if ( HitsOnThisSite[i]->GetIsDirectDamage() ) 
                {DamageType = 1; DirectDamageFlag =true;}
            else
                {DamageType = 2; IndirectDamageFlag =true;}

            int localBasePairID = HitsOnThisSite[i]->GetBasePairID() - DamageSiteStartID +1;
            vector <int > OneHit = {GetStranID(HitsOnThisSite[i]->GetVolumeName()), localBasePairID,  DamageType };

            aDamageSite.push_back(OneHit);
        }

        // ------ sort damage spectrum, make sure damages are recorded from stran1 to stand 2
        sort(aDamageSite.begin(), aDamageSite.end());

        int ChromoLength      = fChromosomeDNAContent[ChromosomeID-1];
        double ChromoPos      = (double)DamageSiteStartID/ChromoLength;
        G4ThreeVector aXYZPos = HitsOnThisSite[0]->GetPosition();

        int aCause =0;
             if (DirectDamageFlag == true  && IndirectDamageFlag == false)  aCause=0;
        else if (DirectDamageFlag == false && IndirectDamageFlag == true )  aCause=1;
        else if (DirectDamageFlag == true  && IndirectDamageFlag == true )  aCause=2;


        XYZPosition.push_back(aXYZPos);            // Field 2
        ChromosomePosition.push_back(ChromoPos);   // Field 4
        Cause.push_back(aCause);                   // Field 5
        DamageSiteDSBflag.push_back(1);            // Field 6, DSB flag
        DamageSites.push_back(aDamageSite);        // Field 7
    }
    
    if(!OnlyIncludeDSBinSDD )
        SeparateDamageSitesOnOneStrand(SSBDamages, XYZPosition, ChromosomePosition, Cause, DamageSiteDSBflag, DamageSites);
    
    // ------ output 
    ofstream outfile;
    outfile.open(filename+"_sdd.txt",std::ios::out|std::ios::app);

    for(int i=0; i<DamageSites.size(); i++)
    {
        vector <vector <int > >  aDamageSite = DamageSites[i];

        G4int newExposoureFlag =0;
        if(EventID != lastOutputEventID)
        {
            newExposoureFlag = 1;
            lastOutputEventID = EventID;
        }
        if (ExposureID !=  lastOutputExposureID)
        {
            newExposoureFlag = 2;
            lastOutputExposureID = ExposureID;
        }

        outfile << newExposoureFlag<<","<< EventID<<"; ";                                                 // Filed 1 
        outfile << XYZPosition[i].x()/um<<","<<XYZPosition[i].y()/um<<","<<XYZPosition[i].z()/um<<"; ";   // Filed 2
        outfile << "1," << ChromosomeID-1<<","  << "1,"<< "0" <<"; ";                                     // Filed 3, strat from 0
        outfile << ChromosomePosition[i]<<"; ";                                                           // Filed 4
        outfile << Cause[i]<<"; ";                                                                        // Filed 5
        outfile <<"0,"<<aDamageSite.size()<<","<<DamageSiteDSBflag[i]<<"; ";                              // Filed 6, num of base damage, num of strand break, DSB flag
          
        for(int j=0; j<aDamageSite.size(); j++)
            outfile <<aDamageSite[j][0]<<","<<aDamageSite[j][1]<<","<<aDamageSite[j][2]<<"/";             //  Filed 7: Full break spec
        
        outfile << "; \n";
    }
    outfile.close();

    XYZPosition.clear();
    ChromosomePosition.clear();
    Cause.clear();
    DamageSiteDSBflag.clear();
    DamageSites.clear();

}



void TsDefineDamage::OutputSDDHeader(string filename)
{
    // Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al. A New Standard DNA Damage (SDD) Data Format[J]. Radiation research, 2018.

    // The header file need to be modifed later

    ofstream outfile;
    outfile.open(filename+"_sdd.txt");
    outfile <<"SDD Version, SDDv1.0; "
            <<"\n"<<"Software,TOPAS-nBio;"
            <<"\n"<<"Author, Hongyu Zhu;"
            <<"\n"<<"Simulation Details, NULL;"
            <<"\n"<<"Source,   Monoenergetic plane-parallel proton beam uniformly exposing nucleus. Energy: 1.0 MeV;"
            <<"\n"<<"Source type,  1;"
            <<"\n"<<"Incident particles, "<<fEventID+1<<";"
            <<"\n"<<"Mean particle energy,  1.0;  # This is an illustrative comment #"
            <<"\n"<<"Energy Distribution,  M, 0;"
            <<"\n"<<"Particle fraction, 1.0; "
            <<"\n"<<"Dose or fluence, 1, 22; "
            <<"\n"<<"Dose rate, 0.0;"
            <<"\n"<<"Irradiation target, Simple spherical cell model, radius 4.65 um;"
            <<"\n"<<"Volumes, 0,5,5,5,0,0,0, 1,4.65,4.65,4.65,0,0,0; "
            <<"\n"<<"Chromosome sizes, "<<fChromosomeDNAContent.size()<<", "; 
            for(int i=0; i<fChromosomeDNAContent.size(); i++) 
            {
                outfile<<(double)fChromosomeDNAContent[i]/1E6;
                if(i<fChromosomeDNAContent.size()-1)
                    outfile<<", "; 
            }
            outfile<<";";
            
    outfile <<"\n"<<"DNA Density, "<<fTotalDNAContent/1E6/421.15<<"; # MBPs/um3 # "
            <<"\n"<<"Cell Cycle Phase, 0; # Simulation without specific cell cycle #"
            <<"\n"<<"DNA Structure, 0, 1; # 0- whole nucleus; 1- wet DNA #"
            <<"\n"<<"In vitro / in vivo, 0; # 0-in vitro # "
            <<"\n"<<"Proliferation status, 0, NULL; # 1- proliferating #"
            <<"\n"<<"Microenvironment, 20, 0.01; # first value defines the temperature in degrees Celsius; second the molar oxygen (O2) concentration in the volume in molarity #"
            <<"\n"<<"Damage definition, 0, 0,"<< fDSBSeparation<<", -1,"<< fDamageThreshold/eV<<"; # direct effects only (0), define in BP(0), DSB seperation,base lesions are not scored(-1), Low energy threshold to induce a SB in eV #"
            <<"\n"<<"Time, 0; "
            <<"\n"<<"Damage and primary count,  1, 1 ;"
            <<"\n"<<"Data entries, " <<"1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0" <<";"
            <<"\n"<<"Additional information, NULL; "
            <<"\n"<<"***EndOfHeader***;"
            <<"\n\n";
    outfile.close(); 


}
