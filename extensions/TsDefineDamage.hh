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

#ifndef TsDefineDamage_hh
#define TsDefineDamage_hh

#include <stdint.h>

#include "TsVNtupleScorer.hh"
#include "G4UnitsTable.hh"
#include "TsHitsRecord.hh"


using namespace CLHEP;
using namespace std;

class TsDefineDamage 
{
public:
    TsDefineDamage();
    ~TsDefineDamage();

    void SeparateHitsOnDifferentDNAStrands(vector<TsHitsRecord*> Hits, vector<TsHitsRecord*> &HitsBack1, G4double IsStrand1);
    void DefineDSBorSSB(vector<TsHitsRecord*> HitsBack1, vector<TsHitsRecord*> HitsBack2,
                        vector <pair<TsHitsRecord*, TsHitsRecord*>> &DSB_pairs,
                        vector<TsHitsRecord*> & SSBonStrand1,vector<TsHitsRecord*> & SSBonStrand2);
    void GetSSBonStrand(vector<TsHitsRecord*> Hits, vector<TsHitsRecord*> & SSBonStrand);
    void DefineDSBType (TsHitsRecord * &hit1, TsHitsRecord * &hit2);
    void ExcludeShortDNAFragments(vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, vector<TsHitsRecord*> SSBonStrand1,vector<TsHitsRecord*> SSBonStrand2, G4int ChromosomeID); 
    void OutputDNAdamageTupleHeader(string filename);
    void OutputDNAdamageTuple(vector<TsHitsRecord*> HitsBack1, vector<TsHitsRecord*> HitsBack2, string filename);
    void OutputDNAdamageSummary(string filename);
    void SingleStrandDamageSiteSeperateration(vector <TsHitsRecord*> &Hits, vector <G4ThreeVector> &XYZPosition, vector <double  > &ChromosomePosition, 
                                              vector <double  > &Cause, vector <int > &DamageSiteDSBflag, vector<vector <int *> > &DamageSites);
    void SeparateDamageSitesOnOneStrand(vector <TsHitsRecord*> Hits, vector <G4ThreeVector> &XYZPosition, vector <double  > &ChromosomePosition, 
                                        vector <double  > &Cause, vector <int > &DamageSiteDSBflag, vector< vector <vector <int > > >  &DamageSites);
    void OutputSDDFile( bool OnlyIncludeDSBinSDD, vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, 
                        vector<TsHitsRecord*> SSBonStrand1,vector<TsHitsRecord*> SSBonStrand2,
                        string filename, G4int EventID, G4int ExposureID, G4int ChromosomeID) ;                       
    void OutputSDDHeader(string filename);
    G4bool CasueDirectDamge(G4double EnergyDeposition);


    // Set method
    inline void SetDamageThreshold(G4double aDamageThreshold) { fDamageThreshold = aDamageThreshold;}
    inline void SetUseLinearProbabilitythreshold(G4bool aUse_Lin_Pro_th)    { fUseLinearProbabilitythreshold = aUse_Lin_Pro_th;}
    inline void SetLinearProbability_lower_limit(G4double aDamageThreshold) { fLinearProbability_lower_limit = aDamageThreshold;}
    inline void SetLinearProbability_upper_limit(G4double aDamageThreshold) { fLinearProbability_upper_limit = aDamageThreshold;}
    inline void SetDSBSeparation(G4int aDSBSeparation) {fDSBSeparation = aDSBSeparation;}
    inline void SetDefineComplexity(G4bool aDefineComplexity) {fDefineComplexity = aDefineComplexity;}
    inline void SetComplexitySeparation(G4int aComplexitySeparation) {fComplexitySeparation = aComplexitySeparation;}
    inline void SetExcludeShortFragment(G4bool aExcludeShortFragment) {fExcludeShortFragment = aExcludeShortFragment;}
    inline void SetLowerFragmentDetectionThreshold(G4int aThreshold) { fLowerFragmentDetectionThreshold = aThreshold;}
    inline void SetUpperFragmentDetectionThreshold(G4int aThreshold) { fUpperFragmentDetectionThreshold = aThreshold;}
    inline void SetEventID(G4int aEventID) {fEventID=aEventID;}
    inline void SetEdep(G4double aEdep) {fEdep=aEdep;}
    inline void SetLET (G4double aLET) {fLET = aLET;}
    inline void SetZetaBeta_sq(G4double aZetaBeta_sq) {fZetaBeta_sq = aZetaBeta_sq;}
    inline void SetNucleusMass(G4double aMass) {fNucleusMass=aMass;}
    inline void SetTotalDNAContent(G4double aTotalContent) {fTotalDNAContent=aTotalContent;}
    inline void SetChromosomeDNAContent(vector<G4int> aDNAContent) { fChromosomeDNAContent = aDNAContent;}
    //inline void SetChromosomeDNAContentSum(vector<G4double> aDNAContentSum) { fChromosomeDNAContentSum = aDNAContentSum;}
    inline void SetOnlyIncludeDSBinSDD (G4bool aInlcude) {fOnlyIncludeDSBinSDD = aInlcude;}


    // Get method
    inline G4int GetNumSSB() {return numSSB;}
    inline G4int GetNumDSB() {return numDSB;}
 
    
private:
    G4double fDamageThreshold;
    G4bool   fUseLinearProbabilitythreshold;
    G4double fLinearProbability_lower_limit;
    G4double fLinearProbability_upper_limit;
    G4int    fDSBSeparation;
    G4bool   fDefineComplexity;
    G4int    fComplexitySeparation;
    G4bool   fExcludeShortFragment;
    G4int    fLowerFragmentDetectionThreshold;
    G4int    fUpperFragmentDetectionThreshold;
    G4int    fDamageSiteSize;
    G4bool   fOnlyIncludeDSBinSDD;

    G4int fEventID;
    G4double fEdep; //MeV
    G4double fLET;  //keV/um
    G4double fNucleusMass; //g
    G4double fTotalDNAContent;  //bp
    G4double fZetaBeta_sq;

    G4int numSB;
    G4int numSB_dir;
    G4int numSB_indir;
    G4int numSSB;
    G4int numSSB_dir;
    G4int numSSB_indir;
    G4int numDSB;
    G4int numDSB_dir;
    G4int numDSB_indir;
    G4int numDSB_hybrid;
    // complex damages
    G4int numSSB_P;
    G4int numDSB_P;
    G4int numDSB_PP;


	G4int Excluded_numSB;
    G4int Excluded_numSB_dir;
    G4int Excluded_numSB_indir;
    G4int Excluded_numSSB;
    G4int Excluded_numSSB_dir;
    G4int Excluded_numSSB_indir;
    G4int Excluded_numDSB;
    G4int Excluded_numDSB_dir;
    G4int Excluded_numDSB_indir;
    G4int Excluded_numDSB_hybrid;
    // complex damages
    G4int Excluded_numSSB_P;
    G4int Excluded_numDSB_P;
    G4int Excluded_numDSB_PP;

    std::vector<G4int> fChromosomeDNAContent; //in unit of bp


    G4int lastOutputEventID;
    G4int lastOutputExposureID;
    G4int lastDSBonlyOutputEventID;
    G4int lastDSBonlyOutputExposureID;

    
};



#endif