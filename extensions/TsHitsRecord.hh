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

#ifndef TsHitsRecord_hh
#define TsHitsRecord_hh


#include "G4UnitsTable.hh"
#include <stdint.h>

using namespace CLHEP;
using namespace std;

class TsHitsRecord 
{
public:
    TsHitsRecord();
    ~TsHitsRecord();

    // Set method
    inline G4ThreeVector GetPosition() {return fPos;}
    inline G4double GetPosX() {return fx;}
    inline G4double GetPosY() {return fy;}
    inline G4double GetPosZ() {return fz;}   
    inline G4double GetEdep() {return fEdep;}
    inline G4double GetParticleEnergy() {return fParticleEnergy;}
    inline G4int GetChromosomeID() {return fChromosomeID;}
    inline G4int GetBasePairID() {return fBasePairID;}
    inline G4int GetEventID() {return fEventID;}
    inline G4int GetRunID() {return fRunID;}
    inline G4int  GetDSBType() {return DSBType;}
    inline G4bool GetIsDirectDamage () {return fIsDirectDamage;}
	inline G4bool GetFlagDSB() {return flagDSB;}
	inline G4bool GetFlagSSB() {return flagSSB;}
    inline G4bool GetFlagSSB_P() {return flagSSB_P;}
    inline G4bool GetFlagDSB_P() {return flagDSB_P;}
    inline G4bool GetFlagDSB_PP() {return flagDSB_PP;}
    inline G4String GetParticleName() {return fParticleName;}
    inline G4String GetVolumeName() {return fVolumeName;}


    // Get method
    inline void SetPosition (G4ThreeVector aPos) {fPos = aPos;}
	inline void SetPosX(G4double aPosX) {fx = aPosX;}
	inline void SetPosY(G4double aPosY) {fy = aPosY;}
	inline void SetPosZ(G4double aPosZ) {fz = aPosZ;}
    inline void SetEdep (G4double aEdep) {fEdep = aEdep;}
    inline void SetParticleEnergy (G4double aParticleEnergy) {fParticleEnergy = aParticleEnergy;}
    inline void SetParticleName (G4String aParticleName) {fParticleName = aParticleName;}
    inline void SetChromosomeID(G4int aID) {fChromosomeID = aID;}
    inline void SetBasePairID(G4int aID) {fBasePairID = aID;}
    inline void SetEventID (G4int aEventID) {fEventID = aEventID;}
    inline void SetRunID (G4int aRunID) {fRunID = aRunID;}
    inline void SetVolumeName (G4String aVolumeName) {fVolumeName = aVolumeName;}
    inline void SetIsDirectDamage (G4bool aIsDirectDamage) {fIsDirectDamage = aIsDirectDamage;}
	inline void SetFlagDSB (G4bool aflag) {flagDSB = aflag;}
	inline void SetFlagSSB (G4bool aflag) {flagSSB = aflag;}
    inline void SetFlagSSB_P (G4bool aflag) {flagSSB_P  = aflag;}
    inline void SetFlagDSB_P (G4bool aflag) {flagDSB_P  = aflag;}
    inline void SetFlagDSB_PP(G4bool aflag) {flagDSB_PP = aflag;}
    inline void SetDSBType (G4int type) {DSBType= type;}


    
protected:

 
    
private:
        G4double fEdep ;  //edep (eV)
        G4double fParticleEnergy;
        G4String fParticleName;
		G4String fVolumeName ;
        G4int fChromosomeID;
        G4int fBasePairID;
        G4int fEventID ; 
        G4int fRunID ;  
		G4double fx ;  //x (nm)
        G4double fy ;  //y (nm)
        G4double fz ;  //z (nm)
        G4ThreeVector fPos;
        G4bool fIsDirectDamage; // true: direct damage; false: indirect damage

		G4bool flagDSB;
		G4bool flagSSB;
        G4int  DSBType;        // -1: not a DSB; 0: dir DSB; 1: indir DSB; 2: hybrid DSB

        G4bool flagComplexity;
        G4bool flagSSB_P;   // SSB+
        G4bool flagDSB_P;   // DSB+
        G4bool flagDSB_PP;  // DSB++
};



#endif