// Scorer for pdb4dna
//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:											      *
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

#include "TsScorePDB4DNA.hh"
#include "TsTrackInformation.hh"

#include "G4SystemOfUnits.hh"

#include <map>

TsScorePDB4DNA::TsScorePDB4DNA(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                               G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
    SetUnit("");
    
    fPDBFileName = fPm->GetStringParameter(GetFullParmName("PDB4DNAFileName"));
	std::fstream in;
	in.open(fPDBFileName,std::ios::in);
	if (!in.is_open() ) {
		G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
		G4cerr << "Input file: " << fPDBFileName << " cannot be opened for Scorer name: " << GetName() << G4endl;
		fPm->AbortSession(1);
	}
	in.close();
	
    fMoleculeList=NULL;
    fBarycenterList=NULL;
    
    G4int verbosity=0;
    unsigned short int isDNA;
    fMoleculeList = fPDBlib.Load(fPDBFileName,isDNA,verbosity);
    
    if (fMoleculeList!=NULL && isDNA==1) {
        fPDBlib.ComputeNbNucleotidsPerStrand(fMoleculeList);
        fBarycenterList=fPDBlib.ComputeNucleotideBarycenters(fMoleculeList );
    }
    
    if (fMoleculeList!=NULL) {
        std::cout << "Read file " << fPDBFileName << std::endl;
    }
    
    // Default parameters
    fThresDistForDSB = 10;
    fThresEdepForSSB = 8.22 * eV;
    
    fNbOfAlgo = 1;
    
    if ( fPm->ParameterExists(GetFullParmName("MinimumDistanceForDSB")) )
        fThresDistForDSB = fPm->GetIntegerParameter(GetFullParmName("MinimumDistanceForDSB"));
    if ( fPm->ParameterExists(GetFullParmName("LowerEnergyForSamplingSSB")) )
        fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("LowerEnergyForSamplingSSB"), "Energy");
    
    // This is for variance reduction
    if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
    	fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
    
    fThresEdepForSSB /= eV;
    
    fNtuple->RegisterColumnI(&fEventID, "Event number");
    fNtuple->RegisterColumnI(&fSSB,     "Single strand breaks");
    fNtuple->RegisterColumnI(&fDSB,     "Double strand breaks");
    
    SuppressStandardOutputHandling();
    
}


TsScorePDB4DNA::~TsScorePDB4DNA() {
}


G4bool TsScorePDB4DNA::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    G4double edep = aStep->GetTotalEnergyDeposit()/eV;
    
    if ( edep > 0 ) {
        G4StepPoint* preStep = aStep->GetPreStepPoint();
        G4ThreeVector pos = preStep->GetPosition();
        G4double x = pos.x()/nanometer;
        G4double y = pos.y()/nanometer;
        G4double z = pos.z()/nanometer;
        
        G4int numStrand = 0;
        G4int numNucl = 0;
        G4int intResidue = -1;
        unsigned short int hit = fPDBlib.ComputeMatchEdepDNA(fBarycenterList, fMoleculeList, x*10, y*10, z*10,
                                                             numStrand, numNucl, intResidue);
        
        if ( 1 == hit ) {
            if ( ( 0 == intResidue ) || ( 1 == intResidue ) ) { //Edep in Phosphate or sugar
                G4int index = 1;
                if ( 1 < fNbOfAlgo ) {
                    TsTrackInformation* trackInformation = (TsTrackInformation*)aStep->GetTrack()->GetUserInformation();
                    index = trackInformation->GetSplitTrackID();
                }
                
                if ( 2 < index ) {
                    if ( numStrand == 1 )
                        fVEdepStrand1[index-3][numNucl] += edep;
                    else
                        fVEdepStrand2[index-3][numNucl] += edep;
                    
                } else {
                    if ( numStrand == 1 )
                        for ( int i = 0; i < fNbOfAlgo; i++ )
                            fVEdepStrand1[i][numNucl] += edep;
                    else
                        for ( int i = 0; i < fNbOfAlgo; i++ )
                            fVEdepStrand2[i][numNucl] += edep;
                }
            }
        }
        return true;
    }
    return false;
}


void TsScorePDB4DNA::UserHookForEndOfEvent() {

    fEventID = GetEventID();
    
    for ( int i = 0; i < fNbOfAlgo; i++ ) {
        G4int sb[2] = {0, 0};
        ComputeStrandBreaks(sb, i);
        fSSB = sb[0];
        fDSB = sb[1];
        fNtuple->Fill();
    }
    
    fVEdepStrand1.erase(fVEdepStrand1.begin(), fVEdepStrand1.end());
    fVEdepStrand2.erase(fVEdepStrand2.begin(), fVEdepStrand2.end());
}


// This class was taken from Geant4/examples/extended/medical/dna/pdb4dna

void TsScorePDB4DNA::ComputeStrandBreaks(G4int* sb, G4int cluster)
{
    // sb quantities
    //
    G4int ssb1=0;
    G4int ssb2=0;
    G4int dsb=0;
    
    // nucleotide id and energy deposit for each strand
    G4int nucl1;
    G4int nucl2;
    G4double edep1;
    G4double edep2;
    
    //Read strand1
    //
    while ( !fVEdepStrand1[cluster].empty() )
    {
        nucl1 = fVEdepStrand1[cluster].begin()->first;
        edep1 = fVEdepStrand1[cluster].begin()->second;
        fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
        
        // SSB in strand1
        //
        if ( edep1 >= fThresEdepForSSB )
        {
            ssb1++;
        }
        
        // Look at strand2
        //
        if ( !fVEdepStrand2[cluster].empty() )
        {
            do
            {
                nucl2 = fVEdepStrand2[cluster].begin()->first;
                edep2 = fVEdepStrand2[cluster].begin()->second;
                if ( edep2 >= fThresEdepForSSB )
                {
                    ssb2++;
                }
                fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
            } while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fVEdepStrand2[cluster].empty()) );
            
            // no dsb
            //
            if ( nucl2-nucl1 > fThresDistForDSB )
            {
                fVEdepStrand2[cluster][nucl2]=edep2;
                if ( edep2 >= fThresEdepForSSB )
                {
                    ssb2--;
                }
            }
            
            // one dsb
            //
            if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
            {
                if ( ( edep2 >= fThresEdepForSSB ) &&
                    ( edep1 >= fThresEdepForSSB ) )
                {
                    ssb1--;
                    ssb2--;
                    dsb++;
                }
            }
        }
    }
    
    // End with not processed data
    //
    while ( !fVEdepStrand1[cluster].empty() )
    {
        nucl1 = fVEdepStrand1[cluster].begin()->first;
        edep1 = fVEdepStrand1[cluster].begin()->second;
        if ( edep1 >= fThresEdepForSSB )
        {
            ssb1++;
        }
        fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
    }
    
    while ( !fVEdepStrand2[cluster].empty() )
    {
        nucl2 = fVEdepStrand2[cluster].begin()->first;
        edep2 = fVEdepStrand2[cluster].begin()->second;
        if ( edep2 >= fThresEdepForSSB )
        {
            ssb2++;
        }
        fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
    }
    sb[0]=ssb1+ssb2;
    sb[1]=dsb;
}
