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

#ifndef  TsChromosome_hh
#define  TsChromosome_hh


#include "G4UnitsTable.hh"
#include "TsHitsRecord.hh"
#include <stdint.h>

using namespace CLHEP;
using namespace std;

class  TsChromosome 
{
public:
     TsChromosome();
    ~ TsChromosome();
    
    G4double               GetTotalDNAContent_MBP() {return totalDNAContent;}
    vector<G4int>          GetDNAcontentofEachChromosome_BP()  {return DNAcontentofEachChromosome_BP;}
    vector<vector<G4int> > GetChromosomeMatrix() { return ChromosomeMatrix;}
    vector<vector<G4int> > SplitChromosome(vector <pair< TsHitsRecord*,  TsHitsRecord*>> DSB_pairs);
    G4int CountDNAFrangmentsWithSize( vector<vector<G4int> > DNAfragments, G4int LowerThreshold, G4int UpperThreshold);

    
protected:

 
    
private:
    
    
    G4double totalDNAContent;
    vector<G4int> DNAcontentofEachChromosome_BP; // In unit of BP
    vector<vector<G4int> > ChromosomeMatrix;     // ChromosomeID, 0, Chromosome DNA content(Mbp) 
};



#endif