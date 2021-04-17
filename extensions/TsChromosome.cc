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


#include <stdint.h>
#include "TsChromosome.hh"
#include "TsHitsRecord.hh"

using namespace CLHEP;
using namespace std;


 TsChromosome:: TsChromosome()
{
    // // Human DNA content
    // DNAcontentofEachChromosome_BP = {252823200,  252823200,  248157000,  248157000,  204040200,  204040200,  
    //                                 195556200,  195556200,  184951200,  184951200,  174770400,  174770400,  
    //                                 162468600,  162468600,  149318400,  149318400,  143379600,  143379600,  
    //                                 138289200,  138289200,  137440800,  137440800,  135319800,  135319800,  
    //                                 116655000,  116655000,  108595200,  108595200,  102656400,  102656400,  
    //                                 90778800, 90778800,  80173800,  80173800,  77628600,  77628600,  
    //                                 65326800,  65326800,  63630000,  63630000,  47934600,  47934600, 
    //                                 50479800,  50479800,  58963800,  158226600};

    // Mouse DNA content
    DNAcontentofEachChromosome_BP = {195471971,  195471971,  182113224,  182113224,  160039680,  160039680,  
                                    156508116,  156508116,  151834684,  151834684,  149736546,  149736546,  
                                    145441459,  145441459,  129401213,  129401213,  124595110,  124595110,  
                                    130694993,  130694993,  122082543,  122082543,  120129022,  120129022,  
                                    120421639,  120421639,  124902244,  124902244,  104043685,  104043685,  
                                    98207768, 98207768,  94987271,  94987271,  90702639,  90702639,  
                                    61431566,  61431566,  171031299,  91744698};
    
    totalDNAContent =0;
    for(G4int i=0; i<DNAcontentofEachChromosome_BP.size(); i++ )
    {
        totalDNAContent += (G4double) DNAcontentofEachChromosome_BP[i]/1E6;
        vector<G4int> oneMat = {0, DNAcontentofEachChromosome_BP[i]};
        ChromosomeMatrix.push_back(oneMat);
    }
        

}
    
TsChromosome::~ TsChromosome()
{}

vector<vector<G4int> > TsChromosome::SplitChromosome(vector <pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs)
{
    // G4cout << "c1" << G4endl;

    vector<vector<G4int> > SplitedChromosome = GetChromosomeMatrix();
    for(G4int i =0; i<SplitedChromosome.size(); i++)
    {
        G4int ChromosomeID = i+1;
        for (G4int iter =0; iter<DSB_pairs.size(); iter++)
        {
            
            G4int DamagedChromosomeID = DSB_pairs[iter].first->GetChromosomeID();
            G4int DamagedBasepairID   = DSB_pairs[iter].first->GetBasePairID();
            if( ChromosomeID == DamagedChromosomeID )
                SplitedChromosome[i].push_back( DamagedBasepairID);
        } 
        sort(SplitedChromosome[i].begin(), SplitedChromosome[i].end());
    }
    return SplitedChromosome;
}

G4int TsChromosome::CountDNAFrangmentsWithSize( vector<vector<G4int> > DNAfragments, G4int LowerFragmentSizeThreshold, G4int UpperFragmentSizeThreshold)
{
    // LowerFragmentSizeThreshold, UpperFragmentSizeThreshold  in the unit of BP

    G4int NumberOfInterestedFragments = 0;
    for(G4int i =0; i<DNAfragments.size(); i++)
    {
        G4double FrangmentSize =0;
        vector<G4int> AllFreeEndofChromosome = DNAfragments[i];
        for(G4int j =1; j<AllFreeEndofChromosome.size(); j++)
        {
            FrangmentSize = AllFreeEndofChromosome[j] - AllFreeEndofChromosome[j-1];
            if(FrangmentSize>=LowerFragmentSizeThreshold && FrangmentSize<=UpperFragmentSizeThreshold )
                NumberOfInterestedFragments++; 
        }
    }
    return NumberOfInterestedFragments;
}

