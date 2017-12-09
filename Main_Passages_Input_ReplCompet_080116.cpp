%Copyright (c) <2017>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory
%Written by Tanya Kostova Vassilevska, kostova@llnl.gov, tan.v.kos@gmail.com
%LLNL Release number LLNL-CODE-727823
%All rights reserved.

%This file is part of <VirESS>, including the main routine. For details, see the comments. 

%Licensed under the Apache License, Version 2.0 (the “Licensee”); you may not use this file except in compliance with the License.  You may obtain a copy of the License at:  http://www.apache.org/licenses/LICENSE-2.0

%Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the license.


/*
 * Author Tanya Kostova Vassilevska
 *
 *
 */

This code simulates the infection and mutation of a virus in a population of cells. The virus is represented by its genetic sequence. 

#include <stdlib.h>
#include <math.h>
#include<vector>
#include <iostream>
#include <fstream>
#include <string>
#include<vector>
#include<set>
#include <cstdio>
#include<list>

#include "Objects_08_10_15.h"
#include "MTRand.h"


using namespace std;

const double MutationRate=0.00005;
const double TransitionProb=0.66;
const double seqFT_basic=1;
const double seqFT_max=5;
const double seqFR5_basic=1;
const double seqFR5_max=5;
const double seqFR3_basic=1;
const double seqFR3_max=5;
const double seqFR_basic=1;
const double seqFR_max=5;
const double seqFE_basic=1;
const double seqFE_max=5;

// REFERENCE sequence; this is GLOBAL

//char* ReferenceSequence= new char[NUMBER_OF_BASES];  

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////TRANSLATE////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double LKq(int L, int k, double q)// calculates binomial coefficient L  over k multiplied by q^k; in the code L is length of genome, q is the probability of one point mutation;
{
  double LoK=1;
  if(k==0){LoK=1; return LoK;}
  if(k>L){cout<<"k>L"; exit(1);}

  for(int i=0; i<k;i++)
    {
      int ii=i;
      double f=(L-ii)*q/(ii+1); 
      LoK=LoK*f; 
    }
     
  return LoK;
  
}


int MutationsNumber(double MutRate, int SeqLength, double rnd)//given a random number rnd, calculates  NM so that rnd<=(L over NM)(1-q)^{L-NM}q^NM; NM is a random number of point mutations happening during one replication event
{
  int NM=0;   //NM is the number of mutations
  int L=SeqLength;
   
  double prob=pow(1-MutRate, SeqLength);
 
  if(prob>=rnd)return NM;
  else
    {
      while(prob<rnd)
	{
	  NM++;
	  prob+=pow(1-MutRate,L-NM)*LKq(L,NM,MutRate); // \sum (L over NM) (1-q)^{L-NM}q^NM probability that less than NM+1 mutations and more than or equal to NM mutations occurred
	} 
    };
  return NM; 
}

void Translate(int SeqL, char* &ntSeq, vector<char> &aaSeq)//translate a given nucleotide sequence into aminoacid sequence; SeqL is the length of the nucleotide sequence ; NOT sure why char* &ntSeq by reference?
{
  int aux=0;
  int l=0;
  char NT[3];

 
  while(aux<SeqL)
    {  
      NT[0]=ntSeq[aux];
      aux++;
      NT[1]=ntSeq[aux];
      aux++;
      NT[2]=ntSeq[aux];
      aux++;

      if(NT[0]=='t'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='c')
	{
	  aaSeq[l]='F'; //cout<<aaSeq[l];
	}
      if(NT[0]=='c'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='g'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='g')
	{
	  aaSeq[l]='L'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='a'&&NT[1]=='t'&&NT[2]=='a')
	{
	  aaSeq[l]='I'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='t'&&NT[2]=='g')
	{
	  aaSeq[l]='M'; //cout<<aaSeq[l];
	}

      if(NT[0]=='g'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='g')
	{
	  aaSeq[l]='V'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='g')
	{
	  aaSeq[l]='S'; //cout<<aaSeq[l];
	}
      if(NT[0]=='c'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='g')
	{
	  aaSeq[l]='P'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='g')
	{
	  aaSeq[l]='T'; //cout<<aaSeq[l];
	}
	
      if(NT[0]=='g'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='g')
	{
	  aaSeq[l]='A'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='a'&&NT[2]=='c')
	{
	  aaSeq[l]='Y'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='a'&&NT[2]=='g')
	{
	  aaSeq[l]='*'; //cout<<aaSeq[l];
	}
      if(NT[0]=='c'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='a'&&NT[2]=='c')
	{
	  aaSeq[l]='H'; //cout<<aaSeq[l];
	}
      if(NT[0]=='c'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='a'&&NT[2]=='g')
	{
	  aaSeq[l]='Q'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='a'&&NT[2]=='c')
	{
	  aaSeq[l]='N'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='a'&&NT[2]=='g')
	{
	  aaSeq[l]='K'; //cout<<aaSeq[l];
	}
      if(NT[0]=='g'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='a'&&NT[2]=='c')
	{
	  aaSeq[l]='D'; //cout<<aaSeq[l];
	}
      if(NT[0]=='g'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='a'&&NT[2]=='g')
	{
	  aaSeq[l]='E'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='g'&&NT[2]=='c')
	{
	  aaSeq[l]='C'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='a')
	{
	  aaSeq[l]='*'; //cout<<aaSeq[l];
	}
      if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='g')
	{
	  aaSeq[l]='W'; //cout<<aaSeq[l];
	}
      if(NT[0]=='c'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='g')
	{
	  aaSeq[l]='R'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='g'&&NT[2]=='c')
	{
	  aaSeq[l]='S'; //cout<<aaSeq[l];
	}
      if(NT[0]=='a'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='g'&&NT[2]=='g')
	{
	  aaSeq[l]='R'; //cout<<aaSeq[l];
	}
      if(NT[0]=='g'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='g')
	{
	  aaSeq[l]='G'; //cout<<aaSeq[l];
	}
            l++;
    }
 	
};

///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////END TRANSLATE ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////CALCULATE FITNESS OF FULL SEQUENCE /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

//

void CalculateFitnessOfFullSequence(int NumberOfHF_5T,int NumberOfHF_5R, vector<int>vHF_5T, vector<int>vHF_5R, vector<char>HF_5T,vector<char>HF_5R, int NumberOfHF_3T,int NumberOfHF_3R, vector<int>vHF_3T,vector<int>vHF_3R, vector<char>HF_3T, vector<char>HF_3R, int NumberOfHF_STR, vector<int>vHF_STR, vector<char>HF_STR,int NumberOfHF_NonSTR,vector<int>vHF_NonSTR,vector<char>HF_NonSTR,  vector<char>Sequence, double& SequenceFitnessOfTranslation,double& SequenceFitnessOfEntry,double& SequenceFitnessOfReplication, int Start_3Prime, vector<double> FPrA,vector<double> FPrC, vector<double> FPrG, vector<double> FPrT,vector<double> TPrA,vector<double> TPrC, vector<double> TPrG, vector<double> TPrT, int fpr, int tpr,int sstrn,int nstrn, vector<map<string,double> > STRCodonFreqs, vector<map<string,double> > NSTRCodonFreqs, int Start_Str_Prot_Region,int Start_NStr_Prot_Region)
{    
  SequenceFitnessOfTranslation=seqFT_basic;//these are default values; they can change to 0 or to 1
  SequenceFitnessOfEntry=seqFE_basic;
  SequenceFitnessOfReplication=seqFR_basic;
  
  //for(int y=0; y<NumberOfHF_5T; y++){cout<<"   "<<y<<"   "<<vHF_5T[y]<<"   "; }
  //for(int y=0; y<NumberOfHF_5R; y++){cout<<"   "<<y<<"   "<<vHF_5R[y]<<"   "; }

  int helper=1;
  int helper1T=0; 
  int helper1R=0;
  
  string str;
  char NT[3];   
  int broT=0;  
  int broR=0;  
  
  vector<double> FPrSequence(fpr);
  vector<double> TPrSequence(tpr);
  
  vector<double> STRSequenceFreq(sstrn);
  vector<double> STRSequenceAAFreq(sstrn);
  
  vector<double> NSTRSequenceFreq(nstrn);
  vector<double> NSTRSequenceAAFreq(nstrn);

  double SequenceFitnessOfReplication5=0;double SequenceFitnessOfReplication3=0; double SequenceFitnessOfReplicationN=0;
  double SequenceFitnessOfTranslation5=0;double SequenceFitnessOfTranslation3=0;

  //cout<<"fpr=  "<<fpr<<endl;  

  for(int i=0; i<fpr;i++)
    {
      if(Sequence[i]=='a')
	{//cout<<i<<"  "<<Sequence[i]<<"           "; 
	  FPrSequence[i]=FPrA[i]; //cout<<FPrA[i]<<endl;
	}
      if(Sequence[i]=='c')
	{//cout<<i<<"  "<<Sequence[i]<<"           ";
	  FPrSequence[i]=FPrC[i];//cout<<FPrC[i]<<endl;
	}
      if(Sequence[i]=='g')
	{//cout<<i<<"  "<<Sequence[i]<<"           ";
	  FPrSequence[i]=FPrG[i];//cout<<FPrG[i]<<endl;
	}
      if(Sequence[i]=='t')
	{//cout<<i<<"  "<<Sequence[i]<<"           ";
	  FPrSequence[i]=FPrT[i];
	  //cout<<FPrT[i]<<endl;
	}
      //the value of FPrSequence is only used in case it is 0, i.e. establisching that the mutation is not observed so far
      
   
      if(!FPrSequence[i])
	{
	  helper=0;
	  //cout<<"at pos. "<< i<<" 0 fitness of repl and transl"<< endl;
	  break; //if helper==0, we will assume that both replication and translation fitness are 0
	}
      else{}
    }
	
  if(NumberOfHF_5T)
    {
      for(int i=0; i<NumberOfHF_5T; i++)
	{
	  if(Sequence[vHF_5T[i]]== HF_5T[i])
	    {
	      helper1T++; //if i=h, then h is x+1, where x is the position in Adam's table
	    }
	  else
	    {}//here calculate how many of the HF for translation position-nt pairs  are actually present in the sequence
	}
    }

  if(NumberOfHF_5R)
    {
      for(int i=0; i<NumberOfHF_5R; i++)
	{
	  if(Sequence[vHF_5R[i]]== HF_5R[i])
	    {
	      helper1R++; //if i=h, then h is x+1, where x is the position in Adam's table
	    }
	  else
	{}	 //here calculate how many of the HF for replication  position-nt pairs  are actually present in the sequence   
	}
    }
    
      //cout<<"helper"<<helper<<endl;
      //cout<<"NumberOfHF_5T  "<<NumberOfHF_5T<<endl;
      //cout<<"NumberOfHF_5R  "<<NumberOfHF_5T<<endl;

  
  if(NumberOfHF_5T&&helper1T==NumberOfHF_5T){SequenceFitnessOfTranslation5=seqFT_max;}// all mutations that make a sequence with high fitness for translation are present!
  else{SequenceFitnessOfTranslation5=seqFT_basic;}

  if(NumberOfHF_5R&&helper1R==NumberOfHF_5R){SequenceFitnessOfReplication5=seqFR_max;}// all mutations that make a sequence with high fitness for repl. are present!
  else{SequenceFitnessOfReplication5=seqFR_basic;}

  if(!helper){SequenceFitnessOfTranslation5=0;SequenceFitnessOfReplication5=0;}

  //cout<<"  SequenceFitnessOfTranslation5  "<< SequenceFitnessOfTranslation5  <<"    SequenceFitnessOfReplication5  "<<SequenceFitnessOfReplication5<<endl;
  //cout<<"helper1T "<<helper1T<<endl;

  
  //=============================================================================================================================

  int ww=Start_Str_Prot_Region; //cout<< "ww=  "<<ww<<endl;
  int  helperE=1;

  int helper2=0;
    
  for(int w=0; w<sstrn; w++)//structural codons numbered from 0 to sstrn-1;
    {
      for(int j=0; j<3; j++)
	{
	  NT[j]=Sequence[ww];	  
	  ww++;
	}     
            
      str="";
      for(int i=0; i<3; i++){str+=NT[i];} //forming a string NOT SURE WHY THIS IS OUTSIDE OF PREV LOOP
      
      STRSequenceFreq[w]=STRCodonFreqs[w][str]; //STRSequenceFreq[w]does not need to be a vector, probably this remained from a previous version of fitness calculation; one could use these frequencies somehow but we are not doing this 
      
      if(!STRSequenceFreq[w])
	{helperE=0;
	//cout<<"helperE=0"<<endl;  
	break;
	}//if the whole codon was never reported, zero fitness; if ANY codon has 0 fitness, stop checking further
    }
  
  //NOTE: in the input data, the full HF codons must be given and they should not have frequency 0!!!

  if(NumberOfHF_STR)
    {
      for(int i=0; i<NumberOfHF_STR; i++)  
	{
	  if(Sequence[vHF_STR[i]]==HF_STR[i])
	    {helper2++;}
	  else{}
	}	    
    }

  if(NumberOfHF_STR&&helper2==NumberOfHF_STR){SequenceFitnessOfEntry=seqFE_max;}
  else{SequenceFitnessOfEntry=seqFE_basic;}
  if(!helperE){SequenceFitnessOfEntry=0;}
   
   

  //==============================================================================

  ww=Start_NStr_Prot_Region;
  
  int helperN=1;

  helper2=0; 
  
  for(int w=0; w<nstrn; w++)
    {
      for(int j=0; j<3; j++)
	{
	  NT[j]=Sequence[ww];
	  ww++; 
	}
         
      str=""; 
      for(int i=0; i<3; i++){str+=NT[i];}//forming a string NOT SURE WHY THIS IS OUTSIDE OF PREV LOOP
      
      NSTRSequenceFreq[w]=NSTRCodonFreqs[w][str];//NSTRSequenceFreq[w]does not need to be a vector, probably this remained from a previous version of fitness calculation; one could use these frequencies somehow but we are not doing this 
      
      if(!NSTRCodonFreqs[w][str])
	{helperN=0;
	//cout<<"helperN=0"<<endl;
	break;
	}//if the whole codon was never reported, zero fitness; if ANY codon has 0 fitness, stop checking further
      
    }
  //NOTE: in the input data, the full HF codons must be given and they should not have frequency 0!!!
  
  if(NumberOfHF_NonSTR)
    {
      for(int i=0; i<NumberOfHF_NonSTR; i++)  
	{
	  if(Sequence[vHF_NonSTR[i]]==HF_NonSTR[i])
	    {helper2++;}
	  else{}
	}	
    }

    if(NumberOfHF_NonSTR&&helper2==NumberOfHF_NonSTR){SequenceFitnessOfReplicationN=seqFR_max;}
    else{SequenceFitnessOfReplicationN=seqFR_basic;}
    if(!helperN){SequenceFitnessOfReplicationN=0;}
   
  
    //==============================================================================================
   
    
    helper=1;
    helper1T=0; 
    helper1R=0;
    
    
    for(int i=0; i<tpr; i++)
      {
	if(Sequence[Start_3Prime+i]=='a')
	  {TPrSequence[i]=TPrA[i];}
	if(Sequence[Start_3Prime+i]=='c')
	  {TPrSequence[i]=TPrC[i];}
	if(Sequence[Start_3Prime+i]=='g')
	  {TPrSequence[i]=TPrG[i];}
	if(Sequence[Start_3Prime+i]=='t')
	  {TPrSequence[i]=TPrT[i];}
	
	
	if(!TPrSequence[i])
	  {
	    helper=0;
	    break;
	  }
	
	else{}
      }
    if(NumberOfHF_3T)
      { 
	for(int i=0; i<NumberOfHF_3T; i++)
	  {
	    if(Sequence[vHF_3T[i]]== HF_3T[i])
	      {
		helper1T++; //if i=h, then h is x+1, where x is the position in Adam's table
	      }
	    else
	      {}//here calculate how many of the HF for translation position-nt pairs  are actually present in the sequence
	  }
      }
    if(NumberOfHF_3R)
      {
	for(int i=0; i<NumberOfHF_3R; i++)
	  {
	    if(Sequence[vHF_3R[i]]== HF_3R[i])
	      {
		helper1R++; //if i=h, then h is x+1, where x is the position in Adam's table
	      }
	    else
	      {}	 //here calculate how many of the HF for replication  position-nt pairs  are actually present in the sequence   
	  }
      }

    //cout<<"helper"<<helper<<endl;
    //cout<<"NumberOfHF_5T  "<<NumberOfHF_5T<<endl;
    //cout<<"NumberOfHF_5R  "<<NumberOfHF_5T<<endl;
    
    
    if(NumberOfHF_3T&&helper1T==NumberOfHF_3T){SequenceFitnessOfTranslation3=seqFT_max;}// all mutations that make a sequence with high fitness for translation are present!
    else{SequenceFitnessOfTranslation3=seqFT_basic;}
    
    if(NumberOfHF_3R&&helper1R==NumberOfHF_3R){SequenceFitnessOfReplication3=seqFR_max ;}// all mutations that make a sequence with high fitness for repl. are present!
    else{SequenceFitnessOfReplication3=seqFR_basic;}
    
    
    if(!helper){SequenceFitnessOfTranslation3=0;SequenceFitnessOfReplication3=0;}
    // ***************************************************************************************	  
    if(NumberOfHF_5T&&NumberOfHF_3T||SequenceFitnessOfTranslation5*SequenceFitnessOfTranslation3==0)
      {
	SequenceFitnessOfTranslation=min(SequenceFitnessOfTranslation5,SequenceFitnessOfTranslation3);
      }
    else
      {
	SequenceFitnessOfTranslation=max(SequenceFitnessOfTranslation5,SequenceFitnessOfTranslation3);
	
      }
    // ***************************************************************************************	  
    if(NumberOfHF_NonSTR&&NumberOfHF_5R&&NumberOfHF_3R||SequenceFitnessOfReplication5*SequenceFitnessOfReplication3*SequenceFitnessOfReplicationN==0)
      {
	SequenceFitnessOfReplication=min(SequenceFitnessOfReplication5, min(SequenceFitnessOfReplication3,SequenceFitnessOfReplicationN));
      }
    else
      {
	if(NumberOfHF_NonSTR&&NumberOfHF_5R){SequenceFitnessOfReplication=min(SequenceFitnessOfReplication5,SequenceFitnessOfReplicationN);}
	if(NumberOfHF_NonSTR&&NumberOfHF_3R){SequenceFitnessOfReplication=min(SequenceFitnessOfReplication3,SequenceFitnessOfReplicationN);}
	if(NumberOfHF_3R&&NumberOfHF_5R){SequenceFitnessOfReplication=min(SequenceFitnessOfReplication5,SequenceFitnessOfReplication3);}
	if(NumberOfHF_3R){SequenceFitnessOfReplication=SequenceFitnessOfReplication3;}
	if(NumberOfHF_5R){SequenceFitnessOfReplication=SequenceFitnessOfReplication5;}
	if(NumberOfHF_NonSTR){SequenceFitnessOfReplication=SequenceFitnessOfReplicationN;}
	if(!NumberOfHF_NonSTR&&!NumberOfHF_5R&&!NumberOfHF_3R){SequenceFitnessOfReplication=seqFR_basic;}
      }
    
	//cout<<"SequenceFitnessOfTranslation  "<<SequenceFitnessOfTranslation<<" SequenceFitnessOfReplication  "<< SequenceFitnessOfReplication<<" SequenceFitnessOfEntry  "<<SequenceFitnessOfEntry<<endl;			      
}

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////END CALCULATE FITNESS OF SEQUENCE////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //***************************************************************************************************************//
    //******************************************** MAIN *************************************************************//
    //***************************************************************************************************************//
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    main()
{
  MTRand_closed(time(NULL));//seed random number generator
  
  ///////////////////////ALPHA and BETA used in fitness calculations///////////////
  double ALPHA=1;
  double BETA=0.;
  string str;
  string ToRead;   //used to jump over a field in the input file  that is not needed for the calculation
 
  double FTR=0; double FE=0; double FREPL=0; double FREPL1=0;
  
 /////////////////////////// READ INPUT PARAMETERS FROM FILE ///////////////////////////////////////////////////
 
  ifstream simInput("SimulationInputWitInoculumOptions1.txt");

  int NumberOfSimulations;
  simInput>>ToRead>>NumberOfSimulations;
  cout<<"NumberOfSimulations "<<NumberOfSimulations<<endl;

  int NumberOfPassages;
  simInput>>ToRead>>NumberOfPassages;
  cout<< "NumberOfPassages "<<NumberOfPassages<<endl;

  int NumberOfDominants;
  simInput>>ToRead>>NumberOfDominants;
  cout<<"NumberOfDominants "<<NumberOfDominants<<endl;//for the output: how many dominants to output

  int NumberOfCellsInCulture;
  simInput>>ToRead>>NumberOfCellsInCulture;
  cout<<"NumberOfCellsInCulture "<<NumberOfCellsInCulture<<endl;

  int MaxRepeatInfection;
  simInput>>ToRead>>MaxRepeatInfection;
  cout<<"MaxRepeatInfection "<<MaxRepeatInfection<<endl;

  double MOI;
  simInput>>ToRead>>MOI;
  cout<<"MOI "<<MOI<<endl;

  int InoculumInput;
  simInput>>ToRead>>InoculumInput;
  cout<<" InoculumInput "<<InoculumInput<<endl;


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// READ INFO about structure of gene sequence ///////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<" READ INFO about structure of gene sequence"<<endl;

  ifstream inRefSeq("ReferenceSeq.txt");
  ifstream inSeqRegions("SequenceRegions.txt");

    
  int NUMBER_OF_BASES; 
  int Start_5Prime; 
  int End_5Prime; 
  int Start_3Prime; 
  int NumberOfCodingRegions; 
  int End_3Prime;
  int Start_Str_Prot_Region; 
  int End_Str_Prot_Region; 
  int Start_NStr_Prot_Region;
  int End_NStr_Prot_Region;
  
  
  if(!inRefSeq){cout<<"NO FILE for Reference Sequence"<<endl; exit(1);}
  if(!inSeqRegions){cout<<"NO FILE for Sequence regions"<<endl; exit(1);}

  inSeqRegions>>ToRead>>NUMBER_OF_BASES; 
  //cout<<NUMBER_OF_BASES<<" "<<endl;
  inSeqRegions>>ToRead>>Start_5Prime;
  //cout<<Start_5Prime<<" "<<endl;
  inSeqRegions>>ToRead>>End_5Prime;
  //cout<<End_5Prime<<" "<<endl;
  inSeqRegions>>ToRead>>Start_3Prime;
  //cout<<Start_3Prime<<" "<<endl;
  inSeqRegions>>ToRead>>NumberOfCodingRegions;
  //cout<<NumberOfCodingRegions<<"   "<<endl;
  End_3Prime=NUMBER_OF_BASES;
  //cout<<End_3Prime<<"   "<<endl;

  int* StartCodingRegion = new int[NumberOfCodingRegions]; //here put the beginnings of each coding region
  int* EndCodingRegion = new int[NumberOfCodingRegions]; //here put the ends of each coding region

   
  for(int i=0;i<NumberOfCodingRegions; i++)
    {
      inSeqRegions>>ToRead>>StartCodingRegion[i];
      inSeqRegions>>ToRead>>EndCodingRegion[i];
    }

 
  inSeqRegions>>ToRead>>Start_Str_Prot_Region>>ToRead>>End_Str_Prot_Region>>ToRead>>Start_NStr_Prot_Region>>ToRead>>End_NStr_Prot_Region;
  
  inSeqRegions.close(); 

 ////////////////////////////////////////////////////////////////////////////////////////
  ///////STRUCTURE SPECIFIC TO FLAVI AND PICORNA (STR BEFORE NONSTR)////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  if(Start_5Prime){cout<<"Wrong Sequence data(1)"<<endl;exit(1);}
  if(End_5Prime!=Start_Str_Prot_Region-1){cout<<"Wrong Sequence data(2)"<<endl;exit(1);}
  if(End_Str_Prot_Region!=Start_NStr_Prot_Region-1){cout<<"Wrong Sequence data(3)"<<endl;exit(1);}
  if(End_NStr_Prot_Region!=Start_3Prime-1){cout<<"Wrong Sequence data(4)"<<endl;exit(1);}



  ////////////////////////////////////////////////////////
  //////// READ REFERENCE SEQUENCE /////////////
  ////////////////////////////////////////////////////////


  
  //ofstream Prot("Proteins.txt");
  
  vector<char*> W;//vector of pointers of character type to keep coding regions
  
  vector<char> ReferenceSequence(NUMBER_OF_BASES);  
  char a;
  
  for(int i=0; i<NUMBER_OF_BASES; i++)
    {
      inRefSeq>>a;
      ReferenceSequence[i]=a; //NOTE that positions in this vector start from 0;
      //cout<<a;
    }
  inRefSeq.close();
    
  //cout<<"end read refseq"<<endl;
  
  /////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////READ and TRANSLATE CODING REGIONS/////////////
  /////////////////////////////////////////////////////////////////////////////////

  char tog;
  vector<vector<char> >Protein;//vector of vectors to keep the virus proteins
    
  for(int i=0;i<NumberOfCodingRegions; i++)
    {
      double y=fmod(EndCodingRegion[i]-StartCodingRegion[i]+1,3); 
      
      if(y)
	{
	  cout<<"unable to translate, remainder not 0"; return 0;
	}
      char* CR=new char[EndCodingRegion[i]-StartCodingRegion[i]+1];
      
      for(int k=0;k<EndCodingRegion[i]-StartCodingRegion[i]+1;k++)
	{
	  CR[k]=ReferenceSequence[k+StartCodingRegion[i]]; 
	  
	}
      
      
      int pp= (EndCodingRegion[i]-StartCodingRegion[i]+1);//pp is the length of the next protein
      
      vector<char> PCR(pp/3);
      
      
      Translate(pp,CR,PCR);//translates coding regions of reference genome

      // Prot<<endl<<"Coding region: "<<i+1<<endl;
     
      //for(int h=0;h<pp/3;h++){Prot<<" "<<h+1<<PCR[h];}

      
      Protein.push_back(PCR);
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// READ VARIABILITY DATA AND CALCULATE CONSENSUS SEQUENCE /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  vector<char> ConsensusSequence(NUMBER_OF_BASES);  

  char CharToSkip;
  int IntToSkip;
  
  int fpr=End_5Prime-Start_5Prime+1;
  int tpr=End_3Prime-Start_3Prime;//because End_3Prime=NumberOfBases, no need to add 1 (which was added to include the 0 element)

  vector<double> FPrA(fpr);//here will keep frequencies of occurrence of nt A, C, G, T in FivePrime region as  observed in gene sequence databases 
  vector<double> FPrC(fpr);
  vector<double> FPrG(fpr);
  vector<double> FPrT(fpr);
  
  vector<double> FPrConsensus(fpr);
  vector<double> FPrReference(fpr);//frequency of nt at position of Reference sequence in 5prime region
  vector<double> FPrMutant; //read mutants from a given file, in this case JMutants.txt

  vector<double> TPrA(tpr);//here will keep frequencies of occurrence of nt A, C, G, T in ThreePrime region as  observed in gene sequence databases 
  vector<double> TPrC(tpr);
  vector<double> TPrG(tpr);
  vector<double> TPrT(tpr);
  
  vector<double> TPrConsensus(tpr);
  vector<double> TPrReference(tpr);//frequency of nt at position of Reference sequence in 3prime region
  vector<double> TPrMutant(tpr);

  double NN;
  double Total;

  ifstream inGENESVdata("GeneSV_heatMAP_codons.matrix");
  if(!inGENESVdata){cout<<"NO FILE"<<endl; exit(1);}
 
  int Number0fpr=0;int Number0tpr=0;int Number0STR=0;int Number0NSTR=0;
  
  inGENESVdata.ignore(1000,'\n');
  inGENESVdata.ignore(1000,'\n');
  inGENESVdata.ignore(1000,'\n');
  inGENESVdata.ignore(1000,'\n');
  
  for(int i=0; i<fpr; i++) //at each position of 5 Prime region calculate frequency of occurrence observed in gene sequence databases of nt A, C, G, T 
    {
      inGENESVdata>>CharToSkip>>IntToSkip>>FPrA[i]>>FPrC[i]>>FPrG[i]>>FPrT[i]>>NN>>Total>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip;
      
      inGENESVdata.ignore(1000,'\n');
      
      FPrA[i]/=(Total-NN); 
      FPrConsensus[i]=FPrA[i]; if(!FPrA[i]){Number0fpr++;}

      FPrC[i]/=(Total-NN); 
      FPrConsensus[i]=max(FPrConsensus[i],FPrC[i]);if(!FPrC[i]){Number0fpr++;}

      FPrG[i]/=(Total-NN); 
      FPrConsensus[i]=max(FPrConsensus[i],FPrG[i]); if(!FPrG[i]){Number0fpr++;}

      FPrT[i]/=(Total-NN); 
      FPrConsensus[i]=max(FPrConsensus[i],FPrT[i]); if(!FPrT[i]){Number0fpr++;}
  
      //cout<<i<<":  " <<FPrA[i]<<" "<<FPrC[i]<<" "<<FPrG[i]<<" "<<FPrT[i]<<" "<<endl;
    }
  
  
  //this version works for viruses encoding polyproteins 

  int sstrn=(End_Str_Prot_Region-Start_Str_Prot_Region+1)/3; // number of codons in str region 
  int nstrn=(End_NStr_Prot_Region-Start_NStr_Prot_Region+1)/3; // number of codons in nstr region
  vector<map<string,double> > STRCodonFreqs(sstrn); 
  vector<map<string,double> > NSTRCodonFreqs(nstrn);
 
  
  double DoubleToRead;

  vector<double> STRConsensusFreq(sstrn);//frequency of codon at position of Consensus sequence in Structural region
  vector<double> NSTRConsensusFreq(nstrn);//frequency of codon at position of Consensus sequence in NonStructural region

  vector<double> STRReferenceFreq(sstrn);//frequency of codon at position of Reference sequence in Structural region
  vector<double> STRReferenceAAFreq(sstrn);//frequency of codon at position of Reference sequence in NonStructural region

  vector<double> NSTRReferenceFreq(nstrn);
  vector<double> NSTRReferenceAAFreq(nstrn);

  vector<double> STRMutantFreq(sstrn);//frequency of codon at position of Mutant sequence in Structural region
  vector<double> STRMutantAAFreq(sstrn);//frequency of codon at position of Mutant sequence in NonStructural region

  vector<double> NSTRMutantFreq(nstrn);
  vector<double> NSTRMutantAAFreq(nstrn);


  vector<char> STRConsensusCodon_NT1(sstrn);//first nt of consensus codon in structural protein 
  vector<char> STRConsensusCodon_NT2(sstrn);//2nd nt of consensus codon in structural protein
  vector<char> STRConsensusCodon_NT3(sstrn);//3rd nt of consensus codon in structural proteins 
  vector<double> STRFrequencyOfConsensusCodon(sstrn);

  vector<char> NSTRConsensusCodon_NT1(nstrn);//first nt of consensus codon in non-structural protein 
  vector<char> NSTRConsensusCodon_NT2(nstrn);//2nd nt of consensus codon in non-structural protein
  vector<char> NSTRConsensusCodon_NT3(nstrn);//3rd nt of consensus codon in non-structural proteins 
  vector<double> NSTRFrequencyOfConsensusCodon(nstrn);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //HERE we find CONSENSUS codon at each position of Structural proteins region/////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i=0; i<sstrn; i++)
    {
      inGENESVdata.ignore(1000,'\n'); 

      double aux=0; double auxold=0;
      string straux="NNN";

      //MUCH OF WHAT IS BELOW IS UNNECESSARY REMNANT FROM A PREVIOUS VERSION; WE USE ONLY THE INFORMATION WHETHER A FREQUENCY IS 0 OR NOT

      inGENESVdata>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip;

      inGENESVdata>>STRCodonFreqs[i]["ttt"];STRCodonFreqs[i]["ttt"]/=100; if(!STRCodonFreqs[i]["ttt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ttt"]);
      if(aux>auxold){straux="ttt"; auxold=aux;}
      
      inGENESVdata>>STRCodonFreqs[i]["ttc"];STRCodonFreqs[i]["ttc"]/=100;if(!STRCodonFreqs[i]["ttc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ttc"]);
      if(aux>auxold){straux="ttc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tta"];STRCodonFreqs[i]["tta"]/=100;if(!STRCodonFreqs[i]["tta"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tta"]);
      if(aux>auxold){straux="tta"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ttg"];STRCodonFreqs[i]["ttg"]/=100;if(!STRCodonFreqs[i]["ttg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ttg"]);
      if(aux>auxold){straux="ttg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ctt"];STRCodonFreqs[i]["ctt"]/=100;if(!STRCodonFreqs[i]["ctt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ctt"]);
      if(aux>auxold){straux="ctt"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ctc"];STRCodonFreqs[i]["ctc"]/=100;if(!STRCodonFreqs[i]["ctc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ctc"]);
      if(aux>auxold){straux="ctc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cta"];STRCodonFreqs[i]["cta"]/=100;if(!STRCodonFreqs[i]["cta"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cta"]);
      if(aux>auxold){straux="cta"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ctg"];STRCodonFreqs[i]["ctg"]/=100;if(!STRCodonFreqs[i]["ctg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ctg"]);
      if(aux>auxold){straux="ctg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["att"];STRCodonFreqs[i]["att"]/=100;if(!STRCodonFreqs[i]["att"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["att"]);
      if(aux>auxold){straux="att"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["atc"];STRCodonFreqs[i]["atc"]/=100;if(!STRCodonFreqs[i]["atc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["atc"]);
      if(aux>auxold){straux="atc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ata"];STRCodonFreqs[i]["ata"]/=100;if(!STRCodonFreqs[i]["ata"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ata"]);
      if(aux>auxold){straux="ata"; auxold=aux;}
      
      inGENESVdata>>STRCodonFreqs[i]["atg"];STRCodonFreqs[i]["atg"]/=100;if(!STRCodonFreqs[i]["atg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["atg"]);
      if(aux>auxold){straux="atg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gtt"];STRCodonFreqs[i]["gtt"]/=100;if(!STRCodonFreqs[i]["gtt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gtt"]);
      if(aux>auxold){straux="gtt"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gtc"];STRCodonFreqs[i]["gtc"]/=100;if(!STRCodonFreqs[i]["gtc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gtc"]);
      if(aux>auxold){straux="gtc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gta"];STRCodonFreqs[i]["gta"]/=100;if(!STRCodonFreqs[i]["gta"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gta"]);
      if(aux>auxold){straux="gta"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gtg"];STRCodonFreqs[i]["gtg"]/=100;if(!STRCodonFreqs[i]["gtg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gtg"]);
      if(aux>auxold){straux="gtg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tct"];STRCodonFreqs[i]["tct"]/=100;if(!STRCodonFreqs[i]["tct"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tct"]);
      if(aux>auxold){straux="tct"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tcc"];STRCodonFreqs[i]["tcc"]/=100;if(!STRCodonFreqs[i]["tcc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tcc"]);
      if(aux>auxold){straux="tcc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tca"];STRCodonFreqs[i]["tca"]/=100;if(!STRCodonFreqs[i]["tca"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tca"]);
      if(aux>auxold){straux="tca"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tcg"];STRCodonFreqs[i]["tcg"]/=100;if(!STRCodonFreqs[i]["tcg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tcg"]);
      if(aux>auxold){straux="tcg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cct"];STRCodonFreqs[i]["cct"]/=100;if(!STRCodonFreqs[i]["cct"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cct"]);
      if(aux>auxold){straux="cct"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ccc"];STRCodonFreqs[i]["ccc"]/=100;if(!STRCodonFreqs[i]["ccc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ccc"]);
      if(aux>auxold){straux="ccc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cca"];STRCodonFreqs[i]["cca"]/=100;if(!STRCodonFreqs[i]["cca"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cca"]);
      if(aux>auxold){straux="cca"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ccg"];STRCodonFreqs[i]["ccg"]/=100;if(!STRCodonFreqs[i]["ccg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ccg"]);
      if(aux>auxold){straux="ccg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["act"];STRCodonFreqs[i]["act"]/=100;if(!STRCodonFreqs[i]["act"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["act"]);
      if(aux>auxold){straux="act"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["acc"];STRCodonFreqs[i]["acc"]/=100;if(!STRCodonFreqs[i]["acc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["acc"]);
      if(aux>auxold){straux="acc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aca"];STRCodonFreqs[i]["aca"]/=100;if(!STRCodonFreqs[i]["aca"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aca"]);
      if(aux>auxold){straux="aca"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["acg"];STRCodonFreqs[i]["acg"]/=100;if(!STRCodonFreqs[i]["acg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["acg"]);
      if(aux>auxold){straux="acg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gct"];STRCodonFreqs[i]["gct"]/=100;if(!STRCodonFreqs[i]["gct"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gct"]);
      if(aux>auxold){straux="gct"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gcc"];STRCodonFreqs[i]["gcc"]/=100;if(!STRCodonFreqs[i]["gcc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gcc"]);
      if(aux>auxold){straux="gcc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gca"];STRCodonFreqs[i]["gca"]/=100;if(!STRCodonFreqs[i]["gca"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gca"]);
      if(aux>auxold){straux="gca"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gcg"];STRCodonFreqs[i]["gcg"]/=100;if(!STRCodonFreqs[i]["gcg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gcg"]);
      if(aux>auxold){straux="gcg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tat"];STRCodonFreqs[i]["tat"]/=100;if(!STRCodonFreqs[i]["tat"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tat"]);
      if(aux>auxold){straux="tat"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tac"];STRCodonFreqs[i]["tac"]/=100;if(!STRCodonFreqs[i]["tac"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tac"]);
      if(aux>auxold){straux="tac"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["taa"];STRCodonFreqs[i]["taa"]/=100;if(!STRCodonFreqs[i]["taa"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["taa"]);
      if(aux>auxold){straux="taa"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tag"];STRCodonFreqs[i]["tag"]/=100;if(!STRCodonFreqs[i]["tag"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tag"]);
      if(aux>auxold){straux="tag"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cat"];STRCodonFreqs[i]["cat"]/=100;if(!STRCodonFreqs[i]["cat"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cat"]);
      if(aux>auxold){straux="cat"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cac"];STRCodonFreqs[i]["cac"]/=100;if(!STRCodonFreqs[i]["cac"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cac"]);
      if(aux>auxold){straux="cac"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["caa"];STRCodonFreqs[i]["caa"]/=100;if(!STRCodonFreqs[i]["caa"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["caa"]);
      if(aux>auxold){straux="caa"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cag"];STRCodonFreqs[i]["cag"]/=100;if(!STRCodonFreqs[i]["cag"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cag"]);
      if(aux>auxold){straux="cag"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aat"];STRCodonFreqs[i]["aat"]/=100;if(!STRCodonFreqs[i]["aat"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aat"]);
      if(aux>auxold){straux="aat"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aac"];STRCodonFreqs[i]["aac"]/=100;if(!STRCodonFreqs[i]["aac"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aac"]);
      if(aux>auxold){straux="aac"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aaa"];STRCodonFreqs[i]["aaa"]/=100;if(!STRCodonFreqs[i]["aaa"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aaa"]);
      if(aux>auxold){straux="aaa"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aag"];STRCodonFreqs[i]["aag"]/=100;if(!STRCodonFreqs[i]["aag"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aag"]);
      if(aux>auxold){straux="aag"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gat"];STRCodonFreqs[i]["gat"]/=100;if(!STRCodonFreqs[i]["gat"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gat"]);
      if(aux>auxold){straux="gat"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gac"];STRCodonFreqs[i]["gac"]/=100;if(!STRCodonFreqs[i]["gac"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gac"]);
      if(aux>auxold){straux="gac"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gaa"];STRCodonFreqs[i]["gaa"]/=100;if(!STRCodonFreqs[i]["gaa"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gaa"]);
      if(aux>auxold){straux="gaa"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gag"];STRCodonFreqs[i]["gag"]/=100;if(!STRCodonFreqs[i]["gag"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gag"]);
      if(aux>auxold){straux="gag"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tgt"];STRCodonFreqs[i]["tgt"]/=100;if(!STRCodonFreqs[i]["tgt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tgt"]);
      if(aux>auxold){straux="tgt"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["tgc"];STRCodonFreqs[i]["tgc"]/=100;if(!STRCodonFreqs[i]["tgc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tgc"]);
      if(aux>auxold){straux="tgc"; auxold=aux;}


      inGENESVdata>>STRCodonFreqs[i]["tga"];STRCodonFreqs[i]["tga"]/=100;if(!STRCodonFreqs[i]["tga"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tga"]);
      if(aux>auxold){straux="tga"; auxold=aux;}


      inGENESVdata>>STRCodonFreqs[i]["tgg"];STRCodonFreqs[i]["tgg"]/=100;if(!STRCodonFreqs[i]["tgg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["tgg"]);
      if(aux>auxold){straux="tgg"; auxold=aux;}


      inGENESVdata>>STRCodonFreqs[i]["cgt"];STRCodonFreqs[i]["cgt"]/=100;if(!STRCodonFreqs[i]["cgt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cgt"]);
      if(aux>auxold){straux="cgt"; auxold=aux;}


      inGENESVdata>>STRCodonFreqs[i]["cgc"];STRCodonFreqs[i]["cgc"]/=100;if(!STRCodonFreqs[i]["cgc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cgc"]);
      if(aux>auxold){straux="cgc"; auxold=aux;}


      inGENESVdata>>STRCodonFreqs[i]["cga"];STRCodonFreqs[i]["cga"]/=100;if(!STRCodonFreqs[i]["cga"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cga"]);
      if(aux>auxold){straux="cga"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["cgg"];STRCodonFreqs[i]["cgg"]/=100;if(!STRCodonFreqs[i]["cgg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["cgg"]);
      if(aux>auxold){straux="cgg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["agt"];STRCodonFreqs[i]["agt"]/=100;if(!STRCodonFreqs[i]["agt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["agt"]);
      if(aux>auxold){straux="agt"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["agc"];STRCodonFreqs[i]["agc"]/=100;if(!STRCodonFreqs[i]["agc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["agc"]);
      if(aux>auxold){straux="agc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["aga"];STRCodonFreqs[i]["aga"]/=100;if(!STRCodonFreqs[i]["aga"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["aga"]);
      if(aux>auxold){straux="aga"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["agg"];STRCodonFreqs[i]["agg"]/=100;if(!STRCodonFreqs[i]["agg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["agg"]);
      if(aux>auxold){straux="agg"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ggt"];STRCodonFreqs[i]["ggt"]/=100;if(!STRCodonFreqs[i]["ggt"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ggt"]);
      if(aux>auxold){straux="ggt"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ggc"];STRCodonFreqs[i]["ggc"]/=100;if(!STRCodonFreqs[i]["ggc"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ggc"]);
      if(aux>auxold){straux="ggc"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["gga"];STRCodonFreqs[i]["gga"]/=100;if(!STRCodonFreqs[i]["gga"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["gga"]);
      if(aux>auxold){straux="gga"; auxold=aux;}

      inGENESVdata>>STRCodonFreqs[i]["ggg"];STRCodonFreqs[i]["ggg"]/=100;if(!STRCodonFreqs[i]["ggg"]){Number0STR++;};
      aux=max(aux,STRCodonFreqs[i]["ggg"]);
      if(aux>auxold){straux="ggg"; auxold=aux;}


      inGENESVdata.ignore(1000,'\n');
      inGENESVdata.ignore(1000,'\n');

    

      //cout<<"     "<<i<<"  Consensus codon "<<straux<< auxold<<"   ";

      STRConsensusCodon_NT1[i]=straux[0];
      STRConsensusCodon_NT2[i]=straux[1];
      STRConsensusCodon_NT3[i]=straux[2];

      STRFrequencyOfConsensusCodon[i]=auxold;//CONSENSUS CODONS NOT USED HERE...

      //cout<<" STRFrequencyNSTRCodonFreqs[i]["ggg"];OfConsensusCodon"<<STRFrequencyOfConsensusCodon[i]<<endl;

    }

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //HERE find consensus codon at each position of NON Structural proteins region so that to calculate distance from mutant to consensus/////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  for(int i=0; i<nstrn; i++)
    {
      double aux=0; double auxold=0;
      string straux="NNN";

      inGENESVdata.ignore(1000,'\n');
      inGENESVdata>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip;

      inGENESVdata>>NSTRCodonFreqs[i]["ttt"]; NSTRCodonFreqs[i]["ttt"]/=100;if(!NSTRCodonFreqs[i]["ttt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ttt"]);
      if(aux>auxold){straux="ttt"; auxold=aux;}
      
      inGENESVdata>>NSTRCodonFreqs[i]["ttc"];NSTRCodonFreqs[i]["ttc"]/=100;if(!NSTRCodonFreqs[i]["ttc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ttc"]); 
      if(aux>auxold){straux="ttc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tta"]; NSTRCodonFreqs[i]["tta"]/=100;if(!NSTRCodonFreqs[i]["tta"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tta"]);
      if(aux>auxold){straux="tta"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ttg"]; NSTRCodonFreqs[i]["ttg"]/=100;if(!NSTRCodonFreqs[i]["ttg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ttg"]);
      if(aux>auxold){straux="ttg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ctt"]; NSTRCodonFreqs[i]["ctt"]/=100;if(!NSTRCodonFreqs[i]["ctt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ctt"]);
      if(aux>auxold){straux="ctt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ctc"];  NSTRCodonFreqs[i]["ctc"]/=100;if(!NSTRCodonFreqs[i]["ctc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ctc"]);
      if(aux>auxold){straux="ctc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cta"]; NSTRCodonFreqs[i]["cta"]/=100;if(!NSTRCodonFreqs[i]["cta"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cta"]);
      if(aux>auxold){straux="cta"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ctg"];  NSTRCodonFreqs[i]["ctg"]/=100;if(!NSTRCodonFreqs[i]["ctg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ctg"]);
      if(aux>auxold){straux="ctg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["att"]; NSTRCodonFreqs[i]["att"]/=100;if(!NSTRCodonFreqs[i]["att"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["att"]);
      if(aux>auxold){straux="att"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["atc"];NSTRCodonFreqs[i]["atc"]/=100;if(!NSTRCodonFreqs[i]["atc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["atc"]);
      if(aux>auxold){straux="atc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ata"]; NSTRCodonFreqs[i]["ata"]/=100;if(!NSTRCodonFreqs[i]["ata"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ata"]);
      if(aux>auxold){straux="ata"; auxold=aux;}
      
      inGENESVdata>>NSTRCodonFreqs[i]["atg"]; NSTRCodonFreqs[i]["atg"]/=100;if(!NSTRCodonFreqs[i]["atg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["atg"]);
      if(aux>auxold){straux="atg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gtt"];  NSTRCodonFreqs[i]["gtt"]/=100;if(!NSTRCodonFreqs[i]["gtt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gtt"]);
      if(aux>auxold){straux="gtt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gtc"];  NSTRCodonFreqs[i]["gtc"]/=100;if(!NSTRCodonFreqs[i]["gtc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gtc"]);
      if(aux>auxold){straux="gtc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gta"];NSTRCodonFreqs[i]["gta"]/=100;if(!NSTRCodonFreqs[i]["gta"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gta"]);
      if(aux>auxold){straux="gta"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gtg"];NSTRCodonFreqs[i]["gtg"]/=100;if(!NSTRCodonFreqs[i]["gtg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gtg"]);
      if(aux>auxold){straux="gtg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tct"]; NSTRCodonFreqs[i]["tct"]/=100;if(!NSTRCodonFreqs[i]["tct"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tct"]);
      if(aux>auxold){straux="tct"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tcc"]; NSTRCodonFreqs[i]["tcc"]/=100;if(!NSTRCodonFreqs[i]["tcc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tcc"]);
      if(aux>auxold){straux="tcc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tca"]; NSTRCodonFreqs[i]["tca"]/=100;if(!NSTRCodonFreqs[i]["tca"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tca"]);
      if(aux>auxold){straux="tca"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tcg"]; NSTRCodonFreqs[i]["tcg"]/=100;if(!NSTRCodonFreqs[i]["tcg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tcg"]);
      if(aux>auxold){straux="tcg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cct"];  NSTRCodonFreqs[i]["cct"]/=100;if(!NSTRCodonFreqs[i]["cct"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cct"]);
      if(aux>auxold){straux="cct"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ccc"];  NSTRCodonFreqs[i]["ccc"]/=100;if(!NSTRCodonFreqs[i]["ccc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ccc"]);
      if(aux>auxold){straux="ccc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cca"];  NSTRCodonFreqs[i]["cca"]/=100;if(!NSTRCodonFreqs[i]["cca"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cca"]);
      if(aux>auxold){straux="cca"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ccg"];NSTRCodonFreqs[i]["ccg"]/=100;if(!NSTRCodonFreqs[i]["ccg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ccg"]);
      if(aux>auxold){straux="ccg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["act"]; NSTRCodonFreqs[i]["act"]/=100;if(!NSTRCodonFreqs[i]["act"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["act"]);
      if(aux>auxold){straux="act"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["acc"];  NSTRCodonFreqs[i]["acc"]/=100;if(!NSTRCodonFreqs[i]["acc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["acc"]);
      if(aux>auxold){straux="acc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aca"];  NSTRCodonFreqs[i]["aca"]/=100;if(!NSTRCodonFreqs[i]["aca"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aca"]);
      if(aux>auxold){straux="aca"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["acg"];  NSTRCodonFreqs[i]["acg"]/=100;if(!NSTRCodonFreqs[i]["acg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["acg"]);
      if(aux>auxold){straux="acg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gct"]; NSTRCodonFreqs[i]["gct"]/=100;if(!NSTRCodonFreqs[i]["gct"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gct"]);
      if(aux>auxold){straux="gct"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gcc"];  NSTRCodonFreqs[i]["gcc"]/=100;if(!NSTRCodonFreqs[i]["gcc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gcc"]);
      if(aux>auxold){straux="gcc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gca"]; NSTRCodonFreqs[i]["gca"]/=100;if(!NSTRCodonFreqs[i]["gca"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gca"]);
      if(aux>auxold){straux="gca"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gcg"];NSTRCodonFreqs[i]["gcg"]/=100;if(!NSTRCodonFreqs[i]["gcg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gcg"]);
      if(aux>auxold){straux="gcg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tat"];NSTRCodonFreqs[i]["tat"]/=100;if(!NSTRCodonFreqs[i]["tat"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tat"]);
      if(aux>auxold){straux="tat"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tac"]; NSTRCodonFreqs[i]["tac"]/=100;if(!NSTRCodonFreqs[i]["tac"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tac"]);
      if(aux>auxold){straux="tac"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["taa"]; NSTRCodonFreqs[i]["taa"]/=100;if(!NSTRCodonFreqs[i]["taa"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["taa"]);
      if(aux>auxold){straux="taa"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tag"];NSTRCodonFreqs[i]["tag"]/=100;if(!NSTRCodonFreqs[i]["tag"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tag"]);
      if(aux>auxold){straux="tag"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cat"];NSTRCodonFreqs[i]["cat"]/=100;if(!NSTRCodonFreqs[i]["cat"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cat"]);
      if(aux>auxold){straux="cat"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cac"];NSTRCodonFreqs[i]["cac"]=100;if(!NSTRCodonFreqs[i]["cac"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cac"]);
      if(aux>auxold){straux="cac"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["caa"];NSTRCodonFreqs[i]["caa"]/=100;if(!NSTRCodonFreqs[i]["caa"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["caa"]);
      if(aux>auxold){straux="caa"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cag"];NSTRCodonFreqs[i]["cag"]/=100;if(!NSTRCodonFreqs[i]["cag"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cag"]);
      if(aux>auxold){straux="cag"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aat"];NSTRCodonFreqs[i]["aat"]/=100;if(!NSTRCodonFreqs[i]["aat"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aat"]);
      if(aux>auxold){straux="aat"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aac"];NSTRCodonFreqs[i]["aac"]/=100;if(!NSTRCodonFreqs[i]["aac"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aac"]);
      if(aux>auxold){straux="aac"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aaa"];NSTRCodonFreqs[i]["aaa"]/=100;if(!NSTRCodonFreqs[i]["aaa"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aaa"]);
      if(aux>auxold){straux="aaa"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aag"];NSTRCodonFreqs[i]["aag"]/=100;if(!NSTRCodonFreqs[i]["aag"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aag"]);
      if(aux>auxold){straux="aag"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gat"];NSTRCodonFreqs[i]["gat"]/=100;if(!NSTRCodonFreqs[i]["gat"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gat"]);
      if(aux>auxold){straux="gat"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gac"];NSTRCodonFreqs[i]["gac"]/=100;if(!NSTRCodonFreqs[i]["gac"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gac"]);
      if(aux>auxold){straux="gac"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gaa"];NSTRCodonFreqs[i]["gaa"]/=100;if(!NSTRCodonFreqs[i]["gaa"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gaa"]);
      if(aux>auxold){straux="gaa"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gag"];NSTRCodonFreqs[i]["gag"]/=100;if(!NSTRCodonFreqs[i]["gag"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gag"]);
      if(aux>auxold){straux="gag"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tgt"];NSTRCodonFreqs[i]["tgt"]/=100;if(!NSTRCodonFreqs[i]["tgt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tgt"]);
      if(aux>auxold){straux="tgt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tgc"];NSTRCodonFreqs[i]["tgc"]/=100;if(!NSTRCodonFreqs[i]["tgc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tgc"]);
      if(aux>auxold){straux="tgc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tga"];NSTRCodonFreqs[i]["tga"]/=100;if(!NSTRCodonFreqs[i]["tga"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tga"]);
      if(aux>auxold){straux="tga"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["tgg"];NSTRCodonFreqs[i]["tgg"]/=100;if(!NSTRCodonFreqs[i]["tgg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["tgg"]);
      if(aux>auxold){straux="tgg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cgt"];NSTRCodonFreqs[i]["cgt"]/=100;if(!NSTRCodonFreqs[i]["cgt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cgt"]);
      if(aux>auxold){straux="cgt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cgc"];NSTRCodonFreqs[i]["cgc"]/=100;if(!NSTRCodonFreqs[i]["cgc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cgc"]);
      if(aux>auxold){straux="cgc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cga"];NSTRCodonFreqs[i]["cga"]/=100;if(!NSTRCodonFreqs[i]["cga"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cga"]);
      if(aux>auxold){straux="cga"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["cgg"];NSTRCodonFreqs[i]["cgg"]/=100;if(!NSTRCodonFreqs[i]["cgg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["cgg"]);
      if(aux>auxold){straux="cgg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["agt"];NSTRCodonFreqs[i]["agt"]/=100;if(!NSTRCodonFreqs[i]["agt"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["agt"]);
      if(aux>auxold){straux="agt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["agc"];NSTRCodonFreqs[i]["agc"]/=100;if(!NSTRCodonFreqs[i]["agc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["agc"]);
      if(aux>auxold){straux="agc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["aga"];NSTRCodonFreqs[i]["aga"]/=100;if(!NSTRCodonFreqs[i]["aga"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["aga"]);
      if(aux>auxold){straux="aga"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["agg"];NSTRCodonFreqs[i]["agg"]/=100;if(!NSTRCodonFreqs[i]["agg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["agg"]);
      if(aux>auxold){straux="agg"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ggt"];NSTRCodonFreqs[i]["ggt"]/=100;if(!NSTRCodonFreqs[i]["ttc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ggt"]);
      if(aux>auxold){straux="ggt"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ggc"];NSTRCodonFreqs[i]["ggc"]/=100;if(!NSTRCodonFreqs[i]["ggc"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ggc"]);
      if(aux>auxold){straux="ggc"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["gga"];NSTRCodonFreqs[i]["gga"]/=100;if(!NSTRCodonFreqs[i]["gga"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["gga"]);
      if(aux>auxold){straux="gga"; auxold=aux;}

      inGENESVdata>>NSTRCodonFreqs[i]["ggg"];NSTRCodonFreqs[i]["ggg"]/=100;if(!NSTRCodonFreqs[i]["ggg"]){Number0NSTR++;};
      aux=max(aux,NSTRCodonFreqs[i]["ggg"]);
      if(aux>auxold){straux="ggg"; auxold=aux;}

      inGENESVdata.ignore(1000,'\n');
      inGENESVdata.ignore(1000,'\n');
   
    
      //cout<<"    "<<i<<"  Consensus codon "<<straux<< auxold<<"   ";


      NSTRConsensusCodon_NT1[i]=straux[0];//these are the nucleotides of the consensus sequence in the nonstr protein
      NSTRConsensusCodon_NT2[i]=straux[1];
      NSTRConsensusCodon_NT3[i]=straux[2];

      NSTRFrequencyOfConsensusCodon[i]=auxold;
 }

 
  
  for(int i=0; i<tpr; i++)//at each position of 3 Prime region calculate frequency of occurrence observed in gene sequence databases of nt A, C, G, T 
    {
      inGENESVdata>>CharToSkip>>IntToSkip>>TPrA[i]>>TPrC[i]>>TPrG[i]>>TPrT[i]>>NN>>Total>>CharToSkip>>IntToSkip>>IntToSkip>>IntToSkip>>CharToSkip;
      
      inGENESVdata.ignore(1000,'\n');
      
      TPrA[i]/=(Total-NN); 
      TPrConsensus[i]=TPrA[i]; if(!TPrA[i]){Number0tpr++;}

      TPrC[i]/=(Total-NN); 
      TPrConsensus[i]=max(TPrConsensus[i],TPrC[i]);if(!TPrC[i]){Number0tpr++;}

      TPrG[i]/=(Total-NN); 
      TPrConsensus[i]=max(TPrConsensus[i],TPrG[i]);if(!TPrG[i]){Number0tpr++;}


      TPrT[i]/=(Total-NN);
      TPrConsensus[i]=max(TPrConsensus[i],TPrT[i]);    if(!TPrT[i]){Number0tpr++;}
      //cout<<"i: "<<i<<TPrA[i]<<" "<<TPrC[i]<<" "<<TPrG[i]<<" "<<TPrT[i]<<" "<<endl;
    }
   
  //cout<<"end read variability info"<<endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////End reading variability Info/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  ////////////////////////////////////////  pick  positions from reference sequence to define highest fitness positions /////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 


  ifstream HFpositions("HFpositions.txt");
  if(! HFpositions){cout<<"NO FILE for HFpositions"<<endl; exit(1);}

  int NumberOfHF_5T;
  int NumberOfHF_3T;
  int NumberOfHF_5R;
  int NumberOfHF_3R; 
  int NumberOfHF_STR; 
  int NumberOfHF_NonSTR;
  string ToReadSkip;


/////////////////////////////////////////////////////////////////////////////////
  /////////////////////////READ HIGHEST FITNESS POSITIONS//////////////////////////
  ///////////////////////////////////////////////////////////////////////////////// 
 
 
  HFpositions>>ToReadSkip>>NumberOfHF_5T>>ToReadSkip;

  vector<int>vHF_5T(NumberOfHF_5T);
  vector<char>HF_5T(NumberOfHF_5T);
  if(NumberOfHF_5T)
    {
      for(int i=0; i<NumberOfHF_5T; i++)
	{
	  HFpositions>>vHF_5T[i]>>HF_5T[i]; 
	  //cout<<" "<<vHF_5T[i]<<"  "<<HF_5T[i]<<endl;
	}
    }
  
  HFpositions>>ToReadSkip>>NumberOfHF_5R>>ToReadSkip;
  vector<int>vHF_5R(NumberOfHF_5R);
  vector<char>HF_5R(NumberOfHF_5R);

  if(NumberOfHF_5R)
    {
      for(int i=0; i<NumberOfHF_5R; i++)
	{
	  HFpositions>>vHF_5R[i]>>HF_5R[i]; 
	  //cout<<" "<<vHF_5R[i]<<"  "<<HF_5R[i]<<endl;
	}
    }

  HFpositions>>ToReadSkip>>NumberOfHF_STR>>ToReadSkip; //NumberOfHF_STR should be a multiple of 3
  
  if(NumberOfHF_STR%3){cout<<'codon not complete in input file HFpositions.txt'; exit(1);}

  vector<int>vHF_STR(NumberOfHF_STR); 
  vector<char>HF_STR(NumberOfHF_STR);
  
  if(NumberOfHF_STR)
    {
      for(int i=0; i<NumberOfHF_STR; i++)
	{
	  HFpositions>>vHF_STR[i]>>HF_STR[i]; 
	  //cout<<" "<<vHF_STR[i]<<"  "<<HF_STR[i]<<endl;
	}
    }

  HFpositions>>ToReadSkip>>NumberOfHF_NonSTR>>ToReadSkip;//NumberOfHF_NonSTR should be a multiple of 3
  if(NumberOfHF_STR%3){cout<<'codon not complete in input file HFpositions.txt'; exit(1);}

  vector<int>vHF_NonSTR(NumberOfHF_NonSTR);
  vector<char>HF_NonSTR(NumberOfHF_NonSTR);
  
  if(NumberOfHF_NonSTR)
    {
      for(int i=0; i<NumberOfHF_NonSTR; i++)
	{
	  HFpositions>>vHF_NonSTR[i]>>HF_NonSTR[i];
	  //cout<<" "<<vHF_NonSTR[i]<<"  "<<HF_NonSTR[i]<<endl;
	}
    }

  HFpositions>>ToReadSkip>>NumberOfHF_3T>>ToReadSkip;
  
  vector<int>vHF_3T(NumberOfHF_3T); 
  vector<char>HF_3T(NumberOfHF_3T); 

  if(NumberOfHF_3T)
    {
      for(int i=0; i<NumberOfHF_3T; i++)
	{
	  HFpositions>>vHF_3T[i]>>HF_3T[i]; 
	  //cout<<" "<<vHF_3T[i]<<"  "<<HF_3T[i]<<endl;
	}
    }

   HFpositions>>ToReadSkip>>NumberOfHF_3R>>ToReadSkip;
  
  vector<int>vHF_3R(NumberOfHF_3R); 
  vector<char>HF_3R(NumberOfHF_3R); 
  
  if(NumberOfHF_3R)
    {
      for(int i=0; i<NumberOfHF_3R; i++)
	{
	  HFpositions>>vHF_3R[i]>>HF_3R[i]; 
	  //cout<<" "<<vHF_3R[i]<<"  "<<HF_3R[i]<<endl;
	}
    }

  // cout<<"end read HF positions"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////CALCULATE FITNESS OF REFERENCE SEQUENCE///////////////////////////////////////////////////
  ////////////////This module also calculates the number of positions at which the minimum  has been reached//////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  char A;
  double GG;
  char NT[3];

      
   double RSFitnessOfTranslation=0;
   double RSFitnessOfEntry=0;
   double RSFitnessOfReplication=0;

   //cout<<"begin calculate fitness of refseq"<<endl;
   
   CalculateFitnessOfFullSequence(NumberOfHF_5T,NumberOfHF_5R, vHF_5T, vHF_5R, HF_5T,HF_5R, NumberOfHF_3T,NumberOfHF_3R, vHF_3T,vHF_3R, HF_3T, HF_3R, NumberOfHF_STR, vHF_STR, HF_STR, NumberOfHF_NonSTR, vHF_NonSTR, HF_NonSTR,  ReferenceSequence, RSFitnessOfTranslation, RSFitnessOfEntry, RSFitnessOfReplication, Start_3Prime, FPrA, FPrC, FPrG, FPrT, TPrA,TPrC, TPrG, TPrT, fpr, tpr, sstrn,nstrn, STRCodonFreqs, NSTRCodonFreqs, Start_Str_Prot_Region,Start_NStr_Prot_Region);


   //cout<<"end calculate fitness of refseq"<<endl;

   cout<<"Reference sequence:"<<"  RSFitnessOfTranslation "<<RSFitnessOfTranslation<<" RSFitnessOfEntry: "<< RSFitnessOfEntry<<" RSFitnessOfReplication: "<<RSFitnessOfReplication<<endl;

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////////END CALCULATE FITNESS OF REFERENCE SEQUENCE////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int NumberOfInfectingGenotypes; //how many viruses will infect cell culture
  NumberOfInfectingGenotypes=floor(MOI*NumberOfCellsInCulture);
  //cout<<"NumberOfInfectingGenotypes "<<NumberOfInfectingGenotypes<<endl;

  int MaxNumberInitialMutations;//to create initial mutants from the reference 
  simInput>>ToRead>>MaxNumberInitialMutations;
  //cout<<"MaxNumberInitialMutations "<<MaxNumberInitialMutations<<endl;


  //********************************************************************************************** 
  Cloud QSCloud;
  QSCloud.InitializeCloud(NUMBER_OF_BASES, 0);

  //********************************************************************************************** 
  vector<char> Mutant(NUMBER_OF_BASES); 

  ofstream MutStat("MutationStats1.txt");
  ofstream Titer("Titer1.txt");
  ofstream TransStat("TransTVStats1.txt");
  ofstream Sequences("CellCloud1.txt");
  ofstream BurstTimes("BurstTimes1.txt");
  ofstream SeqShort("Dominants1.txt");
  ofstream outHM("HeatMap1.txt");
  ofstream outInoculum("Inoculum1.txt");
  ofstream FEDistr("FE_Distr.txt");
  ofstream FRDistr("FR_Distr.txt");
  ofstream FTDistr("FT_Distr.txt");
  ofstream FERDistr("FER_Distr.txt");
  ofstream BTDistr("BurstT_Distr.txt");
  ofstream CellsNonReplic("CellsNonreplVir.txt");
  ofstream Dynamics("Dynamics.txt");
  ifstream mutationsInput("InitialCloud1.txt");
  ifstream thisPassage("thisPassage1.txt");

  int CloudDimension=0;
  
  if(!NumberOfInfectingGenotypes){cout<<"Initial Inoculum=0";exit(1);}
  
  ofstream nextPassage("nextPassage1.txt");
  int NextPassage=0;   

  map<double,int> FEdistribution; 
  map<double,int> FRdistribution; 
  map<double,int> FTdistribution; 
  map<double,int> FERdistribution; 

  map<double,int>BurstTdistr;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// here 3 options for initial innoculum ///////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(InoculumInput==1)
    {
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////INPUT INITIAL INOCULUM FROM DATA////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      NextPassage=0;  
      int position=9;
      char PointMut;
      int r;

      for(int i=0;i<NumberOfInfectingGenotypes; i++) 
	{
	  Genotype G;
	  G.Initialize(NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation,RSFitnessOfEntry,RSFitnessOfReplication);
	  while(position!=0)//this is because in Jonathan's data, position is 0 only at end of file
	    {
	      mutationsInput>>position; //read mutants from a given file, in this case JMutants.txt

	      if(position)
		{
		  mutationsInput>>PointMut; 
		  G.ntMutations[position]=PointMut;
		}
	      else{}
	    }

	  r=G.ntMutations.size();
	  CloudDimension=max(CloudDimension, r); //calculate the max HD from the reference

	  G.PrintToScreen_nt_Mutations();

	  position=9;

	  for(int i=0; i<NUMBER_OF_BASES; i++){Mutant[i]=ReferenceSequence[i];}
	  
	  if(G.ntMutations.size())
	    { 
	      for(map<int,char>::const_iterator ps=G.ntMutations.begin(); ps!=G.ntMutations.end(); ps++)
		{
		  int tt=ps->first;
		  Mutant[tt]=ps->second; 	
		}//now Mutant has the restored sequence of the infecting genotype
	    }
	  else//the infecting sequence is the reference
	    {
	      //Sequences<<"infecting sequence has no mutations"<<endl;
	    }
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////CALCULATE FITNESS OF MUTANT Inoculum SEQUENCE////////////////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  double MutantFitnessOfTranslation=0;
	  double MutantFitnessOfEntry=0;
	  double MutantFitnessOfReplication=0;

	  cout<<"Calculate the fitness of inoculum members.... "<<endl;

	  CalculateFitnessOfFullSequence(NumberOfHF_5T,NumberOfHF_5R, vHF_5T, vHF_5R, HF_5T,HF_5R, NumberOfHF_3T,NumberOfHF_3R, vHF_3T,vHF_3R, HF_3T, HF_3R, NumberOfHF_STR, vHF_STR, HF_STR, NumberOfHF_NonSTR, vHF_NonSTR, HF_NonSTR,  Mutant,MutantFitnessOfTranslation,MutantFitnessOfEntry,MutantFitnessOfReplication, Start_3Prime, FPrA, FPrC, FPrG, FPrT, TPrA,TPrC, TPrG, TPrT, fpr, tpr, sstrn,nstrn, STRCodonFreqs, NSTRCodonFreqs, Start_Str_Prot_Region,Start_NStr_Prot_Region);
	  
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////END CALCULATE FITNESS OF Mutant Inoculum SEQUENCE////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  
	  G.FitnessOfTranslation= MutantFitnessOfTranslation;
	  G.FitnessOfReplication= MutantFitnessOfReplication;
	  G.FitnessOfEntry= MutantFitnessOfEntry;
	  
	  //cout<<"this genotype's fTr " <<G.FitnessOfTranslation<<"this genotype fRepl "<<G.FitnessOfReplication<<"this genotype fEntry "<<G.FitnessOfEntry<<endl;
	  
	  QSCloud.AddGenotypeToCloud(G);
	  
	  //G.Print_nt_Mutations(Sequences);
	  
	  
	}
    }
  
  if(InoculumInput==2) //Read inoculum from saved file from previous passage
    {
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////// INITIAL INOCULUM IS FROM PREVIOUS PASSAGE /////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      int position=9;
      int r;

      thisPassage>>NextPassage;
      cout<<"Passage" <<NextPassage<<endl;

      char PointMut;
      
      for(int i=0;i<NumberOfInfectingGenotypes; i++)  //create NumberOfInfectingGenotypes mutants 
	{
	  Genotype G;
	  G.Initialize(NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation,RSFitnessOfEntry,RSFitnessOfReplication);
	  while(position!=-9)
	    {
	      thisPassage>>position;
	      if(position)
		{
		  thisPassage>>PointMut; 
		  G.ntMutations[position]=PointMut;
		}
	      else{}
	    }
	  r=G.ntMutations.size();

	  CloudDimension=max(CloudDimension, r); //calculate the max HD from the reference

	  //G.PrintToScreen_nt_Mutations();	  
	  
	  position=9;
	  
	  for(int i=0; i<NUMBER_OF_BASES; i++){Mutant[i]=ReferenceSequence[i];}

	  if(G.ntMutations.size())
	    { 
	      for(map<int,char>::const_iterator ps=G.ntMutations.begin(); ps!=G.ntMutations.end(); ps++)
		{
		  int tt=ps->first;
		  Mutant[tt]=ps->second; 	
		}//now Infecting sequence has the restored sequence of the infecting genotype
	    }
	  else
	    {
	      //Sequences<<"infecting sequence has no mutations"<<endl;
	    }

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////CALCULATE FITNESS OF MUTANT Inoculum SEQUENCE////////////////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  double MutantFitnessOfTranslation=0;
	  double MutantFitnessOfEntry=0;
	  double MutantFitnessOfReplication=0;
	  
	  CalculateFitnessOfFullSequence(NumberOfHF_5T,NumberOfHF_5R, vHF_5T, vHF_5R, HF_5T,HF_5R, NumberOfHF_3T,NumberOfHF_3R, vHF_3T,vHF_3R, HF_3T, HF_3R, NumberOfHF_STR, vHF_STR, HF_STR, NumberOfHF_NonSTR, vHF_NonSTR, HF_NonSTR,  Mutant,MutantFitnessOfTranslation,MutantFitnessOfEntry,MutantFitnessOfReplication, Start_3Prime, FPrA, FPrC, FPrG, FPrT, TPrA,TPrC, TPrG, TPrT, fpr, tpr, sstrn,nstrn, STRCodonFreqs, NSTRCodonFreqs, Start_Str_Prot_Region,Start_NStr_Prot_Region);
	  
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////END CALCULATE FITNESS OF Mutant Inoculum SEQUENCE////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  
	  G.FitnessOfTranslation=MutantFitnessOfTranslation;
	  G.FitnessOfReplication=MutantFitnessOfReplication;
	  G.FitnessOfEntry=MutantFitnessOfEntry;
	  
	  //cout<<"this genotype's fTr " <<G.FitnessOfTranslation<<"this genotype fRepl "<<G.FitnessOfReplication<<"this genotype fEntry "<<G.FitnessOfEntry<<endl;
	  
	  QSCloud.AddGenotypeToCloud(G);
	  
	  //G.Print_nt_Mutations(Sequences);
	  
	  
	}
    }
  
  
  if(!InoculumInput)
    { 
      NextPassage=0;  
      for(int i=0;i<NumberOfInfectingGenotypes; i++)  //create NumberOfInfectingGenotypes mutants 
	{
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////GENERATE INITIAL INOCULUM AS random MUTANTS FROM REFERENCE SEQUENCE/////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  MTRand_closed hgrand; 
	  double ss= hgrand(); 
	  int r=floor(ss*MaxNumberInitialMutations); 
	  //cout<<endl<<"Number of mutations in this mutant: "<<r<<endl;//number of mutations in this mutant
	  
	  CloudDimension=max(CloudDimension, r); //calculate the max HD from the reference
	  
	  Genotype G;
	  G.Initialize(NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation,RSFitnessOfEntry,RSFitnessOfReplication);
	  G.CreateMutations(ReferenceSequence,NUMBER_OF_BASES,r,TransitionProb);
	  
	  //G.PrintToScreen_nt_Mutations(); 
	  //G.PrintTwinsToScreen();
	  
	  vector<char> Mutant(NUMBER_OF_BASES); 
	  for(int i=0; i<NUMBER_OF_BASES; i++){Mutant[i]=ReferenceSequence[i];}
	  
	  if(G.ntMutations.size())
	    { 
	      for(map<int,char>::const_iterator ps=G.ntMutations.begin(); ps!=G.ntMutations.end(); ps++)
		{
		  int tt=ps->first;
		  Mutant[tt]=ps->second; 	
		}//now Infecting sequence has the restored sequence of the infecting genotype
	    }
	  else
	    {
	      //Sequences<<"infecting sequence has no mutations"<<endl;
	    }

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////CALCULATE FITNESS OF MUTANT Inoculum SEQUENCE////////////////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  double MutantFitnessOfTranslation=0;
	  double MutantFitnessOfEntry=0;
	  double MutantFitnessOfReplication=0;
	  
	  CalculateFitnessOfFullSequence(NumberOfHF_5T,NumberOfHF_5R, vHF_5T, vHF_5R, HF_5T,HF_5R, NumberOfHF_3T,NumberOfHF_3R, vHF_3T,vHF_3R, HF_3T, HF_3R, NumberOfHF_STR, vHF_STR, HF_STR, NumberOfHF_NonSTR, vHF_NonSTR, HF_NonSTR,  Mutant,MutantFitnessOfTranslation,MutantFitnessOfEntry,MutantFitnessOfReplication, Start_3Prime, FPrA, FPrC, FPrG, FPrT, TPrA,TPrC, TPrG, TPrT, fpr, tpr, sstrn,nstrn, STRCodonFreqs, NSTRCodonFreqs, Start_Str_Prot_Region,Start_NStr_Prot_Region);
	   
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////END CALCULATE FITNESS OF Mutant Inoculum SEQUENCE////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  
	  G.FitnessOfTranslation=MutantFitnessOfTranslation;
	  G.FitnessOfReplication=MutantFitnessOfReplication;
	  G.FitnessOfEntry=MutantFitnessOfEntry;
	  
	  //cout<<"this genotype's fTr " <<G.FitnessOfTranslation<<"this genotype fRepl "<<G.FitnessOfReplication<<"this genotype fEntry "<<G.FitnessOfEntry<<endl;
	  
	  QSCloud.AddGenotypeToCloud(G);
	  
	  //G.Print_nt_Mutations(Sequences);
	  

	}
    }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////END GENERATE INITIAL INNOCULUM ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Initial inoculum generated"<<endl;//now printing its heatmap...

  outHM<<"Passage "<<NextPassage<<endl;
  QSCloud.PrintCloudHeatMap(outHM,ReferenceSequence);//prints a nucleotide frequency matrix of the inoculum
  
  QSCloud.PrintCloudHD(CloudDimension+1,outInoculum); 
  
  QSCloud.CalculateNumberOfMembersInCloud();//now that QSC cloud is populated with inoculum, count its members (to inform cloud member), the total count should be equal to the previously declared NIG

  NumberOfInfectingGenotypes=QSCloud.NumberOfMembers; // counts the number of sequences in the initial inoculum; this is actually done for testing - the  NumberOfInfectingGenotypes should not have changed
  //cout<<"NumberOfInfectingGenotypes "<<NumberOfInfectingGenotypes<<endl;

  CloudDimension=0;
  QSCloud.PrintCloudToFile(Sequences,SeqShort, NumberOfDominants, -1, CloudDimension, NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation, RSFitnessOfEntry, RSFitnessOfReplication);//Sequences and SeqShort are output files; this prints the dominants, everything else was suppressed because files beconme immense; returns CloudDimension being the maximum Hamming distance of the cloud; here -1 is the round of infection 
 
  MutStat<<"Initial Inoculum"<<endl;
  QSCloud.PrintCloudHD(CloudDimension,MutStat);  //not sure why this is printed twice    
  QSCloud.CalculateNumberOfMembersInCloud();

  cout<<"Number of members in cloud  "<<QSCloud.NumberOfMembers<<endl;
   
  Titer<<"passage 0  "; 
  Titer<<QSCloud.NumberOfMembers<<endl;
  
  outInoculum.close();
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////DEFINE PARAMETERS FOR SIMULATION/////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double uyy; 

  int MaxNumberInfectingVirions;
  simInput>>ToRead>>MaxNumberInfectingVirions;
  //cout<<"MaxNumberInfectingVirions "<<MaxNumberInfectingVirions<<endl;

  int NumberReplicationsPosToNeg;
  simInput>>ToRead>>NumberReplicationsPosToNeg;
  //cout<<"NumberReplicationsPosToNeg "<<NumberReplicationsPosToNeg<<endl;

  int NumberReplicationsNegToPos;
  simInput>>ToRead>>NumberReplicationsNegToPos;
  //cout<<"NumberReplicationsNegToPos "<<NumberReplicationsNegToPos<<endl;

  if(!NumberReplicationsPosToNeg){cout<<"NumberReplicationsPosToNeg=0!"; exit(1);}
  if(!NumberReplicationsNegToPos){cout<<"NumberReplicationsNegToPos=0!"; exit(1);}

  int BurstLoad;
  simInput>>ToRead>>BurstLoad;
  //cout<<"BurstLoad "<<BurstLoad<<endl;

  int CellR0;//this is usually a small number   
  simInput>>ToRead>>CellR0;
  //cout<<"CellR0 "<<CellR0<<endl;

  int  MaxNumberVirAttachedToCell;
  simInput>>ToRead>>MaxNumberVirAttachedToCell;
  //cout<<"MaxNumberVirAttachedToCell"<<MaxNumberVirAttachedToCell<<endl;

  int NumberHMtoPrint;
  simInput>>ToRead>>NumberHMtoPrint;

  int HeatMapsToPrint[NumberHMtoPrint];
  simInput>>ToRead;
  //cout<<"Heatmaps to Print ";
   
  for(int g=0; g<NumberHMtoPrint; g++)//Gives option to print heatmaps only for some passages; good for calculating large number of passages. 
    {
      simInput>>HeatMapsToPrint[g];
      //cout<<HeatMapsToPrint[g]<<" "; 
      if(HeatMapsToPrint[g]>=NumberOfPassages)
	{
	  cout<<"HeatMapsToPrint[g]>=NumberOfPassages!";exit(1);
	}
    }
  


  if(!NumberReplicationsPosToNeg){cout<<"NumberReplicationsPosToNeg=0!"; exit(1);}
  if(!NumberReplicationsNegToPos){cout<<"NumberReplicationsNegToPos=0!"; exit(1);}



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////*********************************************************************//////////////////////////////
  ////////////START LOOP FOR A LARGE NUMBER OF SIMULATIONS - CURRENTLY LARGE IS  = 1 NEED TO COMPLETE///////////////////////// 
  //            the problem with this is that often just 1 simulation is too long                    /////////////////////////
  ///////////////////////HERE ALL SIMULATIONS START WITH THE SAME INOCULUM GENERATED PREVIOUSLY///////////////////////////////
  ////////////////////////************************************************************************////////////////////////////
  ////////////////////////////////////////***********************************/////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 
  vector<Cell> NotInfectedCells;
  vector<Cell> InfectedCells; 
  vector<Cell> TempInfectedCells;
  vector<Cell> NewlyInfectedCells; 
  vector<Cell> InfectedNonproducingCells;

  Cell CC;

  int TotalTNum=0;//to count transitions
  int TotalVNum=0;//to count transversions


  double PassageSimTime=0;
  double TimeIncrement=0;//to progress time
  int indexvar;

  double Factor=0; 	
  double Msum=0;
  double Mhelper=0;
  double Mi=0;
  double r1;
  double rho1;
  int Lhelper=0;
  double CommonParam;
  int NRPN;//number replications from positive to negative
  int NRNP;//number replications from negative to positive
  char tt; 
  bool Flag=false;
  int kPM; 
  int kMP;
  
  int spomag;
  int TotalNumberBursts;

  for(int bx=0; bx<NumberOfSimulations; bx++)
    {  
      
      for(int cx = NextPassage; cx<NumberOfPassages; cx++)// simulation can start from a saved passage on
	{
	  TotalNumberBursts=0;
	  BurstTdistr.clear();

	  PassageSimTime=0;
	  cout<<"PASSAGE NUMBER: "<<cx+1<<endl;
	  Sequences<<endl<<"PASSAGE NUMBER: "<<cx+1<<endl;
	  SeqShort<<endl<<"PASSAGE NUMBER: "<<cx+1<<endl;
	  
	  //////////CREATE AN ARRAY OF NON-INFECTED CELLS//////////////////
	  
	  CC.Initialize();

	  NotInfectedCells.clear(); 
	  NewlyInfectedCells.clear();
	  InfectedNonproducingCells.clear();
	  TempInfectedCells.clear();
	  InfectedCells.clear();

	  for(int i=0; i<NumberOfCellsInCulture; i++)
	    {
	      NotInfectedCells.push_back(CC);
	      cout<<"NotInfectedCells.size()"<<NotInfectedCells.size();
	    }

	  InfectedCells.clear();
	  cout<<"InfectedCells.size()"<<InfectedCells.size();
	  
	  cout<<endl<<"INFECT FIRST ROUND OF CELLS WITH innoculum..."<<endl;
	  
	  int RoundOfInfection=0;
	  int RepeatInfection=0;
	  int ICS=0;
	  

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////   FIRST INFECTION ROUND with INOCULUM  --  cells are infected at the same time  ///////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	  CellsNonReplic<<"passage   "<<endl;
	  int CNR=0;

	  RoundOfInfection++; 
	  PassageSimTime+= TimeIncrement;
	  
	  Dynamics<<"Simul time     #Newly inf cells    #Infected with replicating virus     #Infected with non replic virus"<<endl;

	  
	  
	  //cout<<endl<<"Sim. Time is  "<<PassageSimTime<<"Infection with inoculum"<<endl;
	  
	  QSCloud.InfectCellPopulationClumped(NotInfectedCells,NewlyInfectedCells,ICS, MaxNumberInfectingVirions,MaxNumberVirAttachedToCell,CellR0);
	  //in the first round if MOI<1,  the QSCloud will be empty now and there will be two non-empty sets of cells - infected and non infected
	  
	 // cout<<"Infected cells number"<<InfectedCells.size()<<"Not Infected cells number"<<NotInfectedCells.size()<<endl;

	  	  
	  //////////////////////////////INOCULUM HAS INFECTED THE CELL POPULATION//////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  while(NotInfectedCells.size()+NewlyInfectedCells.size()+InfectedCells.size())//all cells were infected BUT virus was so defective that no real infections occurred; THIS IS NOT GOING TO HAPPEN SINCE UNSUCCESSFULLY INFECTED CELLS ARE RETURNED BACK TO THE SUSCEPTIBLE POOL, THIS BRANCH WILL NEVER HAPPEN AM KEEPING IT JUST IN CASE
	    {
	      if(RepeatInfection>=MaxRepeatInfection){cout<<"Failure to infect further!"; break;}
	      
		  if(NewlyInfectedCells.size())//here create new progeny in the newly infected cells and calculate time to burst for each of them
		    {
		      Factor=0; 	
		      Msum=0;
		      Mhelper=0;
		      Mi=0;
		      r1=0;
		      rho1=0;
		      Lhelper=0;
		      
		      
		      RepeatInfection=0;
		      ICS=1;//
		      
		      for(int bz=0; bz<NewlyInfectedCells.size(); bz++)
			{
			  cout<<endl<<endl<<"CELL  "<<bz<<" out of "<<NewlyInfectedCells.size()<<endl;
			  Lhelper=0;
			  Mhelper=0;//use Mhelper to define r1 and rho1
			  Msum=0;
			  Factor=0; 
			  //cout<<"Generating Cell cloud...."<<endl;
			  
			  for(list<Genotype>::iterator kw=NewlyInfectedCells[bz].GenotypesInfecting.begin(); kw!=NewlyInfectedCells[bz].GenotypesInfecting.end(); kw++)//here I am using a formula to calculate the progeny of each type according to fitness of replications and translation and the burst time of each cell
			    {
			      if((*kw).FitnessOfTranslation&&(*kw).FitnessOfReplication)
				{ 
				  Lhelper++;//calculate how many have non-zero rates
				  if(!Mhelper){Mi=0;r1=(*kw).FitnessOfTranslation;rho1=(*kw).FitnessOfReplication;Mhelper++;}//defining the first genotype to apply the formula (4) from the write up Types of virions
				  else{Mi=1/r1-1/(*kw).FitnessOfTranslation+NumberReplicationsPosToNeg*(1/rho1-1/(*kw).FitnessOfReplication);}
				  Msum+=Mi*(*kw).FitnessOfReplication; //this is to calculate the last sum in formula (4)
				  Factor+=(*kw).FitnessOfReplication; 
				}
			    }
			  
			  if(Lhelper)//if Lhelper is not 0, Factor is not 0, so CommonParam has a value; also if Lhelper is not 0, the cell would produce some virus, otherwise, no virus will be produced 
			    {
			      CommonParam = (BurstLoad-NumberReplicationsPosToNeg*Lhelper-Msum)/Factor;
			      CommonParam /= NumberReplicationsPosToNeg; 
			      
			      
			      for(list<Genotype>::iterator kw=NewlyInfectedCells[bz].GenotypesInfecting.begin(); kw!=NewlyInfectedCells[bz].GenotypesInfecting.end(); kw++)
				{
				  
				  ///////////////////////Calculate number of progeny this infecting type will create based on its fitness/////////////
				  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				  
				  
				  if(!(*kw).FitnessOfTranslation||!(*kw).FitnessOfReplication){NRPN=0; NRNP=0;}
				  else
				    {
				      NRPN=NumberReplicationsPosToNeg;
				      Mi=1/r1-1/(*kw).FitnessOfTranslation+NumberReplicationsPosToNeg*(1/rho1-1/(*kw).FitnessOfReplication);
				      NRNP=ceil((*kw).FitnessOfReplication*(Mi/NumberReplicationsPosToNeg+CommonParam));
				      
				      NewlyInfectedCells[bz].BurstTime= max(NewlyInfectedCells[bz].BurstTime, ceil(1/(*kw).FitnessOfTranslation+NRPN/(*kw).FitnessOfReplication*(1+ NRNP))); //here we calculate the max of the times in formula (1). In principle, all these times should be equal but due to round off error and taking ceil for NRNP they may not be equal
				      
				      //cout<<"NRNP= "<<NRNP<<"   ";
				      cout<<"calculating progeny of infecting genotype ...";
				      //cout<<"this genotype's fitness of replication "<<(*kw).FitnessOfReplication<<"  this genotype's fitness of translation "<<(*kw).FitnessOfTranslation<<endl;
				      
				      //It is assumed that each infecting genotype produces its polymerase and is replicated by it; 
				      //The replicative fitness of a genotype is determined by the fitness of the polymerase. 
				      //The same polymerase molecule is assumed to replicate different mutant positive and negative sequences with the same speed. 
				      //Within a single cell the number of genotypes produced would differ; 
				      //MaxRepl is introduced to model the situation that if two cells are infected with a sinle genotype each but one 
				      //is with less replicative fitness than the other, the two cells would produce different number of progeny
				      
				      
				      ///////////////////////////////////////////////////////////////////////////
				      ///END OF Calculate number of progeny this infecting type will create /////
				      ///////////////////////////////////////////////////////////////////////////
				      ///Next calculate the progeny with mutations //////////////////////////////
				      ///////////////////////////////////////////////////////////////////////////
				      
				      vector<char> InfectingSequence(NUMBER_OF_BASES);
				      vector<char> NegativeMutant(NUMBER_OF_BASES); //will be used for primary mutations (from  positive to negative or vice versa)
				      vector<char> PositiveMutant(NUMBER_OF_BASES); //will be used for secondary mutations (from negative to positive or vice versa)
				      
				      for(int i=0; i< NUMBER_OF_BASES; i++)
					{
					  InfectingSequence[i]=ReferenceSequence[i];
					}   
				      
				      //**********************************************
				      
				      
				      
				      //cout<<"now restoring infecting sequence from the mutations map..."<<endl;
				      
				      if((*kw).ntMutations.size())
					{ 
					  for(map<int,char>::const_iterator ps=(*kw).ntMutations.begin(); ps!=(*kw).ntMutations.end(); ps++)
					    {
					      int tt=(*ps).first;
					      InfectingSequence[tt]=(*ps).second; 
					      //cout<<"mutations in inf sec "<<tt<<InfectingSequence[tt]<<endl;	
					    }//now Infecting sequence has the restored sequence of the infecting genotype
					}
				      else
					{
					  //cout<<"infecting sequence has no mutations"<<endl;
					}
				      
				      
				      for(int i=0; i<NUMBER_OF_BASES; i++){NegativeMutant[i]=InfectingSequence[i];} //prepare template infecting sequence 
				      for(int i=0; i<NUMBER_OF_BASES; i++){PositiveMutant[i]=ReferenceSequence[i];} //these initializations need to be here in case there are no mutations in infecting sequence
				      
				      
				      //*************************BEGIN LOOPS OF REPLICATION AND MUTATION***************************
				      
				      kPM=0;
				      
				      ///////////////////replicate original sequence NumberReplicationsPosToNeg times and for each replicate decide whether and how many mutations it will have; 
				      ///////////////////then add it to the appropriate bin in the cloud; NOTE NumberReplicationsPosToNeg>0./////////////////////////////////////////////////////
				      
				      
				      while(kPM<NRPN)  
					{ 
					  Genotype G1;
					  G1.Initialize(NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation,RSFitnessOfEntry,RSFitnessOfReplication);
					  
					  //copy infecting  genotype to G1; G1 will be further mutated to simulate mutations in the transcription (from positive to negative sense, for +ss viruses and vice versa)
					  
					  if((*kw).ntMutations.size()){G1.CopyMutations((*kw).ntMutations);}
					  else{}
					  //next will add option to start with a mutant genotype called by the cell object
					  //next calculate the number of mutations that this replicate genome will receive      
					  
					  //G1.PrintToScreen_nt_Mutations();
					  
					  MTRand_closed prand; 
					  double pr=prand();  
					  
					  int NumM=MutationsNumber(MutationRate,NUMBER_OF_BASES,pr);//generate a random number  of mutations according to a binomial  distribution with probability of mutation MutationRate
					  
					  if(NumM)
					    {
					      G1.CreateMutations(InfectingSequence,NUMBER_OF_BASES,NumM,TransitionProb);//this creates a new ntMutations map from the InfectingSequence
					      //G1.PrintToScreen_nt_Mutations();
					      
					      TotalTNum+=G1.TNumber;//counting total number of transitions and transversions
					      TotalVNum+=G1.VNumber;
					      
					      
					      for(map<int,char>::const_iterator p=G1.ntMutations.begin(); p!=G1.ntMutations.end();p++)
						{
						  int yu=p->first;
						  NegativeMutant[yu]=p->second;//now copy the newly created mutations into the current infecting sequence; this will be further mutated
						}
					    }
					  
					  //G1.PrintToScreen_nt_Mutations();
					  
					  //cout<<"now replicate neg to pos"<<endl;
					  
					  kMP=0;
					  
					  while(kMP<NRNP)
					    {	  
					      Genotype G2; //create new object
					      G2=G1;// and copy the contents of the object G1
					      
					      MTRand_closed qrand;
					      double qr=qrand(); 
					      int PNumM=MutationsNumber(MutationRate,NUMBER_OF_BASES,qr);
					      
					      
					      if(PNumM)
						{
						  G2.CreateMutations(NegativeMutant,NUMBER_OF_BASES,PNumM,TransitionProb);  //... with respect to "new" ReferenceSequence (corresponding to previous mutation)
						  TotalTNum+=G2.TNumber;
						  TotalVNum+=G2.VNumber;
						}
					      
					      //G2.PrintToScreen_nt_Mutations();
					      G2.CompareToReference(ReferenceSequence);//here check if some of the mutations that appeared with respect to the last mutating sequence coincide with values of the reference sequence; if so, empty the mutations map of the object
					      
					      if(!G2.ntMutations.size())
						{
						  //cout<<"  no mutations";
						}
					      
					      //now add all new mutations to calculate fitness of the newly replicated genotype
					      
					      for(int i=0; i<NUMBER_OF_BASES; i++){PositiveMutant[i]=ReferenceSequence[i];}
					      for(map<int,char>::const_iterator p=G2.ntMutations.begin(); p!=G2.ntMutations.end();p++)
						{
						  int ty=p->first;
						  PositiveMutant[ty]=p->second;
						}
					      
					      //G2.PrintToScreen_nt_Mutations();
					      
					      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
					      //////////////////////DEFINE RATES OF ENTRY, REPLICATION and TRANSLATION for newly created genotype PositiveMutant//////////////////////
					      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					      
					      
					      CalculateFitnessOfFullSequence(NumberOfHF_5T,NumberOfHF_5R, vHF_5T, vHF_5R, HF_5T,HF_5R, NumberOfHF_3T,NumberOfHF_3R, vHF_3T,vHF_3R, HF_3T, HF_3R, NumberOfHF_STR, vHF_STR, HF_STR, NumberOfHF_NonSTR, vHF_NonSTR, HF_NonSTR, PositiveMutant,FTR, FE, FREPL, Start_3Prime,  FPrA, FPrC,  FPrG, FPrT, TPrA, TPrC, TPrG, TPrT, fpr, tpr, sstrn,nstrn, STRCodonFreqs, NSTRCodonFreqs, Start_Str_Prot_Region,Start_NStr_Prot_Region);
					      
					      
					      //cout<<endl<<"Mutant Fitness of translation  "<<G2.FitnessOfTranslation<<endl;
					      
					      G2.FitnessOfEntry=FE;
					      
					      //cout<<"Mutant FitnessOfEntry= "<<G2.FitnessOfEntry<<endl;
					      //////////////////////////////////////////////////////////////////	
					      
					      G2.FitnessOfReplication=FREPL;
					      
					      //cout<<"Mutant FitnessOfReplication= "<<G2.FitnessOfReplication<<endl;
					      
					      
					      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
					      ////////////////////////////////////ADDING the new mutant TO THE CELL CLOUD//////////////////////////////////////////////
					      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
					      
					      int hd=G2.ntMutations.size();  
					      
					      NewlyInfectedCells[bz].CellQSC.push_back(G2);
					      //cout<<"now add to cell cloud"<<endl;
					      
					      ////////////////////////////At  the end of the above loop all mutants have been added in CellCloud, 
					      //////////////////////////// and for ecach representative, the number of repeating sequences has been recorded in NumberOfTwins (member of Genotype)
					      
					      
					      //cout<<kMP<<endl;
					      kMP++; 
					      
					    }    
					  
					  kPM++; 
					}  //end of while loop - now the *kw genotype has been replicated
				      
				    }
				  
				  
				  //cout<<endl; 
				}
			      
			    }
			  
			  //NOW add NewlyInfectedCell to the list of infected cells only IF BurstTime is not 0 

			  if(NewlyInfectedCells[bz].BurstTime) //this is non zero if Lhelper is not 0
			    {
			      InfectedCells.push_back(NewlyInfectedCells[bz]);//add newly infected cell to the vector of infected
			      //cout<<"   BurstTime   "<<NewlyInfectedCells[bz].BurstTime;
			    }
			  else{InfectedNonproducingCells.push_back(NewlyInfectedCells[bz]);}
			  
			}// now all NewlyInfectedCells were renamed as either Infected or InfectedNonproducing
		      
		      CellsNonReplic <<"  Number of Infected with non replicating virus  "<<InfectedNonproducingCells.size();
		      
		      
		      Dynamics<<PassageSimTime<<" "<<NewlyInfectedCells.size()<<" "<<InfectedCells.size()<<" "<<InfectedNonproducingCells.size()<<endl;
		      NewlyInfectedCells.clear();
		    }
		  
		  //if newly infected cells were not 0, there are now infected cells ; if newly infected were 0, we next check if there are cells to infect
		  if(!NotInfectedCells.size())//if there are no cells left to infect, just release virus from all infected cells at once to cloud
		    {
		      if(InfectedCells.size())
			{
			  for(int g=0; g<InfectedCells.size(); g++)
				{
				  for(list<Genotype>::const_iterator kg=InfectedCells[g].CellQSC.begin(); kg!= InfectedCells[g].CellQSC.end(); kg++)
				    {
				      QSCloud.AddGenotypeToCloud(*kg);
				      //cout<<"genotype was added"<<endl;
				    }
				  TimeIncrement=InfectedCells[g].BurstTime;
				  BurstTimes<<TimeIncrement<<endl;
				  TotalNumberBursts++; 
				  cout<<"TotalNumberBursts  "<<TotalNumberBursts<<"    ";

				  if(BurstTdistr.count(TimeIncrement)>0)
				    {
				      spomag=BurstTdistr[TimeIncrement];
				      BurstTdistr[TimeIncrement]=spomag+1;
				    }
				  else
				    {
				      BurstTdistr[TimeIncrement]=1;
				    }
				}
			  InfectedCells.clear();
			}//when this is completed, the while loop should end

		    }
		  else//if there are noninfected cells, either there are infectious cells to infect them, or there is virus in cloud, or there is none.
		    {
		      if(!InfectedCells.size())//first, if there are no cells that release virus but there are cells to infect, try to infect them from cloud
			{
			  QSCloud.CalculateNumberOfMembersInCloud();
			  //cout<<"Number of members in cloud  "<<QSCloud.NumberOfMembers<<endl;
			  if(QSCloud.NumberOfMembers)//if there is something left in the cloud ...
			    {			       
			      if(RepeatInfection<MaxRepeatInfection)
				{ 
				  cout<<" Repeat INFECTION"<< RepeatInfection<<endl;
				  // if there were no infected cells in previous round, there is still virus in the cloud; the unsuccessfully attacked cells are still available to infect
				  RepeatInfection++;
				  QSCloud.InfectCellPopulationClumped(NotInfectedCells,NewlyInfectedCells,1, MaxNumberInfectingVirions,MaxNumberVirAttachedToCell,CellR0);//ICS = 1 
				  //cout<<"Infected cells number"<<InfectedCells.size()<<endl; 
				  //cout<<"NonInfected cells number"<<NotInfectedCells.size()<<endl;
				}
			      else{RepeatInfection=MaxRepeatInfection;}
			    }
			  else{RepeatInfection=MaxRepeatInfection;}//if there is nothing in the cloud and if there are no infected cells... this is the end of it- this will be recognized when next time number of cells is checked - in the beginning of loop
			}

		      else//if there are infected cells, and there are non-infected cells, calculate which infected cell will burst first
			{			  
			  TimeIncrement=InfectedCells[0].BurstTime;
			  indexvar=0;
			  RepeatInfection=0;//this is because new virus will be added to cloud to infect cells
			  
			  //now find  the cell with the shortest time to bursting
			  for(int iii=0; iii<InfectedCells.size(); iii++)
			    {
				  if(InfectedCells[iii].BurstTime<TimeIncrement)
				    {
				      indexvar=iii;
				      TimeIncrement=InfectedCells[iii].BurstTime;
				    }
				}
			      //now burst the cell with the shortest time to bursting
			      for(list<Genotype>::const_iterator kg=InfectedCells[indexvar].CellQSC.begin(); kg!= InfectedCells[indexvar].CellQSC.end(); kg++)
				{
				  QSCloud.AddGenotypeToCloud(*kg);
				  
				  //cout<<"genotype was added"<<endl;
				}
			      //cout<<"the virus that first burst was added to cloud"<<endl;

			      BurstTimes<<TimeIncrement<<endl;
			      TotalNumberBursts++;
			      cout<<"TotalNumberBursts  "<<TotalNumberBursts<<"    ";

			      if(BurstTdistr.count(TimeIncrement)>0)
				{
				  spomag=BurstTdistr[TimeIncrement];spomag++;
				  BurstTdistr[TimeIncrement]=spomag;
				}
			      else
				{
				  BurstTdistr[TimeIncrement]=1;
				}

			      PassageSimTime+= TimeIncrement;
			      
			      //cout<<endl<<"Burst Time is  "<<PassageSimTime;
			      
			      //now update InfectedCells list ...
			      TempInfectedCells.clear();
			      
			      for(int iii=0; iii<InfectedCells.size(); iii++) 
				{
				  if(iii!=indexvar)
				    {
				      TempInfectedCells.push_back(InfectedCells[iii]);
				    }
				}
			      
			      InfectedCells.clear();
			      
			      for(int iii=0; iii<TempInfectedCells.size(); iii++)
				{
				  InfectedCells.push_back(TempInfectedCells[iii]);
				}
			      /////////////////////////// 
			      
			      QSCloud.InfectCellPopulationClumped(NotInfectedCells,NewlyInfectedCells,ICS, MaxNumberInfectingVirions,MaxNumberVirAttachedToCell,CellR0);   //ISC=1 still
			} 
		    }
    
	    }
	
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////////////////All cells were infected and now logging final cloud produced by cell culture///////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  cout<<"Now logging cloud to file"<<endl;
	  CloudDimension=0;
	  
	  QSCloud.PrintCloudToFile(Sequences,SeqShort, NumberOfDominants, RoundOfInfection, CloudDimension, NUMBER_OF_BASES,NumberOfCodingRegions, RSFitnessOfTranslation, RSFitnessOfEntry, RSFitnessOfReplication);
	  Sequences<<endl<<"Time is  "<<PassageSimTime;

	  MutStat<<"passage   "<<cx<<endl;
	  QSCloud.PrintCloudHD(CloudDimension,MutStat); 
	  QSCloud.CalculateNumberOfMembersInCloud();
	  cout<<"Number of members in cloud  "<<QSCloud.NumberOfMembers<<endl;
	  
	  Titer<<"passage "<<cx<<"   ";
	  Titer<<QSCloud.NumberOfMembers<<endl;
	  
	  //if passage number is final, log the heat map
	  for(int g=0;g<NumberHMtoPrint; g++)
	    {
	      if(cx==HeatMapsToPrint[g])
		{
		  cout<<"Printing QSCloud's heat map"<<endl;
		  outHM<<"Passage "<<cx<<endl;
		  QSCloud.PrintCloudHeatMap(outHM,ReferenceSequence);
		}
	    }
	  
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  

	  BTDistr<<"Passage"<<cx<<endl;
	  for(map<double,int>::const_iterator sk=BurstTdistr.begin(); sk!=BurstTdistr.end(); sk++)
	    {
	      BTDistr<<"  "<<sk->first<<"  ";	
	    }
	  BTDistr<<endl;

	  for(map<double,int>::const_iterator sk=BurstTdistr.begin(); sk!=BurstTdistr.end(); sk++)
	    {
	      uyy=sk->second;
	      uyy=uyy/TotalNumberBursts;
	      BTDistr<<"  "<<uyy<<"  ";	
	    }
	  BTDistr<<"TotalNumberBursts   "<<TotalNumberBursts<<endl;

	  /////////////////NOW PICK A SAMPLE OF MOI*NUMBER_OF_CELLS FROM THE CLOUD TO INFECT THE NEXT PASSAGE//////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  vector<Genotype> NextCloud(NumberOfInfectingGenotypes);
	  
	  QSCloud.CalculateNumberOfMembersInCloud();
	  cout<<"Number of members in cloud  "<<QSCloud.NumberOfMembers<<endl;
	  if(!QSCloud.NumberOfMembers){cout<<"There is something wrong - no virus in cloud... quitting"<<endl;exit(1);}	 

	  QSCloud.PickSample(NumberOfInfectingGenotypes, NextCloud,FEdistribution,FRdistribution,FTdistribution,FERdistribution);//picks a random sample of size NumberOfInfectingGenotypes

	  //
	  FEDistr<<"Passage"<<cx<<endl;
	  if(QSCloud.NumberOfMembers)
	    {
	      for(map<double,int>::const_iterator sk=FEdistribution.begin(); sk!=FEdistribution.end(); sk++)
		{
		  FEDistr<<"  "<<sk->first<<"  ";	
		}//
	      
	      FEDistr<<endl;
	      
	      for(map<double,int>::const_iterator sk=FEdistribution.begin(); sk!=FEdistribution.end(); sk++)
		{
		  uyy=sk->second;
		  uyy=uyy/QSCloud.NumberOfMembers;
		  FEDistr<<"  "<< uyy<<"  ";	
		}
	      FEDistr<<endl;
	      //
	      
	      FRDistr<<"Passage"<<cx<<endl;
	      for(map<double,int>::const_iterator sk=FRdistribution.begin(); sk!=FRdistribution.end(); sk++)
		{
		  FRDistr<<"  "<<sk->first<<"  ";	
		}//
	      FRDistr<<endl;
	      
	      for(map<double,int>::const_iterator sk=FRdistribution.begin(); sk!=FRdistribution.end(); sk++)
		{
		  uyy=sk->second;
		  uyy=uyy/QSCloud.NumberOfMembers;
		  FRDistr<<"  "<<uyy <<"  ";	
		}
	      FRDistr<<endl;
	      //
	      
	      FTDistr<<"Passage"<<cx<<endl;
	      for(map<double,int>::const_iterator sk=FTdistribution.begin(); sk!=FTdistribution.end(); sk++)
		{
		  FTDistr<<"  "<<sk->first<<"  ";	
		}//
	      FTDistr<<endl;

	      for(map<double,int>::const_iterator sk=FTdistribution.begin(); sk!=FTdistribution.end(); sk++)
		{
		  uyy=sk->second;
		  uyy=uyy/QSCloud.NumberOfMembers;
		  FTDistr<<"  "<<uyy <<"  ";	
		} 
	      FTDistr<<endl;
	      
	      //
	      
	      FERDistr<<"Passage"<<cx<<endl;
	      for(map<double,int>::const_iterator sk=FERdistribution.begin(); sk!=FERdistribution.end(); sk++)
		{
		  FERDistr<<"  "<<sk->first<<"  ";	
		}//
	      FERDistr<<endl;
	      
	      for(map<double,int>::const_iterator sk=FERdistribution.begin(); sk!=FERdistribution.end(); sk++)
		{
		  uyy=sk->second;
		  uyy=uyy/QSCloud.NumberOfMembers;
		  FERDistr<<"  "<< uyy<<"  ";	
		}
	      FERDistr<<endl;
	    }
	  else
	    {
	      FEDistr<<"Cloud is empty"<<endl;
	      FRDistr<<"Cloud is empty"<<endl;
	      FTDistr<<"Cloud is empty"<<endl;
	      FERDistr<<"Cloud is empty"<<endl;
	    }
	  ///////////////////////////////////////////////////////

	  /*for(int i=0; i<NumberOfInfectingGenotypes;i++)
	    {
	    NextCloud[i].Print_nt_Mutations(Sequences);
	    NextCloud[i].PrintTwins(Sequences);
	    }*/
	  
	  QSCloud.InitializeCloud(NUMBER_OF_BASES, 0);//empty previous cloud
	  
	  nextPassage<<cx<<endl;
	  
	  for(int i=0; i<NumberOfInfectingGenotypes;i++)
	    {
	      QSCloud.AddGenotypeToCloud(NextCloud[i]);
	      NextCloud[i].Print_nt_MutationsOnly(nextPassage);
	    }
	  
	  //end of passage 
	}  
      
}  

  cout<<"Number0fpr= "<<Number0fpr<< "  Number0tpr=  "<<Number0tpr<< "  Number0STR=  "<<Number0STR<<"  Number0NSTR=  "<<Number0NSTR<<endl;
  
  Titer.close();   
  outHM.close();
  MutStat.close();     
  Sequences.close();
  TransStat.close();
  simInput.close();
  mutationsInput.close();
  nextPassage.close();
  thisPassage.close();
  BurstTimes.close();
  Dynamics.close();
  FEDistr.close();
  FRDistr.close();
  FTDistr.close();
  FERDistr.close();
  BTDistr.close();
  CellsNonReplic.close();
}
