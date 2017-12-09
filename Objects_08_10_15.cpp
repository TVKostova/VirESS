
%Copyright (c) <2017>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory
%Written by Tanya Kostova Vassilevska, kostova@llnl.gov, tan.v.kos@gmail.com
%LLNL Release number LLNL-CODE-727823
%All rights reserved.

%This file is part of <VirESS>, including definitions of object functions. For details, see the comments. 

%Licensed under the Apache License, Version 2.0 (the “Licensee”); you may not use this file except in compliance with the License.  You may obtain a copy of the License at:  http://www.apache.org/licenses/LICENSE-2.0

%Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the license.
/*
 *
 *
 *
 * Author Tanya Kostova Vassilevska
 */

#include <stdlib.h>
#include <math.h>
#include<vector>
#include <iostream>
#include <fstream>
#include <string>
#include<vector>
#include<map>
#include<set>
#include <cstdio>
#include<list>


#include "Objects_08_10_15.h"
#include "MTRand.h"

//#define NUMBER_OF_BASES 10273
//extern char* ReferenceSequence; 

using namespace std;


////////////////////////////////GENOTYPE///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Genotype::Genotype()
{
  Number_Mutations=0;
  //ntMutations.insert(make_pair(0,'n')); //if only (0,'n') then genotype has no mutation;
};


void Genotype::Initialize(int NB, int NCR, double RS_FTR, double RS_FE, double RS_FREPL)
{
  NumberBases=NB;
  NumberCRegions=NCR;
  
  for(int i=0; i<NCR; i++)
    {
      map<int,char> Map;
      Map.clear();
      aaMutations.push_back(Map);
      }
  
  NumberOfTwins=0;
  Number_Mutations=0;
  ntMutations.clear();
  
  Five_Prime.clear();
  Three_Prime.clear();  

  for(int i=0;i<NumberCRegions;i++){
   aaMutations[i].clear();}
  
  FitnessOfEntry=RS_FE;
  FitnessOfReplication=RS_FREPL;
  FitnessOfTranslation=RS_FTR;
};


Genotype::~Genotype(){};

void Genotype::CopyMutations(map<int,char> MM)
{
  for(map<int,char>::const_iterator p=MM.begin(); p!=MM.end();p++)
    {
      ntMutations[p->first]=p->second;
    }
};

void Genotype::PrintMutations(ofstream &printfile, vector<vector<char> > RefProt, int& Synonymous, int& Nonsynonymous)
{
  //cout<<"Print mutations..."<<endl;
  map<int,char>::const_iterator p;
  
   
  for(p=ntMutations.begin();p!= ntMutations.end(); p++)
    {
      //cout<<p->first<<p->second<<"   ";
      printfile<<p->first<<p->second<<"   ";
    }
  
  if(Five_Prime.size()) 
    {
      //cout<<endl<<"Five_Prime mutations:"<<endl;
      printfile<<"... "<<"Five Prime mutations: ";
      
      for(p=Five_Prime.begin();p!= Five_Prime.end(); p++)
	{
	  //cout<<p->first<<p->second<<"   ";
	  printfile<<p->first<<p->second<<"   ";
	}
   }

 if(Three_Prime.size())
   {
     //cout<<"Three_Prime mutations:"<<endl;
     printfile<<"... "<<"Three_Prime mutations: ";
     
     for(p=Three_Prime.begin();p!= Three_Prime.end(); p++)
       {
	 //cout<<p->first<<p->second<<"   ";
	 printfile<<p->first<<p->second<<"   ";
       }
   }
 
 for(int i=0; i<NumberCRegions; i++)
   {
     if(aaMutations[i].size())
       {
	 //cout<<endl<<"Mutations in Coding region "<<i+1<<endl;
	 printfile<<"... "<<"Mutations in Coding region "<<i+1<<" ";
	 
	 for(p=aaMutations[i].begin();p!= aaMutations[i].end(); p++)
	   {
	     printfile<<p->first<<RefProt[i][p->first-1]<<" -> ";
	     //cout<<p->first<<p->second<<"   ";
	     printfile<<p->first<<p->second<<"   ";
	     if(RefProt[i][p->first-1]==p->second)
	       {
		 Synonymous++;
	       }
	     else
	       {
		 Nonsynonymous++;
	       }
	   }
       }
   }
 
 //cout<<endl<<"Number of Synonimous Mutations: "<<Synonymous<<"  Number of Synonymous Mutations: "<<Nonsynonymous;
 //printfile <<endl<<"Number of Synonymous Mutations: "<<Synonymous<<"  Number of non-Synonymous Mutations: "<<Nonsynonymous;

};

void Genotype::Print_nt_MutationsOnly(ofstream &printfile)
{
  map<int,char>::const_iterator p;
  if(!ntMutations.size()){printfile<<-1<<'n';}
  else
    {
      for(p=ntMutations.begin();p!= ntMutations.end(); p++)
	{
	  printfile<<p->first<<" "<<p->second<<"   ";
	}
    }
  printfile<<-9<<" "<<endl;
};

void Genotype::Print_nt_Mutations(ofstream &printfile)
{
  map<int,char>::const_iterator p;
  if(!ntMutations.size()){printfile<<"NO MUT  ";}
  
  for(p=ntMutations.begin();p!= ntMutations.end(); p++)
    {
      //cout<<p->first<<p->second<<"   ";
      printfile<<p->first<<" "<<p->second<<"   ";
    }
  
  printfile<<"  fTr " <<FitnessOfTranslation<<"  fRepl "<<FitnessOfReplication<<"  fEntry "<<FitnessOfEntry<<"  ";

};

void Genotype::Print_nt_Mutations_No_Fitnesses(ofstream &printfile)
{
  map<int,char>::const_iterator p;
  if(!ntMutations.size()){printfile<<"NO MUT  ";}
  
  for(p=ntMutations.begin();p!= ntMutations.end(); p++)
    {
      //cout<<p->first<<p->second<<"   ";
      printfile<<p->first<<p->second<<"   ";
    }
  
  };


void Genotype::PrintToScreen_nt_Mutations()
{
  cout<<"Print mutations..."<<endl;
  map<int,char>::const_iterator p;
  
  if(!ntMutations.size()){cout<<"no mutations..."<<endl;
  };
  
  for(p=ntMutations.begin();p!= ntMutations.end(); p++)
    {
      cout<<p->first<<" "<<p->second<<"   ";
    }
  cout<<endl;
};


void Genotype::UpdateTwins()
{
  NumberOfTwins=0;
};

void Genotype::PlusTwin()
{
  NumberOfTwins++;
};


void Genotype::PrintTwins(ofstream &printfile)
{
  //cout<<" Repeats:  "<< NumberOfTwins<<" times. "<<endl;
  printfile<<" :  Copies  "<< NumberOfTwins+1<<":";
};



void Genotype::PrintTwinsToScreen()
{
  //cout<<" Repeats:  "<< NumberOfTwins<<" times. "<<endl;

};

void Genotype::PrintFitnesses()
{
  cout<<"Fitness of Replication  "<<FitnessOfReplication<<"Fitness Of Translation  "<<FitnessOfTranslation<<"Fitness Of Entry  "<<FitnessOfEntry<<endl;
};

void Genotype::CreateMutations(const vector<char> &RS,int NB, int r, double TransitionProb) //for a given r and a sequence RS generates r (r>0!) mutations with respect to RS and stores them with their positions in map ntMutations; TNumber and VNumber are the numbers of transitions, transversions
{
 
  TNumber=0;
  VNumber=0;

  char novo;

  int NumberOfBases=NB; 

  int count=0; 
  int prod=1;
  int prod1=1;
  
  vector<int> Positions;//will use this to check whether next position of mutation is new
  bool Transition;
  
  while(count<r)
    {      
      MTRand frand; 
      double ss=frand();
      int ipos=floor(ss*NumberOfBases); 
            
      //this is the candidate mutation position; if it is different from the previous, it becomes a real mutation position
      
      if(count==0){Positions.push_back(ipos);}//this is the first position
      else
	{
	  for(int i=0; i<count; i++){prod1=prod*(Positions[i]-ipos);}
	}
      if(prod1){
	prod=prod1;
	
	Positions.push_back(ipos);
	

	char ref=RS[ipos];//this value will be mutated
	//cout<<"ref "<<RS[ipos]<<endl;
	
	MTRand_closed xrand; 
	MTRand_closed xxrand;
	
	
	double x=xrand(); 
	double xx=xxrand();
	
	double TR=TransitionProb;
	
	if(x<=TransitionProb)
	  {Transition=true;}
	else
	  {Transition=false;}
	
	
	if(ref=='a')
	  {
	    if( Transition==true){novo='g';}
	    else
	      {
		if(xx<=0.5)
		  {novo='c';}
		else
		  {novo='t';}
	      }
	    
	  }
	
	if(ref=='g')
	  {
	    if( Transition==true){novo='a';}
	    else
	      {
		if(xx<=0.5)
		  {novo='c';}
		else
		  {novo='t';}
	      }
	    
	  }
	
	if(ref=='c')
	  {
	    if( Transition==true){novo='t';}
	    else
	      {
		if(xx<=0.5)
		  {novo='a';}
		else
		  {novo='g';}
	      }
	    
	  }
	
	if(ref=='t')
	  {
	    if( Transition==true){novo='c';}
	    else
	      {
		if(xx<=0.5)
		  {novo='a';}
		else
		  {novo='g';}
	      }
	    
	  }
	
	if(ref!='a'&&ref!='c'&&ref!='g'&&ref!='t'){cout<<"UNDEFINED POSITION"<<ipos; exit(1);}
       
	ntMutations[ipos]=novo; //add the new mutation and its position to the map of the object
	
	count++;
	
      }
      else{
	cout<<"Prod=0"<<endl;
      }
      if(Transition){TNumber++;}
      else{VNumber++;}
    }
  
};



void Genotype::CompareToReference(vector<char> RS)//calculates Hamming Distance from a given reference sequence RS and removes revertant mutations from list
{
  for(map<int,char>::const_iterator p=ntMutations.begin(); p!=ntMutations.end();p++)
    {
      if(p->second==RS[p->first]){ntMutations.erase(p->first);}
      else{}
    }
};




void Genotype::WhereMutations(vector<char> SN, int FiveS, int FiveE, int ThreeS, int ThreeE, int NCR, int* SCR, int* ECR)//calculates where in genome the  mutations that happened are ; SCR is an array holding the beginnings of NCR coding regions, ECR similarly for ends; 
{
  vector<char> NT(3);
  vector<char> A(1);


  //create Genotype::aa maps
  for(int i=0; i<NCR; i++)
    {
      map<int, char> AA_map; 
      aaMutations.push_back(AA_map);//create a vector of NRC maps to hold  aaMutants
    }int helper1=1;
			      int RS_helper1=1;
  
  if(ntMutations.size())
    {
      for(map<int,char>::const_iterator iter=ntMutations.begin(); iter!=ntMutations.end(); iter++) 
	{
	  int ntPos=iter->first;//position of mutation; 

	  //cout<<"pos of mutation: "<<ntPos<<endl;

	  if(ntPos>=FiveS&&ntPos<=FiveE){Five_Prime[ntPos]=iter->second; 
	    //cout<<"mut in 5 prime region: "<<ntPos<<" "<<Five_Prime[ntPos]<<" ";
	  }//mutation is in 5prime non-coding region

	  if(ntPos>=ThreeS&&ntPos<=ThreeE){Three_Prime[ntPos]=iter->second;
	    //cout<<"mut in 5 prime region: "<<ntPos<<" "<<Three_Prime[ntPos]<<" ";
	  }

	  else{
	    for(int x=0; x<=NCR; x++)
	      {
		int w;
		
		if(ntPos>=SCR[x]&&ntPos<=ECR[x])
		  {
		    //cout<<"mut in "<<x+1<< "coding region: "<<endl;
		    //calculate aa mutation from nt mutation
		    int y=ntPos-SCR[x]+1;
		    int z=fmod(y,3); 
		    if(z)
		      {
			w=(y-z)/3+1;

			NT[0]=SN[ntPos-z+1];NT[1]=SN[ntPos-z+2];NT[2]=SN[ntPos-z+3];
			//cout<<"nucleotide positions are "<<ntPos-z+1<<SN[ntPos-z+1]<<" "<<ntPos-z+2<<SN[ntPos-z+2]<<" "<<ntPos-z+3<<SN[ntPos-z+3]<<" "<<endl;
		      }
		    else
		      {
			w=y/3;
			//cout<<"codon number is "<<w<<" ";

			NT[0]=SN[ntPos-2];NT[1]=SN[ntPos-1];NT[2]=SN[ntPos];
			//cout<<"nucleotide positions are "<<ntPos-2<<SN[ntPos-2]<<" "<<ntPos-1<<SN[ntPos-1]<<" "<<ntPos<<SN[ntPos]<<" "<<endl;
		      
		      }//this is the codon number in the x protein which has a mutation



		    if(NT[0]=='t'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='c')
		      {
			A[0]='F'; //cout<<A[0];
		      }
		    if(NT[0]=='c'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='t'&&NT[2]=='g'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='t'&&NT[2]=='g')
		      {
			A[0]='L'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='a'&&NT[1]=='t'&&NT[2]=='a')
		      {
			A[0]='I'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='t'&&NT[2]=='g')
		      {
			A[0]='M'; //cout<<A[0];
		      }
		    if(NT[0]=='g'&&NT[1]=='t'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='t'&&NT[2]=='g')
		      {
			A[0]='V'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='c'&&NT[2]=='g')
		      {
			A[0]='S'; //cout<<A[0];
		      }
		    if(NT[0]=='c'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='c'&&NT[2]=='g')
		      {
			A[0]='P'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='c'&&NT[2]=='g')
		      {
			A[0]='T'; //cout<<A[0];
		      }
		    
		    if(NT[0]=='g'&&NT[1]=='c'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='c'&&NT[2]=='g')
		      {
			A[0]='A'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='a'&&NT[2]=='c')
		      {
			A[0]='Y'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='t'&&NT[1]=='a'&&NT[2]=='g')
		      {
			A[0]='*'; //cout<<A[0];
		      }
		    if(NT[0]=='c'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='a'&&NT[2]=='c')
		      {
			A[0]='H'; //cout<<A[0];
		      }
		    if(NT[0]=='c'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='a'&&NT[2]=='g')
		      {
			A[0]='Q'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='a'&&NT[2]=='c')
		      {
			A[0]='N'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='a'&&NT[2]=='g')
		      {
			A[0]='K'; //cout<<A[0];
		      }
		    if(NT[0]=='g'&&NT[1]=='a'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='a'&&NT[2]=='c')
		      {
			A[0]='D'; //cout<<A[0];
		      }
		    if(NT[0]=='g'&&NT[1]=='a'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='a'&&NT[2]=='g')
		      {
			A[0]='E'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='t'&&NT[1]=='g'&&NT[2]=='c')
		      {
			A[0]='C'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='a')
		      {
			A[0]='*'; //cout<<A[0];
		      }
		    if(NT[0]=='t'&&NT[1]=='g'&&NT[2]=='g')
		      {
			A[0]='W'; //cout<<A[0];
		      }
		    if(NT[0]=='c'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='c'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='c'&&NT[1]=='g'&&NT[2]=='g')
		      {
			A[0]='R'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='a'&&NT[1]=='g'&&NT[2]=='c')
		      {
			A[0]='S'; //cout<<A[0];
		      }
		    if(NT[0]=='a'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='a'&&NT[1]=='g'&&NT[2]=='g')
		      {
			A[0]='R'; //cout<<A[0];
		      }
		    if(NT[0]=='g'&&NT[1]=='g'&&NT[2]=='t'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='c'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='a'||NT[0]=='g'&&NT[1]=='g'&&NT[2]=='g')
		      {
			A[0]='G'; //cout<<A[0];
		      }
		    
		    
		    aaMutations[x][w]=A[0];
		  }
	      }
	  }
	}
    }
};


////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CLASS CELL ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


Cell::Cell(){};
Cell::~Cell(){};

void Cell::Initialize()
{
  NumberOfInfectingGenotypes=0;
  NumberOfAttachedGenotypes=0;
  BurstTime=0;

  GenotypesInfecting.clear();
  GenotypesAttaching.clear();

  CellQSC.clear(); 
};


double Cell::LKq(int L, int k, double q)
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



int Cell::MutationsNumber(double MutRate, int SeqLength, double rnd)
{
  int NM=0;
  int L=SeqLength;
   
  double prob=pow(1-MutRate, SeqLength);
  

  if(prob>=rnd)return NM;
  else
    {
      while(prob<rnd)
	{
	  NM++;
	  prob+=pow(1-MutRate,L-NM)*LKq(L,NM,MutRate); 
	} 
    };
  return NM; 
}


void Cell::PrintCellCloud()
{
  list<Genotype>::iterator qq;
  for(qq=CellQSC.begin(); qq!=CellQSC.end();qq++)
    {
	  (*qq).PrintToScreen_nt_Mutations();
	  (*qq).PrintTwinsToScreen();
    } 
};


void Cell::PrintInfectingGenotypes()
{
  list<Genotype>::iterator qq;

  for(qq=GenotypesInfecting.begin(); qq!=GenotypesInfecting.end();qq++)
    {
      (*qq).PrintToScreen_nt_Mutations();
      
    } 
};


/////////////////////////////////CLASS CLOUD//////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

Cloud::Cloud(){};
Cloud::~Cloud(){};


void Cloud::InitializeCloud(int NB, int NM)
{
  NumberOfBases=NB;
  NumberOfMembers=NM;

  QSC.resize(NumberOfBases+1);//will put types without mutations in 0 position of vector, types with mutations in nt =1, in position 1 etc.
  
  for(int i=0;i<NB+1;i++)
    {
       QSC[i].clear();
    };
};
 


void Cloud::AddGenotypeToCloud(Genotype G)//this version is created to reduce the number of comparisons with existing mutants; mutants are added to cloud ordered in sets by the first postion of the mutation; 
{
  char t;
  list<Genotype>::iterator kb;
  //cout<<"First print genotype mutations..."<<endl;
  //G.PrintToScreen_nt_Mutations();
  //G.PrintFitnesses();

  //cout<<" >>..  "<<endl;
  //cout<<"push button"<<endl;
  //cin>>t;
  
  if(G.ntMutations.size())
    {
      //cout<<"G.ntMutations.size()"<<G.ntMutations.size()<<endl;

      bool Flag=false;
      map<int,char>::const_iterator pw=G.ntMutations.begin();
      
      int first=pw->first;
      //cout<<"attempting to add to cloud... first mutant position is    "<<first<<endl;
      //cout<<"push button"<<endl;


      //cout<<"QSC["<<first+1<<"].size()"<<endl;
      //cout<<QSC[first+1].size()<<endl;
      //cout<<"push button"<<endl;
      //cin>>t;
      
      if(QSC[first+1].size()) //genotype with its first mutation in position *first* exists already- now check if it is the same as G
	{

	  //cout<<endl;
	  //cout<<"genotype with its first mutation in position *first* exists already ...."<<endl;

	  //cout<<"push button"<<endl;
	  //cin>>t;

	  for(list<Genotype>::iterator k=QSC[first+1].begin();k!= QSC[first+1].end();k++) 
	    { 		  
	      Flag=true; 
	      map<int,char>::iterator pp,qq;

	      //cout<<"after  map<int,char>::iterator pp,qq"<<endl;
	      //cout<<"push button"<<endl;
	      //cin>>t;

	      if(G.ntMutations.size()!=(*k).ntMutations.size())//this test should make calculations faster
		{
		  Flag=false;
		}
	      else
		{
		  for(pp=G.ntMutations.begin(), qq=(*k).ntMutations.begin(); pp!=G.ntMutations.end()||qq!=(*k).ntMutations.end(); pp++, qq++)//... check if its mutations are the same as G's 
		    {
		      if((*pp).first!=(*qq).first){Flag=false; break;}
		      if((*pp).second!=(*qq).second){Flag=false; break;}
		    }
		}
	      
	      if(Flag==true)//this will be true if all mutations of G and (*k) are identical
		{
		  //(*k).PrintToScreen_nt_Mutations();
		  (*k).PlusTwin();//if G and current member have the same mutations, increase number of repetitions of member and do not add GRef
		  //cout<<"adding twin to existing genotype with mutations.."<<endl;
		  //(*k).PrintTwinsToScreen();
		  break;
		}
	    }
	  
	  if(Flag==false)
	    { //if G and current member do not have the same mutations, add G to list; note here it is assumed that all genotypes do not have twins
	      QSC[first+1].push_back(G);
	      //cout<<"adding this NEW genotype with existing NEW HD  to cloud:"<<endl; 
	      //QSC[first+1].back().PrintToScreen_nt_Mutations();
	    }
	  
	  
	}
      else      
	{ //cout<<"adding this NEW genotype to cloud:"<<endl; 
	  QSC[first+1].push_back(G); 	 
	  
	  //QSC[first+1].back().PrintToScreen_nt_Mutations();
	}
    }
  else//G has no mutations w.r. to Reference 
    {
      if(!QSC[0].size())
	{
	  cout<<"adding genotype with no mutations  to cloud:"<<endl; 
	  QSC[0].push_back(G);
	  QSC[0].back().PrintToScreen_nt_Mutations();
	}
      else
	{
	  //cout<<"adding twin of genotype with no mutations..."<<endl;
	  kb=QSC[0].begin();
	  (*kb).NumberOfTwins++;
	  
	  (*kb).PrintTwinsToScreen();
	}
    }
  //cout<<"end adding genotype to cloud"<<endl;
};



void Cloud::PrintCloudHeatMap(ofstream& outHM, vector<char> RS)
{

  vector<map<char,double> > HeatMap(NumberOfBases);
  char nucl;
    
  NumberOfMembers=0;
  
  for(int i=0; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  for(list<Genotype>::const_iterator kk=QSC[i].begin();kk!=QSC[i].end(); kk++)
	    {
	      NumberOfMembers+=(*kk).NumberOfTwins+1;
	    }
	}
    }



  
      for(int g=0; g<NumberOfBases; g++)
	{
	  nucl= RS[g];
	  
	  HeatMap[g]['a']=0;
	  HeatMap[g]['c']=0;
	  HeatMap[g]['g']=0;
	  HeatMap[g]['t']=0;
	  HeatMap[g][nucl]=NumberOfMembers;//here create initially a heat map consisting only of RS members
	}
  


    
  
  //cout<<"calc HM"<<endl;
  for(int i=1; i<NumberOfBases+1; i++)
    {
      //cout<<"base "<<i<<endl;
      if(QSC[i].size())
	{
	  for(list<Genotype>::const_iterator k=QSC[i].begin();k!= QSC[i].end();k++)
	    {
	      map<int,char>::const_iterator qq;
	      
	      int w=1+(*k).NumberOfTwins;
	      
	      for(qq=(*k).ntMutations.begin();qq!=(*k).ntMutations.end();qq++)
		{
		  nucl=qq->second;
		  HeatMap[qq->first][nucl]+= w;

		  nucl=RS[qq->first];
		  HeatMap[qq->first][nucl]-=w;
		}
      
	    }
	  
	}
      
    }
  


  for(int g=0;g<NumberOfBases;g++)
    {
      outHM<<g<<"    ";

      double rr=HeatMap[g]['a']/NumberOfMembers;
      outHM<<rr<<" ";
      
      rr=HeatMap[g]['c']/NumberOfMembers;
      outHM<<rr<<" ";
      
      rr=HeatMap[g]['g']/NumberOfMembers;
      outHM<<rr<<" ";
      
      rr=HeatMap[g]['t']/NumberOfMembers;
      outHM<<rr<<" ";
      
      outHM<<endl;
      
    }
};


void Cloud::PrintCloudHD(int Dim,ofstream& outfile)//calculates number of different types in the cloud with a given Hamming distance from reference
{
  vector<int> HD(Dim);
  for(int i=0; i< Dim; i++){HD[i]=0;}

  int d;

  if(QSC[0].size())
    {
      for(list<Genotype>::const_iterator k=QSC[0].begin();k!= QSC[0].end();k++)
	{HD[0]+=(*k).NumberOfTwins+1;}
    }
  
  for(int i=1; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  for(list<Genotype>::const_iterator k=QSC[i].begin();k!= QSC[i].end();k++)
	    {
	      d=(*k).ntMutations.size();
	      HD[d]+=(*k).NumberOfTwins+1;
	    }
	}
    }

  //for(int ic=0; ic<Dim; ic++)
  // {
  //  outfile<<endl<<"Number of  DIFFERENT genotypes with "<<ic<<" mutations: "<<HD[ic]<<"    "<<endl;
      
  //}

  for(int ic=0; ic<Dim; ic++)
    {
      outfile<<HD[ic]<<" ";
      
    }
};



void Cloud::CalculateNumberOfMembersInCloud()
{
  NumberOfMembers=0;
  
  for(int i=0; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  for(list<Genotype>::const_iterator kk=QSC[i].begin();kk!=QSC[i].end(); kk++)
	    {
	      NumberOfMembers+=(*kk).NumberOfTwins+1;
	    }
	}
    }
  cout<<"QSCloud has "<<NumberOfMembers<<"members"<<endl;
};



void Cloud::PrintCloudToScreen()
{
  for(int i=0; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  for(list<Genotype>::iterator k=QSC[i].begin();k!= QSC[i].end();k++)
	    {
	      (*k).PrintToScreen_nt_Mutations();
	      (*k).PrintTwinsToScreen();
	    }
	}	  
    }
  
};
     

void Cloud::PrintCloudToFile(ofstream& File,ofstream& FileShort, int NumberOfDominants, int RoundOfInfection, int& MHD, int NB, int NCR, double RSFitnessOfTranslation, double RSFitnessOfEntry, double RSFitnessOfReplication)//suppressed print to file as files become very large
{
  double Freq;
  MHD=0;//this variable calculates maximum Hamming distance
  //calculate number of members in cloud
       
       for(int i=0; i<NumberOfBases+1; i++)
	 {
	   if(QSC[i].size())
	     {
	       int y=0;
	       //cout<<"QSC["<<i<<"].size()"<<QSC[i].size()<<endl;
	       for(list<Genotype>::const_iterator kk=QSC[i].begin();kk!=QSC[i].end(); kk++)
		 {
		   y++;
		   NumberOfMembers+=(*kk).NumberOfTwins+1;
		 }
	       //if(y){cout<<"number of genotypes"<<y<<endl;}
	     }
	 }
  
       Dominant.Initialize(NB, NCR, RSFitnessOfTranslation, RSFitnessOfEntry, RSFitnessOfReplication);//Dominant is a genotype argumant of the class Cloud

  int MaxTwins=-1;
  int MaxTwinsOld=-1;

  int Pomni;
  int yc=0;
  
  //first identify the genotype in the cloud that has the highest frequency

  for(int i=0; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  //File<<"QSC["<<i<<"].size()  "<<QSC[i].size();
	      
	  for(list<Genotype>::iterator k=QSC[i].begin();k!= QSC[i].end();k++)//calculates cloud dimension and dominant genotype
	    {
	      int nts=(*k).ntMutations.size();
	      MHD=max(MHD, nts);

	      yc++;//this calculates the order number of the Genotype in the cloud
	      MaxTwins=max(MaxTwins,(*k).NumberOfTwins);

	      if(MaxTwins>MaxTwinsOld)
		{
		  Dominant=(*k);
		  MaxTwinsOld=MaxTwins;
		  Pomni=yc;//each time a Genotype is a candidate for a dominant, Pomni records its order number
		}
	     

	      //(*k).Print_nt_Mutations(File);
	      //(*k).PrintTwins(File);
	    }
	}	  
    }

  MHD++;//this is done because CloudDimension=max Hamming distance+1 to include the class of genotypes with 0 HD

    //File<<endl<<"DOMINANT GENOTYPE: "<<endl;
  Dominant.Print_nt_Mutations(File);
  Dominant.Print_nt_Mutations_No_Fitnesses(FileShort);

  //Dominant.PrintTwins(File);

  Freq=Dominant.NumberOfTwins+1;
  Freq/=NumberOfMembers;
  File<<Dominant.NumberOfTwins+1<<"/"<<NumberOfMembers<<" Freq "<<Freq<<endl;
  FileShort<<" Freq "<<Freq<<endl;
 
  //File<<"   Fitness of Translation  "<<Dominant.FitnessOfTranslation<<"   Fitness of Replication  "<<Dominant.FitnessOfReplication<<"   Fitness of Entry  "<<Dominant.FitnessOfEntry<<endl;
  

  //Now Pomni has remembered the order number of the dominant  mutant


  vector<int> SubMaxTwins(NumberOfDominants);
  vector<int> SubMaxTwinsOld(NumberOfDominants);
  vector<int> SubPomni(NumberOfDominants);
  vector<int> Sub(NumberOfDominants);
  int Prod;
  int MT;
  
  for(int d=0; d<NumberOfDominants; d++)
    {
      if(d==0){MT=MaxTwins;}
      else{MT=SubMaxTwins[d-1];}

      Sub[d]=0;//Sub now plays the same role as yc 
      SubMaxTwins[d]=0;
      SubMaxTwinsOld[d]=0;
      
      for(int i=0; i<NumberOfBases+1; i++)
	{
	  if(QSC[i].size())
	    {
	      // File<<"QSC["<<i<<"].size()  "<<QSC[i].size();

	      for(list<Genotype>::iterator k=QSC[i].begin();k!= QSC[i].end();k++)
		{
		  Sub[d]++;
		 
		  if((*k).NumberOfTwins<=MT)
		    {
		      Prod=Sub[d]- Pomni;
		      if(d)
			{
			  for(int g=0; g<d;g++ )
			    {
			      Prod*=(Sub[d]-SubPomni[g]);
			    }
			}
		      if(Prod)//here we check if yc[d] is different from all previous yc
			{
			  SubMaxTwins[d]=max(SubMaxTwins[d],(*k).NumberOfTwins);
			  
			  if(SubMaxTwins[d]>SubMaxTwinsOld[d])
			    {
			      Dominant=(*k);
			      SubMaxTwinsOld[d]=Dominant.NumberOfTwins;
			      SubPomni[d]=Sub[d];
			    }
			}
		    }
		  
		  // (*k).Print_nt_Mutations(File);
		  //(*k).PrintTwins(File);
		}
	    }	  
	}
      
         
      
      //if(d+1==1){File<<"st  ";}
      //if(d+1==2){File<<"nd ";};
      //if(d+1==3){File<<"rd ";};
      //if(d+1>3){File<<"th ";}
      //File<<" SUB-DOMINANT: "<<endl;
      Dominant.Print_nt_Mutations(File);
      Dominant.Print_nt_Mutations_No_Fitnesses(FileShort);

      //Dominant.PrintTwins(File);


      Freq=Dominant.NumberOfTwins+1;
        
      Freq/=NumberOfMembers;
      File<<Dominant.NumberOfTwins+1<<"/"<<NumberOfMembers<<" Freq "<<Freq<<endl;
      FileShort<<" Freq "<<Freq<<endl;       
    }
};



void Cloud::InfectCellPopulationClumped(vector<Cell>& NCP,vector<Cell>& ICP, int ICSize, int MIV, int MaxNumberVirAttachedToCell, int CellR0)//MIV=MaxNumberInfectingVirions
{
  NumberOfMembers=0;char t;
  
  //calculate number of members in cloud
  for(int i=0; i<NumberOfBases+1; i++)
    {
      if(QSC[i].size())
	{
	  int y=0;
	  //cout<<"QSC["<<i<<"].size()"<<QSC[i].size()<<endl;
	  for(list<Genotype>::const_iterator kk=QSC[i].begin();kk!=QSC[i].end(); kk++)
	    {
	      y++;
	      NumberOfMembers+=(*kk).NumberOfTwins+1;
	    }
	  //if(y){cout<<"number of genotypes"<<y<<endl;}
	}
    }

  cout<<"QSCloud has "<<NumberOfMembers<<"members"<<endl;
 
  
  int CellN;
  int V;
  int VICP;
  int VNCP;

  cout<<"size of infected cell population"<<ICP.size()<<endl;
 
  //calculate how many viruses will potentially attach to cells 
  if(!ICSize)
    {V=NCP.size()*MaxNumberVirAttachedToCell;}//this should only be true in the beginning of infection; in this case all genotypes from inoculum attach to a cell 
  else
    {
      VICP=ICSize*CellR0*MaxNumberVirAttachedToCell;//this limits the number of attaching genotypes so that infection proceeds more slowly
      VNCP=NCP.size()*MaxNumberVirAttachedToCell;
      V=min(VICP,VNCP);
    }
     
  if(NumberOfMembers<=V)//if viruses are less than non-infected cells*attachmment size... for each virus in cloud find a cell it attaches to; 
    {
      //cout<<"NumberOfMembers<=V...  NOW traversing all genotypes"<<endl;

      for(int i=0; i<NumberOfBases+1; i++)
	{
	  if(QSC[i].size())
	    {
	      //cout<<"now  at base  "<<i<<endl;
	      
 	      for(list<Genotype>::const_iterator kk=QSC[i].begin();kk!=QSC[i].end(); kk++)//go thru list of genotypes in cloud and pick a random cell from the noninfected to attach to
		{
		  //cout<<"number of twins "<<(*kk).NumberOfTwins<<endl;
		  for(int j=0; j<1+(*kk).NumberOfTwins; j++)
		    {
		      
		      MTRand_closed wrand; 
		      double wpr=wrand(); 
		      //cout<<"random number is"<<wpr<<endl;
		      
		      CellN=floor(wpr*NCP.size());//number of cell among first NumberOfCellsToInfect in the noninfected pool to be infected with this genotype

		      cout<<"cell number is "<<CellN<<endl;
		      		      
		      NCP[CellN].NumberOfAttachedGenotypes++; 
		      //cout<<"NumberOfAttachedGenotypes++"<<endl;
		      		      
		      //cout<<"push button"<<endl;
		      //cin>>t;

		      NCP[CellN].GenotypesAttaching.push_back((*kk));
		      //cout<<"GenotypesAttaching.push_back"<<endl;

		      //cout<<"push button"<<endl;
		      //cin>>t;
		    }
		  
		}
	    }
	}
      
      //cout<<"push button to clear QSC"<<endl;
      //cin>>t;

      for(int i=0; i<NumberOfBases+1; i++)
	{
	  QSC[i].clear(); //whole cloud has attached to cells and is empty now 
	}
   }

  else 
    {//create Auxil
      vector<int> Auxil(NumberOfBases+1);  //Auxil is an artificial  structure to enable query of cloud members with fewer traversing steps
      if(!QSC[0].size()){Auxil[0]=0;}//there are no non-mutant types in QSCloud
      else
	{
	  for(list<Genotype>::const_iterator k=QSC[0].begin(); k!=QSC[0].end(); k++)
	    {
	      Auxil[0]+=(*k).NumberOfTwins+1; 
	    }
	}
      
      //cout<<" Auxil[0] "<< Auxil[0]<<endl;
      
      for(int i=1; i<NumberOfBases+1; i++)
	{
	  if(!QSC[i].size())
	    {
	      Auxil[i]=Auxil[i-1];
	    }
	  else
	    {
	      Auxil[i]=Auxil[i-1];
	      for(list<Genotype>::const_iterator k=QSC[i].begin(); k!=QSC[i].end(); k++)
		{
		  Auxil[i]+=(*k).NumberOfTwins+1;
		}
	    }
	  if( Auxil[i]!=Auxil[i-1]){
	    //cout<<" Auxil["<< i<<"] "<< Auxil[i]<<" ";
	  }
	}  
      //end create auxil
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      set<int>SortedTemp;
      set<int>Sorted;//containing a sorted list of non-repeating randomly generated numbers that will be used to pick a random sample of infecting genotypes
      Sorted.clear();
      SortedTemp.clear();
      int GN;
      int VV;

      //pick V randomly chosen DIFFERENT numbers from 0 to NumberOfMembers i.e. pick V randomly chosen members of cloud
      //use a trick because if V is close to NumberOfMembers, finding these random numbers may take much time
      if(NumberOfMembers>=2*V){VV=V;}
      else
	{
	  VV=NumberOfMembers-V;
	}
      while(SortedTemp.size()<VV) 
	{
	  MTRand_closed lrand; 
	  double lpr=lrand(); 
	  GN=floor(lpr*(NumberOfMembers));//floor - to enable using members without mutation 
	  
	  // cout<<"generated infecting genotype number"<<GN<<endl;
	  {
	    SortedTemp.insert(GN);
	  }
	}

      if(NumberOfMembers>=2*V)
	{
	  for(set<int>::const_iterator lm=SortedTemp.begin(); lm!=SortedTemp.end(); lm++) 
	    {
	      if(*lm||!*lm && QSC[0].size())//this condition is needed because if GN=0, and if there are no members without mutations, one less cell will be infected
		{
		  Sorted.insert(*lm);cout<<"  *lm  "<<*lm;
		}
	    }
	}
      else//take the complement
	{
	  int qs=0;
	  for(set<int>::const_iterator lm=SortedTemp.begin(); lm!=SortedTemp.end(); lm++) 
	    {cout<<"  *lm  "<<*lm;
	      while(qs<*lm)
		{
		  if(qs||!qs&&QSC[0].size())//this condition is needed because if GN=0, and if there are no members without mutations, one less cell will be infected
		    {
		      Sorted.insert(qs); cout<<"  qs  "<<qs;
		    }
		  qs++;
		} 
	      qs++;
	    }
	  
	  for(int qq=qs; qq<=NumberOfMembers; qq++){Sorted.insert(qq);cout<<"  qq  "<<qq;}
	}
      ////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////now V sequences have been picked//////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////

      int f;//this is cell number to be randomly picked
      int y=0;
      
      //cout<<"Now attack given number of  not infected cells..."<<endl;
     
 
      for(set<int>::const_iterator lm=Sorted.begin(); lm!=Sorted.end(); lm++)
	{
	  GN=(*lm);
	  //cout<<"Now finding what is genotype"<<GN<<endl;
	  
	  while(GN>Auxil[y])
	    {
	      y++; 
	    }

	  //cout<<"y for Auxil" <<y<<endl;
	  //cout<<"push button"<<endl;
	 // cin>>t;

	  //here GN < Auxil[y];	  
	  
	  int rr=0;
	  
	
	  for(list<Genotype>::iterator k=QSC[y].begin(); k!=QSC[y].end(); k++)//this is executed when GN<=Auxil[y]
	    {
	      rr+=(*k).NumberOfTwins+1;
	      
	      //cout<<"rr= "<<rr<<endl;
	      //cout<<"push button"<<endl;
	     // cin>>t;

	      if(GN<=Auxil[y-1]+rr)//the genotype to infect has been identified!!!
		{
		  //cout<<"genotype to infect has been identified "<<endl;
		  //cout<<"push button"<<endl;
		 // cin>>t;

		  MTRand_closed crand; 
		  double cpr=crand(); 
		  
		  //if(!ICSize||ICSize*CellR0>NCP.size())//In these cases: V=NCP.size()*MaxNumberVirAttachedToCell;this is when only a small number of noninfected cells have remained ; 
		  if(ICSize*CellR0>NCP.size())//removed first option - see above, because if !ICSize then NumbeOfMembers,+V and we are in the previous case
		    {
		      
		      //cout<<"ICSize*CellR0>NCP.size()); 

		      f=floor(cpr*NCP.size());
		      //cout<<"f="<<f<<endl;

		    }
		  else
		    {
		      f=floor(cpr*ICSize*CellR0);//In this case: ICSize*CellR0>NCP.size()) and V=ICP.size()*CellR0*MaxNumberVirAttachedToCell;
		      //cout<<"f="<<f<<endl;
		      
		    }//random cell number to infect 
		  
		  		  
		  NCP[f].GenotypesAttaching.push_back((*k));//attack cell f
		  
		  //NCP[f]-> NumberOfAttachedGenotypes=NCP[f]->GenotypesAttaching.size();
		  		  
		  break;
		}
	    }
	  
	}
      //next infect cells with viable genomes and update two populations - infected and not
      
      
      cout<<"NOW remove infecting genotypes from cloud"<<endl;
      

      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////	  
      ///////////////////////////////////////////////Remove infecting genotypes from cloud/////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
      y=0;
      
      for(set<int>::const_iterator lm=Sorted.begin(); lm!=Sorted.end(); lm++) 
	{	  
	  GN=(*lm);
	  
	  while(GN>Auxil[y]){y++;}
	  
	  
	  int rr=0;
	  
	  for(list<Genotype>::iterator k=QSC[y].begin(); k!=QSC[y].end(); k++)
	    {
	      rr+=(*k).NumberOfTwins+1;
	      
	      if(GN<=Auxil[y-1]+rr)
		{
		  if((*k).NumberOfTwins)
		    {
		      (*k).NumberOfTwins--; 
		      
		    }
		  else
		    {
		      
		      QSC[y].erase(k);
		      k--;
		      
		    }
		  
		  break;
		}
	    }
	}
    } 
  

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

 
 
  cout<<"Counting potentially infected cells "<<endl;   
  double rrr;
  int irr;
  
  cout<<"number of previously non-infected cells "<<NCP.size()<<endl;

  //cout<<"push button"<<endl;
 // cin>>t;

  for(int pp=0; pp<NCP.size(); pp++)
    {
      //cout<<"cell number "<< pp <<"is ";

      if(NCP[pp].GenotypesAttaching.size())
	{
	  multimap<double, Genotype> GASorted; //here store attaching virions with regard to fitness of entry
	 
	  for(int i=0; i<NCP[pp].GenotypesAttaching.size(); i++)
	    {
	      rrr=NCP[pp].GenotypesAttaching[i].FitnessOfEntry;
	      GASorted.insert(make_pair(-rrr,NCP[pp].GenotypesAttaching[i]));//GASorted is in ASCENDING order - this is why we take -rrr
	      //cout<<"fitness of entry of "<<i<<"attaching genotype "<<rrr<<"  ";
	    }

	  //cout<<"GASorted.size "<<GASorted.size()<<endl;


	  /*MTRand_closed drand; 
	  double dpr=drand(); 

	  irr=floor(dpr*(NCP[pp]->GenotypesAttaching.size()));
	  
	  irr=min(irr,MIV);//random number of infecting genotypes from the attached ones no more than MIV
	  */



	  int NGA=NCP[pp].GenotypesAttaching.size();

	  //cout<<"Attaching size of cell "<<pp<<" "<<NGA<<endl;

	  irr = min(MIV, NGA);
	  //cout<<"irr "<<irr<<endl;

	  
	  multimap<double, Genotype>::const_iterator p=GASorted.begin();
	  for(int h=0; h<irr; h++)//now take the irr genotypes with the highest entry rate
	    {
	      //cout<<"p->first  "<<(*p).first<<endl;

	      if(p!=GASorted.end()&&(*p).first)//if at least one genotype has non-zero entry rate
		{
		  NCP[pp].GenotypesInfecting.push_back((*p).second);
		  //cout<<" infected... "<<endl;
		}
	     p++;
	      
	    }
	  //cout<<"cell will be infected by "<<NCP[pp].GenotypesInfecting.size()<<"  genotypes";

	  /////////////// //now add this cell to the population of infected ones
	  if(NCP[pp].GenotypesInfecting.size())
	    {
	      ICP.push_back(NCP[pp]);
	    }
	}
    }	
  
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////Now update the non-infected population//////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  vector<Cell> Temp;
  Temp.clear();
   
  //cout<<"push button"<<endl;
 // cin>>t;


  for(int pp=0; pp<NCP.size(); pp++)
    {
      if(!NCP[pp].GenotypesInfecting.size()) //if this cell will not be infected, put it on the Temp list
	{
	  Temp.push_back(NCP[pp]); 
	}
      else{}
      
    }
  
  NCP.clear();
  //cout<<"<<<<<<<<<<<<<<<<In infectClumped... number of noninfected cells (1)   "<<NCP.size();
  
  //cout<<"push button"<<endl;
 // cin>>t;

  for(int pp=0; pp<Temp.size(); pp++)
    {
      NCP.push_back(Temp[pp]); //some of these Noninfected cells will have viruses attached that have zero fitness for entry though. Is this a problem???
    }

}


/////////////////////////////////////////////////////////////////////////////////////////
///////PickSample also aids the calculation of the distribution of the fitness rates/////

void Cloud::PickSample(int NumberOfInfectingGenotypes,vector<Genotype> &NextCloud,map<double,int> &FEdistribution,map<double,int> &FRdistribution,map<double,int> &FTdistribution,map<double,int> &FERdistribution)
{
 
  double w;  

  FEdistribution.clear();FRdistribution.clear();FTdistribution.clear();FERdistribution.clear();

  NextCloud.clear();

  //cout<<"#################################NextCloud.size"<<NextCloud.size()<<endl;

   //create Auxil
  vector<int> Auxil(NumberOfBases+1);  //Auxil is an artificial  structure to enable query of cloud members with fewer traversing steps
  if(!QSC[0].size())
    {Auxil[0]=0;}//there are no non-mutant types in QSCloud
  else
    {
      for(list<Genotype>::const_iterator k=QSC[0].begin(); k!=QSC[0].end(); k++)
	{
	  Auxil[0]+=(*k).NumberOfTwins+1;

	  if(FEdistribution.count((*k).FitnessOfEntry)>0)
	    {
	      w=FEdistribution[(*k).FitnessOfEntry];
	      FEdistribution[(*k).FitnessOfEntry]=w+(*k).NumberOfTwins+1;
	    }
	  else
	    {
	      FEdistribution[(*k).FitnessOfEntry]=(*k).NumberOfTwins+1;
	    }

	  if(FRdistribution.count((*k).FitnessOfReplication)>0)
	    {
	      w=FRdistribution[(*k).FitnessOfReplication];
	      FRdistribution[(*k).FitnessOfReplication]=w+(*k).NumberOfTwins+1;
	    }
	  else
	    {
	      FRdistribution[(*k).FitnessOfReplication]=(*k).NumberOfTwins+1;
	    }

	  if(FTdistribution.count((*k).FitnessOfTranslation)>0)
	    {
	      w=FTdistribution[(*k).FitnessOfTranslation];
	      FTdistribution[(*k).FitnessOfTranslation]=w+(*k).NumberOfTwins+1;
	    }
	  else
	    {
	      FTdistribution[(*k).FitnessOfTranslation]=(*k).NumberOfTwins+1;
	    }
	  
	  if(FERdistribution.count((*k).FitnessOfEntry)>0)
	    {
	      w=FERdistribution[(*k).FitnessOfReplication*(*k).FitnessOfEntry];
	      FERdistribution[(*k).FitnessOfReplication*(*k).FitnessOfEntry]=w+(*k).NumberOfTwins+1;
	    }
	  else
	    {
	      FERdistribution[(*k).FitnessOfReplication*(*k).FitnessOfEntry]=(*k).NumberOfTwins+1;
	    }
	}
    }

  for(int i=1; i<NumberOfBases+1; i++)
    {
      if(!QSC[i].size())
	{
	  Auxil[i]=Auxil[i-1];
	}
      else
	{
	  Auxil[i]=Auxil[i-1];
	  for(list<Genotype>::const_iterator k=QSC[i].begin(); k!=QSC[i].end(); k++)
	    {
	      Auxil[i]+=(*k).NumberOfTwins+1;

	      w=FEdistribution[(*k).FitnessOfEntry];
	      FEdistribution[(*k).FitnessOfEntry]=w+(*k).NumberOfTwins+1;
	      w=FRdistribution[(*k).FitnessOfReplication];
	      FRdistribution[(*k).FitnessOfReplication]=w+(*k).NumberOfTwins+1;
	      w=FTdistribution[(*k).FitnessOfTranslation];
	      FTdistribution[(*k).FitnessOfTranslation]=w+(*k).NumberOfTwins+1;
	      w=FTdistribution[(*k).FitnessOfTranslation];
	      FTdistribution[(*k).FitnessOfTranslation]=w+(*k).NumberOfTwins+1;
	      w=FERdistribution[(*k).FitnessOfReplication*(*k).FitnessOfEntry];
	      FERdistribution[(*k).FitnessOfReplication*(*k).FitnessOfEntry]=w+(*k).NumberOfTwins+1;
	    }
	}
    }  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  //Auxil created ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Create sorted set Sorted and populate it with different random numbers 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  set<int>Sorted;//containing a sorted list of non-repeating randomly generated numbers that will be used to pick a random sample of infecting genotypes
  int GN;
  
  
  while(Sorted.size()<NumberOfInfectingGenotypes) //pick NumberOfInfectingGenotypes randomly chosen DIFFERENT numbers from 0 to NumberOfMembers i.e. pick V randomly chosen members of cloud
    {
      MTRand_closed lrand; 
      double lpr=lrand(); 
      GN=floor(lpr*(NumberOfMembers));//NOTE: will never pick last member
      
      // cout<<"generated infecting genotype number"<<GN<<endl;
      if(GN||(GN==0&&QSC[0].size()))//this condition is needed because if GN=0, and if there are no members without mutations, one less cell will be infected
	{
	  Sorted.insert(GN);
	}
      
      //cout<<"Sorted.size() "<<Sorted.size()<<endl;
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  //Sorted  created and populated ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int y=0;
  
  for(set<int>::const_iterator lm=Sorted.begin(); lm!=Sorted.end(); lm++)
    {
      GN=(*lm);
      
      //cout<<"Now finding what is genotype with this sequential number (in cloud)"<<GN<<endl;
      
      while(GN>Auxil[y])
	{
	  y++; 
	}
      
       
      int rr=0;
      
      for(list<Genotype>::iterator k=QSC[y].begin(); k!=QSC[y].end(); k++)//this is executed when GN<=Auxil[y]
	{
	  rr+=(*k).NumberOfTwins+1;
	  
	  if(GN<=Auxil[y-1]+rr)//the genotype to add to next cloud has been identified
	    {
	      //cout<<"genotype to infect has been identified "<<endl;
	      
	      NextCloud.push_back((*k)); //add the genotype to list; however, it will inherit the number of "twins"
	      
	      break;
	    }
	}
      
    }

  for(int h=0; h<NumberOfInfectingGenotypes; h++)
    {
      NextCloud[h].NumberOfTwins=0;
    }
  
};








