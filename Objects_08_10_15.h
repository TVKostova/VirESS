%Copyright (c) <2017>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory
%Written by Tanya Kostova Vassilevska, kostova@llnl.gov, tan.v.kos@gmail.com
%LLNL Release number LLNL-CODE-727823
%All rights reserved.

%This file is part of <VirESS>, including the header file. For details, see the comments. 

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
#include <cstdio>
#include<list>
#include<set>

#include "MTRand.h"

using namespace std;

//*****************************

class Genotype
{
 public:
  /*Functions*/
  Genotype();
  ~Genotype();

  void Initialize(int, int, double, double, double);
  void CopyMutations(map<int,char>);
  void InsertMutation(int); //if genotype is destined to mutate (calculated in class CellCloud ) then calculate random mutations

  void CreateMutations(const vector<char> &, int, int, double);
  int InsertOneMutation(vector<char> &, int);

  void WhereMutations(vector<char>, int, int, int, int, int, int*, int*);
  void Print_nt_MutationsOnly(ofstream &);
  void PrintMutations(ofstream &, vector<vector<char> >, int&, int&);
  void PrintToScreen_nt_Mutations();

  void Print_nt_Mutations(ofstream &);    
  void Print_nt_Mutations_No_Fitnesses(ofstream& );                              

  void CompareToReference(vector<char> );
  void UpdateTwins();
  void PlusTwin();
  void PrintTwins(ofstream &);
  void PrintTwinsToScreen();
  void PrintFitnesses();
  

  /* Data*/
  
  map<int,char> ntMutations;//here keep nt mutations WITH RESPECT TO REFERENCE GENOME
  map<int,char> Five_Prime;
  map<int,char> Three_Prime;
  vector<map<int,char> >aaMutations;//here keep aa mutations WITH RESPECT TO REFERENCE GENOME; each map in the vector corresponds to a protein


  int NumberBases;
  int NumberCRegions;

  int NumberOfTwins;//genotype knows how many identical ones are there
  int Number_Mutations;
  int TNumber;
  int VNumber;

  double FitnessOfEntry;
  double FitnessOfReplication;
  double FitnessOfTranslation;


};



/////////////////////////////////////////CELL////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

class Cell
{
 public:

  Cell();
  ~Cell();

  void Initialize();
  
  //void GenerateCellQSCloud(double,double, int, int, int, int,int, int, int, int, int, vector<double>, vector<double>, vector<double>, vector<double>,vector<double>,vector<double>, vector<double>, vector<double>, vector<map<string,double> >, vector<map<string,double> >, vector<double>, vector<double>, vector<char>, vector<double>,  vector<char>, vector<double>, vector<char>, vector<double>,vector<char>,vector<double>,vector<char>, vector<double>, vector<char>,vector<double>, double, double, int, int&, int&, int,int, vector<char>, int*, int*, int&);

  double  LKq(int , int, double);
  int MutationsNumber(double, int, double);
  void PrintInfectingGenotypes();
  void PrintCellCloud();

  double* xCoord;
  double* yCoord;

  int NumberOfInfectingGenotypes;
  int NumberOfAttachedGenotypes;

  double BurstTime;

  vector<Genotype>  GenotypesAttaching;
  list<Genotype>  GenotypesInfecting;
  list<Genotype> CellQSC;

};
          
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////CLOUD//////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

class Cloud
{
 public: 
  Cloud();
  ~Cloud();

  void InitializeCloud(int,int);

  void AddGenotypeToCloud(Genotype);

  void InfectCellPopulationClumped(vector<Cell>&,vector<Cell>&, int, int, int, int);

  void PrintCloudHeatMap(ofstream&,vector<char>);

  void PrintCloudHD(int, ofstream&);

  void PrintCloudToScreen();

  void PrintCloudToFile(ofstream&,ofstream&, int, int, int&, int, int, double, double, double);

  void CalculateNumberOfMembersInCloud();

  void PickSample(int, vector<Genotype> &,map<double,int> &,map<double,int> &,map<double,int> &,map<double,int> &);

  int NumberOfBases;
  int NumberOfMembers;
  
  vector<list<Genotype> >QSC;

  Genotype Dominant;

};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/*
class Sequence
{ 
 public:
  Sequence();
  ~Sequence();
  void Initialize();
  void RestoreSequenceFromMutations();

  int NumberOfBases;
  char* ReferenceSequence;
  map<int,char>MutationsMap;void Genotype::Print_nt_Mutations_No_Fitnesses(ofstream &printfile)


}

*/
