
#include "PTFTTree.hh"
#include <cmath>
#include <vector>
#include "TGraph.h"
#include <iostream>
using namespace std;

void read_txt(){
	
	TTree *t =new TTree("t", "my tree");
	Double_t x,y,z;
	t->Branch("x", &x,"x/D");
	t->Branch("y", &y,"y/D");	
	t->Branch("z", &z,"z/D");  
    //TTree *t = new TTree("t", "my tree");
  t->ReadFile("grid.txt","x/D:y:z" );
  
  //Double_t v_min = t->GetMinimum("x");
  Double_t v_max = t->GetMaximum("x");
  
  cout<<v_max<<endl;
  
//TNtuple * position;
//TNtuple position("position","position","x:y:z");

//position->SetBranchAddress("x",&x);
//position->SetBranchAddress("y",&y);
//position->SetBranchAddress("z",&z);

//position->ReadFile("grid.txt");

  TH2D * pos = new TH2D( "pos", "Position; x (mm); y (mm)",
  		  			    41, -200.0, 200.0, 41, -200.0, 200.0 );	
 int nevents_pos=t->GetEntries();
 cout<<nevents_pos<<endl;
 t->Print()	;
	
  for (unsigned long long iev =0 ; iev < nevents_pos; ++iev ){
	  //	 t->GetEvent( iev );
//	 pos->Fill( t->x[0], t->y[0]);
	 t->GetEntry(iev);
 	 pos->Fill(x, y );
	 cout<<x<<" "<<y<<endl;
 }
pos->Draw()	 ;
 
}
