
#include "PTFTTree.hh"
#include <cmath>
#include <vector>
#include "TGraph.h"
#include "FindCircle.hpp"
#include "Hough.hpp"
#include "HoughDisplay.hpp"
#include "FindCircle.cpp"
#include "Hough.cpp"
#include "HoughDisplay.cpp"
void ptfplot( char * fname = "geant4ptf_100_data_format.root"){
	TTree *t =new TTree("t", "my tree");//Import the grid for the beam incidence
	Double_t x_pos,y_pos,z_pos; 
	t->Branch("x_pos", &x_pos,"x_pos/D");//Assign branch
	t->Branch("y_pos", &y_pos,"y_pos/D");	
	t->Branch("z_pos", &z_pos,"z_pos/D");  
    //TTree *t = new TTree("t", "my tree");
  	t->ReadFile("grid.txt","x/D:y:z" );
   
    // setup input
    TFile * fin = new TFile( fname, "read" );
    PTFTTree * ptf = new PTFTTree;  
    TTree * ptf_tree = (TTree*) fin->Get("ptf_tree");//Get the tree
    ptf->SetBranchAddresses( ptf_tree );
  
    // setup output
    TFile* fout = new TFile( "ptfplot.root", "recreate" );

	double min_bin_x=-300;//Set the bin the match data
	double min_bin_y=-300;
	double max_bin_y=300;
	double max_bin_x=300;
	int bin_x=63;
	int bin_y=59;

 // book histograms using the data
  TH2D * pos = new TH2D( "pos", "Detection efficiency; x (mm); y (mm)",
  		  			    bin_x,  min_bin_x, max_bin_x,  bin_y, min_bin_y, max_bin_y );//Detection efficiency in terms of incident beam	
						
  TH2D * reflection_inc_r = new TH2D( "reflection_inc ", "Beam position/hit position; radial position  (hit) (mm); radial position (incident) (mm)",
					    		  			    200, 0, 500, 200, 0, 500 );	//Radial reflection
  
  TH2D * reflection_inc = new TH2D( "reflection_inc ", "Beam position/hit position; z  (hit) (mm); radial position (incident) (mm)",
											  					    		  			    200, -250, 0, 200, 0, 500 );										
 
  TH2D * hq = new TH2D( "hq", "dig_Q variables; x (mm); y (mm)",
			  		  			      31, -300.0, 0, 31, -300.0, 0  );
  
  TH2D * reflection = new TH2D( "reflection", "Number of hits on photocathode;z (mm);time (nsec)",
			    961, -300,0 ,961, 0,10);
  
  TH2D * reflection_v = new TH2D( "reflection_v", "Number of hits on photocathode; time (nsec); y (mm)",
			  			       961, -1,5, 961, -300, 300.0 );
  
  TH2D * reflection_h = new TH2D( "reflection_h", "Number of hits on photocathode; x (mm); z (mm)",
						  			    200, -325, 325,200, -325, 325);
										
  TH2D * reflection_d = new TH2D( "reflection_d", "Number of hits on photocathode; radial distance (mm) ; z (mm)",
									  						  			    200, 0, 400,200, -150,0 );
							

  TH2D * hscanxy = new TH2D( "hscanxy", "Number of hits on photocathode(not same); x (mm); y (mm)",
			  		  			  bin_x,  min_bin_x, max_bin_x,  bin_y, min_bin_y, max_bin_y );
				
  TH2D * d_eff = new TH2D( "d_eff", "Detection efficiency; x (mm); y (mm)",
			  			    		  			   bin_x,  min_bin_x, max_bin_x,  bin_y,  min_bin_y, max_bin_y  );	
 TH2D * hqscan = new TH2D( "hqscan", "Number of hits on photocathode; x (mm); y (mm)",
			  		  			    bin_x,  min_bin_x,max_bin_x,  bin_y,  min_bin_y, max_bin_y  );
 
  //Testing the 
 TH1D * z_position = new TH1D( "z_position", "Z position; (cm)",
		   						  			  			    100,-254,0 );
  TH1D * dig_Time = new TH1D( "dig_Time", "time distribution; time(nsec)",
						  			  			    200,-300,300 );	
  TH2D * transit_time = new TH2D( "transit_time", "transit time (nsec); x (mm); y (mm)",bin_x,  min_bin_x, max_bin_x,  bin_y, min_bin_y, max_bin_y );	
  
  TH1D * TOF = new TH1D( "TOF", "True_used distribution; true_used",250,-1,5 );
   
   
   TH1D * deposit_energy = new TH1D( "deposit_energy", "Energy of the photon at the hit; ",250,3.00,3.15);
  																    	
 
  
   //Get the entries
  unsigned long long nevents = ptf_tree->GetEntries();
  cout<<nevents<<endl;
  
  vector< double > z_position_list; 
  vector< double > timing_list; 
  
  
  int i=0;
  int n=0;
  
  double qe_cut =0.0;
  double time_cut=0.0;
  for (unsigned long long iev =0 ; iev < nevents; ++iev ){
    
    ptf_tree->GetEvent( iev );
	
	//if (ptf->true_t[iev]>1);{ //Investigation around the value of what dig_T actually represents
	
		//continue;  
		//}
	
    // get average x, y position
    double xav = 0.0, yav = 0.0;
    int    nphotons = ptf->NPhotons;
	
	//double x=ptf->true_x[0];
	//double y=ptf->true_y[0];
	//double z=ptf->true_y[0];
	//TOF->Fill(ptf->true_t[0]-(-0.00337083*ptf->true_z[0]+0.0737));
	
	dig_Time->Fill(ptf->true_x[0]);
	
	reflection_d->Fill( sqrt(ptf->true_x[0]*ptf->true_x[0]+ptf->true_y[0]*ptf->true_y[0]), ptf->true_z[0]);
	
	//reflection_v->Fill( ptf->true_y[0], ptf->true_t[0]);
	reflection_h->Fill( ptf->true_y[0], ptf->true_z[0]);
	//hqscan->Fill(ptf->true_x[0],ptf->true_y[0]);
	//dig_Time->Fill(ptf->true_t[0]);
		
	bool k=ptf->true_used[0];
	if (k==0){ //We only want to display the photons that were detected
		qe_cut++;
		continue;
		};	
	//Timing the cut
	//if (ptf->true_t[0]>2.0){ //Get rid of the reflection that we had
		        //cout<<time_cut<<endl;//ptf->true_t[0]<<endl
				//cout<<ptf->true_used[0]<<endl;
	//			time_cut++;
	//				continue;
	//					};
		
	if (nphotons==1){
	TOF->Fill(ptf->true_used[0]);
	hqscan->Fill(ptf->true_x[0],ptf->true_y[0]);
	transit_time->Fill(ptf->true_x[0],ptf->true_y[0],ptf->true_t[0]);
	reflection->Fill(ptf->true_z[0],ptf->true_t[0]);
	reflection_v->Fill( ptf->true_t[0], ptf->true_y[0]);
						}
						
	
    for ( int iph = 0; iph<nphotons; ++iph ){
		
      xav += ptf->true_x[iph];
      yav += ptf->true_y[iph];
      hscanxy->Fill( ptf->true_x[iph], ptf->true_y[iph] );
      //hscanxy->Reset();//Why reset a plot ?
	 
		//cout<<ptf->true_used[iev] <<endl;
  	//	i++;
    //  dig_Time->Fill( ptf->true_x[iph], ptf->true_y[iph], ptf->dig_T[iph]); //
  }
    if ( nphotons > 0 ){
      xav /= nphotons;
      yav /= nphotons;
    }
    //hq->Fill( xav, yav, ptf->dig_Q[0] );
    int ix = iev / 300;//Switch by the max number of bins 
    int iy = iev % 300;
    
    double xloc = -300.0 + ix * 2.0; //Only a bin that is continued
    double yloc = -300.0 + iy * 2.0;
	
	
    //hqscan->SetBinContent( ix+1, iy+1, ptf->dig_Q[0] );
    //reflection->Fill( ptf->true_z[0], ptf->true_t[0]);

	
	deposit_energy->Fill(ptf->true_e[0]);
	z_position_list.push_back( ptf->true_z[0]);
	timing_list.push_back( ptf->true_t[0]);
	//hqscan->Fill(ptf->true_x[0],ptf->true_y[0]);

	hq->Fill( ptf->true_x[0], ptf->true_y[0]);

    //c.Update();
    //char gif[100];
    //sprintf(gif,"gifsicle/hscanxy%05d.gif",iev);
    //c.SaveAs(gif);
        
    //cout << "iev="<<iev<<" ix="<<ix<<" iy="<<iy
	 //<<" xloc="<<xloc<<" yloc="<<yloc
	 //<<" x="<<xav<<" y="<<yav<<" Q="<<ptf->dig_Q[0]<<endl;
	//cout<<ptf->true_t[iev]<<endl;
  }
  
  
  
  
  
  

long double ratio;
ratio=((double) n)/nevents;
cout<<"The ratio is:"<<ratio<<" "<<n<<"  "<<nevents<<endl;
TGraph* position_correlation=new TGraph(z_position_list.size(),&z_position_list[0],&timing_list[0]);
position_correlation->GetXaxis()->SetTitle("Z_position");//Create normal plot for all the corrected efficiency                                                                     
position_correlation->GetYaxis()->SetTitle("Time");
position_correlation->SetTitle("Correlation between position and time");
position_correlation->SetMarkerStyle(2);
position_correlation->SetMarkerColor(2);
position_correlation->Fit("pol1");


/*
Calculation of the slope for the plot

*/
//Int_t MaxBin=hqscan->GetMaximumBin();
//Int_t x,y,z;
//hqscan->GetBinXYZ(MaxBin, x, y, z);
//int binmax = z_position->GetMaximumBin();
//int binmax = z_position->FindLastBinAbove(1000);
/*double z_max_value=-1;
double timing_max_value=0.1;
double z_min_value=-10;
double timing_min_value=5.0;
for (int i=0;i<=hqscan->GetNbinsX();i++){
	  if( hqscan->GetBinContent(i)<1e-10 ){
		  continue ;
	  };
	  if (z_max_value>hqscan->GetXaxis()->GetBinCenter(i)){
		  z_max_value=z_position->GetXaxis()->GetBinCenter(i);
	  };
      if (z_min_value<z_position->GetXaxis()->GetBinCenter(i)){
 		  z_min_value=z_position->GetXaxis()->GetBinCenter(i);
 	  };
	  
  }
for (int i=0;i<=timing->GetNbinsX();i++){
	   //cout<<timing->GetXaxis()->GetBinCenter(i)<<timing->GetNbinsX()<<endl;

	  if( timing->GetBinContent(i)<1e-10 ){
		  continue ;
	  }
	 
	  bool k;
	  	 
	  if (timing_min_value>timing->GetXaxis()->GetBinCenter(i)){
		  timing_min_value=timing->GetXaxis()->GetBinCenter(i);
 		 //cout<<timing_min_value<<endl;
	  };
      if (timing_max_value<timing->GetXaxis()->GetBinCenter(i)){
		  	  timing_max_value=timing->GetXaxis()->GetBinCenter(i);
			  	  }; 
}

*/
double ratio_de=0;
int weight;
double cut_timing;
int h=0;
for (unsigned long long iev =0 ; iev < nevents; ++iev ){
	   ptf_tree->GetEvent( iev );
	   
 	weight=ptf->true_used[0];
	cut_timing=ptf->true_t[0];
	
	//Apply timing cut
	//if (ptf->true_t[0]>1.0){ //Get rid of the reflection that we had
			//cout<<ptf->true_t[0]<<endl;//ptf->true_t[0]<<endl
			//cout<<ptf->true_used[0]<<endl;
	//			weight=0;
		//cut_timing=-1.0;
	//		}
	
							//};
    //t->GetEvent(h);	.q
	//reflection_inc->Fill(ptf->true_z[0],sqrt(x_pos*x_pos+y_pos*y_pos));	
	//reflection_inc->Fill(sqrt(x_pos*x_pos+y_pos*y_pos),z_pos);	
	//reflection_inc_r->Fill(sqrt(x_pos*x_pos+y_pos*y_pos),sqrt(ptf->true_x[0]*ptf->true_x[0]+ptf->true_y[0]*ptf->true_y[0]));
			
	
	if (iev%100==0 & iev!=0){
	t->GetEvent(iev/100);
	pos->Fill(x_pos,y_pos,((double) ratio_de)/200);	
	h++;
	

	//TOF->Fill(weight);
	//cout<<((double) ratio_de)<<"  "<<x_pos<<"  "<<y_pos<<" "<< ptf->true_used[0]<<"  "<<iev<<endl;
	
	ratio_de=0.0;
//	 cout<<x_pos<<" "<<y_pos<<endl;
	 
	 }
	 	
		
		
    int    nphotons = ptf->NPhotons;
 	if (nphotons==1){
 	 ratio_de=weight+ratio_de;
 }
	
	 //cout<< ratio_de<<" "<<iev<<" "<<ptf->true_used[0]<<endl;
 }
// for (unsigned long long iev =0 ; iev < 10001; ++iev )


 double max_x_data=0.7;//0.6385;//Input coming from the data !
 double min_x_data=0.1;//0.1465;
 double max_y_data=0.7;//0.6185;
 double min_y_data=0.1;//0.1265;
 double x_bin_2;
 double y_bin_2;
 //int bin_sim=31;
 
 TH2D *pos_rescale = new TH2D("pos_rescale","Simulated DE",bin_x,0.1,0.7,bin_y,0.1,0.7);
 //Rescaling the PTF simulation
 for (int i=0;i<=pos->GetNbinsX();i++){
 	x_bin_2=pos->GetXaxis()->GetBinCenter(i);
 	for (int j=0;j<=pos->GetNbinsX();j++){
 	y_bin_2=pos->GetYaxis()->GetBinCenter(j);
 	//if (-250<=round(x_bin_2*100)/100 &&round(x_bin_2*100)/100<=250 && -250<=round(y_bin_2*100)/100 && round(y_bin_2*100)/100<=250){
 		if (x_bin_2<0 ){
 			if (y_bin_2<0){
 				//cout<<pos->GetBinContent(i,j)<<" "<<y_bin_2<<" "<<endl;
 				 pos_rescale->Fill(min_x_data+(x_bin_2+max_bin_x)/1000,min_y_data+(y_bin_2+max_bin_y)/1000,pos->GetBinContent(i,j));
 		}
 			if (y_bin_2>0){
 				//cout<<x_bin_2<<" "<<y_bin_2<<" "<<qe_sim->GetBinContent(i,j)<<endl;
 				 pos_rescale->Fill(min_x_data+(x_bin_2+max_bin_x)/1000,(min_y_data+max_y_data)/2+(y_bin_2)/1000,pos->GetBinContent(i,j));
				
 		}
 	}

 		if (x_bin_2>0 ){
 			if (y_bin_2<0){
 				//cout<<x_bin_2<<" "<<y_bin_2<<" "<<qe_sim->GetBinContent(i,j)<<endl;
 				 pos_rescale->Fill((min_x_data+max_x_data)/2+(x_bin_2)/1000,min_y_data+(y_bin_2+max_bin_y)/1000,pos->GetBinContent(i,j));
				
 			} 
		
 			if (y_bin_2>0){
 				//cout<<x_bin_2<<" "<<y_bin_2<<" "<<qe_sim->GetBinContent(i,j)<<endl;
 				 pos_rescale->Fill((min_x_data+max_x_data)/2+(x_bin_2)/1000,(min_y_data+max_y_data)/2+(y_bin_2)/1000,pos->GetBinContent(i,j));
				
 		}
 	}
 
		
 	}
 }
 
 
 
 
 
cout<<"The timing and detection cut are "<< ((double)time_cut)/nevents<<"   "<<((double)qe_cut)/nevents<<endl;


//Calculation of the detection efficiency

 for(int ix=1; ix<=hqscan->GetNbinsX(); ix++){
	 for(int iy=1; iy<= hqscan->GetNbinsY(); iy++){
		   double weight= ( hqscan->GetBinContent(ix, iy))/200;//Comes from the PMT factor for more statistics
		   if (hqscan->GetBinContent(ix, iy)==0){
		   transit_time->SetBinContent(ix,iy,0);
		}
		else{
		   	double weight_transit_time= ( transit_time->GetBinContent(ix, iy)/hqscan->GetBinContent(ix, iy));
			transit_time->SetBinContent(ix,iy,weight_transit_time);
		}
		   
		   
		    d_eff->SetBinContent(ix,iy,weight);
			//cout<<weight_transit_time<<endl;
	  
		} 
	}
 
 
TH2D* pmt0_qe_corr_grad; //Apply the hough transformation to get something as similar to the data
Circle_st circ = find_circle_max_grad( pos_rescale, pmt0_qe_corr_grad, 0.5 );
zero_outside_circle( pos_rescale, circ );
	
//Save the histogram made
  TCanvas* c = new TCanvas("canvas");
  string plotname;
  hqscan->SetDirectory( fout );
  hqscan->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/hqscan.pdf");
  
  c->SaveAs(plotname.c_str(),"pdf");
  hscanxy->SetStats(0);//Made to not display statistic box
  hscanxy->SetDirectory( fout );
  hscanxy->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/hqscanxy.pdf");
  
  
  pos_rescale->SetStats(0);//Made to not display statistic box
  pos_rescale->SetDirectory( fout );
  //pos->SetMinimum(0.2);
  pos_rescale->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/gun_position.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  
 
  //timing->SetStats(0);//Made to not display statistic box
  position_correlation->Draw("AP");
  gStyle->SetOptFit() ;
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/z_position.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  
  reflection->SetStats(0);
  reflection->SetDirectory( fout );
  reflection->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  reflection_v->SetStats(0);
  reflection_v->SetDirectory( fout );
  reflection_v->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection_vertical.pdf");
  c->SaveAs(plotname.c_str(),"pdf");

  reflection_h->SetStats(0);
  reflection_h->SetDirectory( fout );
  reflection_h->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection_horizontal.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  reflection_d->SetStats(0);
  reflection_d->SetDirectory( fout );
  reflection_d->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection_distance.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  

   reflection_inc->SetStats(0);
   reflection_inc->SetDirectory( fout );
   reflection_inc->Draw("colz0");
   gPad->Modified();
   gPad->Update();
   plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection_incident.pdf");
   c->SaveAs(plotname.c_str(),"pdf");
   
   reflection_inc_r->SetStats(0);
   reflection_inc_r->SetDirectory( fout );
   reflection_inc_r->Draw("colz0");
   gPad->Modified();
   gPad->Update();
   plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/reflection_incident_radial.pdf");
   c->SaveAs(plotname.c_str(),"pdf");
	 
  d_eff->SetStats(0);//Made to not display statistic box
  d_eff->SetDirectory( fout );
  //d_eff->SetMaximum(1.0);
  d_eff->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/DE.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  transit_time->SetStats(0);//Made to not display statistic box
  transit_time->SetDirectory( fout );
  //transit_time->SetMaximum(2.0);
  //transit_time->SetMinimum(0.0);
  transit_time->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/transit_time.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  deposit_energy->SetStats(0);//Made to not display statistic box
  deposit_energy->SetDirectory( fout );
  //d_eff->SetMinimum(0.2);
  deposit_energy->Draw();
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/E_deposit.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  
  /*
  //d_eff_norm->SetMinimum(0.3);
  d_eff_norm->SetMaximum(1);
  d_eff_norm->SetStats(0);//Made to not display statistic box
  d_eff_norm->SetDirectory( fout );
  d_eff_norm->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/De_test_5.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  */

  //TOF->SetMinimum(0.0);
  //dig_Time->SetStats(0);//Made to not display statistic box
  TOF->SetDirectory( fout );
  //c->SetLogy();
  //TOF->SetFillColor( kRed);
  //SetFillStyle( 3001);
  TOF->Draw();
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/TOF.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  //TOF->SetMinimum(0.0);
  //dig_Time->SetStats(0);//Made to not display statistic box
  dig_Time->SetDirectory( fout );
  //c->SetLogy();
  //TOF->SetFillColor( kRed);
  //SetFillStyle( 3001);
  dig_Time->Draw();
  gPad->Modified();
  gPad->Update();
  plotname = string("/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/postprocess/time_distribution.pdf");
  c->SaveAs(plotname.c_str(),"pdf");
  
  
  
  fin->Close();
  
  fout->Write();
  fout->Close();

}