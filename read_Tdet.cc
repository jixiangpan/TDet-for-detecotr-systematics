#include<iostream>
#include<fstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"
#include "TColor.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

/// minuit2
//#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "draw.icc"

///////////////////////////////////////////////////////////////////////////////////////////////////////// TDet
//
// Bootstrapping Method: random sampling with replacement
//
//

class TDet
{
public:
  TDet() {
    num_DetVar = 0;
  }

  void Clear();

  void Exe(TString file_CV, TString file_Var, bool flag_numu, bool flag_FC, int ntoy);
  
  //////////////////////////

  int num_DetVar;
  
  TMatrixD matrix_cov_on_absdiff;

  TH1D *h1_CV_noSel;
  TH1D *h1_CV_wiSel;
  TH1D *h1_Var_noSel;
  TH1D *h1_Var_wiSel;
  
  TH1D *h1_weight;

  map<int, TString>map_total_index_string;
  map<int, double>map_total_index_Erec;
  map<TString, double>map_CV_wiSel_Erec;
  map<TString, double>map_Var_wiSel_Erec;
};


void TDet::Exe(TString file_CV, TString file_Var, bool flag_numu, bool flag_FC, int ntoy)
{  
  TString roostr = "";

  int entry(0), run(0), subrun(0), event(0);
  double weight(0);
  int status(0), post_generic(0), numuCC(0);
  double nueBDT(0);
  int is_FC(0);
  double Erec(0), Evis(0);

  /////////////////////////////////////////////// CV
  
  int line_CV_a0 = 0;
  ifstream InputFile_CV_a0(file_CV, ios::in);
  if(!InputFile_CV_a0) { cerr<<" No input-list"<<endl; exit(1); }    
  while( !InputFile_CV_a0.eof() ) {
    line_CV_a0++;
    InputFile_CV_a0 >> entry >> run >> subrun >> event
		    >> weight >> status >> post_generic >> numuCC
		    >> nueBDT >> is_FC
		    >> Erec >> Evis;    
  }
  line_CV_a0 -= 1;
  cout<<TString::Format(" ---> number of events (CV): %4d", line_CV_a0)<<endl;
     
  /////////////////////////////////////////////// Var
  
  int line_Var_a0 = 0;
  ifstream InputFile_Var_a0(file_Var, ios::in);
  if(!InputFile_Var_a0) { cerr<<" No input-list"<<endl; exit(1); }    
  while( !InputFile_Var_a0.eof() ) {
    line_Var_a0++;
    InputFile_Var_a0 >> entry >> run >> subrun >> event
		    >> weight >> status >> post_generic >> numuCC
		    >> nueBDT >> is_FC
		    >> Erec >> Evis;    
  }
  line_Var_a0 -= 1;
  cout<<TString::Format(" ---> number of events (DetVar): %4d", line_Var_a0)<<endl;
     
  /////////////////////////////////////////////// set event weight
  
  if( line_CV_a0 != line_Var_a0 ) { cerr<<" Error: not common samples between CV and Var"<<endl; exit(1); }

  roostr = "h1_weight";
  h1_weight = new TH1D(roostr, roostr, line_CV_a0, 0.5, line_CV_a0+0.5);
  
  ifstream InputFile_weight(file_CV, ios::in);
  for(int idx=1; idx<=line_CV_a0; idx++) {
    
    InputFile_weight >> entry >> run >> subrun >> event
		     >> weight >> status >> post_generic >> numuCC
		     >> nueBDT >> is_FC >> Erec >> Evis;
    
    h1_weight->SetBinContent( idx, weight );

    roostr = TString::Format("%d_%d_%d", run, subrun, event);
    map_total_index_string[idx] = roostr;
    map_total_index_Erec[idx] = Erec*1000;
      
    // if( weight<0.1 ) {
    //   cout<<" check weight: "
    // 	<< idx <<"\t"<< run <<"\t"<< subrun <<"\t"<< event
    // 	<<"\t"<< weight <<"\t"<< status <<"\t"<< post_generic <<"\t"<< numuCC
    // 	<<"\t"<< nueBDT <<"\t"<< is_FC <<"\t"
    // 	<< Erec <<"\t"<< Evis
    // 	<<endl;
    // }
    
  }// idx

  int total_number_weighted = (int)( h1_weight->Integral() );
  cout<<TString::Format(" ---> Total number of events (weighted): %d", total_number_weighted)<<endl;
  
  ///////////////////////////////////////////////

  int bins_basic = 15;
  double low_basic = 100;
  double hgh_basic = 2350;  
  roostr = "h1_CV_noSel"; h1_CV_noSel = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  roostr = "h1_CV_wiSel"; h1_CV_wiSel = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);  
  roostr = "h1_Var_noSel"; h1_Var_noSel = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  roostr = "h1_Var_wiSel"; h1_Var_wiSel = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);

  //////////////////////// CV
  
  ifstream InputFile_CV_a1(file_CV, ios::in);
  for(int idx=1; idx<=line_CV_a0; idx++) {
    
    InputFile_CV_a1 >> entry >> run >> subrun >> event
		    >> weight >> status >> post_generic >> numuCC >> nueBDT >> is_FC >> Erec >> Evis;
    
    Erec *= 1000;
    roostr = TString::Format("%d_%d_%d", run, subrun, event);

    h1_CV_noSel->Fill( Erec, weight );
      
    if( is_FC==flag_FC ) {

      if( flag_numu ) {// numu candidates
	if( status==1 && post_generic==1 && numuCC==1 ) {
	  h1_CV_wiSel->Fill( Erec, weight );
	  map_CV_wiSel_Erec[roostr] = Erec;
	}
      }
      else {// nue candidates
	if( status==1 && post_generic==1 && nueBDT>6 ) {
	  h1_CV_wiSel->Fill( Erec, weight );
	  map_CV_wiSel_Erec[roostr] = Erec;
	}
      }
      
    }// flag_FC
    
  }// idx
  
  //////////////////////// Var
  
  ifstream InputFile_Var_a1(file_Var, ios::in);
  for(int idx=1; idx<=line_Var_a0; idx++) {
    
    InputFile_Var_a1 >> entry >> run >> subrun >> event
		    >> weight >> status >> post_generic >> numuCC >> nueBDT >> is_FC >> Erec >> Evis;
    
    Erec *= 1000;
    roostr = TString::Format("%d_%d_%d", run, subrun, event);

    h1_Var_noSel->Fill( Erec, weight );
      
    if( is_FC==flag_FC ) {

      if( flag_numu ) {// numu candidates
	if( status==1 && post_generic==1 && numuCC==1 ) {
	  h1_Var_wiSel->Fill( Erec, weight );
	  map_Var_wiSel_Erec[roostr] = Erec;
	}
      }
      else {// nue candidates
	if( status==1 && post_generic==1 && nueBDT>6 ) {
	  h1_Var_wiSel->Fill( Erec, weight );
	  map_Var_wiSel_Erec[roostr] = Erec;
	}
      }
      
    }// flag_FC
    
  }// idx
  
  ///////////////////////////////////////////////

  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = 0;

    content = h1_CV_wiSel->GetBinContent(ibin);
    content = ( (int)(content*10+0.5) )*1./10;
    h1_CV_wiSel->SetBinContent(ibin, content);
    
    content = h1_Var_wiSel->GetBinContent(ibin);
    content = ( (int)(content*10+0.5) )*1./10;
    h1_Var_wiSel->SetBinContent(ibin, content);    
  }

  /////////////////////////////////////////////// Calculate covariance matrix

  TPrincipal principal_test(bins_basic, "ND");
  double *array_test = new double[bins_basic];

  
  for(int itoy=1; itoy<=ntoy; itoy++) {

    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy: %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
    TRandom3 *rand3 = new TRandom3(0);

    for(int ibin=0; ibin<bins_basic; ibin++) array_test[ibin] = 0;
    
    for(int idx=1; idx<=total_number_weighted; idx++) {

      double random = h1_weight->GetRandom();
      int global_index = h1_weight->FindBin( random );
      roostr = map_total_index_string[global_index];
      double Erec = map_total_index_Erec[global_index];
      int bin_index = h1_CV_noSel->FindBin( Erec );

      if( bin_index<1 || bin_index>bins_basic ) {
	//cout<<TString::Format(" ---> check: random %10.2f, global_idx %10d, Erec %8.2f, bin_index %2d", random, global_index, Erec, bin_index )<<endl;
	continue;
      }      
      //cout<<TString::Format(" ---> check: random %10.2f, global_idx %10d, Erec %8.2f, bin_index %2d", random, global_index, Erec, bin_index )<<endl;
      
      int CV_count = 0;
      int Var_count = 0;
      if( map_CV_wiSel_Erec.find(roostr)!=map_CV_wiSel_Erec.end() ) CV_count++;
      if( map_Var_wiSel_Erec.find(roostr)!=map_Var_wiSel_Erec.end() ) Var_count++;

      array_test[bin_index-1] += (Var_count-CV_count);
      
    }// idx

    principal_test.AddRow( array_test );

    delete rand3;
  }// itoy

  
  TMatrixD *tt_matrix_cov_on_absdiff = (TMatrixD *)principal_test.GetCovarianceMatrix();
  int rows = tt_matrix_cov_on_absdiff->GetNrows();
  matrix_cov_on_absdiff.ResizeTo(rows, rows);
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_on_absdiff)(i,j) = (*tt_matrix_cov_on_absdiff)(j,i);
      matrix_cov_on_absdiff(i,j) = (*tt_matrix_cov_on_absdiff)(i,j);
    }
  }

  ////////////////////////
  
  delete[] array_test;
  
}

void TDet::Clear()
{
  num_DetVar++;

  map_total_index_string.clear();
  map_total_index_Erec.clear();  
  map_CV_wiSel_Erec.clear();
  map_Var_wiSel_Erec.clear();
  
  if( num_DetVar==1 ) {
    h1_CV_noSel  = NULL;
    h1_CV_wiSel  = NULL;
    h1_Var_noSel = NULL;
    h1_Var_wiSel = NULL;
    
    h1_weight = NULL;    
  }
  else {      
    delete h1_CV_noSel;
    delete h1_CV_wiSel;
    delete h1_Var_noSel;
    delete h1_Var_wiSel;
    
    delete h1_weight;
  }
  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_Tdet()
{
  TString roostr = "";

  ////////////////////////////////////// Draw style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetEndErrorSize(4);

  //////////////////////////////////////////////////////////////////////////////////

  //
  TDet *det_test = new TDet();

  det_test->Clear();
  
  det_test->Exe("./data_det_syst_total/common_detvar/LYDown/numu_CV", "./data_det_syst_total/common_detvar/LYDown/numu_Var", 1, 1, 100);
  //det_test->Exe("./data_det_syst_total/common_detvar/Recomb2/numu_CV", "./data_det_syst_total/common_detvar/Recomb2/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/SCE/numu_CV", "./data_det_syst_total/common_detvar/SCE/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModThetaXZ/numu_CV", "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModThetaYZ/numu_CV", "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModX/numu_CV", "./data_det_syst_total/common_detvar/WireModX/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModYZ/numu_CV", "./data_det_syst_total/common_detvar/WireModYZ/numu_Var", 1, 1, 1000);

  
  ////////////////////////////////////////////////////////////////////////////////// plotting

  roostr = "canv_h1_weight";
  TCanvas *canv_h1_weight = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_weight, 0.15, 0.2,0.1,0.15);
  det_test->h1_weight->Draw();
  det_test->h1_weight->SetLineColor(kBlack);
  //h1_weight->SetLineStyle(7);
  func_title_size(det_test->h1_weight, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(det_test->h1_weight, "Event index", "Weight");
  det_test->h1_weight->GetXaxis()->SetNdivisions(506);

  ///////////
  
  roostr = "canv_spectra";
  TCanvas *canv_spectra = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_spectra, 0.15, 0.2,0.1,0.15);
  canv_spectra->SetLogy();
  
  det_test->h1_CV_noSel->Draw("hist");
  det_test->h1_CV_noSel->SetMinimum(0.9);
  det_test->h1_CV_noSel->SetTitle("");
  det_test->h1_CV_noSel->SetLineColor(kBlack);
  det_test->h1_CV_noSel->SetLineStyle(7);
  func_title_size(det_test->h1_CV_noSel, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(det_test->h1_CV_noSel, "E_{rec} (MeV) ", "Entries");
  det_test->h1_CV_noSel->GetXaxis()->SetNdivisions(506);
  
  det_test->h1_Var_noSel->Draw("same hist");
  det_test->h1_Var_noSel->SetLineColor(kRed);
  det_test->h1_Var_noSel->SetLineStyle(7);

  det_test->h1_CV_wiSel->Draw("same hist text75");
  det_test->h1_CV_wiSel->SetLineColor(kBlack);
  det_test->h1_CV_wiSel->SetMarkerColor(kBlack);
  
  det_test->h1_Var_wiSel->Draw("same hist text15");
  det_test->h1_Var_wiSel->SetLineColor(kRed);
  det_test->h1_Var_wiSel->SetMarkerColor(kRed);

  det_test->h1_CV_noSel->Draw("same axis");

  ///////////
  
  roostr = "canv_matrix_cov_on_absdiff";
  TCanvas *canv_matrix_cov_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_cov_on_absdiff, 0.15, 0.2,0.1,0.15);
  det_test->matrix_cov_on_absdiff.Draw("colz");   

  int bins_temp = det_test->h1_CV_noSel->GetNbinsX();
  double low_temp = det_test->h1_CV_noSel->GetXaxis()->GetBinLowEdge(1);
  double hgh_temp = det_test->h1_CV_noSel->GetXaxis()->GetBinUpEdge(bins_temp);

  //////
  roostr = "h2_cov_on_absdiff";
  TH2D *h2_cov_on_absdiff = new TH2D(roostr, roostr, bins_temp, low_temp, hgh_temp, bins_temp, low_temp, hgh_temp);
  for(int i=0; i<bins_temp; i++) {
    for(int j=0; j<bins_temp; j++) {
      double content = det_test->matrix_cov_on_absdiff(i,j);
      content = ((int)(content*10 + 0.5))*1./10;
      h2_cov_on_absdiff->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv_h2_cov_on_absdiff";
  TCanvas *canv_h2_cov_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_cov_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_cov_on_absdiff->Draw("colz text");
  h2_cov_on_absdiff->SetTitle("");
  func_title_size(h2_cov_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_cov_on_absdiff, "E_{rec} (MeV) ", "E_{rec} (MeV)");
  h2_cov_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_cov_on_absdiff->GetYaxis()->SetNdivisions(506);

  //////
  roostr = "h1_cov2absdiff";
  TH1D *h1_cov2absdiff = new TH1D(roostr, roostr, bins_temp, low_temp, hgh_temp);
  for(int i=0; i<bins_temp; i++) {
    double val_CV = det_test->h1_CV_wiSel->GetBinContent(i+1);
    double val_Var = det_test->h1_Var_wiSel->GetBinContent(i+1);
    double cov_root = sqrt( det_test->matrix_cov_on_absdiff(i,i) );
    double reldiff = cov_root/(val_Var-val_CV);
    h1_cov2absdiff->SetBinContent( i+1, fabs(reldiff) );
  }
  roostr = "canv_h1_cov2absdiff";
  TCanvas *canv_h1_cov2absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_cov2absdiff, 0.15, 0.2,0.1,0.15);
  h1_cov2absdiff->Draw("hist");
  h1_cov2absdiff->SetTitle("");
  func_title_size(h1_cov2absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_cov2absdiff, "E_{rec} (MeV) ", "|Sqrt(cov) / (DetVar-CV)|");
  h1_cov2absdiff->GetXaxis()->SetNdivisions(506);
  h1_cov2absdiff->GetYaxis()->SetNdivisions(506);
  
  //////
  roostr = "h2_correlation_on_absdiff";
  TH2D *h2_correlation_on_absdiff = new TH2D(roostr, roostr, bins_temp, low_temp, hgh_temp, bins_temp, low_temp, hgh_temp);
  for(int i=0; i<bins_temp; i++) {
    for(int j=0; j<bins_temp; j++) {
      double content = det_test->matrix_cov_on_absdiff(i,j);
      double val_i = sqrt( det_test->matrix_cov_on_absdiff(i,i) );
      double val_j = sqrt( det_test->matrix_cov_on_absdiff(j,j) );
      content = content/val_i/val_j;
      content = ((int)(content*100 + 0.5))*1./100;
      h2_correlation_on_absdiff->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv_h2_correlation_on_absdiff";
  TCanvas *canv_h2_correlation_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_correlation_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_correlation_on_absdiff->Draw("colz text");
  h2_correlation_on_absdiff->SetTitle("");
  func_title_size(h2_correlation_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_correlation_on_absdiff, "E_{rec} (MeV) ", "E_{rec} (MeV)");
  h2_correlation_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_correlation_on_absdiff->GetYaxis()->SetNdivisions(506);

  /////////////////////////////////////////

  TCanvas *canv_sum = new TCanvas("canv_sum", "canv_sum", 1400, 1000);
  canv_sum->Divide(2,2);

  TVirtualPad *pad_spectra = canv_sum->cd(1);
  pad_spectra->SetLogy();
  func_canv_margin(pad_spectra, 0.15, 0.2,0.1,0.15);  
  det_test->h1_CV_noSel->Draw("hist");  
  det_test->h1_Var_noSel->Draw("same hist");
  det_test->h1_CV_wiSel->Draw("same hist text75");
  det_test->h1_Var_wiSel->Draw("same hist text15");
  det_test->h1_CV_noSel->Draw("same axis");

  TVirtualPad *pad_cov = canv_sum->cd(2);
  func_canv_margin(pad_cov, 0.15, 0.2,0.1,0.15);
  h2_cov_on_absdiff->Draw("colz text");

  TVirtualPad *pad_correlation = canv_sum->cd(3);
  func_canv_margin(pad_correlation, 0.15, 0.2,0.1,0.15);
  h2_correlation_on_absdiff->Draw("colz text");
  
  TVirtualPad *pad_cov2absdiff = canv_sum->cd(4);
  func_canv_margin(pad_cov2absdiff, 0.15, 0.2,0.1,0.15);
  h1_cov2absdiff->Draw("hist");

  canv_sum->SaveAs("canv_sum.png");
  
}
