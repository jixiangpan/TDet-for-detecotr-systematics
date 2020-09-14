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
  TH1D *h1_diff_Var2CV;
  TH1D *h1_binomial_err;
  
  TH1D *h1_weight;
  TH1D *h1_sampling;
  TH1D *h1_diff_weight;
  
  map<int, TString>map_total_index_string;
  map<TString, double>map_total_index_CV_Erec;
  map<TString, double>map_total_index_Var_Erec;
  map<TString, double>map_CV_wiSel_Erec;
  map<TString, double>map_Var_wiSel_Erec;

  map<int, TMatrixD>map_matrix_cov_on_absdiff;
  map<int, TH1D*>map_h1_CV_noSel;
  map<int, TH1D*>map_h1_CV_wiSel;
  map<int, TH1D*>map_h1_Var_noSel;
  map<int, TH1D*>map_h1_Var_wiSel;
  map<int, TH1D*>map_h1_diff_Var2CV;
  map<int, TH1D*>map_h1_binomial_err;
  
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
  roostr = "h1_diff_weight";
  h1_diff_weight = new TH1D(roostr, roostr, line_CV_a0, 0.5, line_CV_a0+0.5);

  //// CV
  ifstream InputFile_CV_weight(file_CV, ios::in);
  for(int idx=1; idx<=line_CV_a0; idx++) {
    
    InputFile_CV_weight >> entry >> run >> subrun >> event
			>> weight >> status >> post_generic >> numuCC
			>> nueBDT >> is_FC >> Erec >> Evis;
    
    h1_weight->SetBinContent( idx, weight );

    roostr = TString::Format("%d_%d_%d", run, subrun, event);
    map_total_index_string[idx] = roostr;
    map_total_index_CV_Erec[roostr] = Erec*1000;      
  }// idx

  //// Var
  ifstream InputFile_Var_weight(file_Var, ios::in);
  for(int idx=1; idx<=line_Var_a0; idx++) {
    
    InputFile_Var_weight >> entry >> run >> subrun >> event
			 >> weight >> status >> post_generic >> numuCC
			 >> nueBDT >> is_FC >> Erec >> Evis;

    roostr = TString::Format("%d_%d_%d", run, subrun, event);
    map_total_index_Var_Erec[roostr] = Erec*1000;      
  }// idx

  ////
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
  roostr = "h1_diff_Var2CV"; h1_diff_Var2CV = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  roostr = "h1_binomial_err"; h1_binomial_err = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
    
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

    double content_noSel_CV = 0;
    content_noSel_CV = h1_CV_noSel->GetBinContent(ibin);
    content_noSel_CV = ( (int)(content_noSel_CV*100+0.5) )*1./100;
    h1_CV_noSel->SetBinContent(ibin, content_noSel_CV);
    
    double content_noSel_Var = 0;
    content_noSel_Var = h1_Var_noSel->GetBinContent(ibin);
    content_noSel_Var = ( (int)(content_noSel_Var*100+0.5) )*1./100;
    h1_Var_noSel->SetBinContent(ibin, content_noSel_Var);
    
    double content_wiSel_CV = 0;
    content_wiSel_CV = h1_CV_wiSel->GetBinContent(ibin);
    content_wiSel_CV = ( (int)(content_wiSel_CV*100+0.5) )*1./100;
    h1_CV_wiSel->SetBinContent(ibin, content_wiSel_CV);
    
    double content_wiSel_Var = 0;
    content_wiSel_Var = h1_Var_wiSel->GetBinContent(ibin);
    content_wiSel_Var = ( (int)(content_wiSel_Var*100+0.5) )*1./100;
    h1_Var_wiSel->SetBinContent(ibin, content_wiSel_Var);
    
    h1_diff_Var2CV->SetBinContent( ibin, fabs(content_wiSel_Var-content_wiSel_CV) );

    // double val_p = content_wiSel_CV/content_noSel_CV;
    // double val_err = sqrt(  content_noSel_CV*val_p*(1-val_p) );
    // val_err = ( (int)(val_err*100+0.5) )*1./100;
    // h1_binomial_err->SetBinContent(ibin, val_err);

    double sampleA_total = 0;
    double sampleA_sub   = 0;
    if( content_wiSel_CV>content_wiSel_Var ) {
      sampleA_total = content_wiSel_CV;
      sampleA_sub = (content_wiSel_CV - content_wiSel_Var);
    }
    else {
      sampleA_total = content_wiSel_Var;
      sampleA_sub = (content_wiSel_Var - content_wiSel_CV);
    }
    double val_p = sampleA_sub/sampleA_total;// probility of the non-selected events
    double val_err = sqrt( sampleA_total*val_p*(1-val_p) );// binominal error of the non-selected events
    val_err = ( (int)(val_err*100+0.5) )*1./100;
    h1_binomial_err->SetBinContent(ibin, val_err);

    
  }
  
  /////////////////////////////////////////////// Calculate covariance matrix

  TPrincipal principal_test(bins_basic, "ND");
  double *array_test = new double[bins_basic];
  double *array_CV = new double[bins_basic];
  double *array_Var = new double[bins_basic];

  /// to check
  roostr = "h1_sampling";
  h1_sampling = new TH1D(roostr, roostr, line_CV_a0, 0.5, line_CV_a0+0.5);

  double *array_mean_Var2CV = new double[bins_basic];  
  double *array_total_CV_wiSel = new double[bins_basic];
  double *array_total_Var_wiSel = new double[bins_basic];
  for(int ibin=0; ibin<bins_basic; ibin++) {
    array_mean_Var2CV[ibin] = 0;
    array_total_CV_wiSel[ibin] = 0;
    array_total_Var_wiSel[ibin] = 0;
  }
  
  ///
  for(int itoy=1; itoy<=ntoy; itoy++) {

    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy: %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
    TRandom3 *rand3 = new TRandom3(0);

    for(int ibin=0; ibin<bins_basic; ibin++) {
      array_test[ibin] = 0;
      array_CV[ibin]   = 0;
      array_Var[ibin]  = 0;
    }
    
    for(int idx=1; idx<=total_number_weighted; idx++) {

      /////////////////////////////////////////// common samples
      double random = h1_weight->GetRandom();
      int global_index = h1_weight->FindBin( random );
      roostr = map_total_index_string[global_index];
      h1_sampling->Fill( global_index );
      
      double CV_Erec = map_total_index_CV_Erec[roostr];
      int CV_bin_index = h1_CV_noSel->FindBin( CV_Erec );            
      if( CV_bin_index>=1 && CV_bin_index<=bins_basic ) {
	int CV_count = 0;
	if( map_CV_wiSel_Erec.find(roostr)!=map_CV_wiSel_Erec.end() ) CV_count++;
	array_CV[CV_bin_index-1] += CV_count;
	array_total_CV_wiSel[CV_bin_index-1] += CV_count;
      }      
      
      double Var_Erec = map_total_index_Var_Erec[roostr];
      int Var_bin_index = h1_Var_noSel->FindBin( Var_Erec );            
      if( Var_bin_index>=1 && Var_bin_index<=bins_basic ) {
	int Var_count = 0;
	if( map_Var_wiSel_Erec.find(roostr)!=map_Var_wiSel_Erec.end() ) Var_count++;
	array_Var[Var_bin_index-1] += Var_count;
	array_total_Var_wiSel[Var_bin_index-1] += Var_count;
      }            
      
    }// idx

    for(int ibin=0; ibin<bins_basic; ibin++) {
      array_test[ibin] = array_Var[ibin] - array_CV[ibin];
      
      array_mean_Var2CV[ibin] += ( array_Var[ibin] - array_CV[ibin] )/ntoy;
    }
    
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

  // TVectorD *vector_mean = (TVectorD*)principal_test.GetMeanValues();
  // int size_vector = vector_mean->GetNoElements();
  // for(int idx=0; idx<size_vector; idx++) {
  //   double val_mean = (*vector_mean)(idx);
  //   double val_nominal = h1_Var_wiSel->GetBinContent(idx+1) - h1_CV_wiSel->GetBinContent(idx+1);
  //   double val_diff = val_mean - val_nominal;
  //   cout<<TString::Format(" check mean: bin %2d, nominal %8.2f, toy %8.2f, diff %8.2f",
  // 			  idx+1, val_nominal, val_mean, val_diff
  // 			  )<<endl;    
  // }

  h1_diff_weight->Add( h1_weight, h1_sampling, 1, -1./ntoy );
    
  ////////////////////////
  
  map_matrix_cov_on_absdiff[num_DetVar].ResizeTo(rows, rows);
  map_matrix_cov_on_absdiff[num_DetVar] = matrix_cov_on_absdiff;

  roostr = TString::Format("map_h1_CV_noSel_%02d", num_DetVar); map_h1_CV_noSel[num_DetVar] = (TH1D*)h1_CV_noSel->Clone(roostr);
  roostr = TString::Format("map_h1_CV_wiSel_%02d", num_DetVar); map_h1_CV_wiSel[num_DetVar] = (TH1D*)h1_CV_wiSel->Clone(roostr);
  roostr = TString::Format("map_h1_Var_noSel_%02d", num_DetVar); map_h1_Var_noSel[num_DetVar] = (TH1D*)h1_Var_noSel->Clone(roostr);
  roostr = TString::Format("map_h1_Var_wiSel_%02d", num_DetVar); map_h1_Var_wiSel[num_DetVar] = (TH1D*)h1_Var_wiSel->Clone(roostr);
  roostr = TString::Format("map_h1_diff_Var2CV_%02d", num_DetVar); map_h1_diff_Var2CV[num_DetVar] = (TH1D*)h1_diff_Var2CV->Clone(roostr);
  roostr = TString::Format("map_h1_binomial_err_%02d", num_DetVar); map_h1_binomial_err[num_DetVar] = (TH1D*)h1_binomial_err->Clone(roostr);
  
  ////////////////////////
  
  delete[] array_test;
  delete[] array_CV;
  delete[] array_Var;
  
  delete[] array_mean_Var2CV;
  delete[] array_total_CV_wiSel;
  delete[] array_total_Var_wiSel;
  
}

void TDet::Clear()
{
  num_DetVar++;

  map_total_index_string.clear();
  map_total_index_CV_Erec.clear();
  map_total_index_Var_Erec.clear();  
  map_CV_wiSel_Erec.clear();
  map_Var_wiSel_Erec.clear();
  
  if( num_DetVar==1 ) {
    h1_CV_noSel  = NULL;
    h1_CV_wiSel  = NULL;
    h1_Var_noSel = NULL;
    h1_Var_wiSel = NULL;
    h1_diff_Var2CV = NULL;
    h1_binomial_err= NULL;
      
    h1_weight      = NULL;
    h1_sampling    = NULL;
    h1_diff_weight = NULL;
  }
  else {      
    delete h1_CV_noSel;
    delete h1_CV_wiSel;
    delete h1_Var_noSel;
    delete h1_Var_wiSel;
    delete h1_diff_Var2CV;
    delete h1_binomial_err;
    
    delete h1_weight;
    delete h1_sampling;
    delete h1_diff_weight;
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

  map<int, TString>map_bnb_numu_file_CV;
  map_bnb_numu_file_CV[1] = "./data_det_syst_total/common_detvar/LYDown/numu_CV";
  map_bnb_numu_file_CV[2] = "./data_det_syst_total/common_detvar/Recomb2/numu_CV";
  map_bnb_numu_file_CV[3] = "./data_det_syst_total/common_detvar/SCE/numu_CV";
  map_bnb_numu_file_CV[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_CV";
  map_bnb_numu_file_CV[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_CV";
  map_bnb_numu_file_CV[6] = "./data_det_syst_total/common_detvar/WireModX/numu_CV";
  map_bnb_numu_file_CV[7] = "./data_det_syst_total/common_detvar/WireModYZ/numu_CV";
  
  map<int, TString>map_bnb_numu_file_Var;
  map_bnb_numu_file_Var[1] = "./data_det_syst_total/common_detvar/LYDown/numu_Var";
  map_bnb_numu_file_Var[2] = "./data_det_syst_total/common_detvar/Recomb2/numu_Var";
  map_bnb_numu_file_Var[3] = "./data_det_syst_total/common_detvar/SCE/numu_Var";
  map_bnb_numu_file_Var[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_Var";
  map_bnb_numu_file_Var[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_Var";
  map_bnb_numu_file_Var[6] = "./data_det_syst_total/common_detvar/WireModX/numu_Var";
  map_bnb_numu_file_Var[7] = "./data_det_syst_total/common_detvar/WireModYZ/numu_Var";
  
  //
  TDet *det_test = new TDet();

  //det_test->Clear();  
  //det_test->Exe("./data_det_syst_total/common_detvar/LYDown/numu_CV", "./data_det_syst_total/common_detvar/LYDown/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/Recomb2/numu_CV", "./data_det_syst_total/common_detvar/Recomb2/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/SCE/numu_CV", "./data_det_syst_total/common_detvar/SCE/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModThetaXZ/numu_CV", "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModThetaYZ/numu_CV", "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModX/numu_CV", "./data_det_syst_total/common_detvar/WireModX/numu_Var", 1, 1, 1000);
  //det_test->Exe("./data_det_syst_total/common_detvar/WireModYZ/numu_CV", "./data_det_syst_total/common_detvar/WireModYZ/numu_Var", 1, 1, 1000);


  int nn_DetVar = 1;

  
  for(int idet=1; idet<=nn_DetVar; idet++) {

    cout<<endl<<" ---> processing DetSyst: "<<map_bnb_numu_file_CV[idet]<<endl<<endl;
    
    det_test->Clear(); 
    det_test->Exe( map_bnb_numu_file_CV[idet], map_bnb_numu_file_Var[idet], 1, 1, 500 );
      
    ////////////////////////////////////////////////////////////////////////////////// plotting
    /*
    roostr = "canv_h1_weight";
    TCanvas *canv_h1_weight = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_weight, 0.15, 0.2,0.1,0.15);
    det_test->h1_weight->Draw();
    det_test->h1_weight->SetLineColor(kBlack);
    //h1_weight->SetLineStyle(7);
    func_title_size(det_test->h1_weight, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(det_test->h1_weight, "Event index", "Weight");
    det_test->h1_weight->GetXaxis()->SetNdivisions(506);

    roostr = "canv_h1_sampling";
    TCanvas *canv_h1_sampling = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_sampling, 0.15, 0.2,0.1,0.15);
    det_test->h1_sampling->Draw();
    det_test->h1_sampling->SetLineColor(kBlack);
    //h1_sampling->SetLineStyle(7);
    func_title_size(det_test->h1_sampling, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(det_test->h1_sampling, "Event index", "Weight");
    det_test->h1_sampling->GetXaxis()->SetNdivisions(506);

    roostr = "canv_h1_diff_weight";
    TCanvas *canv_h1_diff_weight = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_diff_weight, 0.15, 0.2,0.1,0.15);
    det_test->h1_diff_weight->Draw();
    det_test->h1_diff_weight->SetLineColor(kBlack);
    //h1_diff_weight->SetLineStyle(7);
    func_title_size(det_test->h1_diff_weight, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(det_test->h1_diff_weight, "Event index", "Diff of Weight");
    det_test->h1_diff_weight->GetXaxis()->SetNdivisions(506);

    ///////////
  
    roostr = "canv_spectra";
    TCanvas *canv_spectra = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_spectra, 0.15, 0.2,0.1,0.15);
    canv_spectra->SetLogy();
  
    det_test->h1_CV_noSel->Draw("hist text75");
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

    det_test->h1_diff_Var2CV->Draw("same hist text0");
    det_test->h1_diff_Var2CV->SetLineColor(kBlue);
    det_test->h1_diff_Var2CV->SetMarkerColor(kBlue);

    det_test->h1_binomial_err->Draw("same hist text75");
    det_test->h1_binomial_err->SetLineColor(kOrange-3);
    det_test->h1_binomial_err->SetMarkerColor(kOrange-3);
    
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
    roostr = "h1_cov2Berr";
    TH1D *h1_cov2Berr = new TH1D(roostr, roostr, bins_temp, low_temp, hgh_temp);
    for(int i=0; i<bins_temp; i++) {
      double cov_root = sqrt( det_test->matrix_cov_on_absdiff(i,i) );
      double Berr = det_test->h1_binomial_err->GetBinContent(i+1);
      double reldiff = cov_root/Berr;
      h1_cov2Berr->SetBinContent( i+1, reldiff );
    }
    double max_h1_cov2absdiff = h1_cov2absdiff->GetMaximum();
    double max_h1_cov2Berr = h1_cov2Berr->GetMaximum();
    if( max_h1_cov2absdiff<max_h1_cov2Berr ) max_h1_cov2absdiff = max_h1_cov2Berr;
  
    roostr = "canv_h1_cov2absdiff";
    TCanvas *canv_h1_cov2absdiff = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_cov2absdiff, 0.15, 0.2,0.1,0.15);
    h1_cov2absdiff->Draw("hist");
    h1_cov2absdiff->SetMaximum(max_h1_cov2absdiff * 1.2);
    h1_cov2absdiff->SetLineColor(kBlue);
    h1_cov2absdiff->SetTitle("");
    func_title_size(h1_cov2absdiff, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_cov2absdiff, "E_{rec} (MeV) ", "Ratio");
    h1_cov2absdiff->GetXaxis()->SetNdivisions(506);
    h1_cov2absdiff->GetYaxis()->SetNdivisions(506);
  
    TF1 *f1_h1_cov2absdiff = new TF1("f1_h1_cov2absdiff", "1", 0, 1e6);
    f1_h1_cov2absdiff->Draw("same");
    f1_h1_cov2absdiff->SetLineColor(kBlack);
    f1_h1_cov2absdiff->SetLineStyle(9);
  
    h1_cov2Berr->Draw("hist same");
    h1_cov2Berr->SetLineColor(kOrange-3);
    
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
    det_test->h1_diff_Var2CV->Draw("same hist text0");
    det_test->h1_binomial_err->Draw("same hist text75");
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
    f1_h1_cov2absdiff->Draw("same");
    f1_h1_cov2absdiff->SetLineStyle(9);
    h1_cov2Berr->Draw("hist same");
    h1_cov2absdiff->Draw("same axis");

    roostr = TString::Format("canv_sum_%02d.png", idet);
    canv_sum->SaveAs(roostr);
    */

    // cout<<TString::Format( " ---> check, CV_wiSel(1) %10.2f, Abs.Diff(1) %10.2f, cov(0,0) %10.2f",
    // 			   det_test->map_h1_CV_wiSel[idet]->GetBinContent(1),
    // 			   det_test->map_h1_diff_Var2CV[idet]->GetBinContent(1),
    // 			   det_test->map_matrix_cov_on_absdiff[idet](0,0) )<<endl;
  }


  // map_matrix_cov_on_absdiff[num_DetVar] = matrix_cov_on_absdiff;
  // roostr = TString::Format("map_h1_CV_noSel_%02d", num_DetVar); map_h1_CV_noSel[num_DetVar] = (TH1D*)h1_CV_noSel->Clone(roostr);
  // roostr = TString::Format("map_h1_CV_wiSel_%02d", num_DetVar); map_h1_CV_wiSel[num_DetVar] = (TH1D*)h1_CV_wiSel->Clone(roostr);
  // roostr = TString::Format("map_h1_Var_noSel_%02d", num_DetVar); map_h1_Var_noSel[num_DetVar] = (TH1D*)h1_Var_noSel->Clone(roostr);
  // roostr = TString::Format("map_h1_Var_wiSel_%02d", num_DetVar); map_h1_Var_wiSel[num_DetVar] = (TH1D*)h1_Var_wiSel->Clone(roostr);
  // roostr = TString::Format("map_h1_diff_Var2CV_%02d", num_DetVar); map_h1_diff_Var2CV[num_DetVar] = (TH1D*)h1_diff_Var2CV->Clone(roostr);
  // roostr = TString::Format("map_h1_binomial_err_%02d", num_DetVar); map_h1_binomial_err[num_DetVar] = (TH1D*)h1_binomial_err->Clone(roostr);

  //////////////////////////////////////////////////// Generating covariance matrix

  int nn_bin_test = nn_DetVar * det_test->map_h1_CV_wiSel[1]->GetNbinsX();
  int line_test = 0;
  
  TPrincipal principal_test(nn_bin_test, "ND");  
  double *array_test = new double[nn_bin_test];

  map<int, double>map_val_CV;
  
  ///  
  int nn_Toy = 1000;

  for( int itoy=1; itoy<=nn_Toy; itoy++ ) {

    ///
    line_test = 0;
    for(int idx=1; idx<=nn_bin_test; idx++) array_test[idx-1] = 0;

    ///
    TRandom *random3 = new TRandom3(0);    
    double rel_random = random3->Gaus(0, 1);

    ///
    for(int idet=1; idet<=nn_DetVar; idet++) {

      int nn_bin = det_test->map_h1_CV_wiSel[idet]->GetNbinsX();



      
      TMatrixDSym DSmatrix_cov_DetVar(nn_bin);
      for(int ibin=0; ibin<nn_bin; ibin++) {
	for(int jbin=0; jbin<nn_bin; jbin++) {
	  DSmatrix_cov_DetVar(ibin, jbin) = det_test->map_matrix_cov_on_absdiff[idet](ibin, jbin);
	}
      }
      TMatrixDSymEigen DSmatrix_eigen_DetVar( DSmatrix_cov_DetVar );
      TMatrixD matrix_eigenvector_DetVar = DSmatrix_eigen_DetVar.GetEigenVectors();
      TMatrixD matrix_eigenvector_DetVar_T(nn_bin, nn_bin);
      matrix_eigenvector_DetVar_T.Transpose( matrix_eigenvector_DetVar );
      TMatrixD matrix_cov_DetVar_diag = matrix_eigenvector_DetVar_T * (det_test->map_matrix_cov_on_absdiff[idet]) * matrix_eigenvector_DetVar;
      for(int i=0; i<nn_bin; i++) {
	for(int j=0; j<nn_bin; j++) {
	  if( matrix_cov_DetVar_diag(i,j)<1e-8 ) matrix_cov_DetVar_diag(i,j) = 0;
	}
      }
      TMatrixD matrix_element(nn_bin, 1);    
      for(int j=0; j<nn_bin; j++) {
	matrix_element(j,0) = random3->Gaus( 0, sqrt( matrix_cov_DetVar_diag(j,j) ) );      
      }
      TMatrixD matrix_variation = matrix_eigenvector_DetVar * matrix_element;      
      map<int, double>user_array;
      for(int idx=0; idx<nn_bin; idx++) user_array[idx] = matrix_variation(idx, 0);

      
      
      
      for(int ibin=1; ibin<=nn_bin; ibin++) {
	double val_CV = det_test->map_h1_CV_wiSel[idet]->GetBinContent(ibin);
	double val_Var = det_test->map_h1_Var_wiSel[idet]->GetBinContent(ibin);
	double val_AbsDiff = val_Var - val_CV;
	val_AbsDiff += user_array[ibin];
	double val_variation = rel_random * val_AbsDiff;
	double val_CV_with_variation = val_CV + val_variation;
	
	///
	line_test++;
	array_test[line_test-1] = val_CV_with_variation;
	
	map_val_CV[line_test-1] = val_CV;
      }// ibin
      
    }// idet

    principal_test.AddRow( array_test );

    ///
    delete random3;
    
  }// itoy
  
  
  ///////////////////
  
  TMatrixD *tt_matrix_cov_on_absdiff = (TMatrixD *)principal_test.GetCovarianceMatrix();
  int rows = tt_matrix_cov_on_absdiff->GetNrows();

  TMatrixD matrix_cov_on_absdiff(2,2);
  matrix_cov_on_absdiff.ResizeTo(rows, rows);
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_on_absdiff)(i,j) = (*tt_matrix_cov_on_absdiff)(j,i);
      matrix_cov_on_absdiff(i,j) = (*tt_matrix_cov_on_absdiff)(i,j);
    }
  }
  
  ////////////////////
  TH2D *h2_matrix_cov_on_absdiff = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  TH2D *h2_matrix_correlation_on_absdiff = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      
      double val_cov = matrix_cov_on_absdiff(i,j);
      double sigma_i = sqrt( matrix_cov_on_absdiff(i,i) );
      double sigma_j = sqrt( matrix_cov_on_absdiff(j,j) );      
      
      h2_matrix_cov_on_absdiff->SetBinContent(i+1,j+1, val_cov/map_val_CV[i]/map_val_CV[j]);

      h2_matrix_correlation_on_absdiff->SetBinContent( i+1,j+1, val_cov/sigma_i/sigma_j );
    }
  }

  
  roostr = "canv_h2_matrix_correlation_on_absdiff";
  TCanvas *canv_h2_matrix_correlation_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_correlation_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_matrix_correlation_on_absdiff->Draw("colz");
  h2_matrix_correlation_on_absdiff->SetTitle("");
  func_title_size(h2_matrix_correlation_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_correlation_on_absdiff, "Index", "Index");
  h2_matrix_correlation_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_matrix_correlation_on_absdiff->GetYaxis()->SetNdivisions(506);

  roostr = "canv_h2_matrix_cov_on_absdiff";
  TCanvas *canv_h2_matrix_cov_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff->Draw("colz");
  h2_matrix_cov_on_absdiff->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff, "Index", "Index");
  h2_matrix_cov_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff->GetYaxis()->SetNdivisions(506);

  
}
