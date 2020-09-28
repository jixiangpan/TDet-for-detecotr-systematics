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
    num_nue_DetVar = 0;
    num_numu_DetVar = 0;
  }

  void Clear();

  void Exe(TString file_CV, TString file_Var, bool flag_numu, int flag_FC, int ntoy);// flag_FC, 0: (patially), 1: fully, 2:no cut

  void Exe_channels_from_bnbnu(TString file_CV, TString file_Var, int idet, int ntoy);
  void Exe_channels_from_intrinsic(TString file_CV, TString file_Var, int idet, int ntoy);
  
  void Exe_obj_cov(int idet, int ntoy);

  void Exe_whole_cov();
  
  //////////////////////////

  int ntoy;
  
  int num_DetVar;
  int num_nue_DetVar;
  int num_numu_DetVar;
  
  TMatrixD matrix_cov_on_absdiff;
  TMatrixD matrix_cov_on_absdiff_check;
  TMatrixD matrix_check_no_CovAbsDiff;
  TMatrixD matrix_cov_obj;
  
  TMatrixD matrix_subset_cov_obj_wiStat;
  TMatrixD matrix_subset_cov_obj_noStat;
  TH1D *h1_subset_CV_wiSel;
  TH1D *h1_subset_Var_wiSel;
  
  int nbins_hist_nue;
  double low_hist_nue;
  double hgh_hist_nue;
  
  int nbins_hist_numu;
  double low_hist_numu;
  double hgh_hist_numu;

  void Set_nue_hist(int nbins, double low, double hgh) {
    nbins_hist_nue = nbins;
    low_hist_nue = low;
    hgh_hist_nue = hgh;
  }
  
  void Set_numu_hist(int nbins, double low, double hgh) {
    nbins_hist_numu = nbins;
    low_hist_numu = low;
    hgh_hist_numu = hgh;
  }
  
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

  /////////////////////////////////////////////////////////// channels from bnbnu
  
  map<TString, double>map_CV_wiSel_Erec_bnbnu_numu_FC;
  map<TString, double>map_Var_wiSel_Erec_bnbnu_numu_FC;
  map<TString, double>map_CV_wiSel_Erec_bnbnu_nue_FC;
  map<TString, double>map_Var_wiSel_Erec_bnbnu_nue_FC;
  map<TString, double>map_CV_wiSel_Erec_bnbnu_numu_PC;
  map<TString, double>map_Var_wiSel_Erec_bnbnu_numu_PC;
  map<TString, double>map_CV_wiSel_Erec_bnbnu_nue_PC;
  map<TString, double>map_Var_wiSel_Erec_bnbnu_nue_PC;

  ///
  map<int, TH1D*>map_bnbnu_h1_CV_wiSel;
  map<int, TH1D*>map_bnbnu_h1_Var_wiSel;
  map<int, TH1D*>map_bnbnu_h1_CV_noSel;
  map<int, TH1D*>map_bnbnu_h1_Var_noSel;
  map<int, TMatrixD>map_bnbnu_matrix_cov_on_AbsDiff_wiStat;
  map<int, TMatrixD>map_bnbnu_matrix_cov_on_AbsDiff_noStat;
  map<int, TMatrixD>map_bnbnu_matrix_cov_obj_wiStat;
  map<int, TMatrixD>map_bnbnu_matrix_cov_obj_noStat;
  
  /////////////////////////////////////////////////////////// channels from intrinsic
  
  map<TString, double>map_CV_wiSel_Erec_intrinsic_LEE_FC;
  map<TString, double>map_Var_wiSel_Erec_intrinsic_LEE_FC;
  map<TString, double>map_CV_wiSel_Erec_intrinsic_nue_FC;
  map<TString, double>map_Var_wiSel_Erec_intrinsic_nue_FC;
  map<TString, double>map_CV_wiSel_Erec_intrinsic_LEE_PC;
  map<TString, double>map_Var_wiSel_Erec_intrinsic_LEE_PC;
  map<TString, double>map_CV_wiSel_Erec_intrinsic_nue_PC;
  map<TString, double>map_Var_wiSel_Erec_intrinsic_nue_PC;

  ///
  map<int, TH1D*>map_intrinsic_h1_CV_wiSel;
  map<int, TH1D*>map_intrinsic_h1_Var_wiSel;
  map<int, TH1D*>map_intrinsic_h1_CV_noSel;
  map<int, TH1D*>map_intrinsic_h1_Var_noSel;
  map<int, TMatrixD>map_intrinsic_matrix_cov_on_AbsDiff_wiStat;
  map<int, TMatrixD>map_intrinsic_matrix_cov_on_AbsDiff_noStat;
  map<int, TMatrixD>map_intrinsic_matrix_cov_obj_wiStat;
  map<int, TMatrixD>map_intrinsic_matrix_cov_obj_noStat;
  
  /////////////////////////////////////////////////////////// subset for all channels
  
  map<int, TMatrixD>map_matrix_subset_cov_obj_wiStat_fractional;
  map<int, TMatrixD>map_matrix_subset_cov_obj_noStat_fractional;

  TMatrixD matrix_whole_cov_obj_wiStat_fractional;
  TMatrixD matrix_whole_cov_obj_noStat_fractional;
  
};
 
/////////
/////////
/////////

void TDet::Exe_whole_cov()
{
  TString roostr = "";
  
  cout<<endl<<" ---> processing whole covariance matrix"<<endl<<endl;

  ////////
  int rows = 0;   
  map<int, TMatrixD>::iterator it_map;  
  for( it_map=map_matrix_subset_cov_obj_wiStat_fractional.begin();  it_map!=map_matrix_subset_cov_obj_wiStat_fractional.end(); it_map++) {
    int idet = it_map->first;
    rows = map_matrix_subset_cov_obj_wiStat_fractional[idet].GetNrows();
  }

  matrix_whole_cov_obj_wiStat_fractional.ResizeTo(rows, rows);
  matrix_whole_cov_obj_noStat_fractional.ResizeTo(rows, rows);  
  
  ////////
  for( it_map=map_matrix_subset_cov_obj_wiStat_fractional.begin();  it_map!=map_matrix_subset_cov_obj_wiStat_fractional.end(); it_map++) {
    int idet = it_map->first;
    cout<<TString::Format(" ---> processing Det.Syst %2d", idet)<<endl;

    matrix_whole_cov_obj_wiStat_fractional += map_matrix_subset_cov_obj_wiStat_fractional[idet];
    matrix_whole_cov_obj_noStat_fractional += map_matrix_subset_cov_obj_noStat_fractional[idet];
  }

  cout<<endl;
}

/////////
/////////
/////////

void TDet::Exe_obj_cov(int idet, int ntoy)
{
  TString roostr = "";
  
  int rows_bnbnu = map_bnbnu_matrix_cov_on_AbsDiff_wiStat[idet].GetNrows();
  int rows_intrinsic = map_intrinsic_matrix_cov_on_AbsDiff_wiStat[idet].GetNrows();
  //rows_intrinsic = 0;
    
  ///////////////////////////////////

  roostr = "h1_CV_noSel"; h1_CV_noSel = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_CV_wiSel"; h1_CV_wiSel = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_Var_noSel"; h1_Var_noSel = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_Var_wiSel"; h1_Var_wiSel = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_diff_Var2CV"; h1_diff_Var2CV = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_binomial_err"; h1_binomial_err = new TH1D(roostr, roostr, 1, 0, 1);

  roostr = "h1_weight"; h1_weight = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_sampling"; h1_sampling = new TH1D(roostr, roostr, 1, 0, 1);
  roostr = "h1_diff_weight"; h1_diff_weight = new TH1D(roostr, roostr, 1, 0, 1);
   
  ///////////////////////////////////
  
  int rows = rows_bnbnu + rows_intrinsic;

  matrix_cov_on_absdiff.ResizeTo(rows, rows);
  int eff_line = 0;

  map<int, double>map_val_CV;
  map<int, double>map_val_Var;

  roostr = "h1_subset_CV_wiSel";
  h1_subset_CV_wiSel = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
  
  roostr = "h1_subset_Var_wiSel";
  h1_subset_Var_wiSel = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
    
  //////
  eff_line = 0;
  
  for(int i=0; i<rows_bnbnu; i++) {    
    for(int j=0; j<rows_bnbnu; j++) {
      matrix_cov_on_absdiff(i+eff_line, j+eff_line) = map_bnbnu_matrix_cov_on_AbsDiff_wiStat[idet](i,j);
    }
  }
  
  for(int ibin=1; ibin<=rows_bnbnu; ibin++) {
    double content = map_bnbnu_h1_CV_wiSel[idet]->GetBinContent( ibin );
    map_val_CV[ ibin-1 +eff_line] = content;
    h1_subset_CV_wiSel->SetBinContent( ibin+eff_line, content );
  }
  
  for(int ibin=1; ibin<=rows_bnbnu; ibin++) {
    double content = map_bnbnu_h1_Var_wiSel[idet]->GetBinContent( ibin );
    map_val_Var[ ibin-1 +eff_line] = content;
    h1_subset_Var_wiSel->SetBinContent( ibin+eff_line, content );
  }
  
  //////
  eff_line += rows_bnbnu;
  
  for(int i=0; i<rows_intrinsic; i++) {    
    for(int j=0; j<rows_intrinsic; j++) {
      matrix_cov_on_absdiff(i+eff_line, j+eff_line) = map_intrinsic_matrix_cov_on_AbsDiff_wiStat[idet](i,j);
    }
  }
  
  for(int ibin=1; ibin<=rows_intrinsic; ibin++) {
    double content = map_intrinsic_h1_CV_wiSel[idet]->GetBinContent( ibin );
    map_val_CV[ ibin-1 +eff_line] = content;
    h1_subset_CV_wiSel->SetBinContent( ibin+eff_line, content );
  }
  
  for(int ibin=1; ibin<=rows_intrinsic; ibin++) {
    double content = map_intrinsic_h1_Var_wiSel[idet]->GetBinContent( ibin );
    map_val_Var[ ibin-1 +eff_line] = content;
    h1_subset_Var_wiSel->SetBinContent( ibin+eff_line, content );
  }
  
  //////
  TMatrixDSym DSmatrix_cov_AbsDiff(rows);
  for(int ibin=0; ibin<rows; ibin++) {
    for(int jbin=0; jbin<rows; jbin++) {
      DSmatrix_cov_AbsDiff(ibin, jbin) = matrix_cov_on_absdiff(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov_AbsDiff );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TMatrixD matrix_eigenvector_T(rows, rows);
  matrix_eigenvector_T.Transpose( matrix_eigenvector );
  TMatrixD matrix_cov_AbsDiff_diag = matrix_eigenvector_T * matrix_cov_on_absdiff * matrix_eigenvector;
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      //if( i!=j )
      if( fabs(matrix_cov_AbsDiff_diag(i,j))<1e-10 ) matrix_cov_AbsDiff_diag(i,j) = 0;
    }
  }

  //////////////////////////////////////////////////
  
  TPrincipal principal_obj(rows, "ND");  
  TPrincipal principal_obj_noStat(rows, "ND");  
    
  for( int itoy=1; itoy<=ntoy; itoy++ ) {
    
    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
    TRandom *random3 = new TRandom3(0);
    double *array_obj = new double[rows];
    double *array_obj_noStat = new double[rows];
      
    TMatrixD matrix_element(rows, 1);    
    for(int j=0; j<rows; j++) {
      matrix_element(j,0) = random3->Gaus( 0, sqrt( matrix_cov_AbsDiff_diag(j,j) ) );      
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;      
    double *array_variation_AbsDiff = new double[rows];    
    
    ///
    for(int idx=0; idx<rows; idx++) {
      array_variation_AbsDiff[idx] = matrix_variation(idx, 0);      
    }
    
    ///
    double rel_random = random3->Gaus(0, 1);

    ///
    for(int ibin=1; ibin<=rows; ibin++) {
      double val_CV = map_val_CV[ibin-1];
      double val_Var = map_val_Var[ibin-1];
      double val_AbsDiff = val_Var - val_CV;

      ////
      array_obj_noStat[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;
      
      ////
      val_AbsDiff += array_variation_AbsDiff[ ibin-1 ];
      array_obj[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;      
    }

    //////
    principal_obj.AddRow( array_obj );
    principal_obj_noStat.AddRow( array_obj_noStat );
    
    //////
    delete random3;
    delete[] array_obj;
    delete[] array_obj_noStat;
    delete[] array_variation_AbsDiff;
      
  }// toy
  
  //////////////////////////////////////////////////////////////////////////////
  
  TMatrixD *tt_matrix_subset_cov_obj_wiStat = (TMatrixD *)principal_obj.GetCovarianceMatrix();
  matrix_subset_cov_obj_wiStat.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_subset_cov_obj_wiStat)(i,j) = (*tt_matrix_subset_cov_obj_wiStat)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_subset_cov_obj_wiStat)(i,j))<1e-10 ) (*tt_matrix_subset_cov_obj_wiStat)(i,j) = 0;
      matrix_subset_cov_obj_wiStat(i,j) = (*tt_matrix_subset_cov_obj_wiStat)(i,j);
    }
  }
  for(int i=0; i<rows; i++) {// for content(bin) = 0, pay attention
    if( h1_subset_CV_wiSel->GetBinContent( i+1 )==0 && h1_subset_Var_wiSel->GetBinContent( i+1 )==0 ) {
      matrix_subset_cov_obj_wiStat(i,i) = 2;
    }
  }
    
  TMatrixD *tt_matrix_subset_cov_obj_noStat = (TMatrixD *)principal_obj_noStat.GetCovarianceMatrix();
  matrix_subset_cov_obj_noStat.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_subset_cov_obj_noStat)(i,j) = (*tt_matrix_subset_cov_obj_noStat)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_subset_cov_obj_noStat)(i,j))<1e-10 ) (*tt_matrix_subset_cov_obj_noStat)(i,j) = 0;
      matrix_subset_cov_obj_noStat(i,j) = (*tt_matrix_subset_cov_obj_noStat)(i,j);
    }
  }
  for(int i=0; i<rows; i++) {// for content(bin) = 0, pay attention
    if( h1_subset_CV_wiSel->GetBinContent( i+1 )==0 && h1_subset_Var_wiSel->GetBinContent( i+1 )==0 ) {
      matrix_subset_cov_obj_noStat(i,i) = 2;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////// map_matrix_subset_cov_obj_wiStat_fractional ttttt

  map_matrix_subset_cov_obj_wiStat_fractional[idet].ResizeTo( rows, rows );
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_subset_cov_obj_wiStat(i,j);
      double val_i = h1_subset_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_subset_CV_wiSel->GetBinContent(j+1);

      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_subset_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_subset_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      map_matrix_subset_cov_obj_wiStat_fractional[idet](i,j) = content;

      //if(i==j) cout<<TString::Format(" ---> check wiStat %3d, %10.6f", i+1, sqrt(content))<<endl;
    }
  }

  //////////
  map_matrix_subset_cov_obj_noStat_fractional[idet].ResizeTo( rows, rows );
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_subset_cov_obj_noStat(i,j);
      double val_i = h1_subset_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_subset_CV_wiSel->GetBinContent(j+1);

      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_subset_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_subset_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      map_matrix_subset_cov_obj_noStat_fractional[idet](i,j) = content;

      //if(i==j) cout<<TString::Format(" ---> check noStat %3d, %10.6f", i+1, sqrt(content))<<endl;
    }
  }
  
  ///////////

  // for(int i=0; i<rows; i++) {
  //   cout<<TString::Format(" ---> check noStat %3d,   %6.4f   %6.4f", i+1,
  // 			  sqrt( map_matrix_subset_cov_obj_wiStat_fractional[idet](i,i) ),
  // 			  sqrt( map_matrix_subset_cov_obj_noStat_fractional[idet](i,i) )
  // 			  )<<endl;
  // }
  
}

/////////
/////////
/////////

void TDet::Exe_channels_from_bnbnu(TString file_CV, TString file_Var, int idet, int ntoy)
{
  // channels:
  // numuCC contained, uncontained
  // nueCC contained, uncontained
  // ---> FlagIndex
  
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

  int user_nbins_hist_nue  = nbins_hist_nue;
  double user_low_hist_nue = low_hist_nue;
  double user_hgh_hist_nue = hgh_hist_nue;
  
  int user_nbins_hist_numu  = nbins_hist_numu;
  double user_low_hist_numu = low_hist_numu;
  double user_hgh_hist_numu = hgh_hist_numu;
  
  roostr = "h1_basic_nue";
  TH1D *h1_basic_nue = new TH1D(roostr, roostr, user_nbins_hist_nue, user_low_hist_nue, user_hgh_hist_nue);  
  roostr = "h1_basic_numu";
  TH1D *h1_basic_numu = new TH1D(roostr, roostr, user_nbins_hist_numu, user_low_hist_numu, user_hgh_hist_numu);


  int bins_basic = user_nbins_hist_numu + user_nbins_hist_numu + user_nbins_hist_nue + user_nbins_hist_nue;// channels
  double low_basic = 0.5;
  double hgh_basic = bins_basic+0.5;

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

    int index_basic_numu = h1_basic_numu->FindBin(Erec);
    int index_basic_nue = h1_basic_nue->FindBin(Erec);

    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {
      int eff_index = index_basic_numu + 0;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
    }
    
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {
      int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
    }
    
    /////////////// numu       
    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {      
      if( status==1 && post_generic==1 && numuCC==1 ) {
        if( is_FC==1 ) {
          map_CV_wiSel_Erec_bnbnu_numu_FC[roostr] = Erec;
          int eff_index = index_basic_numu + 0;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_CV_wiSel_Erec_bnbnu_numu_PC[roostr] = Erec;
          int eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
 
    /////////////// nue       
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_CV_wiSel_Erec_bnbnu_nue_FC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_CV_wiSel_Erec_bnbnu_nue_PC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
           
  }// CV
   
  //////////////////////// Var
  
  ifstream InputFile_Var_a1(file_Var, ios::in);
  for(int idx=1; idx<=line_Var_a0; idx++) {
    
    InputFile_Var_a1 >> entry >> run >> subrun >> event
                    >> weight >> status >> post_generic >> numuCC >> nueBDT >> is_FC >> Erec >> Evis;
    
    Erec *= 1000;
    roostr = TString::Format("%d_%d_%d", run, subrun, event);

    int index_basic_numu = h1_basic_numu->FindBin(Erec);
    int index_basic_nue = h1_basic_nue->FindBin(Erec);

    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {
      int eff_index = index_basic_numu + 0;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
    }
    
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {
      int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
    }
    
    /////////////// numu       
    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {      
      if( status==1 && post_generic==1 && numuCC==1 ) {
        if( is_FC==1 ) {
          map_Var_wiSel_Erec_bnbnu_numu_FC[roostr] = Erec;
          int eff_index = index_basic_numu + 0;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_Var_wiSel_Erec_bnbnu_numu_PC[roostr] = Erec;
          int eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
 
    /////////////// nue       
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_Var_wiSel_Erec_bnbnu_nue_FC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_Var_wiSel_Erec_bnbnu_nue_PC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );   
        }
      } 
    }
           
  }// Var
   
  ///////////////////////////////////////////////

  for(int ibin=1; ibin<=bins_basic; ibin++) {

    double content_noSel_CV = 0;
    content_noSel_CV = h1_CV_noSel->GetBinContent(ibin);
    //content_noSel_CV = ( (int)(content_noSel_CV*100+0.5) )*1./100;
    h1_CV_noSel->SetBinContent(ibin, content_noSel_CV);
    
    double content_noSel_Var = 0;
    content_noSel_Var = h1_Var_noSel->GetBinContent(ibin);
    //content_noSel_Var = ( (int)(content_noSel_Var*100+0.5) )*1./100;
    h1_Var_noSel->SetBinContent(ibin, content_noSel_Var);
    
    double content_wiSel_CV = 0;
    content_wiSel_CV = h1_CV_wiSel->GetBinContent(ibin);
    //content_wiSel_CV = ( (int)(content_wiSel_CV*100+0.5) )*1./100;
    h1_CV_wiSel->SetBinContent(ibin, content_wiSel_CV);
    
    double content_wiSel_Var = 0;
    content_wiSel_Var = h1_Var_wiSel->GetBinContent(ibin);
    //content_wiSel_Var = ( (int)(content_wiSel_Var*100+0.5) )*1./100;
    h1_Var_wiSel->SetBinContent(ibin, content_wiSel_Var);
    
    h1_diff_Var2CV->SetBinContent( ibin, fabs(content_wiSel_Var-content_wiSel_CV) );

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
    
  /////////////////////////////////////////////// Calculate covariance matrix of abs.diff

  TPrincipal principal_test(bins_basic, "ND");
  double *array_test = new double[bins_basic];
  double *array_CV = new double[bins_basic];
  double *array_Var = new double[bins_basic];

  /// to check
  roostr = "h1_sampling";
  h1_sampling = new TH1D(roostr, roostr, line_CV_a0, 0.5, line_CV_a0+0.5);

  double *array_mean_Var2CV = new double[bins_basic];  
  for(int ibin=0; ibin<bins_basic; ibin++) {
    array_mean_Var2CV[ibin] = 0;
  }
  
  ///
  for(int itoy=1; itoy<=ntoy; itoy++) {

    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy ( Cov of Abs.Diff ): %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
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

      /// CV
      double CV_Erec = map_total_index_CV_Erec[roostr];            
      int index_basic_numu_CV = h1_basic_numu->FindBin(CV_Erec);
      int index_basic_nue_CV = h1_basic_nue->FindBin(CV_Erec);
      int eff_index_numu_CV_FC = index_basic_numu_CV + 0;// ---> FlagIndex
      int eff_index_numu_CV_PC = index_basic_numu_CV + user_nbins_hist_numu;// ---> FlagIndex
      int eff_index_nue_CV_FC = index_basic_nue_CV + user_nbins_hist_numu*2;// ---> FlagIndex
      int eff_index_nue_CV_PC = index_basic_nue_CV + user_nbins_hist_nue + user_nbins_hist_numu*2;   // ---> FlagIndex  
      if( index_basic_numu_CV>=1 && index_basic_numu_CV<=user_nbins_hist_numu ) {
        if( map_CV_wiSel_Erec_bnbnu_numu_FC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_numu_FC.end() ) {
          array_CV[eff_index_numu_CV_FC -1] += 1;
        }
        if( map_CV_wiSel_Erec_bnbnu_numu_PC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_numu_PC.end() ) {
          array_CV[eff_index_numu_CV_PC -1] += 1;
        }       
      }
      if( index_basic_nue_CV>=1 && index_basic_nue_CV<=user_nbins_hist_nue ) {
        if( map_CV_wiSel_Erec_bnbnu_nue_FC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_nue_FC.end() ) {
          array_CV[eff_index_nue_CV_FC -1] += 1;
        }
        if( map_CV_wiSel_Erec_bnbnu_nue_PC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_nue_PC.end() ) {
          array_CV[eff_index_nue_CV_PC -1] += 1;
        }       
      }

      /// Var
      double Var_Erec = map_total_index_Var_Erec[roostr];            
      int index_basic_numu_Var = h1_basic_numu->FindBin(Var_Erec);
      int index_basic_nue_Var = h1_basic_nue->FindBin(Var_Erec);
      int eff_index_numu_Var_FC = index_basic_numu_Var + 0;// ---> FlagIndex
      int eff_index_numu_Var_PC = index_basic_numu_Var + user_nbins_hist_numu;// ---> FlagIndex
      int eff_index_nue_Var_FC = index_basic_nue_Var + user_nbins_hist_numu*2;// ---> FlagIndex
      int eff_index_nue_Var_PC = index_basic_nue_Var + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex      
      if( index_basic_numu_Var>=1 && index_basic_numu_Var<=user_nbins_hist_numu ) {
        if( map_Var_wiSel_Erec_bnbnu_numu_FC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_numu_FC.end() ) {
          array_Var[eff_index_numu_Var_FC -1] += 1;
        }
        if( map_Var_wiSel_Erec_bnbnu_numu_PC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_numu_PC.end() ) {
          array_Var[eff_index_numu_Var_PC -1] += 1;
        }       
      }
      if( index_basic_nue_Var>=1 && index_basic_nue_Var<=user_nbins_hist_nue ) {
        if( map_Var_wiSel_Erec_bnbnu_nue_FC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_nue_FC.end() ) {
          array_Var[eff_index_nue_Var_FC -1] += 1;        
        }
        if( map_Var_wiSel_Erec_bnbnu_nue_PC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_nue_PC.end() ) {
          array_Var[eff_index_nue_Var_PC -1] += 1;
        }       
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
      //if(i!=j)
      if( fabs((*tt_matrix_cov_on_absdiff)(i,j))<1e-10 ) (*tt_matrix_cov_on_absdiff)(i,j) = 0;
      matrix_cov_on_absdiff(i,j) = (*tt_matrix_cov_on_absdiff)(i,j);
    }
  }
  for(int i=0; i<rows; i++) {// for content(bin) = 0, pay attention
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) {
      matrix_cov_on_absdiff(i,i) = 2;
    }
  }
  

  h1_diff_weight->Add( h1_weight, h1_sampling, 1, -1./ntoy );
    
  ///////////////////////////////////////////////

  delete h1_basic_nue;
  delete h1_basic_numu;
  
  delete[] array_test;
  delete[] array_CV;
  delete[] array_Var;  
  delete[] array_mean_Var2CV;

  /////////////////////////////////////////////// To handle N selected = 0, ...

  roostr = TString::Format("map_bnbnu_h1_CV_wiSel_%02d", idet);
  map_bnbnu_h1_CV_wiSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_CV_wiSel->GetBinContent( ibin );
    map_bnbnu_h1_CV_wiSel[idet]->SetBinContent( ibin, content );
  }
  
  roostr = TString::Format("map_bnbnu_h1_Var_wiSel_%02d", idet);
  map_bnbnu_h1_Var_wiSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_Var_wiSel->GetBinContent( ibin );
    map_bnbnu_h1_Var_wiSel[idet]->SetBinContent( ibin, content );
  }

  roostr = TString::Format("map_bnbnu_h1_CV_noSel_%02d", idet);
  map_bnbnu_h1_CV_noSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_CV_noSel->GetBinContent( ibin );
    map_bnbnu_h1_CV_noSel[idet]->SetBinContent( ibin, content );
  }
  
  roostr = TString::Format("map_bnbnu_h1_Var_noSel_%02d", idet);
  map_bnbnu_h1_Var_noSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_Var_noSel->GetBinContent( ibin );
    map_bnbnu_h1_Var_noSel[idet]->SetBinContent( ibin, content );
  }

  map_bnbnu_matrix_cov_on_AbsDiff_wiStat[idet].ResizeTo(rows, rows);
  map_bnbnu_matrix_cov_on_AbsDiff_wiStat[idet] = matrix_cov_on_absdiff;

  /////////////////////////////////////////////// generate covarinace matrix: from Abs.Diff + Error of Abs.Diff
  
  TPrincipal principal_obj(rows, "ND");  
  double *array_obj = new double[rows];

  TPrincipal principal_check(rows, "ND");
  
  TPrincipal principal_check_no_CovAbsDiff(rows, "ND"); 
  
  //////
  TMatrixDSym DSmatrix_cov_AbsDiff(rows);
  for(int ibin=0; ibin<rows; ibin++) {
    for(int jbin=0; jbin<rows; jbin++) {
      DSmatrix_cov_AbsDiff(ibin, jbin) = matrix_cov_on_absdiff(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov_AbsDiff );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TMatrixD matrix_eigenvector_T(rows, rows);
  matrix_eigenvector_T.Transpose( matrix_eigenvector );
  TMatrixD matrix_cov_AbsDiff_diag = matrix_eigenvector_T * matrix_cov_on_absdiff * matrix_eigenvector;
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      //if(i!=j)
      if( fabs(matrix_cov_AbsDiff_diag(i,j))<1e-10 ) matrix_cov_AbsDiff_diag(i,j) = 0;      
    }
  }

  for( int itoy=1; itoy<=ntoy; itoy++ ) {
    
    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy ( Covariance ): %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
    TRandom *random3 = new TRandom3(0);
    
    TMatrixD matrix_element(rows, 1);    
    for(int j=0; j<rows; j++) {
      matrix_element(j,0) = random3->Gaus( 0, sqrt( matrix_cov_AbsDiff_diag(j,j) ) );      
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;      
    double *array_variation_AbsDiff = new double[rows];    
    
    ///
    for(int idx=0; idx<rows; idx++) {
      array_variation_AbsDiff[idx] = matrix_variation(idx, 0);      
    }
    principal_check.AddRow( array_variation_AbsDiff );

    ///
    double rel_random = random3->Gaus(0, 1);

    ///
    double *array_check_no_CovAbsDiff = new double[rows];

    ///
    for(int ibin=1; ibin<=rows; ibin++) {
      double val_CV = h1_CV_wiSel->GetBinContent( ibin );
      double val_Var = h1_Var_wiSel->GetBinContent( ibin );
      double val_AbsDiff = val_Var - val_CV;

      ////
      array_check_no_CovAbsDiff[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;

      ////
      val_AbsDiff += array_variation_AbsDiff[ ibin-1 ];
      array_obj[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;
      
    }

    //////
    principal_check_no_CovAbsDiff.AddRow( array_check_no_CovAbsDiff );
    principal_obj.AddRow( array_obj );
    
    //////
    delete random3;
    delete[] array_variation_AbsDiff;
    delete[] array_check_no_CovAbsDiff;    
  }// toy

  ////////
  TMatrixD *tt_matrix_cov_on_absdiff_check = (TMatrixD *)principal_check.GetCovarianceMatrix();
  matrix_cov_on_absdiff_check.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_on_absdiff_check)(i,j) = (*tt_matrix_cov_on_absdiff_check)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_cov_on_absdiff_check)(i,j))<1e-10 ) (*tt_matrix_cov_on_absdiff_check)(i,j) = 0;
      matrix_cov_on_absdiff_check(i,j) = (*tt_matrix_cov_on_absdiff_check)(i,j);
    }
  }
  
  ////////
  TMatrixD *tt_matrix_check_no_CovAbsDiff = (TMatrixD *)principal_check_no_CovAbsDiff.GetCovarianceMatrix();
  matrix_check_no_CovAbsDiff.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_check_no_CovAbsDiff)(i,j) = (*tt_matrix_check_no_CovAbsDiff)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_check_no_CovAbsDiff)(i,j))<1e-10 ) (*tt_matrix_check_no_CovAbsDiff)(i,j) = 0;
      matrix_check_no_CovAbsDiff(i,j) = (*tt_matrix_check_no_CovAbsDiff)(i,j);
    }
  }
  
  for(int i=0; i<rows; i++) {    
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) matrix_check_no_CovAbsDiff(i,i) = 2;
  }
  
  ////////
  TMatrixD *tt_matrix_cov_obj = (TMatrixD *)principal_obj.GetCovarianceMatrix();
  matrix_cov_obj.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_obj)(i,j) = (*tt_matrix_cov_obj)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_cov_obj)(i,j))<1e-10 ) (*tt_matrix_cov_obj)(i,j) = 0;
      matrix_cov_obj(i,j) = (*tt_matrix_cov_obj)(i,j);
    }
  }

  for(int i=0; i<rows; i++) {    
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) matrix_cov_obj(i,i) = 2;
  }
  
  ///////////////////////////////////////////////
  
  map_bnbnu_matrix_cov_obj_wiStat[idet].ResizeTo(rows, rows);
  map_bnbnu_matrix_cov_obj_wiStat[idet] = matrix_cov_obj;
  
  map_bnbnu_matrix_cov_obj_noStat[idet].ResizeTo(rows, rows);
  map_bnbnu_matrix_cov_obj_noStat[idet] = matrix_check_no_CovAbsDiff;
  
  ///////////////////////////////////////////////
  
  delete[] array_obj;
}

/////////
/////////
/////////

void TDet::Exe_channels_from_intrinsic(TString file_CV, TString file_Var, int idet, int ntoy)
{
  // channels:
  // numuCC contained, uncontained
  // nueCC contained, uncontained
  // ---> FlagIndex
  
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

  int user_nbins_hist_nue  = nbins_hist_nue;
  double user_low_hist_nue = low_hist_nue;
  double user_hgh_hist_nue = hgh_hist_nue;
  
  int user_nbins_hist_numu  = nbins_hist_nue;/// LEE
  double user_low_hist_numu = low_hist_nue;
  double user_hgh_hist_numu = hgh_hist_nue;
  
  roostr = "h1_basic_nue";
  TH1D *h1_basic_nue = new TH1D(roostr, roostr, user_nbins_hist_nue, user_low_hist_nue, user_hgh_hist_nue);  
  roostr = "h1_basic_numu";
  TH1D *h1_basic_numu = new TH1D(roostr, roostr, user_nbins_hist_numu, user_low_hist_numu, user_hgh_hist_numu);


  int bins_basic = user_nbins_hist_numu + user_nbins_hist_numu + user_nbins_hist_nue + user_nbins_hist_nue;// channels
  double low_basic = 0.5;
  double hgh_basic = bins_basic+0.5;

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

    int index_basic_numu = h1_basic_numu->FindBin(Erec);
    int index_basic_nue = h1_basic_nue->FindBin(Erec);

    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {
      int eff_index = index_basic_numu + 0;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
    }
    
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {
      int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_CV_noSel->Fill( eff_index-0.5, weight );
    }
    
    /////////////// numu       
    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_CV_wiSel_Erec_bnbnu_numu_FC[roostr] = Erec;
          int eff_index = index_basic_numu + 0;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_CV_wiSel_Erec_bnbnu_numu_PC[roostr] = Erec;
          int eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
 
    /////////////// nue       
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_CV_wiSel_Erec_bnbnu_nue_FC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_CV_wiSel_Erec_bnbnu_nue_PC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_CV_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
           
  }// CV
   
  //////////////////////// Var
  
  ifstream InputFile_Var_a1(file_Var, ios::in);
  for(int idx=1; idx<=line_Var_a0; idx++) {
    
    InputFile_Var_a1 >> entry >> run >> subrun >> event
                    >> weight >> status >> post_generic >> numuCC >> nueBDT >> is_FC >> Erec >> Evis;
    
    Erec *= 1000;
    roostr = TString::Format("%d_%d_%d", run, subrun, event);

    int index_basic_numu = h1_basic_numu->FindBin(Erec);
    int index_basic_nue = h1_basic_nue->FindBin(Erec);

    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {
      int eff_index = index_basic_numu + 0;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
    }
    
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {
      int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
      eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
      h1_Var_noSel->Fill( eff_index-0.5, weight );
    }
    
    /////////////// numu       
    if( index_basic_numu>=1 && index_basic_numu<=user_nbins_hist_numu ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_Var_wiSel_Erec_bnbnu_numu_FC[roostr] = Erec;
          int eff_index = index_basic_numu + 0;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_Var_wiSel_Erec_bnbnu_numu_PC[roostr] = Erec;
          int eff_index = index_basic_numu + user_nbins_hist_numu;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
      } 
    }
 
    /////////////// nue       
    if( index_basic_nue>=1 && index_basic_nue<=user_nbins_hist_nue ) {      
      if( status==1 && post_generic==1 && nueBDT==1 ) {
        if( is_FC==1 ) {
          map_Var_wiSel_Erec_bnbnu_nue_FC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );
        }
        else {// PC
          map_Var_wiSel_Erec_bnbnu_nue_PC[roostr] = Erec;
          int eff_index = index_basic_nue + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex
          h1_Var_wiSel->Fill( eff_index-0.5, weight );   
        }
      } 
    }
           
  }// Var
   
  ///////////////////////////////////////////////

  for(int ibin=1; ibin<=bins_basic; ibin++) {

    double content_noSel_CV = 0;
    content_noSel_CV = h1_CV_noSel->GetBinContent(ibin);
    //content_noSel_CV = ( (int)(content_noSel_CV*100+0.5) )*1./100;
    h1_CV_noSel->SetBinContent(ibin, content_noSel_CV);
    
    double content_noSel_Var = 0;
    content_noSel_Var = h1_Var_noSel->GetBinContent(ibin);
    //content_noSel_Var = ( (int)(content_noSel_Var*100+0.5) )*1./100;
    h1_Var_noSel->SetBinContent(ibin, content_noSel_Var);
    
    double content_wiSel_CV = 0;
    content_wiSel_CV = h1_CV_wiSel->GetBinContent(ibin);
    //content_wiSel_CV = ( (int)(content_wiSel_CV*100+0.5) )*1./100;
    h1_CV_wiSel->SetBinContent(ibin, content_wiSel_CV);
    
    double content_wiSel_Var = 0;
    content_wiSel_Var = h1_Var_wiSel->GetBinContent(ibin);
    //content_wiSel_Var = ( (int)(content_wiSel_Var*100+0.5) )*1./100;
    h1_Var_wiSel->SetBinContent(ibin, content_wiSel_Var);
    
    h1_diff_Var2CV->SetBinContent( ibin, fabs(content_wiSel_Var-content_wiSel_CV) );

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
    
  /////////////////////////////////////////////// Calculate covariance matrix of abs.diff

  TPrincipal principal_test(bins_basic, "ND");
  double *array_test = new double[bins_basic];
  double *array_CV = new double[bins_basic];
  double *array_Var = new double[bins_basic];

  /// to check
  roostr = "h1_sampling";
  h1_sampling = new TH1D(roostr, roostr, line_CV_a0, 0.5, line_CV_a0+0.5);

  double *array_mean_Var2CV = new double[bins_basic];  
  for(int ibin=0; ibin<bins_basic; ibin++) {
    array_mean_Var2CV[ibin] = 0;
  }
  
  ///
  for(int itoy=1; itoy<=ntoy; itoy++) {

    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy ( Cov of Abs.Diff ): %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
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

      /// CV
      double CV_Erec = map_total_index_CV_Erec[roostr];            
      int index_basic_numu_CV = h1_basic_numu->FindBin(CV_Erec);
      int index_basic_nue_CV = h1_basic_nue->FindBin(CV_Erec);
      int eff_index_numu_CV_FC = index_basic_numu_CV + 0;// ---> FlagIndex
      int eff_index_numu_CV_PC = index_basic_numu_CV + user_nbins_hist_numu;// ---> FlagIndex
      int eff_index_nue_CV_FC = index_basic_nue_CV + user_nbins_hist_numu*2;// ---> FlagIndex
      int eff_index_nue_CV_PC = index_basic_nue_CV + user_nbins_hist_nue + user_nbins_hist_numu*2;   // ---> FlagIndex  
      if( index_basic_numu_CV>=1 && index_basic_numu_CV<=user_nbins_hist_numu ) {
        if( map_CV_wiSel_Erec_bnbnu_numu_FC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_numu_FC.end() ) {
          array_CV[eff_index_numu_CV_FC -1] += 1;
        }
        if( map_CV_wiSel_Erec_bnbnu_numu_PC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_numu_PC.end() ) {
          array_CV[eff_index_numu_CV_PC -1] += 1;
        }       
      }
      if( index_basic_nue_CV>=1 && index_basic_nue_CV<=user_nbins_hist_nue ) {
        if( map_CV_wiSel_Erec_bnbnu_nue_FC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_nue_FC.end() ) {
          array_CV[eff_index_nue_CV_FC -1] += 1;
        }
        if( map_CV_wiSel_Erec_bnbnu_nue_PC.find(roostr)!=map_CV_wiSel_Erec_bnbnu_nue_PC.end() ) {
          array_CV[eff_index_nue_CV_PC -1] += 1;
        }       
      }

      /// Var
      double Var_Erec = map_total_index_Var_Erec[roostr];            
      int index_basic_numu_Var = h1_basic_numu->FindBin(Var_Erec);
      int index_basic_nue_Var = h1_basic_nue->FindBin(Var_Erec);
      int eff_index_numu_Var_FC = index_basic_numu_Var + 0;// ---> FlagIndex
      int eff_index_numu_Var_PC = index_basic_numu_Var + user_nbins_hist_numu;// ---> FlagIndex
      int eff_index_nue_Var_FC = index_basic_nue_Var + user_nbins_hist_numu*2;// ---> FlagIndex
      int eff_index_nue_Var_PC = index_basic_nue_Var + user_nbins_hist_nue + user_nbins_hist_numu*2;// ---> FlagIndex      
      if( index_basic_numu_Var>=1 && index_basic_numu_Var<=user_nbins_hist_numu ) {
        if( map_Var_wiSel_Erec_bnbnu_numu_FC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_numu_FC.end() ) {
          array_Var[eff_index_numu_Var_FC -1] += 1;
        }
        if( map_Var_wiSel_Erec_bnbnu_numu_PC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_numu_PC.end() ) {
          array_Var[eff_index_numu_Var_PC -1] += 1;
        }       
      }
      if( index_basic_nue_Var>=1 && index_basic_nue_Var<=user_nbins_hist_nue ) {
        if( map_Var_wiSel_Erec_bnbnu_nue_FC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_nue_FC.end() ) {
          array_Var[eff_index_nue_Var_FC -1] += 1;        
        }
        if( map_Var_wiSel_Erec_bnbnu_nue_PC.find(roostr)!=map_Var_wiSel_Erec_bnbnu_nue_PC.end() ) {
          array_Var[eff_index_nue_Var_PC -1] += 1;
        }       
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
      //if(i!=j)
      if( fabs((*tt_matrix_cov_on_absdiff)(i,j))<1e-10 ) (*tt_matrix_cov_on_absdiff)(i,j) = 0;
      matrix_cov_on_absdiff(i,j) = (*tt_matrix_cov_on_absdiff)(i,j);
    }
  }
  for(int i=0; i<rows; i++) {// for content(bin) = 0, pay attention
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) {
      matrix_cov_on_absdiff(i,i) = 2;
    }
  }
  

  h1_diff_weight->Add( h1_weight, h1_sampling, 1, -1./ntoy );
    
  ///////////////////////////////////////////////

  delete h1_basic_nue;
  delete h1_basic_numu;
  
  delete[] array_test;
  delete[] array_CV;
  delete[] array_Var;  
  delete[] array_mean_Var2CV;

  /////////////////////////////////////////////// To handle N selected = 0, ...

  roostr = TString::Format("map_intrinsic_h1_CV_wiSel_%02d", idet);
  map_intrinsic_h1_CV_wiSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_CV_wiSel->GetBinContent( ibin );
    map_intrinsic_h1_CV_wiSel[idet]->SetBinContent( ibin, content );
  }
  
  roostr = TString::Format("map_intrinsic_h1_Var_wiSel_%02d", idet);
  map_intrinsic_h1_Var_wiSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_Var_wiSel->GetBinContent( ibin );
    map_intrinsic_h1_Var_wiSel[idet]->SetBinContent( ibin, content );
  }

  roostr = TString::Format("map_intrinsic_h1_CV_noSel_%02d", idet);
  map_intrinsic_h1_CV_noSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_CV_noSel->GetBinContent( ibin );
    map_intrinsic_h1_CV_noSel[idet]->SetBinContent( ibin, content );
  }
  
  roostr = TString::Format("map_intrinsic_h1_Var_noSel_%02d", idet);
  map_intrinsic_h1_Var_noSel[idet] = new TH1D(roostr, roostr, bins_basic, low_basic, hgh_basic);
  for(int ibin=1; ibin<=bins_basic; ibin++) {
    double content = h1_Var_noSel->GetBinContent( ibin );
    map_intrinsic_h1_Var_noSel[idet]->SetBinContent( ibin, content );
  }

  map_intrinsic_matrix_cov_on_AbsDiff_wiStat[idet].ResizeTo(rows, rows);
  map_intrinsic_matrix_cov_on_AbsDiff_wiStat[idet] = matrix_cov_on_absdiff;

  /////////////////////////////////////////////// generate covarinace matrix: from Abs.Diff + Error of Abs.Diff
  
  TPrincipal principal_obj(rows, "ND");  
  double *array_obj = new double[rows];

  TPrincipal principal_check(rows, "ND");
  
  TPrincipal principal_check_no_CovAbsDiff(rows, "ND"); 
  
  //////
  TMatrixDSym DSmatrix_cov_AbsDiff(rows);
  for(int ibin=0; ibin<rows; ibin++) {
    for(int jbin=0; jbin<rows; jbin++) {
      DSmatrix_cov_AbsDiff(ibin, jbin) = matrix_cov_on_absdiff(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov_AbsDiff );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TMatrixD matrix_eigenvector_T(rows, rows);
  matrix_eigenvector_T.Transpose( matrix_eigenvector );
  TMatrixD matrix_cov_AbsDiff_diag = matrix_eigenvector_T * matrix_cov_on_absdiff * matrix_eigenvector;
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      //if(i!=j)
      if( fabs(matrix_cov_AbsDiff_diag(i,j))<1e-10 ) matrix_cov_AbsDiff_diag(i,j) = 0;
    }
  }

  for( int itoy=1; itoy<=ntoy; itoy++ ) {
    
    if( itoy%(ntoy/10)==0 ) cout<<TString::Format(" ---> processing toy ( Covariance ): %4.2f, %6d", itoy*1./ntoy, itoy)<<endl;
    
    TRandom *random3 = new TRandom3(0);
    
    TMatrixD matrix_element(rows, 1);    
    for(int j=0; j<rows; j++) {
      matrix_element(j,0) = random3->Gaus( 0, sqrt( matrix_cov_AbsDiff_diag(j,j) ) );      
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;      
    double *array_variation_AbsDiff = new double[rows];    
    
    ///
    for(int idx=0; idx<rows; idx++) {
      array_variation_AbsDiff[idx] = matrix_variation(idx, 0);      
    }
    principal_check.AddRow( array_variation_AbsDiff );

    ///
    double rel_random = random3->Gaus(0, 1);

    ///
    double *array_check_no_CovAbsDiff = new double[rows];

    ///
    for(int ibin=1; ibin<=rows; ibin++) {
      double val_CV = h1_CV_wiSel->GetBinContent( ibin );
      double val_Var = h1_Var_wiSel->GetBinContent( ibin );
      double val_AbsDiff = val_Var - val_CV;

      ////
      array_check_no_CovAbsDiff[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;

      ////
      val_AbsDiff += array_variation_AbsDiff[ ibin-1 ];
      array_obj[ ibin-1 ] = val_CV + rel_random*val_AbsDiff;
    }

    //////
    principal_check_no_CovAbsDiff.AddRow( array_check_no_CovAbsDiff );
    principal_obj.AddRow( array_obj );
    
    //////
    delete random3;
    delete[] array_variation_AbsDiff;
    delete[] array_check_no_CovAbsDiff;    
  }// toy

  ////////
  TMatrixD *tt_matrix_cov_on_absdiff_check = (TMatrixD *)principal_check.GetCovarianceMatrix();
  matrix_cov_on_absdiff_check.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_on_absdiff_check)(i,j) = (*tt_matrix_cov_on_absdiff_check)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_cov_on_absdiff_check)(i,j))<1e-10 ) (*tt_matrix_cov_on_absdiff_check)(i,j) = 0;
      matrix_cov_on_absdiff_check(i,j) = (*tt_matrix_cov_on_absdiff_check)(i,j);
    }
  }
  
  ////////
  TMatrixD *tt_matrix_check_no_CovAbsDiff = (TMatrixD *)principal_check_no_CovAbsDiff.GetCovarianceMatrix();
  matrix_check_no_CovAbsDiff.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_check_no_CovAbsDiff)(i,j) = (*tt_matrix_check_no_CovAbsDiff)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_check_no_CovAbsDiff)(i,j))<1e-10 ) (*tt_matrix_check_no_CovAbsDiff)(i,j) = 0;
      matrix_check_no_CovAbsDiff(i,j) = (*tt_matrix_check_no_CovAbsDiff)(i,j);
    }
  }
  
  for(int i=0; i<rows; i++) {    
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) matrix_check_no_CovAbsDiff(i,i) = 2;
  }
  
  ////////
  TMatrixD *tt_matrix_cov_obj = (TMatrixD *)principal_obj.GetCovarianceMatrix();
  matrix_cov_obj.ResizeTo(rows, rows);  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      if(i<j) (*tt_matrix_cov_obj)(i,j) = (*tt_matrix_cov_obj)(j,i);
      //if(i!=j)
      if( fabs((*tt_matrix_cov_obj)(i,j))<1e-10 ) (*tt_matrix_cov_obj)(i,j) = 0;
      matrix_cov_obj(i,j) = (*tt_matrix_cov_obj)(i,j);
    }
  }

  for(int i=0; i<rows; i++) {    
    if( h1_CV_wiSel->GetBinContent( i+1 )==0 && h1_Var_wiSel->GetBinContent( i+1 )==0 ) matrix_cov_obj(i,i) = 2;
  }
  
  ///////////////////////////////////////////////
  
  map_intrinsic_matrix_cov_obj_wiStat[idet].ResizeTo(rows, rows);
  map_intrinsic_matrix_cov_obj_wiStat[idet] = matrix_cov_obj;
  
  map_intrinsic_matrix_cov_obj_noStat[idet].ResizeTo(rows, rows);
  map_intrinsic_matrix_cov_obj_noStat[idet] = matrix_check_no_CovAbsDiff;
  
  ///////////////////////////////////////////////
  
  delete[] array_obj;
}

///////////
///////////
///////////

void TDet::Exe(TString file_CV, TString file_Var, bool flag_numu, int flag_FC, int ntoy)
{  
  TString roostr = "";

  int entry(0), run(0), subrun(0), event(0);
  double weight(0);
  int status(0), post_generic(0), numuCC(0);
  double nueBDT(0);
  int is_FC(0);
  double Erec(0), Evis(0);

}

void TDet::Clear()
{
  cout<<" ---> begin clear"<<endl;
  
  num_DetVar++;
			     
  map_total_index_string.clear();
  map_total_index_CV_Erec.clear();
  map_total_index_Var_Erec.clear();    
  
  /// channels from bnbnu
  map_CV_wiSel_Erec_bnbnu_numu_FC.clear();
  map_Var_wiSel_Erec_bnbnu_numu_FC.clear();
  map_CV_wiSel_Erec_bnbnu_nue_FC.clear();
  map_Var_wiSel_Erec_bnbnu_nue_FC.clear();
  map_CV_wiSel_Erec_bnbnu_numu_PC.clear();
  map_Var_wiSel_Erec_bnbnu_numu_PC.clear();
  map_CV_wiSel_Erec_bnbnu_nue_PC.clear();
  map_Var_wiSel_Erec_bnbnu_nue_PC.clear();
  
  /// channels from intrinsic
  map_CV_wiSel_Erec_intrinsic_LEE_FC.clear();
  map_Var_wiSel_Erec_intrinsic_LEE_FC.clear();
  map_CV_wiSel_Erec_intrinsic_nue_FC.clear();
  map_Var_wiSel_Erec_intrinsic_nue_FC.clear();
  map_CV_wiSel_Erec_intrinsic_LEE_PC.clear();
  map_Var_wiSel_Erec_intrinsic_LEE_PC.clear();
  map_CV_wiSel_Erec_intrinsic_nue_PC.clear();
  map_Var_wiSel_Erec_intrinsic_nue_PC.clear();
  
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

    h1_subset_CV_wiSel = NULL;
    h1_subset_Var_wiSel = NULL;
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
  cout<<" ---> end clear"<<endl;
					     
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Plot_whole_cov(TMatrixD matrix_cov_obj_wiStat, TMatrixD matrix_cov_obj_noStat)
{
   
  TString roostr = "";
  TString str_suffix = "_zzWhole_";
  TString str_png = TString::Format("_%02d.png", 0);
  
  int rows = matrix_cov_obj_wiStat.GetNrows();
   
  ///////////////////////////////// matrix_cov_obj_noStat
  
  //////
  roostr = "h2_matrix_cov_obj_noStat";
  TH2D *h2_matrix_cov_obj_noStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);      
      h2_matrix_cov_obj_noStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat";
  TCanvas *canv_h2_matrix_cov_obj_noStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat->Draw("colz");
  h2_matrix_cov_obj_noStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat, "Index", "Index");
  h2_matrix_cov_obj_noStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_fractional";
  canv_h2_matrix_cov_obj_noStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_noStat_correlation";
  TH2D *h2_matrix_cov_obj_noStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);
      double val_i = sqrt( matrix_cov_obj_noStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_noStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_noStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_noStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat_correlation->Draw("colz");
  h2_matrix_cov_obj_noStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_noStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  canv_h2_matrix_cov_obj_noStat_correlation->SaveAs(roostr+str_png);
        
  ///////////////////////////////// matrix_cov_obj_wiStat
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat";
  TH2D *h2_matrix_cov_obj_wiStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);      
      h2_matrix_cov_obj_wiStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat";
  TCanvas *canv_h2_matrix_cov_obj_wiStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat->Draw("colz");
  h2_matrix_cov_obj_wiStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat, "Index", "Index");
  h2_matrix_cov_obj_wiStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_fractional";
  canv_h2_matrix_cov_obj_wiStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat_correlation";
  TH2D *h2_matrix_cov_obj_wiStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);
      double val_i = sqrt( matrix_cov_obj_wiStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_wiStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_wiStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_wiStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat_correlation->Draw("colz");
  h2_matrix_cov_obj_wiStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_wiStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  canv_h2_matrix_cov_obj_wiStat_correlation->SaveAs(roostr+str_png);

  
  /////////////////////////////////
  /////////////////////////////////

  roostr = "h1_sqrtCov_wiStat";
  TH1D *h1_sqrtCov_wiStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
  
  roostr = "h1_sqrtCov_noStat";
  TH1D *h1_sqrtCov_noStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);

  for(int ibin=1; ibin<=rows; ibin++) {
    double val_wiStat = sqrt( h2_matrix_cov_obj_wiStat->GetBinContent(ibin,ibin) );
    h1_sqrtCov_wiStat->SetBinContent(ibin, val_wiStat);
    
    // double content = h1_sqrtCov_wiStat->GetBinContent(ibin);
    // content = ((int)(content*1000+0.5))*1./1000;
    // h1_sqrtCov_wiStat->SetBinContent(ibin, content);
    
    double val_noStat = sqrt( h2_matrix_cov_obj_noStat->GetBinContent(ibin,ibin) );
    h1_sqrtCov_noStat->SetBinContent(ibin, val_noStat);      
  }
    
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  TCanvas *canv_h1_sqrtCov = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_sqrtCov, 0.15, 0.1,0.1,0.15);

  h1_sqrtCov_wiStat->Draw("hist");
  h1_sqrtCov_wiStat->SetLineColor(kBlue);
  h1_sqrtCov_wiStat->SetMarkerColor(kBlue);
  h1_sqrtCov_wiStat->SetMarkerSize(0.8);
  h1_sqrtCov_wiStat->SetTitle("");
  func_title_size(h1_sqrtCov_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_sqrtCov_wiStat, "Index", "Rel.Error");
  h1_sqrtCov_wiStat->GetXaxis()->SetNdivisions(506);
  h1_sqrtCov_wiStat->GetYaxis()->SetNdivisions(506);

  h1_sqrtCov_noStat->Draw("same hist");  
  h1_sqrtCov_noStat->SetLineColor(kRed);
  h1_sqrtCov_noStat->SetMarkerColor(kRed);
  
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  canv_h1_sqrtCov->SaveAs(roostr+str_png);
}



void Plot_subset( int idet,
		  TH1D *h1_CV_wiSel, TH1D *h1_Var_wiSel,
		  TMatrixD matrix_cov_on_absdiff,
		  TMatrixD matrix_cov_obj_wiStat, TMatrixD matrix_cov_obj_noStat
		  )
{
  
  TString roostr = "";
  TString str_suffix = "_zDet_";
  TString str_png = TString::Format("_%02d.png", idet);
  
  int rows = matrix_cov_on_absdiff.GetNrows();
  
  ///////////////////////////////// matrix_cov_on_absdiff
  
  //////
  roostr = "h2_matrix_cov_on_absdiff";
  TH2D *h2_matrix_cov_on_absdiff = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);      
      h2_matrix_cov_on_absdiff->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff";
  TCanvas *canv_h2_matrix_cov_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff->Draw("colz");
  h2_matrix_cov_on_absdiff->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff, "Index", "Index");
  h2_matrix_cov_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff";
  canv_h2_matrix_cov_on_absdiff->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_on_absdiff_fractional";
  TH2D *h2_matrix_cov_on_absdiff_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      h2_matrix_cov_on_absdiff_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_fractional";
  TCanvas *canv_h2_matrix_cov_on_absdiff_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff_fractional->Draw("colz");
  h2_matrix_cov_on_absdiff_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff_fractional, "Index", "Index");
  h2_matrix_cov_on_absdiff_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_fractional";
  canv_h2_matrix_cov_on_absdiff_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_on_absdiff_correlation";
  TH2D *h2_matrix_cov_on_absdiff_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);
      double val_i = sqrt( matrix_cov_on_absdiff(i,i) );
      double val_j = sqrt( matrix_cov_on_absdiff(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_on_absdiff_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_correlation";
  TCanvas *canv_h2_matrix_cov_on_absdiff_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff_correlation->Draw("colz");
  h2_matrix_cov_on_absdiff_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff_correlation, "Index", "Index");
  h2_matrix_cov_on_absdiff_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_correlation";
  canv_h2_matrix_cov_on_absdiff_correlation->SaveAs(roostr+str_png);
  
  ///////////////////////////////// matrix_cov_obj_noStat
  
  //////
  roostr = "h2_matrix_cov_obj_noStat";
  TH2D *h2_matrix_cov_obj_noStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);      
      h2_matrix_cov_obj_noStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat";
  TCanvas *canv_h2_matrix_cov_obj_noStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat->Draw("colz");
  h2_matrix_cov_obj_noStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat, "Index", "Index");
  h2_matrix_cov_obj_noStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat";
  canv_h2_matrix_cov_obj_noStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_noStat_fractional";
  TH2D *h2_matrix_cov_obj_noStat_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;
      h2_matrix_cov_obj_noStat_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_fractional";
  TCanvas *canv_h2_matrix_cov_obj_noStat_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat_fractional->Draw("colz");
  h2_matrix_cov_obj_noStat_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat_fractional, "Index", "Index");
  h2_matrix_cov_obj_noStat_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_fractional";
  canv_h2_matrix_cov_obj_noStat_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_noStat_correlation";
  TH2D *h2_matrix_cov_obj_noStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);
      double val_i = sqrt( matrix_cov_obj_noStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_noStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_noStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_noStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat_correlation->Draw("colz");
  h2_matrix_cov_obj_noStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_noStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  canv_h2_matrix_cov_obj_noStat_correlation->SaveAs(roostr+str_png);
      
  ///////////////////////////////// matrix_cov_obj_wiStat
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat";
  TH2D *h2_matrix_cov_obj_wiStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);      
      h2_matrix_cov_obj_wiStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat";
  TCanvas *canv_h2_matrix_cov_obj_wiStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat->Draw("colz");
  h2_matrix_cov_obj_wiStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat, "Index", "Index");
  h2_matrix_cov_obj_wiStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat";
  canv_h2_matrix_cov_obj_wiStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat_fractional";
  TH2D *h2_matrix_cov_obj_wiStat_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      h2_matrix_cov_obj_wiStat_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_fractional";
  TCanvas *canv_h2_matrix_cov_obj_wiStat_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat_fractional->Draw("colz");
  h2_matrix_cov_obj_wiStat_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat_fractional, "Index", "Index");
  h2_matrix_cov_obj_wiStat_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_fractional";
  canv_h2_matrix_cov_obj_wiStat_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat_correlation";
  TH2D *h2_matrix_cov_obj_wiStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);
      double val_i = sqrt( matrix_cov_obj_wiStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_wiStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_wiStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_wiStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat_correlation->Draw("colz");
  h2_matrix_cov_obj_wiStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_wiStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  canv_h2_matrix_cov_obj_wiStat_correlation->SaveAs(roostr+str_png);
    
  /////////////////////////////////
  /////////////////////////////////

  roostr = "h1_sqrtCov_wiStat";
  TH1D *h1_sqrtCov_wiStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
  
  roostr = "h1_sqrtCov_noStat";
  TH1D *h1_sqrtCov_noStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);

  for(int ibin=1; ibin<=rows; ibin++) {
    double val_wiStat = sqrt( h2_matrix_cov_obj_wiStat_fractional->GetBinContent(ibin,ibin) );
    if( h2_matrix_cov_obj_wiStat->GetBinContent(ibin,ibin)==0 ) val_wiStat = 0;
    h1_sqrtCov_wiStat->SetBinContent(ibin, val_wiStat);
    
    double val_noStat = sqrt( h2_matrix_cov_obj_noStat_fractional->GetBinContent(ibin,ibin) );
    if( h2_matrix_cov_obj_noStat->GetBinContent(ibin,ibin)==0 ) val_noStat = 0;
    h1_sqrtCov_noStat->SetBinContent(ibin, val_noStat);      
  }
    
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  TCanvas *canv_h1_sqrtCov = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_sqrtCov, 0.15, 0.1,0.1,0.15);

  h1_sqrtCov_wiStat->Draw("hist");
  h1_sqrtCov_wiStat->SetLineColor(kBlue);
  h1_sqrtCov_wiStat->SetTitle("");
  func_title_size(h1_sqrtCov_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_sqrtCov_wiStat, "Index", "Rel.Error");
  h1_sqrtCov_wiStat->GetXaxis()->SetNdivisions(506);
  h1_sqrtCov_wiStat->GetYaxis()->SetNdivisions(506);

  h1_sqrtCov_noStat->Draw("same hist");
  h1_sqrtCov_noStat->SetLineColor(kRed);
  
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  canv_h1_sqrtCov->SaveAs(roostr+str_png);
}
  

void Plot_samples( int idet, TString str_suffix,
		   TH1D *h1_weight, TH1D *h1_sampling,
		   TH1D *h1_CV_noSel, TH1D *h1_CV_wiSel,
		   TH1D *h1_Var_noSel, TH1D *h1_Var_wiSel,
		   TH1D *h1_diff_Var2CV, TH1D *h1_binomial_err,
		   TMatrixD matrix_cov_on_absdiff,
		   TMatrixD matrix_cov_on_absdiff_check,
		   TMatrixD matrix_cov_obj_wiStat,
		   TMatrixD matrix_cov_obj_noStat
		   )
{

  TString roostr = "";
  TString str_png = TString::Format("_%02d.png", idet);
  
  int rows = matrix_cov_on_absdiff.GetNrows();

  ////////////////////////////// weight
  
  roostr = "canv"+str_suffix+"h1_weight";
  TCanvas *canv_h1_weight = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_weight, 0.15, 0.1,0.1,0.15);
  h1_weight->Draw();
  h1_weight->SetLineColor(kBlack);
  func_title_size(h1_weight, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_weight, "Event index", "Weight");
  h1_weight->GetXaxis()->SetNdivisions(506);
  canv_h1_weight->SaveAs("canv"+str_suffix+"h1_weight"+str_png);
  roostr = "canv"+str_suffix+"h1_weight";
  canv_h1_weight->SaveAs(roostr+str_png);
  
  roostr = "canv"+str_suffix+"h1_sampling";
  TCanvas *canv_h1_sampling = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_sampling, 0.15, 0.1,0.1,0.15);
  h1_sampling->Draw();
  h1_sampling->SetLineColor(kBlack);
  func_title_size(h1_sampling, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_sampling, "Event index", "Weight");
  h1_sampling->GetXaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h1_sampling";
  canv_h1_sampling->SaveAs(roostr+str_png);
  
  ////////////////////////////// spectra
  
  roostr = "canv"+str_suffix+"spectra";
  TCanvas *canv_spectra = new TCanvas(roostr, roostr, 1600, 650);
  func_canv_margin(canv_spectra, 0.1, 0.15,0.1,0.15);
  canv_spectra->SetLogy();

  for(int ibin=1; ibin<=rows; ibin++) {
    double content = 0;

    content = h1_CV_noSel->GetBinContent(ibin);
    content = ( (int)(content*100+0.5) )*1./100;
    h1_CV_noSel->SetBinContent(ibin, content);
    
    content = h1_Var_noSel->GetBinContent(ibin);
    content = ( (int)(content*100+0.5) )*1./100;
    h1_Var_noSel->SetBinContent(ibin, content);
    
    content = h1_CV_wiSel->GetBinContent(ibin);
    content = ( (int)(content*100+0.5) )*1./100;
    h1_CV_wiSel->SetBinContent(ibin, content);
    
    content = h1_Var_wiSel->GetBinContent(ibin);
    content = ( (int)(content*100+0.5) )*1./100;
    h1_Var_wiSel->SetBinContent(ibin, content);
    
    content = h1_diff_Var2CV->GetBinContent(ibin);
    content = ( (int)(content*100+0.5) )*1./100;
    h1_diff_Var2CV->SetBinContent(ibin, content); 
  }
  
  h1_CV_noSel->Draw("hist text75");
  h1_CV_noSel->SetMinimum(0.5);
  h1_CV_noSel->SetTitle("");
  h1_CV_noSel->SetLineColor(kBlack);
  h1_CV_noSel->SetLineStyle(7);
  func_title_size(h1_CV_noSel, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_CV_noSel, "Index", "Entries");
  h1_CV_noSel->GetXaxis()->SetNdivisions(506);
  
  h1_Var_noSel->Draw("same hist");
  h1_Var_noSel->SetLineColor(kRed);
  h1_Var_noSel->SetLineStyle(7);

  h1_CV_wiSel->Draw("same hist text75");
  h1_CV_wiSel->SetLineColor(kBlack);
  h1_CV_wiSel->SetMarkerColor(kBlack);
  
  h1_Var_wiSel->Draw("same hist text15");
  h1_Var_wiSel->SetLineColor(kRed);
  h1_Var_wiSel->SetMarkerColor(kRed);

  h1_diff_Var2CV->Draw("same hist text0");
  h1_diff_Var2CV->SetLineColor(kBlue);
  h1_diff_Var2CV->SetMarkerColor(kBlue);

  h1_binomial_err->Draw("same hist text75");
  h1_binomial_err->SetLineColor(kOrange-3);
  h1_binomial_err->SetMarkerColor(kOrange-3);

  h1_CV_noSel->Draw("same axis");
        
  TLegend *lg = new TLegend(0.81+0.05, 0.5, 0.9+0.1, 0.9);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.04);
  lg->AddEntry( h1_CV_noSel, "CV before sel", "l" );
  lg->AddEntry( h1_CV_wiSel, "CV after sel", "l" );
  lg->AddEntry( "", "", "" );
  lg->AddEntry( h1_Var_noSel, "DV before sel", "l" );
  lg->AddEntry( h1_Var_wiSel, "DV after sel", "l" );
  lg->AddEntry( "", "", "" );
  lg->AddEntry( h1_diff_Var2CV, "Abs.Diff", "l" );
  lg->AddEntry( h1_binomial_err, "Err of Abs.Diff", "l" );    
  lg->Draw();

  roostr = "canv"+str_suffix+"spectra";
  canv_spectra->SaveAs(roostr+str_png);
  
  //////
  roostr = "h1_AbsDiff2CV";
  TH1D *h1_AbsDiff2CV = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
  for(int ibin=1; ibin<=rows; ibin++) {
    double val_diff = h1_diff_Var2CV->GetBinContent(ibin);
    double val_CV = h1_CV_wiSel->GetBinContent(ibin);
    h1_AbsDiff2CV->SetBinContent(ibin, val_diff/val_CV);
    // if( val_CV==0 ) {
    //   for(int idx=ibin; idx>=1; idx--) {
    // 	double content = h1_CV_wiSel->GetBinContent(idx);
    // 	val_CV = content;
    // 	h1_AbsDiff2CV->SetBinContent(ibin, val_diff/val_CV);
    // 	break;
    //   }
    // }
    if( val_CV==0 ) h1_AbsDiff2CV->SetBinContent(ibin, 0);

    double content = h1_AbsDiff2CV->GetBinContent(ibin);
    content = ((int)(content*1000+0.5))*1./1000;
    h1_AbsDiff2CV->SetBinContent(ibin, content);
  }
  roostr = "canv"+str_suffix+"h1_AbsDiff2CV";
  TCanvas *canv_h1_AbsDiff2CV = new TCanvas(roostr, "", 900, 650);
  func_canv_margin(canv_h1_AbsDiff2CV, 0.15, 0.1,0.1,0.15);
  h1_AbsDiff2CV->Draw("hist text75");
  h1_AbsDiff2CV->SetLineColor(kBlack);
  h1_AbsDiff2CV->SetMarkerSize(0.8);
  func_title_size(h1_AbsDiff2CV, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_AbsDiff2CV, "Event index", "Abs.Diff / CV");
  h1_AbsDiff2CV->GetXaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h1_AbsDiff2CV";
  canv_h1_AbsDiff2CV->SaveAs(roostr+str_png);

  ///////////////////////////////// matrix_cov_on_absdiff_check
  
  //////
  roostr = "h2_matrix_cov_on_absdiff_check";
  TH2D *h2_matrix_cov_on_absdiff_check = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff_check(i,j);      
      h2_matrix_cov_on_absdiff_check->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_check";
  TCanvas *canv_h2_matrix_cov_on_absdiff_check = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff_check, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff_check->Draw("colz");
  h2_matrix_cov_on_absdiff_check->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff_check, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff_check, "Index", "Index");
  h2_matrix_cov_on_absdiff_check->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff_check->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_check";
  canv_h2_matrix_cov_on_absdiff_check->SaveAs(roostr+str_png);

  ///////////////////////////////// matrix_cov_on_absdiff
  
  //////
  roostr = "h2_matrix_cov_on_absdiff";
  TH2D *h2_matrix_cov_on_absdiff = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);      
      h2_matrix_cov_on_absdiff->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff";
  TCanvas *canv_h2_matrix_cov_on_absdiff = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff->Draw("colz");
  h2_matrix_cov_on_absdiff->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff, "Index", "Index");
  h2_matrix_cov_on_absdiff->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff";
  canv_h2_matrix_cov_on_absdiff->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_on_absdiff_fractional";
  TH2D *h2_matrix_cov_on_absdiff_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      h2_matrix_cov_on_absdiff_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_fractional";
  TCanvas *canv_h2_matrix_cov_on_absdiff_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff_fractional->Draw("colz");
  h2_matrix_cov_on_absdiff_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff_fractional, "Index", "Index");
  h2_matrix_cov_on_absdiff_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_fractional";
  canv_h2_matrix_cov_on_absdiff_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_on_absdiff_correlation";
  TH2D *h2_matrix_cov_on_absdiff_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_on_absdiff(i,j);
      double val_i = sqrt( matrix_cov_on_absdiff(i,i) );
      double val_j = sqrt( matrix_cov_on_absdiff(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_on_absdiff_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_correlation";
  TCanvas *canv_h2_matrix_cov_on_absdiff_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_on_absdiff_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_on_absdiff_correlation->Draw("colz");
  h2_matrix_cov_on_absdiff_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_on_absdiff_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_on_absdiff_correlation, "Index", "Index");
  h2_matrix_cov_on_absdiff_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_on_absdiff_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_on_absdiff_correlation";
  canv_h2_matrix_cov_on_absdiff_correlation->SaveAs(roostr+str_png);
  
  ///////////////////////////////// comparison

  roostr = "h1_Cov2AbsDiff";
  TH1D *h1_Cov2AbsDiff = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);

  roostr = "h1_Cov2Berr";
  TH1D *h1_Cov2Berr = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);

  for(int ibin=1; ibin<=rows; ibin++) {
    double val_sqrt_cov = sqrt( matrix_cov_on_absdiff(ibin-1, ibin-1) );
    double val_AbsDiff = h1_diff_Var2CV->GetBinContent(ibin);
    double val_Berr = h1_binomial_err->GetBinContent(ibin);

    double val_Cov2AbsDiff = val_sqrt_cov/val_AbsDiff;
    if( val_AbsDiff==0 ) val_Cov2AbsDiff = 0;
    h1_Cov2AbsDiff->SetBinContent(ibin, val_Cov2AbsDiff);
    
    double val_Cov2Berr = val_sqrt_cov/val_Berr;
    if( val_Berr==0 ) val_Cov2Berr = 0;
    h1_Cov2Berr->SetBinContent(ibin, val_Cov2Berr);    
  }

  double max_h1_Cov2AbsDiff = h1_Cov2AbsDiff->GetMaximum();
  double max_h1_Cov2Berr = h1_Cov2Berr->GetMaximum();
  if( max_h1_Cov2AbsDiff<max_h1_Cov2Berr ) max_h1_Cov2AbsDiff = max_h1_Cov2Berr;
  
  roostr = "canv"+str_suffix+"h1_ratio_err";
  TCanvas *canv_h1_ratio_err = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_ratio_err, 0.15, 0.1,0.1,0.15);

  h1_Cov2AbsDiff->Draw("hist");
  h1_Cov2AbsDiff->SetMaximum(max_h1_Cov2AbsDiff * 1.3);
  h1_Cov2AbsDiff->SetLineColor(kBlue);
  h1_Cov2AbsDiff->SetTitle("");
  func_title_size(h1_Cov2AbsDiff, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_Cov2AbsDiff, "Index", "Ratio");
  h1_Cov2AbsDiff->GetXaxis()->SetNdivisions(506);
  h1_Cov2AbsDiff->GetYaxis()->SetNdivisions(506);

  TF1 *f1_ratio_err = new TF1("f1_ratio_err", "1", 0,  1e4);
  f1_ratio_err->Draw("same");
  f1_ratio_err->SetLineColor(kBlack);
  f1_ratio_err->SetLineStyle(7);
  
  //h1_Cov2Berr->Draw("same hist");
  h1_Cov2Berr->SetLineColor(kOrange-3);
          
  TLegend *lg_ratio_err = new TLegend(0.8-0.2, 0.75, 0.9, 0.9);
  lg_ratio_err->SetBorderSize(0);
  lg_ratio_err->SetFillStyle(0);
  lg_ratio_err->SetTextSize(0.05);
  lg_ratio_err->AddEntry( h1_Cov2AbsDiff, "sqrt(cov)/Abs.Diff", "l" );
  //lg_ratio_err->AddEntry( h1_Cov2Berr, "sqrt(cov)/Berr", "l" );  
  lg_ratio_err->Draw();

  roostr = "canv"+str_suffix+"h1_ratio_err";
  canv_h1_ratio_err->SaveAs(roostr+str_png);
  
  ///////////////////////////////// matrix_cov_obj_noStat
  
  //////
  roostr = "h2_matrix_cov_obj_noStat";
  TH2D *h2_matrix_cov_obj_noStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);      
      h2_matrix_cov_obj_noStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat";
  TCanvas *canv_h2_matrix_cov_obj_noStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat->Draw("colz");
  h2_matrix_cov_obj_noStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat, "Index", "Index");
  h2_matrix_cov_obj_noStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat";
  canv_h2_matrix_cov_obj_noStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_noStat_fractional";
  TH2D *h2_matrix_cov_obj_noStat_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      h2_matrix_cov_obj_noStat_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_fractional";
  TCanvas *canv_h2_matrix_cov_obj_noStat_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat_fractional->Draw("colz");
  h2_matrix_cov_obj_noStat_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat_fractional, "Index", "Index");
  h2_matrix_cov_obj_noStat_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_fractional";
  canv_h2_matrix_cov_obj_noStat_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_noStat_correlation";
  TH2D *h2_matrix_cov_obj_noStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_noStat(i,j);
      double val_i = sqrt( matrix_cov_obj_noStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_noStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_noStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_noStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_noStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_noStat_correlation->Draw("colz");
  h2_matrix_cov_obj_noStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_noStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_noStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_noStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_noStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_noStat_correlation";
  canv_h2_matrix_cov_obj_noStat_correlation->SaveAs(roostr+str_png);
      
  ///////////////////////////////// matrix_cov_obj_wiStat
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat";
  TH2D *h2_matrix_cov_obj_wiStat = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);      
      h2_matrix_cov_obj_wiStat->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat";
  TCanvas *canv_h2_matrix_cov_obj_wiStat = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat->Draw("colz");
  h2_matrix_cov_obj_wiStat->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat, "Index", "Index");
  h2_matrix_cov_obj_wiStat->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat";
  canv_h2_matrix_cov_obj_wiStat->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat_fractional";
  TH2D *h2_matrix_cov_obj_wiStat_fractional = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);
      double val_i = h1_CV_wiSel->GetBinContent(i+1);
      double val_j = h1_CV_wiSel->GetBinContent(j+1);
      
      // if( val_i==0 ) {
      // 	for(int idx=i; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_i = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      // if( val_j==0 ) {
      // 	for(int idx=j; idx>=0; idx--) {
      // 	  double content2 = h1_CV_wiSel->GetBinContent(idx+1);
      // 	  if( content2>0 ) {
      // 	    val_j = content2;
      // 	    break;
      // 	  }
      // 	}
      // }
      content = content/val_i/val_j;
      if(val_i==0 || val_j==0) content = 0;      
      h2_matrix_cov_obj_wiStat_fractional->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_fractional";
  TCanvas *canv_h2_matrix_cov_obj_wiStat_fractional = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat_fractional, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat_fractional->Draw("colz");
  h2_matrix_cov_obj_wiStat_fractional->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat_fractional, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat_fractional, "Index", "Index");
  h2_matrix_cov_obj_wiStat_fractional->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat_fractional->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_fractional";
  canv_h2_matrix_cov_obj_wiStat_fractional->SaveAs(roostr+str_png);
  
  //////
  roostr = "h2_matrix_cov_obj_wiStat_correlation";
  TH2D *h2_matrix_cov_obj_wiStat_correlation = new TH2D(roostr, roostr, rows, 0.5, rows+0.5, rows, 0.5, rows+0.5);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<rows; j++) {
      double content = matrix_cov_obj_wiStat(i,j);
      double val_i = sqrt( matrix_cov_obj_wiStat(i,i) );
      double val_j = sqrt( matrix_cov_obj_wiStat(j,j) );
      content = content/val_i/val_j;
      h2_matrix_cov_obj_wiStat_correlation->SetBinContent(i+1, j+1, content);
    }
  }
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  TCanvas *canv_h2_matrix_cov_obj_wiStat_correlation = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h2_matrix_cov_obj_wiStat_correlation, 0.15, 0.2,0.1,0.15);
  h2_matrix_cov_obj_wiStat_correlation->Draw("colz");
  h2_matrix_cov_obj_wiStat_correlation->SetTitle("");
  func_title_size(h2_matrix_cov_obj_wiStat_correlation, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_matrix_cov_obj_wiStat_correlation, "Index", "Index");
  h2_matrix_cov_obj_wiStat_correlation->GetXaxis()->SetNdivisions(506);
  h2_matrix_cov_obj_wiStat_correlation->GetYaxis()->SetNdivisions(506);
  roostr = "canv"+str_suffix+"h2_matrix_cov_obj_wiStat_correlation";
  canv_h2_matrix_cov_obj_wiStat_correlation->SaveAs(roostr+str_png);

  /////////////////////////////////
  /////////////////////////////////

  roostr = "h1_sqrtCov_wiStat";
  TH1D *h1_sqrtCov_wiStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);
  
  roostr = "h1_sqrtCov_noStat";
  TH1D *h1_sqrtCov_noStat = new TH1D(roostr, roostr, rows, 0.5, rows+0.5);

  for(int ibin=1; ibin<=rows; ibin++) {
    double val_wiStat = sqrt( h2_matrix_cov_obj_wiStat_fractional->GetBinContent(ibin,ibin) );
    if( h2_matrix_cov_obj_wiStat->GetBinContent(ibin,ibin)==0 ) val_wiStat = 0;
    h1_sqrtCov_wiStat->SetBinContent(ibin, val_wiStat);
    
    double val_noStat = sqrt( h2_matrix_cov_obj_noStat_fractional->GetBinContent(ibin,ibin) );
    if( h2_matrix_cov_obj_noStat->GetBinContent(ibin,ibin)==0 ) val_noStat = 0;
    h1_sqrtCov_noStat->SetBinContent(ibin, val_noStat);      
  }
    
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  TCanvas *canv_h1_sqrtCov = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_sqrtCov, 0.15, 0.1,0.1,0.15);

  h1_sqrtCov_wiStat->Draw("hist");
  h1_sqrtCov_wiStat->SetLineColor(kBlue);
  h1_sqrtCov_wiStat->SetTitle("");
  func_title_size(h1_sqrtCov_wiStat, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_sqrtCov_wiStat, "Index", "Rel.Error");
  h1_sqrtCov_wiStat->GetXaxis()->SetNdivisions(506);
  h1_sqrtCov_wiStat->GetYaxis()->SetNdivisions(506);

  h1_sqrtCov_noStat->Draw("same hist");
  h1_sqrtCov_noStat->SetLineColor(kRed);
  
  roostr = "canv"+str_suffix+"h1_sqrtCov";
  canv_h1_sqrtCov->SaveAs(roostr+str_png);

  //////////////////////////////////////////////////////////
  
  // TFile *outfile = new TFile("outfile_testBB.root", "recreate");
  // h1_CV_noSel->Write();
  // h1_CV_wiSel->Write();
  // h1_Var_noSel->Write();
  // h1_Var_wiSel->Write();
  // h1_AbsDiff2CV->Write();

  // h1_Cov2AbsDiff->Write();
  // h1_sqrtCov_wiStat->Write();
  
  // outfile->Close();
    
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

  map<int, TString>map_bnbnu_CV;
  map_bnbnu_CV[1] = "./data_det_syst_total/common_detvar/LYDown/numu_CV";
  map_bnbnu_CV[2] = "./data_det_syst_total/common_detvar/Recomb2/numu_CV";
  map_bnbnu_CV[3] = "./data_det_syst_total/common_detvar/SCE/numu_CV";
  map_bnbnu_CV[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_CV";
  map_bnbnu_CV[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_CV";
  map_bnbnu_CV[6] = "./data_det_syst_total/common_detvar/WireModX/numu_CV";
  map_bnbnu_CV[7] = "./data_det_syst_total/common_detvar/WireModYZ/numu_CV";  
  map<int, TString>map_bnbnu_Var;
  map_bnbnu_Var[1] = "./data_det_syst_total/common_detvar/LYDown/numu_Var";
  map_bnbnu_Var[2] = "./data_det_syst_total/common_detvar/Recomb2/numu_Var";
  map_bnbnu_Var[3] = "./data_det_syst_total/common_detvar/SCE/numu_Var";
  map_bnbnu_Var[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/numu_Var";
  map_bnbnu_Var[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/numu_Var";
  map_bnbnu_Var[6] = "./data_det_syst_total/common_detvar/WireModX/numu_Var";
  map_bnbnu_Var[7] = "./data_det_syst_total/common_detvar/WireModYZ/numu_Var";
  
  map<int, TString>map_intrinsic_CV;
  map_intrinsic_CV[1] = "./data_det_syst_total/common_detvar/LYDown/nue_CV";
  map_intrinsic_CV[2] = "./data_det_syst_total/common_detvar/Recomb2/nue_CV";
  map_intrinsic_CV[3] = "./data_det_syst_total/common_detvar/SCE/nue_CV";
  map_intrinsic_CV[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/nue_CV";
  map_intrinsic_CV[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/nue_CV";
  map_intrinsic_CV[6] = "./data_det_syst_total/common_detvar/WireModX/nue_CV";
  map_intrinsic_CV[7] = "./data_det_syst_total/common_detvar/WireModYZ/nue_CV";  
  map<int, TString>map_intrinsic_Var;
  map_intrinsic_Var[1] = "./data_det_syst_total/common_detvar/LYDown/nue_Var";
  map_intrinsic_Var[2] = "./data_det_syst_total/common_detvar/Recomb2/nue_Var";
  map_intrinsic_Var[3] = "./data_det_syst_total/common_detvar/SCE/nue_Var";
  map_intrinsic_Var[4] = "./data_det_syst_total/common_detvar/WireModThetaXZ/nue_Var";
  map_intrinsic_Var[5] = "./data_det_syst_total/common_detvar/WireModThetaYZ/nue_Var";
  map_intrinsic_Var[6] = "./data_det_syst_total/common_detvar/WireModX/nue_Var";
  map_intrinsic_Var[7] = "./data_det_syst_total/common_detvar/WireModYZ/nue_Var";
  
  //
  TDet *det_test = new TDet();
  
  det_test->Set_nue_hist(14, 150, 1550);
  
  det_test->Set_numu_hist(14, 150, 1550);

  det_test->ntoy = 100;

  for(int idet=1; idet<=1; idet++) {

    cout<<endl<<TString::Format(" ------------------> the %2d Detector systematics", idet)<<endl<<endl;
    
    cout<<endl<<" ---> processing DetSyst of bnbnu: "<<map_bnbnu_CV[idet]<<endl<<endl;    
    det_test->Clear();    
    det_test->Exe_channels_from_bnbnu( map_bnbnu_CV[idet], map_bnbnu_Var[idet], idet, det_test->ntoy);
    Plot_samples( idet, "_bnbnu_",
    		  det_test->h1_weight, det_test->h1_sampling,
    		  det_test->h1_CV_noSel, det_test->h1_CV_wiSel,
    		  det_test->h1_Var_noSel, det_test->h1_Var_wiSel,
    		  det_test->h1_diff_Var2CV, det_test->h1_binomial_err,
    		  det_test->matrix_cov_on_absdiff,
    		  det_test->matrix_cov_on_absdiff_check,
    		  det_test->matrix_cov_obj,
    		  det_test->matrix_check_no_CovAbsDiff
    		  );   
    
    
    cout<<endl<<" ---> processing DetSyst of intrinsic: "<<map_intrinsic_CV[idet]<<endl<<endl;
    det_test->Clear();
    det_test->Exe_channels_from_intrinsic( map_intrinsic_CV[idet], map_intrinsic_Var[idet], idet, det_test->ntoy);
    Plot_samples( idet, "_intrinsic_",
    		  det_test->h1_weight, det_test->h1_sampling,
    		  det_test->h1_CV_noSel, det_test->h1_CV_wiSel,
    		  det_test->h1_Var_noSel, det_test->h1_Var_wiSel,
    		  det_test->h1_diff_Var2CV, det_test->h1_binomial_err,
    		  det_test->matrix_cov_on_absdiff,
    		  det_test->matrix_cov_on_absdiff_check,
    		  det_test->matrix_cov_obj,
    		  det_test->matrix_check_no_CovAbsDiff
    		  ); 
    
   
    cout<<endl<<" ---> producing covariance matrix"<<endl<<endl;
    det_test->Clear();
    det_test->Exe_obj_cov(idet, det_test->ntoy);
    Plot_subset( idet,
    		 det_test->h1_subset_CV_wiSel, det_test->h1_subset_Var_wiSel,
    		 det_test->matrix_cov_on_absdiff,
    		 det_test->matrix_subset_cov_obj_wiStat,
    		 det_test->matrix_subset_cov_obj_noStat
    		 );
    
  }// idet 
  
  ////////
  det_test->Exe_whole_cov();
  Plot_whole_cov( det_test->matrix_whole_cov_obj_wiStat_fractional,
		  det_test->matrix_whole_cov_obj_noStat_fractional
		  );

  
  ////////
  TFile *outfile_det = new TFile("outfile_det_wiStat.root", "recreate");
  
  (det_test->matrix_whole_cov_obj_wiStat_fractional).Write("matrix_covariance_fractional_wiStat");
  (det_test->matrix_whole_cov_obj_noStat_fractional).Write("matrix_covariance_fractional_noStat");
    
  for(int idet=1; idet<=7; idet++) {
    (det_test->map_matrix_subset_cov_obj_wiStat_fractional[idet]).Write(TString::Format("matrix_covFrac_det_%02d_wiStat", idet));
    (det_test->map_matrix_subset_cov_obj_noStat_fractional[idet]).Write(TString::Format("matrix_covFrac_det_%02d_noStat", idet));
  }

  outfile_det->Close();
  
}
