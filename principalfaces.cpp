//root
#include "TPrincipal.h"
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"

//C, C++
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

struct face_data_str {
  vector<vector<Double_t>> data_v;
  vector<vector<Double_t>> data_norm_v;
};

void read_face_data(TString name, vector<face_data_str> &face_data_v);
void load_U_S_Vh_data(TString name_U, TString name_S, TString name_Vh);

void fill_2Dhist(vector<TH2D*> &h2_faces_v, vector<TH2D*> &h2_faces_norm_v, const vector<face_data_str> &face_data_v);
void fill_1Dhist(vector<vector<TH1D*>> &h1_faces_v, vector<vector<TH1D*>> &h1_faces_norm_v, const vector<face_data_str> &face_data_v);
void fill_2Dhist_principal_components(vector<TH2D*> &h2_principal_components_v);

void get_mean(const vector<face_data_str> &face_data_v, vector<vector<Double_t>> &face_data_mean);
void get_std(const vector<face_data_str> &face_data_v, const vector<vector<Double_t>> &face_data_mean, vector<vector<Double_t>> &face_data_std);
void get_norm(vector<face_data_str> &face_data_v,
	      const vector<vector<Double_t>> &face_data_mean,
	      const vector<vector<Double_t>> &face_data_std);

void draw_faces( Int_t nn, const vector<TH2D*> &h2_v, TCanvas *c1);

void reshape_data(const vector<face_data_str> &face_data_v, vector<vector<Double_t>> &data);

void get_covariance_data( const vector<vector<Double_t>> &data, vector<vector<Double_t>> &data_sigma, TPrincipal* principal);

void save_to_csv_data(TString file_out_name, const vector<vector<Double_t>> &data);

vector<vector<Double_t>> data_sigma;

const Int_t nn_im = 5000;
const Int_t dd_im = 1024;

Double_t data_U[nn_im][nn_im];
Double_t data_S[dd_im];
Double_t data_Vh[dd_im][dd_im];

Int_t main(){
  
  vector<face_data_str> face_data_v;
  read_face_data("ex7faces.cvs", face_data_v);
  //
  vector<vector<Double_t>> face_data_mean;
  vector<vector<Double_t>> face_data_std;
  //
  get_mean(face_data_v, face_data_mean);
  get_std(face_data_v, face_data_mean, face_data_std);
  get_norm(face_data_v, face_data_mean, face_data_std);
  
  vector<TH2D*> h2_faces_v;
  vector<TH2D*> h2_faces_norm_v;
  fill_2Dhist(h2_faces_v, h2_faces_norm_v, face_data_v);

  vector<vector<TH1D*>> h1_faces_v;
  vector<vector<TH1D*>> h1_faces_norm_v;
  fill_1Dhist(h1_faces_v, h1_faces_norm_v, face_data_v);

  vector<vector<Double_t>> data;
  reshape_data(face_data_v, data);

  save_to_csv_data("data.csv", data);

  load_U_S_Vh_data("U.cvs","S.cvs","Vh.cvs");
  
  //vector<vector<Double_t>> data_sigma(data.size(), vector<Double_t>(data.size()));
  //TPrincipal* principal = new TPrincipal(data.size(),"N");
  //get_covariance_data( data, data_sigma, principal);
  
  // Do the actual analysis
  //principal->MakePrincipals();

  // Print out the result on
  //principal->Print();

  //const double *data_EigenVectors_v = principal->GetEigenVectors()->GetMatrixArray();

  vector<TH2D*> h2_principal_components_v;
  fill_2Dhist_principal_components(h2_principal_components_v);
  
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  draw_faces(10,h2_faces_v,c1);
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1000,1000);
  draw_faces(10,h2_faces_norm_v,c2);
  TCanvas *c3 = new TCanvas("c3","c3",10,10,1000,1000);
  draw_faces(10,h2_principal_components_v,c3);
  //
  
  //----------------------
  TFile* rootFile = new TFile("hist.root", "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<"hist.root"<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<"hist.root"<<endl;

  for(unsigned int i = 0; i < h2_faces_v.size(); i++)
    h2_faces_v.at(i)->Write();
  for(unsigned int i = 0; i < h2_faces_norm_v.size(); i++)
    h2_faces_norm_v.at(i)->Write();
  for(unsigned int i = 0; i < h1_faces_v.size(); i++)
    for(unsigned int j = 0; j < h1_faces_v.at(i).size(); j++)
      h1_faces_v.at(i).at(j)->Write();
  for(unsigned int i = 0; i < h1_faces_norm_v.size(); i++)
    for(unsigned int j = 0; j < h1_faces_norm_v.at(i).size(); j++)
      h1_faces_norm_v.at(i).at(j)->Write();

  for(unsigned int i = 0; i < h2_principal_components_v.size(); i++)
    h2_principal_components_v.at(i)->Write();
  
  c1->Write();
  c2->Write();
  c3->Write();
  
  rootFile->Close();

  
  /*
  Int_t jj = 0;
  
  ifstream file_cvs ();
  if(file_cvs.is_open()){
    for(Int_t i = 0; i < m; i++){
      file_cvs>>evID;
      TString name = "h2";
      name += "_id";
      name += i;
      TH2D *h2 = new TH2D(name.Data(),name.Data(),l,0,l,l,0,l);
      for(Int_t j = 0; j < n; j++){
	file_cvs >> val;
	data[j] = val;
	if(jj == 32)
	jj = 0;
	h2->SetBinContent(32-j/32+1,32-jj+1,val);
	jj++;	
      }
      h2_faces_v.push_back(h2);
      //break;
      principal->AddRow(data);
    }
    file_cvs.close();
  }
  //
  gStyle->SetPalette(kGreyScale);
  h2_faces_v.at(0)->Draw("ZCOLOR");
  */
 
  /*
  //
  // Do the actual analysis
  principal->MakePrincipals();
  
  // Print out the result on
  principal->Print();

  Int_t ii_vec = 1;
  TString name = "h2";
  name += "_EigenVec";
  name += ii_vec ;  
  TH2D *h2 = new TH2D(name.Data(),name.Data(),l,0,l,l,0,l);

  const double *data_v = principal->GetEigenVectors()->GetMatrixArray();

  jj=0;
  for(Int_t j = 0; j < n; j++){
    if(jj == 32)
      jj = 0;
    h2->SetBinContent(32-j/32+1,32-jj+1,data_v[j+ii_vec*n]);
    jj++;	
  }
  
  cout<<principal->GetEigenVectors()->GetNcols()<<endl
      <<principal->GetEigenVectors()->GetNrows()<<endl;
  
  h2->Draw("ZCOLOR");
  */
  //const TMatrixD* eigenV = principal->GetEigenVectors();

  // = new double[4];
  //auto data_v = new double [2][2];
  //double data_v[2][2];
  /*

  const TVectorD* V_v  = principal->GetEigenValues();

  const double *data_v_v = V_v->GetMatrixArray();
  //
  cout<<V_v->GetNrows()<<endl;
  //
  cout<<data_v_v[0]<<endl;
  cout<<data_v[0]<<endl;
  cout<<data_v[1]<<endl;
  cout<<data_v_v[1]<<endl;
  cout<<data_v[2]<<endl;
  cout<<data_v[3]<<endl;

  const double *cov_v = principal->GetCovarianceMatrix()->GetMatrixArray();

  cout<<principal->GetCovarianceMatrix()->GetNcols()<<endl
      <<principal->GetCovarianceMatrix()->GetNrows()<<endl;

  A.Use(m,n,data_data);
  B.Use(m,n,data_data);
  
  cout<<A.GetNcols()<<endl
      <<A.GetNrows()<<endl;

  cout<<B.GetNcols()<<endl
      <<B.GetNrows()<<endl;
    
  //TMatrixD C = A*(B.Transpose());
  TMatrixD C(n,m);
  C.Transpose(B);
  //B.Transpose();

  TMatrixD D(m,m);

  D = A*C;

  TPrincipal* principal_n = new TPrincipal(m,"N");
  
  cout<<C.GetNcols()<<endl
      <<C.GetNrows()<<endl;

  cout<<D.GetNcols()<<endl
      <<D.GetNrows()<<endl;

  
  //const Double_t pData[2][2] = ;
  //(0,0);
  //
  //  
  //GetNcols
  //

  //principal->
  

  // Test the PCA 
  principal->Test();

  // Make some histograms of the orginal, principal, residue, etc data 
  //principal->MakeHistograms();
  
  // Make two functions to map between feature and pattern space 
  //principal->MakeCode();

  // Start a browser, so that we may browse the histograms generated
  // above 
  //TBrowser* b = new TBrowser("principalBrowser", principal);

  gr->Draw("AP");
  */
  
  return 0;
  
}

void fill_2Dhist( vector<TH2D*> &h2_faces_v, vector<TH2D*> &h2_faces_norm_v, const vector<face_data_str> &face_data_v){
  //cout<<"face_data_v.size()              = "<<face_data_v.size()<<endl
  //  <<"face_data_v.at(0).data_v.size() = "<<face_data_v.at(0).data_v.size()<<endl
  //  <<"face_data_v.at(0).data_v.at(0).size() = "<<face_data_v.at(0).data_v.at(0).size()<<endl;
  unsigned int l = face_data_v.at(0).data_v.size();
  for(unsigned int i = 0; i < face_data_v.size(); i++){  
    TString name = "h2";
    name += "_ID";
    name += i;  
    TH2D *h2 = new TH2D(name.Data(),name.Data(),
			face_data_v.at(i).data_v.size(),0,face_data_v.at(i).data_v.size(),
			face_data_v.at(i).data_v.size(),0,face_data_v.at(i).data_v.size());
    //
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++){
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++){ 
	h2->SetBinContent(l-j,l-k,face_data_v.at(i).data_v.at(j).at(k));
      }
    }
    h2_faces_v.push_back(h2);
  }
  for(unsigned int i = 0; i < face_data_v.size(); i++){  
    TString name = "h2_norm";
    name += "_ID";
    name += i;  
    TH2D *h2 = new TH2D(name.Data(),name.Data(),
			face_data_v.at(i).data_norm_v.size(),0,face_data_v.at(i).data_norm_v.size(),
			face_data_v.at(i).data_norm_v.size(),0,face_data_v.at(i).data_norm_v.size());
    //
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++){
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++){ 
	h2->SetBinContent(l-j,l-k,face_data_v.at(i).data_norm_v.at(j).at(k));
      }
    }
    h2_faces_norm_v.push_back(h2);
  }
}

void fill_1Dhist(vector<vector<TH1D*>> &h1_faces_v, vector<vector<TH1D*>> &h1_faces_norm_v, const vector<face_data_str> &face_data_v){
  //
  for(unsigned int j = 0; j < face_data_v.at(0).data_v.size(); j++){
    vector<TH1D*> h1v;
    for(unsigned int k = 0; k < face_data_v.at(0).data_v.at(j).size(); k++){
      TString name = "h1";
      name += "_j";
      name += j;  
      name += "_k";
      name += k;  
      TH1D *h1 = new TH1D(name.Data(),name.Data(),200,-400,400);
      for(unsigned int i = 0; i < face_data_v.size(); i++)
	h1->Fill(face_data_v.at(i).data_v.at(j).at(k));
      h1v.push_back(h1);
    }
    h1_faces_v.push_back(h1v); 
  }
  //
  for(unsigned int j = 0; j < face_data_v.at(0).data_norm_v.size(); j++){
    vector<TH1D*> h1v;
    for(unsigned int k = 0; k < face_data_v.at(0).data_norm_v.at(j).size(); k++){
      TString name = "h1_norm";
      name += "_j";
      name += j;  
      name += "_k";
      name += k;  
      TH1D *h1 = new TH1D(name.Data(),name.Data(),200,-10,10);
      for(unsigned int i = 0; i < face_data_v.size(); i++)
	h1->Fill(face_data_v.at(i).data_norm_v.at(j).at(k));
      h1v.push_back(h1);
    }
    h1_faces_norm_v.push_back(h1v); 
  }
}

void load_U_S_Vh_data(TString name_U, TString name_S, TString name_Vh){
  ifstream file_cvs_U(name_U);
  ifstream file_cvs_S(name_S);
  ifstream file_cvs_Vh(name_Vh);
  if(file_cvs_U.is_open())
    for(Int_t i = 0; i <nn_im; i++)
      for(Int_t j = 0; j < nn_im; j++)
	file_cvs_U>>data_U[i][j];
  file_cvs_U.close();
  //
  if(file_cvs_S.is_open())
    for(Int_t i = 0; i <dd_im; i++)
      file_cvs_S>>data_S[i];
  file_cvs_S.close();
  //
  if(file_cvs_Vh.is_open())
    for(Int_t i = 0; i <dd_im; i++)
      for(Int_t j = 0; j < dd_im; j++)
	file_cvs_Vh>>data_Vh[i][j];
  file_cvs_Vh.close();
}

void read_face_data(TString name, vector<face_data_str> &face_data_v){
  Int_t n=dd_im;
  Int_t l=(Int_t)TMath::Sqrt(n);
  Int_t m=nn_im;
  Int_t evID;
  ifstream file_cvs(name);
  Double_t val;
  if(file_cvs.is_open()){
    for(Int_t i = 0; i < m; i++){
      file_cvs>>evID;
      face_data_str face_str;
      vector<Double_t> val_v;
      vector<Double_t> val_v_zero;
      for(Int_t j = 0; j < n; j++){
	file_cvs>>val;
	val_v.push_back(val);
	val_v_zero.push_back(0);
	if(val_v.size() == (unsigned int)l){
	  face_str.data_v.push_back(val_v);
	  face_str.data_norm_v.push_back(val_v_zero);
	  val_v.clear();
	  val_v_zero.clear();
	}
      }
      face_data_v.push_back(face_str);
    }
    file_cvs.close();
  }
  //
  //gStyle->SetPalette(kGreyScale);
  //h2_faces_v.at(0)->Draw("ZCOLOR"); 
}

void get_mean(const vector<face_data_str> &face_data_v, vector<vector<Double_t>> &face_data_mean){
  for(unsigned int j = 0; j < face_data_v.at(0).data_v.size(); j++){
    vector<Double_t> val;
    for(unsigned int k = 0; k < face_data_v.at(0).data_v.at(j).size(); k++)
      val.push_back(0.0);
    face_data_mean.push_back(val);
  }
  for(unsigned int i = 0; i < face_data_v.size(); i++)
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++)
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++)
	face_data_mean.at(j).at(k) += face_data_v.at(i).data_v.at(j).at(k);
  if(face_data_v.size()>0)
    for(unsigned int j = 0; j < face_data_v.at(0).data_v.size(); j++)
      for(unsigned int k = 0; k < face_data_v.at(0).data_v.at(j).size(); k++)
	face_data_mean.at(j).at(k) /= face_data_v.size();
}

void get_std(const vector<face_data_str> &face_data_v, const vector<vector<Double_t>> &face_data_mean, vector<vector<Double_t>> &face_data_std){
  for(unsigned int j = 0; j < face_data_v.at(0).data_v.size(); j++){
    vector<Double_t> val;
    for(unsigned int k = 0; k < face_data_v.at(0).data_v.at(j).size(); k++)
      val.push_back(0.0);
    face_data_std.push_back(val);
  }
  for(unsigned int i = 0; i < face_data_v.size(); i++)
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++)
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++)
	face_data_std.at(j).at(k) += TMath::Power((face_data_v.at(i).data_v.at(j).at(k) - face_data_mean.at(j).at(k)),2);
  if(face_data_v.size()>0)
    for(unsigned int j = 0; j < face_data_v.at(0).data_v.size(); j++)
      for(unsigned int k = 0; k < face_data_v.at(0).data_v.at(j).size(); k++)
	face_data_std.at(j).at(k) = TMath::Sqrt(face_data_std.at(j).at(k)/face_data_v.size());
}

void get_norm(vector<face_data_str> &face_data_v,
	      const vector<vector<Double_t>> &face_data_mean,
	      const vector<vector<Double_t>> &face_data_std){
  for(unsigned int i = 0; i < face_data_v.size(); i++){
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++){
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++){
	if(face_data_std.at(j).at(k)>0)
	  face_data_v.at(i).data_norm_v.at(j).at(k) = (face_data_v.at(i).data_v.at(j).at(k) - face_data_mean.at(j).at(k))/face_data_std.at(j).at(k);
      }
    }
  }
}

void draw_faces( Int_t nn, const vector<TH2D*> &h2_v, TCanvas *c1){
  c1->cd();
  gStyle->SetPalette(kGreyScale);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  c1->Divide(nn,nn,0.0001,0.0001,0);
  Int_t padID = 0;
  Int_t faceID = 0;
  for(Int_t i = 0;i<nn;i++){
    for(Int_t j = 0;j<nn;j++){
      padID = i*nn+j+1;
      faceID = i*nn+j;
      c1->cd(padID);
      gStyle->SetOptStat(kFALSE);
      gPad->SetLeftMargin(0.0);
      gPad->SetTopMargin(0.0);
      gPad->SetBottomMargin(0.0);
      gPad->SetRightMargin(0.0);
      h2_v.at(faceID)->SetStats(0);
      h2_v.at(faceID)->SetTitle("");
      h2_v.at(faceID)->GetXaxis()->SetLabelSize(0);
      h2_v.at(faceID)->GetYaxis()->SetLabelSize(0);
      h2_v.at(faceID)->GetXaxis()->SetLabelOffset(999);
      h2_v.at(faceID)->GetYaxis()->SetLabelOffset(999);
      h2_v.at(faceID)->GetXaxis()->SetTickLength(0);
      h2_v.at(faceID)->GetYaxis()->SetTickLength(0);
      h2_v.at(faceID)->GetZaxis()->SetTickLength(0);
      h2_v.at(faceID)->Draw("ZCOLOR");
    }
  }
}

void reshape_data(const vector<face_data_str> &face_data_v, vector<vector<Double_t>> &data){
  for(unsigned int i = 0; i < face_data_v.size(); i++){
    vector<Double_t> val_v;
    for(unsigned int j = 0; j < face_data_v.at(i).data_v.size(); j++)
      for(unsigned int k = 0; k < face_data_v.at(i).data_v.at(j).size(); k++)
	val_v.push_back(face_data_v.at(i).data_norm_v.at(j).at(k));
    data.push_back(val_v);
  }
}

void get_covariance_data( const vector<vector<Double_t>> &data, vector<vector<Double_t>> &data_sigma, TPrincipal* principal){
  Double_t data_arr[nn_im][dd_im];
  Double_t data_short_arr[nn_im];
  for(unsigned int i = 0; i < nn_im; i++)
    for(unsigned int j = 0; j < dd_im; j++)
      data_arr[i][j] = data.at(i).at(j);
  Double_t data_sigma_arr[nn_im][nn_im];  
  //cout<<"data.size()       "<<data.size()<<endl
  //<<"data.at(0).size() "<<data.at(0).size()<<endl;
  Double_t val = 0.0;
  for(unsigned int i = 0; i < nn_im; i++){
    if(i%100==0)
      cout<<i<<endl;
    for(unsigned int j = 0; j < nn_im; j++){
      val = 0.0;
      for(unsigned int k = 0; k < dd_im; k++)
	val += data_arr[i][k]*data_arr[j][k];  
      data_sigma_arr[i][j] = val;
      data_short_arr[j]=val;
    }
    principal->AddRow(data_short_arr);
  }
  for(unsigned int i = 0; i < nn_im; i++){
    for(unsigned int j = 0; j < nn_im; j++){
      data_sigma.at(i).at(j) = data_sigma_arr[i][j];
    }
  }
  cout<<"data_sigma.size()       "<<data_sigma.size()<<endl
      <<"data_sigma.at(0).size() "<<data_sigma.at(0).size()<<endl;
}

void fill_2Dhist_principal_components(vector<TH2D*> &h2_principal_components_v){
  Int_t l = TMath::Sqrt(dd_im);
  Int_t i_x, i_y;
  for(unsigned int i = 0; i < dd_im; i++){
    TString name = "h2";
    name += "_principal_components";
    name += i;
    TH2D *h2 = new TH2D( name.Data(), name.Data(), l,0,l, l,0,l);
    //
    for(unsigned int j = 0; j < dd_im; j++){
      i_x = j/l;
      i_y = j-l*i_x;
      h2->SetBinContent(l-i_x,l-i_y,data_Vh[i][j]);
    }
    h2_principal_components_v.push_back(h2);
  }
}

void save_to_csv_data(TString file_out_name, const vector<vector<Double_t>> &data){
  //
  std::ofstream csvfile;
  //
  cout<<"data.size()       = "<<data.size()<<endl
      <<"data.at(0).size() = "<<data.at(0).size()<<endl;
  //
  csvfile.open(file_out_name.Data());
  for(unsigned int i = 0; i <data.size() ; i++){
    for(unsigned int j = 0; j < data.at(i).size(); j++){
      if(j == (data.at(i).size()-1))
	csvfile<<data.at(i).at(j);
      else
	csvfile<<data.at(i).at(j)<<" ";
    }
    csvfile<<endl;
  }
  csvfile.close();  
}
