
#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <TChain.h>

class LikelihoodClass {
public:
  LikelihoodClass();
  ~LikelihoodClass();
  double evaluate_likelihood(std::vector<double> inputvars);
private:
  TH1F * h_fstatekinematics_sig[19];
  TH1F * h_fstatekinematics_bkg[19];
  TFile * file_inputdistributions;
};

LikelihoodClass::LikelihoodClass(){
  file_inputdistributions=new TFile("LikelihoodInput_test_ele.root", "READ");
  h_fstatekinematics_sig[0]= (TH1F*) file_inputdistributions->Get("thq_dipho_pt");
  h_fstatekinematics_sig[1]= (TH1F*) file_inputdistributions->Get("thq_dipho_eta");
  h_fstatekinematics_sig[2]= (TH1F*) file_inputdistributions->Get("thq_dipho_leadPt");
  h_fstatekinematics_sig[3]= (TH1F*) file_inputdistributions->Get("thq_dipho_leadEta");
  h_fstatekinematics_sig[4]= (TH1F*) file_inputdistributions->Get("thq_dipho_subleadPt");
  h_fstatekinematics_sig[5]= (TH1F*) file_inputdistributions->Get("thq_dipho_subleadEta");
//  h_fstatekinematics_sig[6]= (TH1F*) file_inputdistributions->Get("thq_muon1_pt");
//  h_fstatekinematics_sig[7]= (TH1F*) file_inputdistributions->Get("thq_muon1_eta");
  h_fstatekinematics_sig[6]= (TH1F*) file_inputdistributions->Get("thq_ele1_pt");
  h_fstatekinematics_sig[7]= (TH1F*) file_inputdistributions->Get("thq_ele1_eta");
//  h_fstatekinematics_sig[8]= (TH1F*) file_inputdistributions->Get("thq_jet1_pt");
//  h_fstatekinematics_sig[9]= (TH1F*) file_inputdistributions->Get("thq_jet1_eta");
  h_fstatekinematics_sig[8]= (TH1F*) file_inputdistributions->Get("thq_bjet1_pt");
  h_fstatekinematics_sig[9]= (TH1F*) file_inputdistributions->Get("thq_bjet1_eta");
  h_fstatekinematics_sig[10]= (TH1F*) file_inputdistributions->Get("thq_MET_pt");
  h_fstatekinematics_sig[11]= (TH1F*) file_inputdistributions->Get("thq_dr_tHchainfwdjet");
  h_fstatekinematics_sig[12]= (TH1F*) file_inputdistributions->Get("thq_dr_bjetfwdjet");
  h_fstatekinematics_sig[13]= (TH1F*) file_inputdistributions->Get("thq_dr_subleadphobjet");
  h_fstatekinematics_sig[14]= (TH1F*) file_inputdistributions->Get("thq_dr_leadphofwdjet");
  h_fstatekinematics_sig[15]= (TH1F*) file_inputdistributions->Get("thq_dr_subleadphofwdjet");
  h_fstatekinematics_sig[16]= (TH1F*) file_inputdistributions->Get("thq_dr_leptonbjet");
  h_fstatekinematics_sig[17]= (TH1F*) file_inputdistributions->Get("thq_fwdjet1_pt");
  h_fstatekinematics_sig[18]= (TH1F*) file_inputdistributions->Get("thq_fwdjet1_eta");

  h_fstatekinematics_bkg[0]= (TH1F*) file_inputdistributions->Get("tth_dipho_pt");
  h_fstatekinematics_bkg[1]= (TH1F*) file_inputdistributions->Get("tth_dipho_eta");
  h_fstatekinematics_bkg[2]= (TH1F*) file_inputdistributions->Get("tth_dipho_leadPt");
  h_fstatekinematics_bkg[3]= (TH1F*) file_inputdistributions->Get("tth_dipho_leadEta");
  h_fstatekinematics_bkg[4]= (TH1F*) file_inputdistributions->Get("tth_dipho_subleadPt");
  h_fstatekinematics_bkg[5]= (TH1F*) file_inputdistributions->Get("tth_dipho_subleadEta");
//  h_fstatekinematics_bkg[6]= (TH1F*) file_inputdistributions->Get("tth_muon1_pt");
//  h_fstatekinematics_bkg[7]= (TH1F*) file_inputdistributions->Get("tth_muon1_eta");
  h_fstatekinematics_bkg[6]= (TH1F*) file_inputdistributions->Get("tth_ele1_pt");
  h_fstatekinematics_bkg[7]= (TH1F*) file_inputdistributions->Get("tth_ele1_eta");
//  h_fstatekinematics_bkg[8]= (TH1F*) file_inputdistributions->Get("tth_jet1_pt");
//  h_fstatekinematics_bkg[9]= (TH1F*) file_inputdistributions->Get("tth_jet1_eta");
  h_fstatekinematics_bkg[8]= (TH1F*) file_inputdistributions->Get("tth_bjet1_pt");
  h_fstatekinematics_bkg[9]= (TH1F*) file_inputdistributions->Get("tth_bjet1_eta");
  h_fstatekinematics_bkg[10]= (TH1F*) file_inputdistributions->Get("tth_MET_pt");
  h_fstatekinematics_bkg[11]= (TH1F*) file_inputdistributions->Get("tth_dr_tHchainfwdjet");
  h_fstatekinematics_bkg[12]= (TH1F*) file_inputdistributions->Get("tth_dr_bjetfwdjet");
  h_fstatekinematics_bkg[13]= (TH1F*) file_inputdistributions->Get("tth_dr_subleadphobjet");
  h_fstatekinematics_bkg[14]= (TH1F*) file_inputdistributions->Get("tth_dr_leadphofwdjet");
  h_fstatekinematics_bkg[15]= (TH1F*) file_inputdistributions->Get("tth_dr_subleadphofwdjet");
  h_fstatekinematics_bkg[16]= (TH1F*) file_inputdistributions->Get("tth_dr_leptonbjet");
  h_fstatekinematics_bkg[17]= (TH1F*) file_inputdistributions->Get("tth_fwdjet1_pt");
  h_fstatekinematics_bkg[18]= (TH1F*) file_inputdistributions->Get("tth_fwdjet1_eta");

}
 LikelihoodClass::~LikelihoodClass(){
  file_inputdistributions->Close();
}
double LikelihoodClass::evaluate_likelihood(std::vector<double> inputvars){
  if(inputvars.size()!=19){
    std::cout<<"<LikelihoodClass::evaluate_likelihood>: inputvars.size() is "<<inputvars.size()<< " and it is expected to be 8!"<<std::endl;
    std::cout<<"<LikelihoodClass::evaluate_likelihood>: is returning a value of -10.!"<<std::endl;
  return -10.;
  }
  double p_signal=1.;
  double p_background=1.;
  for(unsigned int ielement=0; ielement<inputvars.size(); ielement++){
    int bin_hsignal= (Int_t) h_fstatekinematics_sig[ielement]->FindBin(inputvars.at(ielement));
    int bin_hbackground= (Int_t) h_fstatekinematics_bkg[ielement]->FindBin(inputvars[ielement]);
    p_signal = (p_signal*h_fstatekinematics_sig[ielement]->GetBinContent(bin_hsignal));
    p_background = (p_background*h_fstatekinematics_bkg[ielement]->GetBinContent(bin_hbackground));
//    cout<<"inputvars.at(ielement)"<<inputvars.at(ielement)<<endl;
//    cout<<"bin_hsignal=  "<<bin_hsignal<<endl;; 
//  cout<<"p_signal"<<p_signal<<endl;
  }
//  cout<<"bin_hsignal=  "<<bin_hsignal<<endl;
  double lhood_val_den = (p_signal+p_background);
  return (lhood_val_den != 0. ? (p_signal/lhood_val_den) : -11.);

}


void comparison_tth_thq_test_ele()
   {
   LikelihoodClass *lhoodclass=new LikelihoodClass();
   TChain *tree_thq = new TChain("tagsDumper/trees/THQ_13TeV_THQLeptonicTag");
   tree_thq->Add("output_THQ_0.root");
   tree_thq->Add("output_THQ_1.root");
   tree_thq->Add("output_THQ_2.root");
   tree_thq->Add("output_THQ_3.root");
   tree_thq->Add("output_THQ_4.root");
   tree_thq->Add("output_THQ_5.root");
   tree_thq->Add("output_THQ_6.root");
   tree_thq->Add("output_THQ_7.root");
   tree_thq->Add("output_THQ_8.root");
   tree_thq->Add("output_THQ_9.root");
   tree_thq->Add("output_THQ_10.root");
   tree_thq->Add("output_THQ_11.root");
   tree_thq->Add("output_THQ_12.root");
   tree_thq->Add("output_THQ_13.root");
   tree_thq->Add("output_THQ_14.root");
   tree_thq->Add("output_THQ_15.root");
   tree_thq->Add("output_THQ_16.root");
   tree_thq->Add("output_THQ_17.root");
   tree_thq->Add("output_THQ_18.root");
   tree_thq->Add("output_THQ_19.root");
   tree_thq->Add("output_THQ_20.root");
   tree_thq->Add("output_THQ_21.root");
   tree_thq->Add("output_THQ_22.root");
   tree_thq->Add("output_THQ_23.root");
   tree_thq->Add("output_THQ_24.root");
   tree_thq->Add("output_THQ_25.root");
   tree_thq->Add("output_THQ_26.root");
   tree_thq->Add("output_THQ_27.root");
   tree_thq->Add("output_THQ_28.root");
   tree_thq->Add("output_THQ_29.root");
   tree_thq->Add("output_THQ_30.root");
   tree_thq->Add("output_THQ_31.root");
   tree_thq->Add("output_THQ_32.root");
   tree_thq->Add("output_THQ_33.root");
   tree_thq->Add("output_THQ_34.root");
   tree_thq->Add("output_THQ_35.root");
   tree_thq->Add("output_THQ_36.root");
   tree_thq->Add("output_THQ_37.root");

   TChain *tree_tth = new TChain("tagsDumper/trees/TTHJet_13TeV_THQLeptonicTag");
   tree_tth->Add("output_TTHJet_0.root");
   tree_tth->Add("output_TTHJet_1.root");
   tree_tth->Add("output_TTHJet_2.root");
   tree_tth->Add("output_TTHJet_3.root");
   tree_tth->Add("output_TTHJet_4.root");
   tree_tth->Add("output_TTHJet_5.root");
   tree_tth->Add("output_TTHJet_6.root");
   tree_tth->Add("output_TTHJet_7.root");
   tree_tth->Add("output_TTHJet_8.root");
   tree_tth->Add("output_TTHJet_9.root");
/*   tree_tth->Add("output_TTH_10.root");
   tree_tth->Add("output_TTH_11.root");
   tree_tth->Add("output_TTH_12.root");
   tree_tth->Add("output_TTH_13.root");
   tree_tth->Add("output_TTH_14.root");
   tree_tth->Add("output_TTH_15.root");
   tree_tth->Add("output_TTH_16.root");
   tree_tth->Add("output_TTH_17.root");
   tree_tth->Add("output_TTH_18.root");
   tree_tth->Add("output_TTH_19.root");
*/
   TFile *file_output = new TFile("thq_tth_test_ele.root","RECREATE");


   Float_t qdipho_pt, qdipho_eta, qdipho_leadPt, qdipho_leadEta, qdipho_subleadPt, qdipho_subleadEta;
   Float_t qele1_pt, qele1_eta, qele2_pt, qele2_eta;
   Float_t qmuon1_pt, qmuon2_pt, qmuon1_eta, qmuon2_eta;
   Float_t qjet1_pt, qjet1_eta, qfwdjet1_pt, qfwdjet1_eta, qbjet1_pt, qbjet1_eta, qMET_pt;
   Float_t qdr_tHchainfwdjet, qdr_bjetfwdjet, qdr_subleadphobjet, qdr_leadphofwdjet, qdr_subleadphofwdjet, qdr_leptonbjet;
   Float_t qLeptonType;

//   Float_t rele1_pt, rele1_eta, rele2_pt, rele2_eta;
   Float_t tdipho_pt, tdipho_eta, tdipho_leadPt, tdipho_leadEta, tdipho_subleadPt, tdipho_subleadEta;
   Float_t tele1_pt, tele1_eta, tele2_pt, tele2_eta;
   Float_t tmuon1_pt, tmuon2_pt, tmuon1_eta, tmuon2_eta;
   Float_t tjet1_pt, tjet1_eta, tfwdjet1_pt, tfwdjet1_eta, tbjet1_pt, tbjet1_eta, tMET_pt;
   Float_t tdr_tHchainfwdjet, tdr_bjetfwdjet, tdr_subleadphobjet, tdr_leadphofwdjet, tdr_subleadphofwdjet, tdr_leptonbjet;
   Float_t tLeptonType;
   Float_t tn_jets, tn_fwdjets, tn_bjets, tn_ele, tn_muons; 
   Float_t qn_jets, qn_fwdjets, qn_bjets, qn_ele, qn_muons;

   tree_thq->SetBranchAddress("dipho_pt", &qdipho_pt);
   tree_thq->SetBranchAddress("dipho_eta", &qdipho_eta);
   tree_thq->SetBranchAddress("dipho_leadPt", &qdipho_leadPt);
   tree_thq->SetBranchAddress("dipho_leadEta", &qdipho_leadEta);
   tree_thq->SetBranchAddress("dipho_subleadPt", &qdipho_subleadPt);
   tree_thq->SetBranchAddress("dipho_subleadEta", &qdipho_subleadEta);
   tree_thq->SetBranchAddress("n_muons", &qn_muons);
   tree_thq->SetBranchAddress("muon1_pt", &qmuon1_pt);
   tree_thq->SetBranchAddress("muon1_eta", &qmuon1_eta);
   tree_thq->SetBranchAddress("n_ele", &qn_ele);
   tree_thq->SetBranchAddress("ele1_pt", &qele1_pt);
   tree_thq->SetBranchAddress("ele1_eta", &qele1_eta);
   tree_thq->SetBranchAddress("ele2_pt", &qele2_pt);
   tree_thq->SetBranchAddress("ele2_eta", &qele2_eta);
   tree_thq->SetBranchAddress("n_jets", &qn_jets);
   tree_thq->SetBranchAddress("jet1_pt", &qjet1_pt);
   tree_thq->SetBranchAddress("jet1_eta", &qjet1_eta);
   tree_thq->SetBranchAddress("n_fwdjets", &qn_fwdjets);
   tree_thq->SetBranchAddress("fwdjet1_pt", &qfwdjet1_pt);
   tree_thq->SetBranchAddress("fwdjet1_eta", &qfwdjet1_eta);
   tree_thq->SetBranchAddress("n_bjets", &qn_bjets);
   tree_thq->SetBranchAddress("bjet1_pt", &qbjet1_pt);
   tree_thq->SetBranchAddress("bjet1_eta", &qbjet1_eta);
   tree_thq->SetBranchAddress("dr_tHchainfwdjet", &qdr_tHchainfwdjet);
   tree_thq->SetBranchAddress("dr_bjetfwdjet", &qdr_bjetfwdjet);
   tree_thq->SetBranchAddress("dr_subleadphobjet", &qdr_subleadphobjet);
   tree_thq->SetBranchAddress("dr_leadphofwdjet", &qdr_leadphofwdjet);
   tree_thq->SetBranchAddress("dr_subleadphofwdjet", &qdr_subleadphofwdjet);
   tree_thq->SetBranchAddress("dr_leptonbjet", &qdr_leptonbjet);
   tree_thq->SetBranchAddress("LeptonType", &qLeptonType);
   tree_thq->SetBranchAddress("MET_pt", &qMET_pt);

//   tree_thq->SetBranchStatus("*", 0);
   tree_thq->SetBranchStatus("dipho_pt", 1);
   tree_thq->SetBranchStatus("dipho_eta", 1);
   tree_thq->SetBranchStatus("dipho_leadPt", 1);
   tree_thq->SetBranchStatus("dipho_leadEta", 1);
   tree_thq->SetBranchStatus("dipho_subleadPt", 1);
   tree_thq->SetBranchStatus("dipho_subleadEta", 1);
   tree_thq->SetBranchStatus("n_muons", 1);
   tree_thq->SetBranchStatus("muon1_pt", 1);
   tree_thq->SetBranchStatus("muon1_eta", 1);
   tree_thq->SetBranchStatus("n_ele", 1);
   tree_thq->SetBranchStatus("ele1_pt", 1);
   tree_thq->SetBranchStatus("ele1_eta", 1);
   tree_thq->SetBranchStatus("n_jets", 1);   
   tree_thq->SetBranchStatus("jet1_pt", 1);
   tree_thq->SetBranchStatus("jet1_eta", 1);
   tree_thq->SetBranchStatus("n_fwdjets", 1);
   tree_thq->SetBranchStatus("fwdjet1_pt", 1);
   tree_thq->SetBranchStatus("fwdjet1_eta", 1);
   tree_thq->SetBranchStatus("n_bjets", 1);
   tree_thq->SetBranchStatus("bjet1_pt", 1);
   tree_thq->SetBranchStatus("bjet1_eta", 1);
   tree_thq->SetBranchStatus("dr_tHchainfwdjet", 1);
   tree_thq->SetBranchStatus("dr_bjetfwdjet", 1);
   tree_thq->SetBranchStatus("dr_subleadphobjet", 1);
   tree_thq->SetBranchStatus("dr_leadphofwdjet", 1);
   tree_thq->SetBranchStatus("dr_subleadphofwdjet", 1);
   tree_thq->SetBranchStatus("dr_leptonbjet", 1);
   tree_thq->SetBranchStatus("LeptonType", 1);
   tree_thq->SetBranchStatus("MET_pt", 1); 

   tree_tth->SetBranchAddress("dipho_pt", &tdipho_pt);
   tree_tth->SetBranchAddress("dipho_eta", &tdipho_eta);
   tree_tth->SetBranchAddress("dipho_leadPt", &tdipho_leadPt);
   tree_tth->SetBranchAddress("dipho_leadEta", &tdipho_leadEta);
   tree_tth->SetBranchAddress("dipho_subleadPt", &tdipho_subleadPt);
   tree_tth->SetBranchAddress("dipho_subleadEta", &tdipho_subleadEta);
   tree_tth->SetBranchAddress("n_muons", &tn_muons);
   tree_tth->SetBranchAddress("muon1_pt", &tmuon1_pt);
   tree_tth->SetBranchAddress("muon1_eta", &tmuon1_eta);
   tree_tth->SetBranchAddress("n_ele", &tn_ele);
   tree_tth->SetBranchAddress("ele1_pt", &tele1_pt);
   tree_tth->SetBranchAddress("ele1_eta", &tele1_eta);
   tree_tth->SetBranchAddress("ele2_pt", &tele2_pt);
   tree_tth->SetBranchAddress("ele2_eta", &tele2_eta);
   tree_tth->SetBranchAddress("n_jets", &tn_jets);
   tree_tth->SetBranchAddress("jet1_pt", &tjet1_pt);
   tree_tth->SetBranchAddress("jet1_eta", &tjet1_eta);
   tree_tth->SetBranchAddress("n_fwdjets", &tn_fwdjets);
   tree_tth->SetBranchAddress("fwdjet1_pt", &tfwdjet1_pt);
   tree_tth->SetBranchAddress("fwdjet1_eta", &tfwdjet1_eta);
   tree_tth->SetBranchAddress("n_bjets", &tn_bjets);
   tree_tth->SetBranchAddress("bjet1_pt", &tbjet1_pt);
   tree_tth->SetBranchAddress("bjet1_eta", &tbjet1_eta);
   tree_tth->SetBranchAddress("dr_tHchainfwdjet", &tdr_tHchainfwdjet);
   tree_tth->SetBranchAddress("dr_bjetfwdjet", &tdr_bjetfwdjet);
   tree_tth->SetBranchAddress("dr_subleadphobjet", &tdr_subleadphobjet);
   tree_tth->SetBranchAddress("dr_leadphofwdjet", &tdr_leadphofwdjet);
   tree_tth->SetBranchAddress("dr_subleadphofwdjet", &tdr_subleadphofwdjet);
   tree_tth->SetBranchAddress("dr_leptonbjet", &tdr_leptonbjet);
   tree_tth->SetBranchAddress("LeptonType", &tLeptonType);
   tree_tth->SetBranchAddress("MET_pt", &tMET_pt);

//   tree_tth->SetBranchStatus("*", 0);
   tree_tth->SetBranchStatus("dipho_pt", 1);
   tree_tth->SetBranchStatus("dipho_eta", 1);
   tree_tth->SetBranchStatus("dipho_leadPt", 1);
   tree_tth->SetBranchStatus("dipho_leadEta", 1);
   tree_tth->SetBranchStatus("dipho_subleadPt", 1);
   tree_tth->SetBranchStatus("dipho_subleadEta", 1);
   tree_tth->SetBranchStatus("n_muons", 1);
   tree_tth->SetBranchStatus("muon1_pt", 1);
   tree_tth->SetBranchStatus("muon1_eta", 1);
   tree_tth->SetBranchStatus("n_ele",1);
   tree_tth->SetBranchStatus("ele1_pt", 1);
   tree_tth->SetBranchStatus("ele1_eta", 1);
   tree_tth->SetBranchStatus("n_jets",1);
   tree_tth->SetBranchStatus("jet1_pt", 1);
   tree_tth->SetBranchStatus("jet1_eta", 1);
   tree_tth->SetBranchStatus("n_fwdjets",1);
   tree_tth->SetBranchStatus("fwdjet1_pt", 1);
   tree_tth->SetBranchStatus("fwdjet1_eta", 1);
   tree_tth->SetBranchStatus("n_bjets", 1);
   tree_tth->SetBranchStatus("bjet1_pt", 1);
   tree_tth->SetBranchStatus("bjet1_eta", 1);
   tree_tth->SetBranchStatus("dr_tHchainfwdjet", 1);
   tree_tth->SetBranchStatus("dr_bjetfwdjet", 1);
   tree_tth->SetBranchStatus("dr_subleadphobjet", 1);
   tree_tth->SetBranchStatus("dr_leadphofwdjet", 1);
   tree_tth->SetBranchStatus("dr_subleadphofwdjet", 1);
   tree_tth->SetBranchStatus("dr_leptonbjet", 1);
   tree_tth->SetBranchStatus("LeptonType", 1);
   tree_tth->SetBranchStatus("MET_pt", 1);


   TH1F *hdipho_pt= new TH1F("thq_dipho_pt","dipho P_{T} distribution",50,0,600);       //h->thq, h1->tth
   TH1F *hdipho_eta= new TH1F("thq_dipho_eta","dipho #eta distribution",50,-4,4);
   TH1F *hdipho_leadPt= new TH1F("thq_dipho_leadPt","dipho_lead P_{T} distribution",50,0,600);
   TH1F *hdipho_leadEta= new TH1F("thq_dipho_leadEta","dipho_lead #eta distribution",50,-3,3);
   TH1F *hdipho_subleadPt= new TH1F("thq_dipho_subleadPt","dipho_sublead P_{T} distribution",50,0,300);
   TH1F *hdipho_subleadEta= new TH1F("thq_dipho_subleadEta","dipho_sublead #eta distribution",50,-3,3);
   TH1F *hmuon1_pt= new TH1F("thq_muon1_pt","muon1 P_{T} distribution",50,0,300);
   TH1F *hmuon1_eta= new TH1F("thq_muon1_eta","muon1 #eta distribution",50,-3,3);
   TH1F *hele1_pt= new TH1F("thq_ele1_pt","ele1 P_{T} distribution",50,0,400);
   TH1F *hele1_eta= new TH1F("thq_ele1_eta","ele1 #eta distribution",50,-3,3);
   TH1F *hjet1_pt= new TH1F("thq_jet1_pt","jet1 P_{T} distribution",50,0,600);
   TH1F *hjet1_eta= new TH1F("thq_jet1_eta","jet1 #eta distribution",50,-3,3);
   TH1F *hfwdjet1_pt= new TH1F("thq_fwdjet1_pt","fwdjet1 P_{T} distribution",50,0,600);
   TH1F *hfwdjet1_eta= new TH1F("thq_fwdjet1_eta","fwdjet1 #eta distribution",50,-3,3);
   TH1F *hbjet1_pt= new TH1F("thq_bjet1_pt","bjet1 P_{T} distribution",50,0,600);
   TH1F *hbjet1_eta= new TH1F("thq_bjet1_eta","bjet1 #eta distribution",50,-3,3);
   TH1F *hMET_pt= new TH1F("thq_MET_pt", "MET P_{T} distribution",50,0,600);
   TH1F *hdr_tHchainfwdjet= new TH1F("thq_dr_tHchainfwdjet"," #Delta R tHchain and fwdjet distribution",50,0,7);
   TH1F *hdr_bjetfwdjet= new TH1F("thq_dr_bjetfwdjet","#Delta R bjet and fwdjet distribution",50,0,7);
   TH1F *hdr_subleadphobjet= new TH1F("thq_dr_subleadphobjet","#Delta R subleadpho bjet distribution",50,0,7);
   TH1F *hdr_leadphofwdjet= new TH1F("thq_dr_leadphofwdjet","#Delta R leadpho and fwdjet distribution",50,0,7);
   TH1F *hdr_subleadphofwdjet= new TH1F("thq_dr_subleadphofwdjet","#Delta R subleadpho and fwdjet distribution",50,0,7);
   TH1F *hdr_leptonbjet= new TH1F("thq_dr_leptonbjet","#Delta R lepton and bjet distribution",50,0,7);

   TH1F *h1dipho_pt= new TH1F("tth_dipho_pt","dipho P_{T} distribution",50,0,600);       //h->thq, h1->tth
   TH1F *h1dipho_eta= new TH1F("tth_dipho_eta","dipho #eta distribution",50,-4,4);
   TH1F *h1dipho_leadPt= new TH1F("tth_dipho_leadPt","dipho_lead P_{T} distribution",50,0,600);
   TH1F *h1dipho_leadEta= new TH1F("tth_dipho_leadEta","dipho_lead #eta distribution",50,-3,3);
   TH1F *h1dipho_subleadPt= new TH1F("tth_dipho_subleadPt","dipho_sublead P_{T} distribution",50,0,300);
   TH1F *h1dipho_subleadEta= new TH1F("tth_dipho_subleadEta","dipho_sublead #eta distribution",50,-3,3);
   TH1F *h1muon1_pt= new TH1F("tth_muon1_pt","muon1 P_{T} distribution",50,0,300);
   TH1F *h1muon1_eta= new TH1F("tth_muon1_eta","muon1 #eta distribution",50,-3,3);
   TH1F *h1ele1_pt= new TH1F("tth_ele1_pt","ele1 P_{T} distribution",50,0,400);
   TH1F *h1ele1_eta= new TH1F("tth_ele1_eta","ele1 #eta distribution",50,-3,3);
   TH1F *h1jet1_pt= new TH1F("tth_jet1_pt","jet1 P_{T} distribution",50,0,600);
   TH1F *h1jet1_eta= new TH1F("tth_jet1_eta","jet1 #eta distribution",50,-3,3);
   TH1F *h1fwdjet1_pt= new TH1F("tth_fwdjet1_pt","fwdjet1 P_{T} distribution",50,0,600);
   TH1F *h1fwdjet1_eta= new TH1F("tth_fwdjet1_eta","fwdjet1 #eta distribution",50,-3,3);
   TH1F *h1bjet1_pt= new TH1F("tth_bjet1_pt","bjet1 P_{T} distribution",50,0,600);
   TH1F *h1bjet1_eta= new TH1F("tth_bjet1_eta","bjet1 #eta distribution",50,-3,3);
   TH1F *h1MET_pt= new TH1F("tth_MET_pt", "MET P_{T} distribution",50,0,600);
   TH1F *h1dr_tHchainfwdjet= new TH1F("tth_dr_tHchainfwdjet","#Delta R tHchain and fwdjet distribution",50,0,7);
   TH1F *h1dr_bjetfwdjet= new TH1F("tth_dr_bjetfwdjet","#Delta R bjet and fwdjet distribution",50,0,7);
   TH1F *h1dr_subleadphobjet= new TH1F("tth_dr_subleadphobjet","#Delta R subleadpho and bjet distribution",50,0,7);
   TH1F *h1dr_leadphofwdjet= new TH1F("tth_dr_leadphofwdjet","#Delta R leadpho and fwdjet distribution",50,0,7);
   TH1F *h1dr_subleadphofwdjet= new TH1F("tth_dr_subleadphofwdjet","#Delta R subleadpho and fwdjet distribution",50,0,7);
   TH1F *h1dr_leptonbjet= new TH1F("tth_dr_leptonbjet","#Delta R lepton and bjet distribution",50,0,7);

   TH1F *h_likelihood_thq=new TH1F("likelihood_thq", "Likelihood (e+jets) ", 50, 0., 1.);
   TH1F *h_likelihood_tth=new TH1F("likelihood_tth", "Likelihood (e+jets) ", 50, 0., 1.);   

   Int_t nentries_thq=(Int_t)tree_thq->GetEntries();
   Int_t nentries_tth=(Int_t)tree_tth->GetEntries();
   
   int tth_entries=0;
  int thq_entries=0;
  int tth_sel_evnt=0;
  int thq_sel_evnt=0;
  int tthevents=0;
  int thqevents=0;
   
//thq loop
   for(int ievent=0; ievent<nentries_thq; ievent++)
     {
     tree_thq -> GetEntry( ievent );
     thq_entries++;     
     if( /* qn_jets == 0 && qn_fwdjets == 1 && qn_bjets == 1 &&*/ qn_ele == 1 && qn_muons == 0){     
//     cout<<"qn_bjets_thq="<<qn_bjets<<endl;
//     cout<<"qn_fwdjets_thq="<<qn_fwdjets<<endl;

     hdipho_pt->Fill(qdipho_pt);
     hdipho_eta->Fill(qdipho_eta);
     hdipho_leadPt->Fill(qdipho_leadPt);
     hdipho_leadEta->Fill(qdipho_leadEta);
     hdipho_subleadPt->Fill(qdipho_subleadPt);
     hdipho_subleadEta->Fill(qdipho_subleadEta);
     hmuon1_eta->Fill(qmuon1_eta);
     hmuon1_pt->Fill(qmuon1_pt);
     hele1_eta->Fill(qele1_eta);
     hele1_pt->Fill(qele1_pt);
     hjet1_eta->Fill(qjet1_eta);
     hjet1_pt->Fill(qjet1_pt);
     hfwdjet1_eta->Fill(qfwdjet1_eta);
     hfwdjet1_pt->Fill(qfwdjet1_pt);
     hbjet1_eta->Fill(qbjet1_eta);
     hbjet1_pt->Fill(qbjet1_pt);
     hMET_pt->Fill(qMET_pt);
     hdr_tHchainfwdjet -> Fill( qdr_tHchainfwdjet );
     hdr_bjetfwdjet->Fill( qdr_bjetfwdjet );
     hdr_subleadphobjet->Fill( qdr_subleadphobjet );
     hdr_leadphofwdjet->Fill( qdr_leadphofwdjet );
     hdr_subleadphofwdjet->Fill( qdr_subleadphofwdjet );
     hdr_leptonbjet->Fill( qdr_leptonbjet );
          
     std::vector<double> vec_lhood_calc;
     vec_lhood_calc.push_back(qdipho_pt);
     vec_lhood_calc.push_back(qdipho_eta);
     vec_lhood_calc.push_back(qdipho_leadPt);
     vec_lhood_calc.push_back(qdipho_leadEta);
     vec_lhood_calc.push_back(qdipho_subleadPt);
     vec_lhood_calc.push_back(qdipho_subleadEta);
//     vec_lhood_calc.push_back(qmuon1_pt);
//     vec_lhood_calc.push_back(qmuon1_eta);
     vec_lhood_calc.push_back(qele1_pt);
     vec_lhood_calc.push_back(qele1_eta);
//     vec_lhood_calc.push_back(qjet1_pt);
//     vec_lhood_calc.push_back(qjet1_eta);
     vec_lhood_calc.push_back(qbjet1_pt);
     vec_lhood_calc.push_back(qbjet1_eta);
     vec_lhood_calc.push_back(qMET_pt);
     vec_lhood_calc.push_back(qdr_tHchainfwdjet);
     vec_lhood_calc.push_back(qdr_bjetfwdjet );
     vec_lhood_calc.push_back(qdr_subleadphobjet );
     vec_lhood_calc.push_back(qdr_leadphofwdjet );
     vec_lhood_calc.push_back(qdr_subleadphofwdjet );
     vec_lhood_calc.push_back(qdr_leptonbjet );
     vec_lhood_calc.push_back(qfwdjet1_pt);
     vec_lhood_calc.push_back(qfwdjet1_eta);     

     thq_sel_evnt++;     

     double lhood_value=lhoodclass->evaluate_likelihood(vec_lhood_calc);
//     cout<<"lhood_value"<<lhood_value<<endl;
     h_likelihood_thq->Fill(lhood_value);
     if(lhood_value>0.5){
     thqevents++;
     }

     }//end of thq loop
     }
//tth loop
   for(int ievent=0; ievent<nentries_tth; ievent++)
     {
     tree_tth->GetEntry(ievent);
     tth_entries++;
     if(/*tn_jets == 0 && tn_fwdjets == 1 && tn_bjets ==1 &&*/ tn_ele == 1 && tn_muons == 0 )
     {
//     cout<<"tn_bjets_tth="<<tn_bjets<<endl;
//     cout<<"tn_fwdjets_tth="<<tn_fwdjets<<endl;
     h1dipho_pt->Fill(tdipho_pt);
     h1dipho_eta->Fill(tdipho_eta);
     h1dipho_leadPt->Fill(tdipho_leadPt);
     h1dipho_leadEta->Fill(tdipho_leadEta);
     h1dipho_subleadPt->Fill(tdipho_subleadPt);
     h1dipho_subleadEta->Fill(tdipho_subleadEta);
     h1muon1_eta->Fill(tmuon1_eta);                                     //leptonid=1 (electron), 2 (muon)
     h1muon1_pt->Fill(tmuon1_pt);
     h1ele1_eta->Fill(tele1_eta);
     h1ele1_pt->Fill(tele1_pt);
     h1jet1_eta->Fill(tjet1_eta);
     h1jet1_pt->Fill(tjet1_pt);
     h1fwdjet1_eta->Fill(tfwdjet1_eta);
     h1fwdjet1_pt->Fill(tfwdjet1_pt);
     h1bjet1_eta->Fill(tbjet1_eta);
     h1bjet1_pt->Fill(tbjet1_pt);
     h1MET_pt->Fill(tMET_pt);
     h1dr_tHchainfwdjet -> Fill( tdr_tHchainfwdjet );
     h1dr_bjetfwdjet->Fill( tdr_bjetfwdjet );
     h1dr_subleadphobjet->Fill( tdr_subleadphobjet );
     h1dr_leadphofwdjet->Fill( tdr_leadphofwdjet );
     h1dr_subleadphofwdjet->Fill( tdr_subleadphofwdjet );
     h1dr_leptonbjet->Fill( tdr_leptonbjet );
     
     std::vector<double> vec_lhood_calc;
     vec_lhood_calc.push_back(tdipho_pt);
     vec_lhood_calc.push_back(tdipho_eta);
     vec_lhood_calc.push_back(tdipho_leadPt);
     vec_lhood_calc.push_back(tdipho_leadEta);
     vec_lhood_calc.push_back(tdipho_subleadPt);
     vec_lhood_calc.push_back(tdipho_subleadEta);
//     vec_lhood_calc.push_back(tmuon1_pt);
//     vec_lhood_calc.push_back(tmuon1_eta);
     vec_lhood_calc.push_back(tele1_pt);
     vec_lhood_calc.push_back(tele1_eta);
//     vec_lhood_calc.push_back(tjet1_pt);
//     vec_lhood_calc.push_back(tjet1_eta);
     vec_lhood_calc.push_back(tbjet1_pt);
     vec_lhood_calc.push_back(tbjet1_eta);
     vec_lhood_calc.push_back(tMET_pt);
     vec_lhood_calc.push_back(tdr_tHchainfwdjet);
     vec_lhood_calc.push_back(tdr_bjetfwdjet );
     vec_lhood_calc.push_back(tdr_subleadphobjet );
     vec_lhood_calc.push_back(tdr_leadphofwdjet );
     vec_lhood_calc.push_back(tdr_subleadphofwdjet );
     vec_lhood_calc.push_back(tdr_leptonbjet );
     vec_lhood_calc.push_back(tfwdjet1_pt);
     vec_lhood_calc.push_back(tfwdjet1_eta);

     tth_sel_evnt++;


//     cout<<"vec_lhood_calc.size()"<<vec_lhood_calc.size()<<endl;
     double lhood_value=lhoodclass->evaluate_likelihood(vec_lhood_calc);
//     cout<<"lhood_value=   "<<lhood_value<<endl;
     h_likelihood_tth->Fill(lhood_value);
     if(lhood_value>0.5){
     tthevents++;
     }
     }
     }//end of tth event loop
  cout<<"tth_entries"<<tth_entries<<endl;
  cout<<"thq_entries"<<thq_entries<<endl;
  cout<<"tth_sel_evnt"<<tth_sel_evnt<<endl;
  cout<<"thq_sel_evnt"<<thq_sel_evnt<<endl;
  cout<<"tthevents"<<tthevents<<endl;
  cout<<"thqevents"<<thqevents<<endl;

  hdr_tHchainfwdjet->Scale(1/hdr_tHchainfwdjet->Integral(0, 51));
  hdr_bjetfwdjet->Scale(1/hdr_bjetfwdjet->Integral(0, 51));
  hdr_subleadphobjet->Scale(1/hdr_subleadphobjet->Integral(0, 51));
  hdr_leadphofwdjet->Scale(1/hdr_leadphofwdjet->Integral(0, 51));
  hdr_subleadphofwdjet->Scale(1/hdr_subleadphofwdjet->Integral(0, 51));
  hdr_leptonbjet->Scale(1/hdr_leptonbjet->Integral(0, 51));

  hdipho_eta->Scale(1/hdipho_eta->Integral(0, 51));
  hdipho_leadEta->Scale(1/hdipho_leadEta->Integral(0, 51));
  hdipho_subleadEta->Scale(1/hdipho_subleadEta->Integral(0, 51));
  hmuon1_eta->Scale(1/hmuon1_eta->Integral(0, 51));
  hele1_eta->Scale(1/hele1_eta->Integral(0, 51));
  hfwdjet1_eta->Scale(1/hfwdjet1_eta->Integral(0, 51));
  hbjet1_eta->Scale(1/hbjet1_eta->Integral(0, 51));
 
  h1dipho_eta->Scale(1/h1dipho_eta->Integral(0, 51));
  h1dipho_leadEta->Scale(1/h1dipho_leadEta->Integral(0, 51));
  h1dipho_subleadEta->Scale(1/h1dipho_subleadEta->Integral(0, 51));
  h1muon1_eta->Scale(1/h1muon1_eta->Integral(0, 51));
  h1ele1_eta->Scale(1/h1ele1_eta->Integral(0, 51));
  h1fwdjet1_eta->Scale(1/h1fwdjet1_eta->Integral(0, 51));
  h1bjet1_eta->Scale(1/h1bjet1_eta->Integral(0, 51));

  hdipho_pt->Scale(1/hdipho_pt->Integral(0, 51));
  hdipho_leadPt->Scale(1/hdipho_leadPt->Integral(0, 51));
  hdipho_subleadPt->Scale(1/hdipho_subleadPt->Integral(0, 51));
  hmuon1_pt->Scale(1/hmuon1_pt->Integral(0, 51));
  hele1_pt->Scale(1/hele1_pt->Integral(0, 51));
  hfwdjet1_pt->Scale(1/hfwdjet1_pt->Integral(0, 51));
  hbjet1_pt->Scale(1/hbjet1_pt->Integral(0, 51));
  hMET_pt->Scale(1/hMET_pt->Integral(0, 51));

  h1dipho_pt->Scale(1/h1dipho_pt->Integral(0, 51));
  h1dipho_leadPt->Scale(1/h1dipho_leadPt->Integral(0, 51));
  h1dipho_subleadPt->Scale(1/h1dipho_subleadPt->Integral(0, 51));
  h1muon1_pt->Scale(1/h1muon1_pt->Integral(0, 51));
  h1ele1_pt->Scale(1/h1ele1_pt->Integral(0, 51));
  h1fwdjet1_pt->Scale(1/h1fwdjet1_pt->Integral(0, 51));
  h1bjet1_pt->Scale(1/h1bjet1_pt->Integral(0, 51));
  h1MET_pt->Scale(1/h1MET_pt->Integral(0, 51));

  h1dr_tHchainfwdjet->Scale(1/h1dr_tHchainfwdjet->Integral(0, 51));
  h1dr_bjetfwdjet->Scale(1/h1dr_bjetfwdjet->Integral(0, 51));
  h1dr_subleadphobjet->Scale(1/h1dr_subleadphobjet->Integral(0, 51));
  h1dr_leadphofwdjet->Scale(1/h1dr_leadphofwdjet->Integral(0, 51));
  h1dr_subleadphofwdjet->Scale(1/h1dr_subleadphofwdjet->Integral(0, 51));
  h1dr_leptonbjet->Scale(1/h1dr_leptonbjet->Integral(0, 51));
  
  h_likelihood_tth->Scale(1/h_likelihood_tth->Integral(0, 51));
  h_likelihood_thq->Scale(1/h_likelihood_thq->Integral(0, 51)); 

  hdipho_eta->SetLineColor(4);
  hdipho_eta->SetMarkerColor(4);
  hdipho_eta->SetLineWidth(2);

  hdipho_leadEta->SetLineColor(4);
  hdipho_leadEta->SetMarkerColor(4);
  hdipho_leadEta->SetLineWidth(2);

  hdipho_subleadEta->SetLineColor(4);
  hdipho_subleadEta->SetMarkerColor(4);
  hdipho_subleadEta->SetLineWidth(2);

  hmuon1_eta->SetLineColor(4);
  hmuon1_eta->SetMarkerColor(4);
  hmuon1_eta->SetLineWidth(2);

  hele1_eta->SetLineColor(4);
  hele1_eta->SetMarkerColor(4);
  hele1_eta->SetLineWidth(2);
 
  hfwdjet1_eta->SetLineColor(4);
  hfwdjet1_eta->SetMarkerColor(4);
  hfwdjet1_eta->SetLineWidth(2);

  hbjet1_eta->SetLineColor(4);
  hbjet1_eta->SetMarkerColor(4);
  hbjet1_eta->SetLineWidth(2);

  hdipho_pt->SetLineColor(4);
  hdipho_pt->SetMarkerColor(4);
  hdipho_pt->SetLineWidth(2);

  hdipho_leadPt->SetLineColor(4);
  hdipho_leadPt->SetMarkerColor(4);
  hdipho_leadPt->SetLineWidth(2);

  hdipho_subleadPt->SetLineColor(4);
  hdipho_subleadPt->SetMarkerColor(4);
  hdipho_subleadPt->SetLineWidth(2);

  hmuon1_pt->SetLineColor(4);
  hmuon1_pt->SetMarkerColor(4);
  hmuon1_pt->SetLineWidth(2);

  hele1_pt->SetLineColor(4);
  hele1_pt->SetMarkerColor(4);
  hele1_pt->SetLineWidth(2);

  hfwdjet1_pt->SetLineColor(4);
  hfwdjet1_pt->SetMarkerColor(4);
  hfwdjet1_pt->SetLineWidth(2);

  hbjet1_pt->SetLineColor(4);
  hbjet1_pt->SetMarkerColor(4);
  hbjet1_pt->SetLineWidth(2);

  hMET_pt->SetLineColor(4);
  hMET_pt->SetMarkerColor(4);
  hMET_pt->SetLineWidth(2);

  h1dipho_eta->SetLineColor(2);
  h1dipho_eta->SetMarkerColor(2);
  h1dipho_eta->SetLineWidth(2);

  h1dipho_leadEta->SetLineColor(2);
  h1dipho_leadEta->SetMarkerColor(2);
  h1dipho_leadEta->SetLineWidth(2);

  h1dipho_subleadEta->SetLineColor(2);
  h1dipho_subleadEta->SetMarkerColor(2);
  h1dipho_subleadEta->SetLineWidth(2);

  h1muon1_eta->SetLineColor(2);
  h1muon1_eta->SetMarkerColor(2);
  h1muon1_eta->SetLineWidth(2);

  h1ele1_eta->SetLineColor(2);
  h1ele1_eta->SetMarkerColor(2);
  h1ele1_eta->SetLineWidth(2);

  h1fwdjet1_eta->SetLineColor(2);
  h1fwdjet1_eta->SetMarkerColor(2);
  h1fwdjet1_eta->SetLineWidth(2);

  h1bjet1_eta->SetLineColor(2);
  h1bjet1_eta->SetMarkerColor(2);
  h1bjet1_eta->SetLineWidth(2);

  h1dipho_pt->SetLineColor(2);
  h1dipho_pt->SetMarkerColor(2);
  h1dipho_pt->SetLineWidth(2);

  h1dipho_leadPt->SetLineColor(2);
  h1dipho_leadPt->SetMarkerColor(2);
  h1dipho_leadPt->SetLineWidth(2);

  h1dipho_subleadPt->SetLineColor(2);
  h1dipho_subleadPt->SetMarkerColor(2);
  h1dipho_subleadPt->SetLineWidth(2);

  h1muon1_pt->SetLineColor(2);
  h1muon1_pt->SetMarkerColor(2);
  h1muon1_pt->SetLineWidth(2);

  h1ele1_pt->SetLineColor(2);
  h1ele1_pt->SetMarkerColor(2);
  h1ele1_pt->SetLineWidth(2);
  
  h1fwdjet1_pt->SetLineColor(2);
  h1fwdjet1_pt->SetMarkerColor(2);
  h1fwdjet1_pt->SetLineWidth(2);

  h1bjet1_pt->SetLineColor(2);
  h1bjet1_pt->SetMarkerColor(2);
  h1bjet1_pt->SetLineWidth(2);

  h1MET_pt->SetLineColor(2);
  h1MET_pt->SetMarkerColor(2);
  h1MET_pt->SetLineWidth(2);

  hdr_tHchainfwdjet->SetLineColor(4);
  hdr_tHchainfwdjet->SetMarkerColor(4);
  hdr_tHchainfwdjet->SetLineWidth(2);

  hdr_bjetfwdjet->SetLineColor(4);
  hdr_bjetfwdjet->SetMarkerColor(4);
  hdr_bjetfwdjet->SetLineWidth(2);

  hdr_subleadphobjet->SetLineColor(4);
  hdr_subleadphobjet->SetMarkerColor(4);
  hdr_subleadphobjet->SetLineWidth(2);

  hdr_leadphofwdjet->SetLineColor(4);
  hdr_leadphofwdjet->SetMarkerColor(4);
  hdr_leadphofwdjet->SetLineWidth(2);

  hdr_subleadphofwdjet->SetLineColor(4);
  hdr_subleadphofwdjet->SetMarkerColor(4);
  hdr_subleadphofwdjet->SetLineWidth(2);

  hdr_leptonbjet->SetLineColor(4);
  hdr_leptonbjet->SetMarkerColor(4);
  hdr_leptonbjet->SetLineWidth(2);

  h1dr_tHchainfwdjet->SetLineColor(2);
  h1dr_tHchainfwdjet->SetMarkerColor(2);
  h1dr_tHchainfwdjet->SetLineWidth(2);

  h1dr_bjetfwdjet->SetLineColor(2);
  h1dr_bjetfwdjet->SetMarkerColor(2);
  h1dr_bjetfwdjet->SetLineWidth(2);

  h1dr_subleadphobjet->SetLineColor(2);
  h1dr_subleadphobjet->SetMarkerColor(2);
  h1dr_subleadphobjet->SetLineWidth(2);

  h1dr_leadphofwdjet->SetLineColor(2);
  h1dr_leadphofwdjet->SetMarkerColor(2);
  h1dr_leadphofwdjet->SetLineWidth(2);

  h1dr_subleadphofwdjet->SetLineColor(2);
  h1dr_subleadphofwdjet->SetMarkerColor(2);
  h1dr_subleadphofwdjet->SetLineWidth(2);

  h1dr_leptonbjet->SetLineColor(2);
  h1dr_leptonbjet->SetMarkerColor(2);                      //4(thq)=blue, 2(tth)=red
  h1dr_leptonbjet->SetLineWidth(2);
  
  h_likelihood_thq->SetLineColor(4);
  h_likelihood_thq->SetMarkerColor(4);
  h_likelihood_thq->SetLineWidth(2);

  h_likelihood_tth->SetLineColor(2);
  h_likelihood_tth->SetMarkerColor(2);
  h_likelihood_tth->SetLineWidth(2);

  gStyle->SetOptStat(0);

  TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9, "");
  leg1->SetFillColor(0);
  leg1->AddEntry(h_likelihood_thq,"tHq","l");
  leg1->AddEntry(h_likelihood_tth,"t#bar{t}H","l");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);  

   TCanvas *dr_comperison= new TCanvas("dr_comperison", "dr_comperison",600, 400);
   dr_comperison->Divide(2,3);

   dr_comperison->cd(1);
   h1dr_tHchainfwdjet->Draw();
   hdr_tHchainfwdjet->Draw("SAME");
   h1dr_tHchainfwdjet->GetXaxis()->SetTitle("#Delta R tHchain and fwdjet");
   h1dr_tHchainfwdjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   dr_comperison->cd(2);
   h1dr_bjetfwdjet->Draw();
   hdr_bjetfwdjet->Draw("SAME");
   h1dr_bjetfwdjet->GetXaxis()->SetTitle("#Delta R bjet and fwdjet");
   h1dr_bjetfwdjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   dr_comperison->cd(3);
   h1dr_subleadphobjet->Draw();
   hdr_subleadphobjet->Draw("SAME");
   h1dr_subleadphobjet->GetXaxis()->SetTitle("#Delta R subleadpho and bjet");
   h1dr_subleadphobjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();   

   dr_comperison->cd(4);
   h1dr_leadphofwdjet->Draw();
   hdr_leadphofwdjet->Draw("SAME");
   h1dr_leadphofwdjet->GetXaxis()->SetTitle("#Delta R  leadpho and fwdjet");
   h1dr_leadphofwdjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();
 
   dr_comperison->cd(5);
   hdr_subleadphofwdjet->Draw();
   h1dr_subleadphofwdjet->Draw("SAME");
   hdr_subleadphofwdjet->GetXaxis()->SetTitle("#Delta R  subleadpho and fwdjet");
   hdr_subleadphofwdjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();   

   dr_comperison->cd(6);
   h1dr_leptonbjet->Draw();
   hdr_leptonbjet->Draw("SAME");
   h1dr_leptonbjet->GetXaxis()->SetTitle("#Delta R lepton and bjet");
   h1dr_leptonbjet->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();   

   TCanvas *eta_comperison= new TCanvas("eta_comperison", "eta_comperison",600, 400);
   eta_comperison->Divide(2,3);

   eta_comperison->cd(1);
   h1dipho_eta->Draw();
   hdipho_eta->Draw("SAME");
   h1dipho_eta->GetXaxis()->SetTitle("dipho #eta");
   h1dipho_eta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   eta_comperison->cd(2);
   h1dipho_leadEta->Draw();
   hdipho_leadEta->Draw("SAME");
   h1dipho_leadEta->GetXaxis()->SetTitle("dipho lead #eta");
   h1dipho_leadEta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();
 
   eta_comperison->cd(3);
   h1dipho_subleadEta->Draw();
   hdipho_subleadEta->Draw("SAME");
   h1dipho_subleadEta->GetXaxis()->SetTitle("dipho lead #eta");
   h1dipho_subleadEta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   eta_comperison->cd(4);
   h1fwdjet1_eta->Draw();
   hfwdjet1_eta->Draw("SAME");
   h1fwdjet1_eta->GetXaxis()->SetTitle("fwdjet1 #eta");
   h1fwdjet1_eta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   eta_comperison->cd(5);
   h1bjet1_eta->Draw();         
   hbjet1_eta->Draw("SAME");
   h1bjet1_eta->GetXaxis()->SetTitle("bjet1 #eta");
   h1bjet1_eta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

//   eta_comperison->cd(6);
//   h1muon1_eta->Draw();
//   hmuon1_eta->Draw("SAME");
//   leg1->Draw();

   eta_comperison->cd(6);
   h1ele1_eta->Draw();
   hele1_eta->Draw("SAME");
   h1ele1_eta->GetXaxis()->SetTitle("ele1 #eta");
   h1ele1_eta->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   TCanvas *pt_comperison= new TCanvas("pt_comperison", "pt_comperison",600, 400);
   pt_comperison->Divide(3,3);

   pt_comperison->cd(1);
   h1dipho_pt->Draw();
   hdipho_pt->Draw("SAME");
   h1dipho_pt->GetXaxis()->SetTitle("dipho P_{T}");
   h1dipho_pt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   pt_comperison->cd(2);
   h1dipho_leadPt->Draw();
   hdipho_leadPt->Draw("SAME");
   h1dipho_leadPt->GetXaxis()->SetTitle("dipho lead P_{T}");
   h1dipho_leadPt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   pt_comperison->cd(3);
   h1dipho_subleadPt->Draw();
   hdipho_subleadPt->Draw("SAME");
   h1dipho_subleadPt->GetXaxis()->SetTitle("dipho sublead P_{T}");
   h1dipho_subleadPt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   pt_comperison->cd(4);
   hfwdjet1_pt->Draw();
   h1fwdjet1_pt->Draw("SAME");
   hfwdjet1_pt->GetXaxis()->SetTitle("fwdjet1 P_{T}");
   hfwdjet1_pt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   pt_comperison->cd(5);
   hbjet1_pt->Draw();
   h1bjet1_pt->Draw("SAME");
   hbjet1_pt->GetXaxis()->SetTitle("bjet1 P_{T}");
   hbjet1_pt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();   

//   pt_comperison->cd(6);
//   h1muon1_pt->Draw();
//   hmuon1_pt->Draw("SAME");
//   leg1->Draw();   

   pt_comperison->cd(6);
   hele1_pt->Draw();
   h1ele1_pt->Draw("SAME");
   hele1_pt->GetXaxis()->SetTitle("ele1 P_{T}");
   hele1_pt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   pt_comperison->cd(7);
   hMET_pt->Draw();
   h1MET_pt->Draw("SAME");
   hMET_pt->GetXaxis()->SetTitle("MET P_{T}");
   hMET_pt->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();

   TCanvas *likelihood_comperison= new TCanvas("likelihood_comperison", "likelihood_comperison",600,400);
   likelihood_comperison->Divide(2,2);
   likelihood_comperison->cd(1);
   h_likelihood_tth->Draw();
   h_likelihood_thq->Draw("SAME");
   h_likelihood_tth->GetXaxis()->SetTitle("Likelihood");
   h_likelihood_tth->GetYaxis()->SetTitle("Normalized Probability");
   leg1->Draw();
//   likelihood_comperison->cd(2);
//   h_likelihood_tth->Draw("SAME");

//   likelihood_comperison->cd(3);
//   LikelihoodClass::LikelihoodClass->h_fstatekinematics_bkg[8]->Draw();

   dr_comperison->SaveAs("dr_comperison_thq_tth.pdf");
   eta_comperison->SaveAs("eta_comperison_thq_tth.pdf");
   pt_comperison->SaveAs("pt_comperison_thq_tth.pdf");
   
   file_output->Write();
//   file_output->Close();

}
