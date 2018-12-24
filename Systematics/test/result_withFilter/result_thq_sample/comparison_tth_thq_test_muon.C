#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <TChain.h>
#include "TMath.h"

class LikelihoodClass {
public:
    LikelihoodClass();
    ~LikelihoodClass();
    double evaluate_likelihood(std::vector<double> inputvars);
private:
    TH1F * h_fstatekinematics_sig[21];
    TH1F * h_fstatekinematics_bkg[21];
    TFile * file_inputdistributions;
};

LikelihoodClass::LikelihoodClass() {
    file_inputdistributions=new TFile("LikelihoodInput_test_muon.root", "READ");
    h_fstatekinematics_sig[0]= (TH1F*) file_inputdistributions->Get("thq_dipho_pt");
    h_fstatekinematics_sig[1]= (TH1F*) file_inputdistributions->Get("thq_dipho_eta");
    h_fstatekinematics_sig[2]= (TH1F*) file_inputdistributions->Get("thq_dipho_leadPt");
    h_fstatekinematics_sig[3]= (TH1F*) file_inputdistributions->Get("thq_dipho_leadEta");
    h_fstatekinematics_sig[4]= (TH1F*) file_inputdistributions->Get("thq_dipho_subleadPt");
    h_fstatekinematics_sig[5]= (TH1F*) file_inputdistributions->Get("thq_dipho_subleadEta");
    h_fstatekinematics_sig[6]= (TH1F*) file_inputdistributions->Get("thq_muon1_pt");
    h_fstatekinematics_sig[7]= (TH1F*) file_inputdistributions->Get("thq_muon1_eta");
//  h_fstatekinematics_sig[8]= (TH1F*) file_inputdistributions->Get("thq_ele1_pt");
//  h_fstatekinematics_sig[9]= (TH1F*) file_inputdistributions->Get("thq_ele1_eta");
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
    h_fstatekinematics_bkg[6]= (TH1F*) file_inputdistributions->Get("tth_muon1_pt");
    h_fstatekinematics_bkg[7]= (TH1F*) file_inputdistributions->Get("tth_muon1_eta");
//  h_fstatekinematics_bkg[8]= (TH1F*) file_inputdistributions->Get("tth_ele1_pt");
//  h_fstatekinematics_bkg[9]= (TH1F*) file_inputdistributions->Get("tth_ele1_eta");
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
LikelihoodClass::~LikelihoodClass() {
    file_inputdistributions->Close();
}
double LikelihoodClass::evaluate_likelihood(std::vector<double> inputvars) {
    if(inputvars.size()!=19) {
        std::cout<<"<LikelihoodClass::evaluate_likelihood>: inputvars.size() is "<<inputvars.size()<< " and it is expected to be 8!"<<std::endl;
        std::cout<<"<LikelihoodClass::evaluate_likelihood>: is returning a value of -10.!"<<std::endl;
        return -10.;
    }
    double p_signal=1.;
    double p_background=1.;
    for(unsigned int ielement=0; ielement<inputvars.size(); ielement++) {
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


void comparison_tth_thq_test_muon()
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
    TFile *file_output = new TFile("thq_tth_test_muon.root","RECREATE");


    Float_t qdipho_pt, qdipho_eta, qdipho_leadPt, qdipho_leadEta, qdipho_subleadPt, qdipho_subleadEta;
    Float_t qele1_pt, qele1_eta,qele1_phi, qele1_e, qele2_pt, qele2_eta, qele2_phi, qele2_e;
    Float_t qmuon1_pt, qmuon2_pt, qmuon1_eta, qmuon2_eta, qmuon1_phi, qmuon2_phi, qmuon1_e, qmuon2_e;
    Float_t qjet1_pt, qjet1_eta, qfwdjet1_pt, qfwdjet1_eta, qfwdjet1_phi, qfwdjet1_e, qbjet1_pt, qbjet1_eta, qbjet1_phi, qbjet1_e, qbjet2_pt, qbjet2_eta, qbjet2_phi, qbjet2_e, qbjet3_pt, qbjet3_eta, qbjet3_phi, qbjet3_e, qMET_pt, qMET_phi;
    Float_t qdr_tHchainfwdjet, qdr_bjetfwdjet, qdr_subleadphobjet, qdr_leadphofwdjet, qdr_subleadphofwdjet, qdr_leptonbjet;
    Float_t qLeptonType;

//   Float_t rele1_pt, rele1_eta, rele2_pt, rele2_eta;
    Float_t tdipho_pt, tdipho_eta, tdipho_leadPt, tdipho_leadEta, tdipho_subleadPt, tdipho_subleadEta;
    Float_t tele1_pt, tele1_eta, tele1_phi, tele1_e , tele2_pt, tele2_eta , tele2_phi, tele2_e ;
    Float_t tmuon1_pt, tmuon2_pt, tmuon1_eta, tmuon2_eta, tmuon1_phi, tmuon2_phi, tmuon1_e, tmuon2_e;
    Float_t tjet1_pt, tjet1_eta, tfwdjet1_pt, tfwdjet1_eta, tfwdjet1_phi, tfwdjet1_e, tbjet1_pt, tbjet1_eta, tbjet1_phi, tbjet1_e, tbjet2_pt, tbjet2_eta,  tbjet2_phi, tbjet2_e, tbjet3_pt, tbjet3_eta, tbjet3_phi, tbjet3_e, tMET_pt, tMET_phi;
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
    tree_thq->SetBranchAddress("muon1_phi", &qmuon1_phi);
    tree_thq->SetBranchAddress("muon1_e", &qmuon1_e);
    tree_thq->SetBranchAddress("muon2_pt", &qmuon2_pt);
    tree_thq->SetBranchAddress("muon2_eta", &qmuon2_eta);
    tree_thq->SetBranchAddress("muon2_phi", &qmuon2_phi);
    tree_thq->SetBranchAddress("muon2_e", &qmuon2_e);

    tree_thq->SetBranchAddress("n_ele", &qn_ele);
    tree_thq->SetBranchAddress("ele1_pt", &qele1_pt);
    tree_thq->SetBranchAddress("ele1_eta", &qele1_eta);
    tree_thq->SetBranchAddress("ele1_phi", &qele1_phi);
    tree_thq->SetBranchAddress("ele1_e", &qele1_e);
    tree_thq->SetBranchAddress("ele2_pt", &qele2_pt);
    tree_thq->SetBranchAddress("ele2_eta", &qele2_eta);
    tree_thq->SetBranchAddress("ele2_phi", &qele2_phi);
    tree_thq->SetBranchAddress("ele2_e", &qele2_e);
    tree_thq->SetBranchAddress("n_jets", &qn_jets);
    tree_thq->SetBranchAddress("jet1_pt", &qjet1_pt);
    tree_thq->SetBranchAddress("jet1_eta", &qjet1_eta);
    tree_thq->SetBranchAddress("n_fwdjets", &qn_fwdjets);
    tree_thq->SetBranchAddress("fwdjet1_pt", &qfwdjet1_pt);
    tree_thq->SetBranchAddress("fwdjet1_eta", &qfwdjet1_eta);
    tree_thq->SetBranchAddress("fwdjet1_phi", &qfwdjet1_phi);
    tree_thq->SetBranchAddress("fwdjet1_e", &qfwdjet1_e);

    tree_thq->SetBranchAddress("n_bjets", &qn_bjets);
    tree_thq->SetBranchAddress("bjet1_pt", &qbjet1_pt);
    tree_thq->SetBranchAddress("bjet1_eta", &qbjet1_eta);
    tree_thq->SetBranchAddress("bjet1_phi", &qbjet1_phi);
    tree_thq->SetBranchAddress("bjet1_e", &qbjet1_e);
    tree_thq->SetBranchAddress("bjet2_pt", &qbjet2_pt);
    tree_thq->SetBranchAddress("bjet2_eta", &qbjet2_eta);
    tree_thq->SetBranchAddress("bjet2_phi", &qbjet2_phi);
    tree_thq->SetBranchAddress("bjet2_e", &qbjet2_e);
    tree_thq->SetBranchAddress("bjet3_pt", &qbjet3_pt);
    tree_thq->SetBranchAddress("bjet3_eta", &qbjet3_eta);
    tree_thq->SetBranchAddress("bjet3_phi", &qbjet3_phi);
    tree_thq->SetBranchAddress("bjet3_e", &qbjet3_e);
    tree_thq->SetBranchAddress("dr_tHchainfwdjet", &qdr_tHchainfwdjet);
    tree_thq->SetBranchAddress("dr_bjetfwdjet", &qdr_bjetfwdjet);
    tree_thq->SetBranchAddress("dr_subleadphobjet", &qdr_subleadphobjet);
    tree_thq->SetBranchAddress("dr_leadphofwdjet", &qdr_leadphofwdjet);
    tree_thq->SetBranchAddress("dr_subleadphofwdjet", &qdr_subleadphofwdjet);
    tree_thq->SetBranchAddress("dr_leptonbjet", &qdr_leptonbjet);
    tree_thq->SetBranchAddress("LeptonType", &qLeptonType);
    tree_thq->SetBranchAddress("MET_pt", &qMET_pt);
    tree_thq->SetBranchAddress("MET_phi", &qMET_phi);

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
    tree_thq->SetBranchStatus("muon1_phi", 1);
    tree_thq->SetBranchStatus("muon1_e", 1);
    tree_thq->SetBranchStatus("muon2_pt", 1);
    tree_thq->SetBranchStatus("muon2_eta", 1);
    tree_thq->SetBranchStatus("muon2_phi", 1);
    tree_thq->SetBranchStatus("muon2_e", 1);
    tree_thq->SetBranchStatus("n_ele", 1);
    tree_thq->SetBranchStatus("ele1_pt", 1);
    tree_thq->SetBranchStatus("ele1_eta", 1);
    tree_thq->SetBranchStatus("ele1_phi", 1);
    tree_thq->SetBranchStatus("ele1_e", 1);
    tree_thq->SetBranchStatus("ele2_pt", 1);
    tree_thq->SetBranchStatus("ele2_eta", 1);
    tree_thq->SetBranchStatus("ele2_phi", 1);
    tree_thq->SetBranchStatus("ele2_e", 1);
    tree_thq->SetBranchStatus("n_jets", 1);
    tree_thq->SetBranchStatus("jet1_pt", 1);
    tree_thq->SetBranchStatus("jet1_eta", 1);
    tree_thq->SetBranchStatus("n_fwdjets", 1);
    tree_thq->SetBranchStatus("fwdjet1_pt", 1);
    tree_thq->SetBranchStatus("fwdjet1_eta", 1);
    tree_thq->SetBranchStatus("n_bjets", 1);
    tree_thq->SetBranchStatus("bjet1_pt", 1);
    tree_thq->SetBranchStatus("bjet1_eta", 1);
    tree_thq->SetBranchStatus("bjet1_phi", 1);
    tree_thq->SetBranchStatus("bjet1_e", 1);
    tree_thq->SetBranchStatus("bjet2_pt", 1);
    tree_thq->SetBranchStatus("bjet2_eta", 1);
    tree_thq->SetBranchStatus("bjet2_phi", 1);
    tree_thq->SetBranchStatus("bjet2_e", 1);
    tree_thq->SetBranchStatus("bjet3_pt", 1);
    tree_thq->SetBranchStatus("bjet3_eta", 1);
    tree_thq->SetBranchStatus("bjet3_phi", 1);
    tree_thq->SetBranchStatus("bjet3_e", 1);



    tree_thq->SetBranchStatus("dr_tHchainfwdjet", 1);
    tree_thq->SetBranchStatus("dr_bjetfwdjet", 1);
    tree_thq->SetBranchStatus("dr_subleadphobjet", 1);
    tree_thq->SetBranchStatus("dr_leadphofwdjet", 1);
    tree_thq->SetBranchStatus("dr_subleadphofwdjet", 1);
    tree_thq->SetBranchStatus("dr_leptonbjet", 1);
    tree_thq->SetBranchStatus("LeptonType", 1);
    tree_thq->SetBranchStatus("MET_pt", 1);
    tree_thq->SetBranchStatus("MET_phi", 1);

    tree_tth->SetBranchAddress("dipho_pt", &tdipho_pt);
    tree_tth->SetBranchAddress("dipho_eta", &tdipho_eta);
    tree_tth->SetBranchAddress("dipho_leadPt", &tdipho_leadPt);
    tree_tth->SetBranchAddress("dipho_leadEta", &tdipho_leadEta);
    tree_tth->SetBranchAddress("dipho_subleadPt", &tdipho_subleadPt);
    tree_tth->SetBranchAddress("dipho_subleadEta", &tdipho_subleadEta);
    tree_tth->SetBranchAddress("n_muons", &tn_muons);
    tree_tth->SetBranchAddress("muon1_pt", &tmuon1_pt);
    tree_tth->SetBranchAddress("muon1_eta", &tmuon1_eta);
    tree_tth->SetBranchAddress("muon1_phi", &tmuon1_phi);
    tree_tth->SetBranchAddress("muon1_e", &tmuon1_e);
    tree_tth->SetBranchAddress("muon2_pt", &tmuon2_pt);
    tree_tth->SetBranchAddress("muon2_eta", &tmuon2_eta);
    tree_tth->SetBranchAddress("muon2_phi", &tmuon2_phi);
    tree_tth->SetBranchAddress("muon2_e", &tmuon2_e);
    tree_tth->SetBranchAddress("n_ele", &tn_ele);
    tree_tth->SetBranchAddress("ele1_pt",  &tele1_pt);
    tree_tth->SetBranchAddress("ele1_eta", &tele1_eta);
    tree_tth->SetBranchAddress("ele1_phi", &tele1_phi);
    tree_tth->SetBranchAddress("ele1_e",   &tele1_e);
    tree_tth->SetBranchAddress("ele2_pt", &tele2_pt);
    tree_tth->SetBranchAddress("ele2_eta", &tele2_eta);
    tree_tth->SetBranchAddress("ele2_phi", &tele2_phi);
    tree_tth->SetBranchAddress("ele2_e", &tele2_e);
    tree_tth->SetBranchAddress("n_jets", &tn_jets);
    tree_tth->SetBranchAddress("jet1_pt", &tjet1_pt);
    tree_tth->SetBranchAddress("jet1_eta", &tjet1_eta);
    tree_tth->SetBranchAddress("n_fwdjets", &tn_fwdjets);
    tree_tth->SetBranchAddress("fwdjet1_pt", &tfwdjet1_pt);
    tree_tth->SetBranchAddress("fwdjet1_eta", &tfwdjet1_eta);
    tree_tth->SetBranchAddress("fwdjet1_phi", &tfwdjet1_phi);
    tree_tth->SetBranchAddress("fwdjet1_e", &tfwdjet1_e);
    tree_tth->SetBranchAddress("n_bjets", &tn_bjets);
    tree_tth->SetBranchAddress("bjet1_pt", &tbjet1_pt);
    tree_tth->SetBranchAddress("bjet1_eta", &tbjet1_eta);
    tree_tth->SetBranchAddress("bjet1_phi", &tbjet1_phi);
    tree_tth->SetBranchAddress("bjet1_e", &tbjet1_e);
    tree_tth->SetBranchAddress("bjet2_pt", &tbjet2_pt);
    tree_tth->SetBranchAddress("bjet2_eta", &tbjet2_eta);
    tree_tth->SetBranchAddress("bjet2_phi", &tbjet2_phi);
    tree_tth->SetBranchAddress("bjet2_e", &tbjet2_e);
    tree_tth->SetBranchAddress("bjet3_pt", &tbjet3_pt);
    tree_tth->SetBranchAddress("bjet3_eta", &tbjet3_eta);
    tree_tth->SetBranchAddress("bjet3_phi", &tbjet3_phi);
    tree_tth->SetBranchAddress("bjet3_e", &tbjet3_e);
    tree_tth->SetBranchAddress("dr_tHchainfwdjet", &tdr_tHchainfwdjet);
    tree_tth->SetBranchAddress("dr_bjetfwdjet", &tdr_bjetfwdjet);
    tree_tth->SetBranchAddress("dr_subleadphobjet", &tdr_subleadphobjet);
    tree_tth->SetBranchAddress("dr_leadphofwdjet", &tdr_leadphofwdjet);
    tree_tth->SetBranchAddress("dr_subleadphofwdjet", &tdr_subleadphofwdjet);
    tree_tth->SetBranchAddress("dr_leptonbjet", &tdr_leptonbjet);
    tree_tth->SetBranchAddress("LeptonType", &tLeptonType);
    tree_tth->SetBranchAddress("MET_pt", &tMET_pt);
    tree_tth->SetBranchAddress("MET_phi", &tMET_phi);

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
    tree_tth->SetBranchStatus("muon1_phi", 1);
    tree_tth->SetBranchStatus("muon1_e", 1);
    tree_tth->SetBranchStatus("muon2_pt", 1);
    tree_tth->SetBranchStatus("muon2_eta", 1);
    tree_tth->SetBranchStatus("muon2_phi", 1);
    tree_tth->SetBranchStatus("muon2_e", 1);
    tree_tth->SetBranchStatus("n_ele",1);
    tree_tth->SetBranchStatus("ele1_pt", 1);
    tree_tth->SetBranchStatus("ele1_eta", 1);
    tree_tth->SetBranchStatus("ele1_phi", 1);
    tree_tth->SetBranchStatus("ele1_e", 1);
    tree_tth->SetBranchStatus("ele2_pt", 1);
    tree_tth->SetBranchStatus("ele2_eta", 1);
    tree_tth->SetBranchStatus("ele2_phi", 1);
    tree_tth->SetBranchStatus("ele2_e", 1);
    tree_tth->SetBranchStatus("n_jets",1);
    tree_tth->SetBranchStatus("jet1_pt", 1);
    tree_tth->SetBranchStatus("jet1_eta", 1);
    tree_tth->SetBranchStatus("n_fwdjets",1);
    tree_tth->SetBranchStatus("fwdjet1_pt", 1);
    tree_tth->SetBranchStatus("fwdjet1_eta", 1);
    tree_tth->SetBranchStatus("fwdjet1_phi", 1);
    tree_tth->SetBranchStatus("fwdjet1_e", 1);
    tree_tth->SetBranchStatus("n_bjets", 1);
    tree_tth->SetBranchStatus("bjet1_pt", 1);
    tree_tth->SetBranchStatus("bjet1_eta", 1);
    tree_tth->SetBranchStatus("bjet1_phi", 1);
    tree_tth->SetBranchStatus("bjet1_e", 1);
    tree_tth->SetBranchStatus("bjet2_pt", 1);
    tree_tth->SetBranchStatus("bjet2_eta", 1);
    tree_tth->SetBranchStatus("bjet2_phi", 1);
    tree_tth->SetBranchStatus("bjet2_e", 1);
    tree_tth->SetBranchStatus("bjet3_pt", 1);
    tree_tth->SetBranchStatus("bjet3_eta", 1);
    tree_tth->SetBranchStatus("bjet3_phi", 1);
    tree_tth->SetBranchStatus("bjet3_e", 1);
    tree_tth->SetBranchStatus("dr_tHchainfwdjet", 1);
    tree_tth->SetBranchStatus("dr_bjetfwdjet", 1);
    tree_tth->SetBranchStatus("dr_subleadphobjet", 1);
    tree_tth->SetBranchStatus("dr_leadphofwdjet", 1);
    tree_tth->SetBranchStatus("dr_subleadphofwdjet", 1);
    tree_tth->SetBranchStatus("dr_leptonbjet", 1);
    tree_tth->SetBranchStatus("LeptonType", 1);
    tree_tth->SetBranchStatus("MET_pt", 1);
    tree_tth->SetBranchStatus("MET_phi", 1);


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
    TH1F *hdr_tHchainfwdjet= new TH1F("thq_dr_tHchainfwdjet","#Delta R tHchain and fwdjet distribution",50,0,7);
    TH1F *hdr_bjetfwdjet= new TH1F("thq_dr_bjetfwdjet","#Delta R bjet and fwdjet distribution",50,0,7);
    TH1F *hdr_subleadphobjet= new TH1F("thq_dr_subleadphobjet","#Delta R subleadpho and bjet distribution",50,0,7);
    TH1F *hdr_leadphofwdjet= new TH1F("thq_dr_leadphofwdjet","#Delta R leadpho and fwdjet distribution",50,0,7);
    TH1F *hdr_subleadphofwdjet= new TH1F("thq_dr_subleadphofwdjet","dr_subleadphofwdjet distribution",50,0,7);
    TH1F *hdr_leptonbjet= new TH1F("thq_dr_leptonbjet","#Delta R lepton and bjet distribution",50,0,7);

    TH1F *h1dipho_pt= new TH1F("tth_dipho_pt","dipho P_{T} distribution",50,0,600);       //h->thq, h1->tth
    TH1F *h1dipho_eta= new TH1F("tth_dipho_eta","dipho #eta distribution",50,-4,4);
    TH1F *h1dipho_leadPt= new TH1F("tth_dipho_leadPt","dipho lead P_{T} distribution",50,0,600);
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
    TH1F *h1dr_tHchainfwdjet= new TH1F("tth_dr_tHchainfwdjet"," #Delta R tHchain and fwdjet distribution",50,0,7);
    TH1F *h1dr_bjetfwdjet= new TH1F("tth_dr_bjetfwdjet","#Delta R bjet and fwdjet distribution",50,0,7);
    TH1F *h1dr_subleadphobjet= new TH1F("tth_dr_subleadphobjet","#Delta R subleadpho and bjet distribution",50,0,7);
    TH1F *h1dr_leadphofwdjet= new TH1F("tth_dr_leadphofwdjet","#Delta leadpho and fwdjet distribution",50,0,7);
    TH1F *h1dr_subleadphofwdjet= new TH1F("tth_dr_subleadphofwdjet","#Delta R subleadpho and fwdjet distribution",50,0,7);
    TH1F *h1dr_leptonbjet= new TH1F("tth_dr_leptonbjet","#Delta R lepton and bjet distribution",50,0,7);

    TH1F *h_likelihood_thq=new TH1F("likelihood_thq", "Likelihood (mu+jets)" , 50, 0., 1.);
    TH1F *h_likelihood_tth=new TH1F("likelihood_tth", "Likelihood (mu+jets)" , 50, 0., 1.);

    TH1F *h_topmass11_thq=new TH1F("topmass11_thq", "topmass11_thq", 50, 0. , 800);
    TH1F *h_topmass12_thq=new TH1F("topmass12_thq", "topmass12_thq", 50, 0. , 800);
    TH1F *h_topmass13_thq=new TH1F("topmass13_thq", "topmass13_thq", 50, 0. , 800);
    TH1F *h_topmass21_thq=new TH1F("topmass21_thq", "topmass21_thq", 50, 0. , 800);
    TH1F *h_topmass22_thq=new TH1F("topmass22_thq", "topmass22_thq", 50, 0. , 800);
    TH1F *h_topmass23_thq=new TH1F("topmass23_thq", "topmass23_thq", 50, 0. , 800);

    TH1F *h_topmass11_tth=new TH1F("topmass11_tth", "topmass11_tth", 50, 0. , 800);
    TH1F *h_topmass12_tth=new TH1F("topmass12_tth", "topmass12_tth", 50, 0. , 800);
    TH1F *h_topmass13_tth=new TH1F("topmass13_tth", "topmass13_tth", 50, 0. , 800);
    TH1F *h_topmass21_tth=new TH1F("topmass21_tth", "topmass21_tth", 50, 0. , 800);
    TH1F *h_topmass22_tth=new TH1F("topmass22_tth", "topmass22_tth", 50, 0. , 800);
    TH1F *h_topmass23_tth=new TH1F("topmass23_tth", "topmass23_tth", 50, 0. , 800);

    TH1F *h_top_mt11_thq=new TH1F("top_mt11_thq", "top_mt11_thq", 50, 0. , 800);
    TH1F *h_top_mt12_thq=new TH1F("top_mt12_thq", "top_mt12_thq", 50, 0. , 800);
    TH1F *h_top_mt13_thq=new TH1F("top_mt13_thq", "top_mt13_thq", 50, 0. , 800);
    TH1F *h_top_mt21_thq=new TH1F("top_mt21_thq", "top_mt21_thq", 50, 0. , 800);
    TH1F *h_top_mt22_thq=new TH1F("top_mt22_thq", "top_mt22_thq", 50, 0. , 800);
    TH1F *h_top_mt23_thq=new TH1F("top_mt23_thq", "top_mt23_thq", 50, 0. , 800);

    TH1F *h_top_mt11_tth=new TH1F("top_mt11_tth", "top_mt11_tth", 50, 0. , 800);
    TH1F *h_top_mt12_tth=new TH1F("top_mt12_tth", "top_mt12_tth", 50, 0. , 800);
    TH1F *h_top_mt13_tth=new TH1F("top_mt13_tth", "top_mt13_tth", 50, 0. , 800);
    TH1F *h_top_mt21_tth=new TH1F("top_mt21_tth", "top_mt21_tth", 50, 0. , 800);
    TH1F *h_top_mt22_tth=new TH1F("top_mt22_tth", "top_mt22_tth", 50, 0. , 800);
    TH1F *h_top_mt23_tth=new TH1F("top_mt23_tth", "top_mt23_tth", 50, 0. , 800);


    Int_t nentries_thq=(Int_t)tree_thq->GetEntries();
    Int_t nentries_tth=(Int_t)tree_tth->GetEntries();

    Int_t tth_entries=0;
    Int_t thq_entries=0;
    Int_t tth_sel_evnt=0;
    Int_t thq_sel_evnt=0;
    Int_t tth_after_lhoodcut=0;
    Int_t thq_after_lhoodcut=0;

    TLorentzVector l1_thq, l2_thq, b1_thq, b2_thq, b3_thq, metL_thq;
    TLorentzVector l1b1_thq, l1b2_thq, l1b3_thq, l2b1_thq, l2b2_thq, l2b3_thq;
    TLorentzVector l1b1met_thq, l1b2met_thq, l1b3met_thq, l2b1met_thq, l2b2met_thq, l2b3met_thq;
    Float_t m11_thq=0;
    Float_t m12_thq, m13_thq, m21_thq, m22_thq, m23_thq;
    Float_t top_mt11_thq, top_mt12_thq, top_mt13_thq, top_mt21_thq, top_mt22_thq, top_mt23_thq;
    Float_t met_mt_thq;


//thq loop
    for(int ievent=0; ievent < nentries_thq; ievent++)
    {
        tree_thq -> GetEntry( ievent );
        thq_entries++;

        if(/*qn_jets >= 3  && qn_fwdjets == 1 && qn_bjets == 1 &&*/ qn_ele== 0 && qn_muons >= 1) {
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
            hfwdjet1_pt->Fill(qfwdjet1_pt);
            hfwdjet1_eta->Fill(qfwdjet1_eta);
            hbjet1_eta->Fill(qbjet1_eta);
            hbjet1_pt->Fill(qbjet1_pt);
            hMET_pt->Fill(qMET_pt);
            hdr_tHchainfwdjet -> Fill( qdr_tHchainfwdjet );
            hdr_bjetfwdjet->Fill( qdr_bjetfwdjet );
            hdr_subleadphobjet->Fill( qdr_subleadphobjet );
            hdr_leadphofwdjet->Fill( qdr_leadphofwdjet );
            hdr_subleadphofwdjet->Fill( qdr_subleadphofwdjet );
            hdr_leptonbjet->Fill( qdr_leptonbjet );

            l1b1_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b2_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b3_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b1_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b2_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b3_thq.SetPtEtaPhiE(0., 0., 0., 0.);

            l1b1met_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b2met_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b3met_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b1met_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b2met_thq.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b3met_thq.SetPtEtaPhiE(0., 0., 0., 0.);

            l1_thq.SetPtEtaPhiE(qmuon1_pt, qmuon1_eta,qmuon1_phi, qmuon1_e );
            l2_thq.SetPtEtaPhiE(qmuon2_pt, qmuon2_eta,qmuon2_phi, qmuon2_e );
            b1_thq.SetPtEtaPhiE(qbjet1_pt, qbjet1_eta, qbjet1_phi, qbjet1_e );
            b2_thq.SetPtEtaPhiE(qbjet2_pt, qbjet2_eta, qbjet2_phi, qbjet2_e );
            b3_thq.SetPtEtaPhiE(qbjet3_pt, qbjet3_eta, qbjet3_phi, qbjet3_e );

            met_mt_thq = sqrt(qMET_pt*qMET_pt );

            metL_thq.SetPxPyPzE(qMET_pt * cos(qMET_phi), qMET_pt * sin(qMET_phi), 0, met_mt_thq );


            if( qmuon1_pt > -999 && qbjet1_pt > -999)
            {
                l1b1_thq = l1_thq + b1_thq;
            }
            if( qmuon1_pt > -999 && qbjet2_pt > -999)
            {
                l1b2_thq = l1_thq + b2_thq;
            }
            if( qmuon1_pt > -999 && qbjet3_pt > -999)
            {
                l1b3_thq = l1_thq + b3_thq;
            }
            if( qmuon2_pt > -999 && qbjet1_pt > -999)
            {
                l2b1_thq = l2_thq + b1_thq;
            }
            if( qmuon2_pt > -999 && qbjet2_pt > -999)
            {
                l2b2_thq = l2_thq + b2_thq;
            }
            if( qmuon2_pt > -999 && qbjet3_pt > -999)
            {
                l2b3_thq = l2_thq + b3_thq;
            }
//l+b+met
            if( qmuon1_pt > -999 && qbjet1_pt > -999)
            {
                l1b1met_thq = l1_thq + b1_thq + metL_thq;
            }
            if( qmuon1_pt > -999 && qbjet2_pt > -999)
            {
                l1b2met_thq = l1_thq + b2_thq + metL_thq;
            }
            if( qmuon1_pt > -999 && qbjet3_pt > -999)
            {
                l1b3met_thq = l1_thq + b3_thq + metL_thq;
            }
            if( qmuon2_pt > -999 && qbjet1_pt > -999)
            {
                l2b1met_thq = l2_thq + b1_thq + metL_thq;
            }
            if( qmuon2_pt > -999 && qbjet2_pt > -999)
            {
                l2b2met_thq = l2_thq + b2_thq + metL_thq;
            }
            if( qmuon2_pt > -999 && qbjet3_pt > -999)
            {
                l2b3met_thq = l2_thq + b3_thq + metL_thq;
            }

            top_mt11_thq =sqrt((l1b1met_thq.M() * l1b1met_thq.M()) + (l1b1met_thq.Pt() * l1b1met_thq.Pt()));
            top_mt12_thq =sqrt((l1b2met_thq.M() * l1b2met_thq.M()) + (l1b2met_thq.Pt() * l1b2met_thq.Pt()));
            top_mt13_thq =sqrt((l1b3met_thq.M() * l1b3met_thq.M()) + (l1b3met_thq.Pt() * l1b3met_thq.Pt()));
            top_mt21_thq =sqrt((l2b1met_thq.M() * l2b1met_thq.M()) + (l2b1met_thq.Pt() * l2b1met_thq.Pt()));
            top_mt22_thq =sqrt((l2b2met_thq.M() * l2b2met_thq.M()) + (l2b2met_thq.Pt() * l2b2met_thq.Pt()));
            top_mt23_thq =sqrt((l2b3met_thq.M() * l2b3met_thq.M()) + (l2b3met_thq.Pt() * l2b3met_thq.Pt()));


            std::vector<double> topmass_thq;
            if(l1b1_thq.M()>0 )
            {
                topmass_thq.push_back(l1b1_thq.M());
            }
            if(l1b2_thq.M()>0 )
            {
                topmass_thq.push_back(l1b2_thq.M());
            }
            if(l1b3_thq.M()>0 )
            {
                topmass_thq.push_back(l1b3_thq.M());
            }
            if(l2b1_thq.M()>0 )
            {
                topmass_thq.push_back(l2b1_thq.M());
            }
            if(l2b2_thq.M()>0 )
            {
                topmass_thq.push_back(l2b2_thq.M());
            }
            if(l2b3_thq.M()>0 )
            {
                topmass_thq.push_back(l2b3_thq.M());
            }

//     cout<<"topmass_thq.size()"<<topmass_thq.size()<<endl;

            std::sort(topmass_thq.begin(), topmass_thq.end());


//     cout<< "qbjet3_eta="<<qbjet3_phi<<endl;
//     cout<< "qmuon2_eta="<<qmuon2_phi<<endl;
//     cout<< "m11_thq=" << l1b1_thq.M() << endl;
            /*     cout<< "m11_thq=" <<l1b1_thq.M()<< endl;
                 cout<< "m12_thq=" <<l1b2_thq.M()<< endl;
                 cout<< "m13_thq=" <<l1b3_thq.M()<< endl;
                 cout<< "m21_thq=" <<l2b1_thq.M()<< endl;
                 cout<< "m22_thq=" <<l2b2_thq.M()<< endl;
                 cout<< "m23_thq=" <<l2b3_thq_mass<< endl;
            */
//     cout<<"top_mt11="<<top_mt11_thq<<endl;
//     cout<<"top_mt12="<<top_mt12_thq<<endl;
//     cout<<"top_mt13="<<top_mt13_thq<<endl;
//     cout<<"top_mt21="<<top_mt21_thq<<endl;
//     cout<<"top_mt22="<<top_mt22_thq<<endl;
//     cout<<"top_mt23="<<top_mt23_thq<<endl;

//     cout<<"met_mt="<<met_mt_thq<<endl;

            h_topmass11_thq->Fill(l1b1_thq.M());
            h_topmass12_thq->Fill(l1b2_thq.M());
            h_topmass13_thq->Fill(l1b3_thq.M());
            h_topmass21_thq->Fill(l2b1_thq.M());
            h_topmass22_thq->Fill(l2b2_thq.M());
            h_topmass23_thq->Fill(l2b3_thq.M());

            h_top_mt11_thq->Fill(top_mt11_thq);
            h_top_mt12_thq->Fill(top_mt12_thq);
            h_top_mt13_thq->Fill(top_mt13_thq);
            h_top_mt21_thq->Fill(top_mt21_thq);
            h_top_mt22_thq->Fill(top_mt22_thq);
            h_top_mt23_thq->Fill(top_mt23_thq);

//     cout<< "metL_thq=" << metL_thq.M() << endl;

            std::vector<double> vec_lhood_calc;
            vec_lhood_calc.push_back(qdipho_pt);
            vec_lhood_calc.push_back(qdipho_eta);
            vec_lhood_calc.push_back(qdipho_leadPt);
            vec_lhood_calc.push_back(qdipho_leadEta);
            vec_lhood_calc.push_back(qdipho_subleadPt);
            vec_lhood_calc.push_back(qdipho_subleadEta);
            vec_lhood_calc.push_back(qmuon1_pt);
            vec_lhood_calc.push_back(qmuon1_eta);
//     vec_lhood_calc.push_back(qele1_eta);
//     vec_lhood_calc.push_back(qele1_pt);
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
            if(lhood_value>0.5) {
                thq_after_lhoodcut++;
            }

        }
    }//end of thq loop

//tth loop
//

    TLorentzVector l1_tth, l2_tth, b1_tth, b2_tth, b3_tth, metL_tth;
    TLorentzVector l1b1_tth, l1b2_tth, l1b3_tth, l2b1_tth , l2b2_tth, l2b3_tth;
    TLorentzVector l1b1met_tth, l1b2met_tth, l1b3met_tth, l2b1met_tth , l2b2met_tth, l2b3met_tth;
    Float_t m11_tth, m12_tth, m13_tth, m21_tth, m22_tth, m23_tth;
    Float_t top_mt11_tth, top_mt12_tth, top_mt13_tth, top_mt21_tth, top_mt22_tth, top_mt23_tth;
    Float_t met_mt_tth, top_mt_tth;

    for(int ievent=0; ievent<nentries_tth; ievent++)
    {
        tree_tth->GetEntry(ievent);
        tth_entries++;
        if(/*tn_jets >= 3 && tn_fwdjets == 1 && tn_bjets ==1 &&*/ tn_ele == 0 && tn_muons >=1 )
        {

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
            h1fwdjet1_pt->Fill(tfwdjet1_pt);
            h1fwdjet1_eta->Fill(tfwdjet1_eta);
            h1bjet1_eta->Fill(tbjet1_eta);
            h1bjet1_pt->Fill(tbjet1_pt);
            h1MET_pt->Fill(tMET_pt);

            h1dr_tHchainfwdjet -> Fill( tdr_tHchainfwdjet );
            h1dr_bjetfwdjet->Fill( tdr_bjetfwdjet );
            h1dr_subleadphobjet->Fill( tdr_subleadphobjet );
            h1dr_leadphofwdjet->Fill( tdr_leadphofwdjet );
            h1dr_subleadphofwdjet->Fill( tdr_subleadphofwdjet );
            h1dr_leptonbjet->Fill( tdr_leptonbjet );

            l1b1_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b2_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b3_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b1_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b2_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b3_tth.SetPtEtaPhiE(0., 0., 0., 0.);

            l1b1met_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b2met_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l1b3met_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b1met_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b2met_tth.SetPtEtaPhiE(0., 0., 0., 0.);
            l2b3met_tth.SetPtEtaPhiE(0., 0., 0., 0.);


            l1_tth.SetPtEtaPhiE(tmuon1_pt, tmuon1_eta, tmuon1_phi, tmuon1_e );
            l2_tth.SetPtEtaPhiE(tmuon2_pt, tmuon2_eta, tmuon2_phi, tmuon2_e );
            b1_tth.SetPtEtaPhiE(tbjet1_pt, tbjet1_eta, tbjet1_phi, tbjet1_e );
            b2_tth.SetPtEtaPhiE(tbjet2_pt, tbjet2_eta, tbjet2_phi, tbjet2_e );
            b3_tth.SetPtEtaPhiE(tbjet3_pt, tbjet3_eta, tbjet3_phi, tbjet3_e );

            met_mt_tth = sqrt(tMET_pt*tMET_pt );

            metL_tth.SetPxPyPzE(tMET_pt * cos(tMET_phi), tMET_pt * sin(tMET_phi), 0, met_mt_tth );

            if(tmuon1_pt > 0 && tbjet1_pt > 0 )
            {
                l1b1_tth = l1_tth + b1_tth;
            }
            if(tmuon1_pt > 0 && tbjet2_pt > 0 )
            {
                l1b2_tth = l1_tth + b2_tth;
            }
            if(tmuon1_pt > 0 && tbjet3_pt > 0 )
            {
                l1b3_tth = l1_tth + b3_tth;
            }
            if(tmuon2_pt > 0 && tbjet1_pt > 0 )
            {
                l2b1_tth = l2_tth + b1_tth;
            }
            if(tmuon2_pt > 0 && tbjet2_pt > 0 )
            {
                l2b2_tth = l2_tth + b2_tth;
            }
            if(tmuon2_pt > 0 && tbjet3_pt > 0 )
            {
                l2b3_tth = l2_tth + b3_tth;
            }
//l+b+met

            if(tmuon1_pt > 0 && tbjet1_pt > 0 )
            {
                l1b1met_tth = l1_tth + b1_tth + metL_tth;
            }
            if(tmuon1_pt > 0 && tbjet2_pt > 0 )
            {
                l1b2met_tth = l1_tth + b2_tth + metL_tth;
            }
            if(tmuon1_pt > 0 && tbjet3_pt > 0 )
            {
                l1b3met_tth = l1_tth + b3_tth + metL_tth;
            }
            if(tmuon2_pt > 0 && tbjet1_pt > 0 )
            {
                l2b1met_tth = l2_tth + b1_tth + metL_tth;
            }
            if(tmuon2_pt > 0 && tbjet2_pt > 0 )
            {
                l2b2met_tth = l2_tth + b2_tth + metL_tth;
            }
            if(tmuon2_pt > 0 && tbjet3_pt > 0 )
            {
                l2b3met_tth = l2_tth + b3_tth + metL_tth;
            }

            std::vector<double> topmass_tth;

            if(l1b1_tth.M()>0 )
            {
                topmass_tth.push_back(l1b1_thq.M());
            }
            if(l1b2_tth.M()>0 )
            {
                topmass_tth.push_back(l1b2_thq.M());
            }
            if(l1b3_tth.M()>0 )
            {
                topmass_tth.push_back(l1b3_thq.M());
            }
            if(l2b1_tth.M()>0 )
            {
                topmass_tth.push_back(l2b1_thq.M());
            }
            if(l2b2_tth.M()>0 )
            {
                topmass_tth.push_back(l2b2_thq.M());
            }
            if(l2b3_tth.M()>0 )
            {
                topmass_tth.push_back(l2b3_thq.M());
            }

            /*     cout<< "m11_tth=" <<l1b1_tth.M()<< endl;
                 cout<< "m12_tth=" <<l1b2_tth.M()<< endl;
                 cout<< "m13_tth=" <<l1b3_tth.M()<< endl;
                 cout<< "m21_tth=" <<l2b1_tth.M()<< endl;
                 cout<< "m22_tth=" <<l2b2_tth.M()<< endl;
                 cout<< "m23_tth=" <<l2b3_tth.M()<< endl;
            */
//     top_mt_tth =sqrt((l1b1_tth.M() * l1b1_tth.M()) + (l1b1_tth.Pt() * l1b1_tth.Pt())) + met_mt_tth;

//     cout<< "top_mt_tth="<<top_mt_tth<<endl;

            top_mt11_tth =sqrt((l1b1met_tth.M() * l1b1met_tth.M()) + (l1b1met_tth.Pt() * l1b1met_tth.Pt()));
            top_mt12_tth =sqrt((l1b2met_tth.M() * l1b2met_tth.M()) + (l1b2met_tth.Pt() * l1b2met_tth.Pt()));
            top_mt13_tth =sqrt((l1b3met_tth.M() * l1b3met_tth.M()) + (l1b3met_tth.Pt() * l1b3met_tth.Pt()));
            top_mt21_tth =sqrt((l2b1met_tth.M() * l2b1met_tth.M()) + (l2b1met_tth.Pt() * l2b1met_tth.Pt()));
            top_mt22_tth =sqrt((l2b2met_tth.M() * l2b2met_tth.M()) + (l2b2met_tth.Pt() * l2b2met_tth.Pt()));
            top_mt23_tth =sqrt((l2b3met_tth.M() * l2b3met_tth.M()) + (l2b3met_tth.Pt() * l2b3met_tth.Pt()));


            h_topmass11_tth->Fill(l1b1_tth.M());
            h_topmass12_tth->Fill(l1b2_tth.M());
            h_topmass13_tth->Fill(l1b3_tth.M());
            h_topmass21_tth->Fill(l2b1_tth.M());
            h_topmass22_tth->Fill(l2b2_tth.M());
            h_topmass23_tth->Fill(l2b3_tth.M());

            h_top_mt11_tth->Fill(top_mt11_tth);
            h_top_mt12_tth->Fill(top_mt12_tth);
            h_top_mt13_tth->Fill(top_mt13_tth);
            h_top_mt21_tth->Fill(top_mt21_tth);
            h_top_mt22_tth->Fill(top_mt22_tth);
            h_top_mt23_tth->Fill(top_mt23_tth);
            /*
                 cout<<"top_mt13_tth"<<top_mt13_tth<<endl;
                 cout<<"l1b1met_tth.M()"<<l1b1met_tth.M()<<endl;
                 cout<<"l1b2met_tth.M()"<<l1b2met_tth.M()<<endl;
                 cout<<"l1b3met_tth.M()"<<l1b3met_tth.M()<<endl;
                 cout<<"l2b1met_tth.M()"<<l2b1met_tth.M()<<endl;
                 cout<<"l2b2met_tth.M()"<<l2b2met_tth.M()<<endl;
                 cout<<"l2b3met_tth.M()"<<l2b3met_tth.M()<<endl;
            */
            std::vector<double> vec_lhood_calc;
            vec_lhood_calc.push_back(tdipho_pt);
            vec_lhood_calc.push_back(tdipho_eta);
            vec_lhood_calc.push_back(tdipho_leadPt);
            vec_lhood_calc.push_back(tdipho_leadEta);
            vec_lhood_calc.push_back(tdipho_subleadPt);
            vec_lhood_calc.push_back(tdipho_subleadEta);
            vec_lhood_calc.push_back(tmuon1_pt);
            vec_lhood_calc.push_back(tmuon1_eta);
//     vec_lhood_calc.push_back(tele1_eta);
//     vec_lhood_calc.push_back(tele1_pt);
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

//     cout<<"vec_lhood_calc.size()"<<vec_lhood_calc.size()<<endl;
            tth_sel_evnt++;
            double lhood_value=lhoodclass->evaluate_likelihood(vec_lhood_calc);
//     cout<<"vec_lhood_calc.at(18)"<<vec_lhood_calc.at(19)<<endl;
//     cout<<"lhood_value=   "<<lhood_value<<endl;
            if(lhood_value>0.5) {
                tth_after_lhoodcut++;
            }
            h_likelihood_tth->Fill(lhood_value);
        }
    }//end of tth event loop

    cout<<"tth_entries"<<tth_entries<<endl;
    cout<<"thq_entries"<<thq_entries<<endl;
    cout<<"tth_sel_evnt"<<tth_sel_evnt<<endl;
    cout<<"thq_sel_evnt"<<thq_sel_evnt<<endl;
    cout<<"tth_after_lhoodcut"<<tth_after_lhoodcut<<endl;
    cout<<"thq_after_lhoodcut"<<thq_after_lhoodcut<<endl;

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
    hjet1_eta->Scale(1/hjet1_eta->Integral(0, 51));
    hfwdjet1_eta->Scale(1/hfwdjet1_eta->Integral(0, 51));
    hbjet1_eta->Scale(1/hbjet1_eta->Integral(0, 51));

    h1dipho_eta->Scale(1/h1dipho_eta->Integral(0, 51));
    h1dipho_leadEta->Scale(1/h1dipho_leadEta->Integral(0, 51));
    h1dipho_subleadEta->Scale(1/h1dipho_subleadEta->Integral(0, 51));
    h1muon1_eta->Scale(1/h1muon1_eta->Integral(0, 51));
    h1ele1_eta->Scale(1/h1ele1_eta->Integral(0, 51));
    h1jet1_eta->Scale(1/h1jet1_eta->Integral(0, 51));
    h1fwdjet1_eta->Scale(1/h1fwdjet1_eta->Integral(0, 51));
    h1bjet1_eta->Scale(1/h1bjet1_eta->Integral(0, 51));

    hdipho_pt->Scale(1/hdipho_pt->Integral(0, 51));
    hdipho_leadPt->Scale(1/hdipho_leadPt->Integral(0, 51));
    hdipho_subleadPt->Scale(1/hdipho_subleadPt->Integral(0, 51));
    hmuon1_pt->Scale(1/hmuon1_pt->Integral(0, 51));
    hele1_pt->Scale(1/hele1_pt->Integral(0, 51));
    hjet1_pt->Scale(1/hjet1_pt->Integral(0, 51));
    hfwdjet1_pt->Scale(1/hfwdjet1_pt->Integral(0, 51));
    hbjet1_pt->Scale(1/hbjet1_pt->Integral(0, 51));
    hMET_pt->Scale(1/hMET_pt->Integral(0, 51));

    h1dipho_pt->Scale(1/h1dipho_pt->Integral(0, 51));
    h1dipho_leadPt->Scale(1/h1dipho_leadPt->Integral(0, 51));
    h1dipho_subleadPt->Scale(1/h1dipho_subleadPt->Integral(0, 51));
    h1muon1_pt->Scale(1/h1muon1_pt->Integral(0, 51));
    h1ele1_pt->Scale(1/h1ele1_pt->Integral(0, 51));
    h1jet1_pt->Scale(1/h1jet1_pt->Integral(0, 51));
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

    h_topmass11_thq->Scale(1/h_topmass11_thq->Integral(0, 51));
    h_topmass12_thq->Scale(1/h_topmass12_thq->Integral(0, 51));
    h_topmass13_thq->Scale(1/h_topmass13_thq->Integral(0, 51));
    h_topmass21_thq->Scale(1/h_topmass21_thq->Integral(0, 51));
    h_topmass22_thq->Scale(1/h_topmass22_thq->Integral(0, 51));
    h_topmass23_thq->Scale(1/h_topmass23_thq->Integral(0, 51));

    h_topmass11_tth->Scale(1/h_topmass11_tth->Integral(0, 51));
    h_topmass12_tth->Scale(1/h_topmass12_tth->Integral(0, 51));
    h_topmass13_tth->Scale(1/h_topmass13_tth->Integral(0, 51));
    h_topmass21_tth->Scale(1/h_topmass21_tth->Integral(0, 51));
    h_topmass22_tth->Scale(1/h_topmass22_tth->Integral(0, 51));
    h_topmass23_tth->Scale(1/h_topmass23_tth->Integral(0, 51));

    h_top_mt11_thq->Scale(1/h_top_mt11_thq->Integral(0, 51));
    h_top_mt12_thq->Scale(1/h_top_mt12_thq->Integral(0, 51));
    h_top_mt13_thq->Scale(1/h_top_mt13_thq->Integral(0, 51));
    h_top_mt21_thq->Scale(1/h_top_mt21_thq->Integral(0, 51));
    h_top_mt22_thq->Scale(1/h_top_mt22_thq->Integral(0, 51));
    h_top_mt23_thq->Scale(1/h_top_mt23_thq->Integral(0, 51));

    h_top_mt11_tth->Scale(1/h_top_mt11_tth->Integral(0, 51));
    h_top_mt12_tth->Scale(1/h_top_mt12_tth->Integral(0, 51));
    h_top_mt13_tth->Scale(1/h_top_mt13_tth->Integral(0, 51));
    h_top_mt21_tth->Scale(1/h_top_mt21_tth->Integral(0, 51));
    h_top_mt22_tth->Scale(1/h_top_mt22_tth->Integral(0, 51));
    h_top_mt23_tth->Scale(1/h_top_mt23_tth->Integral(0, 51));


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

    hjet1_eta->SetLineColor(4);
    hjet1_eta->SetMarkerColor(4);
    hjet1_eta->SetLineWidth(2);

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

    hjet1_pt->SetLineColor(4);
    hjet1_pt->SetMarkerColor(4);
    hjet1_pt->SetLineWidth(2);

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

    h1jet1_eta->SetLineColor(2);
    h1jet1_eta->SetMarkerColor(2);
    h1jet1_eta->SetLineWidth(2);

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

    h1jet1_pt->SetLineColor(2);
    h1jet1_pt->SetMarkerColor(2);
    h1jet1_pt->SetLineWidth(2);

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
    h1dr_leptonbjet->SetMarkerColor(2);                      //4=blue, 2=red
    h1dr_leptonbjet->SetLineWidth(2);

    h_likelihood_thq->SetLineColor(4);
    h_likelihood_thq->SetMarkerColor(4);
    h_likelihood_thq->SetLineWidth(2);

    h_likelihood_tth->SetLineColor(2);
    h_likelihood_tth->SetMarkerColor(2);
    h_likelihood_tth->SetLineWidth(2);

    h_topmass11_thq->SetLineColor(4);
    h_topmass11_thq->SetMarkerColor(4);
    h_topmass11_thq->SetLineWidth(2);

    h_topmass11_tth->SetLineColor(2);
    h_topmass11_tth->SetMarkerColor(2);
    h_topmass11_tth->SetLineWidth(2);

    h_topmass12_thq->SetLineColor(4);
    h_topmass12_thq->SetMarkerColor(4);
    h_topmass12_thq->SetLineWidth(2);

    h_topmass12_tth->SetLineColor(2);
    h_topmass12_tth->SetMarkerColor(2);
    h_topmass12_tth->SetLineWidth(2);

    h_topmass13_thq->SetLineColor(4);
    h_topmass13_thq->SetMarkerColor(4);
    h_topmass13_thq->SetLineWidth(2);

    h_topmass13_tth->SetLineColor(2);
    h_topmass13_tth->SetMarkerColor(2);
    h_topmass13_tth->SetLineWidth(2);

    h_topmass21_thq->SetLineColor(4);
    h_topmass21_thq->SetMarkerColor(4);
    h_topmass21_thq->SetLineWidth(2);

    h_topmass21_tth->SetLineColor(2);
    h_topmass21_tth->SetMarkerColor(2);
    h_topmass21_tth->SetLineWidth(2);

    h_topmass22_thq->SetLineColor(4);
    h_topmass22_thq->SetMarkerColor(4);
    h_topmass22_thq->SetLineWidth(2);

    h_topmass22_tth->SetLineColor(2);
    h_topmass22_tth->SetMarkerColor(2);
    h_topmass22_tth->SetLineWidth(2);

    h_topmass23_thq->SetLineColor(4);
    h_topmass23_thq->SetMarkerColor(4);
    h_topmass23_thq->SetLineWidth(2);

    h_topmass23_tth->SetLineColor(2);
    h_topmass23_tth->SetMarkerColor(2);
    h_topmass23_tth->SetLineWidth(2);

    h_top_mt11_thq->SetLineColor(4);
    h_top_mt11_thq->SetMarkerColor(4);
    h_top_mt11_thq->SetLineWidth(2);

    h_top_mt11_tth->SetLineColor(2);
    h_top_mt11_tth->SetMarkerColor(2);
    h_top_mt11_tth->SetLineWidth(2);

    h_top_mt12_thq->SetLineColor(4);
    h_top_mt12_thq->SetMarkerColor(4);
    h_top_mt12_thq->SetLineWidth(2);

    h_top_mt12_tth->SetLineColor(2);
    h_top_mt12_tth->SetMarkerColor(2);
    h_top_mt12_tth->SetLineWidth(2);

    h_top_mt13_thq->SetLineColor(4);
    h_top_mt13_thq->SetMarkerColor(4);
    h_top_mt13_thq->SetLineWidth(2);

    h_top_mt13_tth->SetLineColor(2);
    h_top_mt13_tth->SetMarkerColor(2);
    h_top_mt13_tth->SetLineWidth(2);

    h_top_mt21_thq->SetLineColor(4);
    h_top_mt21_thq->SetMarkerColor(4);
    h_top_mt21_thq->SetLineWidth(2);

    h_top_mt21_tth->SetLineColor(2);
    h_top_mt21_tth->SetMarkerColor(2);
    h_top_mt21_tth->SetLineWidth(2);

    h_top_mt22_thq->SetLineColor(4);
    h_top_mt22_thq->SetMarkerColor(4);
    h_top_mt22_thq->SetLineWidth(2);

    h_top_mt22_tth->SetLineColor(2);
    h_top_mt22_tth->SetMarkerColor(2);
    h_top_mt22_tth->SetLineWidth(2);

    h_top_mt23_thq->SetLineColor(4);
    h_top_mt23_thq->SetMarkerColor(4);
    h_top_mt23_thq->SetLineWidth(2);

    h_top_mt23_tth->SetLineColor(2);
    h_top_mt23_tth->SetMarkerColor(2);
    h_top_mt23_tth->SetLineWidth(2);



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
    hdr_leadphofwdjet->Draw();
    h1dr_leadphofwdjet->Draw("SAME");
    hdr_leadphofwdjet->GetXaxis()->SetTitle("#Delta R  leadpho and fwdjet");
    hdr_leadphofwdjet->GetYaxis()->SetTitle("Normalized Probability");
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

//   eta_comperison->cd(4);
//   h1jet1_eta->Draw();
//   hjet1_eta->Draw("SAME");
//   leg1->Draw();

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

    eta_comperison->cd(6);
    h1muon1_eta->Draw();
    hmuon1_eta->Draw("SAME");
    h1muon1_eta->GetXaxis()->SetTitle("muon1 #eta");
    h1muon1_eta->GetYaxis()->SetTitle("Normalized Probability");
    leg1->Draw();

//   eta_comperison->cd(7);
//   h1ele1_eta->Draw();
//   hele1_eta->Draw("SAME");
//   leg1->Draw();

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

//   pt_comperison->cd(4);
//   hjet1_pt->Draw();
//   h1jet1_pt->Draw("SAME");
//   leg1->Draw();

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

    pt_comperison->cd(6);
    hmuon1_pt->Draw();
    h1muon1_pt->Draw("SAME");
    hmuon1_pt->GetXaxis()->SetTitle("ele1 P_{T}");
    hmuon1_pt->GetYaxis()->SetTitle("Normalized Probability");
    leg1->Draw();

//   pt_comperison->cd(7);
//   h1ele1_pt->Draw();
//   hele1_pt->Draw("SAME");
//   leg1->Draw();

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
    TCanvas *topmass_comperison= new TCanvas("topmass", "topmass",600,400);
    topmass_comperison->Divide(2,3);
    topmass_comperison->cd(1);
    h_topmass11_thq->Draw();
    h_topmass11_tth->Draw("SAME");
    topmass_comperison->cd(2);
    h_topmass12_thq->Draw();
    h_topmass12_tth->Draw("SAME");
    topmass_comperison->cd(3);
    h_topmass13_thq->Draw();
    h_topmass13_tth->Draw("SAME");
    topmass_comperison->cd(4);
    h_topmass21_thq->Draw();
    h_topmass21_tth->Draw("SAME");
    topmass_comperison->cd(5);
    h_topmass22_thq->Draw();
    h_topmass22_tth->Draw("SAME");
    topmass_comperison->cd(6);
    h_topmass23_thq->Draw();
    h_topmass23_tth->Draw("SAME");

    TCanvas *top_mt_comperison= new TCanvas("top_mt", "top_mt",600,400);
    top_mt_comperison->Divide(2,3);
    top_mt_comperison->cd(1);
    h_top_mt11_thq->Draw();
    h_top_mt11_tth->Draw("SAME");
    top_mt_comperison->cd(2);
    h_top_mt12_thq->Draw();
    h_top_mt12_tth->Draw("SAME");
    top_mt_comperison->cd(3);
    h_top_mt13_thq->Draw();
    h_top_mt13_tth->Draw("SAME");
    top_mt_comperison->cd(4);
    h_top_mt21_thq->Draw();
    h_top_mt21_tth->Draw("SAME");
    top_mt_comperison->cd(5);
    h_top_mt22_thq->Draw();
    h_top_mt22_tth->Draw("SAME");
    top_mt_comperison->cd(6);
    h_top_mt23_thq->Draw();
    h_top_mt23_tth->Draw("SAME");


    dr_comperison->SaveAs("dr_comperison_thq_tth_muon.pdf");
    eta_comperison->SaveAs("eta_comperison_thq_tth_muon.pdf");
    pt_comperison->SaveAs("pt_comperison_thq_tth_muon.pdf");
    topmass_comperison->SaveAs("topmass_comperison.pdf");
    top_mt_comperison->SaveAs("top_mt_comperison.pdf");
    file_output->Write();
//  file_output->Close();

}
