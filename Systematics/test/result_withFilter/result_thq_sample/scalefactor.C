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



void scalefactor()

{

/*TFile *file_tHq = new TFile("output_THQ_comb.root","READ");
TTree *tree_thq = (TTree*) file_tHq->Get("tagsDumper/trees/THQ_13TeV_THQLeptonicTag");
*/
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


    Float_t qdipho_pt, qdipho_eta, qdipho_leadPt, qdipho_leadEta, qdipho_subleadPt, qdipho_subleadEta;
    Float_t qele1_pt, qele1_eta,qele1_phi, qele1_e, qele2_pt, qele2_eta, qele2_phi, qele2_e;
    Float_t qmuon1_pt, qmuon2_pt, qmuon1_eta, qmuon2_eta, qmuon1_phi, qmuon2_phi, qmuon1_e, qmuon2_e;
    Float_t qjet1_pt, qjet1_eta, qfwdjet1_pt, qfwdjet1_eta, qfwdjet1_phi, qfwdjet1_e, qbjet1_pt, qbjet1_eta, qbjet1_phi, qbjet1_e, qbjet2_pt, qbjet2_eta, qbjet2_phi, qbjet2_e, qbjet3_pt, qbjet3_eta, qbjet3_phi, qbjet3_e, qMET_pt, qMET_phi;
    Float_t qdr_tHchainfwdjet, qdr_bjetfwdjet, qdr_subleadphobjet, qdr_leadphofwdjet, qdr_subleadphofwdjet, qdr_leptonbjet;
    Float_t qLeptonType;
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

    tree_thq->SetBranchStatus("*", 1);

    TH1F *hbjet1_pt= new TH1F("thq_bjet1_pt","bjet1 P_{T} distribution",50,0,600);

    Int_t nentries_thq=(Int_t)tree_thq->GetEntries();
    Int_t event_passed=0;

    for(int ievent=0; ievent < nentries_thq; ievent++)
    {
        tree_thq -> GetEntry( ievent );
	if(/*abs(qmuon1_eta)<2.4 &&*/ qmuon1_pt>5/* && qdipho_leadPt>30 && qjet1_pt>25 && abs(qjet1_eta)<2.4 && abs(qdipho_leadEta)<2.5*/){
        hbjet1_pt->Fill(qbjet1_pt);
        event_passed++;
	}    
        
    }
cout<<"nentries_thq"<<nentries_thq<<endl;
cout<<"event_passed if loop"<<event_passed<<endl;
TCanvas *c1= new TCanvas("c1", "c1",600, 400);
c1->Divide(1,1);
hbjet1_pt->Draw();
}
