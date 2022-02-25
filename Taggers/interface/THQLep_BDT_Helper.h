#ifndef _THQLEP_BDT_HELPER_
#define _THQLEP_BDT_HELPER_

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
namespace flashgg {

typedef struct
{
        float dipho_leadPtOvermass_;
        float dipho_subleadPtOvermass_;
        float dipho_leadEta_;
        float dipho_subleadEta_;
        float dipho_leadIDMVA_;
        float dipho_subleadIDMVA_;
        float dipho_lead_haspixelseed_;
        float dipho_sublead_haspixelseed_;
        float n_jets_;
        float n_bjets_;
        float n_centraljets_;
        float lepton_ch_;
        float lepton_leadPt_;
        float lepton_leadEta_;
        float fwdJet1_pt_;
        float fwdJet1_eta_;
        float fwdJet1_discr_;
        float top_mt11_;
        float dRtHchainfwdjet_;
        float dRleptonbjet_;
        float dRleptonfwdjet_;
        float dRbjetfwdjet_;
        float dRleadphofwdjet_;
        float dRsubleadphofwdjet_;
        float bjet1_pt_;
        float bjet2_pt_;
        float bjet3_pt_;
        float bjet1_eta_;
        float bjet2_eta_;
        float bjet3_eta_;
        float bjet1_discr_;
        float bjet2_discr_;
        float bjet3_discr_;
        float jet1_pt_;
        float jet2_pt_;
        float jet3_pt_;
        float jet1_eta_;
        float jet2_eta_;
        float jet3_eta_;
        float jet1_discr_;
        float jet2_discr_;
        float jet3_discr_;
        float Xtt0_;
        float top_mass_;
        float Tprime_mass_;
        float Tprime_mt_;
        float ScalarTPtOverAllPt_;
        float recoMET_pt_;
}  InputVariables;


class THQLep_BDT_Helper {

 public:

  THQLep_BDT_Helper(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~THQLep_BDT_Helper() {}
  double evaluate(const std::string& tag, const InputVariables& varList);
  double convert_tmva_to_prob(double score);
  InputVariables varList_;
  std::unique_ptr<TMVA::Reader> reader_;
};

THQLep_BDT_Helper::THQLep_BDT_Helper(const string& mva_algo, const string& xmlfile)
{

  reader_ = std::unique_ptr<TMVA::Reader>(new TMVA::Reader("!Color:!Silent"));

        reader_->AddVariable( "dipho_leadPt/dipho_mass", &varList_.dipho_leadPtOvermass_ );
        reader_->AddVariable( "dipho_subleadPt/dipho_mass", &varList_.dipho_subleadPtOvermass_ );
        reader_->AddVariable( "dipho_leadEta", &varList_.dipho_leadEta_ );
        reader_->AddVariable( "dipho_subleadEta", &varList_.dipho_subleadEta_ );
        reader_->AddVariable( "dipho_leadIDMVA", &varList_.dipho_leadIDMVA_ );
        reader_->AddVariable( "dipho_subleadIDMVA", &varList_.dipho_subleadIDMVA_ );
        reader_->AddVariable( "dipho_lead_haspixelseed", &varList_.dipho_lead_haspixelseed_ );
        reader_->AddVariable( "dipho_sublead_haspixelseed", &varList_.dipho_sublead_haspixelseed_ );
        reader_->AddVariable( "n_jets"  ,              &varList_.n_jets_);
        reader_->AddVariable( "n_bjets",               &varList_.n_bjets_ );
        reader_->AddVariable( "n_centraljets",         &varList_.n_centraljets_);
        reader_->AddVariable( "lepton_charge",         &varList_.lepton_ch_);
        reader_->AddVariable( "lepton_leadPt",         &varList_.lepton_leadPt_);
        reader_->AddVariable( "lepton_leadEta",         &varList_.lepton_leadEta_);
        reader_->AddVariable( "fwdjet1_pt",           &varList_.fwdJet1_pt_);
//        reader_->AddVariable( "fwdjet1_eta",           &varList_.fwdJet1_eta_);
        reader_->AddVariable( "fwdjet1_discr",           &varList_.fwdJet1_discr_);
//        reader_->AddVariable( "top_mt",                 &varList_.top_mt11_  );
        reader_->AddVariable( "dr_tHchainfwdjet",       &varList_.dRtHchainfwdjet_  );
        reader_->AddVariable( "dr_leptonbjet",          &varList_.dRleptonbjet_  );
        reader_->AddVariable( "dr_leptonfwdjet",        &varList_.dRleptonfwdjet_  );
        reader_->AddVariable( "dr_bjetfwdjet",    &varList_.dRbjetfwdjet_);
        reader_->AddVariable( "dr_leadphofwdjet",       &varList_.dRleadphofwdjet_  );
        reader_->AddVariable( "dr_subleadphofwdjet" ,   &varList_.dRsubleadphofwdjet_);
        reader_->AddVariable( "bjet1_pt",              &varList_.bjet1_pt_);
//        reader_->AddVariable( "bjet2_pt",              &varList_.bjet2_pt_);
//        reader_->AddVariable( "bjet3_pt",              &varList_.bjet3_pt_);
        reader_->AddVariable( "bjet1_eta",              &varList_.bjet1_eta_);
//        reader_->AddVariable( "bjet2_eta",              &varList_.bjet2_eta_);
//        reader_->AddVariable( "bjet3_eta",              &varList_.bjet3_eta_);            
        reader_->AddVariable( "bjet1_discr",            &varList_.bjet1_discr_);
//        reader_->AddVariable( "bjet2_discr",            &varList_.bjet2_discr_);
//        reader_->AddVariable( "bjet3_discr",            &varList_.bjet3_discr_);        
        reader_->AddVariable( "jet1_pt",                &varList_.jet1_pt_);
        reader_->AddVariable( "jet2_pt",                &varList_.jet2_pt_);
//        reader_->AddVariable( "jet3_pt",                &varList_.jet3_pt_);
        reader_->AddVariable( "jet1_eta",               &varList_.jet1_eta_);
        reader_->AddVariable( "jet2_eta",               &varList_.jet2_eta_);
//        reader_->AddVariable( "jet3_eta",               &varList_.jet3_eta_);
//        reader_->AddVariable( "jet1_discr",               &varList_.jet1_discr_);
//        reader_->AddVariable( "jet2_discr",               &varList_.jet2_discr_);
//        reader_->AddVariable( "jet3_discr",               &varList_.jet3_discr_);        
//      reader_->AddVariable( "Xtt0",               &varList_.Xtt0_);
//      //      reader_->AddVariable( "top_mass",               &varList_.top_mass_);
//      //      reader_->AddVariable( "Tprime_mass",               &varList_.Tprime_mass_);
//      //      reader_->AddVariable( "Tprime_mt",               &varList_.Tprime_mt_);
        reader_->AddVariable( "(solvedMET_pt + dipho_pt + lepton_leadPt + bjet1_pt)/(HT + recoMET_pt + lepton_leadPt)",               &varList_.ScalarTPtOverAllPt_);
        reader_->AddVariable( "recoMET_pt", &varList_.recoMET_pt_);
        reader_->BookMVA(mva_algo.c_str(), xmlfile.c_str());

}

double THQLep_BDT_Helper::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here                                                                                                      
  return reader_->EvaluateMVA(mva_algo.c_str());
}

double THQLep_BDT_Helper::convert_tmva_to_prob(double score)
{
  // Undo TMVA transformation
  double raw_score = -0.5 * log( (2 / (score + 1)) - 1);
  // Apply logistic (sigmoid) transformation
  double prob = 1 / (1 + exp(-raw_score));
  return prob;
}

}
#endif



