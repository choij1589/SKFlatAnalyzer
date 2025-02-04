#include "Jet.h"

ClassImp(Jet)

Jet::Jet() : Particle() {
  j_area=-999.;
  j_partonFlavour=-999;
  j_hadronFlavour=-999;
  j_GenHFHadronMatcher_flavour=-999;
  j_GenHFHadronMatcher_origin=-999;
  j_DeepCSV=-999.;
  j_DeepCSV_CvsL=-999.;
  j_DeepCSV_CvsB=-999.;
  j_DeepJet=-999;
  j_DeepJet_CvsL=-999;
  j_DeepJet_CvsB=-999;
  j_chargedHadronEnergyFraction=-999.;
  j_neutralHadronEnergyFraction=-999.;
  j_neutralEmEnergyFraction=-999.;
  j_chargedEmEnergyFraction=-999.;
  j_muonEnergyFraction=-999.;
  j_chargedMultiplicity=-999;
  j_neutralMultiplicity=-999;
  j_PileupJetId=-999.;
  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_bJetNN_corr=1.;
  j_bJetNN_res=-999.;
  j_cJetNN_corr=1.;
  j_cJetNN_res=-999.;
  
  j_tightJetID=false;
  j_tightLepVetoJetID=false;
}

Jet::~Jet(){

}

void Jet::SetArea(double area){
  j_area = area;
}
void Jet::SetGenFlavours(int pf, int hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void Jet::SetGenHFHadronMatcher(int flavour, int origin){
  j_GenHFHadronMatcher_flavour = flavour;
  j_GenHFHadronMatcher_origin = origin; 
}
void Jet::SetTaggerResults(std::vector<double> ds){
  j_DeepCSV           = ds.at(0);
  j_DeepCSV_CvsL      = ds.at(1);
  j_DeepCSV_CvsB      = ds.at(2);
  j_DeepJet       = ds.at(3);
  j_DeepJet_CvsL  = ds.at(4);
  j_DeepJet_CvsB  = ds.at(5);
}
void Jet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void Jet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}
void Jet::SetPileupJetId(double v){
  j_PileupJetId = v;
}

void Jet::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void Jet::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}
void Jet::SetBJetNNCorrection(double bJetNN_corr, double bJetNN_res){
  j_bJetNN_corr = bJetNN_corr;
  j_bJetNN_res = bJetNN_res;
}
void Jet::SetCJetNNCorrection(double cJetNN_corr, double cJetNN_res){
  j_cJetNN_corr = cJetNN_corr;
  j_cJetNN_res = cJetNN_res;
}
void Jet::SetTightJetID(double b){
  j_tightJetID = b;
}
void Jet::SetTightLepVetoJetID(double b){
  j_tightLepVetoJetID = b;
}

// Pileup Jet ID
// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetIDUL for more details
bool Jet::Pass_PileupJetID(const TString era, const TString wp) const {
  // Pileup Jet ID is trained for pt < 50 GeV
  if (Pt() < 10.) {
    std::cerr << "[Jet::Pass_PileupJetID] Pt is too small : " << Pt() << std::endl;
    exit(EXIT_FAILURE);
  }
  if (Pt() > 50.) return true;

  struct PUJetIDCuts {
    double eta2p5, eta2p75, eta3p0, eta5p0;
  };
  std::map<std::pair<double, double>, PUJetIDCuts> ptRangeCuts;
  if (era == "2016preVFP" || era == "2016postVFP") {
    if (wp == "tight") {
      ptRangeCuts[make_pair(10., 20.)] = {0.71, -0.32, -0.30, -0.22};
      ptRangeCuts[make_pair(20., 30.)] = {0.87, -0.08, -0.16, -0.12};
      ptRangeCuts[make_pair(30., 40.)] = {0.94, 0.24, 0.05, 0.10};
      ptRangeCuts[make_pair(40., 50.)] = {0.97, 0.48, 0.26, 0.29};
    } else if (wp == "medium") {
      ptRangeCuts[make_pair(10., 20.)] = {0.20, -0.56, -0.43, -0.38};
      ptRangeCuts[make_pair(20., 30.)] = {0.62, -0.39, -0.32, -0.29};
      ptRangeCuts[make_pair(30., 40.)] = {0.86, -0.10, -0.15, -0.08};
      ptRangeCuts[make_pair(40., 50.)] = {0.93, 0.19, 0.04, 0.12};
    } else if (wp == "loose") {
      ptRangeCuts[make_pair(10., 20.)] = {-0.95, -0.70, -0.52, -0.49};
      ptRangeCuts[make_pair(20., 30.)] = {-0.90, -0.57, -0.43, -0.42};
      ptRangeCuts[make_pair(30., 40.)] = {-0.71, -0.36, -0.29, -0.23};
      ptRangeCuts[make_pair(40., 50.)] = {-0.42, -0.09, -0.14, -0.02};
    } else {
      std::cerr << "[Jet::Pass_PileupJetID] Wrong WP : " << wp << std::endl;
      exit(EXIT_FAILURE);
    }
  } else if (era == "2017" || era == "2018") {
    if (wp == "tight") {
      ptRangeCuts[make_pair(10., 20.)] = {0.77, 0.38, -0.31, -0.21};
      ptRangeCuts[make_pair(20., 30.)] = {0.90, 0.60, -0.12, -0.13};
      ptRangeCuts[make_pair(30., 40.)] = {0.96, 0.82, 0.20, 0.09};
      ptRangeCuts[make_pair(40., 50.)] = {0.98, 0.92, 0.47, 0.29};
    } else if (wp == "medium") {
      ptRangeCuts[make_pair(10., 20.)] = {0.26, -0.33, -0.54, -0.37};
      ptRangeCuts[make_pair(20., 30.)] = {0.68, -0.04, -0.43, -0.30};
      ptRangeCuts[make_pair(30., 40.)] = {0.90, 0.36, -0.16, -0.09};
      ptRangeCuts[make_pair(40., 50.)] = {0.96, 0.61, 0.14, 0.12};
    } else if (wp == "loose") {
      ptRangeCuts[make_pair(10., 20.)] = {-0.95, -0.72, -0.68, -0.47};
      ptRangeCuts[make_pair(20., 30.)] = {-0.88, -0.55, -0.60, -0.43};
      ptRangeCuts[make_pair(30., 40.)] = {-0.63, -0.18, -0.43, -0.24};
      ptRangeCuts[make_pair(40., 50.)] = {-0.19, 0.22, -0.13, -0.03};
    } else {
      std::cerr << "[Jet::Pass_PileupJetID] Wrong WP : " << wp << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cerr << "[Jet::Pass_PileupJetID] Wrong era : " << era << std::endl;
    exit(EXIT_FAILURE);
  }

  // Find applicable pt range
  for (const auto& [ptRange, cuts]: ptRangeCuts) {
    if (Pt() > ptRange.first && Pt() < ptRange.second) {
      const double absEta = fabs(Eta());
      if (absEta < 2.5) {
        return j_PileupJetId > cuts.eta2p5;
      } else if (absEta < 2.75) {
        return j_PileupJetId > cuts.eta2p75;
      } else if (absEta < 3.0) {
        return j_PileupJetId > cuts.eta3p0;
      } else if (absEta < 5.0) {
        return j_PileupJetId > cuts.eta5p0;
      } else {
        std::cerr << "[Jet::Pass_PileupJetID] Wrong eta : " << absEta << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  // EOF, should not reach here
  std::cerr << "[Jet::Pass_PileupJetID] Should not reach here" << std::endl;
  exit(EXIT_FAILURE);
}

bool Jet::PassID(TString ID) const {

  if(ID=="tight") return Pass_tightJetID();
  if(ID=="tightLepVeto") return Pass_tightLepVetoJetID();
  if(ID=="tightWithTightPUID16a") return Pass_tightJetID() && Pass_PileupJetID("2016preVFP", "tight");
  if(ID=="tightWithTightPUID16b") return Pass_tightJetID() && Pass_PileupJetID("2016postVFP", "tight");
  if(ID=="tightWithTightPUID17") return Pass_tightJetID() && Pass_PileupJetID("2017", "tight");
  if(ID=="tightWithTightPUID18") return Pass_tightJetID() && Pass_PileupJetID("2018", "tight");
  if(ID=="tightWithMediumPUID16a") return Pass_tightJetID() && Pass_PileupJetID("2016preVFP", "medium");
  if(ID=="tightWithMediumPUID16b") return Pass_tightJetID() && Pass_PileupJetID("2016postVFP", "medium");
  if(ID=="tightWithMediumPUID17") return Pass_tightJetID() && Pass_PileupJetID("2017", "medium");
  if(ID=="tightWithMediumPUID18") return Pass_tightJetID() && Pass_PileupJetID("2018", "medium");
  if(ID=="tightWithLoosePUID16a") return Pass_tightJetID() && Pass_PileupJetID("2016preVFP", "loose");
  if(ID=="tightWithLoosePUID16b") return Pass_tightJetID() && Pass_PileupJetID("2016postVFP", "loose");
  if(ID=="tightWithLoosePUID17") return Pass_tightJetID() && Pass_PileupJetID("2017", "loose");
  if(ID=="tightWithLoosePUID18") return Pass_tightJetID() && Pass_PileupJetID("2018", "loose");

  cout << "[Jet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}

double Jet::GetTaggerResult(JetTagging::Tagger tg) const {

  if(tg==JetTagging::DeepCSV) return j_DeepCSV;
  else if(tg==JetTagging::DeepCSV_CvsL) return j_DeepCSV_CvsL;
  else if(tg==JetTagging::DeepCSV_CvsB) return j_DeepCSV_CvsB;
  else if(tg==JetTagging::DeepJet) return j_DeepJet;
  else if(tg==JetTagging::DeepJet_CvsL) return j_DeepJet_CvsL;
  else if(tg==JetTagging::DeepJet_CvsB) return j_DeepJet_CvsB;
  else{
    cout << "[Jet::GetTaggerResult] ERROR; Wrong tagger : " << tg << endl;
    return -999;
  }
}

