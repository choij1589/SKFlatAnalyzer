#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"

class TriLeptonBase: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // channel flags 
    bool DenseNet, GraphNet;    // network flags
    bool ScaleVar, WeightVar;   // systematics flags
    bool RunSyst, RunTheoryUnc; // systematics flags
    bool FakeStudy;             // for AcceptanceStudy, if FakeStudy is true, reverse prompt matching will be performed
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;
    vector<TString> MASSPOINTs;

    // Trigger Leg Efficiencies
    TH2D *hMu23Leg_Data, *hMu23Leg_MC;
    TH2D *hMu17Leg_Data, *hMu17Leg_MC;
    TH2D *hMu8Leg_Data, *hMu8Leg_MC;
    TH2D *hEl23Leg_Data, *hEl23Leg_MC;
    TH2D *hEl12Leg_Data, *hEl12Leg_MC;
    // ID SF
    TH2D *hMuonIDSF;
    TH2D *hElIDSF;
    TH2D *hMuFR_Central, *hMuFR_PromptNormUp, *hMuFR_PromptNormDown, *hMuFR_MotherJetPtUp, *hMuFR_MotherJetPtDown, *hMuFR_RequireHeavyTag, *hMuFR_MC;
    TH2D *hElFR_Central, *hElFR_PromptNormUp, *hElFR_PromptNormDown, *hElFR_MotherJetPtUp, *hElFR_MotherJetPtDown, *hElFR_RequireHeavyTag, *hElFR_MC;

    TriLeptonBase();
    ~TriLeptonBase();
    void initializeAnalyzer();
    void executeEvent();

    double getMuonRecoSF(const Muon &mu, int sys);
    double getMuonIDSF(const Muon &mu, int sys);
    double getEleIDSF(const Electron &ele, int sys);
    double getTriggerEff(const Muon &mu, TString histkey, bool isDataEff, int sys);
    double getTriggerEff(const Electron &ele, TString histkey, bool isDATA, int sys);
    double getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys);
    double getDblMuTriggerSF(vector<Muon> &muons, int sys);
    double getEMuTriggerEff(vector<Electron> &electrons, vector<Muon> &muons, bool isDATA, int sys);
    double getEMuTriggerSF(vector<Electron> &electrons, vector<Muon> &muons, int sys);
    double getPairwiseFilterEff(TString filter, bool isDATA);
    double getMuonFakeProb(const Muon &mu, const TString &syst);
    double getElectronFakeProb(const Electron &ele, const TString &syst);
    double getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons, const TString &syst);
};

#endif
