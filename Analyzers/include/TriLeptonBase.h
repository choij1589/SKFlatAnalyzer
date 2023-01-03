#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"

class TriLeptonBase: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    TH2D *hMuonIDSF;
    TH2D *hMu17Leg1_Data, *hMu17Leg1_MC;
    TH2D *hMu8Leg2_Data, *hMu8Leg2_MC;
    TH2D *hMuonFR, *hElectronFR;
    TriLeptonBase();
    ~TriLeptonBase();
    void initializeAnalyzer();
    void executeEvent();

    double getMuonIDSF(Muon &mu, int sys);
    double getTriggerEff(Muon &mu, TString histkey, bool isDataEff, int sys);
    double getMuonFakeProb(const Muon &mu);
    double getElectronFakeProb(const Electron &ele);
    double getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons);
};

#endif