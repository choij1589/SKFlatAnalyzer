#ifndef MeasFakeRateV2_h
#define MeasFakeRateV2_h

#include "AnalyzerCore.h"
#include "NonpromptParameter.h"

class IDContainerV2 {
public:
    IDContainerV2() { Tight = Loose = Veto = "";}
    inline void SetIDs(const TString& tight, const TString& loose, const TString& veto) {
        Tight = tight; 
        Loose = loose; 
        Veto = veto; 
    }
    TString GetTightID() { return Tight; }
    TString GetLooseID() { return Loose; }
    TString GetVetoID() { return Veto; }
private:
    TString Tight;
    TString Loose;
    TString Veto;
};

class MeasFakeRateV2: public AnalyzerCore {
public:
    MeasFakeRateV2();
    ~MeasFakeRateV2();

    // Userflags
    bool MeasFakeMu, MeasFakeMu8, MeasFakeMu17;
    bool MeasFakeEl, MeasFakeEl8, MeasFakeEl12, MeasFakeEl23;
    bool RunSyst, RunSystSimple;
    
    // binnings
    vector<double> ptcorr_bins, abseta_bins;

    // ID Definitions
    IDContainerV2 MuID, ElID;

    // Trigger Definitions
    TString isoSglLepTrig;
    double trigSafePtCut;

    // To reuse the containers
    vector<Muon> rawMuons;
    vector<Electron> rawElectrons;
    vector<Jet> rawJets;
    vector<Gen> truth;
    
    // systematic variations
    vector<TString> systematics;
    vector<TString> weightVariations;
    vector<TString> scaleVariations;
    vector<TString> selectionVariations;

    // for nPV reweighting
    TH1D *hNPVData, *hNPVMC;

    // functions to use
    void initializeAnalyzer();
    void executeEvent();
    void executeEventFrom(NonpromptParameter &param);
    void ApplyScaleVariation(vector<Muon> &muons, vector<Electron> &electrons, vector<Jet> &jets, const TString &scale);
    double GetJetPtCut(const TString &selection);
    TString SelectEvent(const vector<Muon> &muons, const vector<Muon> &vetoMuons, const vector<Electron> &electrons, const vector<Electron> &vetoElectrons, const vector<Jet> &jets);
    double GetNPVReweight(unsigned int nPV);
    double GetEventWeight(const NonpromptParameter &param, Event &ev, vector<Muon> &muons, vector<Electron> &electrons, vector<Jet> &jets);
    TString FindBin(const double ptcorr, const double abseta);
    void FillObjects(const TString &channel, const TString &prefix, const vector<Muon> &muons, const vector<Electron> &electrons, const vector<Jet> &jets, const vector<Jet> &bjets, const Particle &METv, const double &weight);

};

#endif
