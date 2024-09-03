#ifndef MeasTrigEff_h
#define MeasTrigEff_h

#include "AnalyzerCore.h"

class TriggerContainer {
public:
    TriggerContainer() {
    SglMuTrigs.clear(); SglElTrigs.clear();
    DblMuTrigs.clear(); DblMuDZTrigs.clear(); DblMuDZMTrings.clear();
    EMuTrigs.clear(); EMuDZTrigs.clear();
    Mu8El23_ElLegFilter = Mu8El23_MuLegFilter = "";
    Mu23El12_ElLegFilter = Mu23El12_MuLegFilter = "";
    }
    inline void SetSglMuTrigs(const vector<TString>& trigs) { SglMuTrigs = trigs; }
    inline void SetSglElTrigs(const vector<TString>& trigs) { SglElTrigs = trigs; }
    inline void SetDblMuTrigs(const vector<TString>& trigs) { DblMuTrigs = trigs; }
    inline void SetDblMuDZTrigs(const vector<TString>& trigs) { DblMuDZTrigs = trigs; }
    inline void SetDblMuDZMTrigs(const vector<TString>& trigs) { DblMuDZMTrings = trigs; }
    inline void SetEMuTrigs(const vector<TString>& trigs) { EMuTrigs = trigs; }
    inline void SetEMuDZTrigs(const vector<TString>& trigs) { EMuDZTrigs = trigs; }
    inline void SetMu8El23_ElLegFilter(const TString& filter) { Mu8El23_ElLegFilter = filter; }
    inline void SetMu8El23_MuLegFilter(const TString& filter) { Mu8El23_MuLegFilter = filter; }
    inline void SetMu23El12_ElLegFilter(const TString& filter) { Mu23El12_ElLegFilter = filter; }
    inline void SetMu23El12_MuLegFilter(const TString& filter) { Mu23El12_MuLegFilter = filter; }
    vector<TString> GetSglMuTrigs() const { return SglMuTrigs; }
    vector<TString> GetSglElTrigs() const { return SglElTrigs; }
    vector<TString> GetDblMuTrigs() const { return DblMuTrigs; }
    vector<TString> GetDblMuDZTrigs() const { return DblMuDZTrigs; }
    vector<TString> GetDblMuDZMTrigs() const { return DblMuDZMTrings; }
    vector<TString> GetEMuTrigs() const { return EMuTrigs; }
    vector<TString> GetEMuDZTrigs() const { return EMuDZTrigs; }
    TString GetMu8El23_ElLegFilter() const { return Mu8El23_ElLegFilter; }
    TString GetMu8El23_MuLegFilter() const { return Mu8El23_MuLegFilter; }
    TString GetMu23El12_ElLegFilter() const { return Mu23El12_ElLegFilter; }
    TString GetMu23El12_MuLegFilter() const { return Mu23El12_MuLegFilter; }
private:
    vector<TString> SglMuTrigs, SglElTrigs;
    vector<TString> DblMuTrigs, DblMuDZTrigs, DblMuDZMTrings;
    vector<TString> EMuTrigs, EMuDZTrigs;
    TString Mu8El23_ElLegFilter, Mu8El23_MuLegFilter;
    TString Mu23El12_ElLegFilter, Mu23El12_MuLegFilter;
};

class MeasTrigEff: public AnalyzerCore {
public:
    MeasTrigEff();
    ~MeasTrigEff();

    // Userflags
    bool MeasMuLegs, MeasElLegs;
    bool MeasDblMuDZ, MeasEMuDZ;
    // MeasMuLegs -> SingleElectron dataset
    // MeasElLegs -> SingleMuon dataset
    // MeasDblMuDZ -> DoubleMuon dataset
    // MeasEMuDZ -> EMu dataset

    // ID and Trigger definitions
    IDContainer MuID, ElID;
    TriggerContainer Triggers;

    // function to use
    void initializeAnalyzer();
    void executeEvent();
    void measEMuTrigEff_ElLeg(vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight);
    void measEMuTrigEff_MuLeg(vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight);
    void measEMuTrigEff_DZ(Event &ev, vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight);
    void measDblMuTrigEff_DZ(Event &ev, vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight);
};

#endif
