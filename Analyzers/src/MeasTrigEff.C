#include "MeasTrigEff.h"

MeasTrigEff::MeasTrigEff() {}
MeasTrigEff::~MeasTrigEff() {}

void MeasTrigEff::initializeAnalyzer() {
    // Userflags
    MeasMuLegs = HasFlag("MeasMuLegs");
    MeasElLegs = HasFlag("MeasElLegs");
    MeasDblMuDZ = HasFlag("MeasDblMuDZ");
    MeasEMuDZ = HasFlag("MeasEMuDZ");

    // ID settings
    if (DataEra == "2016preVFP") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a");
    } else if (DataEra == "2016postVFP") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b");
    } else if (DataEra == "2017") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight17", "HcToWALoose17", "HcToWAVeto17");
    } else if (DataEra == "2018") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight18", "HcToWALoose18", "HcToWAVeto18");
    } else {
        cerr << "[MeasFakeRateV4::initializeAnalzyer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    // Triggers and filters
    if (DataEra == "2016preVFP") {
        // No Mass filter for DblMuTrigs, No DZ filters for EMu
        Triggers.SetSglMuTrigs({"HLT_IsoMu24_v", "HLT_IsoTkMu24_v"});
        Triggers.SetSglElTrigs({"HLT_Ele27_WPTight_Gsf_v"});
        Triggers.SetDblMuTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"});
        Triggers.SetDblMuDZTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"});
        Triggers.SetDblMuDZMTrigs({});
        Triggers.SetEMuTrigs({"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"});
        Triggers.SetEMuDZTrigs({});
        Triggers.SetMu23El12_MuLegFilter("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23");
        Triggers.SetMu8El23_MuLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8");
        Triggers.SetMu23El12_ElLegFilter("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter");
        Triggers.SetMu8El23_ElLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter");
    } else if (DataEra == "2016postVFP") {
        // No mass filter for DblMuTrigs
        Triggers.SetSglMuTrigs({"HLT_IsoMu24_v", "HLT_IsoTkMu24_v"});
        Triggers.SetSglElTrigs({"HLT_Ele27_WPTight_Gsf_v"});
        Triggers.SetDblMuTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"});
        Triggers.SetDblMuDZTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"});
        Triggers.SetDblMuDZMTrigs({});
        Triggers.SetEMuTrigs({"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"});
        Triggers.SetEMuDZTrigs({"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"});
        Triggers.SetMu23El12_MuLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23");
        Triggers.SetMu8El23_MuLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8");
        Triggers.SetMu23El12_ElLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter");
        Triggers.SetMu8El23_ElLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"); 
    } else if (DataEra == "2017") {
        Triggers.SetSglMuTrigs({"HLT_IsoMu27_v"});
        Triggers.SetSglElTrigs({"HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"});
        Triggers.SetDblMuTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"});
        Triggers.SetDblMuDZTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"});
        Triggers.SetDblMuDZMTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"});
        Triggers.SetEMuTrigs({"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"});
        Triggers.SetEMuDZTrigs({"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"});
        Triggers.SetMu23El12_MuLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23");
        Triggers.SetMu8El23_MuLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8");
        Triggers.SetMu23El12_ElLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter");
        Triggers.SetMu8El23_ElLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"); 
    } else if (DataEra == "2018") {
        Triggers.SetSglMuTrigs({"HLT_IsoMu24_v"});
        Triggers.SetSglElTrigs({"HLT_Ele32_WPTight_Gsf_v"});
        Triggers.SetDblMuTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"});
        Triggers.SetDblMuDZTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"});
        Triggers.SetDblMuDZMTrigs({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"});
        Triggers.SetEMuTrigs({"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"});
        Triggers.SetEMuDZTrigs({"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"});
        Triggers.SetMu23El12_MuLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23");
        Triggers.SetMu8El23_MuLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8");
        Triggers.SetMu23El12_ElLegFilter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter");
        Triggers.SetMu8El23_ElLegFilter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"); 
    } else {
        cerr << "[MeasFakeRateV4::initializeAnalzyer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }
}

void MeasTrigEff::executeEvent() {
    if (!PassMETFilter()) return;
    
    Event ev = GetEvent();
    
    bool PassPreTrigs = false;
    if (MeasMuLegs) {
        PassPreTrigs = ev.PassTrigger(Triggers.GetSglElTrigs());
    } else if (MeasElLegs) {
        PassPreTrigs = ev.PassTrigger(Triggers.GetSglMuTrigs());
    } else if (MeasDblMuDZ) {
        PassPreTrigs = ev.PassTrigger(Triggers.GetDblMuTrigs()) or ev.PassTrigger(Triggers.GetDblMuDZTrigs()) or ev.PassTrigger(Triggers.GetDblMuDZMTrigs());
    } else if (MeasEMuDZ) {
        PassPreTrigs = ev.PassTrigger(Triggers.GetEMuTrigs()) or ev.PassTrigger(Triggers.GetEMuDZTrigs());
    } else {
        cerr << "No valid flag set!" << endl;
        exit(EXIT_FAILURE);
    }
    if (!PassPreTrigs) return;

    // Object Definition
    vector<Muon> rawMuons = GetMuons("NOCUT", 5., 2.4);
    vector<Electron> rawElectrons = GetElectrons("NOCUT", 5., 2.5);
    vector<Jet> rawJets = GetJets("tight", 20., 2.4);
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);
    sort(rawElectrons.begin(), rawElectrons.end(), PtComparing);

    vector<Muon> vetoMuons = SelectMuons(rawMuons, MuID.GetVetoID(), 10., 2.4);
    vector<Muon> tightMuons = SelectMuons(rawMuons, MuID.GetTightID(), 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElID.GetVetoID(), 10., 2.5);
    vector<Electron> tightElectrons = SelectElectrons(rawElectrons, ElID.GetTightID(), 10., 2.5);

    // Weights
    double weight = 1.;
    if (!IsDATA) {
        weight = MCweight()*GetKFactor()*ev.GetTriggerLumi("Full");
        weight *= GetPileUpWeight(nPileUp, 0);
        weight *= GetPrefireWeight(0);
    }

    if (MeasMuLegs) measEMuTrigEff_MuLeg(tightMuons, vetoMuons, tightElectrons, vetoElectrons, weight);
    if (MeasElLegs) measEMuTrigEff_ElLeg(tightMuons, vetoMuons, tightElectrons, vetoElectrons, weight);
    if (MeasEMuDZ) measEMuTrigEff_DZ(ev, tightMuons, vetoMuons, tightElectrons, vetoElectrons, weight);
    if (MeasDblMuDZ) measDblMuTrigEff_DZ(ev, tightMuons, vetoMuons, tightElectrons, vetoElectrons, weight);
}

void MeasTrigEff::measEMuTrigEff_MuLeg(vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight) {
    // Event Selection
    if (! (tightMuons.size() == 1 && vetoMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && vetoElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    if (! (el.Pt() > (DataYear == 2016 ? 30: 35))) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    // Define Binning
    vector<double> Mu23PtBins = {0., 10., 20., 23., 25., 30., 40., 50., 100., 200.};
    vector<double> Mu8PtBins = {0., 5., 8., 10., 15., 20., 30., 50., 100., 200.};
    vector<double> EtaBins = {0., 0.8, 1.6, 2.4};

    // Define the HLT objects
    bool passTagHLT = false;
    // Check if the electron matched to the HLT object
    for (const auto &path: Triggers.GetSglElTrigs()) {
        if (el.PassPath(path)) {passTagHLT = true; break; }
    }
    const bool passMu23Leg = mu.PassFilter(Triggers.GetMu23El12_MuLegFilter());
    const bool passMu8Leg = mu.PassFilter(Triggers.GetMu8El23_MuLegFilter());
    if (! passTagHLT) return;

    FillHist("TrigEff_Mu23El12_MuLeg_DENOM/Central/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu23PtBins);
    FillHist("TrigEff_Mu8El23_MuLeg_DENOM/Central/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu8PtBins);
    if (passMu23Leg) FillHist("TrigEff_Mu23El12_MuLeg_NUM/Central/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu23PtBins);
    if (passMu8Leg)  FillHist("TrigEff_Mu8El23_MuLeg_NUM/Central/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu8PtBins);

    // Syst: QCD Contamination
    if (el.Pt() > 40.) {
        FillHist("TrigEff_Mu23El12_MuLeg_DENOM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu23PtBins);
        FillHist("TrigEff_Mu8El23_MuLeg_DENOM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu8PtBins);
        if (passMu23Leg) FillHist("TrigEff_Mu23El12_MuLeg_NUM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu23PtBins);
        if (passMu8Leg)  FillHist("TrigEff_Mu8El23_MuLeg_NUM/AltTag/fEta_Pt", fabs(mu.Eta()), mu.Pt(), weight, EtaBins, Mu8PtBins);
    }
}

void MeasTrigEff::measEMuTrigEff_ElLeg(vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight) {
    // Event Selection
    if (! (tightMuons.size() == 1 && vetoMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && vetoElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    if (! (mu.Pt() > (DataYear == 2017 ? 29: 26))) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    // Define Binning
    vector<double> El23PtBins = {0.,10.,20.,23.,25.,30.,40.,50.,100.,200.};
    vector<double> El12PtBins = {0.,10.,12.,15.,20.,30.,50.,100.,200.};
    vector<double> EtaBins = {0., 0.8, 1.479, 2.5};

    // Define the HLT objects
    bool passTagHLT = false;
    // Check if the electron matched to the HLT object
    for (const auto &path: Triggers.GetSglMuTrigs()) {
        if (mu.PassPath(path)) {passTagHLT = true; break; }
    }
    const bool passEl23Leg = el.PassFilter(Triggers.GetMu8El23_ElLegFilter());
    const bool passEl12Leg = el.PassFilter(Triggers.GetMu23El12_ElLegFilter());
    if (! passTagHLT) return;

    FillHist("TrigEff_Mu8El23_ElLeg_DENOM/Central/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El23PtBins);
    FillHist("TrigEff_Mu23El12_ElLeg_DENOM/Central/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El12PtBins);
    if (passEl23Leg) FillHist("TrigEff_Mu8El23_ElLeg_NUM/Central/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El23PtBins);
    if (passEl12Leg) FillHist("TrigEff_Mu23El12_ElLeg_NUM/Central/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El12PtBins);

    // Syst: QCD Contamination
    if (mu.Pt() > 35. && mu.RelIso() < 0.15) {
        FillHist("TrigEff_Mu8El23_ElLeg_DENOM/AltTag/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El23PtBins);
        FillHist("TrigEff_Mu23El12_ElLeg_DENOM/AltTag/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El12PtBins);
        if (passEl23Leg) FillHist("TrigEff_Mu8El23_ElLeg_NUM/AltTag/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El23PtBins);
        if (passEl12Leg) FillHist("TrigEff_Mu23El12_ElLeg_NUM/AltTag/fEta_Pt", fabs(el.Eta()), el.Pt(), weight, EtaBins, El12PtBins);
    }
}

void MeasTrigEff::measEMuTrigEff_DZ(Event &ev, vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight) {
    // No DZ filter in 2016a
    if (GetEraShort() == "2016a") return;

    // Event Selection
    if (! (tightMuons.size() == 1 && vetoMuons.size() == 1)) return;
    if (! (tightElectrons.size() == 1 && vetoElectrons.size() == 1)) return;
    const Muon &mu = tightMuons.at(0);
    const Electron &el = tightElectrons.at(0);
    if (! (el.Pt() > 15. && mu.Pt() > 10.)) return;
    if (! (el.Pt() > 25. || mu.Pt() > 25.)) return;
    if (! (mu.Charge()+el.Charge() == 0)) return;
    if (! (mu.DeltaR(el) > 0.4)) return;

    // filter settings
    bool fire_isoMu8El23 = false;
    bool fire_isoMu23El8 = false;
    bool fire_isoMu8El23DZ = false;
    bool fire_isoMu23El8DZ = false;

    const TString isoFilterEl23 = "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter";
    const TString isoFilterEl12 = "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"; 
    const TString isoFilterMu23 = "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23";
    const TString isoFilterMu8 = "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8";
    const TString pathMu8El23 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
    const TString pathMu23El12 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
    const TString pathMu8El23DZ = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v";
    const TString pathMu23El12DZ = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

    if (ev.PassTrigger(pathMu8El23) && el.Pt() > 25. && mu.Pt() > 10.) {
        if (el.PassFilter(isoFilterEl23) && mu.PassFilter(isoFilterMu8)) fire_isoMu8El23 = true;
        if (fire_isoMu8El23 && el.PassPath(pathMu8El23DZ) && mu.PassPath(pathMu8El23DZ)) fire_isoMu8El23DZ = true;
    } else if (ev.PassTrigger(pathMu23El12) && el.Pt() > 15. && mu.Pt() > 25.) {
        if (el.PassFilter(isoFilterEl12) && mu.PassFilter(isoFilterMu23)) fire_isoMu23El8 = true;
        if (fire_isoMu23El8 && el.PassPath(pathMu23El12DZ) && mu.PassPath(pathMu23El12DZ)) fire_isoMu23El8DZ = true;
    }

    if (fire_isoMu8El23) {
        FillHist("TrigEff_Mu8El23DZ", 0., weight, 2, 0., 2.);
        if (fire_isoMu8El23DZ) FillHist("TrigEff_Mu8El23DZ", 1., weight, 2, 0., 2.);
    }
    if (fire_isoMu23El8) {
        FillHist("TrigEff_Mu23El8DZ", 0., weight, 2, 0., 2.);
        if (fire_isoMu23El8DZ) FillHist("TrigEff_Mu23El8DZ", 1., weight, 2, 0., 2.);
    }
    if (fire_isoMu8El23 || fire_isoMu23El8) {
        FillHist("TrigEff_EMuDZ", 0., weight, 2, 0., 2.);
        if ((fire_isoMu8El23 && fire_isoMu8El23DZ) || (fire_isoMu23El8 && fire_isoMu23El8DZ)) FillHist("TrigEff_EMuDZ", 1., weight, 2, 0., 2.);
    }
}

void MeasTrigEff::measDblMuTrigEff_DZ(Event &ev, vector<Muon> &tightMuons, vector<Muon> &vetoMuons, vector<Electron> &tightElectrons, vector<Electron> &vetoElectrons, double &weight) {
    // Event Selection
    if (! (tightMuons.size() == 2 && vetoMuons.size() == 2)) return;
    if (! (tightElectrons.size() == 0 && vetoElectrons.size() == 0)) return;
    const Muon &mu1 = tightMuons.at(0);
    const Muon &mu2 = tightMuons.at(1);
    if (! (mu1.Pt() > 20. && mu2.Pt() > 10.)) return;
    if (! (mu1.Charge()+mu2.Charge() == 0)) return;

    double Mll = (mu1+mu2).M();
    double dRll = mu1.DeltaR(mu2);

    bool fire_iso = false;
    bool fire_isoDZ = false;
    bool fire_isoDZM = false;
    bool fire_isoM = false;
    for (unsigned int i=0; i<Triggers.GetDblMuTrigs().size(); i++) {
        const TString path = Triggers.GetDblMuTrigs().at(i);
        const TString pathDZ = Triggers.GetDblMuDZTrigs().at(i);
        const TString pathDZM = (DataYear > 2016) ? Triggers.GetDblMuDZMTrigs().at(i) : "";
        if (! ev.PassTrigger(path)) continue;
        
        unsigned int nMatchIso=0, nMatchIsoDZ=0, nMatchIsoDZM=0, nMatchIsoM=0;
        for (const auto &mu: tightMuons) {
            if (mu.PassPath(path)) nMatchIso++;
            if (mu.PassPath(path) && mu.PassPath(pathDZ)) nMatchIsoDZ++;
            if (mu.PassPath(path) && mu.PassPath(pathDZ) && mu.PassPath(pathDZM)) nMatchIsoDZM++;
            if (mu.PassPath(path) && mu.PassPath(pathDZM)) nMatchIsoM++;
        }
        fire_iso = fire_iso || (nMatchIso == 2);
        fire_isoDZ = fire_isoDZ || (nMatchIsoDZ == 2);
        fire_isoDZM = fire_isoDZM || (nMatchIsoDZM == 2);
        fire_isoM = fire_isoM || (nMatchIsoM == 2);
    }

    bool fire_DZ = false;
    bool fire_DZM = false;
    for (unsigned int i = 0; i < Triggers.GetDblMuDZTrigs().size(); i++) {
        const TString pathDZ = Triggers.GetDblMuDZTrigs().at(i);
        const TString pathDZM = (DataYear > 2016) ? Triggers.GetDblMuDZMTrigs().at(i) : "";
        if (! ev.PassTrigger(pathDZ)) continue;

        unsigned int nMatchDZ=0, nMatchDZM=0;
        for (const auto &mu: tightMuons) {
            if (mu.PassPath(pathDZ)) nMatchDZ++;
            if (mu.PassPath(pathDZ) && mu.PassPath(pathDZM)) nMatchDZM++;   
        }
        fire_DZ = fire_DZ || (nMatchDZ == 2);
        fire_DZM = fire_DZM || (nMatchDZM == 2);
    }

    bool fire_M = false;
    for (unsigned int i = 0; i < Triggers.GetDblMuDZMTrigs().size(); i++) {
        const TString pathDZM = Triggers.GetDblMuDZMTrigs().at(i);
        if (! ev.PassTrigger(pathDZM)) continue;

        unsigned int nMatchDZM=0;
        for (const auto &mu: tightMuons) {
            if (mu.PassPath(pathDZM)) nMatchDZM++;
        }
        fire_M = fire_M || (nMatchDZM == 2);
    }

    vector<TString> SelTagList;
    if(Mll>60 && Mll<120 && dRll>0.3) SelTagList.emplace_back("_ZWin60");
    if(fabs(Mll-91.2)<10 && dRll>0.3) SelTagList.emplace_back("_ZWin20");

    for (const auto &selTag: SelTagList) {
        if (fire_iso) {
            FillHist("TrigEff_Iso"+selTag, 0., weight, 1., 0., 1);
            FillHist("TrigEff_Iso_DZ1"+selTag, fabs(mu1.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_Iso_DZ2"+selTag, fabs(mu2.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_Iso_DZ"+selTag, fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
        }
        if (fire_isoDZ) {
            FillHist("TrigEff_IsoDZ"+selTag, 0., weight, 1., 0., 1);
            FillHist("TrigEff_IsoDZ_DZ1"+selTag, fabs(mu1.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_IsoDZ_DZ2"+selTag, fabs(mu2.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_IsoDZ_DZ"+selTag, fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
        }
        if (fire_DZ) FillHist("TrigEff_DZ"+selTag, 0., weight, 1., 0., 1.);
        if (fire_isoDZM) FillHist("TrigEff_IsoDZM"+selTag, 0., weight, 1., 0., 1.);
        if (fire_isoM) {
            FillHist("TrigEff_IsoM"+selTag, 0., weight, 1., 0., 1);
            FillHist("TrigEff_IsoM_DZ1"+selTag, fabs(mu1.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_IsoM_DZ2"+selTag, fabs(mu2.dZ()), weight, 40, 0., 0.4);
            FillHist("TrigEff_IsoM_DZ"+selTag, fabs(mu1.dZ()-mu2.dZ()), weight, 40, 0., 0.4);
        }
        if (fire_DZM) FillHist("TrigEff_DZM"+selTag, 0., weight, 1., 0., 1.);
        if (fire_M) FillHist("TrigEff_M"+selTag, 0., weight, 1., 0., 1.);
    }
    if (fire_DZ) {
        FillHist("TrigEff_DZ_MllNearCut", Mll, weight, 100, 0., 10.);
        FillHist("TrigEff_DZ_Mll", Mll, weight, 19, 10., 200.);
        if (fire_DZM) {
            FillHist("TrigEff_DZM_MllNearCut", Mll, weight, 100, 0., 10.);
            FillHist("TrigEff_DZM_Mll", Mll, weight, 19, 10., 200.);
        }
    }
}

