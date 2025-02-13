import pandas as pd
from ROOT import gSystem
from ROOT import DataPreprocess
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class DefineTrainData(DataPreprocess):
    def __init__(self):
        super().__init__()
    
    def executeEvent(self):
        mHc = int(str(super().MCSample).split("_")[1].split("-")[1])
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()

        #### Object definition
        vetoMuons = super().SelectMuons(rawMuons, super().MuonIDs[2], 10., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(rawElectrons, super().ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(rawJets, "tight", 15., 2.4)
        # jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            score = jet.GetTaggerResult(3)                      # DeepJet
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if score > wp: bjets.emplace_back(jet)

        # sort objects
        sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(jets, key=lambda x: x.Pt(), reverse=True)
        sorted(bjets, key=lambda x: x.Pt(), reverse=True)

        #### event selection
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        if not (is3Mu or is1E2Mu): return None

        # prompt matching
        if not super().IsDATA:
            truth = super().GetGens()
            promptMuons = vector[Muon]()
            promptElectrons = vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

            if len(promptMuons) != len(tightMuons): return None
            if len(promptElectrons) != len(tightElectrons): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        if is1E2Mu:
            #if not ev.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            #passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            #passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            #passSafeCut = passLeadMu or passLeadEle
            #if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 12.: return None

            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            channel = "SR1E2Mu"
            # orthogonality of SR and CR done by bjet multiplicity
            #if len(bjets) >= 1:
            #    if len(jets) >= 2: channel = "SR1E2Mu"
            #    else: return None
            #else:
            #    mZ = 91.2
            #    isOnZ = abs(pair.M() - mZ) < 10.
            #    if isOnZ: channel = "ZFake1E2Mu"
            #    else:
            #        if abs((mu1+mu2+ele).M() - mZ) < 10.: channel = "ZGamma1E2Mu"
            #        else: return None

        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        else:
            #if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            #if not mu1.Pt() > 20.: return None
            #if not mu2.Pt() > 10.: return None
            #if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None

            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            channel = "SR3Mu"
            # orthogonality of SR and CR done by bjet multiplicity
            #if len(bjets) >= 1:
            #    if len(jets) >= 2: channel = "SR3Mu"
            #    else: return None
            #else:
            #    mZ = 91.2
            #    isOnZ = abs(pair1.M() - mZ) < 10. or abs(pair2.M() - mZ) < 10.
            #    if isOnZ: channel = "ZFake3Mu"
            #    else:
            #        if abs((mu1+mu2+mu3).M() - mZ) < 10.: channel = "ZGamma3Mu"
            #        else: return None

        if not ("1E2Mu" in channel or "3Mu" in channel): return None

        # start matching
        chargedDecays = super().findChargedDecays(truth)
        if not len(chargedDecays) == 3:
            # print([gen.PID() for gen in chargedDecays])
            return None
        
        # Find ACand first
        signalColl, promptColl = [], []
        for mu in tightMuons:
            if super().GetLeptonType(mu, truth) == 2:
                signalColl.append(mu)
            elif super().GetLeptonType(mu, truth) > 0:
                promptColl.append(mu)
            else:
                continue
        if not len(signalColl) == 2: return None
        mu1, mu2 = tuple(signalColl)
        ACand = mu1+mu2

        if 0 < abs(chargedDecays[1].PID()) and abs(chargedDecays[1].PID()) < 6:   # A q q
            super().FillHist("Ajj/CutFlow", 0., 1., 10, 0., 10.)
            # partons should be inside acceptance
            p1, p2 = tuple(chargedDecays[1:])
            if not p1.Pt() > 15.: return None
            if not p2.Pt() > 15.: return None
            if not abs(p1.Eta()) < 2.5: return None
            if not abs(p2.Eta()) < 2.5: return None
            super().FillHist("Ajj/CutFlow", 1., 1., 10, 0., 10.)

            # find nearest jets
            j1 = None; dR1 = 5.
            j2 = None; dR2 = 5.
            for jet in jets:
                if p1.DeltaR(jet) < dR1:
                    j1 = jet; dR1 = p1.DeltaR(jet)
                if p2.DeltaR(jet) < dR2:
                    j2 = jet; dR2 = p2.DeltaR(jet)
            if not j1: return None
            if not j2: return None
            if j1 is j2: return None
            super().FillHist("Ajj/CutFlow", 2., 1., 10, 0., 10.)
            WCand = j1 + j2
            ChargedHiggs = ACand + WCand
            super().FillHist("Ajj/Acceptance/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/Acceptance/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/Acceptance/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Ajj/Acceptance/dRj1", dR1, 1., 500, 0., 5.)
            super().FillHist("Ajj/Acceptance/dRj2", dR2, 1., 500, 0., 5.)

            # charged Higgs mass cut
            mHc = float(str(super().MCSample).split("_")[1].split("-")[1])
            if abs(ChargedHiggs.M() - mHc) > 20.: return None
            super().FillHist("Ajj/CutFlow", 3., 1., 10, 0., 10.)
            super().FillHist("Ajj/MassCut/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/MassCut/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/MassCut/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Ajj/MassCut/dRj1", dR1, 1., 500, 0., 5.)
            super().FillHist("Ajj/MassCut/dRj2", dR2, 1., 500, 0., 5.)

            # matching cut
            if dR1 > 0.3: return None
            if dR2 > 0.3: return None
            super().FillHist("Ajj/CutFlow", 4., 1., 10, 0., 10.)
            super().FillHist("Ajj/Final/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/Final/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Ajj/Final/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Ajj/Final/dRj1", dR1, 1., 500, 0., 5.)
            super().FillHist("Ajj/Final/dRj2", dR2, 1., 500, 0., 5.)

        elif 11 in [abs(gen.PID()) for gen in chargedDecays]:       # A e nu
            super().FillHist("Aenu/CutFlow", 0., 1., 10, 0., 10.)
            if not "1E2Mu" in channel: return None
            super().FillHist("Aenu/CutFlow", 1., 1., 10, 0., 10.)

            eleGen = None
            nuGen = None
            for gen in chargedDecays:
                if abs(gen.PID()) == 11: eleGen = gen
                if abs(gen.PID()) == 12: nuGen = gen
            if not eleGen.Pt() > 8.: return None
            if not abs(eleGen.Eta()) < 2.6: return None
            ele = tightElectrons.at(0)
            WGen = eleGen + nuGen
            WCand = ele + METv
            ChargedHiggs = ACand + WCand
            super().FillHist("Aenu/NoCut/dRele", ele.DeltaR(eleGen), 1., 500, 0., 500.)
            super().FillHist("Aenu/NoCut/mWgen", WGen.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/NoCut/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/NoCut/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/NoCut/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Aenu/NoCut/MT", WCand.Mt(), 1., 200, 0., 200.)

            if not super().GetLeptonType(ele, truth) > 0: return None
            if not ele.DeltaR(eleGen) < 0.1: return None
            super().FillHist("Aenu/CutFlow", 2., 1., 10, 0., 10.)
            super().FillHist("Aenu/Final/dRele", ele.DeltaR(eleGen), 1., 500, 0., 500.)
            super().FillHist("Aenu/Final/mWgen", WGen.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/Final/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/Final/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Aenu/Final/MT", WCand.Mt(), 1., 200, 0., 200.)
            super().FillHist("Aenu/Final/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Aenu/Final/deta", ele.Eta() - nuGen.Eta(), 1., 2000, -10, 10)
            super().FillHist("Aenu/Final/dphi", ele.Phi() - nuGen.Phi(), 1., 2000, -10, 10)
        elif 13 in [abs(gen.PID()) for gen in chargedDecays]:       # A mu nu
            super().FillHist("Amunu/CutFlow", 0., 1., 10, 0., 10.)
            if not "3Mu" in channel: return None
            super().FillHist("Amunu/CutFlow", 1., 1., 10, 0., 10.)
            
            if not len(promptColl) == 1: return None
            super().FillHist("Amunu/CutFlow", 2., 1., 10, 0., 10.)
            muGen = None
            nuGen = None
            for gen in chargedDecays:
                if abs(gen.PID()) == 13: muGen = gen
                if abs(gen.PID()) == 14: nuGen = gen
            
            if not muGen.Pt() > 8.: return None
            if not abs(muGen.Eta()) < 2.5: return None
            promptMu = promptColl[0]
            WGen = muGen + nuGen
            WCand = promptMu + METv
            ChargedHiggs = ACand + WCand
            super().FillHist("Amunu/NoCut/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/NoCut/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/NoCut/mWgen", WGen.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/NoCut/MT", WCand.Mt(), 1., 200, 0., 200.)
            super().FillHist("Amunu/NoCut/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Amunu/NoCut/dR", promptMu.DeltaR(muGen), 1., 500, 0., 5)

            if not promptMu.DeltaR(muGen) < 0.1: return None
            super().FillHist("Amunu/CutFlow", 3., 1., 10, 0., 10.)
            super().FillHist("Amunu/Final/mA", ACand.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/Final/mW", WCand.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/Final/mWgen", WGen.M(), 1., 200, 0., 200.)
            super().FillHist("Amunu/Final/MT", WCand.Mt(), 1., 200, 0., 200.)
            super().FillHist("Amunu/Final/mHc", ChargedHiggs.M(), 1., 500, 0., 500.)
            super().FillHist("Amunu/Final/deta", promptMu.Eta() - nuGen.Eta(), 1., 2000, -10, 10)
            super().FillHist("Amunu/Final/dphi", promptMu.Phi() - nuGen.Phi(), 1., 2000, -10, 10)
        else:       # A tau nu case
            return None
        
if __name__ == "__main__":
    m = DefineTrainData()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
    #if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-160_MA-15_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/221012_081454/0000/SKFlatNtuple_2017_MC_15.root"): exit(1)
    #if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-160_MA-15_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/221012_081454/0000/SKFlatNtuple_2017_MC_1.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializeAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
