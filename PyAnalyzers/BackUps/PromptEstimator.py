from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet, Particle

from itertools import product
from MLTools.helpers import loadModels
from MLTools.helpers import getDenseInput, getDenseScore
from MLTools.helpers import getGraphInput, getGraphScore

class PromptEstimator(TriLeptonBase):
    def __init__(self):
        super().__init__()  # at this point, TriLeptonBase::initializeAnalyzer has not been called
    
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        self.channel = ""
        self.run_syst = False
        if super().Skim1E2Mu:   self.channel = "Skim1E2Mu"
        if super().Skim3Mu:     self.channel = "Skim3Mu"
        if super().RunSyst:     self.run_syst = True
        if not self.channel in ["Skim1E2Mu", "Skim3Mu"]:
            raise KeyError(f"Wrong channel {self.channel}")
            
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        if self.run_syst:
            if self.channel == "Skim1E2Mu":
                self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                          "PileupReweightUp", "PileupReweightDown",
                                          "MuonIDSFUp", "MuonIDSFDown",
                                          "ElectronIDSFUp", "ElectronIDSFDown",
                                          "EMuTrigSFUp", "EMuTrigSFDown",
                                          "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                          "HeavyTagUpCorr", "HeavyTagDownCorr",
                                          "LightTagUpUnCorr", "LightTagDownUnCorr",
                                          "LightTagUpCorr", "LightTagDownCorr"]
            if self.channel == "Skim3Mu":
                self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                          "PileupReweightUp", "PileupReweightDown",
                                          "MuonIDSFUp", "MuonIDSFDown",
                                          "DblMuTrigSFUp", "DblMuTrigSFDown",
                                          "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                          "HeavyTagUpCorr", "HeavyTagDownCorr",
                                          "LightTagUpUnCorr", "LightTagDownUnCorr",
                                          "LightTagUpCorr", "LightTagDownCorr"]
            self.scaleVariations += ["JetResUp", "JetResDown", 
                                     "JetEnUp", "JetEnDown",
                                     "ElectronResUp", "ElectronResDown", 
                                     "ElectronEnUp", "ElectronEnDown",
                                     "MuonEnUp", "MuonEnDown"] 
        self.systematics = self.weightVariations + self.scaleVariations
        self.signalStrings = ["MHc-100_MA-95", "MHc-130_MA-90", "MHc-160_MA-85"]
        self.backgroundStrings = ["nonprompt", "diboson", "ttZ"]
        self.models = loadModels("GraphNeuralNet", self.channel, self.signalStrings, self.backgroundStrings)
        
    def executeEvent(self):
        if not super().PassMETFilter(): return
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        def processEvent(syst, apply_weight_variation=False):
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv) 
            if channel is None: return
            
            # Evaluate scores
            data, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       "data": data,
                       "scores": scores
            }
            if apply_weight_variation:
                for syst in self.weightVariations:
                    weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, syst)
                    self.FillObjects(channel, objects, weight, syst)
            else:
                weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets)
                self.FillObjects(channel, objects, weight, syst)
        
        processEvent("Central", apply_weight_variation=True)
        for syst in self.scaleVariations:
            processEvent(syst, apply_weight_variation=False)
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first hard copy objects
        allMuons = vector[Muon]()
        allElectrons = vector[Electron]()
        allJets = vector[Jet]()

        for mu in rawMuons: allMuons.emplace_back(mu)
        for ele in rawElectrons: allElectrons.emplace_back(ele)
        for jet in rawJets: allJets.emplace_back(jet)

        # check the syst argument
        if not syst in self.systematics:
            raise KeyError(f"[PromptEstimator::defineObjects] Wrong systematics {syst}")
        if syst == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if syst == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if syst == "ElectronResUp":    allElectrons = super().SmearElectrons(allElectrons, +1)
        if syst == "ElectronResDown":  allElectrons = super().SmearElectrons(allElectrons, -1)
        if syst == "ElectronEnUp":     allElectrons = super().ScaleElectrons(allElectrons, +1)
        if syst == "ElectronEnDown":   allElectrons = super().ScaleElectrons(allElectrons, -1)
        if syst == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if syst == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if syst == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if syst == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
        
        vetoMuons = super().SelectMuons(allMuons, super().MuonIDs[2], 10., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 15., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 15., 2.5)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(0)                  # DeepJet score 3
            wp = super().mcCorr.GetJetTaggingCutValue(0, 2)     # DeepJet Medium 3 1
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and len(tightElectrons) == 1 and len(vetoElectrons) == 1)
        
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None

        # for conversion samples
        if "DYJets" in super().MCSample:
            # at least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = vector[Muon]()
            fakeMuons = vector[Muon]()
            convElectrons = vector[Electron]()
            fakeElectrons = vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
                if super().GetLeptonType(mu, truth) in [-1, -2, -3, -4]: fakeMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
                if super().GetLeptonType(ele, truth) in [-1, -2, -3, -4]: fakeElectrons.emplace_back(ele)
            # remove hadronic contribution
            if self.channel == "Skim1E2Mu":
                if not (fakeMuons.size()+fakeElectrons.size() == 0): return None
                if not (convMuons.size()+convElectrons.size() > 0): return None
            if self.channel == "Skim3Mu":
                if not fakeMuons.size() == 0: return None
                if not convMuons.size() > 0: return None

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 12.: return None 
            if not jets.size() >= 2: return None
            if bjets.size() == 0:
                isOnZ = abs(pair.M() - 91.2) < 10.
                if isOnZ: return "ZFake1E2Mu"
                else:     return None
            else:
                return "SR1E2Mu"
        
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets
        if self.channel == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not jets.size() >= 2: return None

            if bjets.size() == 0:
                isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
                if isOnZ: return "ZFake3Mu"
                else:     return None
            else:
                return "SR3Mu"
        raise EOFError       

    def configureChargeOf(self, muons):
        if not muons.size() == 3:
            raise NotImplementedError(f"wrong no. of muons {muons.size()}")

        # The first two is same charged with pt-order, the third one is opposite charged
        mu1, mu2, mu3 = tuple(muons)
        if mu1.Charge() == mu2.Charge():
            return (mu1, mu2, mu3)
        elif mu1.Charge() == mu3.Charge():
            return mu1, mu3, mu2
        elif mu2.Charge() == mu3.Charge():
            return mu2, mu3, mu1
        else:
            raise EOFError()

    #### Get scores for each event
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        data = getGraphInput(muons, electrons, jets, bjets, METv)
        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            scores[f"{sig}_vs_{bkg}"] = getGraphScore(self.models[f"{sig}_vs_{bkg}"], data)
        return data, scores
    
    #### set weight
    def getWeight(self, channel, event, muons, electrons, jets, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            raise KeyError(f"[PromptEstimator::getWeight] Wrong systematic {syst}")

        if not super().IsDATA:
            weight *= super().MCweight() * super().GetKFactor()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            
            if syst == "PileupReweightUp":     w_pileup = super().GetPileUpWeight(super().nPileUp, 1)
            elif syst == "PileupReweightDown": w_pileup = super().GetPileUpWeight(super().nPileUp, -1)
            else:                              w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            
            w_muonIDSF = 1.
            w_eleIDSF = 1.
            w_trigSF = 1.
            if "1E2Mu" in channel:
                for mu in muons:
                    if syst == "MuonIDSFUp":      w_muonIDSF *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
                    elif syst == "MuonIDSFDown":  w_muonIDSF *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
                    else:                         w_muonIDSF *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
                for el in electrons:
                    if syst == "ElectronIDSFUp":     w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), 1) * self.getEleIDSF(el, 1)
                    elif syst == "ElectronIDSFDown": w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), -1) * self.getEleIDSF(el, -1)
                    else:                            w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), 0) * self.getEleIDSF(el, 0) 

                if syst == "EMuTriggerSFUp":     w_trigSF = self.getEMuTriggerSF(electrons, muons, 1)
                elif syst == "EMuTriggerSFDown": w_trigSF = self.getEMuTriggerSF(electrons, muons, -1)
                else:                            w_trigSF = self.getEMuTriggerSF(electrons, muons, 0)
            
            if "3Mu" in channel:
                for mu in muons:
                    if syst == "MuonIDSFUp":      w_muonIDSF *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
                    elif syst == "MuonIDSFDown":  w_muonIDSF *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
                    else:                         w_muonIDSF *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
                if syst == "DblMuTrigSFUp":     w_trigSF = self.getDblMuTriggerSF(muons, 1)
                elif syst == "DblMuTrigSFDown": w_trigSF = self.getDblMuTriggerSF(muons, -1)
                else:                           w_trigSF = self.getDblMuTriggerSF(muons, 0) 
            
            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             # print(f"w_pileup: {w_pileup}")
            weight *= w_muonIDSF           # print(f"muonID: {w_muonIDSF}")
            weight *= w_eleIDSF;           # print(syst, w_eleIDSF)
            weight *= w_trigSF             # print(f"muontrig: {w_dblMuTrigSF}")

            # b-tagging
            jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            if syst == "HeavyTagUpUnCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTag")
            elif syst == "HeavyTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTag")
            elif syst == "HeavyTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTagCorr")
            elif syst == "HeavyTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTagCorr")
            elif syst == "LightTagUpUnCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTag")
            elif syst == "LightTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTag")
            elif syst == "LightTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTagCorr")
            elif syst == "LightTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTagCorr")
            else:                              w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            weight *= w_btag
        
        return weight
    
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        data = objects["data"]
        scores = objects["scores"]
        
        ## fill base observables
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/energy", mu.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/px", mu.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/py", mu.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/pz", mu.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/charge", mu.Charge(), weight, 3, -1, 2)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/energy", ele.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/px", ele.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/py", ele.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/pz", ele.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/charge", ele.Charge(), weight, 3, -1, 2)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/energy", jet.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/px", jet.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/py", jet.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/pz", jet.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/charge", jet.Charge(), weight, 200, -1, 1)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/energy", bjet.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/px", bjet.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/py", bjet.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/pz", bjet.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/charge", bjet.Charge(), weight, 200, -1, 1)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/btagScore", bjet.GetTaggerResult(3), weight, 100, 0., 1.)
        super().FillHist(f"{channel}/{syst}/muons/size", muons.size(), weight, 10, 0., 10.)
        super().FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), weight, 10, 0., 10.)
        super().FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/Central/METv/energy", METv.E(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/METv/px", METv.Px(), weight, 500, -250., 250.)
        super().FillHist(f"{channel}/{syst}/METv/py", METv.Py(), weight, 500, -250., 250.)
        super().FillHist(f"{channel}/Central/METv/pz", METv.Pz(), weight, 500, -250., 250.)
        
        if "1E2Mu" in channel:
            pair = muons.at(0) + muons.at(1)
            super().FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 200, 0., 200.)
        else:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            super().FillHist(f"{channel}/{syst}/stack/pt", pair1.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/stack/eta", pair1.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/stack/phi", pair1.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/stack/mass", pair1.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/stack/pt", pair2.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/stack/eta", pair2.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/stack/phi", pair2.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/stack/mass", pair2.M(), weight, 200, 0., 200.)

        # Fill ZCands
        if "1E2Mu" in channel:
            ZCand = muons.at(0) + muons.at(1)
            nonprompt = electrons.at(0)
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
        else:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            mZ = 91.2
            if abs(pair1.M() - mZ) < abs(pair2.M() - mZ): 
                ZCand, nZCand = pair1, pair2
                nonprompt = mu_ss2
            else:                                         
                ZCand, nZCand = pair2, pair1
                nonprompt = mu_ss1
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
            
        for signal in self.signalStrings:
            if "1E2Mu" in channel:
                ACand = muons.at(0) + muons.at(1)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
            else:
                mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
                pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
                mA = int(signal.split("_")[1].split("-")[1]) 
                if abs(pair1.M() - mA) < abs(pair2.M() - mA): ACand, nACand = pair1, pair2
                else:                                         ACand, nACand = pair2, pair2
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/pt", nACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/eta", nACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/phi", nACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/mass", nACand.M(), weight, 300, 0., 300.)
                
            score_nonprompt = scores[f"{signal}_vs_nonprompt"]
            score_diboson = scores[f"{signal}_vs_diboson"]
            score_ttZ    = scores[f"{signal}_vs_ttZ"]
            super().FillHist(f"{channel}/{syst}/{signal}/score_nonprompt", score_nonprompt, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/score_diboson", score_diboson, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/score_ttZ", score_ttZ, weight, 100, 0., 1.)
