from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT.std import vector
from ROOT import TString
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet

class CR_DiLepton(DiLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        # link flags
        self.channel = None
        self.run_syst = False
        if (super().RunDiMu and super().RunEMu):
            raise ValueError(f"Wrong channel flag {super().RunDiMu} {super().RunEMu}")
        if super().RunDiMu:     self.channel = "RunDiMu"
        if super().RunEMu:      self.channel = "RunEMu"
        if super().RunSyst:     self.run_syst = True
        
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        if self.run_syst:
            if super().RunDiMu:
                self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                          "PileupReweightUp","PileupReweightDown",
                                          "MuonIDSFUp", "MuonIDSFDown",
                                          "DblMuTrigSFUp", "DblMuTrigSFDown"]
            elif super().RunEMu:
                self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                          "PileupReweightUp", "PileupReweightDown",
                                          "MuonIDSFUp", "MuonIDSFDown",
                                          "ElectronIDSFUp", "ElectronIDSFDown",
                                          "EMuTrigSFUp", "EMuTrigSFDown"]
            self.scaleVariations += ["JetResUp", "JetResDown", 
                                     "JetEnUp", "JetEnDown",
                                     "ElectronResUp", "ElectronResDown", 
                                     "ElectronEnUp", "ElectronEnDown",
                                     "MuonEnUp", "MuonEnDown"]
        self.systematics = self.weightVariations + self.scaleVariations
    
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central scale
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
        
        if not channel is None:
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getBaseWeight(channel, ev, jets, truth)
            self.FillObjects(channel, objects, weight, syst="Central_BaseWeight")
            for syst in self.weightVariations:
                weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, truth, syst)
                self.FillObjects(channel, objects, weight, syst)
        
        # Scale variations
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if channel is None: continue
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       }
            weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, truth, syst)
            self.FillObjects(channel, objects, weight, syst)
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        # check the syst argument
        if not syst in self.systematics:
            print(f"[PromptEstimator::defineObjects] Wrong systematics {syst}")
        if syst == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if syst == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if syst == "ElectronResUp":    allElectrons = super().SmearElectrons(allElectrons, +1)
        if syst == "ElectronsResDown": allElectrons = super().SmearElectrons(allElectrons, -1)
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
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)

    def selectEvent(self, event, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        isDiMu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        isEMu = (len(tightMuons) == 1 and len(vetoMuons) == 1 and \
                 len(tightElectrons) == 1 and len(vetoElectrons) == 1)
        
        if self.channel == "RunDiMu":
            if not isDiMu: return None
        if self.channel == "RunEMu":
            if not isEMu: return None
        
        ## DiMu selection
        if self.channel == "RunDiMu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2 = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 50.: return None
            return "DiMu"
        ## EMu selection
        elif self.channel == "RunEMu":
            if not event.PassTrigger(super().EMuTriggers): return None
            mu = tightMuons.at(0)
            ele = tightElectrons.at(0)
            if not ((mu.Pt() > 20. and ele.Pt() > 15.) or (mu.Pt() > 10. and ele.Pt() > 25.)): return None
            if not mu.Charge()+ele.Charge() == 0: return None
            if not mu.DeltaR(ele) > 0.4: return None
            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            return "EMu"
        else:
            print(f"Wrong channel {self.channel}")
    
    def getWeight(self, channel, event, muons, electrons, jets, truth, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            raise ValueError(f"[PromptEstimator::getWeight] Wrong systematic {syst}")

        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            
            if syst == "PileupReweightUp":     w_pileup = super().GetPileUpWeight(super().nPileUp, 1)
            elif syst == "PileupReweightDown": w_pileup = super().GetPileUpWeight(super().nPileUp, -1)
            else:                              w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            
            w_topptweight = 1.
            
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= w_topptweight

            w_muonRecoSF = 1.
            w_muonIDSF = 1.
            w_eleRecoSF = 1.
            w_eleIDSF = 1.
            w_trigSF = 1.
            for mu in muons:
                w_muonRecoSF *= super().getMuonRecoSF(mu, 0)
                if syst == "MuonIDSFUp":     w_muonIDSF *= self.getMuonIDSF(mu, 1)
                elif syst == "MuonIDSFDown": w_muonIDSF *= self.getMuonIDSF(mu, -1)
                else:                        w_muonIDSF *= self.getMuonIDSF(mu, 0)

            for ele in electrons:
                w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0) 
                if syst == "EleIDSFUp":       w_eleIDSF *= self.getEleIDSF(ele, 1);
                elif syst == "EleIDSFDown":   w_eleIDSF *= self.getEleIDSF(ele, -1);
                else:                         w_eleIDSF *= self.getEleIDSF(ele, 0);
            
            if "DiMu" in channel:
                # trigger efficiency
                if syst == "DblMuTrigSFUp":     w_trigSF = self.getDblMuTriggerSF(muons, 1)
                elif syst == "DblMuTrigSFDown": w_trigSF = self.getDblMuTriggerSF(muons, -1)
                else:                           w_trigSF = self.getDblMuTriggerSF(muons, 0)

            if "EMu" in channel:
                if syst == "EMuTrigSFUp":     w_trigSF = self.getEMuTriggerSF(electrons, muons, 1)
                elif syst == "EMuTrigSFDown": w_trigSF = self.getEMuTriggerSF(electrons, muons, -1)
                else:                         w_trigSF = self.getEMuTriggerSF(electrons, muons, 0)

            weight *= w_prefire            #;print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             #;print(f"w_pileup: {w_pileup}")
            weight *= w_muonRecoSF         #;print(f"muonreco: {w_muonRecoSF}")
            weight *= w_muonIDSF           #;print(f"muonid: {w_muonIDSF}")
            weight *= w_eleRecoSF          #;print(f"electronreco: {w_eleRecoSF}")        
            weight *= w_eleIDSF            #;print(f"electronid: {w_eleIDSF}")
            weight *= w_trigSF             #;print(f"trig: {w_trigSF}")

            # b-tagging
            if "EMu" in channel:
                jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
                vjets = vector[Jet]()
                for j in jets: 
                    vjets.emplace_back(j)
                w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
                weight *= w_btag
        
        return weight
    
    # to compare Lepton ID SF / trigger SF w.r.t. baseline
    def getBaseWeight(self, channel, event, jets, truth):
        weight = 1.
        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            weight *= super().GetPrefireWeight(0)
            weight *= super().GetPileUpWeight(super().nPileUp, 0)
            #if "DY" in super().MCSample:
            #    weight *= super().mcCorr.GetOfficialDYReweight(truth, 0)
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                weight *= super().mcCorr.GetTopPtReweight(truth)

            if "EMu" in channel:
                jtp = jParameters(3, 1, 0, 1)
                vjets = vector[Jet]()
                for j in jets: 
                    vjets.emplace_back(j)
                weight *= super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
        
        return weight
    
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
        if "DiMu" in channel:
            pair = muons.at(0) + muons.at(1)
        
        ## fill objects
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        if "DiMu" in channel:
            super().FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 300, 0., 300.)
