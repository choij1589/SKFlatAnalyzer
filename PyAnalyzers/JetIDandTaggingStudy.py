from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet

class JetIDandTaggingStudy(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
        
        # Flags
        self.skim = ""
        if self.Skim1E2Mu: self.skim = "Skim1E2Mu"
        if self.Skim3Mu:   self.skim = "Skim3Mu"
        if self.skim not in ["Skim1E2Mu", "Skim3Mu"]:
            raise ValueError("Invalid skim option")
        
        # Only central systematics for this analyzer
        self.variations = [
            "Central",
            "MediumPUID",
            "LoosePUID",
            "TightBtagging",
            "LooseBtagging"
        ]
        
    def executeEvent(self):
        if not self.PassMETFilter(): return
        ev = self.GetEvent()
        rawMuons = self.GetAllMuons()
        rawElectrons = self.GetAllElectrons()
        rawJets = self.GetAllJets()
        METv = ev.GetMETVector()
        rawObjects = {
            "muons": rawMuons,
            "electrons": rawElectrons,
            "jets": rawJets,
            "METv": METv
        }
        
        def processEvent(variation):
            objects = self.defineObjects(rawObjects, variation)
            channel = self.selectEvent(ev, objects)
            if not channel: return
            
            weight = self.getWeight(ev)
            self.fillObjects(channel, objects, weight, variation)
            
        for variation in self.variations:
            processEvent(variation)
    
    def defineObjects(self, rawObjects, variation):
        jetID = ""
        # Define Jet ID
        if variation == "MediumPUID":
            if self.DataEra == "2016preVFP":
                jetID = "tightWithMediumPUID16a"
            if self.DataEra == "2016postVFP":
                jetID = "tightWithMediumPUID16b"
            if self.DataEra == "2017":
                jetID = "tightWithMediumPUID17"
            if self.DataEra == "2018":
                jetID = "tightWithMediumPUID18"
        elif variation == "LoosePUID":
            if self.DataEra == "2016preVFP":
                jetID = "tightWithLoosePUID16a"
            if self.DataEra == "2016postVFP":
                jetID = "tightWithLoosePUID16b"
            if self.DataEra == "2017":
                jetID = "tightWithLoosePUID17"
            if self.DataEra == "2018":
                jetID = "tightWithLoosePUID18"
        else:
            jetID = "tight"
        
        # first hard copy objects
        allMuons = vector[Muon]()
        allElectrons = vector[Electron]()
        allJets = vector[Jet]()
        
        for mu in rawObjects["muons"]: allMuons.emplace_back(mu)
        for ele in rawObjects["electrons"]: allElectrons.emplace_back(ele)
        for jet in rawObjects["jets"]: allJets.emplace_back(jet)
        
        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs[2], 10., 2.4)
        tightMuons = self.SelectMuons(vetoMuons, self.MuonIDs[0], 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs[2], 15., 2.5)
        tightElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs[0], 15., 2.5)
        jets = self.SelectJets(allJets, jetID, 20., 2.4)
        jets = self.JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)
            if variation == "LooseBtagging":
                wp = self.mcCorr.GetJetTaggingCutValue(3, 0)
            elif variation == "TightBtagging":
                wp = self.mcCorr.GetJetTaggingCutValue(3, 2)
            else:
                wp = self.mcCorr.GetJetTaggingCutValue(3, 1)
            
            if btagScore > wp: bjets.emplace_back(jet)

        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        objects = {
            "vetoMuons": vetoMuons,
            "tightMuons": tightMuons,
            "vetoElectrons": vetoElectrons,
            "tightElectrons": tightElectrons,
            "jets": jets,
            "bjets": bjets,
            "METv": rawObjects["METv"]
        }
        return objects
    
    def selectEvent(self, event, objects):
        vetoMuons = objects["vetoMuons"]
        tightMuons = objects["tightMuons"]
        vetoElectrons = objects["vetoElectrons"]
        tightElectrons = objects["tightElectrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]

        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and len(tightElectrons) == 1 and len(vetoElectrons) == 1)
        
        if self.skim == "Skim1E2Mu" and not is1E2Mu: return
        if self.skim == "Skim3Mu" and not is3Mu: return
        
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.skim == "Skim1E2Mu":
            if not event.PassTrigger(self.EMuTriggers): return
            leptons = vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return
            if not mu1.Charge()+mu2.Charge() == 0: return
            pair = mu1 + mu2
            if not pair.M() > 12.: return  
            if not jets.size() >= 2: return
            if bjets.size() > 0: return "SR1E2Mu"
            
            # No b region
            isOnZ = abs(pair.M() - 91.2) < 10.
            if isOnZ: return "ZFake1E2Mu"
            else:     return
            
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets
        if self.skim == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return 
            if not mu2.Pt() > 10.: return 
            if not mu3.Pt() > 10.: return 
            
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return 
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(tightMuons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            if not pair1.M() > 12.: return 
            if not pair2.M() > 12.: return 
            if not jets.size() >= 2: return 
            if bjets.size() > 0: return "SR3Mu"
            
            # No b region
            isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
            if isOnZ: return "ZFake3Mu"
            else:     return
            
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
            raise EOFError(f"wrong charge configuration {mu1.Charge()} {mu2.Charge()} {mu3.Charge()}")
        
    def getWeight(self, event):
        weight = 1.
        if self.IsDATA: return weight
        
        weight *= self.MCweight() * self.GetKFactor()
        weight *= event.GetTriggerLumi("Full")
        weight *= self.GetPrefireWeight(0)
        weight *= self.GetPileUpWeight(self.nPileUp, 0)
        
        return weight
    
    def fillObjects(self, channel, objects, weight, variation):
        muons = objects["tightMuons"]
        electrons = objects["tightElectrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{variation}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{variation}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{variation}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{variation}/jets/{idx}/pileupScore", jet.PileupJetId(), weight, 200, -1., 1.)
            self.FillHist(f"{channel}/{variation}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
        self.FillHist(f"{channel}/{variation}/muons/size", muons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{variation}/electrons/size", electrons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{variation}/jets/size", jets.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{variation}/bjets/size", bjets.size(), weight, 10, 0., 10.)
