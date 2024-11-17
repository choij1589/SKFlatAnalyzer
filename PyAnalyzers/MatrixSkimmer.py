from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TTree
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet

from array import array
from itertools import product
from MLTools.helpers import loadParticleNet
from MLTools.helpers import getGraphInput, getGraphScore

class MatrixSkimmer(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
        
        ## channels
        if self.Skim1E2Mu: self.skim = "Skim1E2Mu"
        if self.Skim3Mu: self.skim = "Skim3Mu"
        if self.skim not in ["Skim1E2Mu", "Skim3Mu"]:
            raise ValueError("Invalid skim option")
        
        self.network = "ParticleNet"
        self.sigStrings = ["MHc-160_MA-85", "MHc-130_MA-90", "MHc-100_MA-95"]
        self.bkgStrings = ["nonprompt", "diboson", "ttZ"]
        self.models = loadParticleNet("Combined", self.sigStrings, self.bkgStrings)
        self.__prepareTTree()
        
    def executeEvent(self):
        if not self.PassMETFilter(): return
        ev = self.GetEvent()
        rawMuons = self.GetAllMuons()
        rawElectrons = self.GetAllElectrons()
        rawJets = self.GetAllJets()
        METv = ev.GetMETVector()
        
        # initialize contents
        self.__initTreeContents()
        
        # fill contents
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)

        if not channel: return
        pairs = self.makePair(looseMuons)
        _, scores, fold = self.evalScore(looseMuons, looseElectrons, jets, bjets, METv)
        
        if channel == "SR1E2Mu":
            self.mass1["Central"][0] = pairs.M()  
            self.mass2["Central"][0] = -999.
        if channel == "SR3Mu":
            self.mass1["Central"][0] = pairs[0].M()
            self.mass2["Central"][0] = pairs[1].M()
        
        for SIG in self.sigStrings:
            self.scoreX[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_nonprompt"]
            self.scoreY[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_diboson"]
            self.scoreZ[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_ttZ"]
        self.fold["Central"][0] = fold
        self.weight["Central"][0] = self.getFakeWeight(looseMuons, looseElectrons, syst="Central")
        self.tree["Central"].Fill()
    
    def __prepareTTree(self):
        self.tree = {}
        self.mass1 = {}
        self.mass2 = {}
        self.scoreX = {}
        self.scoreY = {}
        self.scoreZ = {}
        self.fold = {}
        self.weight = {}

        tree = TTree(f"Events_Central", "")
        self.mass1["Central"] = array("d", [0.]); tree.Branch("mass1", self.mass1["Central"], "mass1/D")
        self.mass2["Central"] = array("d", [0.]); tree.Branch("mass2", self.mass2["Central"], "mass2/D")
        for SIG in self.sigStrings:
            # vs nonprompt
            self.scoreX[f"{SIG}_Central"] = array("d", [0.])
            tree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_Central"], f"score_{SIG}_vs_nonprompt/D")
            # vs diboson
            self.scoreY[f"{SIG}_Central"] = array("d", [0.])
            tree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_Central"], f"score_{SIG}_vs_diboson/D")
            # vs ttZ
            self.scoreZ[f"{SIG}_Central"] = array("d", [0.])
            tree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_Central"], f"score_{SIG}_vs_ttZ/D")
        self.fold["Central"] = array("i", [0]); tree.Branch("fold", self.fold["Central"], "fold/I")
        self.weight["Central"] = array("d", [0.]); tree.Branch("weight", self.weight["Central"], "weight/D")
        tree.SetDirectory(0)
        self.tree["Central"] = tree

    def __initTreeContents(self):
        self.mass1["Central"][0] = -999.
        self.mass2["Central"][0] = -999.
        for SIG in self.sigStrings:
            self.scoreX[f"{SIG}_Central"][0] = -999.
            self.scoreY[f"{SIG}_Central"][0] = -999.
            self.scoreZ[f"{SIG}_Central"][0] = -999.
        self.fold["Central"][0] = -999
        self.weight["Central"][0] = -999.
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets

        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs[2], 10., 2.4)
        looseMuons = self.SelectMuons(vetoMuons, self.MuonIDs[1], 10., 2.4)
        tightMuons = self.SelectMuons(looseMuons, self.MuonIDs[0], 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs[2], 10., 2.5)
        looseElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs[1], 10., 2.5)
        tightElectrons = self.SelectElectrons(looseElectrons, self.ElectronIDs[0], 10., 2.5)
        jets = self.SelectJets(allJets, "tight", 20., 2.4)
        jets = self.JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = self.mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)

        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        looseMuons = vector[Muon](sorted(looseMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        looseElectrons = vector[Electron](sorted(looseElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))

        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)

    def selectEvent(self, event, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is1E2Mu = (looseMuons.size() == 2 and vetoMuons.size() == 2 and \
                   looseElectrons.size() == 1 and vetoElectrons.size() == 1)
        
        if self.skim == "Skim1E2Mu":
            if not is1E2Mu: return
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return
        if self.skim == "Skim3Mu":
            if not is3Mu: return 
            if tightMuons.size() == looseMuons.size(): return 
            
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets, at least one b-jet
        if self.skim == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return 
            leptons = vector[Lepton]()
            for mu in looseMuons: leptons.emplace_back(mu)
            for ele in looseElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return 
            if not mu1.Charge()+mu2.Charge() == 0: return 
            pair = self.makePair(looseMuons)
            if not pair.M() > 12.: return 
            if not jets.size() >= 2: return 
            if not bjets.size() >= 1: return 
            return "SR1E2Mu"
        
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets, at least one b-jet
        if self.skim == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return 
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return 
            if not mu2.Pt() > 10.: return 
            if not mu3.Pt() > 10.: return 
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return 
            pair1, pair2 = self.makePair(looseMuons)
            if not pair1.M() > 12.: return 
            if not pair2.M() > 12.: return 
            if not jets.size() >= 2: return 
            if not bjets.size() >= 1: return 
            return "SR3Mu"
        
    def makePair(self, muons):
        if muons.size() == 2:
            return (muons[0] + muons[1])
        elif muons.size() == 3:
            mu1, mu2, mu3 = tuple(muons)
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            return (pair1, pair2)
        else:
            raise NotImplementedError(f"Wrong number of muons (muons.size())")
        
    #### Get scores for each event
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        data, fold = getGraphInput(muons, electrons, jets, bjets, METv, self.DataEra)
        for sig, bkg in product(self.sigStrings, self.bkgStrings):
            scores[f"{sig}_vs_{bkg}"] = getGraphScore(self.models[f"{sig}_vs_{bkg}-fold{fold}"], data)
        return data, scores, fold

    def WriteHist(self):
        self.outfile.cd()
        self.tree["Central"].Write()
