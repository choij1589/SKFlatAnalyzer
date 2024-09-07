from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet

from itertools import product
from MLTools.helpers import loadModels
from MLTools.helpers import getGraphInput, getGraphScore

class TriLeptons(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        # Flags
        self.skim = ""
        self.run_syst = False
        if super().Skim1E2Mu: self.skim = "Skim1E2Mu"
        if super().Skim3Mu:   self.skim = "Skim3Mu"
        if super().RunSyst:   self.run_syst = True
        if not self.skim in ["Skim1E2Mu", "Skim3Mu"]:
            raise ValueError(f"Invalid skim {self.skim}")
        
        # Systematics
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        if self.run_syst:
            if self.skim == "Skim1E2Mu":
                # Add PDF!
                self.weightVariations.append(("L1PrefireUp", "L1PrefireDown"))
                self.weightVariations.append(("PileupReweightUp", "PileupReweightDown"))
                self.weightVariations.append(("MuonIDSFUp", "MuonIDSFDown"))
                self.weightVariations.append(("ElectronIDSFUp", "ElectronIDSFDown"))
                self.weightVariations.append(("EMuTrigSFUp", "EMuTrigSFDown"))
                self.weightVariations.append(("HeavyTagUpUnCorr", "HeavyTagDownUnCorr"))
                self.weightVariations.append(("HeavyTagUpCorr", "HeavyTagDownCorr"))
                self.weightVariations.append(("LightTagUpUnCorr", "LightTagDownUnCorr"))
                self.weightVariations.append(("LightTagUpCorr", "LightTagDownCorr")) 
            if self.skim == "Skim3Mu":
                self.weightVariations.append(("L1PrefireUp", "L1PrefireDown"))
                self.weightVariations.append(("PileupReweightUp", "PileupReweightDown"))
                self.weightVariations.append(("MuonIDSFUp", "MuonIDSFDown"))
                self.weightVariations.append(("DblMuTrigSFUp", "DblMuTrigSFDown"))
                self.weightVariations.append(("HeavyTagUpUnCorr", "HeavyTagDownUnCorr"))
                self.weightVariations.append(("HeavyTagUpCorr", "HeavyTagDownCorr"))
                self.weightVariations.append(("LightTagUpUnCorr", "LightTagDownUnCorr"))
                self.weightVariations.append(("LightTagUpCorr", "LightTagDownCorr"))
            self.scaleVariations.append(("JetResUp", "JetResDown"))
            self.scaleVariations.append(("JetEnUp", "JetEnDown"))
            self.scaleVariations.append(("ElectronResUp", "ElectronResDown"))
            self.scaleVariations.append(("ElectronEnUp", "ElectronEnDown"))
            self.scaleVariations.append(("MuonEnUp", "MuonEnDown"))
        self.systematics = self.weightVariations + self.scaleVariations
            
        # Load ML models
        self.signalStrings = ["MHc-100_MA-95", "MHc-130_MA-90", "MHc-160_MA-85"]
        self.backgroundStrings = ["nonprompt", "diboson", "ttZ"]
        self.models = loadModels("GraphNeuralNet", self.channel, self.signalStrings, self.backgroundStrings)
    
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central Scale
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
        if not channel is None:
            data, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            # make objects as a dictionary
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       "data": data,
                       "scores": scores}
            for syst in self.weightVariations:
                weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, syst)
                self.FillObjects(channel, objects, weight, syst)
        
        # Scale Variations
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if not channel is None:
                data, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
                objects = {"muons": tightMuons,
                           "electrons": tightElectrons,
                           "jets": jets,
                           "bjets": bjets,
                           "METv": METv,
                           "data": data,
                           "scores": scores
                           }
                weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, syst)
                self.FillObjects(channel, objects, weight, syst)