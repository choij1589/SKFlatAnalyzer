from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet

import numpy as np
from itertools import product
from MLTools.helpers import loadParticleNet, loadGBDTClassifier
from MLTools.helpers import getGraphInput, getGraphScore

class MatrixSelector(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
       
        # Flags
        self.skim = ""
        if self.Skim1E2Mu: self.skim = "Skim1E2Mu"
        if self.Skim3Mu: self.skim = "Skim3Mu"
        if not self.skim in ["Skim1E2Mu", "Skim3Mu"]:
            raise ValueError(f"Invalid skim flag: {self.skim}")

        self.sigStrings = ["MHc-100_MA-95", "MHc-130_MA-90", "MHc-160_MA-85"]
        self.bkgStrings = ["nonprompt", "diboson", "ttZ"]
        self.models = loadParticleNet("Combined__", self.sigStrings, self.bkgStrings, pilot=False)
        self.classifiers = loadGBDTClassifier(self.DataEra, self.skim.replace("Skim", "SR"), self.sigStrings)
        
    def executeEvent(self):
        if not self.PassMETFilter(): return
        ev = self.GetEvent()
        rawMuons = self.GetAllMuons()
        rawElectrons = self.GetAllElectrons()
        rawJets = self.GetAllJets()
        METv = ev.GetMETVector()
        truth = None
        
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)
        if not channel: return
        
        data, scores, fold = self.evalScore(looseMuons, looseElectrons, jets, bjets, METv)
        objects = {"muons": looseMuons,
                   "electrons": looseElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "data": data,
                   "scores": scores,
                   "fold": fold
        }
        weight = self.getFakeWeight(looseMuons, looseElectrons, syst="Central")
        self.fillObjects(channel, objects, weight, syst="Central")
    
    def defineObjects(self, rawMuons, rawElectrons, rawJets):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs[2], 10., 2.4)
        looseMuons = self.SelectMuons(vetoMuons, self.MuonIDs[1], 10., 2.4)
        tightMuons = self.SelectMuons(looseMuons, self.MuonIDs[0], 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs[2], 15., 2.5)
        looseElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs[1], 15., 2.5)
        tightElectrons = self.SelectElectrons(looseElectrons, self.ElectronIDs[0], 15., 2.5)
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
    
    def selectEvent(self, event, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is1E2Mu = (looseMuons.size() == 2 and vetoMuons.size() == 2 and \
                   looseElectrons.size() == 1 and vetoElectrons.size() == 1)
        
        #### not all leptons tight
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
        ## 4. At least two jets
        if self.skim == "Skim1E2Mu":
            if not event.PassTrigger(self.EMuTriggers): return
            leptons = vector[Lepton]()
            for mu in looseMuons: leptons.emplace_back(mu)
            for ele in looseElectrons: leptons.emplace_back(ele)
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
            if not event.PassTrigger(self.DblMuTriggers): return
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return
            if not mu2.Pt() > 10.: return
            if not mu3.Pt() > 10.: return 
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return 
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(looseMuons)
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
            raise EOFError()
        
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        data, fold = getGraphInput(muons, electrons, jets, bjets, METv, self.DataEra)
        for sig, bkg in product(self.sigStrings, self.bkgStrings):
            scores[f"{sig}_vs_{bkg}"] = getGraphScore(self.models[f"{sig}_vs_{bkg}-fold{fold}"], data)
        
        for sig in self.sigStrings:
            scoreX = scores[f"{sig}_vs_nonprompt"]
            scoreY = scores[f"{sig}_vs_diboson"]
            scoreZ = scores[f"{sig}_vs_ttZ"]
            score = self.classifiers[sig][fold].predict_proba(np.array([[scoreX, scoreY, scoreZ]]))[0][1]
            scores[sig] = score
        
        return data, scores, fold
    
    def fillObjects(self, channel, objects, weight, syst="Central"):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        data = objects["data"]
        scores = objects["scores"]
        
        ## fill base observables
        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/energy", mu.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/px", mu.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/py", mu.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/pz", mu.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/muons/{idx}/charge", mu.Charge(), weight, 3, -1, 2)
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/energy", ele.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/px", ele.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/py", ele.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/pz", ele.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/electrons/{idx}/charge", ele.Charge(), weight, 3, -1, 2)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/energy", jet.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/px", jet.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/py", jet.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/pz", jet.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/charge", jet.Charge(), weight, 200, -1, 1)
            self.FillHist(f"{channel}/{syst}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/energy", bjet.E(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/px", bjet.Px(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/py", bjet.Py(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/pz", bjet.Pz(), weight, 500, -250., 250.)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/charge", bjet.Charge(), weight, 200, -1, 1)
            self.FillHist(f"{channel}/{syst}/bjets/{idx}/btagScore", bjet.GetTaggerResult(3), weight, 100, 0., 1.)
        self.FillHist(f"{channel}/{syst}/muons/size", muons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/electrons/size", electrons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        self.FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        self.FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/{syst}/METv/energy", METv.E(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/{syst}/METv/px", METv.Px(), weight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/py", METv.Py(), weight, 500, -250., 250.)
        self.FillHist(f"{channel}/{syst}/METv/pz", METv.Pz(), weight, 500, -250., 250.)
        
        # Fill discrimination variable
        if "1E2Mu" in channel:
            pair = muons.at(0) + muons.at(1)
            nonprompt = electrons.at(0)
            self.FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 200, 0., 200.)
        else:
            mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
            pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
            self.FillHist(f"{channel}/{syst}/stack/pt", pair1.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/stack/eta", pair1.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/stack/phi", pair1.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/stack/mass", pair1.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/stack/pt", pair2.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/stack/eta", pair2.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/stack/phi", pair2.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/stack/mass", pair2.M(), weight, 200, 0., 200.)
        
        # Fill ZCands
        if "1E2Mu" in channel:
            ZCand = muons.at(0) + muons.at(1)
            nonprompt = electrons.at(0)
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
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
            self.FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            self.FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            self.FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/pt", nonprompt.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/{syst}/nonprompt/eta", nonprompt.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/{syst}/nonprompt/phi", nonprompt.Phi(), weight, 64, -3.2, 3.2)
         
        # Fill signal dependent distributions 
        for signal in self.sigStrings:
            if "1E2Mu" in channel:
                ACand = muons.at(0) + muons.at(1)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
            else:
                mu_ss1, mu_ss2, mu_os = self.configureChargeOf(muons)
                pair1, pair2 = (mu_ss1+mu_os), (mu_ss2+mu_os)
                mA = int(signal.split("_")[1].split("-")[1]) 
                if abs(pair1.M() - mA) < abs(pair2.M() - mA): ACand, nACand = pair1, pair2
                else:                                         ACand, nACand = pair2, pair2
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                self.FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
                self.FillHist(f"{channel}/{syst}/{signal}/nACand/pt", nACand.Pt(), weight, 300, 0., 300.)
                self.FillHist(f"{channel}/{syst}/{signal}/nACand/eta", nACand.Eta(), weight, 100, -5., 5.)
                self.FillHist(f"{channel}/{syst}/{signal}/nACand/phi", nACand.Phi(), weight, 64, -3.2, 3.2)
                self.FillHist(f"{channel}/{syst}/{signal}/nACand/mass", nACand.M(), weight, 300, 0., 300.)
            
            score_nonprompt = scores[f"{signal}_vs_nonprompt"]
            score_diboson   = scores[f"{signal}_vs_diboson"]
            score_ttZ       = scores[f"{signal}_vs_ttZ"]
            score           = scores[signal]
            self.FillHist(f"{channel}/{syst}/{signal}/score_nonprompt", score_nonprompt, weight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/{signal}/score_diboson", score_diboson, weight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/{signal}/score_ttZ", score_ttZ, weight, 100, 0., 1.)
            self.FillHist(f"{channel}/{syst}/{signal}/score", score, weight, 100, 0., 1.)