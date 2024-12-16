from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet, Gen

## Analyzer class for acceptance study
class AcceptanceStudy(TriLeptonBase):
    def __init__(self):
        super().__init__()

    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
        if self.Skim1E2Mu: self.channel = "Skim1E2Mu"
        elif self.Skim3Mu: self.channel = "Skim3Mu"
        else:
            print("Wrong channel")
            exit(1)

    def executeEvent(self):
        # Define objects
        ev = self.GetEvent()
        rawMuons = self.GetAllMuons()
        rawElectrons = self.GetAllElectrons()
        rawJets = self.GetAllJets()
        METv = ev.GetMETVector()
        objects = self.defineObjects(rawMuons, rawElectrons, rawJets)
        objects["METv"] = METv
        
        # set weights
        weight = self.getEvtWeight(ev)
        channel = self.selectEvent(ev, objects, weight)
        if channel is None: return None
        self.fillObjects(channel, objects, weight)

        genObjects = self.getGenObjects(self.GetGens())
        self.fillGenObjects(channel, genObjects, weight)

        recoObjects = self.getRecoObjects(objects)
        self.fillRecoObjects(channel, recoObjects, weight)

    def defineObjects(self, rawMuons, rawElectrons, rawJets):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets

        vetoMuons = super().SelectMuons(allMuons, super().MuonIDs[2], 5., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 5., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 5., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 5., 2.5)
        jets = super().SelectJets(allJets, "tight", 15., 2.4)
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

        # return objects as a dictionary
        return {"vetoMuons": vetoMuons,
                "tightMuons": tightMuons,
                "vetoElectrons": vetoElectrons,
                "tightElectrons": tightElectrons,
                "jets": jets,
                "bjets": bjets}

    def getEvtWeight(self, ev):
        weight = self.MCweight()*ev.GetTriggerLumi("Full")
        weight *= self.GetPrefireWeight(0)
        weight *= self.GetPileUpWeight(self.nPileUp, 0)

        # jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
        # vjets = vector[Jet]()
        # for j in jets: vjets.emplace_back(j)
        # weight *= super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
                
        return weight
    
    def selectEvent(self, ev, objects, weight):
        self.FillHist("cutflow", 0., weight, 10, 0., 10.)

        # MET filter
        if not super().PassMETFilter(): return None
        self.FillHist("cutflow", 1., weight, 10, 0., 10.)

        vetoMuons, tightMuons = objects["vetoMuons"], objects["tightMuons"]
        vetoElectrons, tightElectrons = objects["vetoElectrons"], objects["tightElectrons"]

        # lepton multiplicity
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None

        # prompt matching for fake study
        if self.FakeStudy:
            truth = self.GetGens() if not self.IsDATA else None
            promptMuons = vector[Muon]()
            promptElectrons = vector[Electron]()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)
            isPromptMatched = promptMuons.size() == tightMuons.size() and promptElectrons.size() == tightElectrons.size()
            if isPromptMatched: return None
        self.FillHist("cutflow", 2., weight, 10, 0., 10.)

        if self.channel == "Skim1E2Mu":
            return self._select1E2MuEvent(ev, objects, weight)
        elif self.channel == "Skim3Mu":
            return self._select3MuEvent(ev, objects, weight)
        else:
            raise ValueError(f"[AcceptanceStudy::selectEvent] wrong channel {self.channel}")
    
    def getGenObjects(self, truth: vector[Gen]):
        # find muons and electrons, status 1
        muons = vector[Gen]()
        electrons = vector[Gen]()
        for gen in truth:
            if gen.Status() != 1: continue
            if gen.MotherIndex() < 0: continue
            if abs(gen.PID()) == 11:
                electrons.emplace_back(gen)
            if abs(gen.PID()) == 13:
                muons.emplace_back(gen)
        
        muons_fromA = vector[Gen]()
        muons_fromEW = vector[Gen]()
        muons_fromOffshellW = vector[Gen]()
        muons_fromOther = vector[Gen]()
        for mu in muons:
            leptonType = self.GetLeptonType(mu, truth)
            if leptonType == 1: muons_fromEW.emplace_back(mu)
            elif leptonType == 2: muons_fromA.emplace_back(mu)
            elif leptonType == 6: muons_fromOffshellW.emplace_back(mu)
            else: muons_fromOther.emplace_back(mu)

        electrons_fromEW = vector[Gen]()
        electrons_fromOffshellW = vector[Gen]()
        electrons_fromOther = vector[Gen]()

        for ele in electrons:
            leptonType = self.GetLeptonType(ele, truth)
            if leptonType == 1: electrons_fromEW.emplace_back(ele)
            elif leptonType == 6: electrons_fromOffshellW.emplace_back(ele)
            else: electrons_fromOther.emplace_back(ele)

        # find b parton
        b_justFromTop = vector[Gen]()
        for gen in truth:
            mother = truth[gen.MotherIndex()]
            if abs(mother.PID()) == 6 and abs(gen.PID()) == 5:
                b_justFromTop.emplace_back(gen)

        b_withW = vector[Gen]()
        b_withHc = vector[Gen]()
        
        # b status 23, mother is top
        for b in b_justFromTop:
            motherIdx = b.MotherIndex()
            # find if any charged Higgs daughter from the same mother
            foundHc = False
            for gen in truth:
                if gen is b: continue
                if not gen.MotherIndex() == motherIdx: continue
                if abs(gen.PID()) == 37 and gen.Status() == 22: # resonance decay
                    b_withHc.emplace_back(b)
                    foundHc = True
                    break
            
            if foundHc: continue
            for gen in truth:
                if gen is b: continue
                if not gen.MotherIndex() == motherIdx: continue
                if abs(gen.PID()) == 24 and gen.Status() == 22: # resonance decay
                    b_withW.emplace_back(b)
                    break        

        return {"muons_fromA": muons_fromA,
                "muons_fromEW": muons_fromEW,
                "muons_fromOffshellW": muons_fromOffshellW,
                "muons_fromOther": muons_fromOther,
                "electrons_fromEW": electrons_fromEW,
                "electrons_fromOffshellW": electrons_fromOffshellW,
                "electrons_fromOther": electrons_fromOther,
                "b_withW": b_withW,
                "b_withHc": b_withHc}
    
    def getRecoObjects(self, objects):
        muons, electrons, jets = objects["tightMuons"], objects["tightElectrons"], objects["jets"]
        truth = self.GetGens() if not self.IsDATA else None
        muons_fromA =         vector[Muon]()
        muons_fromEW =        vector[Muon]()
        muons_fromOffshellW = vector[Muon]()
        muons_fromIntConv =   vector[Muon]()
        muons_fromExtConv =   vector[Muon]()
        muons_fromHadron =    vector[Muon]()
        muons_fromOther =     vector[Muon]()

        for mu in muons:
            leptonType = self.GetLeptonType(mu, self.GetGens())
            if leptonType == 1:        muons_fromEW.emplace_back(mu)
            elif leptonType == 2:      muons_fromA.emplace_back(mu)
            elif leptonType == 6:      muons_fromOffshellW.emplace_back(mu)
            elif leptonType in [4, 5]: muons_fromIntConv.emplace_back(mu)
            elif leptonType == -6:     muons_fromExtConv.emplace_back(mu)
            elif leptonType < 0:       muons_fromHadron.emplace_back(mu)
            else:                      muons_fromOther.emplace_back(mu)

        electrons_fromEW = vector[Electron]()
        electrons_fromOffshellW = vector[Electron]()
        electrons_fromIntConv = vector[Electron]()
        electrons_fromExtConv = vector[Electron]()
        electrons_fromHadron = vector[Electron]()
        electrons_fromOther = vector[Electron]()
        for ele in electrons:
            leptonType = self.GetLeptonType(ele, self.GetGens())
            if leptonType == 1: electrons_fromEW.emplace_back(ele)
            elif leptonType == 6: electrons_fromOffshellW.emplace_back(ele)
            elif leptonType in [4, 5]: electrons_fromIntConv.emplace_back(ele)
            elif leptonType == -6: electrons_fromExtConv.emplace_back(ele)
            elif leptonType < 0: electrons_fromHadron.emplace_back(ele)
            else: electrons_fromOther.emplace_back(ele)

        # find b partons
        b_justFromTop = vector[Gen]()
        for gen in truth:
            mother = truth[gen.MotherIndex()]
            if abs(mother.PID()) == 6 and abs(gen.PID()) == 5:
                b_justFromTop.emplace_back(gen)
        
        # find jets matched to b partons
        b_matchedWithHc = vector[Jet]()
        b_matchedWithW = vector[Jet]()
        for jet in jets:
            if not jet.GenHFHadronMatcherFlavour() == 5: continue
            min_dR = 0.3
            matched_parton = None
            for b in b_justFromTop:
                if jet.DeltaR(b) > min_dR: continue
                min_dR = jet.DeltaR(b)
                matched_parton = b
            
            if matched_parton is None: continue
            for gen in truth:
                if gen is matched_parton: continue
                if not gen.MotherIndex() == matched_parton.MotherIndex(): continue
                if abs(gen.PID()) == 37 and gen.Status() == 22: # resonance decay
                    b_matchedWithHc.emplace_back(jet)
                    break
            
            if b_matchedWithHc.size() > 0: continue
            for gen in truth:
                if gen is matched_parton: continue
                if not gen.MotherIndex() == matched_parton.MotherIndex(): continue
                if abs(gen.PID()) == 24 and gen.Status() == 22: # resonance decay
                    b_matchedWithW.emplace_back(jet)
                    break

        return {"muons_fromA": muons_fromA,
                "muons_fromEW": muons_fromEW,
                "muons_fromOffshellW": muons_fromOffshellW,
                "muons_fromIntConv": muons_fromIntConv,
                "muons_fromExtConv": muons_fromExtConv,
                "muons_fromHadron": muons_fromHadron,
                "muons_fromOther": muons_fromOther,
                "electrons_fromEW": electrons_fromEW,
                "electrons_fromOffshellW": electrons_fromOffshellW,
                "electrons_fromIntConv": electrons_fromIntConv,
                "electrons_fromExtConv": electrons_fromExtConv,
                "electrons_fromHadron": electrons_fromHadron,
                "electrons_fromOther": electrons_fromOther,
                "b_matchedWithHc": b_matchedWithHc,
                "b_matchedWithW": b_matchedWithW}
     
    def fillObjects(self, channel, objects, weight):
        muons = objects["tightMuons"]
        electrons = objects["tightElectrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]

        for idx, mu in enumerate(muons, start=1):
            self.FillHist(f"{channel}/RECO/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for idx, ele in enumerate(electrons, start=1):
            self.FillHist(f"{channel}/RECO/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for idx, jet in enumerate(jets, start=1):
            self.FillHist(f"{channel}/RECO/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
        for idx, bjet in enumerate(bjets, start=1):
            self.FillHist(f"{channel}/RECO/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/RECO/muons/size", muons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/RECO/electrons/size", electrons.size(), weight, 10, 0., 10.)
        self.FillHist(f"{channel}/RECO/jets/size", jets.size(), weight, 20, 0., 20.)
        self.FillHist(f"{channel}/RECO/bjets/size", bjets.size(), weight, 15, 0., 15.)
        self.FillHist(f"{channel}/RECO/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/RECO/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)

    def fillGenObjects(self, channel, genObjects, weight):
        muons_fromA, muons_fromEW, muons_fromOffshellW, muons_fromOther = genObjects["muons_fromA"], genObjects["muons_fromEW"], genObjects["muons_fromOffshellW"], genObjects["muons_fromOther"]
        electrons_fromEW, electrons_fromOffshellW, electrons_fromOther = genObjects["electrons_fromEW"], genObjects["electrons_fromOffshellW"], genObjects["electrons_fromOther"]
        b_withW, b_withHc = genObjects["b_withW"], genObjects["b_withHc"]

        try:
            assert muons_fromA.size() == 2
        except:
            print(f"[AcceptanceStudy::fillGenObjects] wrong no. of muons from A {muons_fromA.size()}")
        
        muons_fromA = vector[Gen](sorted(muons_fromA, key=lambda x: x.Pt(), reverse=True))
        self.FillHist(f"{channel}/GEN/muons_fromA/1/pt", muons_fromA.at(0).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/GEN/muons_fromA/1/eta", muons_fromA.at(0).Eta(), weight, 48, -2.4, 2.4)
        self.FillHist(f"{channel}/GEN/muons_fromA/1/phi", muons_fromA.at(0).Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/GEN/muons_fromA/2/pt", muons_fromA.at(1).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/GEN/muons_fromA/2/eta", muons_fromA.at(1).Eta(), weight, 48, -2.4, 2.4)
        self.FillHist(f"{channel}/GEN/muons_fromA/2/phi", muons_fromA.at(1).Phi(), weight, 64, -3.2, 3.2)

        for mu in muons_fromEW:
            self.FillHist(f"{channel}/GEN/muons_fromEW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/muons_fromEW/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/GEN/muons_fromEW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOffshellW:
            self.FillHist(f"{channel}/GEN/muons_fromOffshellW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/muons_fromOffshellW/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/GEN/muons_fromOffshellW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOther:
            self.FillHist(f"{channel}/GEN/muons_fromOther/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/muons_fromOther/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/GEN/muons_fromOther/phi", mu.Phi(), weight, 64, -3.2, 3.2)

        for ele in electrons_fromEW:
            self.FillHist(f"{channel}/GEN/electrons_fromEW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/electrons_fromEW/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/GEN/electrons_fromEW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOffshellW:
            self.FillHist(f"{channel}/GEN/electrons_fromOffshellW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/electrons_fromOffshellW/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/GEN/electrons_fromOffshellW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOther:
            self.FillHist(f"{channel}/GEN/electrons_fromOther/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/GEN/electrons_fromOther/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/GEN/electrons_fromOther/phi", ele.Phi(), weight, 64, -3.2, 3.2)

        for b in b_withW:
            self.FillHist(f"{channel}/LHE/b_withW/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/LHE/b_withW/eta", b.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/LHE/b_withW/phi", b.Phi(), weight, 64, -3.2, 3.2)
        for b in b_withHc:
            self.FillHist(f"{channel}/LHE/b_withHc/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/LHE/b_withHc/eta", b.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/LHE/b_withHc/phi", b.Phi(), weight, 64, -3.2, 3.2)

    def fillRecoObjects(self, channel, recoObjects, weight):
        muons_fromA, muons_fromEW, muons_fromOffshellW, muons_fromIntConv, muons_fromExtConv, muons_fromHadron, muons_fromOther = recoObjects["muons_fromA"], recoObjects["muons_fromEW"], recoObjects["muons_fromOffshellW"], recoObjects["muons_fromIntConv"], recoObjects["muons_fromExtConv"], recoObjects["muons_fromHadron"], recoObjects["muons_fromOther"]
        electrons_fromEW, electrons_fromOffshellW, electrons_fromIntConv, electrons_fromExtConv, electrons_fromHadron, electrons_fromOther = recoObjects["electrons_fromEW"], recoObjects["electrons_fromOffshellW"], recoObjects["electrons_fromIntConv"], recoObjects["electrons_fromExtConv"], recoObjects["electrons_fromHadron"], recoObjects["electrons_fromOther"]
        b_matchedWithHc, b_matchedWithW = recoObjects["b_matchedWithHc"], recoObjects["b_matchedWithW"]
        
        if muons_fromA.size() < 2:
            return
        
        muons_fromA = vector[Muon](sorted(muons_fromA, key=lambda x: x.Pt(), reverse=True))
        self.FillHist(f"{channel}/RECO/muons_fromA/1/pt", muons_fromA.at(0).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/RECO/muons_fromA/1/eta", muons_fromA.at(0).Eta(), weight, 48, -2.4, 2.4)
        self.FillHist(f"{channel}/RECO/muons_fromA/1/phi", muons_fromA.at(0).Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/RECO/muons_fromA/2/pt", muons_fromA.at(1).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/RECO/muons_fromA/2/eta", muons_fromA.at(1).Eta(), weight, 48, -2.4, 2.4)
        self.FillHist(f"{channel}/RECO/muons_fromA/2/phi", muons_fromA.at(1).Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromEW:
            self.FillHist(f"{channel}/RECO/muons_fromEW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromEW/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromEW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOffshellW:
            self.FillHist(f"{channel}/RECO/muons_fromOffshellW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromOffshellW/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromOffshellW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromIntConv:
            self.FillHist(f"{channel}/RECO/muons_fromIntConv/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromIntConv/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromIntConv/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromExtConv:
            self.FillHist(f"{channel}/RECO/muons_fromExtConv/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromExtConv/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromExtConv/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromHadron:
            self.FillHist(f"{channel}/RECO/muons_fromHadron/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromHadron/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromHadron/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOther:
            self.FillHist(f"{channel}/RECO/muons_fromOther/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/muons_fromOther/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/muons_fromOther/phi", mu.Phi(), weight, 64, -3.2, 3.2)

        for ele in electrons_fromEW:
            self.FillHist(f"{channel}/RECO/electrons_fromEW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromEW/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromEW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOffshellW:
            self.FillHist(f"{channel}/RECO/electrons_fromOffshellW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromOffshellW/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromOffshellW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromIntConv:
            self.FillHist(f"{channel}/RECO/electrons_fromIntConv/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromIntConv/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromIntConv/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromExtConv:
            self.FillHist(f"{channel}/RECO/electrons_fromExtConv/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromExtConv/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromExtConv/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromHadron:
            self.FillHist(f"{channel}/RECO/electrons_fromHadron/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromHadron/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromHadron/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOther:
            self.FillHist(f"{channel}/RECO/electrons_fromOther/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/electrons_fromOther/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            self.FillHist(f"{channel}/RECO/electrons_fromOther/phi", ele.Phi(), weight, 64, -3.2, 3.2)

        for idx, b in enumerate(b_matchedWithHc, start=1):
            self.FillHist(f"{channel}/RECO/b_matchedWithHc/{idx}/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/b_matchedWithHc/{idx}/eta", b.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/b_matchedWithHc/{idx}/phi", b.Phi(), weight, 64, -3.2, 3.2)
        for idx, b in enumerate(b_matchedWithW, start=1):
            self.FillHist(f"{channel}/RECO/b_matchedWithW/{idx}/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/RECO/b_matchedWithW/{idx}/eta", b.Eta(), weight, 48, -2.4, 2.4)
            self.FillHist(f"{channel}/RECO/b_matchedWithW/{idx}/phi", b.Phi(), weight, 64, -3.2, 3.2)
    
    def _select1E2MuEvent(self, ev, objects, weight):
        tightMuons, tightElectrons = objects["tightMuons"], objects["tightElectrons"]
        jets, bjets = objects["jets"], objects["bjets"]
        leptons = vector[Lepton]()
        for mu in tightMuons: leptons.emplace_back(mu)
        for ele in tightElectrons: leptons.emplace_back(ele)
        mu1, mu2, ele = tuple(leptons)
        
        # Trigger and kinematic cuts
        if self.CheckTrigger:
            if not ev.PassTrigger(self.EMuTriggers): return None
            self.FillHist("cutflow", 3., weight, 10, 0., 10.)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            if not (passLeadMu or passLeadEle): return None
            self.FillHist("cutflow", 4., weight, 10, 0., 10.)
        
        # charge condition
        if not mu1.Charge()+mu2.Charge() == 0: return None
        self.FillHist("cutflow", 5., weight, 10, 0., 10.)
        
        # mass condition
        pair = self._makePair(tightMuons)
        if not pair.M() > 12.: return None
        self.FillHist("cutflow", 6., weight, 10, 0., 10.)
        
        # jet multiplicity
        if not jets.size() >= 2: return None
        self.FillHist("cutflow", 7., weight, 10, 0., 10.)
        if not bjets.size() >= 1: return None
        self.FillHist("cutflow", 8., weight, 10, 0., 10.)
        
        return "SR1E2Mu"

    def _select3MuEvent(self, ev, objects, weight):

        tightMuons = objects["tightMuons"]
        jets, bjets = objects["jets"], objects["bjets"]
        mu1, mu2, mu3 = tuple(tightMuons)

        if self.CheckTrigger:
            if not ev.PassTrigger(self.DblMuTriggers): return None
            self.FillHist("cutflow", 3., weight, 10, 0., 10.)
            if not (mu1.Pt() > 20. and mu2.Pt() > 10. and mu3.Pt() > 10.): return None
            self.FillHist("cutflow", 4., weight, 10, 0., 10.)

        # charge condition
        if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
        self.FillHist("cutflow", 5., weight, 10, 0., 10.)

        # mass condition
        pair1, pair2 = self._makePair(tightMuons)
        if not (pair1.M() > 12. and pair2.M() > 12.): return None
        self.FillHist("cutflow", 6., weight, 10, 0., 10.)

        # jet multiplicity
        if not jets.size() >= 2: return None
        self.FillHist("cutflow", 7., weight, 10, 0., 10.)
        if not bjets.size() >= 1: return None
        self.FillHist("cutflow", 8., weight, 10, 0., 10.)

        return "SR3Mu"
    
    def _makePair(self, muons):
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
            raise ValueError(f"[AcceptanceStudy::makePair] wrong no. of muons {muons.size}")
