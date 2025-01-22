from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Gen

## Analyzer class for acceptance study
class GenKinematics(TriLeptonBase):
    def __init__(self):
        super().__init__()

    def initializePyAnalyzer(self):
        self.initializeAnalyzer()

    def executeEvent(self):
        # Only Gen level objects are needed
        ev = self.GetEvent()
        weight = self.getEvtWeight(ev)
        genObjects = self.getGenObjects(self.GetGens())
        channel = self.selectEvent(genObjects)
        if channel is None: return
        self.fillGenObjects(channel, genObjects, weight)

    def getEvtWeight(self, ev):
        weight = self.MCweight()*ev.GetTriggerLumi("Full")
        return weight
    
    def getGenObjects(self, truth: vector[Gen]):
        muons = vector[Gen]()
        electrons = vector[Gen]()
        for gen in truth:
            if gen.Status() != 1: continue
            if gen.MotherIndex() < 0: continue #ISR?
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
        b_fromTop = vector[Gen]()
        for gen in truth:
            mother = truth[gen.MotherIndex()]
            if abs(mother.PID()) == 6 and abs(gen.PID()) == 5:
                b_fromTop.emplace_back(gen)

        b_withW = vector[Gen]()
        b_withHc = vector[Gen]()
        
        # b status 23, mother is top
        for b in b_fromTop:
            motherIdx = b.MotherIndex()
            # find if any charged Higgs or W daughter from the same mother
            for gen in truth:
                if gen is b: continue
                if not gen.MotherIndex() == motherIdx: continue
                if gen.Status() != 22: continue # only check resonance decay
                
                if abs(gen.PID()) == 37: # charged Higgs
                    b_withHc.emplace_back(b)
                    break
                elif abs(gen.PID()) == 24: # W boson
                    b_withW.emplace_back(b)
                    break
        
        return {"muons": muons,
                "muons_fromA": muons_fromA,
                "muons_fromEW": muons_fromEW,
                "muons_fromOffshellW": muons_fromOffshellW,
                "muons_fromOther": muons_fromOther,
                "electrons": electrons,
                "electrons_fromEW": electrons_fromEW,
                "electrons_fromOffshellW": electrons_fromOffshellW,
                "electrons_fromOther": electrons_fromOther,
                "b_fromTop": b_fromTop,
                "b_withW": b_withW,
                "b_withHc": b_withHc}

    def selectEvent(self, genObjects):
        muons_fromA = genObjects["muons_fromA"]
        muons_fromEW = genObjects["muons_fromEW"]
        muons_fromOffshellW = genObjects["muons_fromOffshellW"]
        electrons_fromEW = genObjects["electrons_fromEW"]
        electrons_fromOffshellW = genObjects["electrons_fromOffshellW"]
        b_withW = genObjects["b_withW"]
        b_withHc = genObjects["b_withHc"]

        nPromptMuons = muons_fromEW.size() + muons_fromOffshellW.size() + muons_fromA.size()
        nPromptElectrons = electrons_fromEW.size() + electrons_fromOffshellW.size()
        nPromptB = b_withW.size() + b_withHc.size()

        if nPromptMuons == 2 and nPromptElectrons == 1 and nPromptB == 2 and muons_fromA.size() == 2:
            return "GEN1E2Mu"
        elif nPromptMuons == 3 and nPromptElectrons == 0 and nPromptB == 2 and muons_fromA.size() == 2:
            return "GEN3Mu"
        else:
            return None

    def fillGenObjects(self, channel, genObjects, weight):
        muons_fromA = genObjects["muons_fromA"]
        muons_fromEW = genObjects["muons_fromEW"]
        muons_fromOffshellW = genObjects["muons_fromOffshellW"]
        muons_fromOther = genObjects["muons_fromOther"]
        electrons_fromEW = genObjects["electrons_fromEW"]
        electrons_fromOffshellW = genObjects["electrons_fromOffshellW"]
        electrons_fromOther = genObjects["electrons_fromOther"]
        b_withW = genObjects["b_withW"]
        b_withHc = genObjects["b_withHc"]
        
        muons_fromA = vector[Gen](sorted(muons_fromA, key=lambda x: x.Pt(), reverse=True))
        self.FillHist(f"{channel}/muons_fromA/1/pt", muons_fromA.at(0).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/muons_fromA/1/eta", muons_fromA.at(0).Eta(), weight, 100, -5., 5.)
        self.FillHist(f"{channel}/muons_fromA/1/phi", muons_fromA.at(0).Phi(), weight, 64, -3.2, 3.2)
        self.FillHist(f"{channel}/muons_fromA/2/pt", muons_fromA.at(1).Pt(), weight, 300, 0., 300.)
        self.FillHist(f"{channel}/muons_fromA/2/eta", muons_fromA.at(1).Eta(), weight, 100, -5., 5.)
        self.FillHist(f"{channel}/muons_fromA/2/phi", muons_fromA.at(1).Phi(), weight, 64, -3.2, 3.2)

        for mu in muons_fromEW:
            self.FillHist(f"{channel}/muons_fromEW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/muons_fromEW/eta", mu.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/muons_fromEW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOffshellW:
            self.FillHist(f"{channel}/muons_fromOffshellW/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/muons_fromOffshellW/eta", mu.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/muons_fromOffshellW/phi", mu.Phi(), weight, 64, -3.2, 3.2)
        for mu in muons_fromOther:
            self.FillHist(f"{channel}/muons_fromOther/pt", mu.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/muons_fromOther/eta", mu.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/muons_fromOther/phi", mu.Phi(), weight, 64, -3.2, 3.2)

        for ele in electrons_fromEW:
            self.FillHist(f"{channel}/electrons_fromEW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/electrons_fromEW/eta", ele.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/electrons_fromEW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOffshellW:
            self.FillHist(f"{channel}/electrons_fromOffshellW/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/electrons_fromOffshellW/eta", ele.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/electrons_fromOffshellW/phi", ele.Phi(), weight, 64, -3.2, 3.2)
        for ele in electrons_fromOther:
            self.FillHist(f"{channel}/electrons_fromOther/pt", ele.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/electrons_fromOther/eta", ele.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/electrons_fromOther/phi", ele.Phi(), weight, 64, -3.2, 3.2)

        for b in b_withW:
            self.FillHist(f"{channel}/b_withW/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/b_withW/eta", b.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/b_withW/phi", b.Phi(), weight, 64, -3.2, 3.2)
        for b in b_withHc:
            self.FillHist(f"{channel}/b_withHc/pt", b.Pt(), weight, 300, 0., 300.)
            self.FillHist(f"{channel}/b_withHc/eta", b.Eta(), weight, 100, -5., 5.)
            self.FillHist(f"{channel}/b_withHc/phi", b.Phi(), weight, 64, -3.2, 3.2)


