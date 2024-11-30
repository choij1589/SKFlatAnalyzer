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

class PromptSkimmer(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        self.initializeAnalyzer()
        
        ## channels
        if self.Skim1E2Mu: self.skim = "Skim1E2Mu"
        if self.Skim3Mu: self.skim = "Skim3Mu"
        if self.skim not in ["Skim1E2Mu", "Skim3Mu"]:
            raise ValueError("Invalid skim option")
        
        ## ParticleNet Implementation
        self.network = "ParticleNet"
        self.sigStrings = ["MHc-160_MA-85", "MHc-130_MA-90", "MHc-100_MA-95"]
        self.bkgStrings = ["nonprompt", "diboson", "ttZ"]
        self.models = loadParticleNet("Combined__", self.sigStrings, self.bkgStrings, pilot=False)
        
        ## Systematic Sources
        self.systematics = [("Central",)]
        
        # Weight Varitaions
        self.weightVariations = []
        self.weightVariations.append(("L1PrefireUp", "L1PrefireDown"))
        self.weightVariations.append(("PileupReweightUp", "PileupReweightDown"))
        self.weightVariations.append(("MuonIDSFUp", "MuonIDSFDown"))
        self.weightVariations.append(("ElectronIDSFUp", "ElectronIDSFDown"))
        self.weightVariations.append(("TriggerSFUp", "TriggerSFDown"))
        self.weightVariations.append(("HeavyTagUpUnCorr", "HeavyTagDownUnCorr"))
        self.weightVariations.append(("HeavyTagUpCorr", "HeavyTagDownCorr"))
        self.weightVariations.append(("LightTagUpUnCorr", "LightTagDownUnCorr"))
        self.weightVariations.append(("LightTagUpCorr", "LightTagDownCorr"))
        
        # Scale Variations
        self.scaleVariations = []
        self.scaleVariations.append(("JetResUp", "JetResDown"))
        self.scaleVariations.append(("JetEnUp", "JetEnDown"))
        self.scaleVariations.append(("ElectronResUp", "ElectronResDown"))
        self.scaleVariations.append(("ElectronEnUp", "ElectronEnDown"))
        self.scaleVariations.append(("MuonEnUp", "MuonEnDown"))

        # WARNING: Unclustered Energy Varaitions are not integrated in the current version
        # When integrating unclustered energy variation, you should check "Central" METvPt passed to find the fold
        # Since int(METv.Pt)+1 is used to find the fold and should not vary due to the difference scale of METv value
        #self.scaleVariations.append(("UnclusteredEnUp", "UnclusteredEnDown"))
        self.alphas_variations = ["AlpS_down", "AlpS_up", "AlpSfact_down", "AlpSfact_up"]

        if not self.IsDATA:
            self.systematics += self.weightVariations + self.scaleVariations
            
        # For Theory Uncertanties
        self.NPDF = 100
        self.NALPHAS = 2
        self.NALPSFACT = 2
        self.NSCALE = 9
        self.NPSSYST = 4
      
        ## Output Tree
        self.__prepareTTree()
        
    def executeEvent(self):
        if not self.PassMETFilter(): return
        ev = self.GetEvent()
        rawMuons = self.GetAllMuons()
        rawElectrons = self.GetAllElectrons()
        rawJets = self.GetAllJets()
        METv = ev.GetMETVector()
        truth = self.GetGens() if not self.IsDATA else None
        
        ## initialize tree contents
        self.__initTreeContents()
        
        def processEvent(syst, prepare_weight_variations=False):
            ## Define objects
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(
                rawMuons, rawElectrons, rawJets, syst
            )
            channel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if not channel: return
            
            # prepare contents
            pairs = self.makePair(tightMuons)
            _, scores, fold = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            
            if channel == "SR1E2Mu":
                self.mass1[syst][0] = pairs.M()
                self.mass2[syst][0] = -999.
            else:   # SR3Mu
                self.mass1[syst][0] = pairs[0].M()
                self.mass2[syst][0] = pairs[1].M()
                
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_nonprompt"]
                self.scoreY[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_diboson"]
                self.scoreZ[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_ttZ"]
            self.fold[syst][0] = fold
            
            self.weight["Central"][0] = 1.
            if self.IsDATA:
                self.tree[syst].Fill()
                return
            
            # weight settings
            w_norm = self.MCweight() * self.GetKFactor() * ev.GetTriggerLumi("Full")
            w_l1prefire, w_l1prefire_up, w_l1prefire_down = self.L1PrefireWeights(prepare_weight_variations)
            w_pileup, w_pileup_up, w_pileup_down = self.PileupWeights(prepare_weight_variations)    
            sf_muonid, sf_muonid_up, sf_muonid_down = self.MuonIDSF(tightMuons, prepare_weight_variations)
            sf_eleid, sf_eleid_up, sf_eleid_down = self.ElectronIDSF(tightElectrons, prepare_weight_variations)
            sf_trig, sf_trig_up, sf_trig_down = self.TriggerSF(channel, tightElectrons, tightMuons, prepare_weight_variations)
            sf_btag = self.BtaaggingSF(bjets, prepare_weight_variations)
            
            self.weight[syst][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.tree[syst].Fill()
            
            if prepare_weight_variations:
                assert syst == "Central", "Only Central systematic is allowed"
                self.weight["L1PrefireUp"][0] = w_norm * w_l1prefire_up * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
                self.weight["L1PrefireDown"][0] = w_norm * w_l1prefire_down * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
                self.weight["PileupReweightUp"][0] = w_norm * w_l1prefire * w_pileup_up * sf_muonid * sf_eleid * sf_trig * sf_btag
                self.weight["PileupReweightDown"][0] = w_norm * w_l1prefire * w_pileup_down * sf_muonid * sf_eleid * sf_trig * sf_btag
                self.weight["MuonIDSFUp"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid_up * sf_eleid * sf_trig * sf_btag
                self.weight["MuonIDSFDown"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid_down * sf_eleid * sf_trig * sf_btag
                self.weight["ElectronIDSFUp"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid_up * sf_trig * sf_btag
                self.weight["ElectronIDSFDown"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid_down * sf_trig * sf_btag
                self.weight["TriggerSFUp"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_up * sf_btag
                self.weight["TriggerSFDown"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_down * sf_btag
                
                for weight_syst in [var for syst_set in self.weightVariations for var in syst_set]:
                    self.mass1[weight_syst][0] = self.mass1["Central"][0]
                    self.mass2[weight_syst][0] = self.mass2["Central"][0]
                    for SIG in self.sigStrings:
                        self.scoreX[f"{SIG}_{weight_syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                        self.scoreY[f"{SIG}_{weight_syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                        self.scoreZ[f"{SIG}_{weight_syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                    self.tree[weight_syst].Fill()
                
                if not self.RunTheoryUnc: return
                
                # PDF
                for pdf_idx in range(self.NPDF):
                    pdf_syst = f"PDFReweight_{pdf_idx}"
                    self.weight[pdf_syst][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_PDF.at(pdf_idx)
                    self.mass1[pdf_syst][0] = self.mass1["Central"][0]
                    self.mass2[pdf_syst][0] = self.mass2["Central"][0]
                    for SIG in self.sigStrings:
                        self.scoreX[f"{SIG}_{pdf_syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                        self.scoreY[f"{SIG}_{pdf_syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                        self.scoreZ[f"{SIG}_{pdf_syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                    self.tree[pdf_syst].Fill()
                
                # AlphaS
                self.weight["AlpS_down"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_AlphaS.at(0)
                self.weight["AlpS_up"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_AlphaS.at(1)
                self.weight["AlpSfact_down"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_alpsfact.at(0)
                self.weight["AlpSfact_up"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_alpsfact.at(1)
                for alphas_syst in self.alphas_variations:
                    self.mass1[alphas_syst][0] = self.mass1["Central"][0]
                    self.mass2[alphas_syst][0] = self.mass2["Central"][0]
                    for SIG in self.sigStrings:
                        self.scoreX[f"{SIG}_{alphas_syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                        self.scoreY[f"{SIG}_{alphas_syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                        self.scoreZ[f"{SIG}_{alphas_syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                    self.tree[alphas_syst].Fill()
                    
                # Renormalization and Factorization Scale
                for scale_idx in range(self.NSCALE):
                    if scale_idx == 5 or scale_idx == 7: continue
                    scale_syst = f"ScaleVar_{scale_idx}"
                    self.weight[scale_syst][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_Scale.at(scale_idx)
                    self.mass1[scale_syst][0] = self.mass1["Central"][0]
                    self.mass2[scale_syst][0] = self.mass2["Central"][0]
                    for SIG in self.sigStrings:
                        self.scoreX[f"{SIG}_{scale_syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                        self.scoreY[f"{SIG}_{scale_syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                        self.scoreZ[f"{SIG}_{scale_syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                    self.tree[scale_syst].Fill()
                
                # Parton Shower
                for ps_idx in range(self.NPSSYST):
                    ps_syst = f"PSVar_{ps_idx}"
                    self.weight[ps_syst][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag * self.weight_PSSyst.at(ps_idx)
                    self.mass1[ps_syst][0] = self.mass1["Central"][0]
                    self.mass2[ps_syst][0] = self.mass2["Central"][0]
                    for SIG in self.sigStrings:
                        self.scoreX[f"{SIG}_{ps_syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                        self.scoreY[f"{SIG}_{ps_syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                        self.scoreZ[f"{SIG}_{ps_syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                    self.tree[ps_syst].Fill()

        if self.IsDATA:
            processEvent("Central")
        else:
            processEvent("Central", prepare_weight_variations=True)
            for syst_up, syst_down in self.scaleVariations:
                processEvent(syst_up)
                processEvent(syst_down)
        
    def L1PrefireWeights(self, prepare_weight_variations=False):
        w_l1prefire = self.GetPrefireWeight(0)
        w_l1prefire_up, w_l1prefire_down = 1., 1.
        if prepare_weight_variations:
            w_l1prefire_up = self.GetPrefireWeight(1)
            w_l1prefire_down = self.GetPrefireWeight(-1)
        return (w_l1prefire, w_l1prefire_up, w_l1prefire_down)
    
    def PileupWeights(self, prepare_weight_variations=False):
        w_pileup = self.GetPileUpWeight(self.nPileUp, 0)
        w_pileup_up, w_pileup_down = 1., 1.
        if prepare_weight_variations:
            w_pileup_up = self.GetPileUpWeight(self.nPileUp, 1)
            w_pileup_down = self.GetPileUpWeight(self.nPileUp, -1)
        return (w_pileup, w_pileup_up, w_pileup_down)
        
    def MuonIDSF(self, muons, prepare_weight_variations=False): 
        sf_muonid = 1.
        sf_muonid_up = 1.
        sf_muonid_down = 1.
        for mu in muons:
            sf_muonid *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
        if not prepare_weight_variations:
            return (sf_muonid, sf_muonid_up, sf_muonid_down)
            
        for mu in muons:
            sf_muonid_up *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
            sf_muonid_down *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
        return (sf_muonid, sf_muonid_up, sf_muonid_down)
            
    def ElectronIDSF(self, electrons, prepare_weight_variations=False):
        sf_eleid = 1.
        sf_eleid_up = 1.
        sf_eleid_down = 1.
        for ele in electrons:
            sf_eleid *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0) * self.getEleIDSF(ele, 0)
        if not prepare_weight_variations:
            return (sf_eleid, sf_eleid_up, sf_eleid_down)
         
        for ele in electrons:
            sf_eleid_up *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 1) * self.getEleIDSF(ele, 1)
            sf_eleid_down *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), -1) * self.getEleIDSF(ele, -1) 
        return (sf_eleid, sf_eleid_up, sf_eleid_down)
            
    def TriggerSF(self, channel, electrons, muons, prepare_weight_variations=False):
        sf_trig  = 1.
        sf_trig_up = 1.
        sf_trig_down = 1.
        if channel == "SR1E2Mu":
            sf_trig = self.getEMuTriggerSF(electrons, muons, 0)
            if prepare_weight_variations:
                sf_trig_up = self.getEMuTriggerSF(electrons, muons, 1)
                sf_trig_down = self.getEMuTriggerSF(electrons, muons, -1)
        if channel == "SR3Mu":
            sf_trig = self.getDblMuTriggerSF(muons, 0)
            if prepare_weight_variations:
                sf_trig_up = self.getDblMuTriggerSF(muons, 1)
                sf_trig_down = self.getDblMuTriggerSF(muons, -1)
        return (sf_trig, sf_trig_up, sf_trig_down)
        
    def BtaaggingSF(self, jets, prepare_weight_variations=False):
        vjets = vector[Jet]()
        for j in jets: vjets.emplace_back(j)
        jtp = jParameters(3, 1, 0, 1)   # DeepJet, Medium, incl, mujets
        return self.mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            
    
    def __prepareTTree(self):
        self.tree = {}
        self.mass1 = {}
        self.mass2 = {}
        self.scoreX = {}
        self.scoreY = {}
        self.scoreZ = {}
        self.fold = {}
        self.weight = {}
        
        for syst in [var for syst_set in self.systematics for var in syst_set]:
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.sigStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.fold[syst] = array("i", [0]); thisTree.Branch("fold", self.fold[syst], "fold/I")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
            
        if not self.RunTheoryUnc: return
        
        # AlphaS variations
        for syst in self.alphas_variations:
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.sigStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.fold[syst] = array("i", [0]); thisTree.Branch("fold", self.fold[syst], "fold/I")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
        
        # PDF
        for pdf_idx in range(self.NPDF):
            syst = f"PDFReweight_{pdf_idx}"
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.sigStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.fold[syst] = array("i", [0]); thisTree.Branch("fold", self.fold[syst], "fold/I")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
        
        # Renormalization and Factorization Scale
        for scale_idx in range(self.NSCALE):
            if scale_idx == 5 or scale_idx == 7: continue
            syst = f"ScaleVar_{scale_idx}"
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.sigStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.fold[syst] = array("i", [0]); thisTree.Branch("fold", self.fold[syst], "fold/I")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
            
        # Parton Shower
        for ps_idx in range(self.NPSSYST):
            syst = f"PSVar_{ps_idx}"
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.sigStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.fold[syst] = array("i", [0]); thisTree.Branch("fold", self.fold[syst], "fold/I")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
        
    def __initTreeContents(self):
        for syst in [var for syst_set in self.systematics for var in syst_set]:
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.fold[syst][0] = -999
            self.weight[syst][0] = -999.
        
        if not self.RunTheoryUnc: return
        # AlphaS
        for syst in self.alphas_variations:
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.fold[syst][0] = -999
            self.weight[syst][0] = -999.
        
        # PDF
        for pdf_idx in range(self.NPDF):
            syst = f"PDFReweight_{pdf_idx}"
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.fold[syst][0] = -999
            self.weight[syst][0] = -999.
        
        # Renormalization and Factorization Scale
        for scale_idx in range(self.NSCALE):
            if scale_idx == 5 or scale_idx == 7: continue
            syst = f"ScaleVar_{scale_idx}"
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.fold[syst][0] = -999
            self.weight[syst][0] = -999.
        
        # Parton Shower
        for ps_idx in range(self.NPSSYST):
            syst = f"PSVar_{ps_idx}"
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.sigStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.fold[syst][0] = -999
            self.weight[syst][0] = -999.
            
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets

        # systematics assertion
        assert syst in ["Central"] + [syst for syst_set in self.scaleVariations for syst in syst_set], f"Wrong systematic {syst}"
        if syst == "MuonEnUp":         allMuons = self.ScaleMuons(allMuons, +1)
        if syst == "MuonEnDown":       allMuons = self.ScaleMuons(allMuons, -1)
        if syst == "ElectronResUp":    allElectrons = self.SmearElectrons(allElectrons, +1)
        if syst == "ElectronsResDown": allElectrons = self.SmearElectrons(allElectrons, -1)
        if syst == "ElectronEnUp":     allElectrons = self.ScaleElectrons(allElectrons, +1)
        if syst == "ElectronEnDown":   allElectrons = self.ScaleElectrons(allElectrons, -1)
        if syst == "JetResUp":         allJets = self.SmearJets(allJets, +1)
        if syst == "JetResDown":       allJets = self.SmearJets(allJets, -1)
        if syst == "JetEnUp":          allJets = self.ScaleJets(allJets, +1)
        if syst == "JetEnDown":        allJets = self.ScaleJets(allJets, -1)

        vetoMuons = self.SelectMuons(allMuons, self.MuonIDs[2], 10., 2.4)
        tightMuons = self.SelectMuons(vetoMuons, self.MuonIDs[0], 10., 2.4)
        vetoElectrons = self.SelectElectrons(allElectrons, self.ElectronIDs[2], 15., 2.5)
        tightElectrons = self.SelectElectrons(vetoElectrons, self.ElectronIDs[0], 15., 2.5)
        jets = self.SelectJets(allJets, "tight", 20., 2.4)
        jets = self.JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = self.mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)

        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))

        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        # check lepton multiplicity first
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
                  
                  
        if self.skim == "Skim1E2Mu" and not is1E2Mu: return 
        if self.skim == "Skim3Mu" and not is3Mu: return 
        
        # for conversion samples
        if self.MCSample in ["DYJets_MG", "DYJets10to50_MG", "TTG", "WWG"]:
            # at least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = vector[Muon]()
            fakeMuons = vector[Muon]()
            convElectrons = vector[Electron]()
            for mu in tightMuons:
                if self.GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
                if self.GetLeptonType(mu, truth) in [-1, -2, -3, -4]: fakeMuons.emplace_back(mu)
            for ele in tightElectrons:
                if self.GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
            if self.skim == "Skim1E2Mu":
                if not fakeMuons.size() == 0: return
                if not convElectrons.size() == 1: return
            if self.skim == "Skim3Mu":
                if not fakeMuons.size() == 0: return
                if not convMuons.size() == 1: return
                
        # Patching Conversion Sample
        leptons = vector[Lepton]()
        for mu in tightMuons: leptons.emplace_back(mu)
        for ele in tightElectrons: leptons.emplace_back(ele)
        region = ""
        for lep in leptons:
            if lep.Pt() < 15.: region = "LowPT"
        if region == "": region = "HighPT"
        if "DYJets" in self.MCSample and not region == "LowPT": return
        if "ZGToLLG" in self.MCSample and not region == "HighPT": return

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets, at least one b-jet
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
            pair = self.makePair(tightMuons)
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
            if not event.PassTrigger(self.DblMuTriggers): return 
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return 
            if not mu2.Pt() > 10.: return 
            if not mu3.Pt() > 10.: return 
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return 
            pair1, pair2 = self.makePair(tightMuons)
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
        for syst in [var for syst_set in self.systematics for var in syst_set]:
            self.tree[syst].Write()
        
        if not self.RunTheoryUnc: return
        # AlphaS
        for syst in self.alphas_variations:
            self.tree[syst].Write()
            
        # PDF
        for pdf_idx in range(self.NPDF):
            syst = f"PDFReweight_{pdf_idx}"
            self.tree[syst].Write()
        
        # Renormalization and Factorization Scale
        for scale_idx in range(self.NSCALE):
            if scale_idx == 5 or scale_idx == 7: continue
            syst = f"ScaleVar_{scale_idx}"
            self.tree[syst].Write()
        
        # Parton Shower
        for ps_idx in range(self.NPSSYST):
            syst = f"PSVar_{ps_idx}"
            self.tree[syst].Write()
            
        
            
