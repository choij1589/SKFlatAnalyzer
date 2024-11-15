from ROOT import gSystem
from ROOT import TutorialBase

class TutorialRun(TutorialBase):
    def __init__(self):
        super().__init__()

    def initializePyAnalyzer(self):
        super().initializeAnalyzer()

    # Override executeEvent function
    def executeEvent(self):
        allMuons = super().GetAllMuons()
        allJets = super().GetAllJets()

        for muonID in self.MuonIDs:
            super().FillHist(f"{muonID}/Cutflow", 0., 1., 10, 0., 10.)
            
            # MET Filter
            if not super().PassMETFilter(): return None
            ev = super().GetEvent()
            METv = ev.GetMETVector()

            # Trigger
            if not ev.PassTrigger(super().IsoMuTriggerName): return None
            super().FillHist(f"{muonID}/Cutflow", 2., 1., 10, 0., 10.)
            muons = list(super().SelectMuons(allMuons, muonID, 20., 2.4))
            jets = list(super().SelectJets(allJets, "tight", 30., 2.4))
            # Sort in pt order
            muons.sort(key=lambda x: x.Pt(), reverse=True)
            jets.sort(key=lambda x: x.Pt(), reverse=True)

            # b-tagging
            bjets = []
            for jet in jets:
                btagScore = jet.GetTaggerResult(3) # which is DeepJet
                wp = super().mcCorr.GetJetTaggingCutValue(3, 1) # DeepJet Medium
                if (btagScore > wp): bjets.append(jet)
            
            # Event selection
            if not (len(muons) == 2): return None
            if not (muons[0].Charge() + muons[1].Charge() == 0): return None
            if not (muons[0].Pt() > super().TriggerSafePtCut): return None
            ZCand = muons[0]+muons[1]
            MZ = 91.2
            if not (abs(ZCand.M() - MZ) < 15.): return None

            weight = 1.
            if not super().IsDATA:
                genWeight = super().MCweight()
                lumi = ev.GetTriggerLumi("Full")
                weight *= (genWeight * lumi)

                # Muon ID efficiency
                #for mu in muons:
                #    print(super().mcCorr.MuonID_SF("NUM_TopHN_DEN_TrackerMuons", mu.Eta(), mu.Pt()))

            super().FillHist(f"{muonID}/ZCand/mass", ZCand.M(), weight, 40, 70., 110.);
            super().FillHist(f"{muonID}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
            super().FillHist(f"{muonID}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
            super().FillHist(f"{muonID}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
            super().FillHist(f"{muonID}/MET", METv.Pt(), weight, 300, 0., 300.);
            super().FillHist(f"{muonID}/bjets/size", len(bjets), weight, 20, 0., 20.);
