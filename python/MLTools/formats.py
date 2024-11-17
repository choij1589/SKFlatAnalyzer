from ROOT import TLorentzVector
from ROOT import TMath

class NodeParticle(TLorentzVector):
    def __init__(self):
        TLorentzVector.__init__(self)
        self.charge = 0
        self.btagScore = 0.
        self.isMuon = False
        self.isElectron = False
        self.isJet = False

    def Charge(self):
        return self.charge

    def BtagScore(self):
        return self.btagScore

    def MT(self, part):
        dPhi = self.DeltaPhi(part)
        return TMath.Sqrt(2*self.Pt()*part.Pt()*(1-TMath.Cos(dPhi)))

    def IsMuon(self):
        return self.isMuon

    def IsElectron(self):
        return self.isElectron

    def IsJet(self):
        return self.isJet

