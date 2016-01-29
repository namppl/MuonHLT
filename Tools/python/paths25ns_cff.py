class HltPathConfig :
    """Configuration class for paths of which efficiency has to be studied"""

    def __init__(self) :    
        self.TRIGNAME = "DUMMY"
        self.PROBEL1 = "" 
        self.PROBEL2 = "" 
        self.PROBEL3 = "" 
        self.PROBEL3Iso = "" 
        self.NAMEPLOT = "" 

    def __init__(self, TRIGNAME, PROBEL1, PROBEL2, PROBEL3, PROBEL3ISO, NAMEPLOT) :
        self.TRIGNAME = TRIGNAME 
        self.PROBEL1 = PROBEL1 
        self.PROBEL2 = PROBEL2 
        self.PROBEL3 = PROBEL3 
        self.PROBEL3ISO = PROBEL3ISO 
        self.NAMEPLOT = NAMEPLOT
    

# CB extend the path list, may be make paths cff per menu

Mu50 = HltPathConfig("Mu50",
                     "hltL1sL1SingleMu16ORSingleMu25",
                     "hltL2fL1sMu16orMu25L1f0L2Filtered10Q",
                     "hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered50Q",
                     "",
                     "Full")

Mu45 = HltPathConfig("Mu45_eta2p1",
                     "hltL1sL1SingleMu16ORSingleMu25",
                     "hltL2fL1sMu16orMu25L1f0L2Filtered10Q",
                     "hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q",
                     "",
                     "Full")

