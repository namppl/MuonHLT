#ifndef  MuonTree_h
#define  MuonTree_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>



class GenParticleCand {
public:
  Int_t   pdgId; 
  Int_t   status; 
  Float_t energy; 
  Float_t pt; 
  Float_t eta; 
  Float_t phi; 
  std::vector<Int_t>  pdgMother; 
  std::vector<Int_t>  pdgRealMother; 

  GenParticleCand(){};
  virtual ~GenParticleCand(){};
  
  ClassDef(GenParticleCand,1)
};




class MuonCand {
public:

  Float_t pt;  
  Float_t eta; 
  Float_t phi; 
  Int_t   charge;    

  Int_t   isGlobal;
  Int_t   isTracker;

  Int_t   isLoose;
  Int_t   isMedium;
  Int_t   isTight;
  
  Float_t chargedDep_dR03;
  Float_t neutralDep_dR03;
  Float_t photonDep_dR03;
  Float_t puPt_dR03;

  Float_t chargedDep_dR04;
  Float_t neutralDep_dR04;
  Float_t photonDep_dR04;
  Float_t puPt_dR04;

  Float_t ecalPFCluster_dR03;
  Float_t hcalPFCluster_dR03;
  Float_t trkDep_dR03;

  Float_t ecalPFCluster_dR04;
  Float_t hcalPFCluster_dR04;
  Float_t trkDep_dR04;

  MuonCand(){};
  virtual ~MuonCand(){};

  ClassDef(MuonCand,1)
};


class HLTMuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t trkpt;         // pt of the track of the hlt muon [GeV]
  Int_t   charge;         // pt of the track of the hlt muon [GeV]
  
  Float_t ecalDep;
  Float_t hcalDep;
  Float_t trkDep;

  Float_t ecalDep05;
  Float_t hcalDep05;
  Float_t ecalDep1;
  Float_t hcalDep1;

  HLTMuonCand(){};
  virtual ~HLTMuonCand(){};

  ClassDef(HLTMuonCand,1)

};



class HLTObjCand {
public:

  std::string filterTag; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  
  HLTObjCand(){};
  virtual ~HLTObjCand(){};

  ClassDef(HLTObjCand,1)

};





class HLTInfo {
public:
  std::vector<std::string>  triggers;  
  std::vector<HLTObjCand>   objects;   
  double                    rho;

  HLTInfo(){};
  virtual ~HLTInfo(){};
  bool match( const std::string & path ) {
	if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )  return true;
//     if (! iname.compare("HLT_Mu20_v1") == 0) continue;
	return false;
  }

  bool find( const std::string & path ) {
	for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
      if ( it-> compare(path) == 0) return true;
//       if ( it->find ( path ) != std::string::npos ) return true;
	}
	return false;
  }

  ClassDef(HLTInfo,1)

};


class MuonEvent {
public:

  Int_t   runNumber;             
  Int_t   luminosityBlockNumber; 
  Int_t   eventNumber;           

  Int_t   nVtx;                    
  Float_t primaryVertex[3];        
  Float_t cov_primaryVertex[3][3]; 

  Float_t trueNI;   
  Float_t rho; 
  
  Float_t bxId;
  Float_t instLumi; 

  std::vector <GenParticleCand> genParticles; 
  std::vector <MuonCand>        muons;         
  std::vector <HLTMuonCand>     hltmuons;      
  HLTInfo                       hlt;           
  HLTInfo                       hltTag;           

  MuonEvent(){};
  virtual ~MuonEvent(){};

  ClassDef(MuonEvent,1)
};


#endif

