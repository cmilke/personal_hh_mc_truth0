#ifndef mytruth_readtruth_azheavyh_H
#define mytruth_readtruth_azheavyh_H

#include <EventLoop/Algorithm.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <AsgTools/MessageCheck.h>
#include <AthContainers/DataVector.h>
#include <xAODTruth/TruthParticle.h>
#include <xAODTruth/TruthEvent.h>
#include <xAODJet/JetContainer.h>
#include <map>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include "xAODEventInfo/EventInfo.h"

class readtruth_azheavyh : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  bool bbA = false;
  int mH = 300;
  long ctr_posw = 0;
  long ctr_negw = 0;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  const double GeV = 0.001;

  // hist collection
  std::vector<TH1F*> myhist; //!
  TH1F* mkTH1F( const char* n, const char* t, double nb, double x0, double x1);
  const xAOD::TruthParticle* correctedParticle( const xAOD::TruthParticle* part );
  const xAOD::TruthParticle* correctedQuark( const xAOD::TruthParticle* part );
  bool hasParent( const xAOD::TruthParticle* kid, int parentId );
  void printChildren( const xAOD::TruthParticle* part );
  TLorentzVector getJet( const xAOD::TruthParticle* part );

  // my useful functions

  // single particle properties
  
  TH1F* A_pT; //!
  TH1F* A_eta; //!
  TH1F* A_phi; //!
  TH1F* A_E; //!
  TH1F* A_m; //!

  TH1F* H_pT; //!
  TH1F* H_eta; //!
  TH1F* H_phi; //!
  TH1F* H_E; //!
  TH1F* H_m; //!

  TH1F* Z_pT; //!
  TH1F* Z_eta; //!
  TH1F* Z_phi; //!
  TH1F* Z_E; //!
  TH1F* Z_m; //!

  TH1F* ZH_m; //!

  TH1F* l1_pT; //!
  TH1F* l1_eta; //!
  TH1F* l1_phi; //!
  TH1F* l1_E; //!

  TH1F* l2_pT; //!
  TH1F* l2_eta; //!
  TH1F* l2_phi; //!
  TH1F* l2_E; //!

  TH1F* b1_pT; //!
  TH1F* b1_eta; //!
  TH1F* b1_phi; //!
  TH1F* b1_E; //!

  TH1F* b2_pT; //!
  TH1F* b2_eta; //!
  TH1F* b2_phi; //!
  TH1F* b2_E; //!

  // particle system properties
  TH1F* ll_pT; //!
  TH1F* ll_eta; //!
  TH1F* ll_phi; //!
  TH1F* ll_m; //!

  TH1F* bb_pT; //!
  TH1F* bb_eta; //!
  TH1F* bb_phi; //!
  TH1F* bb_m; //!

  TH1F* llbb_pT; //!
  TH1F* llbb_eta; //!
  TH1F* llbb_phi; //!
  TH1F* llbb_m; //!

  // angular distributions
  TH1F* dEta_ll; //!
  TH1F* dPhi_ll; //!
  TH1F* dR_ll; //!
  TH1F* dEta_bb; //!
  TH1F* dPhi_bb; //!
  TH1F* dR_bb; //!

  TH1F* cosTheta_l_Z; //!
  TH1F* cosTheta_b_H; //!


  int lepFlavor; //!
  int ctr_zee; //!
  int ctr_zmm; //!
  int ctr_ztt; //!
  int ctr_tot; //!

  // this is a standard constructor
  readtruth_azheavyh ( bool isbbA = false, int massH = 300 );

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(readtruth_azheavyh, 1);
};

#endif
