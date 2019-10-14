#ifndef mytruth_readtruth_hhbbzzllvv_H
#define mytruth_readtruth_hhbbzzllvv_H

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
#include <TH2F.h>
#include "xAODEventInfo/EventInfo.h"

//enum HHMODEL { non_resonant, scalar, graviton};
//enum HHCHANNEL { bbbb, bbyy, bbtt_lh, bbtt_hh, wwyy, bbww, wwww, bbtt_ll, tttt, bbZZ };

class readtruth_hhbbzzllvv : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  // default values all set for non-resonant
  int mX = -1;
//  HHMODEL model = non_resonant;
//  HHCHANNEL channel = bbZZ;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  const double GeV = 0.001;

  // hist collection
  std::vector<TH1*> myhist; //!
  TH1F* mkTH1F( const char* n, const char* t, double nb, double x0, double x1);
  TH2F* mkTH2F( const char* n, const char* t, double xnb, double x0, double x1, double ynb, double y0, double y1 );
  const xAOD::TruthParticle* correctedParticle( const xAOD::TruthParticle* part );
  const xAOD::TruthParticle* correctedQuark( const xAOD::TruthParticle* part );
  bool hasParent( const xAOD::TruthParticle* kid, int parentId );
  void printChildren( const xAOD::TruthParticle* part );
  TLorentzVector getJet( const xAOD::TruthParticle* part );

  static bool larger_pT( const xAOD::TruthParticle* p1, const xAOD::TruthParticle* p2 );

  // single particle properties
  
  TH1F* H1_pT; //!
  TH1F* H1_eta; //!
  TH1F* H1_phi; //!
  TH1F* H1_E; //!
  TH1F* H1_m; //!

  TH1F* H2_pT; //!
  TH1F* H2_eta; //!
  TH1F* H2_phi; //!
  TH1F* H2_E; //!
  TH1F* H2_m; //!

  // system properties
  TH1F* HH_pT; //!
  TH1F* HH_eta; //!
  TH1F* HH_phi; //!
  TH1F* HH_E; //!
  TH1F* HH_m; //!

  TH2F* HH_pT_weight; //!
  TH1F* hist_weight; //!

  // final state particles in parton level

  //TH1F* b1_pT; //!

  // angular distributions
  TH1F* dEta_HH; //!
  TH1F* dPhi_HH; //!
  TH1F* dR_HH; //!

  //TH1F* cosTheta_l_Z; //!
  //TH1F* cosTheta_b_H; //!

  // b quarks
  TH1F* b1_pT; //!
  TH1F* b1_eta; //!
  TH1F* b1_phi; //!
  TH1F* b1_E; //!

  TH1F* b2_pT; //!
  TH1F* b2_eta; //!
  TH1F* b2_phi; //!
  TH1F* b2_E; //!

  // lepton
  TH1F*  el1_pT ; //!
  TH1F*  el1_eta ; //!
  TH1F*  el1_phi ; //!
  TH1F*  el1_E ; //!

  TH1F*  el2_pT ; //!
  TH1F*  el2_eta ; //!
  TH1F*  el2_phi ; //!
  TH1F*  el2_E ; //!

  TH1F*  mu1_pT ; //!
  TH1F*  mu1_eta ; //!
  TH1F*  mu1_phi ; //!
  TH1F*  mu1_E ; //!

  TH1F*  mu2_pT ; //!
  TH1F*  mu2_eta ; //!
  TH1F*  mu2_phi ; //!
  TH1F*  mu2_E ; //!

  TH1F*  tau1_pT ; //!
  TH1F*  tau1_eta ; //!
  TH1F*  tau1_phi ; //!
  TH1F*  tau1_E ; //!

  TH1F*  tau2_pT ; //!
  TH1F*  tau2_eta ; //!
  TH1F*  tau2_phi ; //!
  TH1F*  tau2_E ; //!

  // neutrino
  TH1F*  nu1_pT ; //!
  TH1F*  nu1_phi ; //!

  TH1F*  nu2_pT ; //!
  TH1F*  nu2_phi ; //!

  // bb for H
  TH1F*  bb_pT ; //!
  TH1F*  bb_eta ; //!
  TH1F*  bb_phi ; //!
  TH1F*  bb_E ; //!
  TH1F*  bb_m ; //!

  // ee/mumu/tautau

  TH1F*  ll_pT ; //!
  TH1F*  ll_eta ; //!
  TH1F*  ll_phi ; //!
  TH1F*  ll_E ; //!
  TH1F*  ll_m ; //!

  // nu nu
  TH1F*  nunu_pT ; //!
  TH1F*  nunu_phi ; //!

  // llvv for H
  TH1F* llvv_pT; //!
  TH1F* llvv_phi; //!
  TH1F* llvv_mT; //!
  TH1F* dPhi_ll_vv; //!

  // bbllvv transver kine for HH
  TH1F*  bbllvv_pT ; //!
  TH1F*  bbllvv_phi ; //!
  TH1F*  bbllvv_mT ; //!

  TH1F*  dPhi_bb_llvv ; //!

  // counters
  int ctr_tot; //!
  int ctr_posw; //!
  int ctr_negw; //!
  int ctr_save; //!

  // this is a standard constructor
  readtruth_hhbbzzllvv ();

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
  ClassDef(readtruth_hhbbzzllvv, 1);
};

#endif
