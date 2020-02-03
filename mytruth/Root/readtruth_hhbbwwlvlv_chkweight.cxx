#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <mytruth/readtruth_hhbbwwlvlv_chkweight.h>

// this is needed to distribute the algorithm to the workers
ClassImp(readtruth_hhbbwwlvlv_chkweight)

readtruth_hhbbwwlvlv_chkweight :: readtruth_hhbbwwlvlv_chkweight ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}

TH1F* readtruth_hhbbwwlvlv_chkweight :: mkTH1F( const char* n, const char* t, double nb, double x0, double x1){
  TH1F* h = new TH1F( n, t, nb, x0, x1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

void readtruth_hhbbwwlvlv_chkweight :: setOption ( std::string opt ){
  option = opt;
}

TH2F* readtruth_hhbbwwlvlv_chkweight::mkTH2F( const char* n, const char* t, double xnb, double x0, double x1, double ynb, double y0, double y1 ){
  TH2F* h = new TH2F( n, t, xnb, x0, x1, ynb, y0, y1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

// in principal we should look at the particlet in the end of self-decay chains
const xAOD::TruthParticle* readtruth_hhbbwwlvlv_chkweight :: correctedParticle( const xAOD::TruthParticle * part ){
  //Info("correctedParticle()", "Check particle with barcode %d, status %d", part->barcode(), part->status() );
  const xAOD::TruthParticle* corrpart = part;
  for( unsigned ikid=0; ikid<part->nChildren(); ++ikid ){
    const xAOD::TruthParticle* kid = part->child( ikid );
	if( kid->pdgId() == part->pdgId() ){
	  //Info( "correctedParticle()", "[pdgId %d] pT %f, eta %f, barcode %d, parent barcode %d, status %d",
	  //      kid->pdgId(), kid->pt(), kid->eta(), kid->barcode(), kid->parent(0)->barcode(), kid->status() );
	  corrpart = correctedParticle( kid );
	  break;
	}
  }
  return corrpart;
}

/*
// only apply to quarks who has hadronization
// get the quark before they get in hadronization (not look good in bb_m)
// does not work ...
const xAOD::TruthParticle* readtruth_hhbbwwlvlv_chkweight :: correctedQuark( const xAOD::TruthParticle * part ){
  //Info("correctedParticle()", "Check particle with barcode %d, status %d", part->barcode(), part->status() );
  const xAOD::TruthParticle* corrpart = part;
  for( unsigned ikid=0; ikid<part->nChildren(); ++ikid ){
    const xAOD::TruthParticle* kid = part->child( ikid );
	if( kid->pdgId() == part->pdgId() ){
	  //Info( "correctedParticle()", "[pdgId %d] pT %f, eta %f, barcode %d, parent barcode %d, status %d",
	  //      kid->pdgId(), kid->pt(), kid->eta(), kid->barcode(), kid->parent(0)->barcode(), kid->status() );
	  if( kid->status()>=71 ) break;
	  corrpart = correctedParticle( kid );
	  break;
	}
  }
  return corrpart;
}
*/


bool readtruth_hhbbwwlvlv_chkweight :: hasParent( const xAOD::TruthParticle* kid, int parentId ){
  for( unsigned int iparent=0; iparent<kid->nParents(); ++iparent ){
    if( kid->parent(iparent)->pdgId() == parentId ) return true;
  }
  return false;
}

void readtruth_hhbbwwlvlv_chkweight :: printChildren( const xAOD::TruthParticle* part ){
  Info( "printChildren()", "Looking at id %d, pt %f, eta %f, status %d", part->pdgId(), part->pt(), part->eta(), part->status() );
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = part->child(ikid);
    Info( "printChildren()", "id %d, pT %f, eta %f, status %d", kid->pdgId(), kid->pt(), kid->eta(), kid->status() );
  }
}

/*
// collect all hadrons after quark undergoing hadroinization
// sort of making a jet
// not working at all !
TLorentzVector readtruth_hhbbwwlvlv_chkweight :: getJet( const xAOD::TruthParticle* part ){
  TLorentzVector jet;
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid ){
    jet += part->child(ikid)->p4();
  }
  return jet;
}
*/

bool readtruth_hhbbwwlvlv_chkweight :: larger_pT( const xAOD::TruthParticle* p1, const xAOD::TruthParticle* p2 ){
  return p1->pt() > p2->pt();
}

EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  // before openning any input files
  job.useXAOD();

  ANA_CHECK_SET_TYPE( EL::StatusCode );
  ANA_CHECK( xAOD::Init() );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  H1_pT = mkTH1F("H1_pT","p_{T}(H1)", 50, 0, 800);
  H1_eta = mkTH1F("H1_eta","#eta(H1)", 40, -5, 5);
  H1_phi = mkTH1F("H1_phi","#phi(H1)", 40, -3.15, 3.15 );
  H1_E = mkTH1F("H1_E","E(H1)", 50, 0, 2000);
  H1_m = mkTH1F("H1_m","m(H1)", 40, 125-0.04, 125+0.04 );
        
  H2_pT = mkTH1F("H2_pT","p_{T}(H2)", 50, 0, 800);
  H2_eta = mkTH1F("H2_eta","#eta(H2)", 40, -5, 5);
  H2_phi = mkTH1F("H2_phi","#phi(H2)", 40, -3.15, 3.15 );
  H2_E = mkTH1F("H2_E","E(H2)", 50, 0, 2000);
  H2_m = mkTH1F("H2_m","m(H2)", 40, 125-0.04, 125+0.04 );

  HH_pT = mkTH1F("HH_pT","p_{T}(HH)", 1000, 0, 1000); //50 bins
  HH_eta = mkTH1F("HH_eta","#eta(HH)", 40, -10, 10);
  HH_phi = mkTH1F("HH_phi","#phi(HH)", 40, -3.15, 3.15 );
  HH_E = mkTH1F("HH_E","E(HH)", 100, 0, 4000);
  HH_m = mkTH1F("HH_m","m(HH)", 3000, 0, 3000);

  HH_pT_weight = mkTH2F( "HH_pT_weight", "p_{T}(HH) vs weight", 40, 0, 400, 25, -0.5, 0.5 );
  hist_weight = mkTH1F("weight","weight",40,-0.5, 0.5);

  dEta_HH = mkTH1F("dEta_HH","#Delta#eta(H,H)", 40, 0, 8);
  dPhi_HH = mkTH1F("dPhi_HH","#Delta#phi(H,H)", 40, 0, 3.15);
  dR_HH = mkTH1F("dR_HH","#DeltaR(H,H)", 40, 0, 10);

  // b quarks
  b1_pT = mkTH1F("b1_pT","p_{T}(b quark 1)", 50, 0, 800);
  b1_eta = mkTH1F("b1_eta","#eta(b quark 1)", 40, -5, 5);
  b1_phi = mkTH1F("b1_phi","#phi(b quark 1)", 40, -3.15, 3.15 );
  b1_E = mkTH1F("b1_E","E(b quark 1)", 50, 0, 2000);

  b2_pT = mkTH1F("b2_pT","p_{T}(b quark 2)", 50, 0, 800);
  b2_eta = mkTH1F("b2_eta","#eta(b quark 2)", 40, -5, 5);
  b2_phi = mkTH1F("b2_phi","#phi(b quark 2)", 40, -3.15, 3.15 );
  b2_E = mkTH1F("b2_E","E(b quark 2)", 50, 0, 2000);

  // lepton
  el1_pT = mkTH1F("el1_pT","p_{T}(electron 1)", 800, 0, 800);
  el1_eta = mkTH1F("el1_eta","#eta(electron 1)", 100, -5, 5);
  el1_phi = mkTH1F("el1_phi","#phi(electron 1)", 40, -3.15, 3.15 );
  el1_E = mkTH1F("el1_E","E(electron 1)", 50, 0, 2000);

  el2_pT = mkTH1F("el2_pT","p_{T}(electron 2)", 800, 0, 800);
  el2_eta = mkTH1F("el2_eta","#eta(electron 2)", 100, -5, 5);
  el2_phi = mkTH1F("el2_phi","#phi(electron 2)", 40, -3.15, 3.15 );
  el2_E = mkTH1F("el2_E","E(electron 2)", 50, 0, 2000);

  el1el2_pT = mkTH2F("el1el2_pT",";p_{T}(electron 1);p_{T}(electron 2);", 800, 0, 800, 800, 0, 800);

  mu1_pT = mkTH1F("mu1_pT","p_{T}(muon 1)", 800, 0, 800);
  mu1_eta = mkTH1F("mu1_eta","#eta(muon 1)", 100, -5, 5);
  mu1_phi = mkTH1F("mu1_phi","#phi(muon 1)", 40, -3.15, 3.15 );
  mu1_E = mkTH1F("mu1_E","E(muon 1)", 50, 0, 2000);

  mu2_pT = mkTH1F("mu2_pT","p_{T}(muon 2)", 800, 0, 800);
  mu2_eta = mkTH1F("mu2_eta","#eta(muon 2)", 100, -5, 5);
  mu2_phi = mkTH1F("mu2_phi","#phi(muon 2)", 40, -3.15, 3.15 );
  mu2_E = mkTH1F("mu2_E","E(muon 2)", 50, 0, 2000);

  mu1mu2_pT = mkTH2F("mu1mu2_pT",";p_{T}(muon 1);p_{T}(muon 2);", 800, 0, 800, 800, 0, 800);

  tau1_pT = mkTH1F("tau1_pT","p_{T}(tau 1)", 800, 0, 800);
  tau1_eta = mkTH1F("tau1_eta","#eta(tau 1)", 100, -5, 5);
  tau1_phi = mkTH1F("tau1_phi","#phi(tau 1)", 40, -3.15, 3.15 );
  tau1_E = mkTH1F("tau1_E","E(tau 1)", 50, 0, 2000);

  tau2_pT = mkTH1F("tau2_pT","p_{T}(tau 2)", 800, 0, 800);
  tau2_eta = mkTH1F("tau2_eta","#eta(tau 2)", 100, -5, 5);
  tau2_phi = mkTH1F("tau2_phi","#phi(tau 2)", 40, -3.15, 3.15 );
  tau2_E = mkTH1F("tau2_E","E(tau 2)", 50, 0, 2000);

  // lepton
  // e mu leptons from e mu tau any of them, universal
  lep1_pT = mkTH1F("lep1_pT","p_{T}(lep 1)",800, 0, 800);
  lep1_eta = mkTH1F("lep1_eta","#eta(lep 1)",100, -5, 5);
  lep1_phi = mkTH1F("lep1_phi","#phi(lep 1)",40, -3.15, 3.15);
  lep1_E = mkTH1F("lep1_E","E(lep 1)",50, 0, 2000);

  lep2_pT = mkTH1F("lep2_pT","p_{T}(lep 2)",800, 0, 800);
  lep2_eta = mkTH1F("lep2_eta","#eta(lep 2)",100, -5, 5);
  lep2_phi = mkTH1F("lep2_phi","#phi(lep 2)",40, -3.15, 3.15);
  lep2_E = mkTH1F("lep2_E","E(lep 2)",50, 0, 2000);

  // e mu leptons from e mu, exlcuding tau
  lep1_viaemu_pT = mkTH1F("lep1_viaemu_pT","p_{T}(lep 1 from emu)",800, 0, 800);
  lep1_viaemu_eta = mkTH1F("lep1_viaemu_eta","#eta(lep 1 from emu)",100, -5, 5);
  lep1_viaemu_phi = mkTH1F("lep1_viaemu_phi","#phi(lep 1 from emu)",40, -3.15, 3.15);
  lep1_viaemu_E = mkTH1F("lep1_viaemu_E","E(lep 1 from emu)",50, 0, 2000);

  lep2_viaemu_pT = mkTH1F("lep2_viaemu_pT","p_{T}(lep 2 from emu)",800, 0, 800);
  lep2_viaemu_eta = mkTH1F("lep2_viaemu_eta","#eta(lep 2 from emu)",100, -5, 5);
  lep2_viaemu_phi = mkTH1F("lep2_viaemu_phi","#phi(lep 2 from emu)",40, -3.15, 3.15);
  lep2_viaemu_E = mkTH1F("lep2_viaemu_E","E(lep 2 from emu)",50, 0, 2000);

  // e mu leptons from tau exclusively
  lep1_viatau_pT = mkTH1F("lep1_viatau_pT","p_{T}(lep 1 from tau)",800, 0, 800);
  lep1_viatau_eta = mkTH1F("lep1_viatau_eta","#eta(lep 1 from tau)",100, -5, 5);
  lep1_viatau_phi = mkTH1F("lep1_viatau_phi","#phi(lep 1 from tau)",40, -3.15, 3.15);
  lep1_viatau_E = mkTH1F("lep1_viatau_E","E(lep 1 from tau)",50, 0, 2000);

  lep2_viatau_pT = mkTH1F("lep2_viatau_pT","p_{T}(lep 2 from tau)",800, 0, 800);
  lep2_viatau_eta = mkTH1F("lep2_viatau_eta","#eta(lep 2 from tau)",100, -5, 5);
  lep2_viatau_phi = mkTH1F("lep2_viatau_phi","#phi(lep 2 from tau)",40, -3.15, 3.15);
  lep2_viatau_E = mkTH1F("lep2_viatau_E","E(lep 2 from tau)",50, 0, 2000);

  // neutrino
  nu1_pT = mkTH1F("nu1_pT","p_{T}(neutrino 1)", 50, 0, 800);
  nu1_phi = mkTH1F("nu1_phi","#phi(neutrino 1)", 40, -3.15, 3.15 );

  nu2_pT = mkTH1F("nu2_pT","p_{T}(neutrino 2)", 50, 0, 800);
  nu2_phi = mkTH1F("nu2_phi","#phi(neutrino 2)", 40, -3.15, 3.15 );

  // bb for H
  bb_pT = mkTH1F("bb_pT","p_{T}(bb)", 50, 0, 800);
  bb_eta = mkTH1F("bb_eta","#eta(bb)", 40, -5, 5);
  bb_phi = mkTH1F("bb_phi","#phi(bb)", 40, -3.15, 3.15 );
  bb_E = mkTH1F("bb_E","E(bb)", 50, 0, 2000);
  bb_m = mkTH1F("bb_m","m(bb)", 40, 125-50, 125+50 );

  // ee/mumu/tautau

  ll_pT = mkTH1F("ll_pT","p_{T}(ll)", 50, 0, 800);
  ll_eta = mkTH1F("ll_eta","#eta(ll)", 40, -5, 5);
  ll_phi = mkTH1F("ll_phi","#phi(ll)", 40, -3.15, 3.15 );
  ll_E = mkTH1F("ll_E","E(ll)", 50, 0, 2000);
  ll_m = mkTH1F("ll_m","m(ll)", 70, 0, 90+50 );

  // nu nu
  nunu_pT = mkTH1F("nunu_pT","p_{T}(#nu#nu)", 50, 0, 800);
  nunu_phi = mkTH1F("nunu_phi","#phi(#nu#nuu)", 40, -3.15, 3.15 );

  // llvv for H
  llvv_pT = mkTH1F("llvv_pT","p_{T}(ll#nu#nu)", 50, 0, 800);
  llvv_phi = mkTH1F("llvv_phi","#phi(ll#nu#nu)", 40, -3.15, 3.15 );
  llvv_mT = mkTH1F("llvv_mT","m_{T}(ll#nu#nu)", 40, 70, 500 );

  dPhi_ll_vv = mkTH1F("dPhi_ll_vv","#Delta#phi(ll,#nu#nu)", 40, 0, 3.15);

  // bbllvv transver kine for HH
  bbllvv_pT = mkTH1F("bbllvv_pT","p_{T}(bbll#nu#nu)", 1000, 0, 1000); //50 bins
  bbllvv_phi = mkTH1F("bbllvv_phi","#phi(bbll#nu#nu)", 40, -3.15, 3.15 );
  bbllvv_mT = mkTH1F("bbllvv_mT","m_{T}(bbll#nu#nu)", 3000, 0, 3000);

  dPhi_bb_llvv = mkTH1F("dPhi_bb_llvv","#Delta#phi(bb,ll#nu#nu)", 40, 0, 3.15);

  // register
  //for_each( myhist.begin(), myhist.end(), std::bind1st( std::mem_fun( &readtruth_hhbbwwlvlv_chkweight::mbook ), this ) );
  // now done in mkTH1F

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  ANA_CHECK_SET_TYPE( EL::StatusCode );

  xAOD::TEvent* event = wk()->xaodEvent();

  Info("initialize()", "Number of events = %lli", event->getEntries() );

  // my init
  ctr_tot = 0;
  ctr_posw = 0;
  ctr_negw = 0;

  ctr_save = 0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  ANA_CHECK_SET_TYPE( EL::StatusCode );

  xAOD::TEvent* event = wk()->xaodEvent();


  // my init
  ++ctr_tot;

  //Info("execute()", "++++++++++++++++++++++++++++++++++++++++++");

  // truth containers
  const DataVector<xAOD::TruthEvent>* mcevt = NULL;
  ANA_CHECK( event->retrieve( mcevt, "TruthEvents" ) );
  //Info("execute()", "number of TruthEvent = %ld", mcevt->size() ); // these always one TruthEvent :)

  const DataVector<xAOD::TruthParticle>* mcpart = NULL;
  ANA_CHECK( event->retrieve( mcpart, "TruthParticles" ) );
  int npart = mcpart->size();
  //Info( "execute()", "number of TruthParticles = %ld", npart );

  // MC weight
  double weight = mcevt->at(0)->weights()[0];
  //Info("execute()", "weight = %f", weight );
  if( weight>0 ) ++ctr_posw;
  else ++ctr_negw;

  // test weight variations
  for( unsigned int iw=0; iw<46; ++iw ){
    Error("execute()", "Print weight [%d] = %f", iw, mcevt->at(0)->weights()[iw] );
  }
  
  //if( abs(weight) > 5 ){
  //  Error( "execute","MC weight large (abs>0.5): %f", weight );
  //  weight = 0; // ignore these events, taking up <<1 per mille
  //}

//  if( weight<0 ){
//  }else{
//    return EL::StatusCode::SUCCESS;
//  }

  // my collection of truth particles
  std::map<std::string, const xAOD::TruthParticle*> mycollection; // historical for all, but now only Higgs and Z before decay
  std::vector<const xAOD::TruthParticle*> mybquarks; // all 2 b quarks ordere by pT
  std::vector<const xAOD::TruthParticle*> myelectrons; // all 2 electrons ordere by pT
  std::vector<const xAOD::TruthParticle*> mymuons; // all 2 muons ordere by pT
  std::vector<const xAOD::TruthParticle*> mytaus; // all 2 taus ordere by pT
  std::vector<const xAOD::TruthParticle*> myleptons; // all e mu from either emu or tau
  std::vector<const xAOD::TruthParticle*> myneutrinos; // all 2 neutrinos ordere by pT

  std::vector<const xAOD::TruthParticle*> myleptons_viaemu; // only e mu leps, only fro emu
  std::vector<const xAOD::TruthParticle*> myleptons_viatau; // only e mu leps, onlu from tau lep decay

  for( int ipart=0; ipart<npart; ++ipart ){

    const xAOD::TruthParticle* particle = mcpart->at( ipart );

	// debug
	//if( particle->pdgId() == 25 ){
	//  Info( "execute()", "ID %d, barcode %d, parent ID %d, pT %f, px %f, py %f, pz %f, eta %f, phi %f",
	//       particle->pdgId(), particle->barcode(), particle->parent(0)->pdgId(),
	//       particle->pt(), particle->px(), particle->py(), particle->pz(), particle->eta(), particle->phi() );
	//}

	// H
	// after shower dressing
	if( particle->pdgId() == 25 ){ //
	  particle = correctedParticle( particle );

    // before shower dressing
	//if( particle->pdgId() == 25 && !hasParent( particle, 25 ) ){ // first particle, i.e. LHE file particle

      // general
	  if( mycollection.count("H1") == 0 ) mycollection["H1"] = particle;
	  else{
	    if( particle != mycollection["H1"] ) mycollection["H2"] = particle;
	  }
	  // sort by pT later

	  // debug
	  //Info( "execute()", "SAVE ID %d, barcode %d, parent ID %d, pT %f, px %f, py %f, pz %f, eta %f, phi %f",
	  //     particle->pdgId(), particle->barcode(), particle->parent(0)->pdgId(),
	  //     particle->pt(), particle->px(), particle->py(), particle->pz(), particle->eta(), particle->phi() );
	}

	// W
	if( abs(particle->pdgId()) == 24 ){
	  particle = correctedParticle( particle ); // after shower dressing

	  if( mycollection.count("W1") == 0 ) mycollection["W1"] = particle;
	  else{
	    if( particle != mycollection["W1"] ) mycollection["W2"] = particle;
	  }
	}

	// once found H1 H2 quit loop (only used without looking into H decays)
	if( mycollection.count("H1")==1 && mycollection.count("H2")==1
	&&  mycollection.count("W1")==1 && mycollection.count("W2")==1 )
	  break;

  }

  // check
  if( mycollection.count("H1")==1 && mycollection.count("H2")==1
  &&  mycollection.count("W1")==1 && mycollection.count("W2")==1 ){
    //++ctr_save; // moved down
  }else{
    Error("execute()","cannot fine a pair of Higgs, skip event"); return EL::StatusCode::SUCCESS;
  }

  // sort by pT
  double _h1_pT = mycollection["H1"]->pt();
  double _h2_pT = mycollection["H2"]->pt();
  if( _h1_pT < _h2_pT ){
    const xAOD::TruthParticle* _h2 = mycollection["H1"];
    const xAOD::TruthParticle* _h1 = mycollection["H2"];
	mycollection["H1"] = _h1;
	mycollection["H2"] = _h2;
  }
  _h1_pT = mycollection["W1"]->pt();
  _h2_pT = mycollection["W2"]->pt();
  if( _h1_pT < _h2_pT ){
    const xAOD::TruthParticle* _h2 = mycollection["W1"];
    const xAOD::TruthParticle* _h1 = mycollection["W2"];
	mycollection["W1"] = _h1;
	mycollection["W2"] = _h2;
  }

  // find bb from H
  for( unsigned int ikid=0; ikid<mycollection["H1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H1"]->child(ikid);
	if( abs(kid->pdgId()) == 24 ) break; // skip H as H->WW
	if( abs(kid->pdgId()) == 5 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mybquarks.begin(), mybquarks.end(), kid) == mybquarks.end() ){
	    mybquarks.push_back(kid);
	  }
	}
	if( mybquarks.size()==2 ) break;
  }
  for( unsigned int ikid=0; ikid<mycollection["H2"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H2"]->child(ikid);
	if( abs(kid->pdgId()) == 24 ) break; // skip H as H->WW
	if( abs(kid->pdgId()) == 5 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mybquarks.begin(), mybquarks.end(), kid) == mybquarks.end() ){
	    mybquarks.push_back(kid);
	  }
	}
	if( mybquarks.size()==2 ) break;
  }
  // all 2 b quarks ordered by pT
  std::sort( mybquarks.begin(), mybquarks.end(), larger_pT );




  // find lvlv from WW
  for( unsigned int ikid=0; ikid<mycollection["W1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["W1"]->child(ikid);
	if( abs(kid->pdgId()) == 11){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myelectrons.begin(), myelectrons.end(), kid) == myelectrons.end() ){
	    myelectrons.push_back(kid);
		myleptons.push_back(kid);
		myleptons_viaemu.push_back(kid);
	  }
	}else if( abs(kid->pdgId()) ==13){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mymuons.begin(), mymuons.end(), kid) == mymuons.end() ){
	    mymuons.push_back(kid);
		myleptons.push_back(kid);
		myleptons_viaemu.push_back(kid);
	  }
	}else if( abs(kid->pdgId()) == 15){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mytaus.begin(), mytaus.end(), kid) == mytaus.end() ){
	    mytaus.push_back(kid);
		// collect tau -> emu
		for( unsigned int _ikid_tau=0; _ikid_tau<kid->nChildren(); ++_ikid_tau ){
		  const xAOD::TruthParticle* _kid_tau = kid->child(_ikid_tau);
		  if( abs(_kid_tau->pdgId()) == 11 or abs(_kid_tau->pdgId()) == 13 ){
		    _kid_tau = correctedParticle( _kid_tau );
		    myleptons.push_back(_kid_tau);
		    myleptons_viatau.push_back(kid);
		  }
		}
	  }
	}else if( abs(kid->pdgId()) == 12 || abs(kid->pdgId()) == 14 || abs(kid->pdgId()) == 16 ){
	  // bbWW lvlv different than bbZZ, v here has all flavours
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myneutrinos.begin(), myneutrinos.end(), kid) == myneutrinos.end() ){
	    myneutrinos.push_back(kid);
	  }
	}
  }
  for( unsigned int ikid=0; ikid<mycollection["W2"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["W2"]->child(ikid);
	if( abs(kid->pdgId()) == 11){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myelectrons.begin(), myelectrons.end(), kid) == myelectrons.end() ){
	    myelectrons.push_back(kid);
		myleptons.push_back(kid);
		myleptons_viaemu.push_back(kid);
	  }
	}else if( abs(kid->pdgId()) == 13){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mymuons.begin(), mymuons.end(), kid) == mymuons.end() ){
	    mymuons.push_back(kid);
		myleptons.push_back(kid);
		myleptons_viaemu.push_back(kid);
	  }
	}else if( abs(kid->pdgId()) == 15){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mytaus.begin(), mytaus.end(), kid) == mytaus.end() ){
	    mytaus.push_back(kid);
		// collect tau -> emu
		for( unsigned int _ikid_tau=0; _ikid_tau<kid->nChildren(); ++_ikid_tau ){
		  const xAOD::TruthParticle* _kid_tau = kid->child(_ikid_tau);
		  if( abs(_kid_tau->pdgId()) == 11 or abs(_kid_tau->pdgId()) == 13 ){
		    _kid_tau = correctedParticle( _kid_tau );
		    myleptons.push_back(_kid_tau);
		    myleptons_viatau.push_back(kid);
		  }
		}
	  }
	}else if( abs(kid->pdgId()) == 12 || abs(kid->pdgId()) == 14 || abs(kid->pdgId()) == 16 ){
	  // bbWW lvlv different than bbZZ, v here has all flavours
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myneutrinos.begin(), myneutrinos.end(), kid) == myneutrinos.end() ){
	    myneutrinos.push_back(kid);
	  }
	}
  }
  // all leptons ordered by pT
  std::sort( myelectrons.begin(), myelectrons.end(), larger_pT );
  std::sort( mymuons.begin(), mymuons.end(), larger_pT );
  std::sort( mytaus.begin(), mytaus.end(), larger_pT );
  std::sort( myleptons.begin(), myleptons.end(), larger_pT );
  std::sort( myneutrinos.begin(), myneutrinos.end(), larger_pT );

  std::sort( myleptons_viaemu.begin(), myleptons_viaemu.end(), larger_pT );
  std::sort( myleptons_viatau.begin(), myleptons_viatau.end(), larger_pT );

  if( myneutrinos.size() < 2 ){
    //Error("execute()","myneutrinos %d !!!", myneutrinos.size());
    //Error("execute()","myelectrons %d !!!", myelectrons.size());
    //Error("execute()","mymuons %d !!!", mymuons.size());
    //Error("execute()","mytaus %d !!!", mytaus.size());
    //Error("execute()","myleptons %d !!!", myleptons.size());
	Error("execute()","SKIP event: less than 2 neutrino is found!!!");
	return EL::StatusCode::SUCCESS;
  }

  // mimic truth filter in JO
  if( option.find("mimicfilter") != std::string::npos ){
    unsigned int _nlepton = 0;
    for( auto _part : myleptons ){
      if( _part->pt()*GeV > 7 && abs(_part->eta()) < 2.8 )
        ++_nlepton;
    }
    if( _nlepton >= 2 ){
      // pass
    }else{
      return EL::StatusCode::SUCCESS;
    }
  }
  ++ctr_save;

  // test
  //std::cout << "myelectrons: " << myelectrons.size() << std::endl;
  //std::cout << "mymuons: " << mymuons.size() << std::endl;
  //std::cout << "mytaus: " << mytaus.size() << std::endl;
  //std::cout << "myleptons: " << myleptons.size() << std::endl;
  //std::cout << "myneutrinos: " << myneutrinos.size() << std::endl;

  // fill the histograms

  // bosons
  H1_pT->Fill( mycollection["H1"]->pt()*GeV, weight );
  H1_eta->Fill( mycollection["H1"]->eta(), weight );
  H1_phi->Fill( mycollection["H1"]->phi(), weight );
  H1_E->Fill( mycollection["H1"]->e()*GeV, weight );
  H1_m->Fill( mycollection["H1"]->m()*GeV, weight );

  H2_pT->Fill( mycollection["H2"]->pt()*GeV, weight );
  H2_eta->Fill( mycollection["H2"]->eta(), weight );
  H2_phi->Fill( mycollection["H2"]->phi(), weight );
  H2_E->Fill( mycollection["H2"]->e()*GeV, weight );
  H2_m->Fill( mycollection["H2"]->m()*GeV, weight );

  // system var
  TLorentzVector v_H1 = mycollection["H1"]->p4();
  TLorentzVector v_H2 = mycollection["H2"]->p4();
  TLorentzVector v_HH = v_H1 + v_H2;


  HH_pT->Fill( v_HH.Pt()*GeV, weight );
//  std::cout << "DDD pTHH: " << v_HH.Pt()*GeV << " weight: " << weight << std::endl;
//  HH_pT->Fill( sqrt(pow(mycollection["H1"]->px()+mycollection["H2"]->px(), 2)+pow(mycollection["H1"]->py()+mycollection["H2"]->py(), 2))*GeV, weight);
  HH_eta->Fill( v_HH.Eta(), weight );
  HH_phi->Fill( v_HH.Phi(), weight );
  HH_E->Fill( v_HH.E()*GeV, weight );
  HH_m->Fill( v_HH.M()*GeV, weight );

  HH_pT_weight->Fill( v_HH.Pt()*GeV, weight, 1 );
  hist_weight->Fill(weight,1);

  //std::cout << "DDD v_HH.Phi(): " << v_HH.Phi() << std::endl;

  // angular var
  dEta_HH->Fill( std::abs( v_H1.Eta() - v_H2.Eta() ), weight);
  dPhi_HH->Fill( v_H1.DeltaPhi(v_H2), weight);
  dR_HH->Fill( v_H1.DeltaR(v_H2), weight);

  // b quarks
  b1_pT->Fill( mybquarks[0]->pt()*GeV, weight );
  b1_eta->Fill( mybquarks[0]->eta(), weight );
  b1_phi->Fill( mybquarks[0]->phi(), weight );
  b1_E->Fill( mybquarks[0]->e()*GeV, weight );

  b2_pT->Fill( mybquarks[1]->pt()*GeV, weight );
  b2_eta->Fill( mybquarks[1]->eta(), weight );
  b2_phi->Fill( mybquarks[1]->phi(), weight );
  b2_E->Fill( mybquarks[1]->e()*GeV, weight );

  // leptons
  if( myelectrons.size() == 2 ){
  el1_pT  -> Fill( myelectrons[0]->pt()*GeV, weight );
  el1_eta -> Fill( myelectrons[0]->eta(), weight );
  el1_phi -> Fill( myelectrons[0]->phi(), weight );
  el1_E   -> Fill( myelectrons[0]->e()*GeV, weight );

  el2_pT  -> Fill( myelectrons[1]->pt()*GeV, weight );
  el2_eta -> Fill( myelectrons[1]->eta(), weight );
  el2_phi -> Fill( myelectrons[1]->phi(), weight );
  el2_E   -> Fill( myelectrons[1]->e()*GeV, weight );

  el1el2_pT -> Fill( myelectrons[0]->pt()*GeV, myelectrons[1]->pt()*GeV, weight);
  }
  else if( mymuons.size() == 2 ){
  mu1_pT  -> Fill( mymuons[0]->pt()*GeV, weight );
  mu1_eta -> Fill( mymuons[0]->eta(), weight );
  mu1_phi -> Fill( mymuons[0]->phi(), weight );
  mu1_E   -> Fill( mymuons[0]->e()*GeV, weight );

  mu2_pT  -> Fill( mymuons[1]->pt()*GeV, weight );
  mu2_eta -> Fill( mymuons[1]->eta(), weight );
  mu2_phi -> Fill( mymuons[1]->phi(), weight );
  mu2_E   -> Fill( mymuons[1]->e()*GeV, weight );

  mu1mu2_pT -> Fill( mymuons[0]->pt()*GeV, mymuons[1]->pt()*GeV, weight);
  }
  else if( mytaus.size() == 2 ){
  tau1_pT  -> Fill( mytaus[0]->pt()*GeV, weight );
  tau1_eta -> Fill( mytaus[0]->eta(), weight );
  tau1_phi -> Fill( mytaus[0]->phi(), weight );
  tau1_E   -> Fill( mytaus[0]->e()*GeV, weight );

  tau2_pT  -> Fill( mytaus[1]->pt()*GeV, weight );
  tau2_eta -> Fill( mytaus[1]->eta(), weight );
  tau2_phi -> Fill( mytaus[1]->phi(), weight );
  tau2_E   -> Fill( mytaus[1]->e()*GeV, weight );
  }

  // unversial leptons e mu that come from e mu and tau !
  if( myleptons.size() >=2 ){
  lep1_pT  -> Fill( myleptons[0]->pt()*GeV, weight );
  lep1_eta -> Fill( myleptons[0]->eta(), weight );
  lep1_phi -> Fill( myleptons[0]->phi(), weight );
  lep1_E   -> Fill( myleptons[0]->e()*GeV, weight );

  lep2_pT  -> Fill( myleptons[1]->pt()*GeV, weight );
  lep2_eta -> Fill( myleptons[1]->eta(), weight );
  lep2_phi -> Fill( myleptons[1]->phi(), weight );
  lep2_E   -> Fill( myleptons[1]->e()*GeV, weight );
  }
  if( myleptons_viaemu.size() >=2 ){
  lep1_viaemu_pT  -> Fill( myleptons_viaemu[0]->pt()*GeV, weight );
  lep1_viaemu_eta -> Fill( myleptons_viaemu[0]->eta(), weight );
  lep1_viaemu_phi -> Fill( myleptons_viaemu[0]->phi(), weight );
  lep1_viaemu_E   -> Fill( myleptons_viaemu[0]->e()*GeV, weight );

  lep2_viaemu_pT  -> Fill( myleptons_viaemu[1]->pt()*GeV, weight );
  lep2_viaemu_eta -> Fill( myleptons_viaemu[1]->eta(), weight );
  lep2_viaemu_phi -> Fill( myleptons_viaemu[1]->phi(), weight );
  lep2_viaemu_E   -> Fill( myleptons_viaemu[1]->e()*GeV, weight );
  }
  if( myleptons_viatau.size() >=2 ){
  lep1_viatau_pT  -> Fill( myleptons_viatau[0]->pt()*GeV, weight );
  lep1_viatau_eta -> Fill( myleptons_viatau[0]->eta(), weight );
  lep1_viatau_phi -> Fill( myleptons_viatau[0]->phi(), weight );
  lep1_viatau_E   -> Fill( myleptons_viatau[0]->e()*GeV, weight );

  lep2_viatau_pT  -> Fill( myleptons_viatau[1]->pt()*GeV, weight );
  lep2_viatau_eta -> Fill( myleptons_viatau[1]->eta(), weight );
  lep2_viatau_phi -> Fill( myleptons_viatau[1]->phi(), weight );
  lep2_viatau_E   -> Fill( myleptons_viatau[1]->e()*GeV, weight );
  }

  // neutrino
  nu1_pT   -> Fill( myneutrinos[0]->pt()*GeV, weight );
  nu1_phi  -> Fill( myneutrinos[0]->phi(), weight );

  nu2_pT   -> Fill( myneutrinos[1]->pt()*GeV, weight );
  nu2_phi  -> Fill( myneutrinos[1]->phi(), weight );

  // bb for H
  TLorentzVector v_bb = mybquarks[0]->p4() + mybquarks[1]->p4();
  bb_pT->Fill( v_bb.Pt()*GeV, weight );
  bb_eta->Fill( v_bb.Eta(), weight );
  bb_phi->Fill( v_bb.Phi(), weight );
  bb_E->Fill( v_bb.E()*GeV, weight );
  bb_m->Fill( v_bb.M()*GeV, weight );

  // below call lvlv as llvv, as the class is derived from ZZ llvv
  // but physics is valid

  // vv for partial WW
  TLorentzVector v_vv = myneutrinos[0]->p4() + myneutrinos[1]->p4();
  nunu_pT  -> Fill( v_vv.Pt()*GeV, weight );
  nunu_phi -> Fill( v_vv.Phi(), weight );

  // ll for partial WW
  if( myleptons.size() >=2 ){
  TLorentzVector v_ll = myleptons[0]->p4() + myleptons[1]->p4();
  ll_pT  -> Fill( v_ll.Pt()*GeV, weight );
  ll_eta -> Fill( v_ll.Eta(), weight );
  ll_phi -> Fill( v_ll.Phi(), weight );
  ll_E   -> Fill( v_ll.E()*GeV, weight );
  ll_m   -> Fill( v_ll.M()*GeV, weight );

  // llvv
  TLorentzVector v_llvv = v_ll + v_vv;
  // mT^2 = (ET_vis+ET_miss)^2 - (vec{pT_vis}+vec{ET_miss})^2
  // ET_vis = sqrt( vec{pT_vis}^2 + m_vis^2 )
  /*
  TLorentzVector v_vis = myleptons[0]->p4()+myleptons[1]->p4()+mybquarks[0]->p4()+mybquarks[1]->p4(); // vec_llbb
  double ET_vis = TMath::Sqrt( v_vis.Pt()*v_vis.Pt() // vec{pT_vis}^2 pT_llbb
                            + v_vis.M()*v_vis.M() ); // m_vis llbb
  double mT2 = TMath::Power( v_vis.Et() + v_vv.Et(), 2 ) // tot ET miss
             - TMath::Power( (v_vis+v_vv).Pt(), 2 );
  */
  llvv_pT  -> Fill( v_llvv.Pt()*GeV, weight );
  llvv_phi -> Fill( v_llvv.Phi(), weight );
  //llvv_mT  -> Fill( TMath::Sqrt(mT2)*GeV, weight );
  llvv_mT  -> Fill( v_llvv.Mt()*GeV, weight );

  dPhi_ll_vv -> Fill( v_ll.DeltaPhi(v_vv), weight );

  // bbllvv for HH
  TLorentzVector v_bbllvv = v_bb + v_llvv;
  bbllvv_pT  -> Fill( v_bbllvv.Pt()*GeV, weight );
  bbllvv_phi -> Fill( v_bbllvv.Phi(), weight );
  bbllvv_mT  -> Fill( v_bbllvv.Mt()*GeV, weight );

  // angular
  dPhi_bb_llvv->Fill( v_bb.DeltaPhi(v_llvv), weight);
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  ANA_CHECK_SET_TYPE( EL::StatusCode );

  xAOD::TEvent* event = wk()->xaodEvent();

  Error("finalize()", "Total events: %d, truth identified (HH+WW) events: %d", ctr_tot, int(ctr_save) );

  Error("finalize()", "Wenu events: %f", el1_pT->GetEntries() );
  Error("finalize()", "Wmunu events: %f", mu1_pT->GetEntries() );
  Error("finalize()", "Wtaunu events: %f", tau1_pT->GetEntries() );
  Error("finalize()", "Neutrino events: %f", nu1_pT->GetEntries() );

  // weight
  Error("finalize()", "# evt with positive weight: %ld", ctr_posw);
  Error("finalize()", "# evt with negative weight: %ld", ctr_negw);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbwwlvlv_chkweight :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
