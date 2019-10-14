#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <mytruth/readtruth_hhbbzz4l.h>

// this is needed to distribute the algorithm to the workers
ClassImp(readtruth_hhbbzz4l)

readtruth_hhbbzz4l :: readtruth_hhbbzz4l ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}

TH1F* readtruth_hhbbzz4l :: mkTH1F( const char* n, const char* t, double nb, double x0, double x1){
  TH1F* h = new TH1F( n, t, nb, x0, x1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

TH2F* readtruth_hhbbzz4l::mkTH2F( const char* n, const char* t, double xnb, double x0, double x1, double ynb, double y0, double y1 ){
  TH2F* h = new TH2F( n, t, xnb, x0, x1, ynb, y0, y1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

// in principal we should look at the particlet in the end of self-decay chains
const xAOD::TruthParticle* readtruth_hhbbzz4l :: correctedParticle( const xAOD::TruthParticle * part ){
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
const xAOD::TruthParticle* readtruth_hhbbzz4l :: correctedQuark( const xAOD::TruthParticle * part ){
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


bool readtruth_hhbbzz4l :: hasParent( const xAOD::TruthParticle* kid, int parentId ){
  for( unsigned int iparent=0; iparent<kid->nParents(); ++iparent ){
    if( kid->parent(iparent)->pdgId() == parentId ) return true;
  }
  return false;
}

void readtruth_hhbbzz4l :: printChildren( const xAOD::TruthParticle* part ){
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
TLorentzVector readtruth_hhbbzz4l :: getJet( const xAOD::TruthParticle* part ){
  TLorentzVector jet;
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid ){
    jet += part->child(ikid)->p4();
  }
  return jet;
}
*/

bool readtruth_hhbbzz4l :: larger_pT( const xAOD::TruthParticle* p1, const xAOD::TruthParticle* p2 ){
  return p1->pt() > p2->pt();
}

EL::StatusCode readtruth_hhbbzz4l :: setupJob (EL::Job& job)
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



EL::StatusCode readtruth_hhbbzz4l :: histInitialize ()
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
  l1_pT = mkTH1F("l1_pT","p_{T}(lepton 1)", 50, 0, 500);
  l1_eta = mkTH1F("l1_eta","#eta(lepton 1)", 40, -5, 5);
  l1_phi = mkTH1F("l1_phi","#phi(lepton 1)", 40, -3.15, 3.15 );
  l1_E = mkTH1F("l1_E","E(lepton 1)", 50, 0, 2000);

  l2_pT = mkTH1F("l2_pT","p_{T}(lepton 2)", 50, 0, 250);
  l2_eta = mkTH1F("l2_eta","#eta(lepton 2)", 40, -5, 5);
  l2_phi = mkTH1F("l2_phi","#phi(lepton 2)", 40, -3.15, 3.15 );
  l2_E = mkTH1F("l2_E","E(lepton 2)", 50, 0, 2000);

  l3_pT = mkTH1F("l3_pT","p_{T}(lepton 3)", 50, 0, 150);
  l3_eta = mkTH1F("l3_eta","#eta(lepton 3)", 40, -5, 5);
  l3_phi = mkTH1F("l3_phi","#phi(lepton 3)", 40, -3.15, 3.15 );
  l3_E = mkTH1F("l3_E","E(lepton 3)", 50, 0, 2000);

  l4_pT = mkTH1F("l4_pT","p_{T}(lepton 4)", 50, 0, 100);
  l4_eta = mkTH1F("l4_eta","#eta(lepton 4)", 40, -5, 5);
  l4_phi = mkTH1F("l4_phi","#phi(lepton 4)", 40, -3.15, 3.15 );
  l4_E = mkTH1F("l4_E","E(lepton 4)", 50, 0, 2000);

  // bb for H
  bb_pT = mkTH1F("bb_pT","p_{T}(bb)", 50, 0, 800);
  bb_eta = mkTH1F("bb_eta","#eta(bb)", 40, -5, 5);
  bb_phi = mkTH1F("bb_phi","#phi(bb)", 40, -3.15, 3.15 );
  bb_E = mkTH1F("bb_E","E(bb)", 50, 0, 2000);
  bb_m = mkTH1F("bb_m","m(bb)", 40, 125-50, 125+50 );

  // ee/mumu/tautau
  ll1_pT = mkTH1F("ll1_pT","p_{T}(ll 1)", 50, 0, 800);
  ll1_eta = mkTH1F("ll1_eta","#eta(ll 1)", 40, -5, 5);
  ll1_phi = mkTH1F("ll1_phi","#phi(ll 1)", 40, -3.15, 3.15 );
  ll1_E = mkTH1F("ll1_E","E(ll 1)", 50, 0, 2000);
  ll1_m = mkTH1F("ll1_m","m(ll 1)", 70, 0, 90+50 );

  ll2_pT = mkTH1F("ll2_pT","p_{T}(ll 2)", 50, 0, 800);
  ll2_eta = mkTH1F("ll2_eta","#eta(ll 2)", 40, -5, 5);
  ll2_phi = mkTH1F("ll2_phi","#phi(ll 2)", 40, -3.15, 3.15 );
  ll2_E = mkTH1F("ll2_E","E(ll 2)", 50, 0, 2000);
  ll2_m = mkTH1F("ll2_m","m(ll 2)", 70, 0, 90+50 );

  // llll for H
  llll_pT = mkTH1F("llll_pT","p_{T}(llll)", 50, 0, 800);
  llll_eta = mkTH1F("llll_eta","#eta(llll)", 40, -5, 5);
  llll_phi = mkTH1F("llll_phi","#phi(llll)", 40, -3.15, 3.15 );
  llll_E = mkTH1F("llll_E","E(llll)", 50, 0, 2000);
  llll_m = mkTH1F("llll_m","m(llll)", 70, 125-50, 125+20 );

  dEta_ll_ll = mkTH1F("dEta_ll_ll","#Delta#eta(ll,ll)", 40, 0, 8);
  dPhi_ll_ll = mkTH1F("dPhi_ll_ll","#Delta#phi(ll,ll)", 40, 0, 3.15);
  dR_ll_ll = mkTH1F("dR_ll_ll","#DeltaR(ll,ll)", 40, 0, 10);

  // bbllll for HH
  bbllll_pT = mkTH1F("bbllll_pT","p_{T}(bbllll)", 1000, 0, 1000);
  bbllll_eta = mkTH1F("bbllll_eta","#eta(bbllll)", 40, -5, 5);
  bbllll_phi = mkTH1F("bbllll_phi","#phi(bbllll)", 40, -3.15, 3.15 );
  bbllll_E = mkTH1F("bbllll_E","E(bbllll)", 50, 0, 2000);
  bbllll_m = mkTH1F("bbllll_m","m(bbllll)", 3000, 0, 3000 );

  dEta_bb_llll = mkTH1F("dEta_bb_llll","#Delta#eta(bb,llll)", 40, 0, 8);
  dPhi_bb_llll = mkTH1F("dPhi_bb_llll","#Delta#phi(bb,llll)", 40, 0, 3.15);
  dR_bb_llll = mkTH1F("dR_bb_llll","#DeltaR(bb,llll)", 40, 0, 10);

  // register
  //for_each( myhist.begin(), myhist.end(), std::bind1st( std::mem_fun( &readtruth_hhbbzz4l::mbook ), this ) );
  // now done in mkTH1F

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: initialize ()
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



EL::StatusCode readtruth_hhbbzz4l :: execute ()
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
  
  if( abs(weight) > 0.5 ){
    Error( "execute","MC weight large (abs>0.5): %f", weight );
	weight = 0; // ignore these events, taking up <<1 per mille
  }

//  if( weight<0 ){
//  }else{
//    return EL::StatusCode::SUCCESS;
//  }

  // my collection of truth particles
  std::map<std::string, const xAOD::TruthParticle*> mycollection; // historical for all, but now only Higgs and Z before decay
  std::vector<const xAOD::TruthParticle*> mybquarks; // all 2 b quarks ordere by pT
  std::vector<const xAOD::TruthParticle*> myleptons_Z1; // all 2 leptons ordere by pT from Z1
  std::vector<const xAOD::TruthParticle*> myleptons_Z2; // all 2 leptons ordere by pT from Z2
  std::vector<const xAOD::TruthParticle*> myleptons; // all 2 leptons ordere by pT

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

	// Z
	if( particle->pdgId() == 23 ){
	  particle = correctedParticle( particle ); // after shower dressing

	  if( mycollection.count("Z1") == 0 ) mycollection["Z1"] = particle;
	  else{
	    if( particle != mycollection["Z1"] ) mycollection["Z2"] = particle;
	  }
	}

	// once found H1 H2 quit loop (only used without looking into H decays)
	if( mycollection.count("H1")==1 && mycollection.count("H2")==1
	&&  mycollection.count("Z1")==1 && mycollection.count("Z2")==1 )
	  break;

  }

  // check
  if( mycollection.count("H1")==1 && mycollection.count("H2")==1
  &&  mycollection.count("Z1")==1 && mycollection.count("Z2")==1 ){
    ++ctr_save;
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
  _h1_pT = mycollection["Z1"]->pt();
  _h2_pT = mycollection["Z2"]->pt();
  if( _h1_pT < _h2_pT ){
    const xAOD::TruthParticle* _h2 = mycollection["Z1"];
    const xAOD::TruthParticle* _h1 = mycollection["Z2"];
	mycollection["Z1"] = _h1;
	mycollection["Z2"] = _h2;
  }

  // find bb from H
  for( unsigned int ikid=0; ikid<mycollection["H1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H1"]->child(ikid);
	if( abs(kid->pdgId() == 23 ) ) break; // skip H as H->ZZ
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
	if( abs(kid->pdgId() == 23 ) ) break; // skip H as H->ZZ
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




  // find llll from ZZ
  for( unsigned int ikid=0; ikid<mycollection["Z1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["Z1"]->child(ikid);
	if( abs(kid->pdgId()) == 11 || abs(kid->pdgId()) == 13 || abs(kid->pdgId()) == 15 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myleptons_Z1.begin(), myleptons_Z1.end(), kid) == myleptons_Z1.end() ){
	    myleptons_Z1.push_back(kid);
		myleptons.push_back(kid);
	  }
	}
  }
  for( unsigned int ikid=0; ikid<mycollection["Z2"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["Z2"]->child(ikid);
	if( abs(kid->pdgId()) == 11 || abs(kid->pdgId()) == 13 || abs(kid->pdgId()) == 15 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(myleptons_Z2.begin(), myleptons_Z2.end(), kid) == myleptons_Z2.end() ){
	    myleptons_Z2.push_back(kid);
		myleptons.push_back(kid);
	  }
	}
  }
  // all leptons ordered by pT
  std::sort( myleptons_Z1.begin(), myleptons_Z1.end(), larger_pT );
  std::sort( myleptons_Z2.begin(), myleptons_Z2.end(), larger_pT );
  std::sort( myleptons.begin(), myleptons.end(), larger_pT );

  // test
  //std::cout << "myleptons_Z1: " << myleptons_Z1.size() << std::endl;
  //std::cout << "myleptons_Z2: " << myleptons_Z2.size() << std::endl;
  //std::cout << "myleptons: " << myleptons.size() << std::endl;

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
  l1_pT  -> Fill( myleptons[0]->pt()*GeV, weight );
  l1_eta -> Fill( myleptons[0]->eta(), weight );
  l1_phi -> Fill( myleptons[0]->phi(), weight );
  l1_E   -> Fill( myleptons[0]->e()*GeV, weight );

  l2_pT  -> Fill( myleptons[1]->pt()*GeV, weight );
  l2_eta -> Fill( myleptons[1]->eta(), weight );
  l2_phi -> Fill( myleptons[1]->phi(), weight );
  l2_E   -> Fill( myleptons[1]->e()*GeV, weight );

  l3_pT  -> Fill( myleptons[2]->pt()*GeV, weight );
  l3_eta -> Fill( myleptons[2]->eta(), weight );
  l3_phi -> Fill( myleptons[2]->phi(), weight );
  l3_E   -> Fill( myleptons[2]->e()*GeV, weight );

  l4_pT  -> Fill( myleptons[3]->pt()*GeV, weight );
  l4_eta -> Fill( myleptons[3]->eta(), weight );
  l4_phi -> Fill( myleptons[3]->phi(), weight );
  l4_E   -> Fill( myleptons[3]->e()*GeV, weight );

  // bb for H
  TLorentzVector v_bb = mybquarks[0]->p4() + mybquarks[1]->p4();
  bb_pT->Fill( v_bb.Pt()*GeV, weight );
  bb_eta->Fill( v_bb.Eta(), weight );
  bb_phi->Fill( v_bb.Phi(), weight );
  bb_E->Fill( v_bb.E()*GeV, weight );
  bb_m->Fill( v_bb.M()*GeV, weight );

  // ll for Z
  TLorentzVector v_ll1 = myleptons_Z1[0]->p4() + myleptons_Z1[1]->p4();
  ll1_pT  -> Fill( v_ll1.Pt()*GeV, weight );
  ll1_eta -> Fill( v_ll1.Eta(), weight );
  ll1_phi -> Fill( v_ll1.Phi(), weight );
  ll1_E   -> Fill( v_ll1.E()*GeV, weight );
  ll1_m   -> Fill( v_ll1.M()*GeV, weight );

  TLorentzVector v_ll2 = myleptons_Z2[0]->p4() + myleptons_Z2[1]->p4();
  ll2_pT  -> Fill( v_ll2.Pt()*GeV, weight );
  ll2_eta -> Fill( v_ll2.Eta(), weight );
  ll2_phi -> Fill( v_ll2.Phi(), weight );
  ll2_E   -> Fill( v_ll2.E()*GeV, weight );
  ll2_m   -> Fill( v_ll2.M()*GeV, weight );

  // llll
  TLorentzVector v_llll = v_ll1 + v_ll2;
  llll_pT  -> Fill( v_llll.Pt()*GeV, weight );
  llll_eta -> Fill( v_llll.Eta(), weight );
  llll_phi -> Fill( v_llll.Phi(), weight );
  llll_E   -> Fill( v_llll.E()*GeV, weight );
  llll_m   -> Fill( v_llll.M()*GeV, weight );

  dEta_ll_ll->Fill( std::abs( v_ll1.Eta() - v_ll2.Eta() ), weight);
  dPhi_ll_ll->Fill( v_ll1.DeltaPhi(v_ll2), weight);
  dR_ll_ll->Fill( v_ll1.DeltaR(v_ll2), weight);

  // bbllll for HH
  TLorentzVector v_bbllll = v_bb + v_llll;
  bbllll_pT  -> Fill( v_bbllll.Pt()*GeV, weight );
  bbllll_eta -> Fill( v_bbllll.Eta(), weight );
  bbllll_phi -> Fill( v_bbllll.Phi(), weight );
  bbllll_E   -> Fill( v_bbllll.E()*GeV, weight );
  bbllll_m   -> Fill( v_bbllll.M()*GeV, weight );

  // angular
  dEta_bb_llll->Fill( std::abs( v_bb.Eta() - v_llll.Eta() ), weight);
  dPhi_bb_llll->Fill( v_bb.DeltaPhi(v_llll), weight);
  dR_bb_llll->Fill( v_bb.DeltaR(v_llll), weight);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: finalize ()
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

  Error("finalize()", "Total events: %d, truth identified (HH+ZZ) events: %d", ctr_tot, int(ctr_save) );

  // weight
  Error("finalize()", "# evt with positive weight: %ld", ctr_posw);
  Error("finalize()", "# evt with negative weight: %ld", ctr_negw);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hhbbzz4l :: histFinalize ()
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
