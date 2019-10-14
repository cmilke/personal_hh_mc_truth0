#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <mytruth/readtruth_azheavyh.h>

// this is needed to distribute the algorithm to the workers
ClassImp(readtruth_azheavyh)

readtruth_azheavyh :: readtruth_azheavyh ( bool isbbA, int massH )
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  bbA = isbbA;
  mH = massH;
}

TH1F* readtruth_azheavyh :: mkTH1F( const char* n, const char* t, double nb, double x0, double x1){
  TH1F* h = new TH1F( n, t, nb, x0, x1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

//
// in principal we should look at the particlet in the end of self-decay chains
// but bb_m distribution looks bad (totally random, no peak ...)
// probably the last particle has pulled anti-particle from vacuum already
//
// then uniformlly only look at the head of chains for all particles
// the following function is not in use
const xAOD::TruthParticle* readtruth_azheavyh :: correctedParticle( const xAOD::TruthParticle * part ){
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

// only apply to quarks who has hadronization
// get the quark before they get in hadronization (not look good in bb_m)
// does not work ...
const xAOD::TruthParticle* readtruth_azheavyh :: correctedQuark( const xAOD::TruthParticle * part ){
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


bool readtruth_azheavyh :: hasParent( const xAOD::TruthParticle* kid, int parentId ){
  for( unsigned int iparent=0; iparent<kid->nParents(); ++iparent ){
    if( kid->parent(iparent)->pdgId() == parentId ) return true;
  }
  return false;
}

void readtruth_azheavyh :: printChildren( const xAOD::TruthParticle* part ){
  Info( "printChildren()", "Looking at id %d, pt %f, eta %f, status %d", part->pdgId(), part->pt(), part->eta(), part->status() );
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = part->child(ikid);
    Info( "printChildren()", "id %d, pT %f, eta %f, status %d", kid->pdgId(), kid->pt(), kid->eta(), kid->status() );
  }
}

// collect all hadrons after quark undergoing hadroinization
// sort of making a jet
// not working at all !
TLorentzVector readtruth_azheavyh :: getJet( const xAOD::TruthParticle* part ){
  TLorentzVector jet;
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid ){
    jet += part->child(ikid)->p4();
  }
  return jet;
}

EL::StatusCode readtruth_azheavyh :: setupJob (EL::Job& job)
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



EL::StatusCode readtruth_azheavyh :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  A_pT = mkTH1F("A_pT","p_{T}(A)", 30, 0, 500);
  A_eta = mkTH1F("A_eta","#eta(A)", 30, -5, 5);
  A_phi = mkTH1F("A_phi","#phi(A)", 30, -3.15, 3.15 );
  A_E = mkTH1F("A_E","E(A)", 50, 0, 2500);
  A_m = mkTH1F("A_m","m(A)", 800, 200, 1000);
        
  H_pT = mkTH1F("H_pT","p_{T}(H)", 30, 0, 500); 
  H_eta = mkTH1F("H_eta","#eta(H)", 30, -5, 5);
  H_phi = mkTH1F("H_phi","#phi(H)", 30, -3.15, 3.15);
  H_E = mkTH1F("H_E","E(H)", 50, 0, 2500);
  H_m = mkTH1F("H_m","m(H)", 900, 100, 1000);
        
  Z_pT = mkTH1F("Z_pT","p_{T}(Z)", 30, 0, 500); 
  Z_eta = mkTH1F("Z_eta","#eta(Z)", 30, -5, 5);
  Z_phi = mkTH1F("Z_phi","#phi(Z)", 30, -3.15, 3.15);
  Z_E = mkTH1F("Z_E","E(Z)", 30, 0, 1500);
  Z_m = mkTH1F("Z_m","m(Z)", 60, 60, 120);

  ZH_m = mkTH1F("ZH_m","m(ZH)", 800, 200, 1000 );
        
  l1_pT = mkTH1F("l1_pT","p_{T}(lepton1)", 30, 0, 800);
  l1_eta = mkTH1F("l1_eta","#eta(lepton1)", 30, -5, 5);
  l1_phi = mkTH1F("l1_phi","#phi(lepton1)", 30, -3.15, 3.15);
  l1_E = mkTH1F("l1_E","E(lepton1)", 30, 0, 800); 
        
  l2_pT = mkTH1F("l2_pT","p_{T}(lepton2)", 30, 0, 800);
  l2_eta = mkTH1F("l2_eta","#eta(lepton2)", 30, -5, 5);
  l2_phi = mkTH1F("l2_phi","#phi(lepton2)", 30, -3.15, 3.15);
  l2_E = mkTH1F("l2_E","E(lepton2)", 30, 0, 800); 

  b1_pT = mkTH1F("b1_pT","p_{T}(bjet1)", 30, 0, 1200);
  b1_eta = mkTH1F("b1_eta","#eta(bjet1)", 30, -5, 5);
  b1_phi = mkTH1F("b1_phi","#phi(bjet1)", 30, -3.15, 3.15);
  b1_E = mkTH1F("b1_E","E(bjet1)", 30, 0, 800); 
        
  b2_pT = mkTH1F("b2_pT","p_{T}(bjet2)", 30, 0, 1200);
  b2_eta = mkTH1F("b2_eta","#eta(bjet2)", 30, -5, 5);
  b2_phi = mkTH1F("b2_phi","#phi(bjet2)", 30, -3.15, 3.15);
  b2_E = mkTH1F("b2_E","E(bjet2)", 30, 0, 800); 

  ll_pT = mkTH1F("ll_pT","p_{T}(ll)", 30, 0, 500);
  ll_eta = mkTH1F("ll_eta","#eta(ll)", 30, -5, 5);
  ll_phi = mkTH1F("ll_phi","#phi(ll)", 30, -3.15, 3.15);
  ll_m = mkTH1F("ll_m","m(ll)", 60, 60, 120);

  bb_pT = mkTH1F("bb_pT","p_{T}(bb)", 30, 0, 500);
  bb_eta = mkTH1F("bb_eta","#eta(bb)", 30, -5, 5);
  bb_phi = mkTH1F("bb_phi","#phi(bb)", 30, -3.15, 3.15);
  bb_m = mkTH1F("bb_m","m(bb)", 1800, 100, 1000);

  llbb_pT = mkTH1F("llbb_pT","p_{T}(llbb)", 30, 0, 500);
  llbb_eta = mkTH1F("llbb_eta","#eta(llbb)", 30, -5, 5);
  llbb_phi = mkTH1F("llbb_phi","#phi(llbb)", 30, -3.15, 3.15);
  llbb_m = mkTH1F("llbb_m","m(llbb)", 1600, 200, 1000);

  dEta_ll = mkTH1F("dEta_ll","#Delta#eta(l,l)", 30, 0, 5);
  dPhi_ll = mkTH1F("dPhi_ll","#Delta#phi(l,l)", 30, 0, 3.15);
  dR_ll = mkTH1F("dR_ll","#DeltaR(l,l)", 30, 0, 10);
  dEta_bb = mkTH1F("dEta_bb","#Delta#eta(b,b)", 30, 0, 5);
  dPhi_bb = mkTH1F("dPhi_bb","#Delta#phi(b,b)",30, 0, 3.15);
  dR_bb = mkTH1F("dR_bb","#DeltaR(b,b)", 30, 0, 10);

  cosTheta_l_Z = mkTH1F("cosTheta_l_Z","cos#theta^{*}(l)", 30, 0, 1);
  cosTheta_b_H = mkTH1F("cosTheta_b_H","cos#theta^{*}(b)", 30, 0, 1);

  // register
  //for_each( myhist.begin(), myhist.end(), std::bind1st( std::mem_fun( &readtruth_azheavyh::mbook ), this ) );
  // now done in mkTH1F

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: initialize ()
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
  ctr_zee = 0;
  ctr_zmm = 0;
  ctr_ztt = 0;
  ctr_tot = 0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  ANA_CHECK_SET_TYPE( EL::StatusCode );

  xAOD::TEvent* event = wk()->xaodEvent();


  // my init
  lepFlavor = 0;
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

//  if( weight<0 ){
//  }else{
//    return EL::StatusCode::SUCCESS;
//  }

  // my collection of truth particles
  std::map<std::string, const xAOD::TruthParticle*> mycollection;
  typedef std::pair<std::string, const xAOD::TruthParticle*> PartType;

  for( unsigned ipart=0; ipart<npart; ++ipart ){

    const xAOD::TruthParticle* particle = mcpart->at( ipart );

    // A
    if( particle->pdgId() == 36 ){
	  if( mycollection.count("A") == 0 ) mycollection["A"] = particle; // correctedParticle( particle );
	}

	// H
	if( particle->pdgId() == 35 && hasParent( particle, 36 ) ){
	  if( mycollection.count("H") == 0 ) mycollection["H"] = particle; // correctedParticle( particle );
	}

	// Z
	if( particle->pdgId() == 23 && hasParent( particle, 36) ){
	  if( mycollection.count("Z") == 0 ) mycollection["Z"] = particle; // correctedParticle( particle );
	}

	// e- and tau->e-
	//if( particle->pdgId() == 11 && ( hasParent( particle, 23) || hasParent( particle, 15 ) ) ){
	// e- only
	if( particle->pdgId() == 11 && hasParent( particle, 23) ){
	  if( mycollection.count("e-") == 0 ) mycollection["e-"] = particle; // correctedParticle( particle );
	}

	// e+ and tau->e+
	//if( particle->pdgId() == -11 && ( hasParent( particle, 23) || hasParent( particle, -15) ) ){
	// e+ only
	if( particle->pdgId() == -11 && hasParent( particle, 23) ){
	  if( mycollection.count("e+") == 0 ) mycollection["e+"] = particle; // correctedParticle( particle );
	}

	// mu- and tau->mu-
	//if( particle->pdgId() == 13 && ( hasParent( particle, 23) || hasParent( particle, 15) ) ){
	// mu- only
	if( particle->pdgId() == 13 && hasParent( particle, 23) ){
	  if( mycollection.count("mu-") == 0 ) mycollection["mu-"] = particle; // correctedParticle( particle );
	}

	// mu+ and tau->mu+
	//if( particle->pdgId() == -13 && ( hasParent( particle, 23) || hasParent( particle, -15) ) ){
	// mu+ only
	if( particle->pdgId() == -13 && hasParent( particle, 23) ){
	  if( mycollection.count("mu+") == 0 ) mycollection["mu+"] = particle; // correctedParticle( particle );
	}

	// b
	if( particle->pdgId() == 5 && hasParent( particle, 35) ){
	  if( mycollection.count("b") == 0 ){
	    mycollection["b"] = particle; // = correctedParticle( particle );
		//printChildren( mycollection["b"] );
	  }
	}

	// b~
	if( particle->pdgId() == -5 && hasParent( particle, 35) ){
	  if( mycollection.count("b~") == 0 ) mycollection["b~"] = particle; //= correctedParticle( particle );
	}

	if( mycollection.count("A")
	 && mycollection.count("H")
	 && mycollection.count("Z")
	 && mycollection.count("e-")
	 && mycollection.count("e+")
	 && mycollection.count("mu-")
	 && mycollection.count("mu+")
	 && mycollection.count("b")
	 && mycollection.count("b~")
	)
	  break;
  }


  // check check check !!!
  // do not fill anything before check!!!
  //
  // check there is one A/H/Z
  if( mycollection.count("A") == 1); else{ Error("execute()","missing A, skip event"); return EL::StatusCode::SUCCESS;}
  if( mycollection.count("H") == 1); else{ Error("execute()","missing H, skip event"); return EL::StatusCode::SUCCESS;}
  if( mycollection.count("Z") == 1); else{ Error("execute()","missing Z, skip event"); return EL::StatusCode::SUCCESS;}
  // following is not a check any more
  // drop counting tautau
  if( ( mycollection.count("e-") + mycollection.count("e+") == 2 ) ){
    lepFlavor = 11;
  }
  if( ( mycollection.count("mu-") + mycollection.count("mu+") == 2 ) ){
    lepFlavor = 13;
  }
  if( lepFlavor != 11 && lepFlavor != 13 ){
    return EL::StatusCode::SUCCESS; // skip tautau
  }

  // counter of Z decay, leptonic tau is incuded, had tau is droped in counting
  ctr_zee += (lepFlavor == 11);
  ctr_zmm += (lepFlavor == 13);
  //ctr_ztt += (lepFlavor == 15);

  // fill the histograms

  // bosons
  A_pT->Fill( mycollection["A"]->pt()*GeV, weight );
  A_eta->Fill( mycollection["A"]->eta(), weight );
  A_phi->Fill( mycollection["A"]->phi(), weight );
  A_E->Fill( mycollection["A"]->e()*GeV, weight );
  A_m->Fill( mycollection["A"]->m()*GeV, weight );

  H_pT->Fill( mycollection["H"]->pt()*GeV, weight );
  H_eta->Fill( mycollection["H"]->eta(), weight );
  H_phi->Fill( mycollection["H"]->phi(), weight );
  H_E->Fill( mycollection["H"]->e()*GeV, weight );
  H_m->Fill( mycollection["H"]->m()*GeV, weight );

  Z_pT->Fill( mycollection["Z"]->pt()*GeV, weight );
  Z_eta->Fill( mycollection["Z"]->eta(), weight );
  Z_phi->Fill( mycollection["Z"]->phi(), weight );
  Z_E->Fill( mycollection["Z"]->e()*GeV, weight );
  Z_m->Fill( mycollection["Z"]->m()*GeV, weight );

  ZH_m->Fill( (mycollection["Z"]->p4() + mycollection["H"]->p4()).M()*GeV, weight );

  // leptons
  std::string lep1 = "", lep2 = "";
  if( lepFlavor == 11 ) { lep1 = "e-"; lep2 = "e+";}
  if( lepFlavor == 13 ) { lep1 = "mu-"; lep2 = "mu+";}
  //if( lepFlavor == 15 ) { lep1 = "tau0"; lep2 = "tau1";}
  // sort lepton by pT
  if( mycollection[lep1]->pt() < mycollection[lep2]->pt() ){
    std::string _lep = lep1;
	lep1 = lep2;
	lep2 = _lep;
  }
  l1_pT->Fill( mycollection[lep1]->pt()*GeV, weight );
  l1_eta->Fill( mycollection[lep1]->eta(), weight );
  l1_phi->Fill( mycollection[lep1]->phi(), weight );
  l1_E->Fill( mycollection[lep1]->e()*GeV, weight );

  l2_pT->Fill( mycollection[lep2]->pt()*GeV, weight );
  l2_eta->Fill( mycollection[lep2]->eta(), weight );
  l2_phi->Fill( mycollection[lep2]->phi(), weight );
  l2_E->Fill( mycollection[lep2]->e()*GeV, weight );

  // b quarks (NOT jets)!
  std::string bjet1 = "b", bjet2 = "b~";
  // sort bjets by pt
  if( mycollection[bjet1]->pt() < mycollection[bjet2]->pt() ){
    std::string _bjet = bjet1;
	bjet1 = bjet2;
	bjet2 = _bjet;
  }
  b1_pT->Fill( mycollection[bjet1]->pt()*GeV, weight );
  b1_eta->Fill( mycollection[bjet1]->eta(), weight );
  b1_phi->Fill( mycollection[bjet1]->phi(), weight );
  b1_E->Fill( mycollection[bjet1]->e()*GeV, weight );

  b2_pT->Fill( mycollection[bjet2]->pt()*GeV, weight );
  b2_eta->Fill( mycollection[bjet2]->eta(), weight );
  b2_phi->Fill( mycollection[bjet2]->phi(), weight );
  b2_E->Fill( mycollection[bjet2]->e()*GeV, weight );

  // system var
  TLorentzVector v_lep1 = mycollection[lep1]->p4();
  TLorentzVector v_lep2 = mycollection[lep2]->p4();
//  TLorentzVector v_bjet1 = mycollection[bjet1]->p4();
//  TLorentzVector v_bjet2 = mycollection[bjet2]->p4();
  TLorentzVector v_bjet1 = getJet(mycollection[bjet1]);
  TLorentzVector v_bjet2 = getJet(mycollection[bjet2]);

  TLorentzVector v_ll = v_lep1 + v_lep2;
  ll_pT->Fill( v_ll.Pt()*GeV, weight );
  ll_eta->Fill( v_ll.Eta(), weight );
  ll_phi->Fill( v_ll.Phi(), weight );
  ll_m->Fill( v_ll.M()*GeV, weight );

  TLorentzVector v_bb = v_bjet1 + v_bjet2;
  bb_pT->Fill( v_bb.Pt()*GeV, weight );
  bb_eta->Fill( v_bb.Eta(), weight );
  bb_phi->Fill( v_bb.Phi(), weight );
  bb_m->Fill( v_bb.M()*GeV, weight );

  TLorentzVector v_llbb = v_lep1 + v_lep2 + v_bjet1 + v_bjet2;
  llbb_pT->Fill( v_llbb.Pt()*GeV, weight );
  llbb_eta->Fill( v_llbb.Eta(), weight );
  llbb_phi->Fill( v_llbb.Phi(), weight );
  llbb_m->Fill( v_llbb.M()*GeV, weight );

  // angular var
  dEta_ll->Fill( std::abs( v_lep1.Eta() - v_lep2.Eta() ), weight);
  dPhi_ll->Fill( v_lep1.DeltaPhi(v_lep2), weight);
  dR_ll->Fill( v_lep1.DeltaR(v_lep2), weight);
  dEta_bb->Fill( std::abs( v_bjet1.Eta() - v_bjet2.Eta() ), weight);
  dPhi_bb->Fill( v_bjet1.DeltaPhi(v_bjet2), weight);
  dR_bb->Fill( v_bjet1.DeltaR(v_bjet2), weight);

  TLorentzVector v_lep1_sysll = v_lep1;
  TLorentzVector v_lep2_sysll = v_lep2;
  v_lep1_sysll.Boost( -(v_ll).BoostVector() );
  v_lep2_sysll.Boost( -(v_ll).BoostVector() );
  double angle_l1 = v_lep1_sysll.Angle(v_ll.Vect());
  double angle_l2 = v_lep2_sysll.Angle(v_ll.Vect());
  double angle_l = angle_l1 < angle_l2 ? angle_l1: angle_l2;
  cosTheta_l_Z->Fill( std::cos(angle_l), weight );

  TLorentzVector v_bjet1_sysbb = v_bjet1;
  TLorentzVector v_bjet2_sysbb = v_bjet2;
  v_bjet1_sysbb.Boost( -(v_bb).BoostVector() );
  v_bjet2_sysbb.Boost( -(v_bb).BoostVector() );
  double angle_b1 = v_bjet1_sysbb.Angle(v_bb.Vect());
  double angle_b2 = v_bjet2_sysbb.Angle(v_bb.Vect());
  double angle_b = angle_b1 < angle_b2 ? angle_b1 : angle_b2;
  cosTheta_b_H->Fill( std::cos(angle_b), weight );


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: finalize ()
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

  // Z decay fractions
  double ctr_sum = ctr_zee+ctr_zmm+ctr_ztt;
  Info("finalize()", "Z decay fractions ee = %f, mm = %f, tt = %f, based on collected",
       ctr_zee / ctr_sum, ctr_zmm/ctr_sum, ctr_ztt/ctr_sum );
  Info("finalize()", "Total events: %d, truth collected events: %d", ctr_tot, int(ctr_sum) );

  // weight
  Info("finalize()", "# evt with positive weight: %ld", ctr_posw);
  Info("finalize()", "# evt with negative weight: %ld", ctr_negw);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_azheavyh :: histFinalize ()
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
