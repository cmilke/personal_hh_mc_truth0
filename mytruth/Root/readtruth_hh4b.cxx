#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <mytruth/readtruth_hh4b.h>

// this is needed to distribute the algorithm to the workers
ClassImp(readtruth_hh4b)

readtruth_hh4b :: readtruth_hh4b ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}

TH1F* readtruth_hh4b :: mkTH1F( const char* n, const char* t, double nb, double x0, double x1){
  TH1F* h = new TH1F( n, t, nb, x0, x1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

TH2F* readtruth_hh4b::mkTH2F( const char* n, const char* t, double xnb, double x0, double x1, double ynb, double y0, double y1 ){
  TH2F* h = new TH2F( n, t, xnb, x0, x1, ynb, y0, y1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

// in principal we should look at the particlet in the end of self-decay chains
const xAOD::TruthParticle* readtruth_hh4b :: correctedParticle( const xAOD::TruthParticle * part ){
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
const xAOD::TruthParticle* readtruth_hh4b :: correctedQuark( const xAOD::TruthParticle * part ){
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


bool readtruth_hh4b :: hasParent( const xAOD::TruthParticle* kid, int parentId ){
  for( unsigned int iparent=0; iparent<kid->nParents(); ++iparent ){
    if( kid->parent(iparent)->pdgId() == parentId ) return true;
  }
  return false;
}

void readtruth_hh4b :: printChildren( const xAOD::TruthParticle* part ){
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
TLorentzVector readtruth_hh4b :: getJet( const xAOD::TruthParticle* part ){
  TLorentzVector jet;
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid ){
    jet += part->child(ikid)->p4();
  }
  return jet;
}
*/

bool readtruth_hh4b :: larger_pT( const xAOD::TruthParticle* p1, const xAOD::TruthParticle* p2 ){
  return p1->pt() > p2->pt();
}

EL::StatusCode readtruth_hh4b :: setupJob (EL::Job& job)
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



EL::StatusCode readtruth_hh4b :: histInitialize ()
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

  b3_pT = mkTH1F("b3_pT","p_{T}(b quark 3)", 50, 0, 800);
  b3_eta = mkTH1F("b3_eta","#eta(b quark 3)", 40, -5, 5);
  b3_phi = mkTH1F("b3_phi","#phi(b quark 3)", 40, -3.15, 3.15 );
  b3_E = mkTH1F("b3_E","E(b quark 3)", 50, 0, 2000);

  b4_pT = mkTH1F("b4_pT","p_{T}(b quark 4)", 50, 0, 800);
  b4_eta = mkTH1F("b4_eta","#eta(b quark 4)", 40, -5, 5);
  b4_phi = mkTH1F("b4_phi","#phi(b quark 4)", 40, -3.15, 3.15 );
  b4_E = mkTH1F("b4_E","E(b quark 4)", 50, 0, 2000);

  // bb1 for H1; bb2 for H2
  bb1_pT = mkTH1F("bb1_pT","p_{T}(bb from H1)", 50, 0, 800);
  bb1_eta = mkTH1F("bb1_eta","#eta(bb from H1)", 40, -5, 5);
  bb1_phi = mkTH1F("bb1_phi","#phi(bb from H1)", 40, -3.15, 3.15 );
  bb1_E = mkTH1F("bb1_E","E(bb from H1)", 50, 0, 2000);
  bb1_m = mkTH1F("bb1_m","m(bb from H1)", 40, 125-50, 125+50 );

  bb2_pT = mkTH1F("bb2_pT","p_{T}(bb from H2)", 50, 0, 800);
  bb2_eta = mkTH1F("bb2_eta","#eta(bb from H2)", 40, -5, 5);
  bb2_phi = mkTH1F("bb2_phi","#phi(bb from H2)", 40, -3.15, 3.15 );
  bb2_E = mkTH1F("bb2_E","E(bb from H2)", 50, 0, 2000);
  bb2_m = mkTH1F("bb2_m","m(bb from H2)", 40, 125-50, 125+50 );

  // bbbb for HH
  bbbb_pT = mkTH1F("bbbb_pT","p_{T}(bbbb)", 1000, 0, 1000); //50 bins
  bbbb_eta = mkTH1F("bbbb_eta","#eta(bbbb)", 40, -10, 10);
  bbbb_phi = mkTH1F("bbbb_phi","#phi(bbbb)", 40, -3.15, 3.15 );
  bbbb_E = mkTH1F("bbbb_E","E(bbbb)", 100, 0, 4000);
  bbbb_m = mkTH1F("bbbb_m","m(bbbb)", 3000, 0, 3000);

  dEta_bb_bb = mkTH1F("dEta_bb_bb","#Delta#eta(bb,bb)", 40, 0, 8);
  dPhi_bb_bb = mkTH1F("dPhi_bb_bb","#Delta#phi(bb,bb)", 40, 0, 3.15);
  dR_bb_bb = mkTH1F("dR_bb_bb","#DeltaR(bb,bb)", 40, 0, 10);

  // register
  //for_each( myhist.begin(), myhist.end(), std::bind1st( std::mem_fun( &readtruth_hh4b::mbook ), this ) );
  // now done in mkTH1F

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: initialize ()
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

  // test/count bbyy/bbZy
  ctr_bbyy = 0;
  ctr_bbZy = 0;

  ctr_save = 0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: execute ()
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

  //std::cout << "DDD MC weight vector size: " << mcevt->at(0)->weights().size() << std::endl;

//  if( weight<0 ){
//  }else{
//    return EL::StatusCode::SUCCESS;
//  }

  // my collection of truth particles
  std::map<std::string, const xAOD::TruthParticle*> mycollection;
  std::vector<const xAOD::TruthParticle*> mybquarks_H1; // 2 b quarks from H1
  std::vector<const xAOD::TruthParticle*> mybquarks_H2; // 2 b quarks from H2
  std::vector<const xAOD::TruthParticle*> mybquarks; // all 4 b quarks ordere by pT
  typedef std::pair<std::string, const xAOD::TruthParticle*> PartType;

  for( unsigned ipart=0; ipart<npart; ++ipart ){

    const xAOD::TruthParticle* particle = mcpart->at( ipart );

	// debug
	if( particle->pdgId() == 25 ){
	  Info( "execute()", "ID %d, barcode %d, parent ID %d, pT %f, px %f, py %f, pz %f, eta %f, phi %f",
	       particle->pdgId(), particle->barcode(), particle->parent(0)->pdgId(),
	       particle->pt(), particle->px(), particle->py(), particle->pz(), particle->eta(), particle->phi() );
	}

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
	  Info( "execute()", "SAVE ID %d, barcode %d, parent ID %d, pT %f, px %f, py %f, pz %f, eta %f, phi %f",
	       particle->pdgId(), particle->barcode(), particle->parent(0)->pdgId(),
	       particle->pt(), particle->px(), particle->py(), particle->pz(), particle->eta(), particle->phi() );
	}

	// once found H1 H2 quit loop (only used without looking into H decays)
	if( mycollection.count("H1")==1 && mycollection.count("H2")==1 )
	  break;

  }

  // check
  if( mycollection.count("H1")==1 && mycollection.count("H2")==1 ){
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

/*
  // test/count bbyy bbZy events
  // from H1
  int _c_b = 0;
  int _c_y = 0;
  int _c_Z = 0;
  for( unsigned int ikid=0; ikid<mycollection["H1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H1"]->child(ikid);
	if( abs(kid->pdgId()) == 5 ){
	  ++_c_b;
	}
	if( abs(kid->pdgId()) == 22 ){
	  ++_c_y;
	}
	if( abs(kid->pdgId()) == 23 ){
	  ++_c_Z;
	}
  }
  // from H2
  for( unsigned int ikid=0; ikid<mycollection["H2"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H2"]->child(ikid);
	if( abs(kid->pdgId()) == 5 ){
	  ++_c_b;
	}
	if( abs(kid->pdgId()) == 22 ){
	  ++_c_y;
	}
	if( abs(kid->pdgId()) == 23 ){
	  ++_c_Z;
	}

  }
  //
  if( (_c_b==2)&&(_c_y==1)&&(_c_Z==1) ){
    ++ctr_bbZy;
  }
  if( (_c_b==2)&&(_c_y==2) ){
    ++ctr_bbyy;
  }
*/


  // original codes for 4b final stats

  // find b quarks originating from these 2 Higgs bosons
  // from H1
  for( unsigned int ikid=0; ikid<mycollection["H1"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H1"]->child(ikid);
	if( abs(kid->pdgId()) == 5 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mybquarks_H1.begin(), mybquarks_H1.end(), kid) == mybquarks_H1.end() ){
	    mybquarks_H1.push_back(kid);
	  }
	}
	if( mybquarks_H1.size()==2 ) break;
  }
  // from H2
  for( unsigned int ikid=0; ikid<mycollection["H2"]->nChildren(); ++ikid){
    const xAOD::TruthParticle* kid = mycollection["H2"]->child(ikid);
	if( abs(kid->pdgId()) == 5 ){
	  kid = correctedParticle( kid ); // after shower dressing
	  if( std::find(mybquarks_H2.begin(), mybquarks_H2.end(), kid) == mybquarks_H2.end() ){
	    mybquarks_H2.push_back(kid);
	  }
	}
	if( mybquarks_H2.size()==2 ) break;
  }
  // all 4 b quarks ordered by pT
  mybquarks.reserve( mybquarks_H1.size() + mybquarks_H2.size() );
  mybquarks.insert( mybquarks.end(), mybquarks_H1.begin(), mybquarks_H1.end() );
  mybquarks.insert( mybquarks.end(), mybquarks_H2.begin(), mybquarks_H2.end() );
  std::sort( mybquarks.begin(), mybquarks.end(), larger_pT );

  //test
  //std::cout << "sorted b quarks" << std::endl;
  //for( auto obj : mybquarks ) std::cout << "b quark: " << obj << " pT: " << obj->pt() << std::endl;

  // check check check !!!
  // do not fill anything before check!!!
  //
  // check there is one A/H/Z
  //if( mycollection.count("A") == 1); else{ Error("execute()","missing A, skip event"); return EL::StatusCode::SUCCESS;}
  //if( mycollection.count("H") == 1); else{ Error("execute()","missing H, skip event"); return EL::StatusCode::SUCCESS;}
  //if( mycollection.count("Z") == 1); else{ Error("execute()","missing Z, skip event"); return EL::StatusCode::SUCCESS;}
  // following is not a check any more
  // drop counting tautau
  //if( ( mycollection.count("e-") + mycollection.count("e+") == 2 ) ){
  //  lepFlavor = 11;
  //}
  //if( ( mycollection.count("mu-") + mycollection.count("mu+") == 2 ) ){
  //  lepFlavor = 13;
  //}
  //if( lepFlavor != 11 && lepFlavor != 13 ){
  //  return EL::StatusCode::SUCCESS; // skip tautau
  //}

  // counter of Z decay, leptonic tau is incuded, had tau is droped in counting
  //ctr_zee += (lepFlavor == 11);
  //ctr_zmm += (lepFlavor == 13);
  //ctr_ztt += (lepFlavor == 15);

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

  b3_pT->Fill( mybquarks[2]->pt()*GeV, weight );
  b3_eta->Fill( mybquarks[2]->eta(), weight );
  b3_phi->Fill( mybquarks[2]->phi(), weight );
  b3_E->Fill( mybquarks[2]->e()*GeV, weight );

  b4_pT->Fill( mybquarks[3]->pt()*GeV, weight );
  b4_eta->Fill( mybquarks[3]->eta(), weight );
  b4_phi->Fill( mybquarks[3]->phi(), weight );
  b4_E->Fill( mybquarks[3]->e()*GeV, weight );

  // bb1 for H1; bb2 for H2
  TLorentzVector v_bb1 = mybquarks_H1[0]->p4() + mybquarks_H1[1]->p4();
  TLorentzVector v_bb2 = mybquarks_H2[0]->p4() + mybquarks_H2[1]->p4();

  bb1_pT->Fill( v_bb1.Pt()*GeV, weight );
  bb1_eta->Fill( v_bb1.Eta(), weight );
  bb1_phi->Fill( v_bb1.Phi(), weight );
  bb1_E->Fill( v_bb1.E()*GeV, weight );
  bb1_m->Fill( v_bb1.M()*GeV, weight );

  bb2_pT->Fill( v_bb2.Pt()*GeV, weight );
  bb2_eta->Fill( v_bb2.Eta(), weight );
  bb2_phi->Fill( v_bb2.Phi(), weight );
  bb2_E->Fill( v_bb2.E()*GeV, weight );
  bb2_m->Fill( v_bb2.M()*GeV, weight );

  // bbbb for HH
  TLorentzVector v_bbbb = v_bb1 + v_bb2;
  bbbb_pT->Fill( v_bbbb.Pt()*GeV, weight );
  bbbb_eta->Fill( v_bbbb.Eta(), weight );
  bbbb_phi->Fill( v_bbbb.Phi(), weight );
  bbbb_E->Fill( v_bbbb.E()*GeV, weight );
  bbbb_m->Fill( v_bbbb.M()*GeV, weight );

  // angular
  dEta_bb_bb->Fill( std::abs( v_bb1.Eta() - v_bb2.Eta() ), weight);
  dPhi_bb_bb->Fill( v_bb1.DeltaPhi(v_bb2), weight);
  dR_bb_bb->Fill( v_bb1.DeltaR(v_bb2), weight);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: finalize ()
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

  Error("finalize()", "Total events: %d, truth collected events: %d", ctr_tot, int(ctr_save) );

  // test/count bbyy/bbZy
  //Error("finalize()", "Only for bbyy/bbZy events: bbyy %d, bbZy %d", ctr_bbyy, ctr_bbZy );

  // weight
  Error("finalize()", "# evt with positive weight: %ld", ctr_posw);
  Error("finalize()", "# evt with negative weight: %ld", ctr_negw);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_hh4b :: histFinalize ()
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
