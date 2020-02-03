#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <mytruth/readtruth_vbfhh.h>

// this is needed to distribute the algorithm to the workers
ClassImp(readtruth_vbfhh)

readtruth_vbfhh :: readtruth_vbfhh ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}

TH1F* readtruth_vbfhh :: mkTH1F( const char* n, const char* t, double nb, double x0, double x1){
  TH1F* h = new TH1F( n, t, nb, x0, x1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

TH2F* readtruth_vbfhh::mkTH2F( const char* n, const char* t, double xnb, double x0, double x1, double ynb, double y0, double y1 ){
  TH2F* h = new TH2F( n, t, xnb, x0, x1, ynb, y0, y1 );
  h->SetStats(0);
  h->Sumw2(1);
  myhist.push_back(h);
  wk()->addOutput(h); // cannot put mkTH1F as inline function due to wk() incomplete class EL::Worker
  return h;
}

// in principal we should look at the particlet in the end of self-decay chains
const xAOD::TruthParticle* readtruth_vbfhh :: correctedParticle( const xAOD::TruthParticle * part ){
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
const xAOD::TruthParticle* readtruth_vbfhh :: correctedQuark( const xAOD::TruthParticle * part ){
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


bool readtruth_vbfhh :: hasParent( const xAOD::TruthParticle* kid, int parentId ){
  for( unsigned int iparent=0; iparent<kid->nParents(); ++iparent ){
    if( kid->parent(iparent)->pdgId() == parentId ) return true;
  }
  return false;
}

void readtruth_vbfhh :: printChildren( const xAOD::TruthParticle* part ){
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
TLorentzVector readtruth_vbfhh :: getJet( const xAOD::TruthParticle* part ){
  TLorentzVector jet;
  for( unsigned int ikid=0; ikid<part->nChildren(); ++ikid ){
    jet += part->child(ikid)->p4();
  }
  return jet;
}
*/

bool readtruth_vbfhh :: larger_pT( const xAOD::TruthParticle* p1, const xAOD::TruthParticle* p2 ){
  return p1->pt() > p2->pt();
}

EL::StatusCode readtruth_vbfhh :: setupJob (EL::Job& job)
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



EL::StatusCode readtruth_vbfhh :: histInitialize ()
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

  dR_qq = mkTH1F("dR_qq","#DeltaR(q,q)", 40, 0, 10);
  Num_q = mkTH1F("Num_q","Number of q", 10, 0, 10);

  // register
  //for_each( myhist.begin(), myhist.end(), std::bind1st( std::mem_fun( &readtruth_vbfhh::mbook ), this ) );
  // now done in mkTH1F

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_vbfhh :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_vbfhh :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_vbfhh :: initialize ()
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



EL::StatusCode readtruth_vbfhh :: execute ()
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
  //    bool flag_q=false;

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
  std::vector<const xAOD::TruthParticle*> myquarks; //addiional quarks from VBF mode
  typedef std::pair<std::string, const xAOD::TruthParticle*> PartType;

  for( unsigned ipart=0; ipart<npart; ++ipart ){

    const xAOD::TruthParticle* particle = mcpart->at( ipart );
//    bool flag_q=false;

	// debug
	if( particle->pdgId() == 25 ){
	  Info( "execute()", "ID %d, barcode %d, parent ID %d, pT %f, px %f, py %f, pz %f, eta %f, phi %f",
	       particle->pdgId(), particle->barcode(), particle->parent(0)->pdgId(),
	       particle->pt(), particle->px(), particle->py(), particle->pz(), particle->eta(), particle->phi() );
	}


	//find VBF quarks
         if(particle->hasProdVtx()==1 && particle->prodVtx()->barcode()==-1 && abs(particle->pdgId())<6){
                particle = correctedParticle( particle );
                myquarks.push_back(particle);
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
  TLorentzVector v_q1, v_q2;
  if(myquarks.size()>0) v_q1.SetPtEtaPhiE(myquarks[0]->pt(),myquarks[0]->eta(),myquarks[0]->phi(),myquarks[0]->e());
  if(myquarks.size()>1) v_q2.SetPtEtaPhiE(myquarks[1]->pt(),myquarks[1]->eta(),myquarks[1]->phi(),myquarks[1]->e());
 

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

  if(myquarks.size()>1) dR_qq->Fill(v_q1.DeltaR(v_q2),weight);
  Num_q->Fill(myquarks.size());

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_vbfhh :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode readtruth_vbfhh :: finalize ()
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



EL::StatusCode readtruth_vbfhh :: histFinalize ()
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
