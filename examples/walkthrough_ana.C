#ifndef __CINT__
#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCGenericObject.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCRelation.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCIterator.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#endif

#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TText.h"
#include "TFile.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>

/** Example script for a simple cut-based selection
 *
 see https://indico.fnal.gov/event/45721 for detailed instructions and explanations
 *
 */


/*
  put this into your .rootlogon.C file

  {
  gInterpreter->AddIncludePath("$LCIO");
  gSystem->Load("${LCIO}/lib/liblcio.so");
  gSystem->Load("${LCIO}/lib/liblcioDict.so");
  }

  for the LCIO API documentation see:
  http://lcio.desy.de/v02-09/doc/doxygen_api/html/index.html
*/


using namespace lcio ;

const TString outPlotName="walkthrough"; // name of output files

// what type of miniDST are we running over?
const int minDstFlavour = 0; // 0=DELPHES, 1=SGV, 2=fullsim
  
// location of input data

// this is for delphes on loginsnowmass21.io
const TString inputDir_SIGNAL = "/collab/project/snowmass21/data/ilc/analysis-walkthrough/signal/miniDST-delphes/";
const TString inputDir_BACKGD = "/collab/project/snowmass21/data/ilc/analysis-walkthrough/backgrounds/";
//const TString inputDir_SIGNAL = "/collab/project/snowmass21/data/ilc/analysis-walkthrough/signal/miniDST-SGV/"; // for the SGV signal files

// on Daniel's machine
//const TString inputDir_SIGNAL = "./signalSamples/miniDST-delphes/";
//const TString inputDir_SIGNAL = "./signalSamples/miniDST-SGV/";
//const TString inputDir_BACKGD = "./";

template<class T>
double inv_mass(T* p1, T* p2){
  double e = p1->getEnergy()+p2->getEnergy() ;
  double px = p1->getMomentum()[0]+p2->getMomentum()[0];
  double py = p1->getMomentum()[1]+p2->getMomentum()[1];
  double pz = p1->getMomentum()[2]+p2->getMomentum()[2];
  return( sqrt( e*e - px*px - py*py - pz*pz  ) );
}

template<class T>
TLorentzVector v4(T* p){
  return TLorentzVector( p->getMomentum()[0],p->getMomentum()[1], p->getMomentum()[2],p->getEnergy());
}

TLorentzVector v4Jet(ReconstructedParticle* p){
  TLorentzVector tlj(0,0,0,0);
  if ( p->getMomentum()[0]==p->getMomentum()[0] ) {
    tlj=v4(p);
  } else { // there's a NAN, need to fudge/recalculate (this is a current "bug" of delphes2LCIO)
    float extraE(0);
    for (size_t j=0; j<p->getParticles().size(); j++) {
      if ( p->getParticles()[j]->getMomentum()[0] == p->getParticles()[j]->getMomentum()[0] ) {
        tlj += v4(p->getParticles()[j]);
      } else { // problematic PFO
        extraE+=p->getParticles()[j]->getEnergy();
      }
    }
    if ( extraE>0 ) { // assume this extra energy is collinear with rest of jet, and that its massless
      TVector3 extraMom=tlj.Vect();
      extraMom*=1/extraMom.Mag();
      extraMom*=extraE;
      tlj+=TLorentzVector(extraMom, extraE);
    }
  }
  return tlj;
}



// maximum events to run over for a given process
int maxEvt; 
TFile* _fout;

std::vector <TH1F*> analyse_process( TString nickname, std::vector <TString> fnames, int minDstFlavour ) {

  // this runs over a list of miniDST file names, 
  //  does some event analysis and selection
  //   fills some histograms (which it returns)
  //
  // "nickname" used in histogram names/titles
  // "fnames" input file names; should all be from the *same* process
  // "minDstFlavour" = 0 (for DELPHES); 1 (SGV) [potentially 2 (fullsim) -> not yet developed, an exercise for the user]

  _fout->mkdir(nickname)->cd();

  std::vector <TH1F*> histos;

  std::string muoColName = "IsolatedMuons";
  std::string jetColName, pfoColName, mcpColName;
  std::string pidName;

  if ( minDstFlavour==0 ) { // delphes-miniDST
    jetColName="Durham2Jets";
    pfoColName="PFOs";
    mcpColName="MCParticles";
    pidName="JetParameters";
  } else if ( minDstFlavour==1 ) { // SGV-miniDST
    jetColName="Refined2Jets";
    pfoColName="PandoraPFOs";
    mcpColName="MCParticlesSkimmed";
    pidName="lcfiplus";
  } else {
    cout << "unknown miniDST flavour!" << endl;
    return histos;
  }

  int nEvents = 0; // event counter

  //---- create some histograms ----
  TH1F* ALL_sel      = new TH1F(nickname+"_ALL_sel", "(ALL) "+nickname+" sel", 20,0,20);
  TH1F* ALL_nmuons   = new TH1F(nickname+"_ALL_nmuons", "(ALL) "+nickname+" nmuons", 5,0,5);

  TH1F* MC_muonmass   = new TH1F(nickname+"_MC_muonmass"  ,"(MC) "+nickname+" inv. mass - muons", 100,  60. , 120. ) ;
  TH1F* MC_recoilm    = new TH1F(nickname+"_MC_recoilm"   ,"(MC) "+nickname+" recoil mass", 100, 70. , 170. ) ;
  TH1F* MC_missm      = new TH1F(nickname+"_MC_missm"     ,"(MC) "+nickname+" missing mass", 300, -150. , 150. ) ;
  TH1F* MC_costhmiss  = new TH1F(nickname+"_MC_costhmiss" ,"(MC) "+nickname+" cosTheta pMiss", 100, -1. , 1. ) ;


  TH1F* PRE_btag        = new TH1F(nickname+"_PRE_btag",       "(PRESEL) "+nickname+" btag", 50,-0.2,1.2);
  TH1F* PRE_njetChg     = new TH1F(nickname+"_PRE_njetchg",    "(PRESEL) "+nickname+" njetchg", 25,0,25);

  TH1F* PRE_muonmass   = new TH1F(nickname+"_PRE_muonmass"  ,"(PRESEL) "+nickname+" inv. mass - muons", 100,  60. , 120. ) ;
  TH1F* PRE_jetmass    = new TH1F(nickname+"_PRE_jetmass"   ,"(PRESEL) "+nickname+" inv. mass - jets", 100, 0. , 150. ) ;
  TH1F* PRE_recoilm    = new TH1F(nickname+"_PRE_recoilm"   ,"(PRESEL) "+nickname+" recoil mass", 100, 70. , 170. ) ;
  TH1F* PRE_missm      = new TH1F(nickname+"_PRE_missm"     ,"(PRESEL) "+nickname+" missing mass", 300, -150. , 150. ) ;
  TH1F* PRE_missm_jjm  = new TH1F(nickname+"_PRE_missm_jjm" ,"(PRESEL) "+nickname+" miss-jj mass difference", 100, -150. , 150. ) ;
  TH1F* PRE_costhmiss  = new TH1F(nickname+"_PRE_costhmiss" ,"(PRESEL) "+nickname+" cosTheta pMiss", 100, -1. , 1. ) ;
  TH1F* PRE_emiss      = new TH1F(nickname+"_PRE_emiss"     ,"(PRESEL) "+nickname+" missing energy", 100, 0. , 150. ) ;
  TH1F* PRE_nBtag       = new TH1F(nickname+"_PRE_nBtag"      ,"(PRESEL) "+nickname+" nBtag", 3, 0. , 3. ) ;

  TH1F* SEL_muonmass   = new TH1F(nickname+"_SEL_muonmass"  ,"(SEL) "+nickname+" inv. mass - muons", 100,  60. , 120. ) ;
  TH1F* SEL_jetmass    = new TH1F(nickname+"_SEL_jetmass"   ,"(SEL) "+nickname+" inv. mass - jets", 100, 0. , 150. ) ;
  TH1F* SEL_recoilm    = new TH1F(nickname+"_SEL_recoilm"   ,"(SEL) "+nickname+" recoil mass", 100, 70. , 170. ) ;
  TH1F* SEL_missm       = new TH1F(nickname+"_SEL_missm"    ,"(SEL) "+nickname+" missing mass", 300, -150. , 150. ) ;
  TH1F* SEL_missm_jjm   = new TH1F(nickname+"_SEL_missm_jjm","(SEL) "+nickname+" miss-jj mass difference", 100, -150. , 150. ) ;
  TH1F* SEL_costhmiss  = new TH1F(nickname+"_SEL_costhmiss" ,"(SEL) "+nickname+" cosTheta pMiss", 100, -1. , 1. ) ;
  TH1F* SEL_emiss      = new TH1F(nickname+"_SEL_emiss"     ,"(SEL) "+nickname+" missing energy", 100, 0. , 150. ) ;
  TH1F* SEL_nBtag       = new TH1F(nickname+"_SEL_nBtag"      ,"(SEL) "+nickname+" nBtag", 3, 0. , 3. ) ;

  histos.push_back(ALL_sel);
  histos.push_back(ALL_nmuons);

  histos.push_back(PRE_btag);
  histos.push_back(PRE_njetChg);

  histos.push_back(MC_muonmass);
  histos.push_back(MC_recoilm);
  histos.push_back(MC_missm);
  histos.push_back(MC_costhmiss);

  histos.push_back(PRE_muonmass);
  histos.push_back(PRE_recoilm);
  histos.push_back(PRE_missm);
  histos.push_back(PRE_costhmiss);
  histos.push_back(PRE_jetmass);
  histos.push_back(PRE_missm_jjm);
  histos.push_back(PRE_emiss);
  histos.push_back(PRE_nBtag);

  histos.push_back(SEL_muonmass);
  histos.push_back(SEL_recoilm);
  histos.push_back(SEL_missm);
  histos.push_back(SEL_costhmiss);
  histos.push_back(SEL_jetmass);
  histos.push_back(SEL_missm_jjm);
  histos.push_back(SEL_emiss);
  histos.push_back(SEL_nBtag);

  for (size_t i=0; i<histos.size(); i++) {
    histos[i]->Sumw2(); // we will be reweighting these histograms
    histos[i]->GetXaxis()->SetLabelSize(0.1);
    histos[i]->GetYaxis()->SetLabelSize(0.1);
  }

  //---------------------------------

  // initial CM energy (nb this ignores the ILC's crossing angle, as do the DELPHES-miniDST files)
  const TLorentzVector ecms(0.,0.,0.,250.) ;

  //  fullsim files have the 14 mrad crossing angle included
  //  const TLorentzVector ecms(2.*125.*sin(0.007), 0., 0.,250.) ; // this is the nominal initial CM with crossing angle


  EVENT::LCEvent* evt ;
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;

  // specify which collections will be read
  lcReader->setReadCollectionNames( { pfoColName, jetColName, muoColName, mcpColName, "RecoMCTruthLink" } );

  // loop over the input files
  for (size_t ifile=0; ifile<fnames.size(); ifile++ ) {

    if ( maxEvt>0 && nEvents>=maxEvt ) break;

    bool firstEvent(true);

    lcReader->open( fnames[ifile].Data() ) ;

    //==================== the event loop ============================================================
    while( ( evt = lcReader->readNextEvent()) != 0  && (maxEvt<=0 || nEvents < maxEvt) ) {

      if ( evt->getEventNumber()==-99 ) continue; // this is a "special" event ontaining summary information

      nEvents++;

      // get xsec and beam polarisations from event parameters
      float xsection = evt->parameters().getFloatVal( "CrossSection_fb");
      float epol     = evt->parameters().getFloatVal( "beamPol1");
      float ppol     = evt->parameters().getFloatVal( "beamPol2");

      if ( xsection<=0. ) {
        if (firstEvent) {
          cout << "WARNING, zero xsec; setting to 1.0" << endl;
        }
        xsection=1.;
      }
      if (firstEvent) {
        std::cout << fnames[ifile] << " xsec " << xsection << " electron pol = " << epol << " positron pol = " << ppol << std::endl  ;
        firstEvent=false;
      }

      int iselStep(0); // to keep track of at which selection stage this event was cut
      ALL_sel->Fill(iselStep++, xsection); // this histogram helps us keep track of the number of selected events

      //-----------------------------------------------------------
      // first look at some MC-level information
      //-----------------------------------------------------------
      LCIterator<MCParticle> mcps( evt, mcpColName ) ;

      // find the initial mu+- and b/b-bar in the event listing
      MCParticle* MCmuplus(0);
      MCParticle* MCmuminus(0);
      MCParticle* MCb(0);
      MCParticle* MCbbar(0);
      auto mcp = mcps.next();
      while ( mcp ) {
        switch ( mcp->getPDG() ) {
        case 13:
          if ( !MCmuminus ) MCmuminus=mcp;
          break;
        case -13:
          if ( !MCmuplus ) MCmuplus=mcp;
          break;
        case 5:
          if ( !MCb ) MCb=mcp;
          break;
        case -5:
          if ( !MCbbar ) MCbbar=mcp;
          break;
        default:
          break;
        }
        mcp = mcps.next();
      }

      if ( MCmuminus && MCmuplus && MCb && MCbbar ) { // we have the 4 fermions we're interested in
        const auto& vMCmuminus = v4( MCmuminus );
        const auto& vMCmuplus  = v4( MCmuplus );
        const auto& vMCb       = v4( MCb );
        const auto& vMCbbar    = v4( MCbbar );
        const auto& vMCmiss    = ecms - (vMCmuminus + vMCmuplus + vMCb + vMCbbar);

        // fill some MC-level histograms
        MC_muonmass ->Fill( (vMCmuminus+vMCmuplus).M() );
        MC_recoilm  ->Fill( (ecms-(vMCmuminus+vMCmuplus)).M() );
        MC_missm    ->Fill( vMCmiss.M() );
        MC_costhmiss->Fill( vMCmiss.CosTheta() );
      }


      //-----------------------------------------------------------
      // then reconstructed-level information
      //-----------------------------------------------------------
      LCIterator<ReconstructedParticle> muons( evt, muoColName ) ;
      LCIterator<ReconstructedParticle> jets( evt, jetColName ) ;

      // particle ID handler for the jets (eg to get b-tag information)
      PIDHandler *pidh = new PIDHandler( evt->getCollection( jetColName ) );
      int ilcfi = pidh->getAlgorithmID( pidName );
      int ibtag = pidh->getParameterIndex(ilcfi, "BTag");

      // links from reco -> MC particles
      LCCollection* linkcol = evt->getCollection( "RecoMCTruthLink" );
      LCRelationNavigator reco_mc_Navi( linkcol );

      // fill some histograms
      ALL_nmuons->Fill(muons.size(), xsection);

      // apply a simple preselection
      if( jets.size() == 2 && muons.size() == 2 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }

      // get the reconstructed muon and jet PFOs
      auto mu1 = muons.next();
      auto mu2 = muons.next();

      auto j1 = jets.next();
      auto j2 = jets.next();


      /*
      // example of using the relations manager
      cout << "number of MC particles linked to the isolated muon PFO = " << reco_mc_Navi.getRelatedToObjects(mu1).size() << endl;
      for (size_t k=0; k<reco_mc_Navi.getRelatedToObjects(mu1).size(); k++) {
	cout << k << " " << reco_mc_Navi.getRelatedToObjects(mu1)[k] << endl;
	if ( reco_mc_Navi.getRelatedToObjects(mu1)[k] ) {
	  MCParticle* mcp = dynamic_cast <MCParticle*> (reco_mc_Navi.getRelatedToObjects(mu1)[k]);
	  float wgt = reco_mc_Navi.getRelatedToWeights(mu1)[k]; // weight of the connection; this is a combination of "track" and "cluster" weights
	  float trackwgt = (int(wgt)%10000)/1000.;  // decomposed into track weight (what fraction of tracker hits in this PFO were created by this particle)
	  float clusterwgt = (int(wgt)/10000)/1000. ; // similar, for calorimeter hit energies
	}
      }
      */

      // count # charged particles in each jet
      int j1_nch(0);
      int j2_nch(0);
      for ( size_t i=0; i<j1->getParticles().size(); i++) {
        if ( j1->getParticles()[i]->getCharge()!=0 ) j1_nch++;
      }
      for ( size_t i=0; i<j2->getParticles().size(); i++) {
        if ( j2->getParticles()[i]->getCharge()!=0 ) j2_nch++;
      }
      PRE_njetChg->Fill( j1_nch, xsection );
      PRE_njetChg->Fill( j2_nch, xsection );

      // btag information for the two jets
      float flvtag1 = pidh->getParticleID(j1, ilcfi).getParameters()[ibtag];
      float flvtag2 = pidh->getParticleID(j2, ilcfi).getParameters()[ibtag];

      int btag1(0);
      int btag2(0);

      if ( minDstFlavour==0 ) { // delphes
        // this is how to get other b/c-tag info from the jets (DELPHES)
        // int btag1_90 = (flvtag1 & ( 1 << 0 )) >> 0; // btag hiEff  90% eff
        // int btag1_70 = (flvtag1 & ( 1 << 1 )) >> 1; // btag medEff 70% eff
        // int btag1_50 = (flvtag1 & ( 1 << 2 )) >> 2; // btag hiPur  50% eff
        // int ctag1_55 = (flvtag1 & ( 1 << 3 )) >> 3; // ctag hiEff  55% eff
        // int ctag1_30 = (flvtag1 & ( 1 << 4 )) >> 4; // ctag medEff 30% eff
        // int ctag1_20 = (flvtag1 & ( 1 << 5 )) >> 5; // ctag hiPur  20% eff
        btag1 = ( int(flvtag1) & ( 1 << 0 )) >> 0;
        btag2 = ( int(flvtag2) & ( 1 << 0 )) >> 0;

	PRE_btag->Fill(btag1, xsection);
	PRE_btag->Fill(btag2, xsection);

      } else if ( minDstFlavour==1 ) { // SGV, cut on lcfi btag output
	PRE_btag->Fill(flvtag1, xsection);
	PRE_btag->Fill(flvtag2, xsection);

        btag1 = int(flvtag1>0.7); // user-defined cut on b-tag score
        btag2 = int(flvtag2>0.7);
      }


      // TLorentzVectors for the objects
      const auto& vm1 = v4(mu1) ;
      const auto& vm2 = v4(mu2) ;
      const auto& vj1 = v4Jet(j1) ;
      const auto& vj2 = v4Jet(j2) ;

      TLorentzVector recoil  = ecms - ( vm1 + vm2 ) ; // recoil against mu mu
      TLorentzVector visible = vm1 + vm2 + vj1 + vj2; // visible 4-mom
      TLorentzVector missing = ecms - visible ; // missing 4-mom

      // calculate kinematic variables
      float mass_mumu = inv_mass( mu1, mu2);
      float mass_jj = (vj1+vj2).M();
      float mass_recoil = recoil.M();
      float mass_miss = missing.M();

      PRE_muonmass->Fill( mass_mumu , xsection) ;
      PRE_jetmass->Fill( mass_jj , xsection ) ;
      PRE_recoilm->Fill( mass_recoil , xsection ) ;
      PRE_missm->Fill( mass_miss , xsection ) ;
      PRE_missm_jjm->Fill( mass_miss-mass_jj , xsection );
      PRE_costhmiss->Fill( missing.CosTheta() , xsection );
      PRE_emiss->Fill( missing.E() , xsection );
      PRE_nBtag->Fill( btag1+btag2 , xsection );

      //------------------------
      // apply full event selection
      //------------------------

      // >3 chg particles per jet to reduce tau BG
      if ( j1_nch>3 && j2_nch>3 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }

      // mu-mu invariant mass, charge
      if ( mu1->getCharge()!=mu2->getCharge() && mass_mumu>81 && mass_mumu<101 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }

      // jet-jet invariant mass
      if ( mass_jj<60 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }

      // mass recoiling against mu-mu system
      if ( mass_recoil > 120 && mass_recoil < 140 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }

      // both jets have b-tag
      if (btag1==1 && btag2==1 ) {
        ALL_sel->Fill(iselStep++, xsection);
      } else {
        continue;
      }


      // some other possible cuts:
      // // missing energy
      // if ( missing.E() > 30 ) {
      //   ALL_sel->Fill(iselStep++, xsection);
      // } else {
      //   continue;
      // }
      // 
      // // polar anagle of missing momentum
      // if ( fabs( missing.CosTheta() ) < 0.95 ) {
      //   ALL_sel->Fill(iselStep++, xsection);
      // } else {
      //   continue;
      // }
      // // difference between missing mass and jet-jet mass
      // if ( mass_miss-mass_jj > -50 && mass_miss-mass_jj < 70 ) {
      //   ALL_sel->Fill(iselStep++, xsection);
      // } else {
      //   continue;
      // }
      // 

      // ----------------- end selection --------------------

      // fill histograms after selection
      SEL_muonmass->Fill( mass_mumu , xsection ) ;
      SEL_jetmass->Fill( mass_jj , xsection );
      SEL_recoilm->Fill( mass_recoil  , xsection) ;
      SEL_missm->Fill( mass_miss  , xsection) ;
      SEL_missm_jjm->Fill( mass_miss-mass_jj , xsection );
      SEL_costhmiss->Fill( missing.CosTheta() , xsection );
      SEL_emiss->Fill( missing.E() , xsection );
      SEL_nBtag->Fill( btag1+btag2 , xsection );

    } // event loop

    lcReader->close() ;

  }

  delete lcReader ;

  // scale according to # MC events used from input sample
  for (size_t i=0; i<histos.size(); i++) {
    histos[i]->Scale(1./nEvents);
  }

  std::cout << nEvents
            << " events read from process : "
            << nickname << std::endl  ;

  return histos;
}




void runall(int _maxevt=0) {

  maxEvt=_maxevt;

  // the different processes to consider
  std::map < TString, std::vector< TString > > allProcs;

  // combine different process in the histograms (eg different polarisations for same final state)

  // the presence of "LR", "RL" etc in the nickname is important: it's used to decide the correct weight
  //   [actually this is a ugly and unnecessary, since the beam polarisations are kept in the event data "beamPol1" and "beamPol2"]

  std::map < TString, std::vector< TString > > hGroups;
  std::map < TString, int > gcols;

  // -- the signal (DELPHES)
  allProcs["sig40_LR"]=std::vector< TString > ();
  allProcs["sig40_LR"].push_back( inputDir_SIGNAL+"eeZH_m40_LR.delphes.slcio" );
  allProcs["sig40_RL"]=std::vector< TString > ();
  allProcs["sig40_RL"].push_back( inputDir_SIGNAL+"eeZH_m40_RL.delphes.slcio" ) ;
  hGroups["Signal40"]=std::vector< TString > (); // this is a group of processes for histogramming
  hGroups["Signal40"].push_back("sig40_LR");     // ...to which we'll add the LR
  hGroups["Signal40"].push_back("sig40_RL");     // ...and RL versions of this final state
  gcols["Signal40"]=4;                           // the colour for this group's histograms

  /*
  allProcs["sig30_LR"]=std::vector< TString > ();
  allProcs["sig30_LR"].push_back( inputDir_SIGNAL+"eeZH_m30_LR.delphes.slcio" );
  allProcs["sig30_RL"]=std::vector< TString > ();
  allProcs["sig30_RL"].push_back( inputDir_SIGNAL+"eeZH_m30_RL.delphes.slcio" ) ;
  hGroups["Signal30"]=std::vector< TString > ();
  hGroups["Signal30"].push_back("sig30_LR");
  hGroups["Signal30"].push_back("sig30_RL");
  gcols["Signal30"]=2;

  allProcs["sig20_LR"]=std::vector< TString > ();
  allProcs["sig20_LR"].push_back( inputDir_SIGNAL+"eeZH_m20_LR.delphes.slcio" );
  allProcs["sig20_RL"]=std::vector< TString > ();
  allProcs["sig20_RL"].push_back( inputDir_SIGNAL+"eeZH_m20_RL.delphes.slcio" ) ;
  hGroups["Signal20"]=std::vector< TString > ();
  hGroups["Signal20"].push_back("sig20_LR");
  hGroups["Signal20"].push_back("sig20_RL");
  gcols["Signal20"]=6;
  */

  /*
  // background (DELPHES)
  // -- 4f_zz_sl
  allProcs["BGzzsl_LR"]=std::vector< TString > ();
  allProcs["BGzzsl_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eL.pR.I106575.001.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eL.pR.I106575.002.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eL.pR.I106575.003.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eL.pR.I106575.004.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eL.pR.I106575.005.stdhep.delphes_card_ILCgen.tcl.slcio");

  allProcs["BGzzsl_RL"]=std::vector< TString > ();
  allProcs["BGzzsl_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eR.pL.I106576.001.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eR.pL.I106576.002.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eR.pL.I106576.003.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eR.pL.I106576.004.stdhep.delphes_card_ILCgen.tcl.slcio");
  allProcs["BGzzsl_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.P4f_zz_sl.Gwhizard-1_95.eR.pL.I106576.005.stdhep.delphes_card_ILCgen.tcl.slcio");

  hGroups["ZZ"]=std::vector< TString > ();
  hGroups["ZZ"].push_back("BGzzsl_LR");
  hGroups["ZZ"].push_back("BGzzsl_RL");
  gcols["ZZ"]=1;

  // -- mu mu H
  allProcs["BGmmh_LR"]=std::vector< TString > ();
  allProcs["BGmmh_LR"].push_back(inputDir_BACKGD+"E250-TDR_ws.Pe2e2h.Gwhizard-1_95.eL.pR.I106479.001.stdhep.delphes_card_ILCgen.tcl.slcio" );
  allProcs["BGmmh_RL"]=std::vector< TString > ();
  allProcs["BGmmh_RL"].push_back(inputDir_BACKGD+"E250-TDR_ws.Pe2e2h.Gwhizard-1_95.eR.pL.I106480.001.stdhep.delphes_card_ILCgen.tcl.slcio" );
  hGroups["mmh"]=std::vector< TString > ();
  hGroups["mmh"].push_back("BGmmh_LR");
  hGroups["mmh"].push_back("BGmmh_RL");
  gcols["mmh"]=8;
  */

  /*
  // -- the signal (SGV) [n.b. you should probably not mix DELPHES and/or SGV and/or FULLSIM samples if you can avoid it]
  allProcs["sig40_LR"]=std::vector< TString > ();
  allProcs["sig40_LR"].push_back( inputDir_SIGNAL+"SGV-mini-DST-eeZH_m40_LR_SGVDST.slcio" );
  allProcs["sig40_RL"]=std::vector< TString > ();
  allProcs["sig40_RL"].push_back( inputDir_SIGNAL+"SGV-mini-DST-eeZH_m40_RL_SGVDST.slcio" );
  hGroups["Signal40"]=std::vector< TString > ();
  hGroups["Signal40"].push_back("sig40_LR");
  hGroups["Signal40"].push_back("sig40_RL");
  gcols["Signal40"]=4;

  allProcs["sig30_LR"]=std::vector< TString > ();
  allProcs["sig30_LR"].push_back( inputDir_SIGNAL+"SGV-mini-DST-eeZH_m30_LR_SGVDST.slcio" );
  allProcs["sig30_RL"]=std::vector< TString > ();
  allProcs["sig30_RL"].push_back( inputDir_SIGNAL+"SGV-mini-DST-eeZH_m30_RL_SGVDST.slcio" );
  hGroups["Signal30"]=std::vector< TString > ();
  hGroups["Signal30"].push_back("sig30_LR");
  hGroups["Signal30"].push_back("sig30_RL");
  gcols["Signal30"]=4;
  */


  _fout = new TFile(outPlotName+".root", "recreate");

  std::map < TString, std::vector<TH1F*> > allHistos;

  // make and fill histograms for each channel
  int icol=1;
  for ( std::map < TString, std::vector< TString > >::iterator jj=allProcs.begin(); jj!=allProcs.end(); jj++) {
    allHistos[jj->first] = analyse_process( jj->first, jj->second , minDstFlavour );
    for ( size_t k=0; k<allHistos[jj->first].size(); k++) {
      allHistos[jj->first][k]->SetLineColor( (icol%9)+1 );
    }
    icol++;
  }

  TCanvas* cc = new TCanvas();
  cc->Print(outPlotName+".pdf[");

  // first draw the raw histograms for each process
  // these are weighted for xsec (for the "pure" sample input polarisation) and for the #MC events used
  //  for (std::map < TString, std::vector<TH1F*> >::iterator jtt = allHistos.begin(); jtt!=allHistos.end(); jtt++) { // the processes
  //    cc->Clear();
  //    cc->Divide(4,6);
  //    int ic=1;
  //    for ( size_t k=0; k<jtt->second.size(); k++) { // the different variable histos
  //      cc->cd(ic++);
  //      TH1F* h = jtt->second[k];
  //      h->GetXaxis()->SetRange(0, h->GetNbinsX() + 1);
  //      h->Draw("hist");
  //    }
  //    cc->Print(outPlotName+".pdf");
  //  }

  _fout->cd();
  gStyle->SetTitleSize(0.08,"");
  gStyle->SetTitleOffset(0.02,"");
  gStyle->SetOptStat(0);

  TText tt;
  tt.SetTextSize(0.03);
  for (int iPolSet = 0; iPolSet<4; iPolSet++) { // the 4 pol sets -+, +-, --, ++

    float epol_target, ppol_target;
    float lumi;
    TString polLab="pol_undef";
    if ( iPolSet==0 ) { // eLpR
      epol_target = -0.8;
      ppol_target = +0.3;
      lumi=900.;
      polLab="pol_eL80_pR30";
    } else if ( iPolSet==1 ) { // eRpL
      epol_target = +0.8;
      ppol_target = -0.3;
      lumi=900.;
      polLab="pol_eR80_pL30";
    } else if ( iPolSet==2 ) { // eLpL
      epol_target = -0.8;
      ppol_target = -0.3;
      lumi=100.;
      polLab="pol_eL80_pL30";
    } else if ( iPolSet==3 ) { // eRpR
      epol_target = +0.8;
      ppol_target = +0.3;
      lumi=100.;
      polLab="pol_eR80_pR30";
    } else if ( iPolSet==4 ) { // nominal unpolarised, 1 ab-1
      epol_target = 0.;
      ppol_target = 0.;
      lumi=1000.;
      polLab="pol_e0_p0";
    } else {
      cout << "unknown beam polarisation set (0=eLpR, 1=eRpL, 2=eLpL, 3=eRpR, 4=unpol)" << endl;
      return;
    }

    _fout->mkdir(polLab)->cd();

    polLab = polLab + "  int lumi = ";
    polLab+=lumi;
    polLab+=" fb^-1";

    // these are the weights to be applied to samples with 100% beam polarisation used in their generation
    float wtLR = (1-epol_target)*(1+ppol_target)/4.;
    float wtRL = (1+epol_target)*(1-ppol_target)/4.;
    float wtLL = (1-epol_target)*(1-ppol_target)/4.;
    float wtRR = (1+epol_target)*(1+ppol_target)/4.;

    TLegend* tl = new TLegend(0.6, 0.3, 0.9, 0.7);

    // now group the histograms, and weight according to polarisation & lumi
    std::map < TString, std::vector<TH1F*> > groupedHistos;
    for ( std::map < TString, std::vector< TString > >::iterator jj = hGroups.begin(); jj!=hGroups.end(); jj++ ) { // loop over groups
      groupedHistos[jj->first]=std::vector<TH1F*> ();
      std::vector< TString > subChanNames = jj->second;

      for ( size_t i=0; i<subChanNames.size(); i++) { // the subchannels within the group

        // choose polarisation weight according to process nickname
        float polweight = 0;
        if ( subChanNames[i].Contains( "LR" ) ) {
          polweight=wtLR;
        } else if ( subChanNames[i].Contains( "RL" ) ) {
          polweight=wtRL;
        } else if ( subChanNames[i].Contains( "RR" ) ) {
          polweight=wtRR;
        } else if ( subChanNames[i].Contains( "LL" ) ) {
          polweight=wtLL;
        } else {
          cout << "ERROR, could not determine sample beam polarisation from process nickname! setting weight to 0." << endl;
        }
        polweight*=lumi; //scale by integrated luminosity


	// add the grouped processes with correct weights
        std::vector<TH1F*> subHists = allHistos[subChanNames[i]];
        for ( size_t k=0; k<subHists.size(); k++) { // the list of histograms
          if (i==0) { // first one, clone and rename histogram
            TString newname=subHists[k]->GetName();
            TString newtitle=subHists[k]->GetTitle();
            newname = newname.ReplaceAll( subChanNames[i], jj->first );
            newtitle = newtitle.ReplaceAll( subChanNames[i], jj->first );
            TH1F* hnew = (TH1F*) subHists[k]->Clone(newname);
            hnew->SetTitle( newtitle );
            hnew->Scale(polweight);
            hnew->SetLineColor(gcols[jj->first]);
            hnew->SetMarkerColor(gcols[jj->first]);
            groupedHistos[jj->first].push_back( hnew );
            if ( k==0 ) {
              tl->AddEntry( hnew, jj->first, "l" );
            }
          } else {
            groupedHistos[jj->first][k]->Add( subHists[k], polweight );
          }
        }


      } // subchannels

    } // groups


    // draw the grouped histograms on top of each other
    cc->Clear();
    cc->Divide(4,6);
    int ic=1;
    for ( size_t k=0; k<groupedHistos.begin()->second.size(); k++) { // the different variable histos
      TVirtualPad* tp = cc->cd(ic++);
      bool first=true;
      TH1F* hh;
      float hmax(0);
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
        if ( jj->second[k]->GetMaximum()>hmax) hmax=jj->second[k]->GetMaximum();
      }
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
        jj->second[k]->SetMaximum(hmax*1.3);
        jj->second[k]->SetMinimum(0.1);
        if (first) {
          first=false;
          jj->second[k]->Draw("hist");
        } else {
          jj->second[k]->Draw("same;hist");
        }
      }
      if ( k==0 ) tl->Draw();
      tp->SetLogy(1);
    }
    cc->cd();
    tt.SetTextColor(2);
    tt.DrawTextNDC(0.01, 0.83, polLab);
    cc->Print(outPlotName+".pdf");


    // count selected events
    cout << endl;
    cout << " --------------------------------------------------" << endl;
    cout << "EVENT cut table " << polLab << endl;
    cout << " --------------------------------------------------" << endl;

    for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
      cout << std::setw(30) << jj->first ;
    }
    cout << endl;

    cout << " ---EVENTS----------------" << endl;
    for (int i=1; i<=10; i++) { // selection steps
      cout << "CUT" << std::setw(2) << i;
      float totBG(0);
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
        TH1F* hsel=jj->second[0];
        // GetBinContent gives the sum of the weights
        // GetBinError gives the sqrt of the sum of the squared weights
        cout << std::setw(15) << hsel->GetBinContent(i) << " '+/-' " << std::setw(15) << hsel->GetBinError(i) ;
        if ( ! jj->first.Contains("Signal") ) {
          totBG+=hsel->GetBinContent(i);
        }
      }
      cout << "  TOTBG " << totBG << endl;
    }

    // calculate the efficiency for each process at each step of the selection
    cout << endl << " ---EFFICIENCY----------------" << endl;
    for (int i=1; i<=10; i++) { // selection steps
      cout << "CUT" << std::setw(2) << i;
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
        TH1F* hsel=jj->second[0];
        // GetBinContent gives the sum of the weights
        // GetBinError gives the sqrt of the sum of the squared weights
        float efficiency = hsel->GetBinContent(i)/hsel->GetBinContent(1);
        float effectiveN = pow(hsel->GetBinContent(1),2)/pow(hsel->GetBinError(1), 2); // (square of sum)/(sum of squares)
        float effErr = sqrt( efficiency*(1.-efficiency)/effectiveN );
        cout << std::setw(15) << efficiency << "  +/-  " << std::setw(15) << effErr ;
      }
      cout << endl;
    }

    cout << endl << " ---SELECTION SUMMARY----------------" << endl;
    for (int icut=1; icut<=10; icut++) { // selection steps
      cout << "CUT" << std::setw(2) << icut;
      float totBG(0);
      // sum the backgrounds
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
	TH1F* hsel=jj->second[0];
	if ( ! jj->first.Contains("Signal") ) {
	  totBG += hsel->GetBinContent(icut);
	}
      }
      cout << " ;  expected SM-BG events " << std::setw(8) << totBG ;

      cout <<  " expected signal events (@ 100% BR) : ";
      // the different signals
      for ( std::map < TString, std::vector<TH1F*> >::iterator jj = groupedHistos.begin(); jj!=groupedHistos.end(); jj++ ) { // loop over groups
	TH1F* hsel=jj->second[0];
	if ( jj->first.Contains("Signal") ) {
	  float nsig = hsel->GetBinContent(icut);
	  cout << jj->first << " " << std::setw(8) << nsig << " ; " ;
	}
      }

      cout << endl;

    }

  }

  cc->Print(outPlotName+".pdf]");


  _fout->Write();
  _fout->Close();

  return;
}
