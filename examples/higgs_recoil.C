#ifndef __CINT__ 
#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCEvent.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCIterator.h"
#endif

#include "TH1F.h"
#include "TLorentzVector.h"

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

using namespace lcio;

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

// Example script for creating invariant mass plots from a ee -> HZ-> X mu mu sample

void higgs_recoil_plots(const char* FILEN, TString outname = "recoil_plots_boost") {

  int nEvents = 0;
  int maxEvt = 10000;  // change as needed

  //---- create some histograms ----

  TH1F* hmuonmass   = new TH1F("hmuonmass",";inv. di-muon mass [GeV]; ; ", 100,  60. , 120. );
  TH1F* hrecoilm    = new TH1F("hrecoilm",";recoil mass [GeV]; ; ", 100, 110. , 170. );
  TH1F* hjetmass    = new TH1F("hjetmass",";inv. di-jet mass [GeV]; ; ", 100, 0. , 150. );
  TH1F* hjetmassno  = new TH1F("hjetmassno",";inv. di-jet mass w/o overlay [GeV]; ; ", 100, 0. , 150. );
  TH1F* hevis       = new TH1F("hevis",";E_{vis} [GeV]; ; ", 300, 0., 300.);
  TH1F* hevisno     = new TH1F("hevisno",";E_{vis} [GeV]; ; ", 300, 0., 300.);

  TH1F* hctags      = new TH1F("hctags",";c-tag values of jets; ; ", 10, 0., 1. );
  TH1F* hbtags      = new TH1F("hbtags",";b-tag values of jets; ; ", 100, 0. , 1. ) ;
  TH1F* hcthbthi    = new TH1F("hcthbthi",";cos #theta of b-tagged jets (high); ; ", 50, 0., 1. );
  TH1F* hcthcthi    = new TH1F("hcthcthi",";cos #theta of c-tagged jets  (high); ; ", 10, 0. , 1. );
  TH1F* hcthbtme    = new TH1F("hcthbtme",";cos #theta of b-tagged jets (medium); ; ", 50, 0., 1. );
  TH1F* hcthctme    = new TH1F("hcthctme",";cos #theta of c-tagged jets  (medium); ; ", 10, 0. , 1. );
  TH1F* hcthbtlo    = new TH1F("hcthbtlo",";cos #theta of b-tagged jets (low); ; ", 50, 0., 1. );
  TH1F* hcthctlo    = new TH1F("hcthctlo",";cos #theta of c-tagged jets  (low); ; ", 10, 0. , 1. );

  //------ collection to use -------
  // DBD & SGV
  const char* pfoColName = "PandoraPFOs";
  const char* jetColName = "Refined2Jets";
  const char* muoColName = "IsolatedMuons";
  const char* mcpColName = "MCParticle";
  //const char* mcpColName = "MCParticlesSkimmed";
  const char* btagAlgoName = "lcfiplus";
  const char* btagParName = "BTag";
  const char* ctagParName = "CTag";
  bool isDelphes = false;
 
  //  // Delphes
  //  const char* pfoColName = "PFOs";
  //  //const char* jetColName = "Jets";
  //  const char* jetColName = "Durham2Jets";
  //  const char* muoColName = "Muons";
  //  const char* mcpColName = "MCParticle";
  //  const char* btagAlgoName = "JetParameters";
  //  const char* btagParName = "BTag";
  //  const char* ctagParName = "BTag";
  //  bool isDelphes = true;

  //---------------------------------
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  lcReader->setReadCollectionNames( { pfoColName, jetColName, muoColName, mcpColName} );
  lcReader->open( FILEN );

  EVENT::LCEvent* evt;
 
  int passHjj = 0;
  int fail2muon = 0;
  int fail2jets = 0;
  int noverlay = 0;

  //==================== the event loop ============================================================
  while( ( evt = lcReader->readNextEvent()) != 0  && nEvents++ < maxEvt ) {

    LCIterator<MCParticle> mcps( evt, mcpColName );
    auto mcp1 = mcps.next();
    if (!mcp1) continue;

    // check Higgs decay mode, only take H -> qq / gg 
    while (abs(mcp1->getPDG())!=25 || abs(mcp1->getDaughters()[0]->getPDG())==25) mcp1 = mcps.next();
    if ( abs(mcp1->getDaughters()[0]->getPDG())>5 && abs(mcp1->getDaughters()[0]->getPDG()) != 21 )  continue;
    passHjj++;
   
    // check if overlay present in event
    bool  nooverlay = true;
    LCIterator<MCParticle> mcpo( evt, mcpColName );
    auto mcp2 = mcpo.next();
    while (mcp2 && !mcp2->isOverlay()) mcp2 = mcps.next();  
    if (mcp2 != 0) nooverlay = false;  
    if (!nooverlay) noverlay++;
   
    double evis = 0.;
    LCIterator<ReconstructedParticle> pfos( evt, pfoColName );
    auto pfo1 = pfos.next();
    while (pfo1) {
      evis += pfo1->getEnergy();
      pfo1 = pfos.next();  
    }
    hevis->Fill( evis );
    if (nooverlay) hevisno->Fill( evis );
   
    LCIterator<ReconstructedParticle> jets( evt, jetColName );
    LCIterator<ReconstructedParticle> muons( evt, muoColName );

    if( muons.size() != 2 ) {fail2muon++; continue;}
    if( jets.size() != 2 ) {fail2jets++; continue;}

    // the invariant di-muon mass
    auto mu1 = muons.next(); 
    auto mu2 = muons.next(); 
    hmuonmass->Fill( inv_mass( mu1, mu2 ) );

    // the recoil mass
    const auto& vm1 = v4(mu1);
    const auto& vm2 = v4(mu2);
    double pxinitial = 0.;
    double Einitial = 250.;
    if (!isDelphes) {  // correct for xing angle boost
      pxinitial = 250.*0.007; 
      Einitial = 250.049;
    }
    TLorentzVector ecms(pxinitial,0.,0.,Einitial);
    TLorentzVector recoil = ecms - ( vm1 + vm2 );
    hrecoilm->Fill( recoil.M() );

    // the invariant di-jet mass
    auto j1 = jets.next(); 
    auto j2 = jets.next(); 
    hjetmass->Fill( inv_mass( j1, j2 ) );
    if (nooverlay) hjetmassno->Fill( inv_mass( j1, j2 ) );
   
    // the b-tag values (only for H->bb events)
    LCCollection* r2j = evt->getCollection( jetColName );
    PIDHandler *pidh = new PIDHandler(r2j);
    int ilcfi = pidh->getAlgorithmID(btagAlgoName);
    int ibtag = pidh->getParameterIndex(ilcfi, btagParName);
    int ictag = pidh->getParameterIndex(ilcfi, ctagParName);
    float btagj1 = pidh->getParticleID(j1, ilcfi).getParameters()[ibtag];
    float btagj2 = pidh->getParticleID(j2, ilcfi).getParameters()[ibtag];
    float ctagj1 = pidh->getParticleID(j1, ilcfi).getParameters()[ictag];
    float ctagj2 = pidh->getParticleID(j2, ilcfi).getParameters()[ictag];

    const auto& vj1 = v4(j1);
    const auto& vj2 = v4(j2);
    double costhj1 = vj1.CosTheta();
    double costhj2 = vj2.CosTheta();
    if (isDelphes) {
      std::bitset<8> ftagsj1 = int(btagj1); 
      std::bitset<8> ftagsj2 = int(btagj2); 
      if ( abs(mcp1->getDaughters()[0]->getPDG()) == 4 ) {
	if (ftagsj1[6]==1) hcthcthi->Fill ( costhj1 );
	if (ftagsj2[6]==1) hcthcthi->Fill ( costhj2 );
	if (ftagsj1[5]==1) hcthctme->Fill ( costhj1 );
	if (ftagsj2[5]==1) hcthctme->Fill ( costhj2 );
	if (ftagsj1[4]==1) hcthctlo->Fill ( costhj1 );
	if (ftagsj2[4]==1) hcthctlo->Fill ( costhj2 );
      }  
      if ( abs(mcp1->getDaughters()[0]->getPDG()) == 5 ) {
	if (ftagsj1[2]==1) hcthbthi->Fill ( costhj1 );
	if (ftagsj2[2]==1) hcthbthi->Fill ( costhj2 );
	if (ftagsj1[1]==1) hcthbtme->Fill ( costhj1 );
	if (ftagsj2[1]==1) hcthbtme->Fill ( costhj2 );
	if (ftagsj1[0]==1) hcthbtlo->Fill ( costhj1 );
	if (ftagsj2[0]==1) hcthbtlo->Fill ( costhj2 );
      }
    }
    else {
      // b/c likeliness
      if ( abs(mcp1->getDaughters()[0]->getPDG()) == 5 ) {
	hbtags->Fill( btagj1 );
	hbtags->Fill( btagj2 );
	if (btagj1 > 0.05) hcthbtlo->Fill ( costhj1 );
	if (btagj2 > 0.05) hcthbtlo->Fill ( costhj2 );
	if (btagj1 > 0.10) hcthbtme->Fill ( costhj1 );
	if (btagj2 > 0.10) hcthbtme->Fill ( costhj2 );
	if (btagj1 > 0.50) hcthbthi->Fill ( costhj1 );
	if (btagj2 > 0.50) hcthbthi->Fill ( costhj2 );
      } 
      if ( abs(mcp1->getDaughters()[0]->getPDG()) == 4 ) {
	hctags->Fill( ctagj1 );
	hctags->Fill( ctagj2 );
	if (ctagj1 > 0.20) hcthctlo->Fill ( costhj1 );
	if (ctagj2 > 0.20) hcthctlo->Fill ( costhj2 );
	if (ctagj1 > 0.40) hcthctme->Fill ( costhj1 );
	if (ctagj2 > 0.40) hcthctme->Fill ( costhj2 );
	if (ctagj1 > 0.70) hcthcthi->Fill ( costhj1 );
	if (ctagj2 > 0.70) hcthcthi->Fill ( costhj2 );
      }  
    }

    delete pidh;
  }

  //===================================================================================================
  gStyle->SetOptStat(10);
  //===================================================================================================
  
  //  TCanvas* c1 = new TCanvas("recoil plots","recoil plots",600,800);
  // 
  //  c1->Divide(2,3,0.001);
  //  c1->cd(1) ;
  //  hmuonmass->Draw() ;
  //  c1->cd(2) ;
  //  hrecoilm->Draw() ;
  //  c1->cd(3) ;
  //  hjetmass->Draw() ;
  //  c1->cd(4) ;
  //  hjetmassno->Draw() ;
  //  c1->cd(5) ;
  //  hevis->Draw();
  //  c1->cd(6) ;
  //  hevisno->Draw();
  //  
  //  TCanvas* c2 = new TCanvas("btag plots","btag plots",600,800);
  //  
  //  c2->Divide(2,4,0.001);
  //  c2->cd(1) ;
  //  hbtags->Draw();
  //  c2->cd(2) ;
  //  hctags->Draw();
  //  c2->cd(3) ;
  //  hcthbthi->Draw();
  //  c2->cd(4) ;
  //  hcthcthi->Draw();
  //  c2->cd(5) ;
  //  hcthbtme->Draw();
  //  c2->cd(6) ;
  //  hcthctme->Draw();
  //  c2->cd(7) ;
  //  hcthbtlo->Draw();
  //  c2->cd(8) ;
  //  hcthctlo->Draw();

  TString outnamepdf = outname + ".pdf";
  TString outnameRoot  = outname + ".root";
  
  TFile outputFile (outnameRoot,"RECREATE");

  hmuonmass->Write();  
  hrecoilm->Write();  
  hjetmass->Write();  
  hjetmassno->Write();  
  hevis->Write();  
  hevisno->Write(); 
  
  hbtags->Write();  
  hctags->Write();  
  hcthbthi->Write();  
  hcthcthi->Write();  
  hcthbtme->Write();  
  hcthctme->Write();  
  hcthbtlo->Write();  
  hcthctlo->Write();  
 
  outputFile.Close();
  //===================================================================================================

  std::cout << std::endl 
	    <<  "  " <<  nEvents 
	    << " events read from file: " 
	    << FILEN << std::endl;
  std::cout <<  " out of which " <<  passHjj  
	    << " events are H->bb/cc/gg " << std::endl;
  std::cout <<  " out of which " <<  fail2muon  
	    << " events don't have two muons " << std::endl;
  std::cout <<  " and " <<  fail2jets  
	    << " additional events don't have two jets " << std::endl;
  std::cout << noverlay << " H->bb/cc/gg events have overlay" << std::endl;         
  
  
  lcReader->close();
  delete lcReader;
}
