#ifndef __CINT__ 
#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
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


using namespace lcio ;



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

/** Example script for creating invariant mass plots from a ee -> HZ-> X mu mu sample
 *
 */
 
void higgs_recoil(const char* FILEN, TString outname = "recoil_plot") {

 int nEvents = 0  ;
 int maxEvt = 10000 ;  // change as needed

//---- create histogram ----

 TH1F* hrecoilm    = new TH1F("hrecoilm",";recoil mass [GeV]; ; ", 100, 110. , 170. ) ;

//------ collection to use: 
//------ ALWAYS read also PandoraPFOs / PFOs when reading IsolatedX collections! -------

 const char* muoColName = "IsolatedMuons";
	
 // DBD & SGV	
 const char* pfoColName = "PandoraPFOs";
 bool isDelphes = false;
 
//  // delphes2lcio
//  const char* pfoColName = "PFOs";
//  bool isDelphes = true;
 
//---------------------------------
 IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
 lcReader->setReadCollectionNames( { pfoColName, muoColName } );
 lcReader->open( FILEN ) ;
  

 EVENT::LCEvent* evt ;
 
 int fail2muon = 0;
 
  
 //==================== the event loop ============================================================
 while( ( evt = lcReader->readNextEvent()) != 0  && nEvents++ < maxEvt ) {

   
   LCIterator<ReconstructedParticle> muons( evt, muoColName ) ;

   // for simplicity only consider events with two muons
   if( muons.size() != 2 ) {fail2muon++; continue;}

   // the two muons and their four-vectors
   auto mu1 = muons.next(); 
   auto mu2 = muons.next();    
   const auto& vm1 = v4(mu1) ;
   const auto& vm2 = v4(mu2) ;
   
   // the recoil mass
   double pxinitial = 0.;
   double Einitial = 250.;
   // in full sim & SGV, correct for crossing angle
   if (!isDelphes) {  
     pxinitial = Einitial*0.007; 
     Einitial = 2.*std::sqrt(std::pow(Einitial/2.,2) + std::pow(pxinitial,2));
   }
   TLorentzVector ecms(pxinitial,0.,0.,Einitial) ;
   TLorentzVector recoil = ecms - ( vm1 + vm2 ) ;
   hrecoilm->Fill( recoil.M() ) ;
   
 }

 //===================================================================================================
 gStyle->SetOptStat(10);
 //===================================================================================================
  
 TCanvas* c1 = new TCanvas("recoil plots","recoil plot",800,600);
 hrecoilm->Draw() ;

 TString outnamepdf = outname + ".pdf";
 c1->SaveAs(outnamepdf);

 TString outnameRoot  = outname + ".root";
 TFile outputFile (outnameRoot,"RECREATE");
 
 hrecoilm->Write();  
 
 outputFile.Close();
 //===================================================================================================

  std::cout << std::endl 
	    <<  "  " <<  nEvents 
	    << " events read from file: " 
	    << FILEN << std::endl  ;
  std::cout <<  " out of which " <<  fail2muon  
	    << " events don't have two muons " << std::endl  ;
  
  
  lcReader->close() ;
  delete lcReader ;
}
