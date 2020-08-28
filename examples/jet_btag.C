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



/** Example script for creating a plot of the b-likeliness of jets in a ee -> HZ-> X mu mu sample
 *
 */
 
void jet_btag(const char* FILEN, TString outname = "btag_plot") {

 int nEvents = 0  ;
 int maxEvt = 10000 ;  // change as needed

//---- create histogram ----

 TH1F* hbtag    = new TH1F("hbtag",";b-likeliness; ; ", 100, 0. , 1. ) ;

//------ collection to use: 

 // DBD & SGV
 const char* jetColName = "Refined2Jets";
 const char* btagAlgoName = "lcfiplus";
 const char* btagParName = "BTag";
 bool isDelphes = false;
 
//  // delphes2lcio
//  const char* jetColName = "Jets";
//  const char* btagAlgoName = "JetParameters";
//  const char* btagParName = "BTag";
//  bool isDelphes = true;
 
//---------------------------------
 IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
 lcReader->setReadCollectionNames( { jetColName } );
 lcReader->open( FILEN ) ;
  

 EVENT::LCEvent* evt ;
 
 int fail2jet = 0;
 
  
 //==================== the event loop ============================================================
 while( ( evt = lcReader->readNextEvent()) != 0  && nEvents++ < maxEvt) {

   LCIterator<ReconstructedParticle> jets( evt, jetColName ) ;
   if (!jets()) continue;

   // for simplicity only consider events with two jets
   if( jets.size() != 2 ) {fail2jet++; continue;}

   // create a PID handler to access th b-tag information
   LCCollection* jetcol = evt->getCollection( jetColName ) ;
   PIDHandler *pidh = new PIDHandler(jetcol);
   int ilcfi = pidh->getAlgorithmID(btagAlgoName);
   int ibtag = pidh->getParameterIndex(ilcfi, btagParName);
   
   // loop over jets
   auto ajet = jets.next(); 
   while (ajet) {
   
     // get their b-tag value
     float btag = pidh->getParticleID(ajet, ilcfi).getParameters()[ibtag];
   
     hbtag->Fill( btag );
     
     ajet = jets.next(); 
   }

   
 }

 //===================================================================================================
 gStyle->SetOptStat(10);
 //===================================================================================================
  
 TCanvas* c1 = new TCanvas("b-tag plot","b-tag plot",800,600);
 hbtag->Draw() ;

 TString outnamepdf = outname + ".pdf";
 c1->SaveAs(outnamepdf); 
	
 TString outnameRoot  = outname + ".root";
 TFile outputFile (outnameRoot,"RECREATE");
 
 hbtag->Write();  
 
 outputFile.Close();
 //===================================================================================================

  std::cout << std::endl 
	    <<  "  " <<  nEvents 
	    << " events read from file: " 
	    << FILEN << std::endl  ;
  std::cout <<  " out of which " <<  fail2jet  
	    << " events don't have two jets " << std::endl  ;
  
  
  lcReader->close() ;
  delete lcReader ;
}
