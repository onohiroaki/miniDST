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
 
void higgs_recoil_with_bkg(const char* DIRNAME = "./", double lumi_target=900., double epol_target=-0.8, double ppol_target=+0.3, TString outname = "recoil_plot_with_bkg") {

 int maxEvt = 10000 ;  // change as needed

//---- create histograms ----

 int nbin = 120;
 float xmin = 60;
 float xmax = 180;

 const int nhistos = 2;
 int icol[nhistos] = {2, 4}; 
 
 TH1F* hrecoilm[nhistos]; 
 hrecoilm[0] = new TH1F("hrecoilm_e2e2h",";recoil mass [GeV]; ; ", nbin, xmin , xmax ) ;
 hrecoilm[1] = new TH1F("hrecoilm_zz_sl",";recoil mass [GeV]; ; ", nbin, xmin , xmax ) ;
 TString legtext[nhistos];
 legtext[0] = "e^{+}e^{-} #rightarrow #mu^{+}#mu^{-} Higgs";
 legtext[1] = "e^{+}e^{-} #rightarrow ZZ #rightarrow #mu^{+}#mu^{-} jj";

 //------ collection to use: 
 //------ ALWAYS read also PandoraPFOs when reading IsolatedX collections! -------

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
 
 const int nfiles = 4;
 TString FILEN[nfiles] = {"rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106479.Pe2e2h.eL.pR-00001-ILDminiDST.slcio",
                          "rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106480.Pe2e2h.eR.pL-00001-ILDminiDST.slcio",
                          "rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106575.P4f_zz_sl.eL.pR-00001-ILDminiDST.slcio",
                          "rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106576.P4f_zz_sl.eR.pL-00001-ILDminiDST.slcio"};

 int nevt_per_file[nfiles]; 
 int histindex_of_file[nfiles] = {0, 0, 1, 1}; 
  
 for (int ifile = 0; ifile < nfiles; ifile++) {
   
   FILEN[ifile] = DIRNAME + FILEN[ifile];
   nevt_per_file[ifile] = 0;
   
   cout << "reading file " << FILEN[ifile] << endl;
   lcReader->open( (const char *)FILEN[ifile] ) ;
  
   EVENT::LCEvent* evt ; 
   
   int nEvents = 0;
  
   //==================== set number of events ============================================================
   nevt_per_file[ifile] = min(maxEvt, lcReader->getNumberOfEvents());   
     
   //==================== the event loop ============================================================
   while( ( evt = lcReader->readNextEvent()) != 0  && nEvents++ <  nevt_per_file[ifile]) {
        
     // CrossSection_fb is the event parameter storing the generator cross-section for this event
     float xsection = evt->parameters().getFloatVal( "CrossSection_fb");
     
     // beamPol1 & beamPol2 are the event parameters storing the generated beam polarisations of this event
     // = 1 = 100% = right-handed  / -1 = -100% = left-handed
     //  polarisation weights for {LR, RL, LL, RR} events, as example for target P(e-,e+)=(-80%,+30%):
     //  LR: polweight = (1-epol_target)*(1+ppol_target)/4.; // -80%,+30% => 1.8 * 1.3 / 4. = 0.585
     //  RL: polweight = (1+epol_target)*(1-ppol_target)/4.; // -80%,+30% => 0.2 * 0.7 / 4. = 0.035
     //  LL: polweight = (1-epol_target)*(1-ppol_target)/4.; // -80%,+30% => 1.8 * 0.7 / 4. = 0.315
     //  RR: polweight = (1+epol_target)*(1+ppol_target)/4.; // -80%,+30% => 0.2 * 1.3 / 4. = 0.065
     float epol = evt->parameters().getFloatVal( "beamPol1");
     float ppol = evt->parameters().getFloatVal( "beamPol2");
     float polweight = (1 + epol*epol_target)*(1 + ppol*ppol_target)/4.;
     
     // assume one file per process-polarisation combination 
     // - in case of several files per process-polarisation combination, need to add up *beforehand* 
     // the numbers of events of all the files belonging to the same process-polarisation combination! 
     float weight = polweight * xsection * lumi_target / nevt_per_file[ifile];
     if (nEvents==1) {
       cout << "weight = " << weight << endl;
     }
   
     LCIterator<ReconstructedParticle> muons( evt, muoColName ) ;

     // for simplicity only consider events with two muons
     if( muons.size() != 2 ) continue;

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
     hrecoilm[histindex_of_file[ifile]]->Fill( recoil.M(), weight ) ;
   
   }  // end of event loop per file
       
   std::cout << std::endl 
	    <<  "  " <<  nevt_per_file[ifile]
	    << " events read from file: " 
	    << FILEN[ifile] << std::endl  ;

 } // end of file loop 

 //===================================================================================================
 gStyle->SetOptStat(0);
 //===================================================================================================
  
 TCanvas* c1 = new TCanvas("recoil plots","recoil plot",800,600);
 
 THStack *hs = new THStack("hs",";recoil mass [GeV]; ; ");

 TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
 leg->SetHeader("ILC at 250 GeV");

 for (int ihist = nhistos-1; ihist > -1; ihist--) {
   hrecoilm[ihist]->SetFillColor(icol[ihist]);
   //hrecoilm[ihist]->SetMarkerColor(icol[ihist]);
   //hrecoilm[ihist]->SetLineColor(icol[ihist]);
   hs->Add(hrecoilm[ihist]) ;
   leg->AddEntry(hrecoilm[ihist], legtext[ihist], "F");
 }  
 
 hs->Draw("hist");
 leg->Draw();
   
 TString outnamepdf = outname + ".pdf";
 c1->SaveAs(outnamepdf);
 
 
 TString outnameRoot  = outname + ".root";
 TFile outputFile (outnameRoot,"RECREATE");
 
 for (int ihist = 0; ihist < nhistos; ihist++) {
   hrecoilm[ihist]->Write();  
 } 
 
 outputFile.Close();
 //===================================================================================================

  
  
  lcReader->close() ;
  delete lcReader ;
}
