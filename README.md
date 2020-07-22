# mini-DST: data format for easy analysis of ILC MC

The mini-DST provides easy access to simulated ILC data in a high-level format, which comprises e.g. ParticleFlow objects, isolated electrons, muons, taus, photons as well as jets. mini-DST files are readable in ROOT (as well as C++, Python, Go, Java and Julia) with the help the corresponding LCIO shared library / dictionary.

mini-DST can be written from the full detector simulation and reconstruction chain, from [SGV](https://inspirehep.net/literature/1091154) and from Delphes via the [delphes2lcio](https://github.com/iLCSoft/LCIO/tree/master/examples/cpp/delphes2lcio) tool.

This page offers a short introduction and a collection of usage examples.

# LCIO
LCIO is an event data model and persistency tool on which all Linear Collider data sets are based, allowing easy exchange of data at various stages of the simulation / reconstruction chain.

For installing LCIO, and for further reading, please consult [https://github.com/iLCSoft/LCIO](https://github.com/iLCSoft/LCIO).
Prerequisite is [cmake](https://cmake.org/) and - if you'd like to read LCIO in ROOT - a ROOT installation.

# Reading LCIO in ROOT
You need a ROOT environment and LCIO (Linear Collider I/O) library.
First, write the following lines in your `.rootlogon.C`.
```
{
 gInterpreter->AddIncludePath("$LCIO");
 gSystem->Load("${LCIO}/lib/liblcio.so");
 gSystem->Load("${LCIO}/lib/liblcioDict.so");
}
```
After this, any [LCIO class](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1LCObject.html) can be used in your ROOT macro.
In context of the mini-DST, you'll only need to deal with the [ReconstructedParticle](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1ReconstructedParticle.html) and the [MCParticle](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1MCParticle.html).

A first example macro, calculating the Higgs recoil mass on Higgsstrahlung events with Z->mumu, can be found in [./examples/higgs_recoil.C](./examples/higgs_recoil.C). An example mini-DST file to try out the example script can be downloaded from [here](https://desycloud.desy.de/index.php/s/KyNEFt3SFFEk6aM). 
For comprehensive SM and Higgs MC data sets in mini-DST format, please consult [http://ilcsnowmass.org/](http://ilcsnowmass.org/). 

Once you have got hold of a mini-DST file FILENAME, you can type

```
.x higgs_recoil.C("FILENAME", "OUTPUTNAME")
```
in your ROOT session to get a plot of the recoil mass.


# mini-DST content

The following information is available on the mini-DST:

```
-------------------------------------------------------------------------------------------------
COLLECTION NAME         COLLECTION NAME       COLLECTION TYPE             EXPLANATION
(SGV / ILD full sim)    (Delphes) 
=================================================================================================
PandoraPFOs             PFOs                  ReconstructedParticle       particle flow objects from the main detector, incl. event shape variables
IsolatedElectrons       Electrons             ReconstructedParticle         
IsolatedMuons           Muons                 ReconstructedParticle   
IsolatedTaus            Taus                  ReconstructedParticle
IsolatedPhotons         Photons               ReconstructedParticle
Refined2Jets            Durham2Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 2 jets (Durham algorithm, plus flavour tag)
Refined3Jets            Durham3Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 3 jets (Durham algorithm, plus flavour tag)
Refined4Jets            Durham4Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 4 jets (Durham algorithm, plus flavour tag)
Refined5Jets            Durham5Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 5 jets (Durham algorithm, plus flavour tag)
Refined6Jets            Durham6Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 6 jets (Durham algorithm, plus flavour tag)
BCalPFOs                N/A                   ReconstructedParticle       particle flow objects from the most forward calorimeter
PrimaryVertex           N/A                   LCVertex                    
PrimaryVertex_RP        N/A                   ReconstructedParticle       "reconstructed particle" representing the primary vertex
MCParticles(Skimmed)    MCParticle            MCParticle                    
MCTruthRecoLink         MCTruthRecoLink       LCRelation                  links from MCParticles to PandoraPFOs                 
RecoMCTruthLink         RecoMCTruthLink       LCRelation                  links from PandoraPFOs to MCParticles
------------------------------------------------------------------------------------------------
```

# usage examples

```
---------------------------------------------------------------------------
PHYSICS QUANTITY             EXAMPLE MACRO                   EXPLANATION  
===========================================================================
(Higgs) recoil mass          ./examples/higgs_recoil.C         
total visible energy
b-tag likelihood
...

---------------------------------------------------------------------------
```


