# mini-DST: data format for easy analysis of ILC MC

The mini-DST provides easy access to simulated ILC data in a high-level format, which comprises e.g. ParticleFlow objects, isolated electrons, muons, taus, photons as well as jets.
The mini-DST files are readable in ROOT (as well as C++, Python, Go, Java and Julia) with the help the corresponding LCIO shared library / dictionary.

The mini-DST can be written from the full detector simulation and reconstruction chain, from [SGV](https://inspirehep.net/literature/1091154) and from Delphes via the [delphes2lcio](https://github.com/iLCSoft/LCIO/tree/master/examples/cpp/delphes2lcio) tool.

This page offers a short introduction and a collection of usage examples. You can download this package with 
```
git clone https://github.com/ILDAnaSoft/miniDST.git
```

# LCIO
LCIO is an event data model and persistency tool on which all Linear Collider data sets are based, allowing easy exchange of data at various stages of the simulation / reconstruction chain.

For installing LCIO, and for further reading, please consult [https://github.com/iLCSoft/LCIO](https://github.com/iLCSoft/LCIO).
Prerequisite is [cmake](https://cmake.org/) and - if you'd like to read LCIO in ROOT - a ROOT installation.

### set-up for OSG (or other systems with access to /cvmfs)
To setup the a consistent set of root, gcc, cmake, python etc for compiling LCIO on the OSG (or e.g. on the German NAF), please source [setenv4LCIO.sh](setenv4LCIO.sh):

```
. ./setenv4LCIO.sh
```
Then follow the instructions on [https://github.com/iLCSoft/LCIO](https://github.com/iLCSoft/LCIO) for downloading and installing LCIO with root dictionaries.

# Reading LCIO in ROOT
You need a ROOT environment and the LCIO library.
The environment variable LCIO is set either by sourcing the setup.sh script in your LCIO installation directory, or when initialising a full iLCSoft installation.

First, write the following lines in your `.rootlogon.C`
```
{
 gInterpreter->AddIncludePath("$LCIO");
 gSystem->Load("${LCIO}/lib/liblcio.so");
 gSystem->Load("${LCIO}/lib/liblcioDict.so");
}
```
and add $LCIO/lib to your LD_LIBRARY_PATH:
```
export LD_LIBRARY_PATH=$LCIO/lib:$LD_LIBRARY_PATH
```

After this, any [LCIO class](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1LCObject.html) can be used in your ROOT macro.
In context of the mini-DST, you'll only need to deal with the [ReconstructedParticle](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1ReconstructedParticle.html) and the [MCParticle](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1MCParticle.html).

# A first example macro
calculating the Higgs recoil mass on Higgsstrahlung events (e+e- ---> ZH) with Z->mumu, can be found in [./examples/higgs_recoil.C](./examples/higgs_recoil.C).
A mini-DST file to try out the example script can be downloaded from [here](https://desycloud.desy.de/index.php/s/5LmrjGWqziQfMe7). 
For comprehensive SM and Higgs MC data sets in mini-DST format, please consult [http://ilcsnowmass.org/](http://ilcsnowmass.org/). 

Once you have got hold of a mini-DST file FILENAME, you can type
```
.x higgs_recoil.C("FILENAME", "OUTPUTNAME")
```
in your ROOT session to get a plot of the recoil mass.

# A second example macro
creates the same plot, but this time with the ZH signal and e+e- --->mumujj background, dominated by e+e- ---> ZZ ---> mumujj, both normalised to given values for the integrated luminosity and the beam polarisations.
By standard, ILC event samples are generated with 100% beam polarisation, for
all allowed sign combinations (usually just the two opposite-sign combinations, P(e-,e+) = (-1,+1) and (+1,-1)).
Distributions for realistic polaristaion values are then created by weighting the events - [./examples/higgs_recoil_with_bkg.C](./examples/higgs_recoil_with_bkg.C) shows you how this works.

It reads four input mini-DST files:

* [ee -> ZH -> mumuH, P(e-,e+) = (-1,+1)](https://desycloud.desy.de/index.php/s/5LmrjGWqziQfMe7)
* [ee -> ZH -> mumuH, P(e-,e+) = (+1,-1)](https://desycloud.desy.de/index.php/s/3ZqPcGPELggW4bP)
* [ee -> ZZ -> mumujj, P(e-,e+) = (-1,+1)](https://desycloud.desy.de/index.php/s/9gKznqtSGcBKBWY)
* [ee -> ZZ -> mumujj, P(e-,e+) = (+1,-1)](https://desycloud.desy.de/index.php/s/3i3tj3adfMPfPaC)

Once you downloaded the input files and the macro, you can type 
```
.x higgs_recoil_bkg.C();
```
in your ROOT session to get the resulting plot.
In this case with quite a lot of background, since no cuts are applied (apart from two muons being present).
The macro optionally takes the following arguments with the following default values:
```
const char* DIRNAME = "./", double lumi_target=900., double epol_target=-0.8, double ppol_target=+0.3, TString outname = "recoil_plot"
```
where DIRNAME is the directory which hosts the input mini-DST files.

# Now it is your turn
Try to improve the signal-to-background ratio by applying a cut on the sum of the b-likeliness values of the two jets.
For this, read in the ```Refined2Jets``` collection, check that it is there and contains 2 jets, and then get the b-likeliness values (MVA output between 0 and 1).
You find an example of how to access jets and b-tag information in [./examples/jet_btag.C](./examples/jet_btag.C).
Take a look at this (of course you can also run it if you like!) and modify your ```higgs_recoil_with_bkg.C``` such that the recoil mass histograms are only filled if the sum of the two b-likeliness values > 1.

# mini-DST content

The following information is available on the mini-DST:

```
-------------------------------------------------------------------------------------------------
COLLECTION NAME         COLLECTION NAME       COLLECTION TYPE             EXPLANATION
(SGV / ILD full sim)    (Delphes) 
=================================================================================================
PandoraPFOs             PFOs                  ReconstructedParticle       particle flow objects from the main detector, incl. event shape variables
IsolatedElectrons       IsolatedElectrons     ReconstructedParticle         
IsolatedMuons           IsolatedMuons         ReconstructedParticle   
IsolatedTaus            IsolatedTaus          ReconstructedParticle
IsolatedPhotons         IsolatedPhotons       ReconstructedParticle
Refined2Jets            Durham2Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 2 jets (Durham algorithm, plus flavour tag)
Refined3Jets            Durham3Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 3 jets (Durham algorithm, plus flavour tag)
Refined4Jets            Durham4Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 4 jets (Durham algorithm, plus flavour tag)
Refined5Jets            Durham5Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 5 jets (Durham algorithm, plus flavour tag)
Refined6Jets            Durham6Jets           ReconstructedParticle       PandoraPFOs minus "IsolatedX" forced into 6 jets (Durham algorithm, plus flavour tag)
BCalPFOs                N/A                   ReconstructedParticle       particle flow objects from the most forward calorimeter
PrimaryVertex           N/A                   LCVertex                    
PrimaryVertex_RP        N/A                   ReconstructedParticle       "reconstructed particle" representing the primary vertex
MCParticlesSkimmed      MCParticles           MCParticle                    
MCTruthRecoLink         MCTruthRecoLink       LCRelation                  links from MCParticles to PandoraPFOs                 
RecoMCTruthLink         RecoMCTruthLink       LCRelation                  links from PandoraPFOs to MCParticles
------------------------------------------------------------------------------------------------
```
#### Useful LCIO utilities for seeing the detailed content of any LCIO file:
To see for each event in the file a list of collections with their number of elements:
```
anajob [your .slcio file]  
```
To print the content of the collections in event with number eventnumber:
```
dumpevent [your .slcio file] [eventnumber + 1]
```

# Usage examples

```
----------------------------------------------------------------------------------
PHYSICS QUANTITY         EXAMPLE MACRO                        EXPLANATION  
==================================================================================
(Higgs) recoil mass      ./examples/higgs_recoil.C           signal only
(Higgs) recoil mass      ./examples/higgs_recoil_with_bkg.C  luminosity and polarisation      
total visible energy
b-tag likelihood         ./examples/jet_btag.C               b-tag MVA output classifier
...

----------------------------------------------------------------------------------
```

