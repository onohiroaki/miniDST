# mini-DST: an LCIO-based format for easy analysis of ILC simulated data

Originally triggered by Snowmass 2021 comminity planning study, but also useful for beginners who usually are not familiar with high energy physics and detector simulation.

# Introduction

The purpose of mini-DST is to provide a "minumum" set of data from fully-simulated (or fast simulation based) MC samples.
When we perform physics analysis, we should use full detector simulation samples to make your analysis as realistic as possible, or, at least, fast simulation-based MC samples (like SGV-based and Delphes-based) are need to be used in the physics analysis.
Usually, these MC samples have tons amount of information.
These information are so important because these are the result of tracking, clustering, calibration, particle flow, and so on.
However, for beginners and newcomers who are typically not familiar with how to handle with it, the full information of simulation/reconstruction are too much and too complex.
The purpose of mini-DST is try to reduce such complexity for beginners.

# How To Use
The mini-DST files are directly readable with ROOT.
You need a ROOT environment and LCIO (Linear Collider I/O) library.



# Detail Description (WIP)

We have fully-simulated MC samples named DST file which contain all information and result of simulation and reconstruction.
In mini-DST file, the complex collections are removed, and useful information are added, e.g.; number of isolated electrons, muons, jets,...

In mini-DST file, the following collections are kept from original DST file.
- PandoraPFOs, BCalRecoParticle (not merged yet)
- MCParticle (MCParticlesSkimmed in next production)
- PrimaryVertex, PrimaryVertex_RP
- RecoMCThuthLink, MCTruthRecoLink (kept full relation so far)

The following collections/variables are added to the mini-DST file.
- event shape variables (used ThrustReconstruction, Sphere, Fox): these are stored at the header of PandoraPFOs
- IsolatedMuons, IsolatedElectrons (used IsolatedLeptonTagging, not tuned)
- IsolatedTaus (used TauFinder, not tuned)
- (IsolatedPhotons) (future)
- RefinedNJets (N = 2, 3, 4, 5, 6) (used LCFIPlus: JetClustering, JetVertexRefiner, FlavorTag, not tuned)
- ErrorFlow is applied to RefinedNJets to calculate covariance matrix for jets.
