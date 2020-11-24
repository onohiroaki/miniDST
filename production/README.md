# mini-DST production

This page will explain how to create mini-DST file. The detailed documentation is work in progress.

# How to create mini-DST
First, you have to initailize your iLCSoft environment, and prepare input files (MC samples of fully-simulated, SGV, ...).
Then specify the place of your input file as INPUTNAME, and output filename as OUTPUTNAME.
If you want to use IsolatedLeptonTagging and/or lcfiplus, you have to check which weight file/place you want to use to create mini-DST file.
Then, type following command and you are done.
'''
Marlin mini-DST-maker.xml
'''
Note that if your input files have very high multiplicity in an event, the job might take a few hours to complete.