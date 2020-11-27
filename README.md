# LLPtriggertiming

Studies for LLP triggers using timing in low-HT events

# Installation
Specify '--recurse-submodules' when cloning, or remember to init+update the submodules
Delphes should show the tag '3.4.2' and Pythia 'pythia8235'

```
cd Pythia
make -j9
cd ../Delphes
patch -p1 < ../Delphes.patch
make HAS_PYTHIA8=true  -j9
cd ..
```

For local runs of Delphes `source setup.sh`, for all the analysis plots `source setup_analysis.sh`

# Generate minbias

Generate individual events to overlay
```
cd run_condor
condor_submit condorSubmit.sub #send 1000 jobs x 6000 events
# wait for jobs
cd ../Delphes
./root2pileup ../LLPtrigger_samples/Minbias.pileup ../LLPtrigger_samples/Minbias.10220229*
#at this point can remove the individual Minbias files
cd ../run_condor
#make sure the delphes card points to the correct Minbias.pileup file within Delphes/cards/delphes_card_ATLAS_PileUp.tcl
condor_submit condorSubmitWithPU.sub #send 1000 jobs x 1000 events
hadd MinbiasWithPU55.root MinbiasWithPU55.10220303.*
#remove the individual files
```

# Check trigger rates
The minbias+PU file can be used to verify the trigger rates. The single jet rates are prety decent but 4J15 is clearly off.
```
cd Analysis
python plot_samples_withPU.py
```

# Shower existing LHE files

This is not strictly needed, and not fully correct since I don't do the proper Pythia matching. It's only useful to have a sample with higher pT spectrum to derive later on the truth2reco calibration
```
cd Delphes
#check the LHE file within examples/Pythia8/configLHE.cmnd
./DelphesPythia8 cards/delphes_card_ATLAS_PileUp.tcl examples/Pythia8/configLHE.cmnd ../LLPtrigger_samples/n3n4_1TeVsquark_150gev_showerLHE_withPU55.20k.root
#set showered file as input for
python truth_to_reco.py
```
The output plots are the truth-reco calibration and relative resolution


# Derive calibrations
```
cd Analysis
python truth_to_reco.py
#check the plot h2d.png, should show something like a -4 GeV offset between truth and reco
python reco_to_l1.py
#check the plot L1Jfit.png, should show a good fit to the turn-on curves
#the four fitted parameters are the calibration constants (reco pT to L1 pT) and the values of the smearing function
#those should be ported to plot_daod.py and convertL1Jet.C
```

#Analysis plots
The analysis works on DAOD_TRUTH3 for any signal sample. Select the signal in the map at the top of the file
```
cd Analysis
python plot_daod.py
```
