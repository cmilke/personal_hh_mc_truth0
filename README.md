# hh_mc_truth0

Originally from https://gitlab.cern.ch/xiaohu/readtruth0

Common tool for HH truth studies with the TRUTH0 format

Contact Xiaohu.Sun@cern.ch

# Make evgen files in TRUTH0

## Event generation

```
setupATLAS
asetup MCProd,20.7.9.9.27 # choose the correct release that your model and JO require
Generate_tf.py --jobConfig yourJO.py --maxEvents 5000 --runNumber 333333 --outputEVNTFile test.EVNT.pool.root --ecmEnergy 13000. —randomSeed=1283
```

## Convert to TRUTH0

```
asetup 21.2.6.0,AthDerivation
Reco_tf.py --inputEVNTFile test.EVNT.pool.root --outputDAODFile test.EVNT.TRUTH0.pool.root --reductionConf TRUTH0
```

# Make histograms with TRUTH0

## For the first time run the package
```
git clone ssh://git@gitlab.cern.ch:7999/xiaohu/readtruth0.git
cd readtruth0
setupATLAS
lsetup 'rcsetup Base,2.4.28'
rc find_packages
rc compile
```

## For every time after only needs below
```
cd readtruth0
setupATLAS
source rcSetup.sh
```

Run example
```
. run_example.sh
```
The histograms are stored in job-bbA_mA800_mH500/hist-*.root

# Develop a new class for your final state

Using 4b as an example, one needs to create or modify the following files accordingly as a full set
```
mytruth/Root/LinkDef.h # for CINT
mytruth/Root/readtruth_hh4b.cxx # for the source code
mytruth/mytruth/readtruth_hh4b.h # for the header
mytruth/util/run_readtruth.cxx # for the binary excutable
```
# Other information

For more information on reading truth, see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TruthDAOD

