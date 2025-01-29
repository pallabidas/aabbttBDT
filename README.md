Installation:
```bash
cmsrel CMSSW_14_1_0_pre5
cd CMSSW_14_1_0_pre5/src
cmsenv
git cms-init
git cms-addpkg PhysicsTools/XGBoost 
git clone https://github.com/pallabidas/aabbttBDT.git
cd $CMSSW_BASE/src
scram b -j 12
```
Run using following command:
```bash
$CMSSW_BASE/bin/$SCRAM_ARCH/bdteval inputFile=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/AN-24-166/skkwan/condorSkim/2024-11-13-14h14m-year-2018-iteration0/TTTo2L2Nu/TTTo2L2Nu_0.root newOutputFile=1.0 newFile=test_skimfile.root
```
