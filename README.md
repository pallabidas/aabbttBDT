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

