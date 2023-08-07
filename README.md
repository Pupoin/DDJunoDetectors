# RUN
Package required: [DD4hep](https://github.com/AIDASoft/DD4hep), 
[Root](https://root.cern.ch/),
[Geant4](https://geant4.web.cern.ch/)

## Environment needed
```bash 
# on cern
source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh
# or using anaconda
conda create -c conda-forge --name <envName> root python geant4 make
```
```bash
# then compile DD4hep and souce its environment
source  DD4hep/install/path/bin/thisdd4hep.sh
# if on lxplus7, just do it.
source /workfs2/bes/yuanchy8/DD4hep_G4/bin/thisdd4hep.sh
```

## How to run me
```bash
git clone https://github.com/Pupoin/DDJunoDetectors.git
cd DDJunoDetectors
mkdir build
cd build
cmake ..
make -j 

cd ../compact/
# show the display
geoDisplay Juno.xml
# check the overlaps
checkOverlaps Juno.xml

```