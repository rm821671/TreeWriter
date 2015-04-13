**TreeWriter** to build a ROOT tree from MiniAOD. Photon Cut- and MVA-IDs are computed.

## Building and Running ##
Tested on lxplus:

Get CMSSW environment 7.2 or later

```
cmsrel CMSSW_7_2_0
cd CMSSW_7_2_0/src
cmsenv
```
Get and build egamma recipes

```
git cms-merge-topic ikrav:egm_id_phys14
git clone https://github.com/ikrav/ElectronWork.git
scram b -j 8
```
Get and build the TreeWriter

```
git clone https://github.com/cms-susy-photon-rwth-1b/TreeWriter.git
cd TreeWriter
scram b
```
Create Pilup Histograms

```
make -C PUreweighting
```
Run the TreeWriter

```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter.py
```

## Configure ##
in the python config, set
- `HT_cut`

## Objects ##
### Photons ###
- based on this recipe for MVA ID: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1552.html)
- this Cut-ID is included manually: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1541.html)
- all photons are used. boolean flags for: loose/medium/tight

### Jets ###
- ak4PFJetsCHS
- only loose id

### Muons ###
- fulfilling loose id
- tight id boolean flag

### Electrons ###
- fulfilling "veto" id
- boolean flags for loose/medium/tight
- recipes: [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_min)
