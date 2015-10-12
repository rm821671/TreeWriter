**TreeWriter** to build a ROOT tree from MiniAOD. Photon Cut- and MVA-IDs are computed.

## Building and Running ##
Get CMSSW environment 74X

```
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src/
cmsenv
git cms-merge-topic ikrav:egm_id_7.4.12_v1
```
Get and build the TreeWriter

```
git clone https://github.com/cms-susy-photon-rwth-1b/TreeWriter.git
scram b
cd TreeWriter
```
Create Pilup Histograms

```
make -C PUreweighting
```
Run the TreeWriter
- locally
```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter.py
```
- on the Grid using CRAB3
```
. /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab
```
for a single dataset
```
crab submit -c crabConfig.py
```
for all datasets
```
python2 crabConfig.py
```

## Configure ##
in the python config, set
- `HT_cut`: minimum HT
- `photon_pT_cut`: minimum leading-photon pT

## Objects ##
### Photons ###
- official cut-based ID and general purpose MVA are taken from [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)
- all photons are used. boolean flags for: loose/medium/tight
- general purpose MVA value is stored

### Jets ###
- ak4PFJetsCHS
- all jets are used
- boolean flag for: loose
- boolean flags for whether a loose electron/photon is found within dR<0.4

### Muons ###
- fulfilling loose id
- tight id boolean flag

### Electrons ###
- fulfilling "veto" id
- boolean flags for loose/medium/tight
- recipes on [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)

### Generated Particles ###
- genJets collection is stored (= full slimmedGenJets)
- gen[Electrons|Photons]: status=1, pT>30
- Decay products (daughters) of W bosons
