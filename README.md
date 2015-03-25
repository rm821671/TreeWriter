*TreeWriter* to build a tree from MiniAOD. Photon Cut-IDs and MVA are computed.

### Photons ### 
- based on this recipe for MVA ID: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1552.html)
  * `git cms-merge-topic …` necessary!
- this Cut-ID is included manually: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1541.html)


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

- generate dictionaries (not necessary, cmssw does so from LinkDef.h already)

    $ rootcint -f ROOT_dicts.cpp -c -I. TreeParticles.hpp LinkDef.h 

