*TreeWriter* to build a tree from MiniAOD. Photon Cut-IDs and MVA are computed.

### Photons ### 
- based on this recipe for MVA ID: https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1552.html
  * "git cms-merge-topic …" necessary!
- this Cut-ID is included manually: https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1541.html


### Jets ###
- ak4PFJetsCHS


- generate dictionaries
$ rootcint -f ROOT_dicts.cpp -c -I. TreeParticles.hpp LinkDef.h 
