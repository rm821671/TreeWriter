#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt

PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 60 --numPileupBins 60 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON)

all: data/puWeights.root

data:
	mkdir -p $@

data/puWeights.root: data createWeights.py nTrueVertexMC.root nTrueVertexData.root
	python2.7 createWeights.py --mc nTrueVertexMC.root --data nTrueVertexData.root -o $@

nTrueVertexMC.root: createMChist.py
	python2.7 $^ $@

# for the pileup uncertainties shift cross section by 20% which is only a rough estimation,
# see https://hypernews.cern.ch/HyperNews/CMS/get/hlt/3453/1.html
nTrueVertexData.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 80000 --pileupHistName pileup tmpHist1.root
	$(PU_CALC) --minBiasXsec 96000 --pileupHistName pileupUp tmpHist2.root
	$(PU_CALC) --minBiasXsec 64000 --pileupHistName pileupDown tmpHist3.root
	hadd -f $@ tmpHist1.root tmpHist2.root tmpHist3.root
	rm tmpHist1.root tmpHist2.root tmpHist3.root




