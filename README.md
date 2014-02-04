HiggsAnalysis-CombinedLimit
===========================

[Manual to run combine](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#How_to_run_the_tool)

### Recipe to check out RazorCMS version of combine
```
cmsrel CMSSW_6_1_2
cd CMSSW_6_1_2/src/
cmsenv
git clone https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git pull origin razor1dpdf
scramv1 b clean; scramv1 b
```
