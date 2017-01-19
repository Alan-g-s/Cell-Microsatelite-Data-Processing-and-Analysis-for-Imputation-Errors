from MSTreeDataImputationWorkflow import *

#parameters to select dataset to run
frumkin2005treeadataset="Frumkin2005_processedandreformateddata_Tue May 10 2016"
frumkin2005treeadatasetv2="Frumkin2005_processedandreformateddata_Sun Aug 14 2016" # with only leaves
excelfilename="preprocessMSdataproperlyformateddatafortables2a.xlsx"
currentsheet='tables2a'

currentdataset=frumkin2005treeadatasetv2

saveddic=readinsaveddatastructs(currentdataset)
print("read in saved data structs")

impdic=readinexceltablesasdict(saveddic)
print("read in excel files with original datasets")

producetableswithmissingentries(impdic,200,[2**exp for exp in xrange(7)],rauvfm)
producetableswithmissingentries(impdic,200,[2**exp for exp in xrange(7)],rauvfmlbiaroc)
print("finished producing tables with missing entries ")

listofimpmethod=[(onennimputationuniform,None),
                 (linearregressionendswithcutoff,None),
                 (Fitchforimputation,None)]

for impmethod in listofimpmethod:imputeandsavetableusingspecifiedimpmethod(impdic,impmethod[0],impmethod[1])
print("finished imputing blanked out tables")

statdic=makedicofallstats(impdic,[ratioofdiff,ratioofmagnitudediff])
print("finish making stat dictionary")

saveimputedata(currentdataset,impdic, statdic)
print("save impute dictionary and maybe statdic")
del impdic

print("display plots for imputed table data")
currentimpdic,statdic=readinsavedimpdatastructs(currentdataset)
turnstatdicsintoplots(statdic,colorder=["linearregressionendswithcutoff_imputedtables","onennimputationuniform_imputedtables","Fitchforimputation_imputedtables"])
turnstatdicsintoplots_strippedaxis(statdic,colorder=["linearregressionendswithcutoff_imputedtables","onennimputationuniform_imputedtables","Fitchforimputation_imputedtables"],coltitle=["Linear Regression","1NN Imputation","Fitch's Imputation"])

seefilter,colcounts=makefiltereddata(currentimpdic,basicfilter)
savefilterdata(currentdataset,seefilter,colcounts)
del seefilter
del colcounts

print("finished filtering tables to remove duplicate cols")

filterdata,coldic=readinsavedfilterdatastructs(currentdataset)
#this function produces trees from the tables and prepares them for analysis
readintreecmpdata,titles=readinandsafetreecmps(filterdata,coldic,runallweightnjfortree)
dataseperateddic={}
print titles[4:len(titles)]

treestatdic=makedicofallstatsfortrees(readintreecmpdata,[ratioofdiff,ratioofmagnitudediff],titles,4,len(titles))

print("finished making tree data")

print("display plots for tree data")
turntreestatdicsintoplots(treestatdic,"R-F",colorder=["linearregressionendswithcutoff_imputedtables","onennimputationuniform_imputedtables","Fitchforimputation_imputedtables"])
turntreestatdicsintoplots_strippedaxis(treestatdic,"R-F",colorder=["linearregressionendswithcutoff_imputedtables","onennimputationuniform_imputedtables","Fitchforimputation_imputedtables"])
turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlog(treestatdic,"R-F",["linearregressionendswithcutoff_imputedtables","onennimputationuniform_imputedtables","Fitchforimputation_imputedtables"],hardyrange=None,coltitle=["Linear Regression","1NN Imputation","Fitch's Imputation"])
savetreedata(currentdataset,readintreecmpdata,treestatdic)
