from MSTreeDataImputationWorkflow import *

#parameters to select dataset to run
#tree a
#frumkin2005treeadataset="Frumkin2005_processedandreformateddata_Tue May 10 2016"
#frumkin2005treeadatasetv2="Frumkin2005_processedandreformateddata_Sun Aug 14 2016" # with only leaves
frumkin2005treeadatasettreea="Frumkin2005_treea_processedandreformateddata_Fri Feb 10 2017"
excelfilename="preprocessMSdataproperlyformateddatafortables2a.xlsx"
currentsheet='tables2a'

#tree b
frumkin2005treeadatasettreeb="Frumkin2005_treeb_processedandreformateddata_Tue Feb  16:05:27 2017"
excelfilenametreeb="preprocessMSdataproperlyformateddatafortables2b.xlsx"
currentsheettreeb='tables2b'

#tree c
frumkin2005treeadatasettreec="Frumkin2005_treec_processedandreformateddata_Tue Feb  16:06:23 2017"
excelfilenametreec="preprocessMSdataproperlyformateddatafortables2c.xlsx"
currentsheettreec='tables2c'


currentdataset=frumkin2005treeadatasettreea

saveddic=readinsaveddatastructs(currentdataset)
print("read in saved data structs")

impdic=readinexceltablesasdict(saveddic)
print("read in excel files with original datasets")

producetableswithmissingentries(impdic,200,[2**exp for exp in xrange(7)],rauvfm)
producetableswithmissingentries(impdic,200,[2**exp for exp in xrange(7)],rauvfmlbiaroc)
print("finished producing tables with missing entries ")

oldlistofimpmethod=[(onennimputationuniform,None),
                 (linearregressionendswithcutoff,None),
                 (Fitchforimputation,None)]

listofimpmethod=[(OneNN,None),
                 (LinearReg,None),
                 (FitchImp,None)]

#assert(False)

for impmethod in listofimpmethod:imputeandsavetableusingspecifiedimpmethod(impdic,impmethod[0],impmethod[1])
print("finished imputing blanked out tables")

statdic=makedicofallstats(impdic,[ratioofdiff,ratioofmagnitudediff])
print("finish making stat dictionary")

saveimputedata(currentdataset,impdic, statdic)
print("save impute dictionary and maybe statdic")
del impdic

print("display plots for imputed table data")
currentimpdic,statdic=readinsavedimpdatastructs(currentdataset)
turnstatdicsintoplots(statdic,colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"])
turnstatdicsintoplots_strippedaxis(statdic,colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"],coltitle=["Linear Regression","1NN Imputation","Fitch's Imputation"])

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
turntreestatdicsintoplots(treestatdic,"R-F",colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"])
turntreestatdicsintoplots_strippedaxis(treestatdic,"R-F",colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"])
turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlog(treestatdic,"R-F",colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"],hardyrange=None,coltitle=["Linear Regression","1NN Imputation","Fitch's Imputation"])
savetreedata(currentdataset,readintreecmpdata,treestatdic)

#ttu=impdic["data_for_tables2a"]["missingentriestablesdata_madeby_rauvfm"]["tableswith 64 missingentries"]["tables"][0];
#for i in ttu:print i;
#disttable,labels=makedisttable(ttu);nodedist,searchabletree=doneighborjoiningtreebuilding(disttable,labels,{});impttuassignmissingvaluesusingmodifiedfitch(ttu,searchabletree,labels);
