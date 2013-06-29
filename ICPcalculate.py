from __future__ import division as division
from pandas import DataFrame, read_csv
from numpy import random
import easygui as eg
import numpy as np
import os
import datetime
#import matplotlib.pyplot as plt
import pandas as pd


cancel = lambda x: exit() if x == None else None
PATH = os.getcwd()
configPATH = os.path.join(PATH, "\config"
inputLogin = eg.fileopenbox("", "Select the Sample Login file...", "*.csv")
cancel(inputLogin)
inputICPdata = eg.fileopenbox("", "Select the ICP export file...", "*.csv")
cancel(inputICPdata)
inputRawDataPDF = eg.fileopenbox("", "Select the PDF printout...", "*.pdf")
cancel(inputRawDataPDF)
outputFolder = eg.diropenbox("",'Choose the folder you want to save to...', '*.csv')
cancel(outputFolder)
#inputAnalytes = os.path.join(PATH, "Analytes.csv") #list of analytes to calculate

#LOQ/conversion data
inputICPsettings = os.path.join(configPATH, "ICPsettings.csv")
inputConversions = os.path.join(configPATH, "Conversions.csv") 

#with open("output.txt", "w") as myOutput:
	#myOutput.write("inputLogin\n")
	#myOutput.write(inputLogin)
	#myOutput.write("\n\ninputICPdata\n")
	#myOutput.write(inputICPdata)
	#myOutput.write("\n\ninputRawDataPDF\n")
	#myOutput.write(inputRawDataPDF)
	#myOutput.write("\n\ninputICPsettings\n")
	#myOutput.write(inputICPsettings)
	#myOutput.write("\n\ninputAnalytes\n")
	#myOutput.write(inputAnalytes)

#####Parse raw ICP export ddata #####
#Samples
#Sample Blank Correction lookup Table
#QC recovery lookup table [QCAverage]	
#High STD lookup table [HSTD]
#Low STD lookup table [LSD]
conversions = read_csv(inputConversions, sep=",")

dfExportICP = read_csv(inputICPdata, sep=",")
dfICPsorted = dfExportICP.sort_index()

dfICPsorted = dfICPsorted.applymap(lambda x: np.nan if str(x).isspace() else x)
dfICPsorted['Diluted To Vol.'].fillna(1, inplace=True)
dfICPsorted['Sample Prep Vol.'].fillna(1, inplace=True)
#dfICPsorted['Conc (Calib)'].fillna(99999, inplace=True)
dfICPsorted['Diluted To Vol.'] = dfICPsorted['Diluted To Vol.'].astype('int')
dfICPsorted['Sample Prep Vol.'] = dfICPsorted['Sample Prep Vol.'].astype('int')
dfICPsorted['Conc (Calib)'] = dfICPsorted['Conc (Calib)'].astype('float')
dfICPsorted['Analyte Name'] = dfICPsorted['Analyte Name'].apply(lambda x: x[0:2])
dfICPsorted['Analyte Name'] = dfICPsorted['Analyte Name'].apply(lambda x: x.strip(' '))
dfICPsorted['Analyte Name'] = dfICPsorted['Analyte Name'].apply(lambda x: x.strip('.'))
dfICPsorted['Total ug'] = dfICPsorted['Diluted To Vol.'] * dfICPsorted['Sample Prep Vol.'] * dfICPsorted['Conc (Calib)']



SamplesCriterion = dfICPsorted['Sample ID'].map(lambda x: str(x)[0:1].isalpha() == False)
HIGHcriterion = dfICPsorted['Sample ID'].map(lambda x: str(x)[0:4] == 'HIGH')
QBcriterion = dfICPsorted['Sample ID'].map(lambda x: str(x)[0:3] == 'QCB')
QC1criterion = dfICPsorted['Sample ID'].map(lambda x: str(x)[0:2] == 'QC')
LSDcriterion = dfICPsorted['Sample ID'].map(lambda x: str(x)[0:4] == 'LSD-')

#Create Grouped Dataframe chunks
VSTD = dfICPsorted[dfICPsorted['Sample ID'].isin(['ICV1','CCVF'])]
HSTD = dfICPsorted[HIGHcriterion]

QB = dfICPsorted[QBcriterion]
QC = dfICPsorted[QC1criterion]
QC2criterion = QC['Sample ID'].map(lambda x: str(x)[0:3] != 'QCB')
QC = QC[QC2criterion]



#QC Acceptance Criteria
#potentially make this based on QC historical data
defaultHSTD = 10
QClow = 0.9 
QChigh = 1.1
approveQC = lambda x, y, z: np.nan if ((x < y) or (x > z)) else x

###Create High STD lookup Table [HSTD]###
HSTD = dfICPsorted[HIGHcriterion]
HSTD['Description'] = HSTD['Description'].astype('int')

#Determine % recovery, then drop any that fails
HSTD['Recovery'] = HSTD['Total ug'] / HSTD['Description']
HSTD['Recovery'] = HSTD['Recovery'].map(lambda x: approveQC(x, QClow, QChigh))
HSTD.ix[HSTD['Recovery'].isnull(),['Description']] = defaultHSTD
#keep only highest passing STD
###HSTD###
HSTD = HSTD.groupby(['Analyte Name']).apply(lambda x: x.ix[x['Description'].idxmax()]).convert_objects(convert_numeric=True)



dfLogin = read_csv(inputLogin, sep=",")
dfLogin = dfLogin.sort_index()
#dfLogin.dtypes



dfLogin['OrderID'] = dfLogin['OrderID'].astype(str)
dfLogin.rename(columns={'SampleDetails.Method' : 'Method'}, inplace=True)
dfLogin = dfLogin.drop('Modified', axis=1)
#dfLogin.dtypes



Samples = dfICPsorted[SamplesCriterion]
Samples.rename(columns={'Sample ID' : 'SampleNumber'}, inplace=True)
Samples = Samples.drop('Description', axis=1)
#Samples



Samples['LookupMagic'] = 'Description'
Samples['High STD'] = HSTD.lookup(Samples['Analyte Name'], Samples['LookupMagic'])
Samples = Samples.drop('LookupMagic', axis=1)
Samples['TotalPrepVolume'] = Samples['Sample Prep Vol.'] * Samples['Diluted To Vol.']



#Check only Samples with dilutions to remove any with high concs
#This check will fail if there are multiple dilutions with conc all > HSTD
#But will leave a value if no dilutions were run even if sample conc > HSTD

while True:
    Samples = Samples.dropna()
    dupCheckFirst = Samples.duplicated(cols=['SampleNumber','Analyte Name'])
    dupCheckLast = Samples.duplicated(cols=['SampleNumber','Analyte Name'],take_last=True)
    dupCheck = dupCheckFirst|dupCheckLast
    dilutionsCheck = dupCheck & (Samples['Conc (Calib)'] > Samples['High STD'])
    Samples.loc[dilutionsCheck, 'Conc (Calib)'] = np.nan
    duplicateSize = Samples.loc[dilutionsCheck].index.size
    if duplicateSize == 0: break



Samples.sort()
Samples = Samples.groupby(['Analyte Name', 'SampleNumber']).apply(lambda x: x.ix[x['TotalPrepVolume'].idxmin()]).convert_objects(convert_numeric=True)



Samples = pd.merge(dfLogin, Samples, on='SampleNumber')
#Samples.dtypes



def mediaType(mediaShorthand):
    if mediaShorthand == 'MCE':
        return '0.8u MCE Filter'
    elif mediaShorthand == 'PVC':
        return '5u PVC Filter'
    elif mediaShorthand == 'IPE':
        return 'Wipe'
    else:
        return mediaShorthand
    
matrixExtraction = lambda x: mediaType(str(x)[len(str(x))-3:])




QC['Description'] = QC['Description'].astype('int')

#Parse media type from Sample ID
QC['media'] = QC['Sample ID'].map(lambda x: matrixExtraction(x))
QB['media'] = QB['Sample ID'].map(lambda x: matrixExtraction(x))
QB.rename(columns={'media' : 'Matrix'}, inplace=True)
QC.rename(columns={'media' : 'Matrix'}, inplace=True)
QB = QB.drop(['Diluted To Vol.', 'Sample Prep Vol.', 'Sample ID', 'Conc (Calib)', 'Wavelength', 'Description'], axis = 1)
#QB.dtypes



QB = QB.pivot(index='Analyte Name', columns='Matrix', values='Total ug')
QC['QCB'] = QB.lookup(QC['Analyte Name'], QC['Matrix'])
QC['Recovery'] = (QC['Total ug'] - QC['QCB'])/QC['Description']
QC['QC Factor'] = QC['Recovery'].apply(lambda x: approveQC(x, QClow, QChigh))

QCAverage = QC.pivot_table(values='QC Factor', rows='Matrix', cols='Analyte Name', aggfunc='mean')
QCstdev = QC.pivot_table(values='QC Factor', rows='Matrix', cols='Analyte Name', aggfunc='std')
#QCAverage.head


###STUB - Eventually call historical QC data to get bias and uncertainty
Samples['Recovery'] = QCAverage.lookup(Samples['Matrix'],Samples['Analyte Name'])
QC_stdev = QCstdev.lookup(Samples['Matrix'],Samples['Analyte Name'])
Samples['QC_uncertainty'] = 2 * QC_stdev / Samples['Recovery']
Samples['QC_bias'] = Samples['Recovery'] - 1
Samples['QC_Passed'] = Samples['Recovery'].apply(lambda x: 'False' if np.isnan(x) else 'True')
Samples['Recovery'] = Samples['Recovery'].apply(lambda x: 1 if np.isnan(x) else x)


LSD = dfICPsorted[LSDcriterion]
LOQs = read_csv(inputICPsettings)
LOQs.set_index('Element', inplace=True)
LOQs.sort_index(inplace=True)
LOQs['LOQ_filter'] = LOQs['LOQ_filter']/10
LOQs['LOQ_wipe'] = LOQs['LOQ_wipe']/20
#LOQs.dtypes

LSD['Matrix'] = 'LOQ_filter'
LSD['LOQ'] = LOQs.lookup(LSD['Analyte Name'],LSD['Matrix'])
LSD['Recovery'] = LSD['Conc (Calib)'] / LSD['Description']

LOQupper = 1.35
LOQlower = 0.65
acceptableLSD = (LSD['Description'] <= LSD['LOQ']) & (LSD['Recovery'] <= LOQupper) & (LSD['Recovery'] >= LOQlower)

LSD['Recovery'] = LSD['Recovery'].apply(lambda x: approveQC(x, 0.65, 1.35))
LSD = LSD.dropna()


LSDtable = LSD.loc[acceptableLSD]
goodLSDrows = LSDtable.groupby('Analyte Name')['Description'].idxmax()
goodLSDs = LSDtable.ix[goodLSDrows]



def matrixType(strCategory = 'LOQ|TLV', matrix = 'nothing'):
    if matrix == '5u PVC Filter' or matrix == '0.8u MCE Filter':
        return strCategory + '_filter'
    elif matrix == 'Wipe':
        return strCategory + '_wipe'
    else: return strCategory + '_filter'



goodLSDs.set_index('Analyte Name')


Samples['TestName'] = 'Test'
Samples['ConversionFactor'] = 1
#Samples['ConversionFactor'] = LOQs.lookup(Samples['Analyte Name'], Samples['ConversionFactor'])
Samples['TestName'] = LOQs.lookup(Samples['Analyte Name'], Samples['TestName'])
Samples['LOQ_Type'] = Samples['Matrix'].apply(lambda x: matrixType('LOQ',str(x)))
Samples['TLV_Type'] = Samples['Matrix'].apply(lambda x: matrixType('TLV',str(x)))
Samples['TLV'] = LOQs.lookup(Samples['Analyte Name'], Samples['TLV_Type']) #TLV in mg/m3 or wipe TLV in ug/s
Samples['LOQ'] = LOQs.lookup(Samples['Analyte Name'], Samples['LOQ_Type']) #ug/mL
Samples['LOQperSample'] = Samples['LOQ'] * Samples['Diluted To Vol.'] * Samples['Sample Prep Vol.']
notDetected = Samples['Total ug'] < Samples['LOQperSample']
Samples['notDetected'] = notDetected


###Split samples to handle blanks separately
isBlanks = Samples.groupby(Samples['IsBlank'] == 1)
notBlanks, Blanks = [Samples.ix[isBlanks.indices[tF]] for tF in (False, True)]



#Setup Blank lookup table with Blanks that pass QC criteria and are the smallest value per set/media
BlankFails = (Blanks['Total ug'] > (3 * Blanks['LOQ'] * Blanks['Sample Prep Vol.']) * Blanks['Diluted To Vol.']) #times prop volume
BlankTable = Blanks[['OrderID', 'Matrix', 'SampleNumber', 'Analyte Name', 'LOQ', 'LOQperSample','Total ug']]
BlankTable.loc[BlankFails,'Total ug'] = 0
BlankTable['Fails'] = BlankFails

negativeBlanks = Blanks['Total ug'] < 0
BlankTable.loc[negativeBlanks,'Total ug'] = 0
#BlankTable.head()

BlankTable = BlankTable.groupby(['OrderID','Matrix','Analyte Name']).apply(lambda x: x.ix[x['Total ug'].idxmin()]).convert_objects(convert_numeric=True)
#BlankTable.dtypes

Blanks.reset_index(drop=True,inplace=True)

Blanks['Results'] = Blanks['Total ug']
Blanks['Results'] = np.where(Blanks['notDetected'],Blanks['LOQperSample'],Blanks['Total ug'])
Blanks['Results'] = Blanks['Results'] / 1000 #convert to mg/sample

blankUnits = 'mg/sample' #in case this ever needs to be changed
Blanks['Units'] = blankUnits


#Assign blanks to non blanks
notBlanks.set_index(['OrderID','Matrix','Analyte Name'],inplace=True)
notBlanks['Blank_ug'] = BlankTable['Total ug']
notBlanks['BlankFailed'] = BlankTable['Fails']
notBlanks['LOQperSample'] = 0
notBlanks['LOQperSample'] = BlankTable['LOQperSample']
notBlanks.reset_index(inplace=True) #STUB:need to add try() only if multiindex

#Blank correct, QC correct, then
#Assign LOQ to sample if not detected
notBlanks['BlankCorrected_ug'] = (notBlanks['Total ug'] - notBlanks['Blank_ug'])
notBlanks['QC_Corrected_ug'] = notBlanks['BlankCorrected_ug'] / notBlanks['Recovery']
notDetected = notBlanks['QC_Corrected_ug'] < notBlanks['LOQperSample']
notBlanks['notDetected'] = notDetected
notBlanks['QC_Corrected_ug'] = np.where(notBlanks['notDetected'],notBlanks['LOQperSample'],notBlanks['QC_Corrected_ug'])

#Bias and uncertainty
notBlanks['SampleBias'] = notBlanks['QC_Corrected_ug'] * notBlanks['QC_bias']
notBlanks['SampleUncertainty'] = notBlanks['QC_Corrected_ug'] * notBlanks['QC_uncertainty']

#units in unicode messes up in excel pre2013
mg_m3 = 'mg/m'+u'\xb3'
ug_100 = u'\xb5'+'g/100cm'+u'\xb2'
isFilter = notBlanks['LOQ_Type'] == 'LOQ_filter'
isWipe = notBlanks['LOQ_Type'] == 'LOQ_wipe'

notBlanks['Units'] = ug_100
notBlanks['Results'] = notBlanks['QC_Corrected_ug']

notBlanks.loc[isFilter,'Units'] = mg_m3
filterResults = notBlanks.loc[isFilter,'QC_Corrected_ug'] / notBlanks.loc[isFilter,'Volume']
notBlanks.loc[isFilter,'Results'] = [filterResults]

notBlanks['PercentTLV'] = notBlanks['Results'] / notBlanks['TLV']



Blanks.set_index(['SampleNumber','Analyte Name'])
notBlanks.set_index(['SampleNumber','Analyte Name'])
Blanks.sort()
notBlanks.sort()

from math import log10, floor
def rnd_sigfigs(num,sig_figs=2):
    #assert isNumeric(num), 'x must be numeric!'
    if num != 0:
        rounded = round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
        if rounded >= 10: return str(int(rounded))
        else: return str(rounded)
    else:
        return str(0)  # Can't take the log of 0

simpleResults = pd.merge(Blanks,notBlanks,how='outer')

#Stoichiometric Conversions
convertedResults = DataFrame()
conversions.set_index('ToAnalyte', inplace=True)
for compound in conversions.index:
    element = conversions.loc[compound,'FromAnalyte']
    copyTest = simpleResults.loc[simpleResults['TestName'] == element]
    copyTest['TestName'] = compound
    copyTest['ConversionFactor'] = conversions.loc[compound,'Multiplier']
	#copyTest['TLV'] = conversions.loc[compound,'TLV']
    convertedResults = pd.concat([convertedResults, copyTest])

#potentially do this calc after combining simpleResults and convertedResults
#
convertedResults['LOQ'] = convertedResults['LOQ'] * convertedResults['ConversionFactor']
convertedResults['LOQperSample'] = convertedResults['LOQperSample'] * convertedResults['ConversionFactor']
convertedResults['Results'] = convertedResults['Results'] * convertedResults['ConversionFactor']

	
results = pd.concat([simpleResults,convertedResults])
results['percentTLV'] = results['Results'] / results['TLV']
	
	#copy base element group from results to memor
	#lookup multiplier
	#Convert result
	#append to converted dF

today = datetime.datetime.now()
resultsFile = 'Results-' + today.strftime('%Y%m%d') + '.csv'
results.to_csv(os.path.join(outputFolder,resultsFile),
							encoding='utf-8' ,float_format='%F',
							 index=False, index_label=False) #STUB:need to add try() while loop if blocked by permissions




