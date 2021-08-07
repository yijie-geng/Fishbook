
import csv, os

os.chdir(r'') # set working directory, store files for analysis (the Fishbook script outputs) in this folder 

outputFile = open(r'', 'w', newline='') # set output file name and path
outputWriter = csv.writer(outputFile)

for csvFilename in os.listdir('.'):
        listLaneNo = ['Lane Number']
        listSocialScore = ['Social Score']
        
        if not csvFilename.endswith('.csv'):
            continue # skip non-csv files
        csvFileObj = open(csvFilename)
        readerObj = csv.reader(csvFileObj)
        Data = list(readerObj)

        for i in range(0, len(Data), 5):
                listLaneNo.append(Data[i][0])
        for i in range(1, len(Data), 5):
                listSocialScore.append(Data[i][1])

        csvFileObj.close()

        outputWriter.writerow([str(csvFilename)])
        outputWriter.writerow(listLaneNo)
        outputWriter.writerow(listSocialScore)
        outputWriter.writerow(['************************'])

outputFile.close()


