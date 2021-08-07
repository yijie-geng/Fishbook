
import csv, os

os.chdir(r'') # set working directory, store files for analysis in this folder

outputFile = open(r'', 'w', newline='') # set output file name and path
outputWriter = csv.writer(outputFile)

global filecount, invalidFilecount, SocialScoreTotal
filecount = 0
invalidFilecount = 0
SocialScoreTotal = 0

for csvFilename in os.listdir('.'):
    pilotFrames = 0
    totalFrames = 0
    Ytotal = 0
    SocialScore = 0
    nancount = 0
    zerocount = 0
    if not csvFilename.endswith('.csv'):
            continue # skip non-csv files
    filecount += 1
    csvFileObj = open(csvFilename)
    readerObj = csv.reader(csvFileObj)
    for row in readerObj:
        pilotFrames += 1
        v = row

        if (pilotFrames > 4500): # videos were recorded at 450 frames per min; the 4500 frames setting will analyze 10 mins of recording 
            break
        totalFrames += 1
        if (v[0] == 'NaN'): 
            nancount += 1
        elif (v[0] == '0'):
            zerocount += 1
        else:
            y = float(v[0])
            Ytotal += ((325-y)/325) # total length of each test arena is 650 pixels in the recording

    csvFileObj.close()

    if totalFrames-nancount-zerocount == 0:
        invalidFilecount += 1
    else:
        SocialScore = Ytotal/(totalFrames-nancount-zerocount)

    SocialScoreTotal += SocialScore
    
    print(csvFilename)
    print('Social Score = ' + str(SocialScore))
    print('************************')

    outputWriter.writerow([str(csvFilename)])
    outputWriter.writerow(['Social Score = ', str(SocialScore)])
    outputWriter.writerow(['NaN = ', str(nancount)])
    outputWriter.writerow(['zero = ', str(zerocount)])
    outputWriter.writerow(['************************'])

SocialScoreAve = SocialScoreTotal/(filecount-invalidFilecount)

print('Social Score Ave = ', str(SocialScoreAve))

outputWriter.writerow(['Social Score Ave = ', str(SocialScoreAve)])


outputFile.close()


