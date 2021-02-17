"""Code for calculating rPSA low,medium,high vulnerability categories.
Compares results to sPSA and outputs in figure.
Prints number of species categorized in low, medium, high for sPSA, then rPSA.
Prints proportion and number of species changing categories with rPSA compared to sPSA.
"""

from matplotlib import pyplot as plt
import numpy as np
import csv
import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import norm
from scipy.stats import rankdata

SMALL_SIZE = 36
MEDIUM_SIZE = 48
BIGGER_SIZE = 72

def mscatter(x,y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    ax = ax or plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc





def geo_mean(iterable): # calculate a geometric mean
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

def inputFile(path1, path2): # input file from paths and read in data for productivity and susceptibility
    with open(path1, newline='') as csvfile:  # read in susceptibility data and store in data
        data = list(csv.reader(csvfile))

    with open(path2, newline='') as csvfile:  # read in productivity data and store in data2
        data2 = list(csv.reader(csvfile))

    return data, data2

def variance(prodPercentiles, suscPercentiles): # output associated variances related to underlying analysis
    """find common denominator to calculate multipliers for variance calculations"""
    tempProd = max([1 / eval(prodPercentiles[0]), 1 / (eval(prodPercentiles[1]) - eval(prodPercentiles[0])),
                    1 / (1 - eval(prodPercentiles[1]))])
    tempSusc = max([1 / eval(suscPercentiles[0]), 1 / (eval(suscPercentiles[1]) - eval(suscPercentiles[0])),
                    1 / (1 - eval(suscPercentiles[1]))])
    tempBoth = [tempProd, tempSusc]
    permBoth = [0, 0]
    counted = [0, 0]
    for i in range(2):
        counting = 1
        flag = True
        while flag == True:
            integerCheck = tempBoth[i] * counting
            counting += 1
            if integerCheck % 1 == 0:
                flag = False
                permBoth[i] = integerCheck
                counted[i] = counting

    multipliersProd = [permBoth[0] / eval(prodPercentiles[0]),
                       permBoth[0] / (eval(prodPercentiles[1]) - eval(prodPercentiles[0])),
                       permBoth[0] / (1 - eval(prodPercentiles[1]))]
    multipliersSusc = [permBoth[1] / eval(suscPercentiles[0]),
                       permBoth[1] / (eval(suscPercentiles[1]) - eval(suscPercentiles[0])),
                       permBoth[1] / (1 - eval(suscPercentiles[1]))]

    addVarianceProd = (multipliersProd[0] + multipliersProd[2]) / sum(multipliersProd)
    addVarianceSusc = (multipliersSusc[0] + multipliersSusc[2]) / sum(multipliersSusc)

    multMean = np.log(6) / 3
    multVarianceProd = (multipliersProd[0] * (multMean) ** 2 + multipliersProd[1] * (multMean - np.log(2)) ** 2 +
                        multipliersProd[2] * (multMean - np.log(3)) ** 2) / sum(multipliersProd)
    multVarianceSusc = (multipliersSusc[0] * (multMean) ** 2 + multipliersSusc[1] * (multMean - np.log(2)) ** 2 +
                        multipliersSusc[2] * (multMean - np.log(3)) ** 2) / sum(multipliersSusc)

    return [addVarianceProd, addVarianceSusc, multVarianceProd, multVarianceSusc]

def addAnalysis(data, num, var, type): # additive model output scores and associated standard errors
    error = []
    scores = []
    for i in range(len(data)):
        prescore=[]
        weight=[]
        for j in range(1, 2+num):
            if data[i][j] != '':
                weight.append(float(data[i][j + num + 2]))
                for b in range(int(data[i][j+num+2])):
                    if type == 'prod':
                        prescore.append(4-float(data[i][j]))
                    else:
                        prescore.append(float(data[i][j]))
        if prescore != []:
            scored = np.average(prescore)
            scores.append(scored)
            totalWeight = np.sum(weight)

            for z in range(len(weight)):
                weight[z] = (weight[z]/totalWeight)**2
            error.append((var*np.sum(weight))**(1/2))
    return scores, error

def multAnalysis(data, num, var, type): # multiplicative model output scores and associated standard errors
    error = []
    scores = []
    logscores = []
    for i in range(len(data)):
        prescore=[]
        weight=[]
        for j in range(1, 2+num):
            if data[i][j] != '':
                weight.append(float(data[i][j + num + 2]))
                for b in range(int(data[i][j+num+2])):
                    if type == 'prod':
                        prescore.append(4-float(data[i][j]))
                    else:
                        prescore.append(float(data[i][j]))
        if prescore != []:
            totalWeight = np.sum(weight)
            scored = geo_mean(prescore)
            logscores.append(np.log(scored))
            scores.append(scored)


            for z in range(len(weight)):
                weight[z] = (weight[z]/totalWeight)**2
            error.append((var*np.sum(weight))**(1/2))
    return scores, error, logscores


def main():

    abs_file_path = "C:/Users/Richard/Desktop/desktop/PSAReanalysis/NOAA_UnitedStates/NOAA_Susceptibility.csv" # path for susceptibility spreadsheet
    abs_file_path2 = "C:/Users/Richard/Desktop/desktop/PSAReanalysis/NOAA_UnitedStates/NOAA_Productivity.csv" # path for productivity spreadsheet


    data, data2 = inputFile(abs_file_path, abs_file_path2)

    prodScaling = data2[0][1] # store scaling choice (multiplicative or additive) for productivity attributes
    suscScaling = data[0][1] # store scaling choice (multiplicative or additive) for susceptibility attributes
    scaling = [prodScaling, suscScaling]
    prodPercentiles = data2[1][1:3] # store percentile cut-offs for attribute scores (1,2,3) for productivity attributes
    suscPercentiles = data[1][1:3] # store percentile cut-offs for attribute scores (1,2,3) for susceptibility attributes
    prodThresholds = data2[1][1:3] # store productivity thresholds
    suscThresholds = data[1][1:3] # store susceptibility thresholds
    if prodThresholds != suscThresholds:
        print('Warning: Thresholds in Productivity spreadsheet do not match thresholds in Susceptibility spreadsheet!')

    varianceList = variance(prodPercentiles, suscPercentiles)

    numProd = int(data2[2][1]) # store number of productivity attributes used in analysis
    numSusc = int(data[2][1]) # store number of susceptibility attributes used in analysis


    data = data[6:]  # exclude first row (data header)
    data2 = data2[6:]

    newdata = []
    newdata2 =[]

    for x in range(len(data)): # exclude first column
        newdata.append(data[x][1:])
        newdata2.append(data2[x][1:])

    mean = [2,2]

    """define scaling model based on first letter in scaling input"""

    if scaling[0][0] == 'a' or scaling[0][0] == 'A':
        productivity, producError = addAnalysis(newdata2, numProd, varianceList[0], 'prod')
        mean[0] = 2
        choiceProd = productivity
    else:
        productivity, producError, logproductivity = multAnalysis(newdata2, numProd, varianceList[2], 'prod')
        mean[0] = np.log(6)/3
        choiceProd = logproductivity

    if scaling[1][0] == 'a' or scaling[1][0] == 'A':
        susceptibility, susceptError = addAnalysis(newdata, numSusc, varianceList[1], 'susc')
        mean[1] = 2
        choiceSus = susceptibility
    else:
        susceptibility, susceptError, logsusceptibility = multAnalysis(newdata, numSusc, varianceList[3], 'susc')
        mean[1] = np.log(6)/3
        choiceSus = logsusceptibility



    SEp = producError
    SEs = susceptError
    SEps = []
    riskVector = []
    projectionMatrix =[]
    mean = np.matrix(mean).transpose()
    transformation=[]
    projection = []
    projection_p = []
    projection_s = []
    distanceMetric = []
    revisedCategory = []
    low = 0
    medium = 0
    high = 0

    lowerThresh = norm.ppf(eval(prodThresholds[0]))
    upperThresh = norm.ppf(eval(prodThresholds[1]))
    print(lowerThresh, upperThresh)
    newVuln = []


    """calculate vulnerabilities based on projections to risk axis and assign categories low, medium, high"""
    for i in range(len(SEp)):
        SEps.append(np.sqrt(2) * SEp[i] * SEs[i] / np.sqrt(SEp[i] ** 2 + SEs[i] ** 2))
        riskVector.append(np.matrix([SEs[i], SEp[i]]).transpose())
        projectionMatrix.append(riskVector[i] * (riskVector[i].transpose() * riskVector[i]) ** -1 * riskVector[i].transpose())
        transformation.append(mean - projectionMatrix[i] * mean)

        projection.append(projectionMatrix[i] * np.matrix([choiceProd[i], choiceSus[i]]).transpose())
        projection_p.append(projection[i][0, 0] + transformation[i][0, 0])
        projection_s.append(projection[i][1, 0] + transformation[i][1, 0])
        diffp = projection_p[i] - mean[0, 0]
        diffs = projection_s[i] - mean[1, 0]
        distanceMetric.append(np.sign(diffp) * np.sqrt((diffp) ** 2 + (diffs) ** 2))
        newVuln.append(norm.cdf(distanceMetric[i]/SEps[i]))
        if distanceMetric[i] < lowerThresh * SEps[i]:
            low += 1
            revisedCategory.append('low')
            print('low', distanceMetric[i], lowerThresh*SEps[i], 'Vs = '+ str(newVuln[i]))
        elif lowerThresh * SEps[i] <= distanceMetric[i] < upperThresh * SEps[i]:
            medium += 1
            revisedCategory.append('medium')
            print('medium', distanceMetric[i], lowerThresh*SEps[i], upperThresh * SEps[i], 'Vs = '+ str(newVuln[i]))
        else:
            high += 1
            revisedCategory.append('high')
            print('high', distanceMetric[i], upperThresh * SEps[i], 'Vs = '+ str(newVuln[i]))

    """calculate sPSA categories for comparison"""

    vulnerability = []
    oldCategory = []
    for g in range(len(productivity)):
        vulnerability.append(np.sqrt(productivity[g]**2+susceptibility[g]**2))
        if vulnerability[g] < 2.64:
            #print('low')
            oldCategory.append('low')
        elif 2.64 <= vulnerability[g] <= 3.18:
            #print('medium')
            oldCategory.append('medium')
        else:
            #print('high')
            oldCategory.append('high')

    """output figure with colors representing changes in category from sPSA to rPSA, size gives num overlapping points"""
    markerColor = []
    markertype = []
    lowerCat = 0
    higherCat = 0
    sameCat = 0
    print(oldCategory.count('low'), oldCategory.count('medium'), oldCategory.count('high'))
    print(revisedCategory.count('low'), revisedCategory.count('medium'), revisedCategory.count('high'))
    for w in range(len(revisedCategory)):
        if revisedCategory[w] == 'low':
            markerColor.append('b')
        elif revisedCategory[w] == 'medium':
            markerColor.append('y')
        else:
            markerColor.append('r')

    for v in range(len(oldCategory)):
        if oldCategory[v] == revisedCategory[v]:
            markertype.append("o")
            sameCat+=1
        elif (oldCategory[v] == 'low' and revisedCategory[v] == 'medium') or (oldCategory[v] == 'medium' and revisedCategory[v] == 'high'):
            markertype.append(6)
            higherCat+=1
        else:
            markertype.append(7)
            lowerCat+=1

    markerArea = []
    totalCat = lowerCat+higherCat+sameCat
    print(lowerCat/totalCat, higherCat/totalCat, sameCat/totalCat)
    print(lowerCat, higherCat, sameCat, totalCat)
    for i in range(len(oldCategory)):
        count = 0
        uniqueMeasure = [productivity[i], susceptibility[i]]
        for x in range(len(oldCategory)):
            if [productivity[x], susceptibility[x]] == uniqueMeasure:
                count += 1
        markerArea.append(200 * count)



    projPlot = plt.figure(figsize=(10,10))
    #plt.scatter(productivity, susceptibility, alpha=0.4, c=markerColor, s=markerArea, marker=markertype)
    mscatter(productivity, susceptibility, m=markertype, c=markerColor, s=markerArea, alpha=0.4)
    plt.axis((1, 3, 1, 3))
    plt.xlabel('Productivity')
    plt.ylabel('Susceptibility')
    circle1=plt.Circle((0, 0), 2.64, color='k', fill=None, linestyle='--', linewidth='2')
    circle2=plt.Circle((0, 0), 3.18, color='k', fill=None, linestyle='--', linewidth='2')
    currentAxis = plt.gca()
    currentAxis.add_artist(circle1)
    currentAxis.add_artist(circle2)
    plt.tight_layout()
    plt.xticks([1,2,3])
    plt.yticks([1,2,3])
    plt.show()
    projPlot.savefig("figure.png")

    vulnerability = np.array(vulnerability)/(3*np.sqrt(2))
    a = rankdata(vulnerability)
    b = rankdata(newVuln)
    print(sum(i != j for i, j in zip(a, b)))


main()