"""Code for simulating NOAA study PSA
Number of Productivity Attributes: 10
Model for Productivity: additive
Number of Susceptibility Attributes: 12
Model for Susceptibility: multiplicative
Weighting: attributes weighted in csv file
Other notes:
"""

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
import random
import pandas as pd
import numpy as np
import csv
import seaborn as sns; sns.set(style="white", color_codes=True)

SMALL_SIZE = 36
MEDIUM_SIZE = 48
BIGGER_SIZE = 72

def geo_mean(iterable):
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
def main():

    abs_file_path = "C:/Users/Richard/Desktop/PSAReanalysis/NOAA_UnitedStates/NOAA_Susceptibility.csv"
    abs_file_path2 = "C:/Users/Richard/Desktop/PSAReanalysis/NOAA_UnitedStates/NOAA_Productivity.csv"

    with open(abs_file_path, newline='') as csvfile:  # read in mosquito data and store in totalData
        data = list(csv.reader(csvfile))

    with open(abs_file_path2, newline='') as csvfile:  # read in mosquito data and store in totalData
        data2 = list(csv.reader(csvfile))

    data = data[1:]  # exclude first row (data header)
    data2 = data2[1:]

    newdata = []
    newdata2 =[]

    for x in range(len(data)):
        newdata.append(data[x][1:])
        newdata2.append(data2[x][1:])




    numProd = 10
    numSusc = 12

    susceptibility = []
    logsusceptibility = []
    productivity = []
    susceptError = []
    producError = []

    for i in range(len(newdata)):
        presusceptibility = []
        susceptweight = []
        for j in range(1, 2 + numSusc):
            if newdata[i][j] != '':
                susceptweight.append(float(newdata[i][j + numSusc + 2]))
                for b in range(int(newdata[i][j + numSusc + 2])):
                    presusceptibility.append(float(newdata[i][j]))
        if presusceptibility != []:
            suscep = geo_mean(presusceptibility)
            logsuscep = np.log(suscep)
            susceptibility.append(suscep)
            logsusceptibility.append(logsuscep)
            totalWeight = np.sum(susceptweight)
            for y in range(len(susceptweight)):
                susceptweight[y] = (susceptweight[y] / totalWeight) ** 2
            susceptError.append((.2058 * np.sum(susceptweight)) ** (1 / 2))



    for i in range(len(newdata2)):
        preproductivity=[]
        producweight=[]
        for j in range(1, 2+numProd):
            if newdata2[i][j] != '':
                producweight.append(float(newdata2[i][j + numProd + 2]))
                for b in range(int(newdata2[i][j+numProd+2])):
                    preproductivity.append(4-float(newdata2[i][j]))
        if preproductivity != []:
            prod = np.average(preproductivity)
            productivity.append(prod)
            totalWeight = np.sum(producweight)

            for z in range(len(producweight)):
                producweight[z] = (producweight[z]/totalWeight)**2
            producError.append((.6667*np.sum(producweight))**(1/2))


    SEp = producError
    SEs = susceptError
    SEps = []
    riskVector = []
    projectionMatrix =[]
    mean = np.matrix([2, np.log(2)]).transpose()
    transformation=[]
    projection = []
    projection_p = []
    projection_s = []
    distanceMetric = []
    revisedCategory = []
    low = 0
    medium = 0
    high = 0


    for i in range(len(producError)):
        SEps.append(np.sqrt(2) * SEp[i] * SEs[i] / np.sqrt(SEp[i] ** 2 + SEs[i] ** 2))
        riskVector.append(np.matrix([SEs[i], SEp[i]]).transpose())
        projectionMatrix.append(riskVector[i] * (riskVector[i].transpose() * riskVector[i]) ** -1 * riskVector[i].transpose())
        transformation.append(mean - projectionMatrix[i] * mean)

        projection.append(projectionMatrix[i] * np.matrix([productivity[i], logsusceptibility[i]]).transpose())
        projection_p.append(projection[i][0, 0] + transformation[i][0, 0])
        projection_s.append(projection[i][1, 0] + transformation[i][1, 0])
        diffp = projection_p[i] - mean[0, 0]
        diffs = projection_s[i] - mean[1, 0]
        distanceMetric.append(np.sign(diffp) * np.sqrt((diffp) ** 2 + (diffs) ** 2))
        if distanceMetric[i] < -.431 * SEps[i]:
            low += 1
            revisedCategory.append('low')
            #print('low')
        elif -.431 * SEps[i] <= distanceMetric[i] < .431 * SEps[i]:
            medium += 1
            revisedCategory.append('medium')
            #print('medium')
        else:
            high += 1
            revisedCategory.append('high')
            #print('high')


    vulnerability = []
    oldCategory = []
    for g in range(len(productivity)):
        vulnerability.append(np.sqrt(productivity[g]**2+susceptibility[g]**2))
        if vulnerability[g] < 2.64:
            print('low')
            oldCategory.append('low')
        elif 2.64 <= vulnerability[g] <= 3.18:
            print('medium')
            oldCategory.append('medium')
        else:
            print('high')
            oldCategory.append('high')

    markerColor = []
    markerStyle = []
    lowerCat = 0
    higherCat = 0
    sameCat = 0
    print(oldCategory.count('low'), oldCategory.count('medium'), oldCategory.count('high'))
    print(revisedCategory.count('low'), revisedCategory.count('medium'), revisedCategory.count('high'))
    for v in range(len(oldCategory)):
        if oldCategory[v] == revisedCategory[v]:
            markerColor.append('k')
            markerStyle.append('o')
            sameCat+=1
        elif (oldCategory[v] == 'low' and revisedCategory[v] == 'medium') or (oldCategory[v] == 'medium' and revisedCategory[v] == 'high'):
            markerColor.append('r')
            markerStyle.append('^')
            higherCat+=1
        else:
            markerColor.append('b')
            markerStyle.append('v')
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
    plt.scatter(productivity, susceptibility, alpha=0.4, c=markerColor, s=markerArea)
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
    projPlot.savefig("NOAAFigure1.tiff")



main()