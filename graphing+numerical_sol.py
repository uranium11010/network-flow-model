import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy import optimize as opt
import xlrd
import warnings

warnings.filterwarnings("ignore")

#order species in decreasing order of degree to make the topology look as nested as possible (used to assess degree of nestedness and detect complete topological nestedness)
def makeNested(Fns):
    rowNonzero = np.count_nonzero(Fns, axis=1)
    rNZ = []
    for i in range(SH):
        rNZ.append((rowNonzero[i], i))
    colNonzero = np.count_nonzero(Fns, axis=0)
    cNZ = []
    for j in range(SP):
        cNZ.append((colNonzero[j], j))
    rNZ.sort(key=lambda tup: tup[0], reverse=True)
    cNZ.sort(key=lambda tup: tup[0], reverse=True)

    F = np.empty(Fns.shape)
    for i in range(Fns.shape[0]):
        for j in range(Fns.shape[1]):
            F[i,j] = Fns[rNZ[i][1],cNZ[j][1]]

    return F

#determine whether the topology is completely nested
def nested(A): #A must be an output from makeNested(Fns)
    B = A.astype(bool)
    return bool(np.sum(np.append(np.full((1, B.shape[1]), True), B, 0) != np.append(B, np.full((1, B.shape[1]), False), 0)) == B.shape[1])

#generate N random networks with the same topology and aggregate flow rates of Fin
def generateRandom(Fin, dx, ey, SH, SP, fact, N):
    #multiply by fact to make dx, ey integers
    dx = np.rint(dx * fact).astype(int)
    ey = np.rint(ey * fact).astype(int)

    #generate stubs on resources and consumers (for later generation of multigraph)
    Hstubs = []
    for x in range(SH):
        for i in range(dx[x]):
            Hstubs.append(x)
    Pstubs = []
    for y in range(SP):
        for i in range(ey[y]):
            Pstubs.append(y)

    #create random networks
    Farrays = []
    for k in range(N):
        print("random network #" + str(k))
        length = min(len(Hstubs), len(Pstubs))

        #create a multigraph by joining stubs
        success = False
        while not success:
            F = np.zeros((SH, SP))
            copy = [Pstubs[i] for i in range(length)]
            for s in range(length-1,-1,-1):
                good = False
                for t in range(s+1):
                    if Fin[Hstubs[s],copy[t]]:
                        good = True
                        break
                if not good:
                    break #no available stub on a predator to connect to prey; we're stuck so restart process of connecting stubs
                t = np.random.randint(s+1)
                
                #find a random available stub on a predator to connect to prey
                while not Fin[Hstubs[s],copy[t]]:
                    t = np.random.randint(s+1)
                F[Hstubs[s],copy[t]] += 1
                copy.pop(t)
            
            success = s == 0
        
        F /= fact #flow network from multigraph by dividing by fact
        Farrays.append(F)

    return Farrays


#IMPORTING DATA

#Glen Canyon (Cross et al. 2011)
Fglen = np.array([[1.558,	2.849,	1.388,	0.074,	44.595,	28.359,	0.818,	13.466],
[3.385,	6.191,	3.716,	0.316,	70.47,	33.792,	3.048,	22.246],
[0.278,	0.508,	0.297,	0.004,	0.918,	2.409,	0.273,	5.849],
[0.198,	0.362,	0.322,	0.002,	3.387,	3.661,	0.344,	0],
[0.03,	0.056,	0.049,	0.002,	1.158,	0.719,	0,	0],
[0.017,	0.032,	0.106,	0,	0,	0,	0,	0]])
factglen = 1000 #factor to multiply to generate the multigraph when applying the modified configuration model to generate a random network

#Forest streams (Bumpers et al. 2017); values are *per gut* values
sheet = xlrd.open_workbook("potential energy flow data/forest streams energy flows.xlsx").sheet_by_index(1)
shape = (sheet.nrows-1, sheet.ncols-1)
Fforest = np.empty(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        Fforest[i,j] = sheet.cell(i+1, j+1).value
i = 0
while i < Fforest.shape[0]:
    if np.count_nonzero(Fforest[i]) <= 1:
        Fforest = np.delete(Fforest, i, 0) #remove resources with only 1 consumer (since their flow is predicted with 100% accuracy and would thus overestimate the accuracy of our model)
    else:
        i += 1
factforest = 10

#Coral reef (Cocheret de la Morinière et al. 2003) bay; values are *percentages*
sheet = xlrd.open_workbook("potential energy flow data/coral reef energy flows.xlsx").sheet_by_index(0)
shape = (sheet.nrows-1, sheet.ncols-1)
Fcoralbay = np.empty(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        Fcoralbay[i,j] = sheet.cell(i+1, j+1).value
Fcoralbay = np.transpose(Fcoralbay)
factcoralbay = 1

#Coral reef (Cocheret de la Morinière et al. 2003) reef; values are *percentages*
sheet = xlrd.open_workbook("potential energy flow data/coral reef energy flows.xlsx").sheet_by_index(1)
shape = (sheet.nrows-1, sheet.ncols-1)
Fcoralreef = np.empty(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        Fcoralreef[i,j] = sheet.cell(i+1, j+1).value
j = 0
while j < Fcoralreef.shape[1]:
    if np.count_nonzero(Fcoralreef[:,j]) == 1:
        Fcoralreef = np.delete(Fcoralreef, j, 1) #remove resources with only 1 consumer
    else:
        j += 1
Fcoralreef = np.transpose(Fcoralreef)
factcoralreef = 1

#Crustacean (Rudnick & Resh 2005); values are *percentage abundances*
sheet = xlrd.open_workbook("potential energy flow data/crustacean energy flows.xlsx").sheet_by_index(0)
shape = (sheet.nrows-1, sheet.ncols-1)
Fcrust = np.empty(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        Fcrust[i,j] = sheet.cell(i+1, j+1).value
factcrust = 10

#Earthworm (Judas et al. 1992); values are median % dry weight in all guts of a consumer species
Fworm = np.array([[24, 39, 61, 60, 56, 71],
[35, 58, 37, 38, 44, 29],
[38, 0, 0, 0, 0, 0]])
Fworm = Fworm[:2]
factworm = 1

#Chesapeake
Fchesa = np.array([[2.7, 12.3, 17.2],
[2.6, 8, 10.6]])
factchesa = 10

#Caribbean
FcarA = np.array([[0.08, 4.8],
[1.13, 6.3]])
factcarA = 100

FcarB = np.array([[0.224, 0.28, 0.48, 0.18],
[2.02, 2.54, 4.34, 1.54],
[0.269, 0, 0.58, 0.20]])
factcarB = 1000

FcarC = np.array([[6.42, 102.85, 10.28, 40],
[8.57, 77.14, 7.71, 30]])
factcarC = 100

FcarD = np.array([[1200, 900],
[1200, 900]])
factcarD = 1/300

FcarE = np.array([[15, 610],
[15, 610]])
factcarE = 0.2

FseaA = np.array([[10, 0, 1, 18],
[37.5, 37.5, 48, 25],
[40, 50, 49, 50]])
factseaA = 2

FseaB = np.array([[1, 9, 33, 40],
[2, 4, 30, 27],
[1, 2, 17, 12.5],
[1, 1, 3, 2],
[7, 0, 0, 2],
[7, 16, 0, 0]])
factseaB = 2

#IMPORTING DATA ENDS


Flows = [Fglen, Fforest, Fcoralbay, Fcoralreef, Fcrust, Fworm, Fchesa, FcarA, FcarB, FcarC, FcarD, FcarE, FseaA, FseaB]
factors = [factglen, factforest, factcoralbay, factcoralreef, factcrust, factworm, factchesa, factcarA, factcarB, factcarC, factcarD, factcarE, factseaA, factseaB]
graphing = [False, False, False, False, False, False, False, False, False, False, False, False, False, False] #whether to graph log10(data/model)
recordDataModel = False #record values of log10(data/model)
recordGlenRandomModel = False #record random model output for Diatom-Simuliid pair for Glen Canyon food web

print("# networks:", len(Flows))

if recordDataModel:
    dataRes = np.empty((0,))

for k in range(len(Flows)):
    Fns = Flows[k]
    SH, SP = Fns.shape
    F = makeNested(Fns) #order species in decreasing degree
    print("# resources:", F.shape[0], "\t# consumers:", F.shape[1], "\t# linkages:", np.count_nonzero(F), "\tnested?", nested(F))
    
    #calculate model flow rates by numerically solving Eq. 6ab of Main Text
    L = F != 0
    dx = np.sum(F, axis=1)
    ey = np.sum(F, axis=0)

    def eqs(coef): #[a2, a3, ..., b1, b2, ...]
        D = coef[:(SH-1)] * np.dot(L, coef[(SH-1):])[1:] - dx[1:]
        E = coef[(SH-1):] * np.dot(np.transpose(L), np.append(1, coef[:(SH-1)])) - ey
        return np.append(D, E)

    sol = opt.fsolve(eqs, np.append(dx[1:], ey))

    A = np.append(1, sol[:(SH-1)])
    B = sol[(SH-1):]

    F0 = np.dot(np.transpose(np.array([A,])), np.array([B,])) * L #MODEL FLOW RATES

    #random model for glen canyon diatoms --> simuliids
    if k == 0 and recordGlenRandomModel:
        for i in range(SH):
            print("data/model:", i, F[i,3]/F0[i,3]) #simuliids
        print("data:", F[1,3], "model:", F0[1,3]) #diatoms --> simuliids

        diatSimulFlow = []
        for random in generateRandom(F, dx, ey, SH, SP, factors[k], 1000):
            diatSimulFlow.append(random[1,3])
        
        with open("glen_diat_simul_random.txt", "w") as fout:
            for f in diatSimulFlow:
                fout.write(str(f)+'\n')
    
    #record values of log10(data/model)
    if recordDataModel:
        dataDivide = (np.log10(F/F0)).flatten()
        dataRes = np.append(dataRes, dataDivide)

    #graph log10(data/model) of a single network in colored map format
    if graphing[k]:
        current_cmap = cm.get_cmap()
        current_cmap.set_bad(color="black")
        plt.axis("off")
        plt.imshow(np.log10(F/F0), origin="upper")
        plt.colorbar()
        plt.show()

#record log10(data/model)
if recordDataModel:
    with open("log10_data_model.txt", "w") as fout:
        for i in range(len(dataRes)):
            if not np.isnan(dataRes[i]):
                fout.write(str(dataRes[i])+'\n')
