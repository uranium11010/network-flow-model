import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import patches

Sarr = []
SHarr = []
SParr = []
sigmaarr = []
tauarr = []

Sarr.append(6)
SHarr.append(9)
SParr.append(13)
sigmaarr.append(np.array([0, 2, 5, 6, 7, 8, SHarr[-1]]))
tauarr.append(np.array([0, 1, 2, 5, 11, 12, 13, SParr[-1]]))

Sarr.append(11)
SHarr.append(20)
SParr.append(15)
sigmaarr.append(np.array([0, 6, 8, 9, 10, 13, 14, 16, 17, 18, 19, SHarr[-1]]))
tauarr.append(np.array([0, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, SParr[-1]]))

Sarr.append(5)
SHarr.append(20)
SParr.append(15)
sigmaarr.append(np.array([0, 13, 15, 18, 19, SHarr[-1]]))
tauarr.append(np.array([0, 9, 12, 13, 14, SParr[-1]]))

Sarr.append(10)
SHarr.append(20)
SParr.append(15)
sigmaarr.append(np.array([0, 6, 9, 10, 11, 12, 13, 17, 18, 19, SHarr[-1]]))
tauarr.append(np.array([0, 3, 4, 6, 7, 8, 9, 12, 13, 14, SParr[-1]]))

for bruh in range(4):
    S = Sarr[bruh]
    SH = SHarr[bruh]
    SP = SParr[bruh]
    sigma = sigmaarr[bruh]
    tau = tauarr[bruh]

    d = np.array([1/SH for k in range(SH)])
    e = np.array([1/SP for k in range(SP)])

    A = np.empty((SH,), dtype=int)
    B = np.empty((SP,), dtype=int)
    for k in range(S):
        A[sigma[k]:sigma[k+1]] = k+1
        B[tau[k]:tau[k+1]] = k+1

    F = np.zeros((SH+1, SP+1))
    for x in range(SH):
        for y in range(tau[S+1-A[x]]):
            F[x+1,y+1] = d[x]*e[y] \
            * np.prod([np.sum(d[:sigma[S-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+1-A[x])]) \
            / np.prod([np.sum(d[:sigma[S+1-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+2-A[x])])

    #current_cmap = cm.get_cmap()
    #current_cmap.set_bad(color="black")
    #plt.axis("off")
    fig, ax = plt.subplots()
    plt.imshow(F, origin="upper")
    #plt.clim(0, 0.02)
    plt.colorbar()
    ax.xaxis.set_label_position("top")
    plt.xlabel("j")
    plt.ylabel("i")
    ax.xaxis.tick_top()
    plt.xticks(np.arange(1, SP+1, 1))
    plt.yticks(np.arange(1, SH+1, 1))
    plt.xlim([0.5, SP+0.5])
    plt.ylim([SH+0.5, 0.5])
    plt.plot([0.5, SP+0.5], [SH+0.5, 0.5], c="C1")
    if bruh == 0:
        ax.add_patch(patches.Circle((2, 8), radius=1.7, fill=False, color="C3", lw=2))
    elif bruh == 3:
        ax.add_patch(patches.Circle((9, 13), radius=2.3, fill=False, color="C3", lw=2))
    plt.show()
