import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as opt


def AB(S, sigma, tau):
    """
    get A(i) and B(j) from sigma and tau
    """
    A = np.empty((sigma[-1],), dtype=int)
    B = np.empty((tau[-1],), dtype=int)
    for k in range(S):
        A[sigma[k]:sigma[k+1]] = k+1
        B[tau[k]:tau[k+1]] = k+1

    return A, B


def get_E(S, sigma, tau):
    """
    get number of edges from sigma and tau
    """
    ans = 0
    for i in range(S):
        ans += sigma[S-i] * (tau[i+1] - tau[i])
    return ans


def get_random_L(R, C, E):
    """
    get random topology with given # edges
    """
    L = np.zeros((R,C), dtype=bool)
    while True:
        ind = np.random.choice(np.arange(R*C, dtype=int), (E,), False)
        for i in ind:
            L[i//C, i%C] = True
        if np.all(np.any(L, axis=1)) and np.all(np.any(L, axis=0)):
            break
    return L


def get_random_topology(R, C, get_L=True):
    """
    get random topology given R,C
    """
    S = np.random.choice(np.arange(2, min(R,C)+1, dtype=int))
    sigma = np.concatenate((np.concatenate(([0], np.sort(np.random.choice(np.arange(1,R), size=(S-1,), replace=False)))), [R]))
    tau = np.concatenate((np.concatenate(([0], np.sort(np.random.choice(np.arange(1,C), size=(S-1,), replace=False)))), [C]))
    E = get_E(S, sigma, tau)
    if get_L:
        L = get_random_L(R, C, E)
        return S, sigma, tau, E, L
    else:
        return S, sigma, tau, E


def max_ent_model(SH, SP, S, d, e, sigma, tau, L=None):
    """
    get MaxEnt flow rates
    """
    if not (L is None): #numerical solution
        def eqs(coef): #[a2, a3, ..., b1, b2, ...]
            D = coef[:(SH-1)] * np.dot(L, coef[(SH-1):])[1:] - d[1:]
            E = coef[(SH-1):] * np.dot(np.transpose(L), np.append(1, coef[:(SH-1)])) - e
            return np.append(D, E)

        sol = opt.fsolve(eqs, np.append(d[1:], e))

        A = np.append(1, sol[:(SH-1)])
        B = sol[(SH-1):]

        return np.dot(np.transpose(np.array([A,])), np.array([B,])) * L #MODEL FLOW RATES

    #analytic solution
    A, B = AB(S, sigma, tau)
    F = np.zeros((SH, SP))
    for x in range(SH):
        for y in range(tau[S+1-A[x]]):
            F[x,y] = d[x]*e[y] \
            * np.prod([np.sum(d[:sigma[S-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+1-A[x])]) \
            / np.prod([np.sum(d[:sigma[S+1-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+2-A[x])])

    return F


def get_random_F(R, C, S, sigma, tau, F_dist, e_d_dist, nested=True, max_ent=True, L=None):
    """
    generate random flow rate matrix
    """
    F = np.zeros((R, C))
    if max_ent:
        d = e_d_dist(size=(R,))*C
        e = e_d_dist(size=(C,))*R
        F = max_ent_model(R, C, S, d, e, sigma, tau, L)
        if np.any(F < 0):
            return None
    else:
        if nested:
            A, _ = AB(S, sigma, tau)
            linkages = 0
            for i in range(R):
                linkages += tau[S+1-A[i]]
            for i in range(R):
                for j in range(tau[S+1-A[i]]):
                    F[i,j] = F_dist() * R*C / linkages
        else:
            F = F_dist(size=(R,C))
    return F


def community(R, C, S, sigma, tau, F_dist, e_d_dist, eta_dist, dFdB_dist, g_dist, other_dist, sub_dist, nested=True, max_ent=True, L=None):
    """
    generate random community matrix
    """
    F = get_random_F(R, C, S, sigma, F_dist, e_d_dist, tau, nested, max_ent, L)
    if F is None:
        return
    dFdB = np.zeros((R, C))
    eta = eta_dist(size=(R, C))
    M = np.zeros((R+C, R+C))
    
    M[:R,:R] = other_dist(size=(R,R))
    dFdB = F * dFdB_dist(size=(R,C))
    M[:R,R:] -= F * sub_dist(size=(R,C))
    M[R:,:R] = (eta * dFdB * sub_dist(size=(R,C))).T
    for i in range(R):
        M[i,i] = -g_dist()
    for j in range(R,R+C):
        M[j,j] = -g_dist()/10

    return M, np.linalg.eigvals(M), F


def analytic_eigen(R, C, F_dist, e_d_dist, nested=True, max_ent=True):
    """
    Analytic approach: test z.T dM z on eigenvectors z of dM where dM is symmetric and M0 is anti-symmetric
    """
    F = None
    while F is None:
        S, sigma, tau, E, L = get_random_topology(R, C)
        F = get_random_F(R, C, S, sigma, tau, F_dist, e_d_dist, nested, max_ent, None if nested and max_ent else L)
    
    M0 = np.vstack((np.hstack((np.zeros((R, R)), -F)), np.hstack((F.T, np.zeros((C, C))))))
    # dM = np.vstack((np.hstack((np.zeros((R, R)), F)), np.hstack((F.T, np.zeros((C, C))))))
    dM = np.vstack((np.zeros((R, R+C)), np.hstack((F.T*np.random.uniform(size=(C, R)), np.zeros((C, C))))))

    eig_val, eig_vec = np.linalg.eig(M0)
    all_pos = True
    for i in range(R+C):
        v = eig_vec[:,i]
        print(v.conj().T)
        print(dM @ v)
        c = v.conj().T @ dM @ v
        print(c)
        all_pos = all_pos and c
    return all_pos


#probability distributions
def F_dist(size=None):
    return np.random.exponential(size=size)
def e_d_dist(size=None):
    return np.sqrt(np.random.chisquare(2, size=size) * (2/np.pi))
def eta_dist(size=None):
    return np.random.beta(2, 18, size=size)
def dFdB_dist(size=None):
    return np.random.uniform(size=size)
def g_dist(size=None):
    return np.random.exponential(size=size)
def other_dist(size=None):
    return np.zeros(size)
def sub_dist(size=None):
    return np.ones(size)


def generate_evaluate(min_S, max_R, max_C, K):
    """
    generate random community matrices and assess stability
    """
    max_l_rand_all = np.empty((max_R, max_C))
    max_l_rand_nest_all = np.empty((max_R, max_C))
    max_l_nest_all = np.empty((max_R, max_C))
    max_l_max_all = np.empty((max_R, max_C))
    for R in range(min_S, max_R):
        for C in range(min_S, max_C):
            print(R,C)
            max_l_rand = np.empty((K,))
            max_l_rand_nest = np.empty((K,))
            max_l_nest = np.empty((K,))
            max_l_max = np.empty((K,))

            k = 0
            while k < K:
                S, sigma, tau, E, L = get_random_topology(R,C)

                l_max = community(R, C, S, sigma, tau, F_dist, e_d_dist, eta_dist, dFdB_dist, g_dist, other_dist, sub_dist, True, True)
                l_nest = community(R, C, S, sigma, tau, F_dist, e_d_dist, eta_dist, dFdB_dist, g_dist, other_dist, sub_dist, False, True, L)
                if l_max != None and l_nest != None:
                    l_rand = community(R, C, S, sigma, tau, F_dist, e_d_dist, eta_dist, dFdB_dist, g_dist, other_dist, sub_dist, False, False)
                    l_rand_nest = community(R, C, S, sigma, tau, F_dist, e_d_dist, eta_dist, dFdB_dist, g_dist, other_dist, sub_dist, True, False)
                    max_l_rand[k] = np.max(np.real(l_rand[1]))
                    max_l_rand_nest[k] = np.max(np.real(l_rand_nest[1]))
                    max_l_nest[k] = np.max(np.real(l_nest[1]))
                    max_l_max[k] = np.max(np.real(l_max[1]))
                    
                    k += 1
                
            
            max_l_rand_all[R,C] = np.sum(max_l_rand < 0)/K
            max_l_rand_nest_all[R,C] = np.sum(max_l_rand_nest < 0)/K
            max_l_nest_all[R,C] = np.sum(max_l_nest < 0)/K
            max_l_max_all[R,C] = np.sum(max_l_max < 0)/K

    #record data and graph
    all_df = [pd.DataFrame(max_l_rand_all[min_S:,min_S:]), pd.DataFrame(max_l_rand_nest_all[min_S:,min_S:]), pd.DataFrame(max_l_nest_all[min_S:,min_S:]), pd.DataFrame(max_l_max_all[min_S:,min_S:])]
    writer = pd.ExcelWriter("2022-03-24 stabilities.xlsx", engine="xlsxwriter")
    for k in range(len(all_df)):
        all_df[k].to_excel(writer, sheet_name="Sheet"+str(k+1), header=False, index=False)
    writer.save()

    fig, ax = plt.subplots(2, 2)
    fig.tight_layout()

    ax[0,0].imshow(max_l_rand_all, vmin=0, vmax=1)
    ax[0,1].imshow(max_l_rand_nest_all, vmin=0, vmax=1)
    ax[1,0].imshow(max_l_nest_all, vmin=0, vmax=1)
    im = ax[1,1].imshow(max_l_max_all, vmin=0, vmax=1)

    ax[0,0].set_title("a")
    ax[0,1].set_title("b")
    ax[1,0].set_title("c")
    ax[1,1].set_title("d")

    for j in [(0,0), (0,1), (1,0), (1,1)]:
        ax[j].set_xlim([min_S-0.5, max_C-0.5])
        ax[j].set_ylim([max_R-0.5, min_S-0.5])
        ax[j].set_xticks(np.arange(min_S, max_C, 1))
        ax[j].set_yticks(np.arange(min_S, max_R, 1))

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.6])
    plt.colorbar(im, cax=cbar_ax)

    plt.show()

    fig.savefig("2022-03-24 stabilities.png")
    fig.savefig("2022-03-24 stabilities.pdf")


if __name__ == "__main__":
    # generate_evaluate(min_S = 2, max_R = 51, max_C = 51, K = 10)
    analytic_eigen(5, 5, F_dist, e_d_dist, False, False)
