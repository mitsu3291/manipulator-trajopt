import numpy as np
import matplotlib.pyplot as plt
import time

# set parameters
dimx = 2         # number of state
dimu = 1         # number of input
tf = 20          # simulation time
dt = 1e-2        # time step
N = int(tf/dt)   # number of time step
eps = 1e-2       # terminal condition (stop if error is smaller than eps)
loop= 100        # maximun number of input updates
dalpha = 1e-4    # number of alpha step
lsloop = 100     # maximun number of alpha updates (linear search loop)

# set problem-specific parameters
beta1 = 3        # increase rate of x1
beta2 = 1        # increase rate of x2
gamma = 0.1      # x1-u1 coefficient
rho = 2          # x1-x2 coefficient (food chain effect)
delta = 0.1      # regulate food chain effect
kappa1 = 100     # capacity of x1
kappa2 = 1       # capacity of x2
umin = 0         # minimun value of u

# set weightings
sf1 = 1          # terminal x1
sf2 = 1          # terminal x2
q1 = 1           # stage x1
q2 = 1           # stage x2
r1 = 10          # stage u1
w1 = 10000       # constraint umin

# set variables
t = np.linspace(0, tf, N+1)     # time
x = np.zeros((N+1, dimx))       # state
lmd = np.zeros((N+1, dimx))     # costate
u = np.zeros((N, dimu))         # input
s = np.zeros((N, dimu))         # decent direction
x_tmp = np.zeros((N+1, dimx))   # state in linear search
u_tmp = np.zeros((N, dimu))     # input in linear search
x_es = np.zeros((N+1, dimx))    # state without input (only ecosystem)
printcode1 = 0                  # print "error did not converge"
printcode2 = 0                  # print "dalpha is too large"
printcode3 = 0                  # print "dalpha or lsloop is too small"
err = []                        # history of error

# initial conditon and initial guess
x0 = [100., 10.]
u0 = np.zeros((N, dimu))

# functions
# 1. equation of state: dx/dt = f
# f: dynamics
def f(x, u, t):
    f1 = x[0]*(beta1*(1 - x[0]/kappa1) - rho*x[1]/(x[0]+delta)) + gamma*u[0]
    f2 = x[1]*beta2*(1 - kappa2*x[1]/x[0])
    return np.array([f1, f2])

# 2. constraint
# C1: u>umin
def C1(x, u, t):
    if u[0] < umin:
        return 0.5*w1*(u[0] - umin)**2
    else:
        return 0.

# 3. cost function: J = phi + int L
# phi: terminal cost
def phi(x):
    phi = -0.5*(sf1*x[0]**2 + sf2*x[1]**2)
    return phi
# L: stage cost
def L(x, u, t):
    L = 0.5*(-q1*x[0]**2 - q2*x[1]**2 + r1*u[0]**2)
    return L
# J: cost function
def J(x, u, t):
    J = phi(x[N])
    for i in range(N):
        J += (L(x[i], u[i], t[i]) + C1(x[i], u[i], t[i]))*dt
    return J

# 4. derivatives
# Hamiltonian: H = L + lmd^{T} f
# dphi/dx
def dphidx(x):
    dphidx1 = -sf1*x[0]
    dphidx2 = -sf2*x[1]
    return np.array([dphidx1, dphidx2])
# dH/dx
def dHdx(x, u, lmd, t):
    dHdx1 = -q1*x[0] + lmd[0]*(beta1 - 2*beta1*x[0]/kappa1 - rho*delta*x[1]/(x[0]+delta)**2) + lmd[1]*beta2*kappa2*x[1]**2/x[0]**2
    dHdx2 = -q2*x[1] + lmd[0]*(-rho*x[0]/(x[0]+delta)) + lmd[1]*(beta2 - 2*beta2*kappa2*x[1]/x[0])
    return np.array([dHdx1, dHdx2])
# dH/du
def dHdu(x, u, lmd, t):
    if u[0] < umin:
        dHdu1 = r1*u[0] + lmd[0]*gamma + w1*(u[0] - umin)
    else:
        dHdu1 = r1*u[0] + lmd[0]*gamma
    return np.array([dHdu1])

# calculate x without input
# 1. initial state x[0] is given
x_es[0] = x0
# 2. calculate next x by equation of state (forward Euler method)
for i in range(N):
    x_es[i+1] = x_es[i] + f(x_es[i], u[i], t[i])*dt

# start clock
t0 = time.time()

# main
u = u0
for i in range(loop):
    # 1. calculate x
    # 1-1. initial state x[0] is given
    x[0] = x0
    # 1-2. calculate next x by equation of state (forward Euler method)
    for j in range(N):
        x[j+1] = x[j] + f(x[j], u[j], t[j])*dt

    # 2. calculate lmd
    # 2-1. terminal condition
    lmd[N] = dphidx(x[N])
    # 2-2. calculate previous lmd by adjoint equation (backward Euler method)
    for j in reversed(range(N)):
        lmd[j] = lmd[j+1] - dHdx(x[j], u[j], lmd[j+1], t[j])*(-dt)

    # 3. evaluate error
    # 3-1. calculate s=dH/du
    for j in range(N):
        s[j] = - dHdu(x[j], u[j], lmd[j+1], t[j])
    # 3-2. error is sqrt(sum((dH/du)**2))
    err_tmp = np.sqrt(np.sum(s**2))
    err.append(err_tmp)
    # 3-3. evaluate error (decide or update input)
    if err_tmp < eps:
        loop = i+1
        break
    if i == loop-1:
        printcode1 = 1

    # 4. update input
    #4-1. initial alpha and J
    alpha = 0.
    J_ = J(x, u, t)
    # 4-2. linear search of stepsiza
    for j in range(lsloop):
        # 4-2-1. update alpha
        alpha_tmp = alpha + dalpha
        # 4-2-2. update input
        for k in range(N):
            u_tmp[k] = u[k] + alpha_tmp*s[k]
        # 4-2-3. update state
        x_tmp[0] = x0
        for k in range(N):
            x_tmp[k+1] = x_tmp[k] + f(x_tmp[k], u_tmp[k], t[k])*dt
        # 4-2-4. update J
        J_tmp = J(x_tmp, u_tmp, t)
        # 4-2-5. evaluate J (decide or update alpha)
        if J_tmp > J_:
            if j == 0:
                printcode2 = 1
            break
        if j == lsloop-1:
            printcode3 = 1
            alpha = alpha_tmp
        J_ = J_tmp
        alpha = alpha_tmp
    # 4-3. update input
    for j in range(N):
        u[j] += alpha*s[j]

# stop clock and print time for calculate optimal input
tf = time.time()
print ("time:%f" % (tf-t0))

# print converge error
if printcode1 == 1:
    print("error did not converge")
if printcode2 == 1:
    print("dalpha is too large")
if printcode3 == 1:
    print("dalpha or lsloop is too small")

# plot
# fig1 state
fig1 = plt.figure()
plt.plot(t, x[:,0], label = "x1(t)")
plt.plot(t, x[:,1], label = "x2(t)")
plt.plot(t, x_es[:,0], linestyle = "dashed", label = "x1(t) (without input)")
plt.plot(t, x_es[:,1], linestyle = "dashed", label = "x2(t) (without input)")
plt.xlabel("time")
plt.ylabel("state")
plt.legend(loc = "upper right")
# fig2 input
fig2 = plt.figure()
plt.axhline(y = umin, color = "red")
plt.plot(t[:-1], u[:,0], label = "u(t)")
plt.xlabel("time")
plt.ylabel("input")
plt.legend(loc = "upper right")
# fig3 error
fig3 = plt.figure()
plt.plot(np.arange(0, loop, 1), err, label = "||F||")
plt.xlabel("iteration")
plt.ylabel("error")
plt.legend(loc = "upper right")
plt.show()