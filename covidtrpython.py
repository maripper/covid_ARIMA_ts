# %%

#Abaixo definimos os parâmetros o valor é o "value":

delta = 0.4
eta = 0.02
mi = 0.0000040849
alfa = 0.02
theta = 0.03
beta = 4.3e-6
#  beta variável? beta_w = [0.00001503443,0.00002822209,0.00000207524]
gamma = 0.97
omega = 0.0000058
alfa_2 = [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
       40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
       40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
       40, 40, 40, 40, 40, 40, 40, 40, 40, 70, 70, 70, 70, 70, 70, 70, 70,
       70, 70, 70, 70, 70, 70, 70, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45]
beta_c = 4.3e-7
len(alfa_2)

# %%
#Definimos os valores incicais das cidades 1: Araraquara, 2: Ribeirão e 3: Sertãozinho
N_C1 = 538339 
N_C2 = 720116 
N_C3 = 428432 

#Abaixo os valores iniciais de cada compartimento C: Confinados, I: Infectados, R: Recuperados, D: Mortos e S: Sucetíveis
C0_C1 = 0
I0_C1 = 1
R0_C1 = 0
D0_C1 = 0
S0_C1 = N_C1 - (C0_C1+I0_C1+R0_C1+D0_C1)
C0_C2 = 0
I0_C2 = 1
R0_C2 = 0
D0_C2 = 0
S0_C2 = N_C2 - (C0_C2+I0_C2+R0_C2+D0_C2)
C0_C3 = 0
I0_C3 = 1
R0_C3 = 0
D0_C3 = 0
S0_C3 = N_C3 - (C0_C3+I0_C3+R0_C3+D0_C3)
#Colocamos os valores iniciais em uma array:
y0 = [S0_C1,C0_C1,I0_C1,R0_C1,D0_C1,S0_C2,C0_C2,I0_C2,R0_C2,D0_C2,S0_C3,C0_C3,I0_C3,R0_C3,D0_C3]
f1 = 3.847
f2 = 3.829
f3 = 5.215
f4 = 5.343
f5 = 4.550
f6 = 4.500

# %%
#Criamos uma função com o modelo:
def deriv(y,t,beta, gamma, delta, eta, mi, alfa, theta, alfa_2,beta_c,f1,f2,f3,f4,f5,f6):
    S_C1,C_C1,I_C1,R_C1,D_C1,S_C2,C_C2,I_C2,R_C2,D_C2,S_C3,C_C3,I_C3,R_C3,D_C3= y
    dSdt_C1 = mi * (S_C1+C_C1+I_C1+R_C1) - alfa_2 * S_C1 + delta * C_C1 - beta * S_C1 * I_C1 + eta * R_C1 - mi * S_C1 - f1*beta*I_C2*S_C1 - f6*beta*I_C3*S_C1
    dCdt_C1 = alfa_2 * S_C1 - delta * C_C1 - mi * C_C1 - beta_c * C_C1 * I_C1
    dIdt_C1 = beta * S_C1 * I_C1 - theta * I_C1 - gamma * I_C1 - mi * I_C1 + beta_c * C_C1 * I_C1 + f2*beta*I_C1* S_C2 + f5*beta*I_C1* S_C3
    dRdt_C1 = gamma * I_C1 - eta * R_C1 - mi * R_C1
    dDdt_C1 = theta * I_C1
    dSdt_C2 = mi * (S_C2+C_C2+I_C2+R_C2) - alfa * S_C2 + delta * C_C2 - beta * S_C2 * I_C2 + eta * R_C2 - mi * S_C2 - f2*beta*I_C1* S_C2 - f3*beta*I_C3* S_C2
    dCdt_C2 = alfa * S_C2 - delta * C_C2 - mi * C_C2 - beta_c * C_C2 * I_C2
    dIdt_C2 = beta * S_C2 * I_C2 - theta * I_C2 - gamma * I_C2 - mi * I_C2 + beta_c * C_C2 * I_C2 +  f1*beta*I_C2*S_C1  + f4*beta*I_C2*S_C3
    dRdt_C2 = gamma * I_C2 - eta * R_C2 - mi * R_C2
    dDdt_C2 = theta * I_C2
    dSdt_C3 = mi * (S_C3+C_C3+I_C3+R_C3) - alfa * S_C3 + delta * C_C3 - beta * S_C3 * I_C3 + eta * R_C3 - mi * S_C3 - f4*beta*I_C2*S_C3 - f5*beta*I_C1* S_C3
    dCdt_C3 = alfa * S_C3 - delta * C_C3 - mi * C_C3 - beta_c * C_C3 * I_C3
    dIdt_C3 = beta * S_C3 * I_C3 - theta * I_C3 - gamma * I_C3 - mi * I_C3 + beta_c * C_C3 * I_C3 + f3*beta*I_C3* S_C2 + f6*beta*I_C3*S_C1
    dRdt_C3 = gamma * I_C3 - eta * R_C3 - mi * R_C3
    dDdt_C3 = theta * I_C3
    return [dSdt_C1,dCdt_C1,dIdt_C1,dRdt_C1,dDdt_C1,dSdt_C2,dCdt_C2,dIdt_C2,dRdt_C2,dDdt_C2,dSdt_C3,dCdt_C3,dIdt_C3,dRdt_C3,dDdt_C3]
    

# %%
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Integrate the equations over the time grid, t.
tspan = np.linspace(0,60,60)
#Integramos a função no tempo para confinamento de Araraquara de 35% em 60 dias
sol1 = odeint(deriv, y0, tspan, args=( beta, gamma, delta, eta, mi, alfa, theta, 0.35,beta_c, f1,f2,f3,f4,f5,f6))

tspan2 = np.linspace(60,75,15)
tspan21 = np.linspace(0,15,15)

y02 = [sol1.T[0][59],sol1.T[1][59],sol1.T[2][59],sol1.T[3][59],sol1.T[4][59],sol1.T[5][59],sol1.T[6][59],sol1.T[7][59],
        sol1.T[8][59],sol1.T[9][59],sol1.T[10][59],sol1.T[11][59],sol1.T[12][59],sol1.T[13][59],sol1.T[14][59]]
#Integramos a função no tempo para confinamento de Araraquara de 65% em 15 dias
sol2 = odeint(deriv, y02, tspan21, args=( beta, gamma, 0, eta, mi, alfa, theta, 0.60,beta_c,0,0,f3,f4,0,0))

tspan3 = np.linspace(75,100,25)
tspan31 = np.linspace(0,25,25)

y03 = [sol2.T[0][14],sol2.T[1][14],sol2.T[2][14],sol2.T[3][14],sol2.T[4][14],sol2.T[5][14],sol2.T[6][14],sol2.T[7][14],
        sol2.T[8][14],sol2.T[9][14],sol2.T[10][14],sol2.T[11][14],sol2.T[12][14],sol2.T[13][14],sol2.T[14][14]]
#Integramos a função no tempo para confinamento de Araraquara de 40% em 25 dias
sol3 = odeint(deriv, y03, tspan31, args=( beta, gamma, delta, eta, mi, alfa, theta, 0.4,beta_c,f1,f2,f3,f4,f5,f6))

#Junção dos tempos
tspan = tspan.tolist() + tspan2.tolist() + tspan3.tolist()

# %%
#Junção dos dados em cada período de confinamento:
dSdt_C1 = sol1.T[0].tolist() + sol2.T[0].tolist() + sol3.T[0].tolist()
dCdt_C1 = sol1.T[1].tolist() + sol2.T[1].tolist() + sol3.T[1].tolist()
dIdt_C1 = sol1.T[2].tolist() + sol2.T[2].tolist() + sol3.T[2].tolist()
dRdt_C1 = sol1.T[3].tolist() + sol2.T[3].tolist() + sol3.T[3].tolist()
dDdt_C1 = sol1.T[4].tolist() + sol2.T[4].tolist() + sol3.T[4].tolist()
dSdt_C2 = sol1.T[5].tolist() + sol2.T[5].tolist() + sol3.T[5].tolist()
dCdt_C2 = sol1.T[6].tolist() + sol2.T[6].tolist() + sol3.T[6].tolist()
dIdt_C2 = sol1.T[7].tolist() + sol2.T[7].tolist() + sol3.T[7].tolist()
dRdt_C2 = sol1.T[8].tolist() + sol2.T[8].tolist() + sol3.T[8].tolist()
dDdt_C2 = sol1.T[9].tolist() + sol2.T[9].tolist() + sol3.T[9].tolist()
dSdt_C3 = sol1.T[10].tolist() + sol2.T[10].tolist() + sol3.T[10].tolist()
dCdt_C3 = sol1.T[11].tolist() + sol2.T[11].tolist() + sol3.T[11].tolist()
dIdt_C3 = sol1.T[12].tolist() + sol2.T[12].tolist() + sol3.T[12].tolist()
dRdt_C3 = sol1.T[13].tolist() + sol2.T[13].tolist() + sol3.T[13].tolist()
dDdt_C3 = sol1.T[14].tolist() + sol2.T[14].tolist() + sol3.T[14].tolist()

# %%
#Plotamos todos os compartimentos:
plt.rcParams["figure.figsize"] = (20,10)
plt.subplot(1, 2, 1)
plt.plot(tspan, dSdt_C1, 'b', label='Sucetíveis Cidade1')
plt.plot(tspan, dSdt_C2, 'g', label='Sucetíveis Cidade2')
plt.plot(tspan, dSdt_C3, 'r', label='Sucetíveis Cidade3')

plt.legend(loc='best')
plt.xlabel('t')
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(tspan, dCdt_C1 , 'b', label='Confinados Cidade1')
plt.plot(tspan, dCdt_C2 , 'g', label='Confinados Cidade2')
plt.plot(tspan, dCdt_C3 , 'r', label='Confinados Cidade3')

plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()


# %%
plt.rcParams["figure.figsize"] = (20,10)
plt.subplot(1, 2, 1)
plt.plot(tspan, dIdt_C1 , 'b', label='Infectados Cidade1')
plt.plot(tspan, dIdt_C2 , 'g', label='Infectados Cidade2')
plt.plot(tspan, dIdt_C3 , 'r', label='Infectados Cidade3')

plt.legend(loc='best')
plt.xlabel('t')
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(tspan, dRdt_C1 , 'b', label='Recuperados Cidade1')
plt.plot(tspan, dRdt_C2 , 'g', label='Recuperados Cidade2')
plt.plot(tspan, dRdt_C3 , 'r', label='Recuperados Cidade3')

plt.legend(loc='best')
plt.xlabel('t')
plt.grid()

plt.show()

# %%

plt.plot(tspan, dDdt_C1 , 'b', label='Dead Cidade1')
plt.plot(tspan, dDdt_C2 , 'g', label='Dead Cidade2')
plt.plot(tspan, dDdt_C3 , 'r', label='Dead Cidade3')

plt.legend(loc='best')
plt.xlabel('t')
plt.grid()

# %%


# %%



