
# coding: utf-8

# In[5]:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# P1 Tarea 9

np.random.seed(100)


def Bootstrap(datos, alfa):
    '''Recibe un arreglo de datos y retorna un
    intervalo de confianza al 100*(1-alfa)%
    usando el metodo de bootstrap'''
    numdatos = len(datos)
    numiter = int(numdatos * np.log10(numdatos)**2)
    promedios = np.zeros(numiter)
    for i in range(numiter):
        indice = np.random.randint(low=0, high=numdatos, size=numdatos)
        promedios[i] = np.mean(datos[indice])
    promedios = np.sort(promedios)
    liminf = promedios[int(numiter * (alfa / 2))]
    limsup = promedios[int(numiter * (1 - alfa / 2))]
    return liminf, limsup

# Main Setup

datos = np.loadtxt('data/hubble_original.dat')
D = datos[:, 0]  # Primera columna data
V = datos[:, 1]  # Segunda columna data
H0ec1 = np.sum(V * D) / np.sum(D**2)  # H0 en función de V y D
H0ec2 = np.sum(V**2) / np.sum(V * D)
H0mean = (H0ec1 + H0ec2) / 2
recta = D * H0mean
fig1 = plt.figure(1)
fig1.clf()
plt.plot(D, V, 'r*')
plt.plot(D, recta, 'b-')
plt.xlabel('Distancia ($Mpc$)')
plt.ylabel(r'Velocidad ($\frac{km}{s}$)')
plt.title('Grafico ley de hubble con $H_0$ obtenida por Hubble')
plt.grid(True)
fig1.savefig('hubbleoriginal')
plt.show()

# Intervalo de confianza

alfa = 0.05
liminf, limsup = Bootstrap(V / D, alfa)
print 'H0 obtenido por Hubble = ' + str(H0mean) + ' [km/(s Mpc)]'
print ('El limite inferior para el intervalo de confianza para H0 es ' +
       str(liminf) +
       ' [km/(s Mpc)] y el limite superior es ' + str(limsup) +
       ' [km/(s Mpc)]')

# P2 Tarea 9

# Main Setup

datos = np.loadtxt('data/SNIa.dat', usecols=(1, 2))
V = datos[:, 0]  # Primera columna data
D = datos[:, 1]  # Segunda columna data
H0ec1 = np.sum(V * D) / np.sum(D**2)  # H0 en función de V y D
H0ec2 = np.sum(V**2) / np.sum(V * D)
H0mean = (H0ec1 + H0ec2) / 2
recta = D * H0mean
fig2 = plt.figure(2)
fig2.clf()
plt.plot(D, V, 'r*')
plt.plot(D, recta, 'g-')
plt.xlabel('Distancia ($Mpc$)')
plt.ylabel(r'Velocidad ($\frac{km}{s}$)')
plt.title(r'Grafico ley de hubble con $H_0$ corregida por Friedman')
plt.grid(True)
fig2.savefig('hubblefriedman')
plt.show()

# Intervalo de confianza

alfa = 0.05
liminf, limsup = Bootstrap(V / D, alfa)
print 'H0 obtenido por Friedman = ' + str(H0mean) + ' [km/(s Mpc)]'
print ('El limite inferior para el intervalo de confianza para H0 es ' +
       str(liminf) +
       ' [km/(s Mpc)] y el limite superior es ' + str(limsup) +
       ' [km/(s Mpc)]')

# P3 Tarea 9


def Montecarlo(datos1,  datos2, alfa, numejecuciones):
    '''Recibe un arreglo de datos y retorna un
    intervalo de confianza al 100*(1-alfa)% usando
    una simulacion de montecarlo'''
    flujo1 = datos1[0]
    error1 = datos1[1]
    flujo2 = datos2[0]
    error2 = datos2[1]
    a = np.zeros(numejecuciones)  # constante
    b = np.zeros(numejecuciones)  # pendiente
    for i in range(numejecuciones):
        r = np.random.normal(0, 1, size=len(flujo1))
        muestra1 = flujo1 + error1 * r
        muestra2 = flujo2 + error2 * r
        b[i], a[i] = np.polyfit(muestra1, muestra2, 1)
    b = np.sort(b)
    a = np.sort(a)
    liminfa = a[int(numejecuciones * alfa)]
    limsupa = a[int(numejecuciones * (1 - alfa))]
    liminfb = b[int(numejecuciones * alfa)]
    limsupb = b[int(numejecuciones * (1 - alfa))]
    return liminfa, limsupa, liminfb, limsupb

# Main Setup

datos = (3.631 *
         np.loadtxt('data/DR9Q.dat', usecols=(80, 81, 82, 83)))
a, b, c, d = datos[:, 0], datos[:, 1], datos[:, 2], datos[:, 3]
datosi = np.array([a, b])
datosz = np.array([c, d])
polyfit = np.polyfit(datosi[0], datosz[0], 1)
y = np.zeros(len(datosi[0]))
k = 0
for i in datosi[0]:
    y[k] = polyfit[0] * i + polyfit[1]
    k += 1
fig3 = plt.figure(3)
fig3.clf()
plt.plot(datosi[0], y, 'g-')
plt.plot(datosi[0], datosz[0], 'r*')
plt.xlabel('Flujo banda i ($1e-6 \ Jy$)')
plt.ylabel(r'Flujo banda z ($1e-6 \ Jy$)')
plt.title(r'Grafico flujo banda i vs flujo banda z')
plt.grid(True)
fig3.savefig('sdss')
plt.show()

# Intervalo de confianza

alfa = 0.05
numejecuciones = 10000
liminfa, limsupa, liminfb, limsupb = Montecarlo(datosi, datosz,
                                                alfa,
                                                numejecuciones)
print 'Parametro A = ' + str(polyfit[0])
print 'Parametro B = ' + str(polyfit[1])
print ('El limite inferior para el intervalo de ' +
       'confianza para el parametro A es ' +
       str(liminfb) +
       ' [1e-6 Jy] y el limite superior es ' + str(limsupb) +
       ' [1e-6 Jy]')
print ('El limite inferior para el intervalo de ' +
       'confianza para el parametro B es ' +
       str(liminfa) +
       ' [1e-6 Jy] y el limite superior es ' + str(limsupa) +
       ' [1e-6 Jy]')


# In[ ]:



