import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from configparser import ConfigParser


#LECTURA DE PARAMETROS DE ENTRADA

parser = ConfigParser()
parser.read('Entradas.ini') 

K = parser.getfloat('Parametros','K')
n = parser.getint('Parametros','n')

t0 = parser.getfloat('CondicionesIniciales','t0')
p0 = parser.getfloat('CondicionesIniciales','p0')


#defino arreglos de dimensión n, siendo n el número de iteraciones
tt = np.zeros(n) 
pp = np.zeros(n)
ee = np.zeros(n) 
xx = np.zeros(n)
delta = np.zeros(n)
Lyap = np.zeros(n)

#condiciones iniciales

t = t0
p = p0

#Separación inicial entre las trayectorias
d0=1.e-9
e=d0
x=d0
#e = eta = p' - p
#x = xi = t' - t

K=K/(2*pi)


for i in range(1,n):		
	#con mod(,1) me aseguro que p y t estén entre 0 y 1
	p = np.mod(p + K * np.sin(2.*pi*t),1.)
	t = np.mod(t + p,1.)
	
	e = e + K*np.cos(2.*pi*t)*x
	x = x + e
	
	pp[i] = p
	tt[i] = t
	
	xx[i] = x
	ee[i] = e
	
	delta[i]=np.sqrt(e**2+x**2)
	Lyap[i]=(1./i)*np.log(delta[i]/d0)

x=np.linspace(1,n,n)

#Defino un string que contiene los valores de las condiciones iniciales		
nombre='t0='+str(t0)+'. p0='+str(p0)



#para graficar el punto
tp0=[t0,p0]
plt.figure()
plt.plot(tt,pp,'.',markersize='5',color='k')	
plt.plot(tp0[0],tp0[1],'.',markersize='10',color='r')
plt.xlim(0,1)
plt.ylim(0,1)
plt.grid()
plt.xlabel('t')
plt.ylabel('p')
plt.title('MapaStandart.'+nombre)		
plt.savefig('MapaStandart.'+nombre+'.png')		


plt.show()

plt.figure()
plt.plot(x,Lyap,color='k',label='Lyapunov')	
plt.xlim(1,n)
plt.xlabel('n')
plt.ylabel(r'$\sigma$(n)')
plt.title('Exponente de Lyapunov.'+nombre)
plt.savefig('Lyapunov.'+nombre+'.png')
	
plt.figure()
plt.plot(x,delta,color='k',label=r'$\delta$(n)')	
plt.xlabel('n')
plt.ylabel(r'$\delta$(n)')
plt.xlim(1,n)
plt.title('Delta.'+nombre)
plt.savefig('Delta.'+nombre+'.png')
