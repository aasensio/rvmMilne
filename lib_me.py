#/usr/local/bin/ipython

import numpy as np
import pylab as pl
import pdb
from scipy.special import wofz

###########################################
# Returns the Voigt function for an axis of wavelengths l and damping parameter a
###########################################
def fvoigt(l, a):
        z = l + 1j*a
        return [wofz(z).real, wofz(z).imag/2.]

def frontl(x1,x2,x3,y1,y2,y3):
	l1=int(x1+x2-x3)
	l2=int(x2+x3-x1)
	l3=int(x3+x1-x2)
	l4=int(x1+x2+x3+1)
	l5=int(x1+y1)
	l6=int(x1-y1)
	l7=int(x2+y2)
	l8=int(x2-y2)
	l9=int(x3+y3)
	l10=int(x3-y3)

	dum = float(np.math.factorial(l1)) * float(np.math.factorial(l2)) * float(np.math.factorial(l3)) * float(np.math.factorial(l5)) / float(np.math.factorial(l4))
	dum = dum * float(np.math.factorial(l6)) * float(np.math.factorial(l7)) * float(np.math.factorial(l8)) * float(np.math.factorial(l9)) * float(np.math.factorial(l10))
	dum = np.sqrt(dum)

	return dum

def tresj(j1,j2,j3,m1,m2,m3):
	if (abs(m1) > j1): return 0.
	if (abs(m2) > j2): return 0.
	if (abs(m3) > j3): return 0.
	
	if (j1+j2-j3 < 0.): return 0.
	if (j2+j3-j1 < 0.): return 0.
	if (j3+j1-j2 < 0.): return 0.

	if (m1+m2+m3 != 0.): return 0.

	kmin1=j3-j1-m2
	kmin2=j3-j2+m1

	kmin = - np.array([kmin1,kmin2]).min()
	if (kmin < 0.): kmin = 0.
	kmax1=j1+j2-j3
	kmax2=j1-m1
	kmax3=j2+m2
	kmax = np.array([kmax1,kmax2,kmax3]).min()
	kmax=min([kmax1,kmax2,kmax3])

	if (kmin > kmax): return 0.
	term1 = frontl(j1,j2,j3,m1,m2,m3)
	msign = (-1)**(j1-j2-m3)
	suma = 0.
	for j in range(int(kmin), int(kmax)+1):
		term2 = np.math.factorial(j) * np.math.factorial(kmin1+j) * np.math.factorial(kmin2+j)
		term2 = term2 * np.math.factorial(kmax1-j) * np.math.factorial(kmax2-j) * np.math.factorial(kmax3-j)
		term = (-1)**j * msign * term1 / term2
		suma = suma + term

	return suma

def land_f(s, l, j):
	if (j == 0):
		return 0.
	else:
		return 3./2. + (s * (s + 1.) - l * (l + 1.)) / (2. * j * (j + 1.))

def lamb_vec(lamb0, n):
	deltal = 1. + (lamb0 - 5000.) * 1./10000.
	return np.linspace(lamb0-deltal, lamb0+deltal, num = n, endpoint=True)
	#return (np.arange(0,201)-100)*0.01 + lamb0

def at_elem(elem):
	if (elem.strip() == 'H'):
		return 1.
	elif (elem.strip() == 'HE'):
		return 4.
	elif (elem.strip() == 'LI'):
		return 7.
	elif (elem.strip() == 'BE'):
		return 9.
	elif (elem.strip() == 'B'):
		return 11.
	elif (elem.strip() == 'C'):
		return 12.
	elif (elem.strip() == 'N'):
		return 14.
	elif (elem.strip() == 'O'):
		return 16.
	elif (elem.strip() == 'F'):
		return 19.
	elif (elem.strip() == 'NE'):
		return 20.
	elif (elem.strip() == 'NA'):
		return 23.
	elif (elem.strip() == 'MG'):
		return 24.
	elif (elem.strip() == 'AL'):
		return 27.
	elif (elem.strip() == 'SI'):
		return 28.
	elif (elem.strip() == 'P'):
		return 31.
	elif (elem.strip() == 'S'):
		return 32.
	elif (elem.strip() == 'CL'):
		return 35.
	elif (elem.strip() == 'AR'):
		return 40.
	elif (elem.strip() == 'K'):
		return 39.
	elif (elem.strip() == 'CA'):
		return 40.
	elif (elem.strip() == 'SC'):
		return 45.
	elif (elem.strip() == 'TI'):
		return 48.
	elif (elem.strip() == 'V'):
		return 51.
	elif (elem.strip() == 'CR'):
		return 52.
	elif (elem.strip() == 'MN'):
		return 55.
	elif (elem.strip() == 'FE' or elem.strip() == 'XX'):
		return 56.
	elif (elem.strip() == 'CO'):
		return 59.
	elif (elem.strip() == 'NI'):
		return 59.
	elif (elem.strip() == 'CU'):
		return 64.
	elif (elem.strip() == 'ZN'):
		return 65.
	elif (elem.strip() == 'GA'):
		return 70.
        elif (elem.strip() == 'GE'):
                return 73.
        elif (elem.strip() == 'AS'):
                return 75.
        elif (elem.strip() == 'SE'):
                return 79.
        elif (elem.strip() == 'BR'):
                return 80.
        elif (elem.strip() == 'KR'):
                return 84.
	else:
		print 'UNKNOWN ATOMIC ELEMENT'
		return None


def ang_mom(mom):
	if (mom.strip() == 'S'):
		return 0.
	elif (mom.strip() == 'P'):
		return 1.
	elif (mom.strip() == 'D'):
		return 2.
	elif (mom.strip() == 'F'):
		return 3.
	elif (mom.strip() == 'G'):
		return 4.
	elif (mom.strip() == 'H'):
		return 5.
	else:
		print 'MOMENTO ANGULAR DESCONOCIDO'
		return None

def lineas(filename, ind_in):
	fdatos = open(filename, 'r')
	lineas = fdatos.readlines()
	for index in lineas:
		if (index.split('=')[0].strip() == ind_in):
			dum = index.split('=')[1]
			dum = dum.split()
			elem = dum[0]
			ions = dum[1]
			lamb = dum[2]
			vandw_factor = dum[3]
			potexc = dum[4]
			fuerza_osc = dum[5]
			dum2 = dum[6]
			s1 = dum2[0]
			l1 = dum2[1]
			j1 = dum[7].split('-')[0]
			dum2 = dum[8]
			s2 = dum2[0]
			l2 = dum2[1]
			j2 = dum[9].split('-')[0]
			try:
				param1 = dum[10]
				param2 = dum[11]
				param1 = float(param1)
				param2 = float(param2)
			except:
				print 'No aditional atomic parameters considered'
				param1 = None
				param2 = None
			l1 = ang_mom(l1)
			l2 = ang_mom(l2)
			mn = at_elem(elem)
			lamb = float(lamb)
			vandw_factor = float(vandw_factor)
			potexc = float(potexc)
			s1 = (float(s1) - 1.) / 2.
			l1 = float(l1)
			j1 = float(j1)
			s2 = (float(s2) - 1.) / 2.
			l2 = float(l2)
			j2 = float(j2)
			mn = float(mn)
			return elem, ions, lamb, vandw_factor, potexc, s1, l1, j1, s2, l2, j2, param1, param2, mn
	return None

def sintetizador(Bmag, deltalambdaD, omegam, theta, chi, br, adamp, eta0, lamb0, lambv, s2, l2, j2, s1, l1, j1, mn, noise=None):
	npoints = np.shape(lambv)
	#DETERMINAMOS LOS FACTORES DE LANDE DE CADA NIVEL
	g1 = land_f(s1, l1, j1)
	g2 = land_f(s2, l2, j2)
	#CONSTANTES UTILES:
	c = 2.99792458e+8
	mas = 1.66053886e-27*mn
	
	v = (lambv - lamb0) / deltalambdaD
	vm = omegam / (c) * lamb0 / deltalambdaD
	#INICIALIZAMOS VALORES PARA EL MOMENTO ANGULAR MAGNETICO DE CADA NIVEL:
	m1 = np.linspace(-j1, j1, num = 2*j1+1, endpoint=True)
	m2 = np.linspace(-j2, j2, num = 2*j2+1, endpoint=True)
	#INICIALIZAMOS LOS PERFILES DE ABSORCION Y EMISION??
	phir = np.zeros(npoints)
	phip = np.zeros(npoints)
	phib = np.zeros(npoints)
	psir = np.zeros(npoints)
	psip = np.zeros(npoints)
	psib = np.zeros(npoints)

	for index in range(2*int(j1)+1):
		#SIGMA ROJA:
		mdum1 = m1[index]
		caso = np.where(m2 - mdum1 == -1)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			phi_r0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phir = phir + phi_r0
			psi_r0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psir = psir + psi_r0
		#PI:
		caso = np.where(m2 - mdum1 == 0)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			phi_p0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phip = phip + phi_p0
			psi_p0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psip = psip + psi_p0
		#SIGMA AZUL:
		caso = np.where(m2 - mdum1 == 1)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			phi_b0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phib = phib + phi_b0
			psi_b0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psib = psib + psi_b0

	etai = 1. + eta0 / 2. * (phip * np.sin(theta)**2. + 1. / 2. * (phib + phir) * (1. + np.cos(theta)**2.))
	etaq = eta0 / 2. * (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.cos(2. * chi)
	etau = eta0 / 2. * (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.sin(2. * chi)
	etav = eta0 / 2. * (phir - phib) * np.cos(theta)
	rhoq = eta0 / 2. * (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.cos(2. * chi)
	rhou = eta0 / 2. * (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.sin(2. * chi)
	rhov = eta0 / 2. * (psir - psib) * np.cos(theta)

	deltaM = etai**2. * (etai**2. - etaq**2. - etau**2. - etav**2. + rhoq**2. + rhou**2. + rhov**2.) - (etaq * rhoq + etau * rhou + etav * rhov)**2.

	stoki = (1. + br * deltaM**(-1.) * etai * (etai**2. + rhoq**2. + rhou**2. + rhov**2.)) / (1. + br)
	stokq = - br / (1. + br) * deltaM**(-1.) * (etai**2. * etaq + etai * (etav * rhou - etau * rhov) + rhoq * (etaq * rhoq + etau * rhou + etav * rhov))
	stoku = - br / (1. + br) * deltaM**(-1.) * (etai**2. * etau + etai * (etaq * rhov - etav * rhoq) + rhou * (etaq * rhoq + etau * rhou + etav * rhov))
	stokv = - br / (1. + br) * deltaM**(-1.) * (etai**2. * etav + etai * (etau * rhoq - etaq * rhou) + rhov * (etaq * rhoq + etau * rhou + etav * rhov))

	if (noise != None):
		stoki = stoki + np.random.randn(npoints[0]) * noise
		stokq = stokq + np.random.randn(npoints[0]) * noise
		stoku = stoku + np.random.randn(npoints[0]) * noise
		stokv = stokv + np.random.randn(npoints[0]) * noise

	return stoki, stokq, stoku, stokv

#RESPECTO DE Br, eta0 y theta

def derivador(Bmag, deltalambdaD, omegam, theta, chi, br, adamp, eta0, lamb0, lambv, s2, l2, j2, s1, l1, j1, mn):
	npoints = np.shape(lambv)
	#DETERMINAMOS LOS FACTORES DE LANDE DE CADA NIVEL
	g1 = land_f(s1, l1, j1)
	g2 = land_f(s2, l2, j2)
	#CONSTANTES UTILES:
	c = 2.99792458e+8
	mas = 1.66053886e-27*mn

	v = (lambv - lamb0) / deltalambdaD
	vm = omegam / (c) * lamb0 / deltalambdaD
	#INICIALIZAMOS VALORES PARA EL MOMENTO ANGULAR MAGNETICO DE CADA NIVEL:
	m1 = np.linspace(-j1, j1, num = 2*j1+1, endpoint=True)
	m2 = np.linspace(-j2, j2, num = 2*j2+1, endpoint=True)
	#INICIALIZAMOS LOS PERFILES DE ABSORCION Y EMISION??
	phir = np.zeros(npoints)
	phip = np.zeros(npoints)
	phib = np.zeros(npoints)
	psir = np.zeros(npoints)
	psip = np.zeros(npoints)
	psib = np.zeros(npoints)
	der_phir_db = np.zeros(npoints)
	der_psir_db = np.zeros(npoints)
	der_phip_db = np.zeros(npoints)
	der_psip_db = np.zeros(npoints)
	der_phib_db = np.zeros(npoints)
	der_psib_db = np.zeros(npoints)
	der_phir_dadamp = np.zeros(npoints)
	der_psir_dadamp = np.zeros(npoints)
	der_phip_dadamp = np.zeros(npoints)
	der_psip_dadamp = np.zeros(npoints)
	der_phib_dadamp = np.zeros(npoints)
	der_psib_dadamp = np.zeros(npoints)
	der_phir_ddeltalambdaD = np.zeros(npoints)
	der_psir_ddeltalambdaD = np.zeros(npoints)
	der_phip_ddeltalambdaD = np.zeros(npoints)
	der_psip_ddeltalambdaD = np.zeros(npoints)
	der_phib_ddeltalambdaD = np.zeros(npoints)
	der_psib_ddeltalambdaD = np.zeros(npoints)
	der_phir_dvm = np.zeros(npoints)
	der_psir_dvm = np.zeros(npoints)
	der_phip_dvm = np.zeros(npoints)
	der_psip_dvm = np.zeros(npoints)
	der_phib_dvm = np.zeros(npoints)
	der_psib_dvm = np.zeros(npoints)

	der_vprime_dvm = -lamb0 / (c) / deltalambdaD
	for index in range(2*int(j1)+1):
		#SIGMA ROJA:
		mdum1 = m1[index]
		caso = np.where(m2 - mdum1 == -1)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			#*********************************************************************
			#DERIVADA RESPECTO DE B
			der_vprime_db = -(g1 * mdum1 - g2 * mdum2) * 4.67e-13 * lamb0**2. / deltalambdaD
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc
			der_hfunc_db = der_hfunc_dvprime * der_vprime_db
			der_lfunc_db = der_lfunc_dvprime * der_vprime_db

			der_phir_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_db
			der_psir_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_db

			der_phir_db = der_phir_db + der_phir_db0
			der_psir_db = der_psir_db + der_psir_db0
			#RESPECTO DE ADAMP
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc

			der_hfunc_dadamp = -2. * der_lfunc_dvprime
			der_lfunc_dadamp = 1. / 2. * der_hfunc_dvprime

			der_phir_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dadamp
			der_psir_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dadamp

			der_phir_dadamp = der_phir_dadamp + der_phir_dadamp0
			der_psir_dadamp = der_psir_dadamp + der_psir_dadamp0
			#RESPECTO DE DELTALAMBDAD
			der_vprime_ddeltalambdaD = -vprime / deltalambdaD

			der_hfunc_ddeltalambdaD = der_hfunc_dvprime * der_vprime_ddeltalambdaD
			der_lfunc_ddeltalambdaD = der_lfunc_dvprime * der_vprime_ddeltalambdaD

			der_phir_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_hfunc_ddeltalambdaD - 1. / deltalambdaD**2. * hfunc)
			der_psir_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_lfunc_ddeltalambdaD - 1. / deltalambdaD**2. * lfunc)

			der_phir_ddeltalambdaD = der_phir_ddeltalambdaD + der_phir_ddeltalambdaD0
			der_psir_ddeltalambdaD = der_psir_ddeltalambdaD + der_psir_ddeltalambdaD0
			#RESPECTO DE vm
			der_hfunc_dvm = der_hfunc_dvprime * der_vprime_dvm
			der_lfunc_dvm = der_lfunc_dvprime * der_vprime_dvm

			der_phir_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dvm
			der_psir_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dvm

			der_phir_dvm = der_phir_dvm + der_phir_dvm0
			der_psir_dvm = der_psir_dvm + der_psir_dvm0
			#*********************************************************************
			phi_r0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phir = phir + phi_r0
			psi_r0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psir = psir + psi_r0

		#PI:
		caso = np.where(m2 - mdum1 == 0)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			#*********************************************************************
			#DERIVADA RESPECTO DE B
			der_vprime_db = 0.
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc
			der_hfunc_db = der_hfunc_dvprime * der_vprime_db
			der_lfunc_db = der_lfunc_dvprime * der_vprime_db

			der_phip_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_db
			der_psip_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_db

			der_phip_db = der_phip_db + der_phip_db0
			der_psip_db = der_psip_db + der_psip_db0
			#DERIVADA RESPECTO DE ADAMP
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc

			der_hfunc_dadamp = -2. * der_lfunc_dvprime
			der_lfunc_dadamp = 1. / 2. * der_hfunc_dvprime

			der_phip_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dadamp
			der_psip_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dadamp

			der_phip_dadamp = der_phip_dadamp + der_phip_dadamp0
			der_psip_dadamp = der_psip_dadamp + der_psip_dadamp0
			#RESPECTO DE DELTALAMBDAD
			der_vprime_ddeltalambdaD = -vprime / deltalambdaD

			der_hfunc_ddeltalambdaD = der_hfunc_dvprime * der_vprime_ddeltalambdaD
			der_lfunc_ddeltalambdaD = der_lfunc_dvprime * der_vprime_ddeltalambdaD

			der_phip_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_hfunc_ddeltalambdaD - 1. / deltalambdaD**2. * hfunc)
			der_psip_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_lfunc_ddeltalambdaD - 1. / deltalambdaD**2. * lfunc)

			der_phip_ddeltalambdaD = der_phip_ddeltalambdaD + der_phip_ddeltalambdaD0
			der_psip_ddeltalambdaD = der_psip_ddeltalambdaD + der_psip_ddeltalambdaD0
			#RESPECTO DE vm
			der_hfunc_dvm = der_hfunc_dvprime * der_vprime_dvm
			der_lfunc_dvm = der_lfunc_dvprime * der_vprime_dvm

			der_phip_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dvm
			der_psip_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dvm

			der_phip_dvm = der_phip_dvm + der_phip_dvm0
			der_psip_dvm = der_psip_dvm + der_psip_dvm0

			#*********************************************************************
			phi_p0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phip = phip + phi_p0
			psi_p0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psip = psip + psi_p0

		#SIGMA AZUL:
		caso = np.where(m2 - mdum1 == 1)
		if (np.size(caso) != 0):
			mdum2 = m2[caso]
			mdum2 = mdum2[0]
			delta_lm1m2 = (g1 * mdum1 - g2 * mdum2) * 4.67e-13 * Bmag * lamb0**2.
			vb = delta_lm1m2 / deltalambdaD
			vprime = v - vm - vb
			ff=3.e0*tresj(j1,j2,1.e0,mdum1,-mdum2,mdum2-mdum1)**2.
			hfunc, lfunc = fvoigt(vprime, adamp)
			#*********************************************************************
			#DERIVADA RESPECTO DE B
			der_vprime_db = -(g1 * mdum1 - g2 * mdum2) * 4.67e-13 * lamb0**2. / deltalambdaD
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc
			der_hfunc_db = der_hfunc_dvprime * der_vprime_db
			der_lfunc_db = der_lfunc_dvprime * der_vprime_db

			der_phib_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_db
			der_psib_db0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_db

			der_phib_db = der_phib_db + der_phib_db0
			der_psib_db = der_psib_db + der_psib_db0
			#DERIVADA RESPECTO DE ADAMP
			der_hfunc_dvprime = -2.*vprime*hfunc+4.*adamp*lfunc
			der_lfunc_dvprime = 1. / np.sqrt(np.pi) - adamp * hfunc - 2. * vprime * lfunc

			der_hfunc_dadamp = -2. * der_lfunc_dvprime
			der_lfunc_dadamp = 1. / 2. * der_hfunc_dvprime

			der_phib_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dadamp
			der_psib_dadamp0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dadamp

			der_phib_dadamp = der_phib_dadamp + der_phib_dadamp0
			der_psib_dadamp = der_psib_dadamp + der_psib_dadamp0
			#DERIVADA RESPECTO DE DELTALAMBDAD
			der_vprime_ddeltalambdaD = -vprime / deltalambdaD

			der_hfunc_ddeltalambdaD = der_hfunc_dvprime * der_vprime_ddeltalambdaD
			der_lfunc_ddeltalambdaD = der_lfunc_dvprime * der_vprime_ddeltalambdaD

			der_phib_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_hfunc_ddeltalambdaD - 1. / deltalambdaD**2. * hfunc)
			der_psib_ddeltalambdaD0 = 1. / np.sqrt(np.pi) * ff * (1. / deltalambdaD * der_lfunc_ddeltalambdaD - 1. / deltalambdaD**2. * lfunc)

			der_phib_ddeltalambdaD = der_phib_ddeltalambdaD + der_phib_ddeltalambdaD0
			der_psib_ddeltalambdaD = der_psib_ddeltalambdaD + der_psib_ddeltalambdaD0
			#RESPECTO DE vm
			der_hfunc_dvm = der_hfunc_dvprime * der_vprime_dvm
			der_lfunc_dvm = der_lfunc_dvprime * der_vprime_dvm

			der_phib_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_hfunc_dvm
			der_psib_dvm0 = 1. / np.sqrt(np.pi) * 1. / deltalambdaD * ff * der_lfunc_dvm

			der_phib_dvm = der_phib_dvm + der_phib_dvm0
			der_psib_dvm = der_psib_dvm + der_psib_dvm0
			#*********************************************************************
			phi_b0 = ff * hfunc / deltalambdaD / np.sqrt(np.pi)
			phib = phib + phi_b0
			psi_b0 = ff * lfunc / deltalambdaD / np.sqrt(np.pi)
			psib = psib + psi_b0

	#****************************************************************************************************
	#DERIVADAS DE LOS ETA Y RHO RESPECTO DE LOS PHI Y PSI
	der_etai_dphir = eta0 / 4. * (1. + np.cos(theta)**2.)
	der_etai_dphip = eta0 / 2. * np.sin(theta)**2.
	der_etai_dphib = eta0 / 4. * (1. + np.cos(theta)**2.)
	der_etaq_dphir = -eta0 / 4. * np.sin(theta)**2. * np.cos(2. * chi)
	der_etaq_dphip = eta0 / 2. * np.sin(theta)**2. * np.cos(2. * chi)
	der_etaq_dphib = -eta0 / 4. * np.sin(theta)**2. * np.cos(2. * chi)
	der_etau_dphir = -eta0 / 4. * np.sin(theta)**2. * np.sin(2. * chi)
	der_etau_dphip = eta0 / 2. * np.sin(theta)**2. * np.sin(2. * chi)
	der_etau_dphib = -eta0 / 4. * np.sin(theta)**2. * np.sin(2. * chi)
	der_etav_dphir = eta0 / 2. * np.cos(theta)
	der_etav_dphip = 0.
	der_etav_dphib = -eta0 / 2. * np.cos(theta)
	der_rhoq_dpsir = -eta0 / 4. * np.sin(theta)**2. * np.cos(2. * chi)
	der_rhoq_dpsip = eta0 / 2. * np.sin(theta)**2. * np.cos(2. * chi)
	der_rhoq_dpsib = -eta0 / 4. * np.sin(theta)**2. * np.cos(2. * chi)
	der_rhou_dpsir = -eta0 / 4. * np.sin(theta)**2. * np.sin(2. * chi)
	der_rhou_dpsip = eta0 / 2. * np.sin(theta)**2. * np.sin(2. * chi)
	der_rhou_dpsib = -eta0 / 4. * np.sin(theta)**2. * np.sin(2. * chi)
	der_rhov_dpsir = eta0 / 2. * np.cos(theta)
	der_rhov_dpsip = 0.
	der_rhov_dpsib = -eta0 / 2. * np.cos(theta)
	#****************************************************************************************************
	#DERIVADAS RESPECTO RESPECTO DE LOS ETA Y LOS RHO
	der_etai_deta0 = 1. / 2. * (phip * np.sin(theta)**2. + 1. / 2. * (phib + phir) * (1. + np.cos(theta)**2.))
	der_etaq_deta0 = (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.cos(2. * chi) / 2.
	der_etau_deta0 = (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.sin(2. * chi) / 2.
	der_etav_deta0 = (phir - phib) * np.cos(theta) / 2.
	der_rhoq_deta0 = (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.cos(2. * chi) / 2.
	der_rhou_deta0 = (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.sin(2. * chi) / 2.
	der_rhov_deta0 = (psir - psib) * np.cos(theta) / 2.

	der_etai_dtheta = eta0 * np.sin(theta) * np.cos(theta) * (phip - (phir + phib) / 2.)
	der_etaq_dtheta = eta0 * np.cos(theta) * np.sin(theta) * np.cos(2. * chi) * (phip - (phir + phib) / 2.)
	der_etau_dtheta = eta0 * np.cos(theta) * np.sin(theta) * np.sin(2. * chi) * (phip - (phir + phib) / 2.)
	der_etav_dtheta = -eta0 / 2. * np.sin(theta) * (phir - phib)
	der_rhoq_dtheta = eta0 * np.sin(theta) * np.cos(theta) * np.cos(2. * chi) * (psip - (psir + psib) / 2.)
	der_rhou_dtheta = eta0 * np.sin(theta) * np.cos(theta) * np.sin(2. * chi) * (psip - (psir + psib) / 2.)
	der_rhov_dtheta = -eta0 / 2. * np.sin(theta) * (psir - psib)

	der_etai_dchi = np.zeros(len(lambv))
	der_etaq_dchi = -eta0 * np.sin(theta)**2. * np.sin(2. * chi) * (phip - (phir + phib) / 2.)
	der_etau_dchi = eta0 * np.sin(theta)**2. * np.cos(2. * chi) * (phip - (phir + phib) / 2.)
	der_etav_dchi = np.zeros(len(lambv))
	der_rhoq_dchi = -eta0 * np.sin(theta)**2. * np.sin(2. * chi) * (psip - (psir + psib) / 2.)
	der_rhou_dchi = eta0 * np.sin(theta)**2. * np.cos(2. * chi) * (psip - (psir + psib) / 2.)
	der_rhov_dchi = np.zeros(len(lambv))

	der_etai_db = der_etai_dphir * der_phir_db + der_etai_dphip * der_phip_db + der_etai_dphib * der_phib_db
	der_etaq_db = der_etaq_dphir * der_phir_db + der_etaq_dphip * der_phip_db + der_etaq_dphib * der_phib_db
	der_etau_db = der_etau_dphir * der_phir_db + der_etau_dphip * der_phip_db + der_etau_dphib * der_phib_db
	der_etav_db = der_etav_dphir * der_phir_db + der_etav_dphip * der_phip_db + der_etav_dphib * der_phib_db
	der_rhoq_db = der_rhoq_dpsir * der_psir_db + der_rhoq_dpsip * der_psip_db + der_rhoq_dpsib * der_psib_db
	der_rhou_db = der_rhou_dpsir * der_psir_db + der_rhou_dpsip * der_psip_db + der_rhou_dpsib * der_psib_db
	der_rhov_db = der_rhov_dpsir * der_psir_db + der_rhov_dpsip * der_psip_db + der_rhov_dpsib * der_psib_db

	der_etai_dadamp = der_etai_dphir * der_phir_dadamp + der_etai_dphip * der_phip_dadamp + der_etai_dphib * der_phib_dadamp
	der_etaq_dadamp = der_etaq_dphir * der_phir_dadamp + der_etaq_dphip * der_phip_dadamp + der_etaq_dphib * der_phib_dadamp
	der_etau_dadamp = der_etau_dphir * der_phir_dadamp + der_etau_dphip * der_phip_dadamp + der_etau_dphib * der_phib_dadamp
	der_etav_dadamp = der_etav_dphir * der_phir_dadamp + der_etav_dphip * der_phip_dadamp + der_etav_dphib * der_phib_dadamp
	der_rhoq_dadamp = der_rhoq_dpsir * der_psir_dadamp + der_rhoq_dpsip * der_psip_dadamp + der_rhoq_dpsib * der_psib_dadamp
	der_rhou_dadamp = der_rhou_dpsir * der_psir_dadamp + der_rhou_dpsip * der_psip_dadamp + der_rhou_dpsib * der_psib_dadamp
	der_rhov_dadamp = der_rhov_dpsir * der_psir_dadamp + der_rhov_dpsip * der_psip_dadamp + der_rhov_dpsib * der_psib_dadamp

	der_etai_ddeltalambdaD = der_etai_dphir * der_phir_ddeltalambdaD + der_etai_dphip * der_phip_ddeltalambdaD + der_etai_dphib * der_phib_ddeltalambdaD
	der_etaq_ddeltalambdaD = der_etaq_dphir * der_phir_ddeltalambdaD + der_etaq_dphip * der_phip_ddeltalambdaD + der_etaq_dphib * der_phib_ddeltalambdaD
	der_etau_ddeltalambdaD = der_etau_dphir * der_phir_ddeltalambdaD + der_etau_dphip * der_phip_ddeltalambdaD + der_etau_dphib * der_phib_ddeltalambdaD
	der_etav_ddeltalambdaD = der_etav_dphir * der_phir_ddeltalambdaD + der_etav_dphip * der_phip_ddeltalambdaD + der_etav_dphib * der_phib_ddeltalambdaD
	der_rhoq_ddeltalambdaD = der_rhoq_dpsir * der_psir_ddeltalambdaD + der_rhoq_dpsip * der_psip_ddeltalambdaD + der_rhoq_dpsib * der_psib_ddeltalambdaD
	der_rhou_ddeltalambdaD = der_rhou_dpsir * der_psir_ddeltalambdaD + der_rhou_dpsip * der_psip_ddeltalambdaD + der_rhou_dpsib * der_psib_ddeltalambdaD
	der_rhov_ddeltalambdaD = der_rhov_dpsir * der_psir_ddeltalambdaD + der_rhov_dpsip * der_psip_ddeltalambdaD + der_rhov_dpsib * der_psib_ddeltalambdaD

	der_etai_dvm = der_etai_dphir * der_phir_dvm + der_etai_dphip * der_phip_dvm + der_etai_dphib * der_phib_dvm
	der_etaq_dvm = der_etaq_dphir * der_phir_dvm + der_etaq_dphip * der_phip_dvm + der_etaq_dphib * der_phib_dvm
	der_etau_dvm = der_etau_dphir * der_phir_dvm + der_etau_dphip * der_phip_dvm + der_etau_dphib * der_phib_dvm
	der_etav_dvm = der_etav_dphir * der_phir_dvm + der_etav_dphip * der_phip_dvm + der_etav_dphib * der_phib_dvm
	der_rhoq_dvm = der_rhoq_dpsir * der_psir_dvm + der_rhoq_dpsip * der_psip_dvm + der_rhoq_dpsib * der_psib_dvm
	der_rhou_dvm = der_rhou_dpsir * der_psir_dvm + der_rhou_dpsip * der_psip_dvm + der_rhou_dpsib * der_psib_dvm
	der_rhov_dvm = der_rhov_dpsir * der_psir_dvm + der_rhov_dpsip * der_psip_dvm + der_rhov_dpsib * der_psib_dvm
	#****************************************************************************************************

	etai = 1. + eta0 / 2. * (phip * np.sin(theta)**2. + 1. / 2. * (phib + phir) * (1. + np.cos(theta)**2.))
	etaq = eta0 / 2. * (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.cos(2. * chi)
	etau = eta0 / 2. * (phip - 1. / 2. * (phib + phir)) * np.sin(theta)**2. * np.sin(2. * chi)
	etav = eta0 / 2. * (phir - phib) * np.cos(theta)
	rhoq = eta0 / 2. * (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.cos(2. * chi)
	rhou = eta0 / 2. * (psip - 1. / 2. * (psib + psir)) * np.sin(theta)**2. * np.sin(2. * chi)
	rhov = eta0 / 2. * (psir - psib) * np.cos(theta)

	#****************************************************************************************************
	#DERIVADAS DE DELTA RESPECTO DE ETA Y RHO
	der_delta_detai = 4. * etai**3. + 2. * etai * (rhoq**2. + rhou**2. + rhov**2. - etaq**2. - etau**2. - etav**2.)
	der_delta_detaq = -2. * etai**2. * etaq - 2. * rhoq * (etaq * rhoq + etau * rhou + etav * rhov)
	der_delta_detau = -2. * etai**2. * etau - 2. * rhou * (etaq * rhoq + etau * rhou + etav * rhov)
	der_delta_detav = -2. * etai**2. * etav - 2. * rhov * (etaq * rhoq + etau * rhou + etav * rhov)
	der_delta_drhoq = 2. * etai**2. * rhoq - 2. * etaq * (etaq * rhoq + etau * rhou + etav * rhov)
	der_delta_drhou = 2. * etai**2. * rhou - 2. * etau * (etaq * rhoq + etau * rhou + etav * rhov)
	der_delta_drhov = 2. * etai**2. * rhov - 2. * etav * (etaq * rhoq + etau * rhou + etav * rhov)
	#****************************************************************************************************
	#DERIVADA DE DELTA RESPECTO DE ETA0
	der_delta_deta0 = der_delta_detai * der_etai_deta0 + \
		der_delta_detaq * der_etaq_deta0 + \
		der_delta_detau * der_etau_deta0 + \
		der_delta_detav * der_etav_deta0 + \
		der_delta_drhoq * der_rhoq_deta0 + \
		der_delta_drhou * der_rhou_deta0 + \
		der_delta_drhov * der_rhov_deta0

	der_delta_dtheta = der_delta_detai * der_etai_dtheta + \
		der_delta_detaq * der_etaq_dtheta + \
		der_delta_detau * der_etau_dtheta + \
		der_delta_detav * der_etav_dtheta + \
		der_delta_drhoq * der_rhoq_dtheta + \
		der_delta_drhou * der_rhou_dtheta + \
		der_delta_drhov * der_rhov_dtheta

	der_delta_dchi = der_delta_detai * der_etai_dchi + \
		der_delta_detaq * der_etaq_dchi + \
		der_delta_detau * der_etau_dchi + \
		der_delta_detav * der_etav_dchi + \
		der_delta_drhoq * der_rhoq_dchi + \
		der_delta_drhou * der_rhou_dchi + \
		der_delta_drhov * der_rhov_dchi

	der_delta_db = der_delta_detai * der_etai_db + \
		der_delta_detaq * der_etaq_db + \
		der_delta_detau * der_etau_db + \
		der_delta_detav * der_etav_db + \
		der_delta_drhoq * der_rhoq_db + \
		der_delta_drhou * der_rhou_db + \
		der_delta_drhov * der_rhov_db

	der_delta_dadamp = der_delta_detai * der_etai_dadamp + \
		der_delta_detaq * der_etaq_dadamp + \
		der_delta_detau * der_etau_dadamp + \
		der_delta_detav * der_etav_dadamp + \
		der_delta_drhoq * der_rhoq_dadamp + \
		der_delta_drhou * der_rhou_dadamp + \
		der_delta_drhov * der_rhov_dadamp

	der_delta_ddeltalambdaD = der_delta_detai * der_etai_ddeltalambdaD + \
		der_delta_detaq * der_etaq_ddeltalambdaD + \
		der_delta_detau * der_etau_ddeltalambdaD + \
		der_delta_detav * der_etav_ddeltalambdaD + \
		der_delta_drhoq * der_rhoq_ddeltalambdaD + \
		der_delta_drhou * der_rhou_ddeltalambdaD + \
		der_delta_drhov * der_rhov_ddeltalambdaD

	der_delta_dvm = der_delta_detai * der_etai_dvm + \
		der_delta_detaq * der_etaq_dvm + \
		der_delta_detau * der_etau_dvm + \
		der_delta_detav * der_etav_dvm + \
		der_delta_drhoq * der_rhoq_dvm + \
		der_delta_drhou * der_rhou_dvm + \
		der_delta_drhov * der_rhov_dvm

	#****************************************************************************************************
	#DERIVADA DE LAS "CONSTANTES" RESPECTO DE LOS ETA Y RHO:
	der_csi_detai = 3. * etai**2. + rhoq**2. + rhou**2. + rhov**2.
	der_csi_detaq = np.zeros(len(der_csi_detai))
	der_csi_detau = np.zeros(len(der_csi_detai))
	der_csi_detav = np.zeros(len(der_csi_detai))
	der_csi_drhoq = 2. * etai * rhoq
	der_csi_drhou = 2. * etai * rhou
	der_csi_drhov = 2. * etai * rhov

	der_csq_detai = 2. * etai * etaq + (etav * rhou - etau * rhov)
	der_csq_detaq = etai**2. + etaq**2.
	der_csq_detau = rhoq * rhou - rhov * etai
	der_csq_detav = etai * rhou + rhoq * rhov
	der_csq_drhoq = 2. * etaq * rhoq + etau * rhou + etav * rhov
	der_csq_drhou = etai * etav + rhoq * etau
	der_csq_drhov = rhoq * etav - etai * etau

	der_csu_detai = 2. * etai * etau + (etaq * rhov - etav * rhoq)
	der_csu_detaq = etai * rhov + rhou * rhoq
	der_csu_detau = etai**2. + rhou**2.
	der_csu_detav = rhou * rhov - etai * rhoq
	der_csu_drhoq = rhou * etaq - etai * etav
	der_csu_drhou = etaq * rhoq + 2. * etau * rhou + etav * rhov
	der_csu_drhov = etai * etaq + rhou * rhov

	der_csv_detai = 2. * etai * etav + (etau * rhoq - etaq * rhou)
	der_csv_detaq = rhov * rhoq - etai * rhou
	der_csv_detau = etai * rhoq + rhov * rhou
	der_csv_detav = etai**2. + rhov**2.
	der_csv_drhoq = etai * etau + rhov * etaq
	der_csv_drhou = rhov * etau - etai * etaq
	der_csv_drhov = etaq * rhoq + etau * rhou + 2. * etav * rhov

	#****************************************************************************************************
	#DERIVADAS DE LAS "CONSTANTES" RESPECTO DE ETA0
	der_csi_deta0 = der_csi_detai * der_etai_deta0 + \
		der_csi_detaq * der_etaq_deta0 + \
		der_csi_detau * der_etau_deta0 + \
		der_csi_detav * der_etav_deta0 + \
		der_csi_drhoq * der_rhoq_deta0 + \
		der_csi_drhou * der_rhou_deta0 + \
		der_csi_drhov * der_rhov_deta0
	der_csi_dtheta = der_csi_detai * der_etai_dtheta + \
		der_csi_detaq * der_etaq_dtheta + \
		der_csi_detau * der_etau_dtheta + \
		der_csi_detav * der_etav_dtheta + \
		der_csi_drhoq * der_rhoq_dtheta + \
		der_csi_drhou * der_rhou_dtheta + \
		der_csi_drhov * der_rhov_dtheta
	der_csi_dchi = der_csi_detai * der_etai_dchi + \
		der_csi_detaq * der_etaq_dchi + \
		der_csi_detau * der_etau_dchi + \
		der_csi_detav * der_etav_dchi + \
		der_csi_drhoq * der_rhoq_dchi + \
		der_csi_drhou * der_rhou_dchi + \
		der_csi_drhov * der_rhov_dchi
	der_csi_db = der_csi_detai * der_etai_db + \
		der_csi_detaq * der_etaq_db + \
		der_csi_detau * der_etau_db + \
		der_csi_detav * der_etav_db + \
		der_csi_drhoq * der_rhoq_db + \
		der_csi_drhou * der_rhou_db + \
		der_csi_drhov * der_rhov_db
	der_csi_dadamp = der_csi_detai * der_etai_dadamp + \
		der_csi_detaq * der_etaq_dadamp + \
		der_csi_detau * der_etau_dadamp + \
		der_csi_detav * der_etav_dadamp + \
		der_csi_drhoq * der_rhoq_dadamp + \
		der_csi_drhou * der_rhou_dadamp + \
		der_csi_drhov * der_rhov_dadamp
	der_csi_ddeltalambdaD = der_csi_detai * der_etai_ddeltalambdaD + \
		der_csi_detaq * der_etaq_ddeltalambdaD + \
		der_csi_detau * der_etau_ddeltalambdaD + \
		der_csi_detav * der_etav_ddeltalambdaD + \
		der_csi_drhoq * der_rhoq_ddeltalambdaD + \
		der_csi_drhou * der_rhou_ddeltalambdaD + \
		der_csi_drhov * der_rhov_ddeltalambdaD
	der_csi_dvm = der_csi_detai * der_etai_dvm + \
		der_csi_detaq * der_etaq_dvm + \
		der_csi_detau * der_etau_dvm + \
		der_csi_detav * der_etav_dvm + \
		der_csi_drhoq * der_rhoq_dvm + \
		der_csi_drhou * der_rhou_dvm + \
		der_csi_drhov * der_rhov_dvm

	der_csq_deta0 = der_csq_detai * der_etai_deta0 + \
		der_csq_detaq * der_etaq_deta0 + \
		der_csq_detau * der_etau_deta0 + \
		der_csq_detav * der_etav_deta0 + \
		der_csq_drhoq * der_rhoq_deta0 + \
		der_csq_drhou * der_rhou_deta0 + \
		der_csq_drhov * der_rhov_deta0
	der_csq_dtheta = der_csq_detai * der_etai_dtheta + \
		der_csq_detaq * der_etaq_dtheta + \
		der_csq_detau * der_etau_dtheta + \
		der_csq_detav * der_etav_dtheta + \
		der_csq_drhoq * der_rhoq_dtheta + \
		der_csq_drhou * der_rhou_dtheta + \
		der_csq_drhov * der_rhov_dtheta
	der_csq_dchi = der_csq_detai * der_etai_dchi + \
		der_csq_detaq * der_etaq_dchi + \
		der_csq_detau * der_etau_dchi + \
		der_csq_detav * der_etav_dchi + \
		der_csq_drhoq * der_rhoq_dchi + \
		der_csq_drhou * der_rhou_dchi + \
		der_csq_drhov * der_rhov_dchi
	der_csq_db = der_csq_detai * der_etai_db + \
		der_csq_detaq * der_etaq_db + \
		der_csq_detau * der_etau_db + \
		der_csq_detav * der_etav_db + \
		der_csq_drhoq * der_rhoq_db + \
		der_csq_drhou * der_rhou_db + \
		der_csq_drhov * der_rhov_db
	der_csq_dadamp = der_csq_detai * der_etai_dadamp + \
		der_csq_detaq * der_etaq_dadamp + \
		der_csq_detau * der_etau_dadamp + \
		der_csq_detav * der_etav_dadamp + \
		der_csq_drhoq * der_rhoq_dadamp + \
		der_csq_drhou * der_rhou_dadamp + \
		der_csq_drhov * der_rhov_dadamp
	der_csq_ddeltalambdaD = der_csq_detai * der_etai_ddeltalambdaD + \
		der_csq_detaq * der_etaq_ddeltalambdaD + \
		der_csq_detau * der_etau_ddeltalambdaD + \
		der_csq_detav * der_etav_ddeltalambdaD + \
		der_csq_drhoq * der_rhoq_ddeltalambdaD + \
		der_csq_drhou * der_rhou_ddeltalambdaD + \
		der_csq_drhov * der_rhov_ddeltalambdaD
	der_csq_dvm = der_csq_detai * der_etai_dvm + \
		der_csq_detaq * der_etaq_dvm + \
		der_csq_detau * der_etau_dvm + \
		der_csq_detav * der_etav_dvm + \
		der_csq_drhoq * der_rhoq_dvm + \
		der_csq_drhou * der_rhou_dvm + \
		der_csq_drhov * der_rhov_dvm

	der_csu_deta0 = der_csu_detai * der_etai_deta0 + \
		der_csu_detaq * der_etaq_deta0 + \
		der_csu_detau * der_etau_deta0 + \
		der_csu_detav * der_etav_deta0 + \
		der_csu_drhoq * der_rhoq_deta0 + \
		der_csu_drhou * der_rhou_deta0 + \
		der_csu_drhov * der_rhov_deta0
	der_csu_dtheta = der_csu_detai * der_etai_dtheta + \
		der_csu_detaq * der_etaq_dtheta + \
		der_csu_detau * der_etau_dtheta + \
		der_csu_detav * der_etav_dtheta + \
		der_csu_drhoq * der_rhoq_dtheta + \
		der_csu_drhou * der_rhou_dtheta + \
		der_csu_drhov * der_rhov_dtheta
	der_csu_dchi = der_csu_detai * der_etai_dchi + \
		der_csu_detaq * der_etaq_dchi + \
		der_csu_detau * der_etau_dchi + \
		der_csu_detav * der_etav_dchi + \
		der_csu_drhoq * der_rhoq_dchi + \
		der_csu_drhou * der_rhou_dchi + \
		der_csu_drhov * der_rhov_dchi
	der_csu_db = der_csu_detai * der_etai_db + \
		der_csu_detaq * der_etaq_db + \
		der_csu_detau * der_etau_db + \
		der_csu_detav * der_etav_db + \
		der_csu_drhoq * der_rhoq_db + \
		der_csu_drhou * der_rhou_db + \
		der_csu_drhov * der_rhov_db
	der_csu_dadamp = der_csu_detai * der_etai_dadamp + \
		der_csu_detaq * der_etaq_dadamp + \
		der_csu_detau * der_etau_dadamp + \
		der_csu_detav * der_etav_dadamp + \
		der_csu_drhoq * der_rhoq_dadamp + \
		der_csu_drhou * der_rhou_dadamp + \
		der_csu_drhov * der_rhov_dadamp
	der_csu_ddeltalambdaD = der_csu_detai * der_etai_ddeltalambdaD + \
		der_csu_detaq * der_etaq_ddeltalambdaD + \
		der_csu_detau * der_etau_ddeltalambdaD + \
		der_csu_detav * der_etav_ddeltalambdaD + \
		der_csu_drhoq * der_rhoq_ddeltalambdaD + \
		der_csu_drhou * der_rhou_ddeltalambdaD + \
		der_csu_drhov * der_rhov_ddeltalambdaD
	der_csu_dvm = der_csu_detai * der_etai_dvm + \
		der_csu_detaq * der_etaq_dvm + \
		der_csu_detau * der_etau_dvm + \
		der_csu_detav * der_etav_dvm + \
		der_csu_drhoq * der_rhoq_dvm + \
		der_csu_drhou * der_rhou_dvm + \
		der_csu_drhov * der_rhov_dvm

	der_csv_deta0 = der_csv_detai * der_etai_deta0 + \
		der_csv_detaq * der_etaq_deta0 + \
		der_csv_detau * der_etau_deta0 + \
		der_csv_detav * der_etav_deta0 + \
		der_csv_drhoq * der_rhoq_deta0 + \
		der_csv_drhou * der_rhou_deta0 + \
		der_csv_drhov * der_rhov_deta0
	der_csv_dtheta = der_csv_detai * der_etai_dtheta + \
		der_csv_detaq * der_etaq_dtheta + \
		der_csv_detau * der_etau_dtheta + \
		der_csv_detav * der_etav_dtheta + \
		der_csv_drhoq * der_rhoq_dtheta + \
		der_csv_drhou * der_rhou_dtheta + \
		der_csv_drhov * der_rhov_dtheta
	der_csv_dchi = der_csv_detai * der_etai_dchi + \
		der_csv_detaq * der_etaq_dchi + \
		der_csv_detau * der_etau_dchi + \
		der_csv_detav * der_etav_dchi + \
		der_csv_drhoq * der_rhoq_dchi + \
		der_csv_drhou * der_rhou_dchi + \
		der_csv_drhov * der_rhov_dchi
	der_csv_db = der_csv_detai * der_etai_db + \
		der_csv_detaq * der_etaq_db + \
		der_csv_detau * der_etau_db + \
		der_csv_detav * der_etav_db + \
		der_csv_drhoq * der_rhoq_db + \
		der_csv_drhou * der_rhou_db + \
		der_csv_drhov * der_rhov_db
	der_csv_dadamp = der_csv_detai * der_etai_dadamp + \
		der_csv_detaq * der_etaq_dadamp + \
		der_csv_detau * der_etau_dadamp + \
		der_csv_detav * der_etav_dadamp + \
		der_csv_drhoq * der_rhoq_dadamp + \
		der_csv_drhou * der_rhou_dadamp + \
		der_csv_drhov * der_rhov_dadamp
	der_csv_ddeltalambdaD = der_csv_detai * der_etai_ddeltalambdaD + \
		der_csv_detaq * der_etaq_ddeltalambdaD + \
		der_csv_detau * der_etau_ddeltalambdaD + \
		der_csv_detav * der_etav_ddeltalambdaD + \
		der_csv_drhoq * der_rhoq_ddeltalambdaD + \
		der_csv_drhou * der_rhou_ddeltalambdaD + \
		der_csv_drhov * der_rhov_ddeltalambdaD
	der_csv_dvm = der_csv_detai * der_etai_dvm + \
		der_csv_detaq * der_etaq_dvm + \
		der_csv_detau * der_etau_dvm + \
		der_csv_detav * der_etav_dvm + \
		der_csv_drhoq * der_rhoq_dvm + \
		der_csv_drhou * der_rhou_dvm + \
		der_csv_drhov * der_rhov_dvm
	#****************************************************************************************************


	delta = etai**2. * (etai**2. - etaq**2. - etau**2. - etav**2. + rhoq**2. + rhou**2. + rhov**2.) - (etaq * rhoq + etau * rhou + etav * rhov)**2.

	stoki = (1. + br * delta**(-1.) * etai * (etai**2. + rhoq**2. + rhou**2. + rhov**2.)) / (1. + br)
	stokq = - br / (1. + br) * delta**(-1.) * (etai**2. * etaq + etai * (etav * rhou - etau * rhov) + rhoq * (etaq * rhoq + etau * rhou + etav * rhov))
	stoku = - br / (1. + br) * delta**(-1.) * (etai**2. * etau + etai * (etaq * rhov - etav * rhoq) + rhou * (etaq * rhoq + etau * rhou + etav * rhov))
	stokv = - br / (1. + br) * delta**(-1.) * (etai**2. * etav + etai * (etau * rhoq - etaq * rhou) + rhov * (etaq * rhoq + etau * rhou + etav * rhov))

	csi = etai * (etai**2. + rhoq**2. + rhou**2. + rhov**2.)
	csq = etai**2. * etaq + etai * (etav * rhou - etau * rhov) + rhoq * (etaq * rhoq + etau * rhou + etav * rhov)
	csu = etai**2. * etau + etai * (etaq * rhov - etav * rhoq) + rhou * (etaq * rhoq + etau * rhou + etav * rhov)
	csv = etai**2. * etav + etai * (etau * rhoq - etaq * rhou) + rhov * (etaq * rhoq + etau * rhou + etav * rhov)
	#****************************************************************************************************
	#DERIVADAS DE LOS PARAMETROS DE STOKES RESPECTO DE DELTA Y CONSTANTES
	der_stokesi_ddelta = -br * csi / (1. + br) / delta**2.
	der_stokesi_dcsi = br / (1. + br) / delta

	der_stokesq_ddelta = br * csq / (1. + br) / delta**2.
	der_stokesq_dcsq = -br / (1. + br) / delta

	der_stokesu_ddelta = br * csu / (1. + br) / delta**2.
	der_stokesu_dcsu = -br / (1. + br) / delta

	der_stokesv_ddelta = br * csv / (1. + br) / delta**2.
	der_stokesv_dcsv = -br / (1. + br) / delta


	#stokes I:
	der_stokesi_dbr = (csi / delta - 1.) / (1. + br)**2.
	der_stokesi_deta0 = der_stokesi_ddelta * der_delta_deta0 + der_stokesi_dcsi * der_csi_deta0
	der_stokesi_dtheta = der_stokesi_ddelta * der_delta_dtheta + der_stokesi_dcsi * der_csi_dtheta
	der_stokesi_dchi = der_stokesi_ddelta * der_delta_dchi + der_stokesi_dcsi * der_csi_dchi
	der_stokesi_db = der_stokesi_ddelta * der_delta_db + der_stokesi_dcsi * der_csi_db
	der_stokesi_dadamp = der_stokesi_ddelta * der_delta_dadamp + der_stokesi_dcsi * der_csi_dadamp
	der_stokesi_ddeltalambdad = der_stokesi_ddelta * der_delta_ddeltalambdaD + der_stokesi_dcsi * der_csi_ddeltalambdaD
	der_stokesi_dvm = der_stokesi_ddelta * der_delta_dvm + der_stokesi_dcsi * der_csi_dvm
	#stokes Q:
	der_stokesq_dbr = -csq / delta / (1. + br)**2.
	der_stokesq_deta0 = der_stokesq_ddelta * der_delta_deta0 + der_stokesq_dcsq * der_csq_deta0
	der_stokesq_dtheta = der_stokesq_ddelta * der_delta_dtheta + der_stokesq_dcsq * der_csq_dtheta
	der_stokesq_dchi = der_stokesq_ddelta * der_delta_dchi + der_stokesq_dcsq * der_csq_dchi
	der_stokesq_db = der_stokesq_ddelta * der_delta_db + der_stokesq_dcsq * der_csq_db
	der_stokesq_dadamp = der_stokesq_ddelta * der_delta_dadamp + der_stokesq_dcsq * der_csq_dadamp
	der_stokesq_ddeltalambdad = der_stokesq_ddelta * der_delta_ddeltalambdaD + der_stokesq_dcsq * der_csq_ddeltalambdaD
	der_stokesq_dvm = der_stokesq_ddelta * der_delta_dvm + der_stokesq_dcsq * der_csq_dvm
	#stokes U:
	der_stokesu_dbr = -csu / delta / (1. + br)**2.
	der_stokesu_deta0 = der_stokesu_ddelta * der_delta_deta0 + der_stokesu_dcsu * der_csu_deta0
	der_stokesu_dtheta = der_stokesu_ddelta * der_delta_dtheta + der_stokesu_dcsu * der_csu_dtheta
	der_stokesu_dchi = der_stokesu_ddelta * der_delta_dchi + der_stokesu_dcsu * der_csu_dchi
	der_stokesu_db = der_stokesu_ddelta * der_delta_db + der_stokesu_dcsu * der_csu_db
	der_stokesu_dadamp = der_stokesu_ddelta * der_delta_dadamp + der_stokesu_dcsu * der_csu_dadamp
	der_stokesu_ddeltalambdad = der_stokesu_ddelta * der_delta_ddeltalambdaD + der_stokesu_dcsu * der_csu_ddeltalambdaD
	der_stokesu_dvm = der_stokesu_ddelta * der_delta_dvm + der_stokesu_dcsu * der_csu_dvm
	#stokes V:
	der_stokesv_dbr = -csv / delta / (1. + br)**2.
	der_stokesv_deta0 = der_stokesv_ddelta * der_delta_deta0 + der_stokesv_dcsv * der_csv_deta0
	der_stokesv_dtheta = der_stokesv_ddelta * der_delta_dtheta + der_stokesv_dcsv * der_csv_dtheta
	der_stokesv_dchi = der_stokesv_ddelta * der_delta_dchi + der_stokesv_dcsv * der_csv_dchi
	der_stokesv_db = der_stokesv_ddelta * der_delta_db + der_stokesv_dcsv * der_csv_db
	der_stokesv_dadamp = der_stokesv_ddelta * der_delta_dadamp + der_stokesv_dcsv * der_csv_dadamp
	der_stokesv_ddeltalambdad = der_stokesv_ddelta * der_delta_ddeltalambdaD + der_stokesv_dcsv * der_csv_ddeltalambdaD
	der_stokesv_dvm = der_stokesv_ddelta * der_delta_dvm + der_stokesv_dcsv * der_csv_dvm

	si = np.vstack([der_stokesi_dbr, der_stokesi_deta0, der_stokesi_dtheta, der_stokesi_dchi, der_stokesi_db, der_stokesi_dadamp, der_stokesi_ddeltalambdad, der_stokesi_dvm])
	sq = np.vstack([der_stokesq_dbr, der_stokesq_deta0, der_stokesq_dtheta, der_stokesq_dchi, der_stokesq_db, der_stokesq_dadamp, der_stokesq_ddeltalambdad, der_stokesq_dvm])
	su = np.vstack([der_stokesu_dbr, der_stokesu_deta0, der_stokesu_dtheta, der_stokesu_dchi, der_stokesu_db, der_stokesu_dadamp, der_stokesu_ddeltalambdad, der_stokesu_dvm])
	sv = np.vstack([der_stokesv_dbr, der_stokesv_deta0, der_stokesv_dtheta, der_stokesv_dchi, der_stokesv_db, der_stokesv_dadamp, der_stokesv_ddeltalambdad, der_stokesv_dvm])

	return si, sq, su, sv

def marquadt(stokIn, stokQn, stokUn, stokVn, Bmag, deltalambdaD, omegam, theta, chi, br, adamp, eta0, lamb, lambv, s1, l1, j1, s2, l2, j2, mn, itmax=10, lambdap=10., noise=1.,Mneg=None,fact=10.):
	nparam = 8l
	itm = 0l
	condit = 1l
	npoint = 4l * len(lambv)

	stokIp, stokQp, stokUp, stokVp = sintetizador(Bmag, deltalambdaD, omegam, theta, chi, br, adamp, eta0, lamb, lambv, s1, l1, j1, s2, l2, j2, mn)

	while (condit == 1):

        	chi2 = (((stokIn - stokIp) / noise)**2. + \
                	((stokQn - stokQp) / noise)**2. + \
                	((stokUn - stokUp) / noise)**2. + \
                	((stokVn - stokVp) / noise)**2.).sum()

        	deris, derqs, derus, dervs = derivador(Bmag, deltalambdaD, omegam, theta, chi, br, adamp, eta0, lamb, lambv, s1, l1, j1, s2, l2, j2, mn)

        	beta = np.zeros((nparam))
        	alphap = np.zeros((nparam, nparam))

		if (Mneg == None):
#AQUI HABRA QUE HACER MODIFICACIONES PARA LA ESTIMACION DE ENTRADA
        		Mneg = np.diag(np.repeat(1./noise**2., npoint))

		s = np.concatenate([(stokIn - stokIp), (stokQn - stokQp), (stokUn - stokUp), (stokVn - stokVp)])
       		mat = np.dot(Mneg,s)

        	for index in range(nparam):
			r = np.concatenate([deris[index,:], derqs[index,:], derus[index,:], dervs[index,:]])
                	beta[index] = np.dot(r, mat)

        	for index1 in range(nparam):
                	for index2 in range(nparam):
				r1 = np.concatenate([deris[index1,:], derqs[index1,:], derus[index1,:], dervs[index1,:]])
				r2 = np.concatenate([deris[index2,:], derqs[index2,:], derus[index2,:], dervs[index2,:]])
				alphap[index1, index2] = np.dot(r1, np.dot(Mneg,r2))
                        	if (index1 == index2):
                                	alphap[index1, index1] = alphap[index1, index1] * (1. + lambdap)

        	#for index1 in range(nparam):
                	#for index2 in range(nparam):
                        	#alphap[index1, index2] = (1. / noise**2. * (deris[index1,:] * deris[index2,:]) + \
                                	#1. / noise**2. * (derqs[index1,:] * derqs[index2,:]) + \
                                	#1. / noise**2. * (derus[index1,:] * derus[index2,:]) + \
                                	#1. / noise**2. * (dervs[index1,:] * dervs[index2,:])).sum()
				#if (index1 == index2):
					#alphap[index1, index1] = alphap[index1, index1] * (1. + lambdap)
                			#beta[index1] = ((stokIn - stokIp) / noise**2. * deris[index1, :] + \
                        			#(stokQn - stokQp) / noise**2. * derqs[index1, :] + \
                        			#(stokUn - stokUp) / noise**2. * derus[index1, :] + \
                        			#(stokVn - stokVp) / noise**2. * dervs[index1, :]).sum()

        	perturb, _, _, _ = np.linalg.lstsq(alphap, beta)#, rcond = 1.e-5)

        	sBmag = np.abs(Bmag + perturb[4])
        	sadamp = adamp + perturb[5]
        	somegam = omegam + perturb[7]
        	sdeltalambdaD = deltalambdaD + perturb[6]
        	stheta = theta + perturb[2]
        	schi = chi + perturb[3]
        	seta0 = eta0 + perturb[1]
        	sbr = br + perturb[0]

        	stokIs, stokQs, stokUs, stokVs = sintetizador(sBmag, sdeltalambdaD, somegam, stheta, schi, sbr, sadamp, seta0, lamb, lambv, s1, l1, j1, s2, l2, j2, mn)

        	chi2s = (((stokIn - stokIs) / noise)**2. + \
                	((stokQn - stokQs) / noise)**2. + \
                	((stokUn - stokUs) / noise)**2. + \
                	((stokVn - stokVs) / noise)**2.).sum()

        	if (chi2s > chi2):
                	lambdap = lambdap * fact
        	else:
                	Bmag = sBmag
                	adamp = sadamp
                	omegam = somegam
                	deltalambdaD = sdeltalambdaD
                	theta = stheta
                	chi = schi
                	eta0 = seta0
                	br = sbr
                	lambdap = lambdap / fact
                	if (lambdap < 1.e-5): lambdap = 1.e-5

        	stokIp = stokIs.copy()
        	stokQp = stokQs.copy()
        	stokUp = stokUs.copy()
        	stokVp = stokVs.copy()

        	itm += 1
		
        	if (itm > itmax): condit = 0

	return Bmag, omegam, adamp, br, theta * 180. / np.pi, chi * 180. / np.pi, eta0, deltalambdaD

