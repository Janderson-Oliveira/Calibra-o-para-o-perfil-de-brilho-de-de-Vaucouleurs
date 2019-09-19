#!/home/j/anaconda2/bin/python


'''
    Setembro 2019
    
    Perfil de brilho de uma galaxia
    
    @author:  Janderson Oliveira

    Simple usage example: 

./atv_1.py --priorsfile=NGC0426_SDSSr.txt --titulo=NGC0426 --mag_ref=11.02 --cts=4752. --galaxia_mag_ref=12.99


./atv_1.py --priorsfile=NGC0430_SDSSr.txt --titulo=NGC0430 --mag_ref=11.02 --cts=4752. --galaxia_mag_ref=12.60

'''


import sys

from optparse import OptionParser

from uncertainties.umath import * 

from uncertainties import ufloat

import uncertainties.umath as umath

from scipy.stats import linregress

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from matplotlib import pyplot as pl


import time

import numpy as np

import math

from math import sqrt


parser = OptionParser()
parser.add_option("-i", "--priorsfile", dest="priorsfile", help='Prior parameters file',type='string',default="")
parser.add_option("-t", "--titulo", dest="titulo", help='Titulo dos graficos',type='string',default="Ajuste")
parser.add_option("-m", "--mag_ref", dest="mag_ref", help='Mag estrela de referencia',type='string',default="")
parser.add_option("-M", "--galaxia_mag_ref", dest="galaxia_mag_ref", help='Mag estrela de referencia',type='string',default="")
parser.add_option("-c", "--cts", dest="cts", help='Contagens obtidas pelas isofotas',type='string',default="")
parser.add_option("-v",action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with eclipsefit.py -h ";sys.exit(1);

if options.verbose:
    print 'Prior parameters file: ', options.priorsfile
    print 'Title: ', options.titulo


t0 = time.time()
t0_str = time.ctime()


# load data

SMA, RSMA, INTENS, RMS, MAG, MAG_LERR, MAG_UERR = np.loadtxt(options.priorsfile, unpack=True, comments='#')

titulo = "Ajuste para " + str(options.titulo)

mag_ref = float(options.mag_ref)
cts_ref = float(options.cts)
galaxia_mag_ref = float(options.galaxia_mag_ref)



################# biblioteca de erros #####################

def convert_ufloat(var, error):
    list_var_ufloat = []
    for i in range(len(var)):
        var_ufloat = ufloat(var[i], error[i])
        list_var_ufloat.append(var_ufloat)
    return list_var_ufloat


def vec_nominal_value(vec):
    nominal = []
    for i in range(len(vec)):
        value = vec[i].n
        nominal.append(value)
    return np.array(nominal)

def vec_erro_value(vec_err):
    erro = []
    for i in range(len(vec_err)):
        value = vec_err[i].s
        erro.append(value)
    return np.array(erro)



##########################################################


# utilizar valores com erros 
Vec_MAG = np.array(MAG)
Vec_MAG_LERR = np.array(MAG_LERR)

MAG_err_ufloat = convert_ufloat(Vec_MAG, Vec_MAG_LERR)


def convert_mag(mag_ref, cts_ref, cts):

    #mag_ref = mag da estrela de referencia
    #cts = valor das contagens do objeto obtido pelo iraf
    #cts_ref = valor das contagens da estrela de referencia obtida pelo iraf
    
    pix  = 0.39597 #arcsec
    A = log(pix**2)
    B = cts/cts_ref

    #mag_arcsec = mag_ref - 2.5*log(B) - 2.5*A

    mag_arcsec = mag_ref + cts  + 2.5*log(cts_ref) - 2.5*A
   
    mag_arcsec_2 = 22.5 - 2.5*log(cts)

    return mag_arcsec, mag_arcsec_2


mag_arcsec_list = []
mag_arcsec_list_2 = []



for i in range(len(MAG)):

    mag_arcsec, mag_arcsec_2 = convert_mag(mag_ref, cts_ref, MAG_err_ufloat[i])
    #print mag_arcsec, ' ', mag_arcsec_2
    print mag_arcsec
    mag_arcsec_list.append(mag_arcsec)
    mag_arcsec_list_2.append(mag_arcsec_2)


'''
for i in range(len(mag_arcsec_list)):
    print mag_arcsec_erro_list[i], mag_arcsec_erro_list_2[i]
'''






def plot_(x, y,  y_err, titulo):
    list_ = []

    for i in range(len(x)):
        list_.append(0.)
 
    plt.errorbar(x, y,  y_err, label=u"Dados")
    #plt.errorbar(x, y2,  y2_err, label=u"Dados")
    #pl.plot(Potencial,list_, '-')

    plt.title(titulo) 
    plt.grid(True)
    pl.plot(x,y, 'o')

    ay = plt.gca()
    ay.invert_yaxis()

    #pl.plot(x,y2, 'o')
    plt.xlabel('Raio')
    plt.ylabel('Magnitude')
    plt.legend()
    plt.legend(loc=2)
    plt.savefig(titulo+"_sem conversao")
    plt.show()


x = RSMA

y = vec_nominal_value(mag_arcsec_list)
y_err = vec_erro_value(mag_arcsec_list)

#y = vec_nominal_value(mag_arcsec_list_2)
#y_err = vec_erro_value(mag_arcsec_list_2)


pl_ = plot_(x, y, y_err, titulo)


def determine_regressao(vec_x, vec_y, vecy_err, titulo):


    x = vec_x
    y = vec_y
 
 
    m, b, R, p, SEm = linregress(x, y)

    #m> declive
    #b> ordenada na origem
    #R> coeficiente de correlacao (de Pearson)
    #p> p-value do teste F em que H0: y = const, independente de x
    #SEm> erro padrao do declive


    def lin_regression(x, y):
        """Simple linear regression (y = m * x + b + error)."""
        m, b, R, p, SEm = linregress(x, y)

        # need to compute SEb, linregress only computes SEm
        n = len(x)
        SSx = np.var(x, ddof=1) * (n-1)  # this is sum( (x - mean(x))**2 )
        SEb2 = SEm**2 * (SSx/n + np.mean(x)**2)
        SEb = SEb2**0.5

        return m, b, SEm, SEb, R, p

    m, b, Sm, Sb, R, p = lin_regression(x, y)

    m0_e = ufloat(b, Sb)
    b_e = ufloat(m, Sm)
    #v_x0 = (-b_e)/m_e

    print('b>>b_err = {:>.9g} +- {:6.9f}'.format(m, Sm))
    print('m0>>b_err = {:>.11g} +- {:6.9f}\n'.format(b, Sb))
    #print ('v para c_0 = {:6.9g}'.format(v_x0))
    print ''

    #print('R2 = {:7.5f}'.format(R**2))
    #print('p of test F : {:<8.6f}'.format(p))
 
    plt.subplot(211)
    plt.errorbar(x, y, y_err, linestyle='-.',marker="o", label=u"Dados")
    x2 = np.array([min(x), max(x)])
    pl.plot(x2, m * x2 + b, linestyle='-', color='g', label=u"Regressao")

    # Anotacao sobre o grafico:
    #ptxt = 'm = {:>.4g} +- {:6.4f}\nb = {:>.4g} +- {:6.4f}\nR2 = {:7.5f}'
    #plt.text(4.5, 4, ptxt.format(m, Sm, b, Sb, R**2), fontsize=14)

    ay = plt.gca()
    ay.invert_yaxis()


    plt.grid(True)
    plt.xlabel('Raio')
    plt.ylabel('Magnitude')
    #plt.ylabel(r'$Y(m)$   e   $V_y(\frac{m}{s})$')
    plt.title(titulo)
    plt.legend()
    plt.legend(loc=2)
  

    modelo  = m * x + b
    plt.subplot(212)
    residuos = y - modelo
    plt.plot(y, residuos, label='Residuos', lw=0.5, marker='.', ms=2, drawstyle='default')
    #plt.plot(y, residuos)
    plt.grid(True)
    plt.xlabel('Magnitude pelo raio',fontsize=16)
    plt.ylabel('Residuos',fontsize=16)
    #plt.xlim((x[0],x[-1]))

    plt.savefig(titulo)
    pl.show()

   

    per = ufloat(m, Sm)
    t0 = ufloat(b, Sb)
    return per, t0

b, m0 = determine_regressao(x, y, y_err, titulo)

#print b, m0
# calculo do raio efetivo e da luminosidade

def calculo_raio_lum(b,m0, mag_ref):
    print "b=", b
    print "m0=", m0
    constante = 1.0857362
    r = ( constante / b )**4

    f=40320.
    i0 = exp(m0/(-2.5)) 
    l  = f*3.1415*i0*(mag_ref**2)
        
    return r, l


raio, l = calculo_raio_lum(b,m0, galaxia_mag_ref)
print "raio=", raio
print "luminosidade=", l 


tf = time.time()
tf_str = time.ctime()
dt = tf - t0
print "Inicio:",t0_str," Fim:",tf_str
print "Tempo total: ", dt, "s"
