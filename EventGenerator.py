#!/usr/bin/env python3


import random
import numpy as np
import csv
from scipy import constants
from scipy import optimize
from scipy import interpolate
import unum
from unum import IncompatibleUnitsError
from scipy.integrate import simps

def estado_inicial(eP, m1, m2):
    
    """
	Esta funcao e responsavel por gerar o quadrimomento das particulas que irao ser colididas.
	:param eP: Energia do proton utilizado na colisao.
	:param m1: Massa da particula 1.
	:param m2: Massa da particula 2.
	:return: Os quadrimomentos das particulas 1 e 2.
    """
        
#Calculando a energia e o quadrimomento da particula 1:
    e1 = eP*gerador_prob()
    p_Parton1 = np.array([e1,0,0,np.sqrt(e1**2-m1**2)])
    
    
#Calculando a energia e o quadrimomento da particula 2:    
    e2 = eP*gerador_prob()
    p_Parton2 = np.array([e2,0,0,-np.sqrt(e2**2-m2**2)])
    
    
    return p_Parton1,p_Parton2



def colisao(p1,p2,mX,mY):
    
    """
    Esta funcao realiza a colisao das particulas 1 e 2.
    :param p1: Quadrimomento da particula 1.
    :param p2: Quadrimomento da particula 2.
    :param mX: Massa da particula X (Boson).
    :param my: Massa da particula Y.
    :return: Quadrimomento das particulas X e Y.
    """
    
    
    #Velocidade do boost 
    v = velocidade_boost(p1,p2)
       
    #Encontrando os quadrimomentos no referencial de centro de massa do sistema:
    p1_cm = transformacaoDeLorentz(p1,v)
    p2_cm = transformacaoDeLorentz(p2,v)
    
    #Quadrimomento total no centro de massa:
    pT_cm = p1_cm+p2_cm
        
    
    pT_cm[3] = np.around(pT_cm[3], decimals = 4)
    
#Depois da colisao:
    #Gerando os angulos que a particula X foi emitida no referencial do centro de massa:
    theta = np.arccos(1 - (2 * np.random.random()))
    phi = random.uniform(0,2*np.pi)
    
    #Calculando a energia da particula X:
    eX = ((pT_cm[0]**2)+(mX**2)-(mY**2))/(2*pT_cm[0])
    
    
    #Calculando o modulo do momento tridimensional de X:
    if ((eX**2-mX**2) < 0):
        eX = np.around(eX, decimals = 4)
        
    mod_pX = np.sqrt(eX**2-mX**2)
    
    
    #Calculando o quadrimomento de X:
    pX_cm = np.array([eX,mod_pX*(np.cos(phi))*(np.sin(theta)),mod_pX*(np.sin(phi))*(np.sin(theta)),mod_pX*(np.cos(theta))])
    
    #Calculando a energia da particula Y:
    eY = ((pT_cm[0]**2)+(mY**2)-(mX**2))/(2*pT_cm[0])
    
    #Calculando o quadrimomento de Y:
    pY_cm = np.array([eY,-mod_pX*(np.cos(phi))*(np.sin(theta)),-mod_pX*(np.sin(phi))*(np.sin(theta)),-mod_pX*(np.cos(theta))])
    
    #Retornando ao referencial do laboratorio:
    pX = transformacaoDeLorentz(pX_cm,-v)
    
    pY = transformacaoDeLorentz(pY_cm,-v)
    
  
    return pX,pY



def decaimento(mX, mA, mB, pX):
    
    """
    Esta fucao gera as componentes do quadrimomento das particulas A e B.
    :param mX: Massa da particula X (GeV).
    :param mA: Massa da particula A (GeV).
    :param mB: Massa da particula B (GeV).
    :param pX: Quadrimomento da particula mae.
    :param v: Vetor com as componentes da velocidade da particula mae (m/s).
    :param theta: Angulo que a particula A foi emitida em relacao ao eixo z (rad).
    :param phi: Angulo que a particula A foi emitida em relacao ao eixo x (rad).
    :return: Vetor contendo as componentes do quadrimomento das particulas A e B.
    """
    
    #Vetor velocidade da particula mae:
    v = np.array([pX[1],pX[2],pX[3]])/pX[0]
    
    #Gerando os angulos em que a particula A foi emitida
    theta = np.arccos(1 - 2 * np.random.random())
    phi = random.uniform(0,2*np.pi)

    #Calculando a Energia de A e B no referencial do CM de X
    eA = ((mX**2)+(mA**2)-(mB**2))/(2*mX)
    eB = ((mX**2)+(mB**2)-(mA**2))/(2*mX)

    #Calculando o modulo do 3-momento de A e B no referencial do CM de X
    pAB = ((eA**2)-(mA**2))**0.5

    #Calculando as componentes do 4-momento de A no referencial do CM de X
    pA = [eA,pAB*(np.cos(phi))*(np.sin(theta)),pAB*(np.sin(phi))*(np.sin(theta)),pAB*(np.cos(theta))]
    
    #Calculando as componentes do 4-momento de A no referencial do CM de X
    pB = [eB,(-1)*pA[1],(-1)*pA[2],(-1)*pA[3]]
  
    
    #Fazendo o Boost em A
    pA_Lab = transformacaoDeLorentz(pA,v)
    
    #Fazendo o Boost em B
    pB_Lab = transformacaoDeLorentz(pB,v)
        
       
    return pA_Lab,pB_Lab



def detector(pD,desv_pad = 0.1):
    
    """
    Esta funcao simula o efeito detector.
    :param pD: Quadrimomento medido pelo detector.
    :param alpha: Fator de correcao.
    :return: O valor real do quadrimomento.
    """
    
#Sortear com numpy ou utilizar acceptance-rejection?
    alpha  = np.random.normal(1,desv_pad)
#Valor com a correcao do efeito detector
    return alpha*pD
    
	
	
	
def probabilidade_up(x,x_pdf,y_pdf,const):
    
    """
    Esta funcao e a densidade de probabilidade de o parton ter uma determinada fracao da energia do proton.
    :param x: Variavel no intervalo [0,1] que denota a fracao de energia do parton.
    :param x_pdf: Valores do dominio da pdd provenientes da tabela.
    :param y_pdf: Valores da imagem da pdf provenientes da tabela.
    :param a: Limite inferior do dominio da funcao.
    :param b: Limite superior do dominio da funcao.
    :return: O valor da funcao densidade no ponto x
    """
    
    #pdf
    f = interpolate.interp1d(x_pdf, y_pdf)    

    #Limites do dominio:
    a = 0.001
    b = 1
    if (x >= a and x <= b):
        return const*f(x)
    else:
        return 0

	

def gerador_prob():
    
    """
    Esta funcao sorteia numeros aleatorios utilizando o metodo acceptance-rejection.
    :param g(x): Funcao maior que p(x) para todo x.
    :param u: Numero aleatorio distribuido uniformemente entre 0 e 1.
    :param a: Numero aleatorio pertencendo ao dominio de g(x) e p(x).
    :return: Numero aleatorio que se encontra em p(x) ou abaixo dela.
    """
	
	#Criando parametros para a funcao interpolate:
	
    pdf_table = np.loadtxt('upPDF.dat',delimiter=',',skiprows=1) 
    pdf_table = np.array(pdf_table)
    x_pdf = []
    y_pdf = []

    for i in pdf_table:
    	x_pdf.append(i[0])
    	y_pdf.append(i[1])
		
    x_pdf = np.array(x_pdf)
    y_pdf = np.array(y_pdf)
	
    f = interpolate.interp1d(x_pdf, y_pdf)


	#Constante de normalizacao
    domain = np.arange(0.001,1.0,0.01)
    const = 1/simps(f(domain), domain)
    
    #g = optimize.fmin(lambda x: (-1)*probabilidade_up(x,x_pdf,y_pdf), 1, disp = 0)
    g = probabilidade_up(0.001,x_pdf,y_pdf,const)
    u = np.random.random()
    a = np.random.uniform(0,1)
    if (u*g <= probabilidade_up(a,x_pdf,y_pdf,const)):
        return a
    else:
        return gerador_prob()
	
	
	
def transformacaoDeLorentz(quadrivetor,v):
   
    """
    Esta funcao recebe um quadrivetor e aplica as transformacoes de Lorentz nele a fim de representar esse quadrivetor em outro referencial.
    :param quadrivetor: Vetor com as quatro componentes do quadrivetor a ser transformado.
    :param v: Vetor com as componentes da velocidade com que o referencial se move (m/s).
    :param vel: Norma do vetor v.
    :param gamma: Fator de Lorentz.
    :param matriz_Lambda: Forma matricial das Transformacoes de Lorentz.
    :return: Vetor com as componentes do quadrivetor no novo referencial.
    """
    
    #Calculando o modulo da velocidade
    vel = np.linalg.norm(v)
    if (vel == 0):
        return quadrivetor
    
    #Calculando o fator de Lorentz
    gamma = (1/(1-(vel**2)))**0.5
    
    #Definindo a matriz da transformacao
    
    matriz_Lambda = np.zeros((4,4))
    matriz_Lambda[0][0] = gamma
    for i in range(1,4):
        matriz_Lambda[i][0] = matriz_Lambda[0][i] = -gamma*v[i-1]
    
    m3x3 = ((gamma-1)/vel**2)*np.einsum('i,j->ij',v,v) + np.identity(3)
    m3x3 = np.insert(m3x3,0,[0,0,0], axis = 1)
    m3x3 = np.insert(m3x3,0,[0,0,0,0], axis = 0)
    
    matriz_Lambda = matriz_Lambda + m3x3
    
    #Definindo a matriz a ser retornada
    quadrivetor_Lab = []
    
    #Produto de matrizes
    quadrivetor_Lab = np.matmul(matriz_Lambda, quadrivetor)
    
    return quadrivetor_Lab



def velocidade_boost(p1,p2=np.array([0,0,0,0])):
    
    """
   Esta funcao encontra a velocidade do boost para o centro de massa de um sistema de ate duas particulas.
    :param p1: Quadrimomento da particula 1.
    :param p2: Quadrimomento da particula 2.
    :param pT: Quadrimomento total do sistema.
    :return: Vetor contendo as 3 componentes da velocidade do boost
    """
    
    pT = np.array([p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2],p1[3]+p2[3]])
    v_boost = np.array([pT[1],pT[2],pT[3]])/pT[0]
    
    return v_boost



def check_E(p1,p2,mX,mY):
    v = velocidade_boost(p1,p2)
       
    #Encontrando os quadrimomentos no referencial de centro de massa do sistema:
    p1_cm = transformacaoDeLorentz(p1,v)
    p2_cm = transformacaoDeLorentz(p2,v)
    
    #Quadrimomento total no centro de massa:
    pT_cm = p1_cm+p2_cm
        
    if (pT_cm[0] < mX+mY):
        return False
    else:
        return True
	
	
	
def main(n,eP,mX,mY,m1,m2,mA,mB):
	


	#Declaracao das unidades de medida

	unum.Unum.reset()
	unum.Unum.VALUE_FORMAT = "%0.2E" 
	unum.Unum.UNIT_HIDE_EMPTY = True

	#Unidades de massa e energia
	eV = unum.Unum.unit('eV')
	keV = unum.Unum.unit('keV', 10 ** 3 * eV)
	MeV = unum.Unum.unit('MeV', 10 ** 6 * eV)
	GeV = unum.Unum.unit('GeV', 10 ** 9 * eV)
	TeV = unum.Unum.unit('TeV', 10 ** 12 * eV)
	kg = unum.Unum.unit('kg', (constants.value('speed of light in vacuum') ** 2) * (6.2415 * 10**18) * eV)
	g = unum.Unum.unit('g', (constants.value('speed of light in vacuum') ** 2) * (6.2415 * 10**21) * eV)
	J = unum.Unum.unit('J', (6.2415 * 10**18) * eV)

	#Unidades de distancia
	m = unum.Unum.unit('m')

	#Unidade de tempo
	s = unum.Unum.unit('s')

	#Unidades de velocidade
	c = unum.Unum.unit('c', constants.value('speed of light in vacuum') * (m/s))
	
	#Dados de entrada

	eP = eP * TeV #Energia dos protons. (kg, g, eV, keV, MeV, GeV,TeV, J)
	eP = eP.asNumber(GeV)

	m1 = m1 * MeV #Massa do parton 1. (kg, g, eV, keV, MeV, GeV,TeV, J)
	m1 = m1.asNumber(GeV)

	m2 = m2 * MeV #Massa do parton 2. (kg, g, eV, keV, MeV, GeV,TeV, J)
	m2 = m2.asNumber(GeV)


	mX = mX * GeV #Massa da particula X. (kg, g, eV, keV, MeV, GeV,TeV, J)
	mX = mX.asNumber(GeV)

	mY= mY * MeV #Massa da particula Y. (kg, g, eV, keV, MeV, GeV,TeV, J)
	mY= mY.asNumber(GeV)

	mA = mA * GeV #Massa da particula A. (kg, g, eV, keV, MeV, GeV,TeV, J)
	mA = mA.asNumber(GeV)

	mB = mB * GeV #Massa da particula B. (kg, g, eV, keV, MeV, GeV,TeV, J)
	mB = mB.asNumber(GeV)


	n = n #Numero de eventos a serem gerados.
	
	#Criando o arquivo onde os dados serao gravados
	f = open('Quadrimomentos.csv', 'w')
	writer = csv.writer(f)
	writer.writerow( ('p0A','p1A','p2A','p3A','p0B','p1B','p2B','p3B','p0X','p1X','p2X','p3X','p0Y','p1Y','p2Y','p3Y','p01','p11','p21','p31','p02','p12','p22','p32') )
	f.close()
	
	#Gerando os eventos

	for i in range(n):
	#Bloco estado-inicial:
		p1, p2 = estado_inicial(eP, m1, m2)

		#Verificando se a energia e maior ou igual a minima:
		while (check_E(p1,p2,mX,mY) == False):
			p1, p2 = estado_inicial(eP, m1, m2)


	#Bloco colisao:
		pX, pY = colisao(p1, p2, mX, mY)

	#Bloco decaimento:
		pA_Lab, pB_Lab = decaimento(mX, mA, mB, pX)

	#Bloco detector:
		#pY = detector(pY)
		pA_Lab = detector(pA_Lab)
		pB_Lab = detector(pB_Lab)

	#Escrevendo no arquivo csv    
		f = open('Quadrimomentos.csv', 'a')
		writer = csv.writer(f)
		writer.writerow( (pA_Lab[0],pA_Lab[1],pA_Lab[2],pA_Lab[3],pB_Lab[0],pB_Lab[1],pB_Lab[2],pB_Lab[3],pX[0],pX[1],pX[2],pX[3],pY[0],pY[1],pY[2],pY[3],p1[0],p1[1],p1[2],p1[3],p2[0],p2[1],p2[2],p2[3]) )
		f.close()
		
		
		
if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the event generator." )
    ap.add_argument('-n','--nevents',default=1000, type = int,
				   help ='Number of events to be generated')
    ap.add_argument('-e', '--eP', default=6.5,type = float,
            help='half of center of mass energy in TeV')
    ap.add_argument('-x', '--mX', default=125.5, type = float,
            help='Mass of particle X in GeV')
    ap.add_argument('-p', '--m1', default=2.4, type = float,
            help='Mass of parton 1 in MeV')
    ap.add_argument('-q', '--m2', default=2.4, type = float,
            help='Mass of parton 2 in MeV')
    ap.add_argument('-y', '--mY', default=2.4, type = float,
            help='Mass of particle Y in MeV')
    ap.add_argument('-a', '--mA', default=0.0, type = float,
            help='Mass of particle A in GeV')
    ap.add_argument('-b', '--mB', default=0.0, type = float,
            help='Mass of particle B in GeV')
    


    

    args = ap.parse_args()
main(args.nevents,args.eP,args.mX,args.mY,args.m1,args.m2,args.mA,args.mB)
