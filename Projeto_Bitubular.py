#Comentário: Nesse projeto as propriedades referentes a região anular serão representadas
#por um subscrito "o", enquanto a região interna terá um subscrito "i". E os subscitos "h" e "f",
#representam, frio e quente, respectivamente. As propriedades referentes ao material dos tubos serão representadas
#por letra maiuscula


#Bibliotecas usadas no projeto
import CoolProp.CoolProp as CP #Para obter as propriedades dos fluidos
import math as mt

#Parametros de Entrada (ENTRADA)
Fluido_i = str(input('Qual é o fluído que passa na região interna:    '))
Fluido_o = str(input('Qual é o fluído que passa na região anular:    '))
Th_i = float(input('Qual a MAIOR temperatura do fluido na região interna [°C]:  '))
Tc_i = float(input('Qual a MENOR temperatura do fluido na região interna [°C]: '))
Th_o = float(input('Qual a MAIOR temperatura do fluido na região anular [°C]:  '))
Tc_o = float(input('Qual a MENOR temperatura do fluido na região anular [°C]: '))
di = float(input('Qual o diâmetro interno do tubo interno [m]:   ')) #USAR JAYA
do = float(input('Qual o diâmetro externo do tubo interno [m]:   ')) #USAR JAYA
Do = float(input('Qual o diâmetro interno do tubo externo [m]:   ')) #USAR JAYA
m_o = float(input('Qual o fluxo de massa do fluido na região anular [kg/s]: ')) #PODE SER ALTERADO
l = 3.5 #Definir

#Material do tubo interno
K = 54 #W/mK

#Temperatura de Bulk
Tb_i = ((Th_i + Tc_i)/2)+273 #Em K
Tb_o = ((Th_o + Tc_o)/2)+273 #Em K

#Determinando propriedades dos fluidos
#Parte interna
Q = 0 #Liquido Saturadp pois é monofásico
dens_i = CP.PropsSI('D', 'T', Tb_i, 'Q', Q, Fluido_i)      # Densidade [kg/m³]
visc_i = CP.PropsSI('V', 'T',Tb_i , 'Q', Q, Fluido_i)       # Viscosidade [Pa·s]
k_i = CP.PropsSI('L', 'T', Tb_i, 'Q', Q, Fluido_i)        # Condutividade térmica [W/m·K]
Cp_i = CP.PropsSI('C', 'T', Tb_i, 'Q', Q, Fluido_i)       # Calor específico [J/kg·K]
Pdt_i = CP.PropsSI('Prandtl', 'T', Tb_i, 'Q', Q, Fluido_i) # Número de Prandtl [-]

#ii) parte anular
Q = 0
dens_o = CP.PropsSI('D', 'T', Tb_o, 'Q', Q, Fluido_o)      # Densidade [kg/m³]
visc_o =CP.PropsSI('V', 'T', Tb_o, 'Q', Q, Fluido_o)       # Viscosidade [Pa·s]
k_o =CP.PropsSI('L', 'T', Tb_o, 'Q', Q, Fluido_o)        # Condutividade térmica [W/m·K]
Cp_o = CP.PropsSI('C', 'T', Tb_o, 'Q', Q, Fluido_o)    # Calor específico [J/kg·K]
Pdt_o =CP.PropsSI('Prandtl', 'T', Tb_o, 'Q', Q, Fluido_o) # Número de Prandtl [-]

#Determinando Área Interna e Anular, Diametro Hidraulico, Diametro Equivalente
A_i = (mt.pi/4)*(di**2) #m2
A_o = (mt.pi/4)*((Do**2)-(do**2))
Dh_i = di #m
Dh_o = Do-do #m
De_i = di #m
De_o = ((Do**2)-(do**2))/do #m

#DETERMINANDO FLUXO DE MASSA E VELOCIDADES MÉDIA DO ESCOAMENTO (AJUSTAR PRA LEVAR EM CONSIDERAÇÃO TROCA DO ESCOAMENTO)
m_i = (m_o*Cp_o*(Th_o-Tc_o))/(Cp_i*(Th_i-Tc_i)) #Balanço de Energia
u_i = m_i/(dens_i*A_i) #m/s
u_o = m_o/(dens_o*A_o) #m/s

#REYNOLDS
Re_i = (4*m_i)/(mt.pi*visc_i*Dh_i)
Re_o = (dens_o*u_o*Dh_o)/visc_o

#INCRUSTAÇÕES
Rf_o = 0.000352
Rf_i = 0.000176

#NUSSELT PARA A PARTE INTERNA (COLOCAR OUTRAS CORRELAÇÕES)
if 2300 < Re_i < 10**4:  #Gnielinski (TURBULENTO: MAS APLICAVEL NA REGIÃO DE TRANSIÇÃO)
    f_i = (1.58 * mt.log(Re_i, mt.e) - 3.28) ** -2
    Nu_i = ((f_i / 2)*(Re_i-1000)*Pdt_i)/(1+12.7*mt.sqrt(f_i/2)*(Pdt_i**(2/3)-1))

elif 10**4 < Re_i < 5*10**6 and 0.5<Pdt_i<2000: #PETUKHOV (TURBULENTO)
    f_i = (1.58*mt.log(Re_i,mt.e)-3.28)**-2
    Nu_i = ((f_i/2)*Re_i*Pdt_i)/(1.07+12.7*mt.sqrt(f_i/2)*(Pdt_i**(2/3)-1))

#NUSSELT PARA A PARTE ANULAR
if 2300 < Re_o < 10**4:  #Gnielinski (TURBULENTO: MAS APLICAVEL NA REGIÃO DE TRANSIÇÃO)
    f_o = (1.58 * mt.log(Re_o, mt.e) - 3.28) ** -2
    Nu_i = ((f_o / 2)*(Re_o-1000)*Pdt_o)/(1+12.7*mt.sqrt(f_o/2)*(Pdt_o**(2/3)-1))

elif 10**4 < Re_o < 5*10**6 and 0.5<Pdt_o <2000: #PETUKHOV (TURBULENTO)
    f_o = (1.58*mt.log(Re_o,mt.e)-3.28)**-2
    Nu_o = ((f_o / 2) * Re_o * Pdt_o) / (1.07 + 12.7 * mt.sqrt(f_o / 2) * (Pdt_o ** (2 / 3) - 1))

#COEFICIENTES DE CONVECÇÃO
h_i = (Nu_i*k_i)/De_i
h_o = (Nu_o*k_o)/De_o

#Coeficiente Global de transferencia de calor com incrustações
Resistencias = (do/(h_i*di))+(do*Rf_i/(di))+(1/h_o)+(Rf_o)+((do*mt.log(do/di,mt.e))/(2*K)) #Em função de do
U_o_f = 1/Resistencias

#Diferença de Temperatura Média Logaritmica (DTML)
DeltaT1, DeltaT2 = (Th_i-Th_o), (Tc_i-Tc_o)
if DeltaT1 == DeltaT2:
    DTML = (DeltaT1 + DeltaT2)/2
else:
    DTML = (DeltaT1-DeltaT2)/(mt.log(DeltaT1/DeltaT2,mt.e))

#Taxa de troca de Calor
Q = m_i*Cp_i*(Th_i-Tc_i)

#Area
A_o = Q/(U_o_f*DTML)

#Número de Grampos
A_gp = 2*mt.pi*l*do
N = A_o/A_gp #número de grampos
N = round(N) #Arredonda

#Coeficiente Global de transferencia de calor sem incrustações
Resistencias = (do/(h_i*di))+(1/h_o)+((do*mt.log(do/di,mt.e))/(2*K)) #Em função de do
U_o_c = 1/Resistencias

#FATOR DE LIMPEZA
CF = U_o_f/U_o_c

#Porcentagem sobre a superficie OS
Rft = (1-CF)/(U_o_c*CF)
OS = 100*U_o_c*Rft

#QUEDA DE PRESSÃO NA PARTE INTERNA E ANULAR
f_i = (1.58*mt.log(Re_i,mt.e)-3.28)**-2
f_o = (1.58*mt.log(Re_o,mt.e)-3.28)**-2

Deltapi = 4*f_i*((2*l)/Dh_i)*N*dens_i*(u_i**2)/2
Deltapo = 4*f_o*((2*l)/Dh_o)*N*dens_o*(u_o**2)/2


#POTÊNCIA DA BOMBA EXIGIDA
Nb = 0.8 #Eficiencia de bomba
Pi = (Deltapi*m_i)/(Nb*dens_i)
Po = (Deltapo*m_o)/(Nb*dens_o)