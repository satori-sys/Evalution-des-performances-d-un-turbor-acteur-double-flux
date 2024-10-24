# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 21:08:45 2024

"""

#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
from thermo import *
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import quad
#
#---------------------------------------------------------------------------------------#		
# Programme principal
#---------------------------------------------------------------------------------------#
print("Sample code for UPMC M2 5AE10 projet. V0")
#
species_names = ["C10H22","O2","CO2","H2O","N2"]
ne = len(species_names)
#---------------------------------------------------------------------------------------#		
# Construction des instances de la classe ThermoSpecies()
#---------------------------------------------------------------------------------------#
species = []
for spec in species_names:
  ie = species_names.index(spec)
  species.append(ThermoSpecies())	
  species[ie].name = spec
  species[ie].num  = ie
#	
#---------------------------------------------------------------------------------------#		
# Lecture du fichier xml
#---------------------------------------------------------------------------------------#
thermofn = "thermo.xml"
thermoreadxml(species, thermofn)
#
#---------------------------------------------------------------------------------------#		
# Vérification des données thermodynamiques d'espèces lues 
#---------------------------------------------------------------------------------------#
for spec in species:
  print("Species name: ", spec.name)
  print("mm : ", spec.mm, "kg/mol","/ Href : ", spec.Href, "J/mol")
#		
#---------------------------------------------------------------------------------------#		
# Calcul de données thermodynamiques dérivées (exemple pour deux espèces, O2 et N2)
#---------------------------------------------------------------------------------------#
iO2 = species_names.index("O2")
print("O2 Cp for T = 300 K: ", species[iO2].Cp(300), "J/mol/K")
print("O2 Cp for T = 1400 K: ", species[iO2].Cp(1400), "J/mol/K")
print("O2 cp for T = 300 K: ", species[iO2].cp(300), "J/kg/K")
print("O2 cp for T = 1400 K: ", species[iO2].cp(1400), "J/kg/K")
#
iN2 = species_names.index("N2") 
print("N2 Cp for T = 300 K: ", species[iN2].Cp(300), "J/mol/K")
print("N2 Cp for T = 1400 K: ", species[iN2].Cp(300), "J/mol/K")
print("N2 cp for T = 300 K: ", species[iN2].cp(300), "J/kg/K")
print("N2 cp for T = 1400 K: ", species[iN2].cp(300), "J/kg/K")
#---------------------------------------------------------------------------------------#		
# Calcul de données thermodynamiques dérivées 
#---------------------------------------------------------------------------------------#
iCO2 = species_names.index("CO2")
print("CO2 Cp for T = 300 K: ", species[iCO2].Cp(300), "J/mol/K")
print("CO2 Cp for T = 1400 K: ", species[iCO2].Cp(1400), "J/mol/K")
print("CO2 cp for T = 300 K: ", species[iCO2].cp(300), "J/kg/K")
print("CO2 cp for T = 1400 K: ", species[iCO2].cp(1400), "J/kg/K")
#
iH2O = species_names.index("H2O") 
print("H2 Cp for T = 300 K: ", species[iH2O].Cp(300), "J/mol/K")
print("H2 Cp for T = 1400 K: ", species[iH2O].Cp(300), "J/mol/K")
print("H2 cp for T = 300 K: ", species[iH2O].cp(300), "J/kg/K")
print("H2 cp for T = 1400 K: ", species[iH2O].cp(300), "J/kg/K")   

iC10H22 = species_names.index("C10H22") 
print("C10H22 Cp for T = 300 K: ", species[iC10H22].Cp(300), "J/mol/K")
print("C10H22 Cp for T = 1400 K: ", species[iC10H22].Cp(300), "J/mol/K")
print("C10H22 cp for T = 300 K: ", species[iC10H22].cp(300), "J/kg/K")
print("C10H22 cp for T = 1400 K: ", species[iC10H22].cp(300), "J/kg/K")
#---------------------------------------------------------------------------------------#		
# Calcul des propriétés du mélange
#---------------------------------------------------------------------------------------#
p = 101325
T = 300
y = [0.0, 0.233, 0.0, 0.0, 0.767]
#
mix = ThermoMix()                      # Création d'une instance de la classe ThermoMix()
#
thermomixprop(species, y, p, T, mix)   # Appel à la fonction thermomixprop
#
print("Mix properties (T=1400 K, p=101325 Pa)") 
print("Molecular weight: ", mix.mm, "kg/mol") 
print("cp: ", mix.cp, "J/kg/K") 
print("h: ", mix.h, "J/kg")
print("rho: ", mix.rho, "kg/mol") 
print('---------------------------------------------------------------------------------------')
#---------------------------------------------------------------------------------------#	
# Partie 2 - "Description de l’algorithme mis en œuvre" du rapport
#---------------------------------------------------------------------------------------#		
# Données du problème
#---------------------------------------------------------------------------------------#
M             = 0.85 
Tta           = 259      # K
Pta           = 46000    # Pa
rr_compFan    = 1.6
taux_compComp = 25
eff_isen      = 0.9
gamma         = 1.4
PCI_kero      = 43e6     # J.kg-1
cp_air        = 1005     # J K-1 kg-1
Tt4           = 1450     # K
#---------------------------------------------------------------------------------------#		
# Calcul de la pression atmosphérique
#---------------------------------------------------------------------------------------#
def Pression_atmospherique(Pta, gamma):
  Pa = Pta / ( (1 + ((gamma-1)/2) * M**2)**(gamma/(gamma-1)) )
  return Pa
Pa = Pression_atmospherique(Pta, gamma)
#print("Pa=", Pa,"Pa")
#---------------------------------------------------------------------------------------#		
# Calcul de la température atmosphérique
#---------------------------------------------------------------------------------------#
def Temperature_atmospherique(Tta, gamma, M):
  Ta = Tta / (1 + ((gamma-1)/2) * M**2)
  return Ta
Ta = Temperature_atmospherique(Tta, gamma, M)
#print("Ta=", Ta,"K")
#---------------------------------------------------------------------------------------#		
# Calcul de la vitesse de l'aéronef V
#---------------------------------------------------------------------------------------#
def vitesse_aeronef(Ta, M, gamma):
  V = M * math.sqrt(gamma * cp_air * (gamma-1) / gamma*Ta)
  return V
V = vitesse_aeronef(Ta, M, gamma)
#print("V=",V,"m/s") 
#---------------------------------------------------------------------------------------#	
# Calcul des grandeurs caractéristiques de l'écoulement jusqu'à la sortie de la turbine HP, qui sont indépendantes du BPR
#---------------------------------------------------------------------------------------#
i   = 0
P2  = Pta  #Pa
Tt2 = Tta  #K
T4  = Tt4  #Pour un mach faible
#---------------------------------------------------------------------------------------#
# Partie 3 - "Fonctions employées en python" du rapport
#---------------------------------------------------------------------------------------#
def grandeurs_cara(gamma, taux_compComp, eff_isen, Tt2, P2, i):
  Pt23    = 1.6 * P2
  Tt23_is = Tta * (P2/Pt23)**((1-gamma)/gamma)
  Tt23    = ((Tt23_is-Tt2)/eff_isen) + Tt2
  #Travail fourni par le fan dans le flux primaire
  Wf_FP   = cp_air * (Tt23-Tt2)
  Pt3     = taux_compComp * Pt23
  Tt3_is  = Tt23 * (Pt23/Pt3)**((1-gamma)/gamma)
  Tt3     = ((Tt3_is-Tt23)/eff_isen) + Tt23
  #Travail fourni par le compresseur HP dans le flux primaire
  Wc_HP   = cp_air * (Tt3-Tt23)
  Pt4     = Pt3 
  #Travail fourni par la turbine HP dans le flux primaire
  Wt_HP   = - Wc_HP
  Tt45    = Tt4 - Wc_HP/cp_air
  Tt45_is = Tt4 + (Tt45 - Tt4)/eff_isen
  Pt45    = Pt4 * (Tt45_is/Tt4)**(gamma/(gamma-1))
  #Rendement de la turbine HP
  rt_HP   = Pt45/Pt4 
  #Initialisation du travail rajouté à la turbine BP
  travail = 0
  #Boucle pour déterminer le BPR souhaité
  for i in range(1000000):
    #Incrémentation du travail de la turbine BP
    travail += 1
    #Travail nécessaire fourni par la turbine BP  
    Wt_BP   = -Wf_FP - travail
    Tt5     = Tt45 + Wt_BP/cp_air
    Tt5_is  = Tt45 + (Tt5 - Tt45)/eff_isen
    Pt5     = Pt45 * (Tt5_is/Tt45)**(gamma/(gamma-1))
    T9      = Tt5 * ( Pt5/Pa ) ** ( ( 1 - gamma )/gamma)
    T19     = Ta
    Tt13    = Tt5 - T9 + T19
    #Travail fourni par la compression du flux secondaire
    Wc_2    = cp_air * (Tt13-Tt2)
    #Formule du By-pass Ratio 
    BPR     = (Tt45 - Tt5 - Tt23 + Tt2)/(Tt13 - Tt2)
    #Condition pour déterminer le BPR voulu
    #print(BPR)
    if BPR >= 0: # Faire le choix d'un BPR
      break
  #La vitesse en sortie des deux flux     
  Vj      = math.sqrt(max(0, 2*cp_air*(Tt5-T9)))
  #Débit de combustible injecté
  m_f     = (cp_air*(Tt4-Tt3))/PCI_kero
  #Condition imposée pour un taux de compression égal à 0 pour le BPR=0
  if BPR > 0 and BPR < 0.0001:
    Pt13 = 0
  else:
    Tt13_is = eff_isen * (Tt13-Tt2) + Tt2
    Pt13    = P2 * (Tt2/Tt13_is)**(gamma/(1-gamma))
  #Poussée brute par unité de masse d’air
  Fg = Vj * (BPR+1)
  #Poussée nette par unité de masse d’air
  Fn = (1+BPR) * (Vj - V)
  #Specific Fuel Consumption
  SFC = (m_f/Fn)
  #Rendement propulsif
  rt_prop = 2 * V/(Vj+V)
  #Rendement global
  rt_global = 1/SFC * V/PCI_kero
  #Le taux de compression du fan pour flux secondaire
  Tauxcomp2 = Pt13/P2
  return Wt_BP, BPR, Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2
#---------------------------------------------------------------------------------------#
Wt_BP, BPR,Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2 = grandeurs_cara(gamma, taux_compComp, eff_isen, Tt2, P2, i)
#---------------------------------------------------------------------------------------#
# Partie 4 - "Résultats obtenus" du rapport
#---------------------------------------------------------------------------------------#
print("Travail spécifique consommé par la turbine BP =", Wt_BP,"kJ.kg-1")
print("BPR =",BPR)
print("Vj =",Vj,"m/s") 
print("Poussée brute=",Fg/1000,"kN.s−1.kg−1") # /1000 pour obtenir du kN.s−1.kg−1
print("Poussée nette=",Fn/1000,"kN.s−1.kg−1") # /1000 pour obtenir du kN.s−1.kg−1
#print("SFC=",SFC,"s m-1")
print("SFC=",SFC*(1000/1)*1000,"g.kN−1.s−1") # *1000 vient du Fn non diviser par 1000
print("SFC=",SFC*((3600*1000*0.001)/(100*1))*1000,"kg.h−1.DaN−1") # *1000 vient du Fn non diviser par 1000
print("Rendement propulsif=",rt_prop,"Rendement global=",rt_global)
print("Taux de compression du flux secondaire=",Tauxcomp2)
#---------------------------------------------------------------------------------------#	
# Liste des BPR allant de 0 à 14
#---------------------------------------------------------------------------------------#
bpr_values   = []
#---------------------------------------------------------------------------------------#
# Liste des valeurs de SFC pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
sfc_values_0 = [2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05, 2.2306883540441598e-05]
sfc_values  = [2.2306883540441598e-05, 1.8121494685125777e-05, 1.6308049725206232e-05, 1.5250650383894466e-05, 
1.454576372852258e-05, 1.4037752890283356e-05, 1.3652225702156564e-05, 1.3348645448300347e-05, 1.3102872371368161e-05, 
1.289947592934648e-05, 1.272824423795458e-05, 1.2581957896392301e-05, 1.2455365935847979e-05, 1.2344844534646727e-05,
1.2247374172733221e-05]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de Fg pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
fg_values_0 = [933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409, 933.442127193409]
fg_values   = [933.442127193409, 1346.1948491753878, 1695.250393129337, 2015.835729949671, 2320.2008946794804,
2614.150510635545, 2900.9132697229193, 3182.4638900828704, 3460.0644425228847, 3734.6595574133185, 4006.759215621649,
4276.882965442282, 4545.584883538227, 4812.727859034212, 5078.828056895296]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de Fn pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
fn_values_0 = [677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047, 677.0675164378047]
fn_values   = [677.0675164378047, 833.4448399883361, 926.123386467853, 990.3358780123913, 1038.327482838069, 
1075.9034124790946, 1106.2860054978698, 1131.4456059747397, 1152.6683470715284, 1170.8433986712596, 1186.5946281229767, 
1200.3907788091378, 1212.5911286737226, 1223.4472613896185, 1233.1840299138646]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de SFC/SFC(BPR=0), Fg/Fg(BPR=0), Fn/Fn(BPR=0)
#---------------------------------------------------------------------------------------#
sfc_sfc_0 = [a / b for a, b in zip(sfc_values, sfc_values_0)]
fg_fg_0   = [a / b for a, b in zip(fg_values, fg_values_0)]
fn_fn_0   = [a / b for a, b in zip(fn_values, fn_values_0)]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de vj pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
Vj_0 = [933.4398097645088, 673.0947195131911,565.07972921447,503.95698894568983, 464.03889725801326, 435.69083051346036, 414.4140004857358, 397.8028483547519, 384.4465517129795, 373.45482435817297, 364.24426736418263, 356.4063929913336, 349.64651831020257, 343.76196975845437, 338.5855072399162]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de taux de compression pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
taux_compComp_0 = [ 6.01394082084749, 3.5693002700956264, 2.680770504955299, 2.2369289654334255, 1.9744460371801067, 1.8021778428816453, 1.6808747668918782, 1.5910341137153363, 1.5218999904410477, 1.4671237138010875, 1.4226722300967265, 1.3858680367693086, 1.3549468928735815, 1.328581025078975]
x_taux_compComp_0 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de rendement_propulsif pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
rendement_0 = [0.43094806549374776, 0.5516570401542062, 0.624195796250013, 0.6743746780491904, 0.7117417925115068, 0.7408958597314832, 0.7643964528685019, 0.783806351428026, 0.8001428289830327, 0.814106864670751, 0.8261889743693008, 0.8367564888253609, 0.8460901154475806, 0.8543863330181377, 0.8618199465106985]
#---------------------------------------------------------------------------------------#
# Fonction pour tracer les courbes en fonction du BPR
#---------------------------------------------------------------------------------------#
def tracer_sfc_bpr():
    bpr_values = np.arange(0, 15, 1)
    for bpr_value in np.arange(0, 15, 1):
        Wt_BP, BPR, Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2 = grandeurs_cara(gamma, taux_compComp, eff_isen, Tt2, P2, i)
        #sfc_values.append(SFC)
        #fg_values.append(Fg)
        #fn_values.append(Fn)
    #print(sfc_values)
    #print(fg_values)
    #print(fn_values)
    plt.figure(1)
    plt.plot(bpr_values,sfc_sfc_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('SFC/SFC(BPR=0)')
    plt.title('Évolution de la consommation spécifique réduite par rapport au taux de dilution.')
    plt.grid(True)
    plt.figure(2)
    plt.plot(bpr_values,fg_fg_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Fg/Fg(BPR=0)')
    plt.title('Évolution des poussées brutes par rapport au taux de dilution.')
    plt.grid(True)   
    plt.figure(3)
    plt.plot(bpr_values,fn_fn_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Fn/Fn(BPR=0)')
    plt.title('Évolution des poussées nettes par rapport au taux de dilution.')
    plt.grid(True)
    plt.figure(4)
    plt.plot(bpr_values, Vj_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Vj')
    plt.title('Évolution de la vitesse du jet par rapport au taux de dilution.')
    plt.grid(True)
    plt.figure(5)
    plt.plot(bpr_values, marker=' ', linestyle='-')
    plt.xlabel('')
    plt.ylabel('Taux de dilution(BPR)')
    plt.title('Évolution du BPR.')
    plt.grid(True)
    plt.figure(6)
    plt.plot(x_taux_compComp_0, taux_compComp_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Taux de compression du fan')
    plt.title('Évolution du taux de compression du fan par rapport au taux de dilution.')
    plt.grid(True)
    plt.figure(7)
    plt.plot(bpr_values,rendement_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Rendement propulsif')
    plt.title('Évolution du rendement propulsif par rapport au taux de dilution')
    plt.grid(True)  
    plt.show()
#---------------------------------------------------------------------------------------#
courbes=tracer_sfc_bpr()
#---------------------------------------------------------------------------------------#
# Partie 4.2 "Application de l’exercice vu en cours" du rapport
#---------------------------------------------------------------------------------------#
print('---------------------------------------------------------------------------------------')
#Calcul de la chaleur de formation
Q      = species[iC10H22].Href - 10 * species[iCO2].Href - 11 * species[iH2O].Href
print('Chaleur de formation :',Q,'J.mol-1')
# Ti_mel visée = 4849.4 K 
Ti_mel = 300 #Prendre la température 300K puis les bornes [4800K, 4887.1K]
# Fonction de l'intégrale de Cp dT pour le C10H22 
def cp_function1(T):    
    Hs_C10H22 = species[iC10H22].Cp(Ti_mel) * T
    return Hs_C10H22  
# Fonction de l'intégrale de Cp dT pour le O2     
def cp_function2(T): 
    Hs_O2     = species[iO2].Cp(Ti_mel) * T    
    return Hs_O2 
# Fonction de l'intégrale de Cp dT pour le N2 
def cp_function3(T):
    Hs_N2     = species[iN2].Cp(Ti_mel) * T
    return Hs_N2       
# Bornes de l'intégrale
T_min = 300   # Température minimale
T_max = 1400  # Température maximale
# Calcul de l'intégrale Cp dT
result_C10H22, err = quad(cp_function1, T_min, T_max)
result_O2, err     = quad(cp_function2, T_min, T_max)
result_N2, err     = quad(cp_function3, T_min, T_max)
print('Enthalpie sensible du C10H22 :',result_C10H22,'J.mol-1','Enthalpie sensible du O2 :', result_O2,'J.mol-1','Enthalpie sensible du N2 :', result_N2,'J.mol-1')
#Calcul de l'enthalpie sensible d'un mélange air-kérosène 
Hs_Tu = result_C10H22 + 31/2 * result_O2 + 31/2*3.76 * result_N2
print(Hs_Tu)
#Calcul de l'enthalpie du mélange frais de l'air-kérosène 
Hs_Tb = Hs_Tu + Q
print('Enthalpie du mélange frais :',Hs_Tb,'J.mol-1')
Hs_Tb_300  = 2238231007.7887344       # J.mol-1
Tbrulee    = (Ti_mel*Hs_Tb_300)/Hs_Tb # Déterminer la température brûlée
print(Tbrulee)                        # valeurs calculées de ci-dessous
Tb_4800K   = 4435.9000567329695       # K
Tb_4887_1K = 5262.550823647796        # K
Tb_moy     = (Tb_4800K+Tb_4887_1K)/2
print('La température brûlée moyenne est de',Tb_moy,'K')
#---------------------------------------------------------------------------------------#
Tb4 = Tb_moy
def grandeurs_cara2(gamma, taux_compComp, eff_isen, Tt2, P2, i):
  Pt23    = 1.6 * P2
  Tt23_is = Tta * (P2/Pt23)**((1-gamma)/gamma)
  Tt23    = ((Tt23_is-Tt2)/eff_isen) + Tt2
  #Travail fourni par le fan dans le flux primaire
  Wf_FP   = mix.cp * (Tt23-Tt2)
  Pt3     = taux_compComp * Pt23
  Tt3_is  = Tt23 * (Pt23/Pt3)**((1-gamma)/gamma)
  Tt3     = ((Tt3_is-Tt23)/eff_isen) + Tt23
  #Travail fourni par le compresseur HP dans le flux primaire
  Wc_HP   = mix.cp * (Tt3-Tt23)
  Pt4     = Pt3 
  #Travail fourni par la turbine HP dans le flux primaire
  Wt_HP   = - Wc_HP
  Tt45    = Tb4 - Wc_HP/mix.cp
  Tt45_is = Tb4 + (Tt45 - Tb4)/eff_isen
  Pt45    = Pt4 * (Tt45_is/Tb4)**(gamma/(gamma-1))
  #Rendement de la turbine HP
  rt_HP   = Pt45/Pt4 
  #Initialisation du travail rajouté à la turbine BP
  travail = 0
  #Boucle pour déterminer le BPR souhaité
  for i in range(100000000000):
    #Incrémentation du travail de la turbine BP
    travail += 1
    #Travail nécessaire fourni par la turbine BP  
    Wt_BP   = -Wf_FP - travail
    Tt5     = Tt45 + Wt_BP/mix.cp
    Tt5_is  = Tt45 + (Tt5 - Tt45)/eff_isen
    Pt5     = Pt45 * (Tt5_is/Tt45)**(gamma/(gamma-1))
    T9      = Tt5 * ( Pt5/Pa ) ** ( ( 1 - gamma )/gamma)
    T19     = Ta
    Tt13    = Tt5 - T9 + T19
    #Travail fourni par la compression du flux secondaire
    Wc_2    = mix.cp * (Tt13-Tt2)
    #Formule du By-pass Ratio 
    BPR     = (Tt45 - Tt5 - Tt23 + Tt2)/(Tt13 - Tt2)
    #Condition pour déterminer le BPR voulu
    #print(BPR)
    if BPR >= 0: # Faire le choix d'un BPR
      break
  #La vitesse en sortie des deux flux     
  Vj      = math.sqrt(max(0, 2*mix.cp*(Tt5-T9)))
  #Débit de combustible injecté
  m_f     = (mix.cp*(Tb4-Tt3))/PCI_kero
  #Condition imposée pour un taux de compression égal à 0 pour le BPR=0
  if BPR > 0 and BPR < 0.0001:
    Pt13 = 0
  else:
    Tt13_is = eff_isen * (Tt13-Tt2) + Tt2
    Pt13    = P2 * (Tt2/Tt13_is)**(gamma/(1-gamma))
  #Poussée brute par unité de masse d’air
  Fg = Vj * (BPR+1)
  #Poussée nette par unité de masse d’air
  Fn = (1+BPR) * (Vj - V)
  #Specific Fuel Consumption
  SFC = (m_f/Fn)
  #Rendement propulsif
  rt_prop = 2 * V/(Vj+V)
  #Rendement global
  rt_global = 1/SFC * V/PCI_kero
  #Le taux de compression du fan pour flux secondaire
  Tauxcomp2 = Pt13/P2
  return Wt_BP, BPR, Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2
#---------------------------------------------------------------------------------------#
Wt_BP, BPR,Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2 = grandeurs_cara2(gamma, taux_compComp, eff_isen, Tt2, P2, i)
#---------------------------------------------------------------------------------------#
# Partie 4 - "Résultats obtenus" du rapport
#---------------------------------------------------------------------------------------#
print('---------------------------------------------------------------------------------------')
print("Travail spécifique consommé par la turbine BP =", Wt_BP,"kJ.kg-1")
print("BPR =",BPR)
print("Vj =",Vj,"m/s") 
print("Poussée brute=",Fg/1000,"kN.s−1.kg−1") # /1000 pour obtenir du kN.s−1.kg−1
print("Poussée nette=",Fn/1000,"kN.s−1.kg−1") # /1000 pour obtenir du kN.s−1.kg−1
#print("SFC=",SFC,"s m-1")
print("SFC=",SFC*(1000/1)*1000,"g.kN−1.s−1") # *1000 vient du Fn non diviser par 1000
print("SFC=",SFC*((3600*1000*0.001)/(100*1))*1000,"kg.h−1.DaN−1") # *1000 vient du Fn non diviser par 1000
print("Rendement propulsif=",rt_prop,"Rendement global=",rt_global)
print("Taux de compression du flux secondaire=",Tauxcomp2)
#---------------------------------------------------------------------------------------#
# Liste des valeurs de SFC pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
sfc2_values_0 = [4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05, 4.469399292027858e-05]
sfc2_values  = [4.469399292027858e-05, 3.360792841167953e-05, 2.869416006408007e-05, 2.575516233709142e-05, 2.3744730661086466e-05, 2.2259238556487012e-05, 2.110495066889216e-05, 2.0175602924917216e-05, 1.9407213283019673e-05, 1.875873438372404e-05, 1.8202414682066107e-05, 1.771866035194745e-05, 1.729332811024071e-05, 1.6915724476612172e-05, 1.6577833463455425e-05]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de Fg pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
fg2_values_0 = [2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154, 2379.458398032154]
fg2_values = [2379.458398032154, 3336.1625628860984, 4076.036692601508, 4709.7729009993045, 5278.088493655373, 5801.1536818850245, 6290.6818418129515, 6754.153840884414, 7196.745047102426, 7622.147714699239, 8033.117260015606, 8431.825962879671, 8819.903071682525, 9198.780794509039, 9569.478801071918]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de Fn pour 0 et allant de 0 à 14
#---------------------------------------------------------------------------------------#
fn2_values_0 = [2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867, 2123.0843321401867]
fn2_values = [2123.0843321401867, 2823.414610608719, 3306.9138771764224, 3684.2755975633227, 3996.217833093169, 4262.909347461678, 4496.059602247286, 4703.161360924594, 4889.373591459987, 5058.396487140312, 5212.996064929627, 5355.321126148156, 5487.03612774436, 5609.521261771696, 5723.85506942169]
#---------------------------------------------------------------------------------------#
# Liste des valeurs de SFC/SFC(BPR=0), Fg/Fg(BPR=0), Fn/Fn(BPR=0)
#---------------------------------------------------------------------------------------#
sfc_sfc2_0 = [a / b for a, b in zip(sfc_values, sfc_values_0)]
fg_fg2_0   = [a / b for a, b in zip(fg_values, fg_values_0)]
fn_fn2_0   = [a / b for a, b in zip(fn_values, fn_values_0)]
#---------------------------------------------------------------------------------------#
# Fonction pour tracer les courbes en fonction du BPR
#---------------------------------------------------------------------------------------#
def tracer_sfc_bpr():
    bpr_values = np.arange(0, 15, 1)
    for bpr_value in np.arange(0, 15, 1):
        Wt_BP, BPR, Vj, Fg, Fn, SFC, rt_prop, rt_global, Tauxcomp2 = grandeurs_cara2(gamma, taux_compComp, eff_isen, Tt2, P2, i)
    #     sfc2_values.append(SFC)
    #     fg2_values.append(Fg)
    #     fn2_values.append(Fn)
    # print(sfc2_values)
    # print(fg2_values)
    # print(fn2_values)
    plt.figure(8)
    plt.plot(bpr_values,sfc_sfc2_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('SFC/SFC(BPR=0)')
    plt.title('Évolution de la consommation spécifique réduite (T=4849.23K) par rapport au taux de dilution.')
    plt.grid(True)
    plt.figure(9)
    plt.plot(bpr_values,fg_fg2_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Fg/Fg(BPR=0)')
    plt.title('Évolution des poussées brutes (T=4849.23K) par rapport au taux de dilution.')
    plt.grid(True)   
    plt.figure(10)
    plt.plot(bpr_values,fn_fn2_0, marker=' ', linestyle='-')
    plt.xlabel('BPR')
    plt.ylabel('Fn/Fn(BPR=0)')
    plt.title('Évolution des poussées nettes (T=4849.23K) par rapport au taux de dilution.')
    plt.grid(True)  
    plt.show()
#---------------------------------------------------------------------------------------#
courbes=tracer_sfc_bpr()
#---------------------------------------------------------------------------------------#