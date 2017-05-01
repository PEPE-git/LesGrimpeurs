#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""


import sys, os
import matplotlib.pyplot as plt

import Conformation_analysis as conf_anal

#-----------------------------------------------------------------------
# GRAPHIQUES
#-----------------------------------------------------------------------

def plotRes(l_dict) :
	'''
	Affiche et enregistre dans le fichier correspondant les plot de l'anayle globale, puis ceux de l'analyse locale
	'''
	centreMasse = conf_anal.choixMeth()
	global type_analyse
	type_analyse=centreMasse+"_"+str(len(l_dict[1]["liste_conformations"])-1)


	global decision
	decision = raw_input("\nVoulez-vous afficher les plot un par un maintenant ? O/N\n")
	
	while ((decision != 'O') & (decision != 'o') & (decision != 'N')  & (decision != 'n')) :
		print "Erreur : Repondre O pour oui ou N pour non" 
		decision = raw_input("\nVoulez-vous afficher les plot un par un maintenant ? O/N\n")

	plotGlobal(l_dict[1])
	plotLocal(l_dict[0]) 

	print "Les images des plot ont ete enregistrees dans Barstar_Results_"+type_analyse+"/\n"


def plotGlobal(d_conf) :
	plotGiration(d_conf)
	plotDistance(d_conf)
	plotGlobalRMSD(d_conf)
	plotFlexibiteEnfouissement(d_conf)
	
def plotGiration(d_conf) :

	plt.title('Rayon de Giration en fonction des conformations')
	plt.plot(d_conf["rayonGiration"])
	plt.axhline(y=d_conf["rayonGiration"][0],ls='--',color='black')

	plt.xlabel('Conformations')
	plt.ylabel('Rayon de Giration')


	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/Giration_conf.png")

	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()

def plotDistance(d_conf) :

	moy = d_conf["distance_moy"]
	sd = d_conf["distance_sd"]
	moy_s = [x+y for (x,y) in zip(moy,sd)]
	moy_i = [x-y for (x,y) in zip(moy,sd)]

	plt.plot(moy, "b", label = "Distance moyenne")
	plt.plot(moy_s, "r--", label = "+/- ecart-type")
	plt.plot(moy_i, "r--",)
	plt.title('Distance moyenne des residus(+/- ecart-type) \npar rapport au CdM, en fonction des conformations')
	plt.xlabel('Conformations')
	plt.ylabel('Distance moyenne')
	
	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/Distance_conf.png")
	
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()

def plotGlobalRMSD(d_conf) :

	moy = d_conf["RMSDmoy"]
	sd = d_conf["RMSDmoy_sd"]
	moy_s = [x+y for (x,y) in zip(moy,sd)]
	moy_i = [x-y for (x,y) in zip(moy,sd)]
	
	plt.plot(moy, "b", label = "RMSD moyen")
	plt.plot(moy_s, "r--", label = "+/- ecart-type")
	plt.plot(moy_i, "r--",)
	plt.title('RMSD moyen des residus(+/- ecart-type) \npar rapport a la ref, en fonction des conformations')
	plt.xlabel('Conformations')
	plt.ylabel('RMSD moyen')

	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/RMSD_conf.png")
		
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()
	


def plotFlexibiteEnfouissement(d_conf) :
	plt.subplot(211)
	plt.plot(d_conf["corEnfFlexi_conf"][0])
	plt.title('Correlation de la Flexibilite moyenne des residus\n et leur Enfouissement moyen, en fonction des conformations')
	plt.xlabel('Conformations')
	plt.ylabel('Correlation')
	
	plt.subplot(212)
	plt.plot(d_conf["corEnfFlexi_conf"][1])
	plt.title('p-value de la Correlation Flexibilite/Enfouissement moyens \n en fonction des conformations')
	plt.xlabel('Conformations')
	plt.ylabel('p-valeur')
	
	plt.tight_layout()

	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/Correlation_conf.png")
		
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()

#-----------------------------------------------------------------------
def plotLocal(d_ref) :
	#Appel des fonctions de plot pour les graphes de l'analyse locale
	plotDistanceLocal(d_ref)
	plotRMSDLocal(d_ref)
	plotFlexibiteEnfouissement_residus(d_ref)

def plotDistanceLocal(d_ref) :

	moy = d_ref["enfRes_mean"]
	sd = d_ref["enfRes_sd"]
	moy_s = [x+y for (x,y) in zip(moy,sd)]
	moy_i = [x-y for (x,y) in zip(moy,sd)]

	plt.plot(moy, "b", label = "Distance moyenne")
	plt.plot(moy_s, "r--", label = "+/- ecart-type")
	plt.plot(moy_i, "r--",)
	plt.xlabel('Residus')
	plt.ylabel('Distance')
	plt.title('Distance moyenne (+/- ecart-type) des residus par rapport au CdM \npour chaque conformation, en fonction du residu')

	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/Distance_residu.png")
		
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()



def plotRMSDLocal(d_ref) :
	moy = d_ref["RMSDres_mean"]
	sd = d_ref["RMSDres_sd"]
	moy_s = [x+y for (x,y) in zip(moy,sd)]
	moy_i = [x-y for (x,y) in zip(moy,sd)]
	
	plt.plot(moy, "b", label = "RMSD moyen")
	plt.plot(moy_s, "r--", label = "+/- ecart-type")
	plt.plot(moy_i, "r--",)
	plt.xlabel('Residus')
	plt.ylabel('RMSD')
	plt.title('RMSD moyen (+/- ecart-type), par rapport a la reference, de la position \ndes residus pour chaque conformation, en fonction du residu')

	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/RMSD_residu.png")
		
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()


def plotFlexibiteEnfouissement_residus(d_ref) :
	plt.subplot(211)
	plt.plot(d_ref["corEnfFlexi_ref"][0])
	plt.title("Correlation entre la Flexibilite et l'Enfouissement de chaque residu \nde chaque conformation, en fonction des residus")
	plt.xlabel('Residus')
	plt.ylabel('Correlation')
	
	plt.subplot(212)
	plt.plot(d_ref["corEnfFlexi_ref"][1])
	plt.axhline(y=0.05,ls='--',color='black')
	plt.title('p-value de la Correlation Flexibilite/Enfouissement \n en fonction des residus')
	plt.xlabel('Residus')
	plt.ylabel('p-valeur')
	
	plt.tight_layout()

	global type_analyse
	plt.savefig("Barstar_Results_"+type_analyse+"/Correlation_residu.png")
		
	global decision
	if ((decision == 'O') | (decision == 'o')) :
		plt.show()