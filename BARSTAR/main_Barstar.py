#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""

# usage : python2 main_Barstar.py start_prot_only.pdb md_prot_only_skip100.pdb all

#Import des modules
import sys, os
from math import sqrt
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import csv

#Import des fonctions
import Parsing_dico as parse
import Conformation_analysis as conf_anal
import Ecriture as write
import Graphes as graph

'''
Fonction MAIN du projet BARSTAR.
On regroupe ici les 4 importantes etapes de notre programme :
-Verification du bon usage du programme et Parsing des fichiers pdb donnes en arguments
-Calcul des parametres permettant l'analyse de ces fichiers pdb, en prenant en compte la methode de calcul du centre de masse choisie en argument
-Creation du dossier (avec nom pertinent) dans lequel seront enregistres les resultats et ecriture des resultats (tableaux excel des valeurs des parametres calcules et images des graphes)
-Representations graphiques pertinentes des parametres calcules afin d'analyser les changemements conformationnels globaux et locaux.
'''

#-----------------------------------------------------------------------
# MAIN
#-----------------------------------------------------------------------

if __name__ == '__main__':
	'''
	Creation d'une liste : [dictionnaire_des_conformations,dictionnaire_de_la_reference]
	Chaque dictionnaire possede la position de chaque atome de chaque residu de chaque conformation.
	'''
	liste_dictionnaire = parse.dictionnaire()


	'''
	Choix de la methode de calcul du centre de masse.
	Centre de Masse des residus calcule a partir de la position du carbone alpha "CM_CA" ou de la distance moyenne separant les atomes "CM_moyall"
	'''			
	centreMasse = conf_anal.choixMeth()			


	'''
	Calcul des valeurs des variables d'interet (Centre de Masse, Distance au Centre de Masse, RMSD, Correlation...) et ajout dans les dictionnaires, associees a des cles pertinentes
	'''
	conf_anal.conformation_analysis(liste_dictionnaire)	



	#string globale servant a donner le bon nom aux dossiers et fichiers en fonction du nombre de conformations analysees et de la methode pour le calcul du centre de masse
	#Pour l'analyse d'un fichier pdb de 200 conformations avec la methode d'approximation du centre de masse des residus par la position du carbone alpha : type_analyse= "CA_200"
	global type_analyse
	type_analyse=centreMasse+"_"+str(len(liste_dictionnaire[1]["liste_conformations"])-1)	


	'''
	Ecriture des resultats dans les tableaux excel de sortie, dans les dossiers correspondants, en fonction du nombre de conformations analysees et de la methode de calcul du centre de masse
	'''
	write.ecriture(liste_dictionnaire,type_analyse)			


	'''
	Representations graphiques des resultats et enregistrement des images dans les dossiers correspondants.
	'''
	graph.plotRes(liste_dictionnaire) 						


	#Suppression des fichiers .pyc une fois le programme execute.
	cwd = os.getcwd()
	for file in os.listdir(cwd):
		if file.endswith('.pyc'):
			print "suppression de "+file
			rm = 'rm '+file
			os.popen(rm)
	
	print "\nMerci d'avoir utilise notre programme.\nBonne analyse !\n"
