#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""

# usage : python2 main_Barstar.py start_prot_only.pdb md_prot_only_skip100.pdb all

import sys, os
from math import sqrt
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import csv

import Parsing_dico as parse
import Conformation_analysis as conf_anal
import Ecriture as write
import Graphes as graph


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
	Centre de Masse des residus calcule a partir de la position du carbone alpha "CM_CA" ou de la distance moyenne separant les atomes "CM_moyall"
	'''			
	centreMasse = conf_anal.choixMeth()			


	'''
	Calcul des valeurs des variables d'interet (CdM, Distance au CdM, RMSD, Correlation...) et ajout dans les dictionnaires a des cles pertinentes
	'''
	conf_anal.conformation_analysis(liste_dictionnaire)	



	#string globale servant a donner le bon nom aux dossiers et fichiers en fonction du nombre de conformations analysees et de la methode pour le calcul du centre de masse
	global type_analyse
	type_analyse=centreMasse+"_"+str(len(liste_dictionnaire[1]["liste_conformations"])-1)	

	'''
	Ecriture des resultats dans les tableaux excel de sortie, en fonction de la methode de calcul des CdM
	'''
	write.ecriture(liste_dictionnaire,type_analyse)			
	graph.plotRes(liste_dictionnaire) 						# representations graphiques et enregistrement des images


	cwd = os.getcwd()
	for file in os.listdir(cwd):
		if file.endswith('.pyc'):
			print "suppression de "+file
			rm = 'rm '+file
			os.popen(rm)