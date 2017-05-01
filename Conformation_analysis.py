#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""

import sys, os
from math import sqrt
from scipy.stats.stats import pearsonr
import numpy as np


def choixMeth() :
	'''
	Deux methodes de calcul sont proposees :
	- soit en utilisant la position des carbones alpha (hyp. representatif des residus)
	- soit en calculant la position des residus a partir du centre de masse de l'ensemble des atomes les consituant
	'''
	global centreMasse
	
	centreMasse = sys.argv[3]
	
	if sys.argv[3] == "CA" :
		centreMasse = "CM_CA"
		return "CM_CA"
	else :
		centreMasse = "CM_moyAll"
		return "CM_moyAll"


def conformation_analysis(l_dict) :
	'''
	Analyse conformationnelle basee sur le Rayon de Giration, le RMSD
	et l'eventuelle correlation entre l'enfouissement et la flexibilite des residus
	'''
	centreMasseCalc(l_dict)
	RMSD(l_dict)						
	distance(l_dict)
	rayonGiration(l_dict)
	corEnfouissementFlexibilite_conf(l_dict[1])
	corEnfouissementFlexibilite_res(l_dict[0])



#-----------------------------------------------------------------------
# PARTIE 2 CALCUL DU CENTRE DE MASSE
#-----------------------------------------------------------------------

def centreMasseCalc(l_dict) :
	print "Calcul des Centres de Masse des Residus et des proteines"
	__centreMasseResidus(l_dict[0])
	__centreMasseResidus(l_dict[1])
	__centreMasseProteine(l_dict[0])
	__centreMasseProteine(l_dict[1])


	
def __centreMasseProteine(d_prot) :
	'''
		Fonction qui calcule le centre de masse des proteines.
		Cela correspond a la position moyenne des residus la constituant.
		On cree des nouvelles cles au dictionnaire donne en argument, dont les valeurs
		correspondent 
	'''
	global centreMasse
	
	l_CM = []
	l_CMx = [] # stockage des coordonnees des CM des residus sur l'axe x
	l_CMy = [] # stockage des coordonnees des CM des residus sur l'axe y
	l_CMz = [] # stockage des coordonnees des CM des residus sur l'axe z
		
	for conf in d_prot["liste_conformations"] :
		l_CM = [list(),list(),list()] 	# liste des coordonnees (x,y,z) des CM pour chaque conformation
		
		for resid in d_prot[conf]["liste_n_residus"] :
			l_CM[0].append(d_prot[conf][resid][centreMasse]["x"])
			l_CM[1].append(d_prot[conf][resid][centreMasse]["y"])
			l_CM[2].append(d_prot[conf][resid][centreMasse]["z"])
		l_CMx.append(moyenne(l_CM[0]))
		l_CMy.append(moyenne(l_CM[1]))
		l_CMz.append(moyenne(l_CM[2]))
	d_prot["liste_CM_x"] = l_CMx
	d_prot["liste_CM_y"] = l_CMy
	d_prot["liste_CM_z"] = l_CMz
	
	return d_prot


def __centreMasseResidus(d_prot) :
	'''
	Le centre d'un masse d'un residus peut etre considere comme etant :
		- la position moyenne de ses atomes
		- la position de ses carbones alpha
	'''
	global centreMasse
	
	if centreMasse == "CM_CA" :
		return __centreMasseResCa(d_prot)
	return __centreMasseResAll(d_prot)


def __centreMasseResCa(d_prot) :
	'''		
	Fonction qui calcule le centre de masse des residus selon la position de leur carbone alpha.
	
	@param dictionnaire de la proteine avec une ou plusieurs conformations
	'''
	
	# Pour toutes les conformations de la proteine
	for conf in d_prot["liste_conformations"] :
		
		
		d_prot[conf]["CM_res"] = dict()
		d_prot[conf]["CM_res"]["x"] = list()
		d_prot[conf]["CM_res"]["y"] = list()
		d_prot[conf]["CM_res"]["z"] = list()

		# Pour tous les residus de chaque conformation de la proteine
		for resid in d_prot[conf]["liste_n_residus"] :
			d_prot[conf][resid]["CM_CA"] = dict()
			d_prot[conf][resid]["CM_CA"]["x"] = d_prot[conf][resid]["CA"]["x"]
			d_prot[conf][resid]["CM_CA"]["y"] = d_prot[conf][resid]["CA"]["y"]
			d_prot[conf][resid]["CM_CA"]["z"] = d_prot[conf][resid]["CA"]["z"]
			
			d_prot[conf]["CM_res"]["x"].append(d_prot[conf][resid]["CA"]["x"])
			d_prot[conf]["CM_res"]["y"].append(d_prot[conf][resid]["CA"]["y"])
			d_prot[conf]["CM_res"]["z"].append(d_prot[conf][resid]["CA"]["z"])

def __centreMasseResAll(d_prot) :
	'''		
	Fonction qui calcule le centre de masse des residus
	Le centre de masse d'un residu correspond a la position moyenne
	des atomes le constituant
	
	@param dictionnaire de la proteine avec une ou plusieurs conformations
	'''
	
	# Pour toutes les conformations de la proteine
	for conf in d_prot["liste_conformations"] :
		# Pour tous les residus de chaque conformation de la proteine
		for resid in d_prot[conf]["liste_n_residus"] :
			lcoord = [list(),list(),list()] # stockage coordonnees [x,y,z] pour tous les atomes de chaque residus
			
			
			d_prot[conf]["CM_res"] = dict()
			d_prot[conf]["CM_res"]["x"] = list()
			d_prot[conf]["CM_res"]["y"] = list()
			d_prot[conf]["CM_res"]["z"] = list()

			for atom in d_prot[conf][resid]["liste_atomes"] :
				lcoord[0].append(d_prot[conf][resid][atom]["x"])
				lcoord[1].append(d_prot[conf][resid][atom]["y"])
				lcoord[2].append(d_prot[conf][resid][atom]["z"])

			xmoy = moyenne(lcoord[0])
			ymoy = moyenne(lcoord[1])
			zmoy = moyenne(lcoord[2])
			
			d_prot[conf][resid]["CM_moyAll"] = dict()
			d_prot[conf][resid]["CM_moyAll"]["x"] = xmoy
			d_prot[conf][resid]["CM_moyAll"]["y"] = ymoy
			d_prot[conf][resid]["CM_moyAll"]["z"] = zmoy

			d_prot[conf]["CM_res"]["x"].append(xmoy)
			d_prot[conf]["CM_res"]["y"].append(ymoy)
			d_prot[conf]["CM_res"]["z"].append(zmoy)

#-----------------------------------------------------------------------
# PARTIE 3 : CALCUL DU RMSD
#-----------------------------------------------------------------------
'''
Le RMSD est considere comme la distance moyenne entre les residus
	
'''

def RMSD(l_dict) :
	'''	
	Fonction qui calcule le RMSD
	Fonction permettant de realiser l'ensemble des calculs relatifs au RMSD
	C'est la seule qui peut etre utilisee par l'utilisateur
	'''

	print "Calcul des RMSD"
	__RMSDresidus(l_dict[0],l_dict[1])
	__RMSDconf(l_dict[1])
	__RMSDres(l_dict[0], l_dict[1])


def __RMSDresidus(d_ref, d_conf) :
	'''
	Calcul du RMSD entre la position de chaque residu de chaque conformation et
	sa position dans la reference
	'''
	global centreMasse
	
	REF = d_ref["liste_conformations"][0] # il n'y a qu'une seule conformation de reference
	
	for conf in d_conf["liste_conformations"] :
		d_conf[conf]["RMSD"] = list()
		for resid in d_conf[conf]["liste_n_residus"] :
			xref = d_ref[REF][resid][centreMasse]["x"]
			yref = d_ref[REF][resid][centreMasse]["y"]
			zref = d_ref[REF][resid][centreMasse]["z"]

			xconf = d_conf[conf][resid][centreMasse]["x"]
			yconf = d_conf[conf][resid][centreMasse]["y"]
			zconf = d_conf[conf][resid][centreMasse]["z"]
			
			d_conf[conf]["RMSD"].append(sqrt((xref-xconf)**2 + (yref-yconf)**2 + (zref-zconf)**2))
	

def __RMSDconf(d_conf) :
	'''
	Calcul du RMSD de maniere globale a la conformation.
	Moyenne du RMSD des residus.
	'''
	
	d_conf["RMSDmoy"] = list()
	d_conf["RMSDmoy_sd"] = list()
	
	for conf in d_conf["liste_conformations"] :
		d_conf["RMSDmoy"].append(moyenne(d_conf[conf]["RMSD"]))
		d_conf["RMSDmoy_sd"].append(ecart_type(d_conf[conf]["RMSD"]))


def __RMSDres(d_ref,d_conf) :
	'''
	Calcul du RMSD de maniere locale.
	Moyenne du RMSD en chaque position a partir de l'ensemble des conformations.
	'''
	
	l_res = d_ref[d_ref["liste_conformations"][0]]["liste_n_residus"]

	d_ref["RMSDres_mean"] = list()
	d_ref["RMSDres_sd"] = list()
	d_ref["list_RMSDres"] = list()

	for i in range(len(l_res)) :
		l_rmsd = list() # liste contenant les valeurs des RMSD d'un residu pour toutes ses conformations
		for conf in d_conf["liste_conformations"] :
			l_rmsd.append(d_conf[conf]["RMSD"][i])
		d_ref["list_RMSDres"].append(l_rmsd)
		d_ref["RMSDres_mean"].append(moyenne(l_rmsd))
		d_ref["RMSDres_sd"].append(ecart_type(l_rmsd))


#-----------------------------------------------------------------------
# PARTIE 4 : DISTANCE
#-----------------------------------------------------------------------
'''
Calcul de la distance entre chaque residu et le centre de Masse de la conformation.
La distance calculee pour un residu represente son enfouissement au sein de la conformation.
Plus la distance est grande, plus ce residu est en peripherie (faible enfouissement relatif)
'''

def distance(l_dict) :

	print "Calcul des distances Residus-Centre de Masse"
	__distanceConf(l_dict[0])
	__distanceConf(l_dict[1])
	__distanceRes(l_dict[0], l_dict[1])
	
	
def __distanceConf(d_conf) :
	'''
	Calcul de la distance entre le CdM de chaque residu et le CdM de la conformation.
	Creation des cles :
		-enfouissement : a pour valeur la liste de l'enfouissement de chaque residu de chaque conformation
		-distance_moy : a pour valeur la liste de l'enfouissement moyen des residus de chaque conformation
		-distance_sd : a pour valeur la liste de l'ecart-type de l'enfouissement des residus de chaque conformation
	'''
	global centreMasse
	# pour chaque conformation, son centre de masse :
	xconf = d_conf["liste_CM_x"]
	yconf = d_conf["liste_CM_y"]
	zconf = d_conf["liste_CM_z"]
	
	d_conf["distance_moy"] = list()
	d_conf["distance_sd"] = list()

	for conf in d_conf["liste_conformations"] :
		l_dist = []
		
		for resid in d_conf[conf]["liste_n_residus"] :
			i=0
			xres = d_conf[conf][resid][centreMasse]["x"]
			yres = d_conf[conf][resid][centreMasse]["y"]
			zres = d_conf[conf][resid][centreMasse]["z"]
			
			distance = sqrt((xres-xconf[i])**2 + (yres-yconf[i])**2 + (zres-zconf[i])**2)
			l_dist.append(distance)
			
			i += 1
		
		d_conf[conf]["enfouissement"] = l_dist
		d_conf["distance_moy"].append(moyenne(l_dist))
		d_conf["distance_sd"].append(ecart_type(l_dist))


def __distanceRes(d_ref, d_conf) :
	'''
	Calcul de la distance moyenne de chaque residu au centre de masse de la proteine.
	Distance locale moyenne.
	'''

	l_res = d_ref[d_ref["liste_conformations"][0]]["liste_n_residus"]

	d_ref["enfRes_mean"] = list()
	d_ref["enfRes_sd"] = list()
	d_ref["list_enfRes"] = list()

	for i in range(len(l_res)) :
		l_enf = list()
		for conf in d_conf["liste_conformations"] :
			l_enf.append(d_conf[conf]["enfouissement"][i])
		d_ref["list_enfRes"].append(l_enf)
		d_ref["enfRes_mean"].append(moyenne(l_enf))
		d_ref["enfRes_sd"].append(ecart_type(l_enf))




#-----------------------------------------------------------------------
# PARTIE 5 : RAYON DE GIRATION
#-----------------------------------------------------------------------
'''
Le rayon de giration est calcule comme etant la distance maximale entre 
un residu et le centre de masse d'une conformation.
C'est la distance associee au residu qui a l'enfouissement miniimal.
'''


def rayonGiration(l_dict) :
	'''
	Creation de la cle ratio_giration qui a pour valeur le rapport entre le rayon de Giration des conformations
	et celui de la reference.
	'''
	print "Calcul des Rayons de Giration"
	d_ref = l_dict[0]
	d_conf = l_dict[1]
		
	rayon = __maxDistance(d_conf)
	ref = __maxDistance(d_ref)[0]
	
	d_conf["ratio_giration"] = [x/ref for x in rayon]


def __maxDistance(d_prot) :
	'''
	Recherche du Rayon de Giration parmi toutes les distances calculees precedemment.
	Creation de la cle rayonGiration qui a pour valeur la distance entre entre le CdM
	et le residu le plus eloigne du CdM, pour chaque conformation
	'''
	d_prot["rayonGiration"] = list()
	for conf in d_prot["liste_conformations"] :
		d_prot["rayonGiration"].append(max(d_prot[conf]["enfouissement"]))
	return d_prot["rayonGiration"]




#-----------------------------------------------------------------------
def corEnfouissementFlexibilite_conf(d_conf) :
	'''
	Calcul de la correlation entre l'enfouissement des residus et la flexibilite des regions pour chaque conformation.
	L'enfouissement est inversement proportionnel a la distance des residus aux CdM.
	La flexibilite augmente avec le RMSD entre le residu et sa position dans la reference.

	'''
	d_conf["corEnfFlexi_conf"] = [list(), list()] # correlation et pvaleur
	for conf in d_conf["liste_conformations"] :
		if(d_conf[conf]["RMSD"] == [0] * len(d_conf[conf]["RMSD"])) :
			cor = [1,0]
		else :
			cor = pearsonr(d_conf[conf]["enfouissement"],d_conf[conf]["RMSD"])
		d_conf["corEnfFlexi_conf"][0].append(cor[0])
		d_conf["corEnfFlexi_conf"][1].append(cor[1])
	return d_conf["corEnfFlexi_conf"]


def corEnfouissementFlexibilite_res(d_ref) :
	#Calcul, pour chaque residu de chaque conformation, de la correlation entre l'enfouissement du residu et sa flexibilite.
	d_ref["corEnfFlexi_ref"] = [list(), list()] # correlation et pvaleur
	for i in range(len(d_ref["list_enfRes"])) :
		cor = pearsonr(d_ref["list_enfRes"][i],d_ref["list_RMSDres"][i])
		d_ref["corEnfFlexi_ref"][0].append(cor[0])
		d_ref["corEnfFlexi_ref"][1].append(cor[1])
	return d_ref["corEnfFlexi_ref"]




#-----------------------------------------------------------------------
# FONCTIONS DE CALCUL STATISTIQUE
#-----------------------------------------------------------------------

def moyenne(liste) :
	return sum(liste)/len(liste)

def variance(liste) :
	n = len(liste)
	if (n != 0) :
		m = moyenne(liste)**2
		s = sum([x**2 for x in liste])
		return s/n-m
	print "ERREUR : Dictionnaire sans conformation associee"
	sys.exit(1)

def ecart_type(liste) :
	return sqrt(variance(liste))