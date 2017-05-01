#!/usr/bin/env python2

# usage : python2 Barstar.py start_prot_only.pdb md_prot_only_skip100.pdb CA

import sys, os
from math import sqrt
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import numpy as np
import csv

#PLAN
	#Changements conformationnels globaux

		#Calcul rayon de Giration (distance entre CdM et residu le plus eloigne du CdM)
		#et RMSD de toutes les conformations par rapport a la ref

		#Representation de la variation Giration/RMSD, Giration/Temps, RMSD/Temps.

	#Changements conformationnels locaux
	
		#RMSD de chaque residu de chaque conformation par rapport a sa position dans la ref
			#Region flexible = region dont residus ont un grand RMSD

		#Enfouissement : Calcul distance entre CdM de chaque residu et CdM du centre de la prot

		#Comparer enfouissement des residus et RMSD des residus
			#Regions flexibles enfouies ou en surface ?

		#Calcul RMSD et Enfouissement moyen pour chaque residu de chaque conformation vs reference
			#Representation Graphique

#-----------------------------------------------------------------------
def	__usage(arguments) :
	# Verifie que l'utilisateur a fourni le bon nombre d'arguments
	# et que le 3e argument correspond bien a la methode de calcul du Centre de Masse. (CA ou all)

	if len(arguments) != 4 :
		__erreurMes()
	if ((arguments[3] != 'CA') & (arguments[3] != 'all')) :
		__erreurMes()

def __erreurMes() :
	print "ATTENTION : Mauvais usage de Barstar \nArguments necessaires :\t<fileRef.pdb> <fileConf.pdb> <nom_atome>\n"
	print "fileRef.pdb :\t\tfichier pdb contenant la conformation de reference"
	print "fileConf.pdb :\t\tfichiers pdb contenant les conformations issues de la dynamique moleculaire"
	print "nom_atome :\t\tCA ou all (Methode de calcul du centre de masse)"
	sys.exit(0)

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
def dictionnaire() :
	global argv
	
	# Verification : les fichiers ont ete fournis dans le bon ordre
	__usage(sys.argv)
	
	# Transformation des fichiers pdb en dictionnaires
	dico1 = __parsePDBMultiConf(sys.argv[1])
	dico2 = __parsePDBMultiConf(sys.argv[2])

	# Verification ordre des fichiers
	if (__verifData(dico1, dico2) == 1) :
		__inversion(dico1, dico2)
	
	return [dico1, dico2]


def __parsePDBMultiConf(infile) :
	"""
	Cette fonction permet de parser un fichier pdb. Conversion en 
	dictionnaire directement utilisable dans le script.
	"""

	try:
		f = open(infile, "r")
		lines = f.readlines()
		f.close()
	except:
		print("Erreur chargement fichier : verifiez existence fichier et relancer\n")
		sys.exit(0) # Arret execution programme

	print "Parsing de "+infile
	d_PDB = dict()
	d_PDB["liste_conformations"] = list()

	for element in lines :		
		if element[0:5] == "MODEL" :
			conformation = element[10:14]
			d_PDB["liste_conformations"].append(conformation)
			d_PDB[conformation] = dict()
			d_PDB[conformation]["liste_n_residus"] = list()
			d_PDB[conformation]["liste_seq_residus"] = list()
			d_PDB[conformation]["liste_residus"] = list()
        
		if element[0:4] == "ATOM" :
			residus = element[17:20]
			n_res = element[23:26].strip() 
			atome = element[13:16].strip() 
			identifiant = element[7:11].strip() 
			x = float(element[30:38])
			y = float(element[38:46])
			z = float(element[46:54])
			
			if not n_res in d_PDB[conformation]["liste_n_residus"] :
				d_PDB[conformation]["liste_n_residus"].append(n_res)
				d_PDB[conformation]["liste_seq_residus"].append(residus)
				
				if not residus in d_PDB[conformation]["liste_residus"] :
					d_PDB[conformation]["liste_residus"].append(residus)
				
				d_PDB[conformation][n_res] = dict()
				d_PDB[conformation][n_res]["liste_atomes"] = list()
			
			d_PDB[conformation][n_res]["liste_atomes"].append(atome)
			d_PDB[conformation][n_res][atome] = dict()
			d_PDB[conformation][n_res][atome]["x"] = x
			d_PDB[conformation][n_res][atome]["y"] = y
			d_PDB[conformation][n_res][atome]["z"] = z
			d_PDB[conformation][n_res][atome]["id"] = identifiant

	return d_PDB


#-----------------------------------------------------------------------
def __verifData(dico1, dico2) :
	'''
	Verification des fichiers : 
		- dico1 contient la conformation de reference
		- dico2 contient l'ensemble des conformations pour tous les temps consideres
	S'ils ont ete renseignes dans le bon ordre, la fonction retourne 0
	sinon elle retourne 1 : necessite d'inverser les dictionnaires
	'''
	if (dico1 == dico2) :
		print "Erreur : les dictionnaires sont identiques"
		sys.exit(0)
	if ((len(dico1["liste_conformations"]) != 1) | (len(dico2["liste_conformations"]) == 1)) :
		return 1

def __inversion(dico1, dico2) :
		d_tmp = dico2
		dico2 = dico1
		dico1 = d_tmp





#-----------------------------------------------------------------------
# PARTIE 2 CALCUL DU CENTRE DE MASSE
#-----------------------------------------------------------------------

def centreMasseCalc(l_dict) :
	print "Calcul des Centres de Masse des Residus et des proteines"
	__centreMasseResidus(l_dict[0])
	__centreMasseResidus(l_dict[1])
	__centreMasseProteine(l_dict[0])
	__centreMasseProteine(l_dict[1])


def choixMeth() :
	'''
	Deux methodes de calcul sont proposees :
	- soit en utilisant la position des carbones alpha (hyp. representatif des residus)
	- soit en calculant la position des residus a partir du centre de masse de l'ensemble des atomes les consituant
	'''
	global argv
	
	arg = sys.argv[3]
	
	if arg == "CA" :
		return "CM_CA"
	else :
		return "CM_moyAll"

	
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
        


#-----------------------------------------------------------------------
# REPRESENTATION DES DONNEES
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
	#Calcul de la correlation entre l'enfouissement des residus et leur flexibilite en fonction du residu.
	d_ref["corEnfFlexi_ref"] = [list(), list()] # correlation et pvaleur
	for i in range(len(d_ref["list_enfRes"])) :
		cor = pearsonr(d_ref["list_enfRes"][i],d_ref["list_RMSDres"][i])
		d_ref["corEnfFlexi_ref"][0].append(cor[0])
		d_ref["corEnfFlexi_ref"][1].append(cor[1])
	return d_ref["corEnfFlexi_ref"]


#-----------------------------------------------------------------------
# PARTIE 6 : ECRITURE DES RESULTATS
#-----------------------------------------------------------------------

def ecriture(l_dict,methode) :

	print "Ecriture des resultats dans les fichiers de sortie 'res_barstar_globaux_"+methode+".csv' et 'res_barstar_locaux_"+methode+".csv'"
	output1 = __verificationfFichier("res_barstar_globaux_"+methode+".csv") # nom de fichier par defaut, mais on ne veut pas ecraser des resultats precedents
	output2 = __verificationfFichier("res_barstar_locaux_"+methode+".csv")
	
	x = input("Indiquez le nombre de decimales souhaitees :\t")
	__outputGlobaux(output1, l_dict[0], l_dict[1], x)
	__outputLocaux(output2, l_dict[0], x)



def __verificationfFichier(output) :
	#Si le fichier existe deja, demande si l'utilisateur veut le remplacer ou non

	if (os.path.exists(output)) :
		decision = raw_input("Fichier de sortie "+str(output)+" deja existant : Voulez vous l'ecraser ? O/N\n")
	
		while ((decision != 'O') & (decision != 'o') & (decision != 'N')  & (decision != 'n')) :
			print "Erreur : Repondre O pour oui ou N pour non" 
			decision = raw_input("Fichier de sortie deja existant : Voulez vous l'ecraser ? O/N\n")
	
		if ((decision == 'N') | (decision == 'n')) :
			while os.path.exists(output) :
				output = raw_input("Nom de fichier de sortie deja existant : entrez un nouveau nom de fichier : \n")
	return output



def __outputGlobaux(output, d_ref, d_conf, x) :
	#Ecriture des resultats dans un tableau excel

	try:
		with open(output, "w") as f:
			fieldnames = ['Conf', 'Rayon Giration','Distance','ecart-type Distance','RMSD','ecart-type RMSD','Ratio Giration','Correlation','p-value Correlation']
			writer = csv.DictWriter(f, fieldnames=fieldnames)

			writer.writeheader()
			writer.writerow({'Conf':'REF','Rayon Giration': round(d_ref["rayonGiration"][0],x), 'Distance': 0})

			for i in range(len(d_conf["liste_conformations"])) :

				num = d_conf["liste_conformations"][i].strip()
				rayonG = d_conf["rayonGiration"][i]
				d_moy = d_conf["distance_moy"][i]
				d_sd = d_conf["distance_sd"][i]
				rmsd = d_conf["RMSDmoy"][i]
				rmsd_sd = d_conf["RMSDmoy_sd"][i]
				ratio_gir = d_conf["ratio_giration"][i]
				cor = d_conf["corEnfFlexi_conf"][0][i]
				pvalue = d_conf["corEnfFlexi_conf"][1][i]
				writer.writerow({'Conf': num, 'Rayon Giration': round(rayonG,x), 'Distance': round(d_moy,x), 'ecart-type Distance': round(d_sd,x), 'RMSD': round(rmsd,x), 'ecart-type RMSD': round(rmsd_sd,x), 'Ratio Giration': round(ratio_gir,x), 'Correlation': round(cor,x),'p-value Correlation': round(pvalue,x)})	

	except:
		print("Erreur chargement fichier"+output+"\n")
		sys.exit(0)


def __outputLocaux(output, d_ref, x) :
	#Ecriture des resultats dans un tableau excel

	try:
		with open(output, "w") as f:
			fieldnames = ['Residus', 'Nom','RMSD', 'ecart-type RMSD','Distance residu/CdM','ecart-type Distance residu/CdM']
			writer = csv.DictWriter(f, fieldnames=fieldnames)
		
			writer.writeheader()

			conf_ref = d_ref["liste_conformations"][0]
			lres = d_ref[conf_ref]["liste_n_residus"]
		
			for i in range(len(lres)) :
				res = lres[i]
				nom = d_ref[d_ref["liste_conformations"][0]]["liste_seq_residus"][i]
				rmsd = d_ref["RMSDres_mean"][i]
				rmsd_sd = d_ref["RMSDres_sd"][i]
				d = d_ref["enfRes_mean"][i]
				d_sd = d_ref["enfRes_sd"][i]
			
				writer.writerow({'Residus': res, 'Nom': nom, 'RMSD': round(rmsd,x),'ecart-type RMSD': round(rmsd_sd,x), 'Distance residu/CdM': round(d,x), 'ecart-type Distance residu/CdM' : round(d_sd,x)})	

				
	except:
		print("Erreur chargement fichier"+output+"\n")
		sys.exit(0)


#-----------------------------------------------------------------------
# GRAPHIQUES
#-----------------------------------------------------------------------

def plotRes(l_dict) :
	plotGlobal(l_dict[1])
	plotLocal(l_dict[0]) 

def plotGlobal(d_conf) :
	plotGiration(d_conf)
	plotDistance(d_conf)
	plotGlobalRMSD(d_conf)
	plotFlexibiteEnfouissement(d_conf)
	
def plotGiration(d_conf) :
	# plt.subplot(211)
	plt.title('Rayon de Giration en fonction des conformations')
	plt.plot(d_conf["rayonGiration"])
	plt.axhline(y=d_conf["rayonGiration"][0],ls='--',color='black')
	# plt.plot(d_conf["rayonGiration"][0], color="black")
	plt.xlabel('Conformations')
	plt.ylabel('Rayon de Giration')


	# plt.subplot(212)
	# plt.title('Ratio du Rayon de Giration des conformations\n sur celui de la reference en fonction des conformations')
	# plt.plot(d_conf["ratio_giration"])
	# plt.xlabel('Conformation')
	# plt.ylabel('Ratio Rayon Giration\n Conformation/Reference')

	# plt.tight_layout()
	plt.show()

def plotDistance(d_conf) :
	# plt.subplot(311)
	# plt.plot(d_conf["distance_moy"])
	# plt.xlabel('Conformation')
	# plt.ylabel('Distance moyenne')
	
	# plt.subplot(312)
	# plt.plot(d_conf["distance_sd"])
	# plt.xlabel('Conformation')
	# plt.ylabel('Ecart type')
	
	# plt.subplot(313)
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
	
	plt.show()

def plotGlobalRMSD(d_conf) :
	# plt.subplot(311)
	# plt.plot(d_conf["RMSDmoy"])
	# plt.xlabel('Conformation')
	# plt.ylabel('RMSD moyen')

	# plt.subplot(312)
	# plt.plot(d_conf["RMSDmoy_sd"])
	# plt.xlabel('Conformation')
	# plt.ylabel('Ecart type')

	# plt.subplot(313)
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
	plt.show()

#-----------------------------------------------------------------------
def plotLocal(d_ref) :
	#Appel des fonctions de plot pour les graphes de l'analyse locale
	plotDistanceLocal(d_ref)
	plotRMSDLocal(d_ref)
	plotFlexibiteEnfouissement_residus(d_ref)
	

def plotDistanceLocal(d_ref) :
	# plt.subplot(311)
	# plt.plot(d_ref["enfRes_mean"])
	# plt.xlabel('Conformation')
	# plt.ylabel('Distance')
	
	# plt.subplot(312)
	# plt.plot(d_ref["enfRes_sd"])
	# plt.xlabel('Conformation')
	# plt.ylabel('ecart type')
	
	# plt.subplot(313)
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
	plt.show()



def plotRMSDLocal(d_ref) :
	# plt.subplot(311)
	# plt.plot(d_ref["RMSDres_mean"])
	# plt.xlabel('Conformation')
	# plt.ylabel('RMSD')
	
	# plt.subplot(312)
	# plt.plot(d_ref["RMSDres_sd"])
	# plt.xlabel('Conformation')
	# plt.ylabel('ecart type')

	# plt.subplot(313)
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
	plt.show()

#-----------------------------------------------------------------------
# MAIN

if __name__ == '__main__':
	liste_dictionnaire = dictionnaire()					# Creation d'une liste : [dictionnaire_des_conformations,dictionnaire_de_la_reference]
	centreMasse = choixMeth()							# "CM_CA" ou "CM_moyall" en fonction de la methode de calcul choisie
	conformation_analysis(liste_dictionnaire) 			# calcul des variables d'interet et ajout dans les dictionnaires 
	ecriture(liste_dictionnaire,centreMasse) 			# ecriture des resultats dans les tableaux excel de sortie, en fonction de la methode de calcul des CdM
	plotRes(liste_dictionnaire) 						# representations graphiques

	print ""
