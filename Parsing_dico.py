#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""

import sys, os

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