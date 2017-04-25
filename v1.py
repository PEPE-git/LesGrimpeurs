

# usage : python v1.py start_prot_only.pdb md_prot_only_skip10.pdb

import sys, os
from math import sqrt
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import numpy as np

#PLAN
	#Changements conformationnels globaux

		#Calcul rayon de Giration (distance entre CdM et residu le plus eloigne du CdM)
		#et RMSD de toutes les conformations par rapport a la ref

		#Representation de la variation Giration/RMSD, Giration/Temps, RMSD/Temps.

	#Changements conformationnels locaux
	
		#RMSD de chaque residu de chaque conformation par rapport à sa position dans la ref
			#Region flexible = region dont residus ont un grand RMSD

		#Enfouissement : Calcul distance entre CdM de chaque residu et CdM du centre de la prot

		#Comparer enfouissement des residus et RMSD des residus
			#Regions flexibles enfouies ou en surface ?

		#Calcul RMSD et Enfouissement moyen pour chaque residu de chaque conformation vs reference
			#Representation Graphique

#-----------------------------------------------------------------------
def	usage(arguments) :
	
	if len(arguments) != 4 :
		print "Mauvais usage de Barstar \n Arguments necessaires  <fileRef.pdb> <fileConf.pdb> <nom_atome>"
		print "fileRef.pdb : fichier pdb contenant la conformation de reference"
		print "fileConf.pdb : fichiers pdb contenant les conformations issues de la dynamique moleculaire"
		print "nom_atome : CA ou all --> calcul du RMSD"
		sys.exit(0)
	if ((arguments[3] != 'CA') & (arguments[3] != 'all')) :
		print "Mauvais usage de Barstar \n Arguments necessaires  <fileRef.pdb> <fileConf.pdb> <nom_atome>"
		print "fileRef.pdb : fichier pdb contenant la conformation de reference"
		print "fileConf.pdb : fichiers pdb contenant les conformations issues de la dynamique moleculaire"
		print "nom_atome : CA ou all --> calcul du RMSD"
		sys.exit(0)


#-----------------------------------------------------------------------
def parsePDBMultiConf(infile) :
	"""
		Cette fonction permet de parser un fichier pdb. Conversion en 
		dictionnaire directement utilisable dans le script.
	"""

	try:
		f = open(infile, "r")
		lines = f.readlines()
		f.close()
	except:
		print("Erreur chargement fichier : verifiez existance fichier et relancer\n")
		sys.exit(0) # Arret execution programme

	d_PDB = dict()
	d_PDB["liste_conformations"] = list()

	for element in lines :		
		if element[0:5] == "MODEL" :
			conformation = element[13:17]
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
# Choix methode de calcul

def choixMeth(arg) :
	if arg == "CA" :
		return "CM_CA"
	else :
		return "CM_moyAll"

#-----------------------------------------------------------------------
def verifData(dico1, dico2) :
	# Fonction qui permet de verifier si les fichiers ont ete fournis dans le bon sens
	# sinon retourne 1 --> indique la necessite d'inverser les dictionnaires correspondants
	if (dico1 == dico2) :
		print "Erreur : les dictionnaires sont identiques"
		sys.exit(0)
	if ((len(dico1["liste_conformations"]) != 1) | (len(dico2["liste_conformations"]) == 1)) :
		return 1



#~ class DataException(Exception):
    #~ def __init__(self,raison):
        #~ self.raison = raison
    
    #~ def __str__(self):
        #~ return self.raison

#~ def exactData(dico1,dico2):
    
    #~ if :
		#~ raise DataException("Erreur : les conformations n'ont pas toutes les memes residus")
	
	#~ elif (dico1 == dico2) :
		#~ raise DataException("Erreur : les dictionnaires sont identiques")
    #~ else :
		#~ # Fonction qui permet de verifier si les fichiers ont ete fournis dans le bon sens
		#~ # sinon retourne 1 --> indique la necessite d'inverser les dictionnaires correspondants
		#~ if ((len(dico1["liste_conformations"]) != 1) | (len(dico2["liste_conformations"]) == 1)) :
			#~ return 1

def inversion(dico1, dico2) :
		d_tmp = dico2
		dico2 = dico1
		dico1 = d_tmp

#-----------------------------------------------------------------------
def timeList (infile): #Recupere et met dans un liste le temps des differentes conformations
    time = []
    file = open(infile,"r")
    lines = file.readlines()
    for line in lines :
        if line[0:5] == "TITLE":
            timet=line[65:80].strip()
            time.append(timet)
    #print temps
    file.close()
    return(time)
#-----------------------------------------------------------------------
# CENTRE DE MASSE

def centreMasseProteine(d_prot, centreMasse) :
	'''
		Fonction qui calcule le centre de masse des proteines. Cela
		correspond a la position moyenne des residus la constituant.
		
		@param : indique la methode du calcul de centre de masse CM_CA ou CM_moyAll
	'''
	
	l_CMx = [] # stockage des coordonnees des CM des residus sur l'axe x
	l_CMy = [] # stockage des coordonnees des CM des residus sur l'axe y
	l_CMz = [] # stockage des coordonnees des CM des residus sur l'axe z
		
	for conf in d_prot["liste_conformations"] :
		l_CM = [0,0,0] # liste de la somme des coordonnees des CM
		cpt = 0 # cpt du nombre de residus
		
		for resid in d_prot[conf]["liste_n_residus"] :
			cpt += 1
			l_CM[0] += d_prot[conf][resid][centreMasse]["x"]
			l_CM[1] += d_prot[conf][resid][centreMasse]["y"]
			l_CM[2] += d_prot[conf][resid][centreMasse]["z"]
		l_CMx.append(l_CM[0]/cpt)
		l_CMy.append(l_CM[1]/cpt)
		l_CMz.append(l_CM[2]/cpt)
	d_prot["liste_CM_x"] = l_CMx
	d_prot["liste_CM_y"] = l_CMy
	d_prot["liste_CM_z"] = l_CMz
	
	return d_prot

def centreMasseResidus(d_prot, meth) :
	'''
		Le centre d'un masse d'un residus peut etre considere comme etant :
			- la position moyenne de ses atomes
			- la position de ses carbones alpha
	'''
	if meth == "CM_CA" :
		return centreMasseResCa(d_prot)
	return centreMasseResAll(d_prot)


def centreMasseResCa(d_prot) :
	'''		
		Fonction qui calcule le centre de masse des residus
		Le centre de masse d'un residus correspond a la position de son
		carbone alpha
		
		@param diconnaire de la proteine avec une ou plusieurs 
		conformations
		@return dictionnaire de la proteine avec la valeur des centres
		de masse "CA" des residus ajoutes
	'''
	
	# Pour toutes les conformations de la proteine
	for conf in d_prot["liste_conformations"] :
		# Pour tous les residus de chaque conformation de la proteine
		
		d_prot[conf]["CM_res"] = dict()
		d_prot[conf]["CM_res"]["x"] = list()
		d_prot[conf]["CM_res"]["y"] = list()
		d_prot[conf]["CM_res"]["z"] = list()
		for resid in d_prot[conf]["liste_n_residus"] :
			d_prot[conf][resid]["CM_CA"] = dict()
			d_prot[conf][resid]["CM_CA"]["x"] = d_prot[conf][resid]["CA"]["x"]
			d_prot[conf][resid]["CM_CA"]["y"] = d_prot[conf][resid]["CA"]["y"]
			d_prot[conf][resid]["CM_CA"]["z"] = d_prot[conf][resid]["CA"]["z"]
			
			d_prot[conf]["CM_res"]["x"].append(d_prot[conf][resid]["CA"]["x"])
			d_prot[conf]["CM_res"]["y"].append(d_prot[conf][resid]["CA"]["y"])
			d_prot[conf]["CM_res"]["z"].append(d_prot[conf][resid]["CA"]["z"])
	return d_prot

def centreMasseResAll(d_prot) :
	'''		
		Fonction qui calcule le centre de masse des residus
		Le centre de masse d'un residus correspond a la position moyenne
		des atomes le constituant
		
		@param dictionnaire de la proteine avec une ou plusieurs 
		conformations
		@return dictionnaire de la proteine avec la valeur des centres
		de masse "all" des residus ajoutes
	'''
	
	# Pour toutes les conformations de la proteine
	for conf in d_prot["liste_conformations"] :
		# Pour tous les residus de chaque conformation de la proteine
		for resid in d_PDB[conf]["liste_n_residus"] :
			cpt = 0
			sumx = 0
			sumy = 0
			sumz = 0
			# Pour tous les atomes de chaque residus

			d_prot[conf]["CM_res"]["x"] = list()
			d_prot[conf]["CM_res"]["y"] = list()
			d_prot[conf]["CM_res"]["z"] = list()

			for atom in d_PDB[conformation][resid]["liste_atomes"] :
				cpt += 1
				sumx += d_PDB[conf][resid][atom]["x"]
				sumy += d_PDB[conf][resid][atom]["y"]
				sumz += d_PDB[conf][resid][atom]["z"]

			d_PDB[conformation][resid]["CM_moyAll"]["x"] = sumx/cpt
			d_PDB[conformation][resid]["CM_moyAll"]["y"] = sumy/cpt
			d_PDB[conformation][resid]["CM_moyAll"]["z"] = sumz/cpt

			d_prot[conf]["CM_res"]["x"].append(sumx/cpt)
			d_prot[conf]["CM_res"]["y"].append(sumy/cpt)
			d_prot[conf]["CM_res"]["z"].append(sumz/cpt)

	return d_prot


#-----------------------------------------------------------------------
# RMSD

'''
	Le RMSD est considere comme la distance moyenne entre les residus
	
'''

def RMSDresidus(d_ref, d_conf, centreMasse) :
	
	REF = d_ref["liste_conformations"][0] # il n'y a qu'une conformation : celle de reference
	
	for conf in d_conf["liste_conformations"] :
		d_conf[conf]["RMSD"] = list()
		for resid in d_conf[conf]["liste_n_residus"] :
			xref = d_ref[REF][resid][centreMasse]["x"]
			yref = d_ref[REF][resid][centreMasse]["y"]
			zref = d_ref[REF][resid][centreMasse]["z"]

			xconf = d_conf[conf][resid][centreMasse]["x"]
			yconf = d_conf[conf][resid][centreMasse]["y"]
			zconf = d_conf[conf][resid][centreMasse]["z"]
			
			d_conf[conf]["RMSD"].append((xref-xconf)**2 + (yref-yconf)**2 + (zref-zconf)**2)

	return d_conf


#-----------------------------------------------------------------------
# DISTANCE --> Calcul de la distance de chaque residu au centre de Masse
# de la conformation  --> enfouissement des residus au sein de la conformation

def distance(d_prot, centreMasse) :
	
	# pour chaque conformation, son centre de masse :
	xconf = d_prot["liste_CM_x"]
	yconf = d_prot["liste_CM_y"]
	zconf = d_prot["liste_CM_z"]
	
	for conf in d_prot["liste_conformations"] :
		l_dist = []
		
		for resid in d_prot[conf]["liste_n_residus"] :
			i=0
			xres = d_prot[conf][resid][centreMasse]["x"]
			yres = d_prot[conf][resid][centreMasse]["y"]
			zres = d_prot[conf][resid][centreMasse]["z"]
			
			distance = sqrt((xres-xconf[i])**2 + (yres-yconf[i])**2 + (zres-zconf[i])**2)
			l_dist.append(distance)
			
			i += 1
		d_prot[conf]["enfouissement"] = l_dist

	return d_prot

#-----------------------------------------------------------------------
# RAYON DE GIRATION --> distance maximale entre un residus et le centre de masse
# d'une conformation <-> residus a l'enfouissement maximal

def rayonGiration(d_prot) :
	d_prot["rayonGiration"] = list()
	for conf in d_prot["liste_conformations"] :
		d_prot["rayonGiration"].append(max(d_prot[conf]["enfouissement"]))
	return d_prot["rayonGiration"]


#-----------------------------------------------------------------------
# REPRESENTATION DES DONNEES

def corEnfouissementFlexibilite(d_prot) :
	'''
		Fonction qui permet de visualiser la correlation entre l'enfouissement
		des residus et la flexibilite des regions
		La flexibilite augmente avec la distance des residus aux CdM
		--> calcul de la correlation
	'''
	d_prot["corEnfFlexi"] = list()
	for conf in d_prot["liste_conformations"] :
		d_prot["corEnfFlexi"].append(pearsonr(d_prot[conf]["enfouissement"],d_prot[conf]["RMSD"]))
	return d_prot["corEnfFlexi"]

def plotRMSD_Giration(listRMSD, listGiration):
	y=np.array(listRMSD)
	x=np.array(listGiration)
	plt.scatter(x,y,c='red')
	axes = plt.gca()
	axes.set_xlim(-30, 2100)
	axes.set_ylim(-1,25)
	plt.title('RMSD en fonction de la Giration')
	plt.legend(['RMSD','Giration'])
	plt.show()

def plotRMSD_Temps(listRMSD, listTemps):

def plotGiration_Temps(listGiration, listTemps):

#-----------------------------------------------------------------------
# Ecriture des resultats

def verificationfFichier(output) :
	if (os.path.exists(output)) :
		decision = raw_input("Fichier de sortie deja existant : Voulez vous l'ecraser ? O/N")
		while (decision != "O" | decision != "o" | decision != "N"  | decision != "n") :
			print "Erreur : Repondre O pour oui ou N pour non" 
			decision = raw_input("Fichier de sortie deja existant : Voulez vous l'ecraser ? O/N")
		if ((decision == 'N') | (decision == 'n')) :
			while os.path.exists(output) :
				output = raw_input("Nom de fichier de sortie deja existant : entrez un nouveau nom de fichier : ")
	return output

def ecritureResultat(output) :
	output = verificationfFichier(output)
	try:
		f = open(output, "w")
		f.write()
		f.close()
	except:
		print("Erreur chargement fichier\n")
		sys.exit(0)

def outputGlobaux() :
	outfile = "res_barstar_globaux.txt"
	ecritureResultat(outfile)

def outputLocaux() :
	outfile = "res_barstar_locaux.txt"
	ecritureResultat(outfile)


#-----------------------------------------------------------------------
# MAIN

if __name__ == '__main__':
	usage(sys.argv)
	
	# Extraction de la conformation de reference
	d_ref = parsePDBMultiConf(sys.argv[1])
	
	# Extraction des autres conformations
	d_conf = parsePDBMultiConf(sys.argv[2])
	
	# Verification ordre des fichiers
	if (verifData(d_ref, d_conf) == 1) :
		inversion(d_ref, d_conf)

	# Il faut verifier que toutes les conformations ont le meme nombre de residus que 
	# la conformation de reference
	# stocke ce nombre de residus --> evitera d'avoir a mettre un compteur dans les 
	# fonction
	
	# idem pour les residus : meme nombre d'atome
	
	# pour cela comparaison des listes de residus et listes atomes !
	
	# --> sinon exit avec erreur
	
	# Methode de calcul du centre de masse
	meth = choixMeth(sys.argv[3]) # variable globale ??
	
	# Centre Masse Residus
	d_ref = centreMasseResidus(d_ref, meth)
	d_conf = centreMasseResidus(d_conf, meth)
	
	# Centre Masse Proteine
	d_ref = centreMasseProteine(d_ref, meth)
	d_conf = centreMasseProteine(d_conf, meth)
	
	# RMSD
	d_conf = RMSDresidus(d_ref, d_conf, meth)
	
	# distance
	d_conf = distance(d_conf, meth)
	d_ref = distance(d_ref, meth)

	# rayon de Giration :
	rayGiration = rayonGiration(d_conf)
	rayGirationRef = rayonGiration(d_ref)
	#~ print d_conf["rayonGiration"]
	#~ rayGiration = [i/rayGirationRef[0] for i in rayGiration]
	#~ plt.plot(rayGiration)
	
	# pour chaque conformation --> calcul de la correlation entre
	correlation = corEnfouissementFlexibilite(d_conf)
	#~ print correlation
	plt.plot(correlation)
	plt.show()
