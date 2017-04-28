

# usage : python v1.py start_prot_only.pdb md_prot_only_skip100.pdb CA

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
	
		#RMSD de chaque residu de chaque conformation par rapport a sa position dans la ref
			#Region flexible = region dont residus ont un grand RMSD

		#Enfouissement : Calcul distance entre CdM de chaque residu et CdM du centre de la prot

		#Comparer enfouissement des residus et RMSD des residus
			#Regions flexibles enfouies ou en surface ?

		#Calcul RMSD et Enfouissement moyen pour chaque residu de chaque conformation vs reference
			#Representation Graphique

#-----------------------------------------------------------------------
def	usage(arguments) :
	if len(arguments) != 4 :
		erreurMes()
	if ((arguments[3] != 'CA') & (arguments[3] != 'all')) :
		erreurMes()

def erreurMes() :
	print "ATTENTION : Mauvais usage de Barstar \nArguments necessaires :\t<fileRef.pdb> <fileConf.pdb> <nom_atome>\n"
	print "fileRef.pdb :\t\tfichier pdb contenant la conformation de reference"
	print "fileConf.pdb :\t\tfichiers pdb contenant les conformations issues de la dynamique moleculaire"
	print "nom_atome :\t\tCA ou all (Methode de calcul du centre de masse)"
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

def centreMasseProteine(d_prot) :
	'''
		Fonction qui calcule le centre de masse des proteines. Cela
		correspond a la position moyenne des residus la constituant.
		
		@param : indique la methode du calcul de centre de masse CM_CA ou CM_moyAll
	'''
	
	global centreMasse
	
	l_CM = []
	l_CMx = [] # stockage des coordonnees des CM des residus sur l'axe x
	l_CMy = [] # stockage des coordonnees des CM des residus sur l'axe y
	l_CMz = [] # stockage des coordonnees des CM des residus sur l'axe z
		
	for conf in d_prot["liste_conformations"] :
		l_CM = [list(),list(),list()] # liste de la somme des coordonnees des CM
		cpt = 0 # cpt du nombre de residus
		
		for resid in d_prot[conf]["liste_n_residus"] :
			cpt += 1
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

def centreMasseResidus(d_prot) :
	'''
	Le centre d'un masse d'un residus peut etre considere comme etant :
		- la position moyenne de ses atomes
		- la position de ses carbones alpha
	'''
	global centreMasse
	
	if centreMasse == "CM_CA" :
		return centreMasseResCa(d_prot)
	return centreMasseResAll(d_prot)


def centreMasseResCa(d_prot) :
	'''		
	Fonction qui calcule le centre de masse des residus
	Le centre de masse d'un residus correspond a la position de son
	carbone alpha
	
	@param diconnaire de la proteine avec une ou plusieurs conformations
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
	
	@param dictionnaire de la proteine avec une ou plusieurs conformations
	@return dictionnaire de la proteine avec la valeur des centres de masse 
	"all" des residus ajoutes
	'''
	
	# Pour toutes les conformations de la proteine
	for conf in d_prot["liste_conformations"] :
		# Pour tous les residus de chaque conformation de la proteine
		for resid in d_PDB[conf]["liste_n_residus"] :
			lcoord = list()
			lcoord[0] = list() # stockage coordonnees x
			lcoord[1] = list() # stockage coordonnees y
			lcoord[2] = list() # stockage coordonnees z
			# Pour tous les atomes de chaque residus

			d_prot[conf]["CM_res"]["x"] = list()
			d_prot[conf]["CM_res"]["y"] = list()
			d_prot[conf]["CM_res"]["z"] = list()

			for atom in d_PDB[conformation][resid]["liste_atomes"] :
				lcoord[0].append(d_PDB[conf][resid][atom]["x"])
				lcoord[1].append(d_PDB[conf][resid][atom]["y"])
				lcoord[2].append(d_PDB[conf][resid][atom]["z"])

			xmoy = moyenne(lcoord[0])
			ymoy = moyenne(lcoord[1])
			zmoy = moyenne(lcoord[2])

			d_PDB[conformation][resid]["CM_moyAll"]["x"] = xmoy
			d_PDB[conformation][resid]["CM_moyAll"]["y"] = ymoy
			d_PDB[conformation][resid]["CM_moyAll"]["z"] = zmoy

			d_prot[conf]["CM_res"]["x"].append(xmoy)
			d_prot[conf]["CM_res"]["y"].append(ymoy)
			d_prot[conf]["CM_res"]["z"].append(zmoy)

	return d_prot


#-----------------------------------------------------------------------
# RMSD

'''
	Le RMSD est considere comme la distance moyenne entre les residus
	
'''

def RMSD(d_ref, d_conf) :
	RMSDresidus(d_ref, d_conf)
	
	RMSDconf(d_ref)
	RMSDconf(d_conf)
	
	RMSDres(d_ref, d_conf)
	
	RMSDratio(d_ref, d_conf)


def RMSDratio(dico1, dico2) :
	ref = dico1["RMSDres_mean"][0]
	dico1["ratio_RMSD"] = [x/ref for x in dico2["RMSDres_mean"]]
	return dico1["ratio_RMSD"]
	
	
	
	
def RMSDresidus(d_ref, d_conf) :
	
	global centreMasse
	
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
			
			d_conf[conf]["RMSD"].append(sqrt((xref-xconf)**2 + (yref-yconf)**2 + (zref-zconf)**2))

	return d_conf
	
def RMSDconf(dico) :
	'''
	Calcul du RMSD de maniere globale a la conformation
	Moyenne du RMSD des residus
	'''
	
	dico["RMSDmoy"] = list()
	dico["RMSDmoy_sd"] = list()
	
	for conf in d_conf["liste_conformations"] :
		dico["RMSDmoy"].append(moyenne(d_conf[conf]["RMSD"]))
		dico["RMSDmoy_sd"].append(ecart_type(d_conf[conf]["RMSD"]))
	
	return dico




def RMSDres(dico1,dico2) :
	'''
	Calcul du RMSD de maniere locale
	Moyenne du RMSD en chaque position a partir de l'ensemble des conformations
	'''
	
	l_res = dico1[dico1["liste_conformations"][0]]["liste_n_residus"]

	dico1["RMSDres_mean"] = list()
	dico1["RMSDres_sd"] = list()

	for i in range(len(l_res)) :
		l_rmsd = list() # liste contenant les valeurs des RMSD d'un residu pour toutes ses conformations
		for conf in dico2["liste_conformations"] :
			l_rmsd.append(dico2[conf]["RMSD"][i])
		
		dico1["RMSDres_mean"].append(moyenne(l_rmsd))
		dico1["RMSDres_sd"].append(ecart_type(l_rmsd))
	return dico1["RMSDres_mean"]

#-----------------------------------------------------------------------
# Fonction de calcul statistique
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
# DISTANCE --> Calcul de la distance de chaque residu au centre de Masse
# de la conformation  --> enfouissement des residus au sein de la conformation

def distance(d_prot) :
	
	global centreMasse
	# pour chaque conformation, son centre de masse :
	xconf = d_prot["liste_CM_x"]
	yconf = d_prot["liste_CM_y"]
	zconf = d_prot["liste_CM_z"]
	
	d_prot["distance_moy"] = list()
	d_prot["distance_sd"] = list()

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
		d_prot["distance_moy"].append(moyenne(l_dist))
		d_prot["distance_sd"].append(ecart_type(l_dist))
	return d_prot



def distanceRes(dico1, dico2) :
	'''
	Calcul de la distance moyenne de chaque residu a centre de masse de la proteine
	'''

	l_res = dico1[dico1["liste_conformations"][0]]["liste_n_residus"]

	dico1["enfRes_mean"] = list()
	dico1["enfRes_sd"] = list()

	for i in range(len(l_res)) :
		l_enf = list()
		for conf in dico2["liste_conformations"] :
			l_enf.append(dico2[conf]["enfouissement"][i])
		
		dico1["enfRes_mean"].append(moyenne(l_enf))
		dico1["enfRes_sd"].append(ecart_type(l_enf))
	return dico1["enfRes_mean"]

#-----------------------------------------------------------------------
# RAYON DE GIRATION --> distance maximale entre un residus et le centre de masse
# d'une conformation <-> residus a l'enfouissement maximal

def max_distance(d_prot) :
	d_prot["rayonGiration"] = list()
	for conf in d_prot["liste_conformations"] :
		d_prot["rayonGiration"].append(max(d_prot[conf]["enfouissement"]))
	return d_prot["rayonGiration"]

def rayonGiration(d_prot1,d_prot2) :
	rayon = max_distance(d_prot2)
	ref = max_distance(d_prot1)[0]
	d_prot2["ratio_giration"] = [x/ref for x in rayon]
	return d_prot2["ratio_giration"]

#-----------------------------------------------------------------------
# REPRESENTATION DES DONNEES

def corEnfouissementFlexibilite(d_prot) :
	'''
	Fonction qui permet de visualiser la correlation entre l'enfouissement
	des residus et la flexibilite des regions
	La flexibilite augmente avec la distance des residus aux CdM
	--> calcul de la correlation
	'''
	d_prot["corEnfFlexi"] = [list(), list()] # correlation et pvaleur
	for conf in d_prot["liste_conformations"] :
		if(d_prot[conf]["RMSD"] == [0] * len(d_prot[conf]["RMSD"])) :
			cor = [1,0]
		else :
			cor = pearsonr(d_prot[conf]["enfouissement"],d_prot[conf]["RMSD"])
		d_prot["corEnfFlexi"][0].append(cor[0])
		d_prot["corEnfFlexi"][1].append(cor[1])
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

#~ def plotRMSD_Temps(listRMSD, listTemps):

#~ def plotGiration_Temps(listGiration, listTemps):

#-----------------------------------------------------------------------
# Ecriture des resultats

def verificationfFichier(output) :
	if (os.path.exists(output)) :
		decision = raw_input("Fichier de sortie "+str(output)+" deja existant : Voulez vous l'ecraser ? O/N\n")
		while ((decision != 'O') & (decision != 'o') & (decision != 'N')  & (decision != 'n')) :
			print "Erreur : Repondre O pour oui ou N pour non" 
			decision = raw_input("Fichier de sortie deja existant : Voulez vous l'ecraser ? O/N\n")
		if ((decision == 'N') | (decision == 'n')) :
			while os.path.exists(output) :
				output = raw_input("Nom de fichier de sortie deja existant : entrez un nouveau nom de fichier : \n")
	return output



def outputGlobaux(output, dico1, dico2,x) :
	try:
		f = open(output, "w")
		
		texte = "Conformation\t|\tRayon Giration\t|\tRMSD\t(+/- ecart-type)\n"
		#~ texte = "Conformation\t|\tRayon Giration\t|\tRMSD\t(+/- ecart-type)\tratio Giration\t ratio RMSD\n"
		
		rayonG_ref = dico1["rayonGiration"][0]
		rmsd_ref = dico1["RMSDmoy"][0]

		texte += "\tREF\t\t\t|\t\t"+str(round(rayonG_ref,x))+"\t\t|\t"+str(round(rmsd_ref,x))+"\n"
		
		for i in range(len(dico2["liste_conformations"])) :
			
			num = dico2["liste_conformations"][i].strip()
			rayonG = dico2["rayonGiration"][i]

			rmsd = dico2["RMSDmoy"][i]

			rmsd_sd = dico2["RMSDmoy_sd"][i]

			#~ ratio_rmsd = rmsd/rmsd_ref # a calculer dans fonction
			#~ ratio_gir = rayonG/rayonG_ref

			texte += "\t"+str(num)+"\t\t\t|\t\t"+str(round(rayonG,x))+"\t|\t"+str(round(rmsd,x))+"\t(+/-"+str(round(rmsd_sd,x))+")\n"
			#~ texte += "\t"+str(num)+"\t\t\t|\t\t"+str(round(rayonG,x))+"\t|\t"+str(round(rmsd,x))+"\t(+/-"+str(round(rmsd_sd,x))+str(round(ratio_gir,x))+str(round(ratio_rmsd,x))+")\n"

		f.write(texte)
		f.close()		
	except:
		print("Erreur chargement fichier global\n")
		sys.exit(0)

def outputLocaux(output, dico1,x) :
	
	try:
		
		f = open(output, "w")

		texte = "Residus\t|\t Nom \t|\tRMSD (+/- ecart-type)\t|\tDistance residu/CdM\t(+/- ecart-type)\n"
		
		conf_ref = dico1["liste_conformations"][0]
		lres = dico1[conf_ref]["liste_n_residus"]

		for i in range(len(lres)) :

			res = lres[i]

			nom = dico1[dico1["liste_conformations"][0]]["liste_seq_residus"][i]
			
			rmsd = dico1["RMSDres_mean"][i]
			rmsd_sd = dico1["RMSDres_sd"][i]
			
			d = dico1["enfRes_mean"][i]
			d_sd = dico1["enfRes_sd"][i]
			#~ texte += "\t"+str(res)+"\t|\t"+nom+"\t|\t"+str(round(rmsd,x))+"\t(+/-"+")\t"+str(round(d,x))+"\t(+/-"+")\n"
			texte += "\t"+str(res)+"\t|\t"+nom+"\t|\t"+str(round(rmsd,x))+"\t(+/-"+str(round(rmsd_sd,x))+")\t"+str(round(d,x))+"\t(+/-"+str(round(d_sd,x))+")\n"
			
		f.write(texte)
		f.close()
		
	except:
		print("Erreur chargement fichier local\n")
		sys.exit(0)

def ecriture(dico1, dico2) :
	output1 = verificationfFichier("res_barstar_globaux.txt") # nom de fichier par defaut, mais on ne veut pas ecraser des resultats precedents
	output2 = verificationfFichier("res_barstar_locaux.txt")
	x = input("Indiquez le nombre de decimales souhaitees :\t")

	outputGlobaux(output1, dico1, dico2,x)
	outputLocaux(output2, dico1,x)
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
	centreMasse = choixMeth(sys.argv[3]) # variable globale ??
	
	# Centre Masse Residus
	d_ref = centreMasseResidus(d_ref)
	d_conf = centreMasseResidus(d_conf)
	
	# Centre Masse Proteine
	d_ref = centreMasseProteine(d_ref)
	d_conf = centreMasseProteine(d_conf)
	
	# RMSD
	d_conf = RMSDresidus(d_ref, d_conf)
	
	# distance
	d_conf = distance(d_conf)
	d_ref = distance(d_ref)


	#~ # AU NIVEAU GLOBAL
	#~ #-----------------
	
	#~ # rayon de Giration
	#~ rayonGiration(d_ref, d_conf)
	#~ plt.plot(d_conf["rayonGiration"])
	#~ plt.show()
	#~ plt.plot(d_conf["ratio_giration"])
	#~ plt.show()
	
	#~ # Distance
	#~ plt.plot(d_conf["distance_moy"])
	#~ plt.show()
	#~ plt.plot(d_conf["distance_sd"])
	#~ plt.show()
	
	# RMSD
	RMSD(d_ref, d_conf)
	plt.plot(d_conf["RMSDmoy"])
	plt.show()
	plt.plot(d_conf["RMSDmoy_sd"])
	plt.show()
	plt.plot(d_conf["ratio_RMSD"])
	plt.show()
	
	#~ # LIENS
	#~ #------
	#~ corEnfouissementFlexibilite(d_conf)
	#~ plt.plot(d_conf["corEnfFlexi"][0])
	#~ plt.show()
	#~ plt.plot(d_conf["corEnfFlexi"][1])
	#~ plt.show()
	
	
	
	
	
	
	#~ plt.plot()
	
	#~ plt.plot(RMSDconf(d_ref))
	#~ plt.show()
	
	#~ plt.plot(RMSDres(d_ref,d_conf))
	#~ plt.show()
	
	#~ distanceRes(d_ref, d_conf)
	#~ ecriture(d_ref, d_conf)

		
	
