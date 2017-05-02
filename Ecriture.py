#!/usr/bin/env python2

"""
Author: LesGrimpeurs
Date: 02/05/2017
Description: Projet Barstar
"""

import sys, os
import csv

import Conformation_analysis as conf_anal

'''
Fonction Ecriture du projet BARSTAR
Creation du dossier (nom fonction du nombre de conformations analysees et de la methode pour le centre de masse chosie)
dans lequel seront ecrits et enregistres les resultats (tableaux excel des valeurs des parametres calcules et images des graphes)

'''


#-----------------------------------------------------------------------
# ECRITURE DES RESULTATS
#-----------------------------------------------------------------------


def ecriture(l_dict,methode) :

	'''
	Creation du dossier de sortie s'il n'existe pas.
	Verification que les fichiers de sortie n'existent pas deja dans le dossier.
	
	'''

	if not os.path.exists('Barstar_Results_'+methode+'/'):
		os.makedirs('Barstar_Results_'+methode)
		
	output1 = __verificationfFichier("Barstar_Results_"+methode+"/res_barstar_globaux_"+methode+".csv") # nom de fichier par defaut, mais on ne veut pas ecraser des resultats precedents
	output2 = __verificationfFichier("Barstar_Results_"+methode+"/res_barstar_locaux_"+methode+".csv")
	print "Ecriture dans le dossier Barstar_Results_"+methode+"/ les fichiers de sortie 'res_barstar_globaux_"+methode+".csv' et 'res_barstar_locaux_"+methode+".csv'"
	x = raw_input("\nIndiquez le nombre de decimales souhaitees :\t")
	while not x.isdigit() :
		print "Indiquez un chiffre svp"
		x = raw_input("Indiquez le nombre de decimales souhaitees :\t")

	__outputGlobaux(output1, l_dict[0], l_dict[1], int(x))
	__outputLocaux(output2, l_dict[0], int(x))


def __verificationfFichier(output) :
	#Si le fichier existe deja, demande si l'utilisateur veut le remplacer ou non
	#si non, entrer un nouveau non de fichier.

	if (os.path.exists(output)) :
		decision = raw_input("\nFichier de sortie "+str(output)+" deja existant : Voulez vous l'ecraser ? O/N\n")
	
		while ((decision != 'O') & (decision != 'o') & (decision != 'N')  & (decision != 'n')) :
			print "Erreur : Repondre O pour oui ou N pour non" 
			decision = raw_input("Fichier de sortie deja existant : Voulez vous l'ecraser ? O/N\n")
		
		if ((decision == 'N') | (decision == 'n')) :
			while os.path.exists(output) :
				output = raw_input("Nom de fichier de sortie deja existant : entrez un nouveau nom de fichier : \n")
		elif ((decision == 'O') | (decision == 'o')) :
			rm = 'rm '+output
			os.popen(rm)	
	return output



def __outputGlobaux(output, d_ref, d_conf, x) :
	#Ecriture des valeurs des parametres de l'analyse locale dans un tableau excel

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
	#Ecriture des valeurs des parametres de l'analyse globale dans un tableau excel

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
