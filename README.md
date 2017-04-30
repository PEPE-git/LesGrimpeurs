 LesGrimpeurs
 ------------
Projet BARstar de DA SILVA OPHELIE et MERCKAERT PIERRE

 Synopsis
 --------

Une étude de dynamique moléculaire du complexe barstar-barnase a montré l'importance de deux résidus clés appartenant à l'interface de la barstar dans l'interaction avec la barnase (Kimura et al, Biophysical Journal, 2001). Ces deux résidus subissent un changement conformationnel lors de l'interaction de la barstar avec la barnase et jouent un rôle clé dans cette interaction. Dans ce projet, nous avons réalisé une dynamique moléculaire de la barstar afin d'étudier en fonction des conformations la stabilité de la structure et les changements conformationnels subis.
En particulier, on s'intéresse aux changements conformationnels globaux (compaction, dépliement,mouvements de boucles de grande amplitude, mouvements allostériques...) et locaux (petits mouvements de boucles, mouvements de chaîne latérales...).


 Code Example
 ------------

Utilisez Barstar.py pour obtenir toutes les informations pertinentes quant aux changements conformationnels de votre protéine d'intérêt par rapport à sa référence.
code :: python
	#Parsing des fichiers .pdb en arguments sous forme de dictionnaire
	barstar_conf= __parsePDBM(barstar_conf.pdb)
	barstar_ref = __parsePDBM(barstar_ref.pdb)

	#Calcul des changements conformationnels vis-à-vis de la référence et ajout de ces valeurs dans le dictionnaire avec des clés pertinentes.
	l_dict=[barstar_conf, barstar_ref]

	#Calcul du centre de masse des résidus et de chaque conformation
	centreMasseCalc(l_dict)	
	#Calcul du RMSD des residus de chaque conformation par rapport à la référence			
	RMSD(l_dict)			
	#Calcul de la distance entre les résidus et le CdM pour chaque conformation			
	distance(l_dict)	
	#Calcul du Rayon de Giration pour chaque conformation				
	rayonGiration(l_dict)				
	#Calcul de la corrélation entre le Rayon de Giration (Enfouissement) et la Flexibilité (RMSD) pour chaque conformation
	corEnfouissementFlexibilite_conf(l_dict)		


 API
 ---
??????????????
Depending on the size of the project, if it is small and simple enough the reference docs can be added to the README. For medium size to larger projects it is important to at least provide a link to where the API reference docs live.
??????????????

 License
 -------
??????????????
A short snippet describing the license (MIT, Apache, etc.)
??????????????