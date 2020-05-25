# -*- coding: utf-8 -*-
#Marcheix Benjamin & Bourdeau Emie & Sebile Wanda & Broyez Vianney

class Smith_Waterman(object):
    #############################################################
    #Constructeur
    #############################################################
    #score - la table des scores
    #seuil_match_minimum - proportion de similarité minimum entre les séquences
    #seq_ref - séquence de reference
    #print_details - si l'utilisateur souhaite afficher les détails de l'alignement
    def __init__(self,seq_ref,pos_etudiees={},score={},seuil_match_minimum=0.95,print_details=False, get_align_seq=True):
        self.__seuil_match_minimum=seuil_match_minimum#Proportion de similarité minimum entre les séquences
        self.__seq_ref=seq_ref#La séquence de référence
        self.print_details=print_details#Si l'utilisateur souhaite afficher les détails des résultats de chaque étape
        self.__first_seq=True#Pour éviter les duplications dans le dictionnaire des positions étudiées lors des appels récursifs
        self.__get_align_seq=get_align_seq

        #Le dictionnaire des postiions étudiées est initialisé
        self.__pos_etudiees={}
        for pos in pos_etudiees:
            self.__pos_etudiees[pos]={"A":0,"C":0,"G":0,"T":0,"D":{"A":0,"C":0,"G":0,"T":0},"I":{"A":0,"C":0,"G":0,"T":0},"total":0}

        #Matrice des scores
        if(score=={}):#Matrice par défaut
            self.score={}
            self.score["A"]={"A":3.0,"C":-1.0,"G":-1.0,"T":-1.0}
            self.score["C"]={"A":-1.0,"C":3.0,"G":-1.0,"T":-1.0}
            self.score["G"]={"A":-1.0,"C":-1.0,"G":3.0,"T":-1.0}
            self.score["T"]={"A":-1.0,"C":-1.0,"G":-1.0,"T":3.0}
            self.score["-"]=-2.0
        else:#Si __matrice spécifiée
            self.score=score
        #Moyenne des scores pour un match de nucléotides
        self.moy_match_score=(self.score["A"]["A"]+self.score["T"]["T"]+self.score["G"]["G"]+self.score["C"]["C"])/4.0
        self.__reset_align()#Déclare les attributs qui seront réinitialisés après chaque alignement de read

    #############################################################
    #Initialise une nouvelle matrice
    #############################################################
    #read - le read à aligner
    def __initialisation_matrice(self,read):
        larg = len(read)+1
        haut = len(self.__seq_ref)+1
        for lig in range(haut): # Pour chaques lignes de la matrice
            ligne = [] # Création d'une liste pour chaques lignes de la matrice 
            for col in range(larg): # Parcours des éléments dans chaques lignes
                ligne.append(0.0) # Ajout de 0 dans chaques lignes
            self.__matrice.append(ligne) # Ajout de chaques lignes dans la matrice

    #############################################################
    #Actualise la valeur max de la matrice en fonction du point
    #spécifié
    #############################################################
    #i, j- les coordonnées du point à analyser
    def __update_max(self,i,j):
        if(self.__matrice[i][j]>self.__max): #Si un nouveau max est trouvé
            self.__max=self.__matrice[i][j]#Actualise la valeur du max
            self.__liste_pos_max=[]#La liste est réinitialisée
        if(self.__matrice[i][j]==self.__max):#Si une valeur équivalente au max est trouvée
            self.__liste_pos_max.append([i,j])#Ajoute ses coordonnées dans la liste

    #############################################################
    #Crée une nouvelle matrice
    #############################################################
    #read - les deux séquences à aligner
    #score - la matrice des scores
    def __creation_matrice(self, read):
        self.__initialisation_matrice(read)
        for i in range (1,len(self.__matrice)): # Parcours des lignes
            for j in range (1,len(self.__matrice[0])): # Parcours dans les lignes
                self.__matrice[i][j] = max(self.score[read[j-1]][self.__seq_ref[i-1]]+self.__matrice[i-1][j-1],self.score["-"]+self.__matrice[i][j-1],self.score["-"]+self.__matrice[i-1][j],0.0)
                self.__update_max(i,j)#Modifie le max pour qu'il corresponde à la valeur maximum trouvée
        if(self.print_details):
            print("matrice créee : "+str(len(self.__liste_pos_max))+" max trouvés pour une valeur de "+str(self.__max))

    #############################################################
    #Affiche la __matrice
    #############################################################
    #__matrice - la __matrice à afficher
    def print_matrice(self, maxNumbreSize=6):
        print("\n", end="")
        for i in self.__matrice:
            for n in i:
                nbrSpaces=maxNumbreSize-len(str(n))
                print(str(n)+(nbrSpaces*" "), end="")
            print("")
        print("\n", end="")
    #############################################################
    #Trouve les meilleurs alignements entre les deux séquences s1
    #et s2. Retourne une liste avec les résultats.
    #############################################################
    #lig, col - les coordonnées d'où lancer l'alignement
    #s1, s2 - les deux séquences à aligner
    #align1, align2 - les deux séquences de l'alignement en cours
    #score - la matrice des scores
    #tmp_liste_alignements - la liste des alignements
    def __aligner(self,lig,col,align1,align2,s1,s2,tmp_liste_alignements) :#affiche l'ensemble des alignements
        if(lig in self.__pos_etudiees and len(align1)!=0):#Si la position est contenue dans le dictionnaire des positions à étudier
            if(align2[0]=="-"):#En cas d'insertion
                self.__pos_etudiees[lig]["I"][align1[0]]+=1#Ajoute la base insérée dans le dictionnaire
            elif(align1[0]=="-"):#En cas de déletion
                self.__pos_etudiees[lig]["D"][align2[0]]+=1#Ajoute la base délétée dans le dictionnaire
                if(self.__first_seq):
                    self.__pos_etudiees[lig]["total"]+=1
            elif(self.__first_seq):#Cas normal
                self.__pos_etudiees[lig][align1[0]]+=1#Ajoute le résultat de l'alignement au dictionnaire
                self.__pos_etudiees[lig]["total"]+=1

        if((lig==0 and col==0) or self.__matrice[lig][col]==0): # On enregistre dans la liste l'alignement car on est arrivé au bout
            tmp_liste_alignements.append([align1,align2])
            self.__first_seq=False
        else:
            #test de la case au dessus :
            if lig>0 :#si il y a une case au dessus
                if self.__matrice[lig-1][col]+self.score["-"]==self.__matrice[lig][col] :#si on peut venir de cette case 
                    self.__aligner(lig-1,col,'-'+align1,s2[lig-1]+align2,s1,s2,tmp_liste_alignements)
            #test de la case de gauche :
            if col>0 :#si il y a une case a gauche :
                if self.__matrice[lig][col-1]+self.score["-"]==self.__matrice[lig][col] :#si on peut venir de cette case : 
                    self.__aligner(lig,col-1,s1[col-1]+align1,'-'+align2,s1,s2,tmp_liste_alignements)
            #test de la diagonale :
            if lig>0 and col>0 :#si il existe une case en diagonale :
                if self.__matrice[lig-1][col-1]+self.score[s1[col-1]][s2[lig-1]]==self.__matrice[lig][col] :#si on peut venir de la diagonale :
                    self.__aligner(lig-1,col-1,s1[col-1]+align1,s2[lig-1]+align2,s1,s2,tmp_liste_alignements)
        return(tmp_liste_alignements)

    #############################################################
    #Réinitialise les attributs du dernier alignement
    #############################################################
    def __reset_align(self):
        self.__matrice=[]#La matrice de Smith et Waterman
        self.__liste_alignements=[]#Liste de tous les alignements
        self.__liste_pos_max=[]#Liste de toutes les positions correspondant au maximum
        self.__max=-1#La valeur maximum dans la matrice

    #############################################################
    #Aligne le read passé en paramètre. Retourne le tableau des 
    #positions etudiees
    #############################################################
    #read - le read à aligner
    def align(self, read):
        self.__reset_align()
        self.__creation_matrice(read)#Crée la matrice
        maxScoreAlignement=self.moy_match_score*float(len(self.__matrice[0]))-3.0#Le score théorique d'un alignement parfait
        seuilMatriceMinimum=maxScoreAlignement-(1.0-self.__seuil_match_minimum)*maxScoreAlignement#Le seuil minimum de non rejet du max
        if(self.__max>seuilMatriceMinimum and self.__get_align_seq):#Si le max a une valeur acceptable
            for i in self.__liste_pos_max:#Lance l'alignement en partant de chaque maximum
                tmp_liste_alignements=[]
                self.__liste_alignements.append(self.__aligner(i[0],i[1],"","",read,self.__seq_ref,tmp_liste_alignements))
                self.__first_seq=True
        #Si le mode affichage est activé, affiche les détails
        if(self.print_details):
            print("maximum alignment score : "+str(int(maxScoreAlignement)))
            print("seuil : "+str(int(seuilMatriceMinimum)))
            if(self.__max<=seuilMatriceMinimum):
                print("read exclu car inférieur au seuil")
            else:
                print("meilleurs alignements trouvés : ")
                print(self.__liste_alignements)
        return(self.__pos_etudiees)

    #############################################################
    #Aligne toutes les séquences de la liste et retourne le 
    #tableau des positions étudiées
    #############################################################
    #read - le read à aligner
    def multiple_align(self, reads):
        from tqdm import tqdm
        for i in tqdm(reads,unit="read",bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}|remaining:{remaining}|{rate_fmt}"):
            self.align(i)
        return(self.__pos_etudiees)

    #############################################################
    #Getter pour le dictionnaire 
    #############################################################
    def get_dict_sw(self):
        return(self.__pos_etudiees)

    #############################################################
    #Pour retourner la liste des alignements dans un print
    #############################################################
    def __repr__(self):
        tmpStr=""
        for pos in self.__pos_etudiees:
            tmpStr+=str(pos)+" -> "+str(self.__pos_etudiees[pos])+"\n"
        return(tmpStr)