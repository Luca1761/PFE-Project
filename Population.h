/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <time.h>
#include "Node.h"
#include "Individual.h"

using namespace std ;

struct SousPop
{
	// individus de la sous-population
	vector<Individual*> individus ;

	// nombre de ces individus
	unsigned int nbIndiv ;
};

class Population
{
   private:

   // Acces aux parametres de l'instance et du genetique
   Params* params;

   unsigned int nbScenario;

   // liste qui repertorie si les XXX derniers individus produits etaient valides en terme de charge ou non
   list<bool> listeValiditeCharge ;

   // auxiliary data structure to apply the local search
   Individual * trainer;

   void education_scenario(Individual * indiv);

   // fonction booleenne verifiant si le fitness n'existe pas deja
   bool fitExist(SousPop * pop, Individual * indiv);

   // place un individu donne dans le tableau
   // retourne une erreur si echec sinon sa position dans la population
   unsigned int placeIndividu(SousPop * pop, Individual * indiv);

   public:

   // clock time when the best individual was found
   clock_t timeBest ;

   clock_t totalTime ;

   // calcule le fitness etendu des elements de la sous-population
   void evalExtFit(SousPop * pop);

   // ajoute un individu a la population
   // la population se debrouille pour le placer ou il lui semble bon
   // updateNbValides est mis a true si on veut mettre jour la proportion de valides aussi 
   // retourne une erreur si echec, sinon sa nouvelle position dans la population
   unsigned int addIndividu (Individual * indiv);

   // enleve un individu de la population en fonction de la replacement policy
   void removeIndividu(SousPop * pop, unsigned int p);

   // met a jour les individus les plus proches d'une population
   // en fonction de l'arrivant
   void updateProximity (SousPop * pop, Individual * indiv);

   // procede de redemarrage avec remplacement d'une partie de la population
   // modele tres simplifie
   // on remplace la moitie individus de fitness situes sous la moyenne par des individus aleatoires
   void diversify ();

   // recopie un Individu dans un autre
   // ATTENTION !!! ne recopie que le chromT et les attributs du fitness
   void recopieIndividu (Individual * destination , Individual * source);
   
   // differents individus valides presents dans la population
   SousPop * valides;

   // differents individus invalides presents dans la population
   SousPop * invalides;

   // getter de 1 individu par binary tournament
   Individual * getIndividuBinT (double & rangRelatif);

   // getter du meilleur individu valide
   // retourne NULL si il n'y a pas de valides
   Individual * getIndividuBestValide ();

   // getter du meilleur individu invalide
   // retourne NULL si il n'y a pas de invalides
   Individual * getIndividuBestInvalide ();

   // compromis entre fitness et diversite
   unsigned int selectCompromis (SousPop * souspop) ;

   // recalcule l'evaluation des individus a partir des violations
   // puis effectue un tri a bulles de la population
   // operateur de tri peu efficace mais fonction appellee tres rarement
   void validatePen ();

   // getter simple d'un individu
   Individual * getIndividu(unsigned int p);

   // exporte la meilleure solution
   void ExportPop(string nomFichier, bool add, std::vector<double> deliveries, double &totalCost) ;

   // retourne la fraction d'individus valides en terme de charge
   double fractionValidesCharge() ;

   // diversite de la population
   double getDiversity(SousPop * pop);

   // retourne la moyenne des solutions valides
   double getMoyenneValides();

    // retourne la moyenne des solutions invalides
   double getMoyenneInvalides();

   // affiche l'etat de la population
   void afficheEtat(unsigned int NbIter);

   // met a jour le compte des valides
   void updateNbValides (Individual * indiv);

   void measureAndUpdateQuantities(std::vector<double> &deliveries, double &totalCost);

   //constructeur
   Population(Params* params);

   //destructeur
   ~Population();
};

#endif
