

/****************************************************************/
/*								*/
/*          programme servant a corréler 2 séquences            */
/*								*/
/****************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <string.h>
//#include <malloc.h>
//#include <values.h>
#include <sys/time.h>
#include <assert.h>


#include "types.h"
#include "LoadAndSaveImages.h"


/*************************************************************************/
/*                                                                       */
/*                     sauvegarde en floats                              */
/*                                                                       */
/*************************************************************************/
void
save_tab_flt (float *tab_flt, int nb_lignes, int nb_colonnes, char *nom_fic)
{
	FILE *f = NULL;
	int nb_total = nb_lignes * nb_colonnes;

	if ((f = fopen (nom_fic, "wb")) == NULL)
	{
		printf ("Impossible d'ouvrir <%s> en mode wb\n", nom_fic);
		exit (-1);
	}
	fwrite (&nb_lignes, sizeof (int), 1, f);
	fwrite (&nb_colonnes, sizeof (int), 1, f);
	fwrite (tab_flt, sizeof (float), nb_total, f);
	fclose (f);
}

/*************************************************************************/
/*                                                                       */
/*                        lecture de floats                              */
/*                                                                       */
/*************************************************************************/
float *
load_tab_flt (int *nb_lignes, int *nb_colonnes, char *nom_fic)
{
	FILE *f = NULL;
	int nb_total;
	float *tab_tmp = NULL;

	if ((f = fopen (nom_fic, "rb")) == NULL)
	{
		printf ("Impossible d'ouvrir <%s> en mode rb\n", nom_fic);
		exit (-1);
	}
	fread (nb_lignes, sizeof (int), 1, f);
	fread (nb_colonnes, sizeof (int), 1, f);
	nb_total = (*nb_lignes) * (*nb_colonnes);
	if ((tab_tmp = (float *) malloc (sizeof (float) * nb_total)) == NULL)
	{
		printf ("Impossible d'allouer la memoire\n");
		exit (-1);
	}
	fread (tab_tmp, sizeof (float), nb_total, f);
	fclose (f);
	return (tab_tmp);
}

/*************************************************************************/
/*                                                                       */
/*                sauvegarde dans un fichier .txt                        */
/*                                                                       */
/*************************************************************************/
void
enreg_txt (double *donnees, int nb_bins, char nom[256])
{
	FILE *fptr = NULL;
	int a, i;
	double buff;
	char nom_fic[256];

	sprintf (nom_fic, "%s.txt", nom);
	printf ("sauvegarde dans %s\n", nom_fic);

	if ((fptr = fopen (nom_fic, "w")) == NULL)
	{
		printf ("Impossible de creer le fichier %s\n", nom_fic);
		exit (1);
	}

	for (i = 0; i <= nb_bins; i++)
	{
//  if (i==0)
//    {
//      fprintf(fptr,"%d\n",nb_bins);
//    }
//  else
//    {
		buff = donnees[i - 1];
		a = fprintf (fptr, "%f\n", buff);
		if (a < 0)
			printf ("PROBLEME\n");
//    }
	}
	fclose (fptr);
}

/*************************************************************************/
/*                                                                       */
/*                sauvegarde dans un fichier .txt                        */
/*                                                                       */
/*************************************************************************/

double
XCorr (double *x, double *y, int n, int delay)
{
	int i, k;
	double sum = 0, mean1 = 0, mean2 = 0;

	for (i = 0; i < n; i++)
	{
		mean1 += x[i];
		mean2 += y[i];
	}
	mean1 /= n;
	mean2 /= n;

	for (i = 0; i < n; i++)
	{
		k = (i + delay) % n;
		sum += (x[i] - mean1) * (y[k] - mean2);
	}

	return (sum);
}

/*******************************************************************/
/*                                                                 */
/*                      programme main                             */
/*                                                                 */
/*******************************************************************/
int
main (int argc, char *argv[])
{
/*******************************************************/
/*                     variables                       */
/*******************************************************/

	FILE *fptr;
	char nom_X[256];
	char nom_Y[256];
	char TypeCh[256];

	double mx = 0.0, my = 0.0;
	double diff_x = 0.0, diff_y = 0.0;
	double sx, sy, sxy;
	double denom, r;
	double *CrossCorr;
	double *Sigmas;

	int i, j, d, k, MaxD, Type;
	entierns nb_lignes1, nb_colonnes1;
	entierns nb_lignes2, nb_colonnes2;
	int nb_bins1, nb_bins2, nb_bins;
	float max1 = 0, max2 = 0, min1 = 255, min2 = 255;
	float *X, *Y;
	float *X2, *Y2;
	int option;
        int Incomp=0;

	unsigned char *ImageX, *ImageY;

	float *FX, *FY;
	float chrono;
	float a,Force;
	char CHForce[256];


	sprintf (TypeCh, "%s", argv[1]);
	Type = atoi (TypeCh);

	if ((!strcmp (argv[1], "TXT")) || (!strcmp (argv[1], "T")))
		option = 0;
	if ((!strcmp (argv[1], "RAS")) || (!strcmp (argv[1], "R")))
		option = 1;

    if (argv[5] != NULL)
       sprintf (CHForce, "%s", argv[5]);

    Force = atof(CHForce);

    if (argv[5] != NULL)
      printf("\n\n\tForce : %f\n\n",Force);

	chrono = (float) clock ();

	switch (option)
	{
	case 0:

/*******************************************************/
/*            lecture du premier fichier               */
/*******************************************************/
		sprintf (nom_X, "%s.txt", argv[2]);
//  sprintf(nom_X,"/home/fautruss/extraction/extraction/%s.txt",argv[1]);
//  sprintf(nom_X,"%s.txt",argv[1]);
		printf ("ouverture de %s\t\t", nom_X);

		if ((fptr = fopen (nom_X, "r")) == NULL)
		{
			printf ("Impossible d'ouvrir le fichier %s\n", nom_X);
			exit (1);
		}

		fscanf (fptr, "%d\n", &nb_bins1);
		printf ("%d bits\n", nb_bins1);

		X = (float *) calloc (2 * nb_bins1, sizeof (float));

		for (i = 0; i < nb_bins1; i++)
		{
			a = fscanf (fptr, "%f\n", &X[i]);
			if (a < 0)
			{
				printf ("PROBLEME 1\n");
				exit (1);
			}
		}
		fclose (fptr);
/*******************************************************/


/*******************************************************/
/*             lecture du second fichier               */
/*******************************************************/
		sprintf (nom_Y, "%s.txt", argv[3]);
//  sprintf(nom_Y,"/home/fautruss/extraction/extraction/%s.txt",argv[2]);
//  sprintf(nom_Y,"%s.txt",argv[2]);
		printf ("ouverture de %s\t\t", nom_Y);

		if ((fptr = fopen (nom_Y, "r")) == NULL)
		{
			printf ("Impossible d'ouvrir le fichier %s\n", nom_Y);
			exit (1);
		}

		fscanf (fptr, "%d\n", &nb_bins2);
		printf ("%d bits\n", nb_bins2);

		Y = (float *) calloc (2 * nb_bins2, sizeof (float));

		for (i = 0; i < nb_bins2; i++)
		{
			a = fscanf (fptr, "%f\n", &Y[i]);
			if (a < 0)
			{
				printf ("PROBLEME 1\n");
				exit (1);
			}
		}
		fclose (fptr);
/*******************************************************/

        Incomp=0;
		/* vérification compatibilité de la longueur des vecteurs : */
		if (nb_bins1 != nb_bins2)
		{
			printf ("\ndimensions incompatibles\n\n");
            //exit (0);
            if (nb_bins1 > nb_bins2)
            {
               Incomp=1;
               Y2 = (float *) calloc (2 * nb_bins1, sizeof (float));
               for (i=0;i<nb_bins2;i++)
                   Y2[i] = Y[i];
               for (i=nb_bins2;i<nb_bins1;i++)
                   Y2[i] = 0;
               nb_bins = nb_bins1;
            }
            else if (nb_bins2 > nb_bins1)
            {
               Incomp=2;
               X2 = (float *) calloc (2 * nb_bins2, sizeof (float));
               for (i=0;i<nb_bins1;i++)
                   X2[i] = X[i];
               for (i=nb_bins1;i<nb_bins2;i++)
                   X2[i] = 0;
               nb_bins = nb_bins2;
            }
		}
		else if (nb_bins1 == nb_bins2)
			nb_bins = nb_bins1;


		FX = (float *) calloc (2 * nb_bins, sizeof (float));
		FY = (float *) calloc (2 * nb_bins, sizeof (float));

        if (Incomp == 0)
        {
		   for (i = 0; i < nb_bins; i++)
  		   {
			  FX[i] = (float) X[i];
			  FY[i] = (float) Y[i];
		   }
        }
        else if (Incomp == 1)
        {
		   for (i = 0; i < nb_bins; i++)
  		   {
			  FX[i] = (float) X[i];
			  FY[i] = (float) Y2[i];
		   }
        }
        else if (Incomp == 2)
        {
		   for (i = 0; i < nb_bins; i++)
  		   {
			  FX[i] = (float) X2[i];
			  FY[i] = (float) Y[i];
		   }
        }

        if (argv[5] != NULL)
          for (i = 0; i < nb_bins; i++)
            FY[i] *= Force;


		for (i = 0; i < nb_bins; i++)
		{
			if (FX[i] >= max1)
				max1 = FX[i];
			if (FX[i] <= min1)
				min1 = FX[i];

			if (FY[i] >= max2)
				max2 = FY[i];
			if (FY[i] <= min2)
				min2 = FY[i];
		}

		printf ("\npour %s : max = %f\tmin = %f\n", argv[2], max1,
			min1);
		printf ("pour %s : max = %f\tmin = %f\n\n", argv[3], max2,
			min2);



		break;

	case 1:

		sprintf (nom_X, "%s", argv[2]);
		sprintf (nom_Y, "%s", argv[3]);

		lit_raster (&ImageX, nom_X, &nb_lignes1, &nb_colonnes1);
		lit_raster (&ImageY, nom_Y, &nb_lignes2, &nb_colonnes2);

		if ((nb_lignes1 != nb_lignes2)
		    || (nb_colonnes1 != nb_colonnes2))
		{
			printf ("\n\n\t\tdimensions des images incompatibles...\n\n");
			exit (0);
		}

		nb_bins = nb_lignes1 * nb_colonnes1;

		FX = (float *) calloc (2 * nb_bins, sizeof (float));
		FY = (float *) calloc (2 * nb_bins, sizeof (float));

		for (i = 0; i < nb_bins; i++)
		{
			FX[i] = (float) ImageX[i];
			FY[i] = (float) ImageY[i];
		}

		break;
	}			// fin du switch.


/****************************************/
	printf ("\n\tCross-Correlation : \n\n");
/****************************************/

	// INITIALISATIONS :
	mx = 0.0;
	my = 0.0;
	diff_x = 0.0;
	diff_y = 0.0;

	//  calcul des Moyennes :
	for (i = 0; i < nb_bins; i++)
		mx += (double) FX[i];
	mx = mx / (nb_bins);

	for (i = 0; i < nb_bins; i++)
		my += (double) FY[i];
	my = my / (nb_bins);

	printf ("moyX = %f\t", mx);
	printf ("moyY = %f\n\n", my);

	CrossCorr = (double *) calloc (3 * nb_bins, sizeof (double));
	//Sigmas = (double *) calloc (3*nb_bins,sizeof(double));


	sx = 0;
	sy = 0;
	for (i = 0; i < nb_bins; i++)
	{
		sx += (FX[i] - mx) * (FX[i] - mx);
		sy += (FY[i] - my) * (FY[i] - my);
	}
	denom = sqrt (sx * sy);
	printf ("denom = %f\n\n", denom);


	k = 0;
	MaxD = (int) nb_bins;
	for (d = -MaxD; d < MaxD; d++)
	{
		sxy = 0;
		for (i = 0; i < nb_bins; i++)
		{
			j = i + d;
			if (j < 0 || j >= nb_bins)
				continue;
			sxy += (FX[i] - mx) * (FY[j] - my);
		}
		r = sxy / denom;
		CrossCorr[k] = r;
		k++;
	}

	enreg_txt (CrossCorr, 2 * nb_bins, argv[4]);


	printf ("\nTemps : %6.3f sec.\n",
		(clock () - chrono) / CLOCKS_PER_SEC);
	printf ("Temps : %6.3f min.\n",
		(clock () - chrono) / (CLOCKS_PER_SEC * 60));

	return (0);

}
