////////////////////////////////////////////////////////////////////////////////////
//		       	types
////////////////////////////////////////////////////////////////////////////////////
typedef long 	         		entier;     // entier doit etre sur 32 bits (int ou long)
typedef unsigned  long 	 		entierns;	// entierns doit etre sur 32 bits (uint ou ulong)
typedef float	         		reel;	    // reel  (float ou double)
typedef unsigned  short       	motns;		// mot 16 bits non signé
typedef unsigned  char        	octetns;	// octet 8 bits non signé

////////////////////////////////////////////////////////////////////////////////////
//    		        define
////////////////////////////////////////////////////////////////////////////////////
#define	       SAVEDCP          	1      			// 1:sauver toutes les infos de la DCP
#define	       SAVEALL	       		0      			// 1:sauver toutes les infos (sauf infos DCP,image du parcours & image de contraste)
#define	       DEBUGINFO        	0      			// affiche des infos de debug
#define	       DEBUGMODE       		0      			// génère des données de debug (tableaux, fichiers...)
#define        BMPMAGIC				19778  			// "BM"
#define	       PIx2	     			6.28319   		// 2*pi
#define	       PIon4	     		0.785398  		// pi/4
//#define	    PI	     						3.14159   	// pi

////////////////////////////////////////////////////////////////
//    		       LUM defines
////////////////////////////////////////////////////////////////
#define         luminancemax    	90     	// luminance max de l'écran
#define         luminanceRmax   	30     	// luminance max de l'écran pour la composante R
#define         luminanceGmax   	30     	// luminance max de l'écran pour la composante G
#define         luminanceBmax   	30     	// luminance max de l'écran pour la composante B
#define         gammaR          	1.0    	// gamma pour R
#define         gammaG          	1.0    	// gamma pour G
#define         gammaB          	1.0    	// gamma pour B
#define         Rmax	        	255    	// 256 niveaux
#define         Gmax	       		255    	// 256 niveaux
#define         Bmax	       		255    	// 256 niveaux


////////////////////////////////////////////////////////////////
//    		       FFT defines
////////////////////////////////////////////////////////////////
#define	       OUVERT	       		1      				// from pgr_fft.h  (definitions FFT)
#define	       FERME	       			0      				// from pgr_fft.h  (definitions FFT)
#define	       W1_3r           		-0.5    				// from pgr_fft.h  (definitions FFT)
#define	       W1_3i	  				0.86602540   	// from pgr_fft.h  (definitions FFT)
#define	       W1_5r	  				0.30901699   	// from pgr_fft.h  (definitions FFT)
#define	       W1_5i	  				0.95105652   	// from pgr_fft.h  (definitions FFT)
#define	       W2_5r	 				-0.80901699  	// from pgr_fft.h  (definitions FFT)
#define	       W2_5i	  				0.58778525   	// from pgr_fft.h  (definitions FFT)
#define	       DIRECTE	       		1.0    				// from constant.h (transformee de Fourier
#define	       INVERSE	       		-1.0    				// directe ou inverse)
#define	       irint(x)	   				(((x)>0) ? ((entier)((x)+0.5)) : ((entier)((x)-0.5))) // from ???
#define	       SQUARELF(x) 		((reel)(x) * (reel)(x))			// from macro.h
// #define	   M_PI						3.14159265359  // from Borland C++ 4.5: constante M_PI


////////////////////////////////////////////////////////////////
//    		       RAS defines
////////////////////////////////////////////////////////////////
#define  RAS_MAGIC              			0x59a66a95
#define  RAS_MAGIC_SWAPPED      	0x956aa659
#define  RT_MODERN      					1
#define  RMT_NONE       						0
#define  RMT_EQUAL_RGB					1



////////////////////////////////////////////////////////////////
//    		   GLOBAL defines
////////////////////////////////////////////////////////////////

#define 			BMPorRAS      				0  			// == 1 si on veut traiter du BMP,
															// == 0 si on veut traiter du RAS.

////////////////////////////////////////////////////////////////
//    	SWAP d'octets pour
//		compatibilité Windows / Unix.
////////////////////////////////////////////////////////////////
void swapLong (char *p) {
    char t;
    t = p[0];
    p[0] = p[3];
    p[3] = t;
    t = p[1];
    p[1] = p[2];
    p[2] = t; }

////////////////////////////////////////////////////////////////////////////////////
//	  		struct
////////////////////////////////////////////////////////////////////////////////////

typedef struct
{

  int             ras_magic;   /* = RAS_MAGIC */
  int             ras_width;
  int             ras_height;
  int             ras_depth;  /* souvent 8bpp */
  int             ras_length; /* nb bytes of image datas */
  int             ras_type; /* RT_OLD or RT_EXPERIMENTAL */
  int             ras_maptype;  /* RMT_NONE */
  int             ras_maplength;  /* 0 */

}	raster_header;

typedef struct
{
  motns	       	Signature;      // = "BM"
  entierns     	FileSize;       // taille du fichier en octets
  entierns     	Reserved;       // non utilisé (=0)
  entierns     	DataOffset;    	// offset vers les données RASTER

}	BitmapFileHeaderType;

typedef struct
{
  entierns     	Size;	        // Taille de l'InfoHeader(=40)
  entierns     	Width;	        // Largeur du Bitmap
  entierns     	Height;	        // Hauteur du Bitmap
  motns	       	Planes;  	// Nombre de plans (=1)
  motns	       	BitCount;       // Nombre de bits par pixel : (1:NumColors=1, 4:NumColors=16,
                                        //8:NumColors=256, 16:NumColors=65536, 24:NumColors=16M)

  entierns     	Compression;   	// Type de compression (0:pas de compression, 1:8bits RLE, 2:4bits RLE)
  entierns     	ImageSize;     	// Taille compressée de l'image (=0 si Compression=0)
  //entierns   	XpixelsPerM;  	// 0
  entierns     	XpixelsPerM;   	// résolution horizontale (pixels/metre)
  entierns     	YpixelsPerM;   	// résolution verticale (pixels/metre)
  entierns     	ColorsUsed;    	// nombre de couleurs utilisées
  entierns     	ColorsImportant;// nombre de couleurs importants (0: toutes)
}	BitmapInfoHeaderType;

typedef struct
{
  octetns       Blue;	       	// bleu
  octetns       Green;	       	// vert
  octetns		Red;	       	// rouge
}	RGBTri; // (NumColors) fois (puis BitmapInfoHeaderType.ImageSize octets)

////////////////////////////////////////////////////////////////////////////////////
//									vars
////////////////////////////////////////////////////////////////////////////////////
