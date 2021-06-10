/////////////////////////////////////////////////////////////////////////////
// myAssert
// - signale une erreur a l'aide d'un message si !condition
/////////////////////////////////////////////////////////////////////////////
void myAssert(entier condition, char *message)
{
  if (!condition)
    {
      printf("myAssert: %s \n", message);
      exit(-1);
    }
}

///////////////////////////////////////////////////////////////////////////////
// concat
// concatenation de 2 chaines
///////////////////////////////////////////////////////////////////////////////
char* concat(char *a, char *b)

{
  char 	        *dummy;	       	// utilise pour le resultat
  entierns      size=0;

  size = (a !=NULL) ? strlen(a) : 0 ;
  size+= (b !=NULL) ? strlen(b) : 0 ;
  dummy=(char*)malloc(size+1);	       	// Allocation de la chaine resultat
  myAssert( (dummy!=NULL), "concat: Out of Memory!");
  strcpy(dummy,"");
  if (a !=NULL) strcat(dummy,a);       	// Concatenation
  if (b !=NULL) strcat(dummy,b);        // Concatenation
  return dummy;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// VerticalFlip
// Flip vertical de l'image pour compenser le format BMP
/////////////////////////////////////////////////////////////////////////////////////////////
void VerticalFlip(unsigned char *image, entierns SizeX, entierns SizeY, char mode)
{

  char		dummy;
  entierns	iX,iY;
  entierns	datatowrite;

  // printf("VerticalFLip: BEFORE: R(0,0) = %d\n",image[0]);
  // printf("VerticalFLip: BEFORE: R(H,0) = %d\n",image[(SizeY-1)*SizeX]);
  datatowrite=SizeX*SizeY;
  for (iY=0; iY<(SizeY/2); iY++)
    {
      for (iX=0; iX<SizeX; iX++)	
	{
	  dummy=image[iY*SizeX+iX];
	  image[iY*SizeX+iX]=image[(SizeY-1-iY)*SizeX+iX];
	  image[(SizeY-1-iY)*SizeX+iX]=dummy;
	  if (mode==3)
	    {
	      dummy=image[iY*SizeX+iX+datatowrite];
	      image[iY*SizeX+iX+datatowrite]=image[(SizeY-1-iY)*SizeX+iX+datatowrite];
	      image[(SizeY-1-iY)*SizeX+iX+datatowrite]=dummy;
	      dummy=image[iY*SizeX+iX+2*datatowrite];
	      image[iY*SizeX+iX+2*datatowrite]=image[(SizeY-1-iY)*SizeX+iX+2*datatowrite];
	      image[(SizeY-1-iY)*SizeX+iX+2*datatowrite]=dummy;
	    }
	}
    }
  // printf("VerticalFLip: AFTER: R(0,0) = %d\n",image[0]);
  // printf("VerticalFLip: AFTER: R(H,0) = %d\n",image[(SizeY-1)*SizeX]);
}

/************************************************************/
/*                                                          */
/*             lecture d'une image raster                   */
/*                                                          */
/************************************************************/
void lit_raster(unsigned char **image,char *chemin,entierns *lig,entierns *col)
{
    FILE             *textu;
    raster_header     entete;
    int               i;
    unsigned char     rasterLut[768];
    int               colorNumber;
    char cheminras[256];


    sprintf(cheminras,"%s.ras",chemin);


    if ((textu=fopen(cheminras,"r")) == NULL)
      printf("%s inexistant\n",cheminras);

    fread(&entete,sizeof(entete),1,textu);

    if (entete.ras_magic==RAS_MAGIC_SWAPPED)
      {
	swapLong ((char *) &entete.ras_magic);
	swapLong ((char *) &entete.ras_height);
	swapLong ((char *) &entete.ras_width);
	swapLong ((char *) &entete.ras_depth);
	swapLong ((char *) &entete.ras_length);
	swapLong ((char *) &entete.ras_type);
	swapLong ((char *) &entete.ras_maptype);
	swapLong ((char *) &entete.ras_maplength);
      }

    if (entete.ras_magic!=RAS_MAGIC)
      {
	printf("%s n'est pas du raster!\n",cheminras);
	exit(1);
      }

    *col = entete.ras_width;
    *lig = entete.ras_height;

    printf ("\npour le fichier %s :\n\tnb_lignes=%ld  nb_colonnes=%ld\n",cheminras,*lig,*col);

    if (entete.ras_depth != 8)
	{
	    printf("%s n'est pas du 8bpp!\n",cheminras);
	    exit(1);
	}

    if (entete.ras_maplength > 0)
      {
	assert (entete.ras_maplength <= sizeof (rasterLut));
	if (fread (rasterLut,1,entete.ras_maplength,textu) != entete.ras_maplength)
	  printf("pb pour lecture de lut!\n");
      }

    *image = (unsigned char *)malloc((*lig)*(*col));
    fread((*image),1,(*lig)*(*col),textu);

    if (entete.ras_maplength > 0)
      {
	colorNumber = entete.ras_maplength/3;
	for (i = 0; i < colorNumber; i++)
	  rasterLut [i] = (rasterLut [i] + rasterLut [i + colorNumber] + rasterLut [i + colorNumber*2])/3;

	for (i = 0; i < *lig * *col; i++)
	  (*image) [i] = rasterLut [(*image) [i]];
      }

    fclose(textu);
}

/************************************************************/
/*                                                          */
/*           sauvegarde d'une image raster                  */
/*                                                          */
/************************************************************/
void save_raster(unsigned char *image,char *chemin,int Lignes,int Colonnes)
{
    FILE         *textu;
    raster_header entete;

    entete.ras_magic=RAS_MAGIC;
    entete.ras_width=Colonnes;
    entete.ras_height=Lignes;
    entete.ras_depth=8;
    entete.ras_length=Lignes*Colonnes;
    entete.ras_type=RT_MODERN;
    entete.ras_maplength= 0;
    entete.ras_maptype=RMT_NONE;

#ifdef LITTLE_ENDIAN
    swapLong ((char *) &entete.ras_magic);
    swapLong ((char *) &entete.ras_height);
    swapLong ((char *) &entete.ras_width);
    swapLong ((char *) &entete.ras_depth);
    swapLong ((char *) &entete.ras_length);
    swapLong ((char *) &entete.ras_type);
    swapLong ((char *) &entete.ras_maptype);
    swapLong ((char *) &entete.ras_maplength);
#endif

    fflush(stdout);

    sprintf(chemin,"%s.ras",chemin);
//    printf("\n   sauvegarde dans %s\n",chemin);

    textu=fopen(chemin,"w");
    fwrite(&entete, sizeof(raster_header), 1,               textu);
    fwrite(image,   1,                     Lignes*Colonnes, textu);

    fclose(textu);
}


/////////////////////////////////////////////////////////////////////////////////////////////
// LoadBMP
// Charge une image en BMP 24 bits
/////////////////////////////////////////////////////////////////////////////////////////////
entierns LoadBMP(char *chemin, char *thefilename, unsigned char **image, entierns *SizeX, entierns *SizeY)
{
  FILE 					*bmpfile;

  BitmapFileHeaderType 	BitmapFileHeader;	// file header
  BitmapInfoHeaderType 	BitmapInfoHeader;	// info header
  entierns				datatoread; 		// quantite de donnees a lire
  entierns				readdata;   		// quantite de donnees lue
  entierns				i,iX,iY,iPadding;
  RGBTri					RGBt;
  entierns				PaddingNum;
  octetns					Padding;
  char					BMPCase;

  fpos_t					GrosBugPos;



  bmpfile=fopen(concat(concat(chemin,thefilename),".bmp"),"rb");
  if (bmpfile==NULL)
    {
      printf("LoadBMP: Can't open file %s.bmp\n",thefilename) ;
      exit(-1) ;
    }

  // lecture du file header
  fread(&BitmapFileHeader.Signature, sizeof(BitmapFileHeader.Signature),1,bmpfile);
  fread(&BitmapFileHeader.FileSize, sizeof(BitmapFileHeader.FileSize),1,bmpfile);
  fread(&BitmapFileHeader.Reserved, sizeof(BitmapFileHeader.Reserved),1,bmpfile);
  fread(&BitmapFileHeader.DataOffset, sizeof(BitmapFileHeader.DataOffset),1,bmpfile);

  if (BitmapFileHeader.Signature!=BMPMAGIC)
    {
      printf("LoadBMP: pas un fichier BMP (BMPMAGIC)\n");
      printf("LoadBMP: BMPMAGIC==%d\n",BitmapFileHeader.Signature);
      exit(-1);
    }

  fread(&BitmapInfoHeader.Size,sizeof(BitmapInfoHeader.Size),1,bmpfile); 	// lecture du info header
  fread(&BitmapInfoHeader.Width,sizeof(BitmapInfoHeader.Width),1,bmpfile); 	// lecture du info header
  fread(&BitmapInfoHeader.Height,sizeof(BitmapInfoHeader.Height),1,bmpfile); 	// lecture du info header
  fread(&BitmapInfoHeader.Planes,sizeof(BitmapInfoHeader.Planes),1,bmpfile); 	// lecture du info header
  fread(&BitmapInfoHeader.BitCount,sizeof(BitmapInfoHeader.BitCount),1,bmpfile);// lecture du info header
  fread(&BitmapInfoHeader.Compression,sizeof(BitmapInfoHeader.Compression),1,bmpfile);
  fread(&BitmapInfoHeader.ImageSize,sizeof(BitmapInfoHeader.ImageSize),1,bmpfile);
  fread(&BitmapInfoHeader.XpixelsPerM,sizeof(BitmapInfoHeader.XpixelsPerM),1,bmpfile);
  fread(&BitmapInfoHeader.YpixelsPerM,sizeof(BitmapInfoHeader.YpixelsPerM),1,bmpfile);
  fread(&BitmapInfoHeader.ColorsUsed,sizeof(BitmapInfoHeader.ColorsUsed),1,bmpfile);
  fread(&BitmapInfoHeader.ColorsImportant,sizeof(BitmapInfoHeader.ColorsImportant),1,bmpfile);

  *SizeX=BitmapInfoHeader.Width;
  *SizeY=BitmapInfoHeader.Height;

  if (BitmapInfoHeader.BitCount!=24)
    {
      printf("LoadBMP: %s n'est pas du 24 bits!\n",chemin);
      printf("LoadBMP: BMPInfoHeader.BitCount = %d \n",BitmapInfoHeader.BitCount);
      exit(-1);
    }
  if (BitmapInfoHeader.Compression!=0)
    {
      printf("LoadBMP: %s est compressé!\n",chemin);
      printf("LoadBMP: BMPInfoHeader.Compression = %ld \n",BitmapInfoHeader.Compression);
      exit(-1);
    }

  printf ("\npour le fichier %s.bmp :\t\tnb_lignes=%d  nb_colonnes=%d\n",
	  thefilename,(int )*SizeY,(int )*SizeX);


  datatoread= ((*SizeY)*(*SizeX));

  *image = (unsigned char *)calloc((3*datatoread),sizeof(unsigned char));   // allocation
  if (*image==NULL)
    {
      printf("LoadBMP: erreur d'allocation pour %s!\n",chemin);
      exit(-1);
    }

  for (i=0; i<3*datatoread; i++)
    (*image)[i]=0;    // effacement (declenche une erreur si pb d'allocation)
  RGBt.Red	= 0;
  RGBt.Green	= 0;
  RGBt.Blue	= 0;
  PaddingNum	= (*SizeX) % 4;
  if (BitmapInfoHeader.ImageSize==datatoread*3)
    BMPCase=0;
  else if (BitmapInfoHeader.ImageSize==datatoread*3+PaddingNum*(*SizeY))
    BMPCase=1;
  else if (BitmapInfoHeader.ImageSize==datatoread*4)
    BMPCase=2; // format n'existant psa a priori
  else
    {
      printf("LoadBMP: format BMP non reconnu\n");
      exit(-1);
    }

  for (iY=0; iY<(*SizeY); iY++)
    {
      for (iX=0; iX<(*SizeX); iX++)
	{
	  // lecture de la map (datatoread octets)
	  readdata= fread(&RGBt, sizeof(RGBt), (size_t) 1, bmpfile);
	  if (readdata!=1)
	    {
	      printf("LoadBMP: Erreur de lecture du fichier %s\n",chemin);
	      printf("LoadBMP: total pixels num = %ld\n", 3*datatoread );
	      printf("LoadBMP: readdata = %ld  bytes\n",readdata);
	      printf("LoadBMP: fgetpos = %d\n",fgetpos(bmpfile,&GrosBugPos));
	      //printf("LoadBMP: fgetpos returned = %d\n", GrosBugPos);
	      fclose(bmpfile);
	      exit(-1);
	    }
	  (*image)[iY*(*SizeX)+iX]=RGBt.Red;
	  (*image)[iY*(*SizeX)+iX+datatoread]=RGBt.Green;
	  (*image)[iY*(*SizeX)+iX+2*datatoread]=RGBt.Blue;
	}
      if (BMPCase==1)
	{
	  for (iPadding=0; iPadding<PaddingNum; iPadding++)
	    readdata= fread(&Padding, sizeof(Padding), (size_t) 1, bmpfile);	// lecture des 0 de padding
	}
      else if (BMPCase==2)
	{
	  for (iPadding=0; iPadding<(*SizeX); iPadding++)
	    readdata= fread(&Padding, sizeof(Padding), (size_t) 1, bmpfile);	// lecture des 0 de padding
	}
    }
  fclose(bmpfile);
  VerticalFlip(*image,(*SizeX),(*SizeY),3);
  if (DEBUGINFO>0) printf("LoadBMP: %s.bmp (%ldx%ld)\n", chemin, *SizeX, *SizeY);
  return(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// SaveBMP
// Sauve une image en BMP 24 bits
/////////////////////////////////////////////////////////////////////////////////////////////
entierns SaveBMP(char *chemin, char *thefilename, unsigned char *image, entierns SizeX, entierns SizeY, unsigned char mode)
{
  FILE 	       *bmpfile;

  BitmapFileHeaderType 	BitmapFileHeader;   // file header
  BitmapInfoHeaderType 	BitmapInfoHeader;	// info header
  entierns     	datatowrite; 	    // quantite de donnees a ecrire
  entierns     	writtendata;       	// quantite de donnees ecrites
  entierns     	iX,iY,iPadding;
  RGBTri       	RGBt;
  octetns      	Padding;
  entierns    	PaddingNum;

  VerticalFlip(image,SizeX,SizeY,mode);
  BitmapFileHeader.Signature=BMPMAGIC; 	       	// ecriture du file header
  BitmapFileHeader.FileSize=54+3*SizeX*SizeY;  	// ecriture du file header
  BitmapFileHeader.Reserved=0; 		       	// ecriture du file header
  BitmapFileHeader.DataOffset=54; 	       	// ecriture du file header

  BitmapInfoHeader.Size=40;		       	// ecriture du file header
  BitmapInfoHeader.Width=SizeX;	 	       	// ecriture du file header
  BitmapInfoHeader.Height=SizeY;	       	// ecriture du file header
  BitmapInfoHeader.Planes=1; 		       	// ecriture du file header
  BitmapInfoHeader.BitCount=24; 	       	// ecriture du file header
  BitmapInfoHeader.Compression=0; 	       	// ecriture du file header
  BitmapInfoHeader.ImageSize=3*SizeX*SizeY;    	// ecriture du file header
  BitmapInfoHeader.XpixelsPerM=3780; 	       	// ecriture du file header
  BitmapInfoHeader.YpixelsPerM=3780; 	       	// ecriture du file header
  BitmapInfoHeader.ColorsUsed=0; 	       	// ecriture du file header
  BitmapInfoHeader.ColorsImportant=0; 	        // ecriture du file header

  datatowrite = SizeX*SizeY;

  bmpfile=fopen(concat(concat(chemin,thefilename),".bmp"),"wb");
  myAssert( (bmpfile!=NULL), concat("SaveBMP: can't write to %s",chemin) );

  // ecriture du file header
  fwrite(&BitmapFileHeader.Signature, sizeof(BitmapFileHeader.Signature),1,bmpfile);
  fwrite(&BitmapFileHeader.FileSize, sizeof(BitmapFileHeader.FileSize),1,bmpfile);
  fwrite(&BitmapFileHeader.Reserved, sizeof(BitmapFileHeader.Reserved),1,bmpfile);
  fwrite(&BitmapFileHeader.DataOffset, sizeof(BitmapFileHeader.DataOffset),1,bmpfile);

  // ecriture du info header
  fwrite(&BitmapInfoHeader.Size,sizeof(BitmapInfoHeader.Size),1,bmpfile);
  fwrite(&BitmapInfoHeader.Width,sizeof(BitmapInfoHeader.Width),1,bmpfile);
  fwrite(&BitmapInfoHeader.Height,sizeof(BitmapInfoHeader.Height),1,bmpfile);
  fwrite(&BitmapInfoHeader.Planes,sizeof(BitmapInfoHeader.Planes),1,bmpfile);
  fwrite(&BitmapInfoHeader.BitCount,sizeof(BitmapInfoHeader.BitCount),1,bmpfile);
  fwrite(&BitmapInfoHeader.Compression,sizeof(BitmapInfoHeader.Compression),1,bmpfile);
  fwrite(&BitmapInfoHeader.ImageSize,sizeof(BitmapInfoHeader.ImageSize),1,bmpfile);
  fwrite(&BitmapInfoHeader.XpixelsPerM,sizeof(BitmapInfoHeader.XpixelsPerM),1,bmpfile);
  fwrite(&BitmapInfoHeader.YpixelsPerM,sizeof(BitmapInfoHeader.YpixelsPerM),1,bmpfile);
  fwrite(&BitmapInfoHeader.ColorsUsed,sizeof(BitmapInfoHeader.ColorsUsed),1,bmpfile);
  fwrite(&BitmapInfoHeader.ColorsImportant,sizeof(BitmapInfoHeader.ColorsImportant),1,bmpfile);

  Padding		= 0;
  PaddingNum	= SizeX % 4;
  for (iY=0; iY<SizeY; iY++)
    {
      for (iX=0; iX<SizeX; iX++)
	{
	  if (mode==1)
	    {
	      RGBt.Red	= (image)[iY*SizeX+iX];
	      RGBt.Green	= RGBt.Red;
	      RGBt.Blue	= RGBt.Red;
	    }
	  if (mode==3)
	    {
	      RGBt.Red	= (image)[iY*SizeX+iX];
	      RGBt.Green	= (image)[iY*SizeX+iX+datatowrite];
	      RGBt.Blue	= (image)[iY*SizeX+iX+2*datatowrite];
	    }
	  writtendata= fwrite(&RGBt, sizeof(RGBt), (size_t) 1, bmpfile);
	  // => ecriture de la map (datatoread octets)

	  if (writtendata!=1)
	    {
	      printf("SaveBMP: Erreur d'ecriture du fichier %s\n",chemin);
	      printf("SaveBMP: pixel= %ld / %ld\n", iY*SizeX+iX, datatowrite );
	      printf("SaveBMP: writtendata = %ld  bytes\n",writtendata);
	      fclose(bmpfile);
	      return(-1);
	    }
	}
      for (iPadding=0; iPadding<PaddingNum; iPadding++)
	writtendata= fwrite(&Padding, sizeof(Padding), (size_t) 1, bmpfile);
      // => lecture de la map (datatoread octets)

    }

	// Affichage :
	//printf("\nSauvegarde dans %s.bmp \n",thefilename);

  fclose(bmpfile);
  if (DEBUGINFO>0) printf("SaveBMP: %s.bmp (%ldx%ld)\n", chemin, SizeX, SizeY);
  VerticalFlip(image,SizeX,SizeY,mode);
  return(0);
}
