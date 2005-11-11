/*
*  
* $Id: delsolve.c,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
*/
/*
    fxsolv.c - C Part of Ambiguitiy Processor by Daniel Wicke (BUGH Wuppertal)

    Modified version for BRAHMS simulation package, based on FXSOLV 253.

                         Markus Elsing

    History:
    0.2   26.11.1999 - Included simplifications of FXSOLV 2.54. DW


    ME - Markus Elsing
    DW - Daniel Wicke

    ToDo:

    check SCORE/MAXSCORE
          tanagra()
	  TKincomplete stuff !!!
*/

#define VERSION       "0.2"
#define VERSIONDATE   "26.11.99"

/********************************************
 * Most important parameters
 *
 * STANDALONE   - When this is set, a standalone program is created
 *                by the compiler. Useful for testing single events.
 * N2LEVEL      - If level>N2LEVEL the simple n^2 algorithm is used.
 *                this wasn't as succesfull as I hoped, so default
 *                is 50 (=infinity)
 * RESIGNLEVEL  - If level gets equal RESIGNLEVEL the algorithm acts
 *                as if an timeout occured. This should be carfully
 *                tuned with respect to the actual TIMEOUT value, not
 *                to throw away events that would have gone through.
 * REFIT        - #define this to get substrings refitted (default).
 * TIMEOUT      - Timeout value in seconds.
 * USENDONE     - #define this to switch to counter based 'timeout'.
 *                This uses TIMEOUT*NDONEPERSEC as limit for the
 *                counter.
 ******************************************** */
/*#define STANDALONE*/
/*#define NOSUBTKS*/
#define MAXSCORECUT
/*#define MAXSCORE2*/
/*#define MAXTKCUT*/
#define N2LEVEL 50
/*#define RESIGNLEVEL 9*/
#define RESIGNLEVEL 8
#define REFIT
/*#define TIMEOUT 60*/
#define TIMEOUT 1
#define USENDONE
#define TKBONUS  100
#define NDONEPERSEC 60000

/*#define FULLNFACTORIAL*/

/* The natural logarithm of 10 */
#define LN10 (2.302585093)


/********************************************
 * Parameters for the different array sizes
 ******************************************** */
#define MAXTE 30
#define MAXTK 1000
#define MAXTKCON 1000
#define MAXTKINLIST 5000
#define MAXSUB MAXTK
#define MAXEXCL 50
#define MAXLINK 50
#define MAXEXLI 1000
#define MAXTEINLIST ((MAXTK*MAXTE)/4)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

/********************************************
 * Parameter CUTOFFNAME
 ******************************************** */
#ifdef MAXSCORECUT
#define CUTOFFNAME "Score cutoff (MAXSCORE1)"
#elif defined(MAXTKCUT)
#define CUTOFFNAME "N TK cutoff"
#else
#define CUTOFFNAME "None"
#endif

/********************************************
 * The Macro FORTRANNAME defines how 'name' has to be written,
 * so it can be seen as 'name' from FORTRAN
 ******************************************** */
#ifdef vms
#  define FORTRANNAME(name) name
#else
#  ifdef _AMIGA
#     define FORTRANNAME(name) name
#  else
#     define FORTRANNAME(name) name##_
#  endif
#endif

/*********************************************
 * The type INTEGER should reflect, what is FORTRANs idea of integer
 * The size of INTEGER should be 4 Bytes!
 ******************************************** */
#if   (USHRT_MAX==4294967295U)
#   define INTEGER unsigned short
#elif (UINT_MAX==4294967295U)
#   define INTEGER unsigned int
#elif (ULONG_MAX==4294967295U)
#   define INTEGER unsigned long
#endif

/********************************************
 * The type REAL should reflect, what is FORTRANS idea of real
 * The size of REAL should be 4 Bytes!
 ******************************************** */
#ifndef REAL
#define REAL float
#endif

/********************************************
 * Computer type specific stuff
 ******************************************** */
#ifdef _AMIGA
  /********************************************
   * AMIGA specific stuff
   ******************************************** */
#  include <dos.h>
   char const VersionStr[]="$VER: FXAMBI_FXC " VERSION " (" VERSIONDATE ")";
   long __near __stack=500000;
#  define iqsort(badlist,nbad) lqsort((long *)badlist,nbad)
   void FORTRANNAME(timex)(REAL *timeval)
   {  *timeval=(REAL)time(NULL);
   }
#  define EXIT_ERROR    30
#  undef  MAXTKCON
#  define MAXTKCON     500
#  undef  MAXTKINLIST
#  define MAXTKINLIST 1000
#  undef MAXEXCL
#  define MAXEXCL 1
#  undef MAXLINK
#  define MAXLINK 1
/*#  undef  USENDONE*/
#else
  /********************************************
   * Default stuff
   ******************************************** */
#  define chkabort()
   int
   intcmp(int const *i1,int const *i2)
   {  return *i1-*i2;
   }
#  define iqsort(l,n) qsort(l,n,sizeof(int),(cmpfkt)intcmp)
#  define EXIT_ERROR 1
#endif

/********************************************
 * Some compilers dont define the following
 ******************************************** */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/********************************************
 * Special changes for STANDALONE version
 ******************************************** */
#ifdef STANDALONE
REAL FORTRANNAME(fxprob)(REAL *a,INTEGER *b)
{  return ((REAL)1);
}
#endif

/********************************************
 * typedef's
 ******************************************** */
typedef int (*cmpfkt)(void const *,void const *);

/********************************************
 * Internal TE structure
 ******************************************** */
struct TE {
      int id;                    /* The Tanagra ID of this TE */
      int detcode;               /* Detector which made this TE */
      int nexcl;                 /* Number of excluded TEs */
      struct TE *excl[MAXEXCL];  /* Pointer to excluded TEs */
      int exclid[MAXEXCL];       /* Tanagra IDs to excluded TEs */
      int nlink;                 /* Number of linked TEs */
      struct TE *link[MAXLINK];  /* Pointer to linked TEs */
      int linkid[MAXLINK];       /* Tanagra IDs to linked TEs */
      int used;
      int labl;                  /* Tanagra Labl. */
   };

/********************************************
 *  Tanagra TKR structure
 *  (This is how a TKR bank looks like in C)
 ******************************************** */
struct TKR {
      INTEGER software_module_id;
      INTEGER det_used;
      INTEGER measurement_code;
      INTEGER track_type;
      INTEGER n_points_fitted;
      INTEGER charge;
      INTEGER mass_id;
      INTEGER ndf;
      REAL    chi2;
      REAL    length;
      REAL    x_start;
      REAL    y_start;
      REAL    z_start;
      REAL    x_end;
      REAL    y_end;
      REAL    z_end;
      REAL    xR_ref;
      REAL    yRPhi_ref;
      REAL    z_ref;
      REAL    theta;
      REAL    phi;
      REAL    Over_p;
      REAL    Covariance[15];
   };


/********************************************
 * Internal TK structure
 ******************************************** */
struct TK {
      int id;               /* The Tanagra ID of this TK or ID given by CreateSubTK() */
      int detcode;          /* Detectors used */
      unsigned int hash;    /* Hash Code: TKs with different hash are different */
      double score;         /* The TKs score */
      int nte;              /* Vector of Pointers to the TEs used in this TK */
      struct TE *te[MAXTE]; /* Pointer to the TEs used in this TK */
      struct TK *subtk[MAXTE];/* Pointer to the substring which is obtained
                               by removing te[i]. NULL if not yet calculated */
      struct TKR tkr;
   };

static struct TK EmptyTK={0};

/********************************************
 * Deklaration of external functions
 ******************************************** */

extern REAL FORTRANNAME(fxprob)(REAL *chi2,INTEGER *ndf);
extern void FORTRANNAME(timex)(REAL *);
extern void FORTRANNAME(fxrefit)(struct TKR *,struct TKR *, INTEGER *, INTEGER *, INTEGER *, INTEGER *);


/********************************************
 * Global Variables
 * (I should really make a structure from this.)
 ******************************************** */

static struct TK TKList[MAXTKINLIST];  /* The global List of all TKs we operate on */
static int TKListUsed;                 /* Number of valid entries in TKList */
static struct TE TEList[MAXTEINLIST];  /* The global List of all TEs we operate on */
static int TEListUsed;                 /* Number of valid entries in TKList */
static int LastTKID=0;                 /* LastTKID I gave to newly created TKs */

static REAL StartClockTime;            /* Time when Clock was started */
static int  NDone;                     /* Counter to gues CPU time used */

/*  2 dim array. Indices are indices into TKList.
 *  The Value:
 *        if non-negativ:   Index of the first TE in 1st TK
 *                          excluded by the 2nd TK.
 *        TKC_UNCONNECTED  if TKs are not connected
 *        TKC_UNKNOWN      if relationship is unknown.
 */
#define TKC_UNCONNECTED  (-2)
#define TKC_UNKNOWN      (-1)
static signed char TKConnection[MAXTKCON][MAXTKCON];

/* Some statistics */
static int ncombtested;
static int maxlevel;
static int fxdebug;
static int thissetmaxlevel;
static int nfxscorecalled;
static int largestsubsetsize;
static int errorbits;

/********************************************
 * ERROR codes, Detector codes and FLAGS.
 ******************************************** */

#define ERR_TIMEOUT  (1<<0)
#define ERR_RESIGNED (1<<1)
#define ERR_NDONE    (1<<2)

#define ERR_UNKNOWN_DETCODE (1<<16)
#define ERR_WRONGARGS       (1<<17)
#define ERR_TOMANYEXLI      (1<<18)
#define ERR_EMPTYTK_IN_1VD  (1<<19)

/* Explanation of FLAGs see fxsolv() */
#define FLAG_DEBUGLEVEL  ((1<<0)|(1<<1))
#define FLAG_1VDREMOVAL  (1<<2)
#define FLAG_USEN2       (1<<3)
#define FLAG_NOSUBTK     (1<<4)
/* Not really a FLAG, but used for internal steering */
#define FLAG_RECURPROT   (1<<15)

/* The detector types as given to us by fxambi */
/* ME new BRAHMS detectors, see COMMON ~ DETID */
#define TPC_CODE    5
#define FTD_CODE    3
#define ITC_CODE    4
#define VTX_CODE    6
#define FCH_CODE    7
#define SIT1_CODE   8
#define SIT2_CODE   9
#define VTX1_CODE  10
#define VTX2_CODE  11
#define VTX3_CODE  12
#define VTX4_CODE  13
#define VTX5_CODE  14
#define FTD1_CODE  15
#define FTD2_CODE  16
#define FTD3_CODE  17
#define FTD4_CODE  18
#define FTD5_CODE  19
#define FTD6_CODE  20
#define FTD7_CODE  21

/********************************************
 *  TK comparison by score as used by qsort().
 *
 *  Factor 1024 to avoid floating arithmetic
 *  in calling function and still have some
 *  digits of the probability in the comparison.
 ******************************************** */

int
tkcmp(struct TK **a, struct TK **b)
{
   return (int)(1024.0*((*a)->score-(*b)->score));
}


/********************************************
 *
 ******************************************** */

static
void startclock(void)
{   FORTRANNAME(timex)(&StartClockTime);

    if(fxdebug>0)
      printf("StartClock: %f\n",StartClockTime);
}

/********************************************
 *
 ******************************************** */
static
int timeout(void)
{   static REAL CurrentTime;
    static int count=0;

    if(count++==0)
      FORTRANNAME(timex)(&CurrentTime);

    count%=128;

    /*printf("Clock  %f  Current Time %f\n",(double)StartClockTime,(double)CurrentTime);
     */
    return (CurrentTime-StartClockTime)>TIMEOUT;
}

/********************************************
 *
 ******************************************** */
static
void
PrintTK(struct TK *tk)
{  int i,j;

   if(tk)
   {  printf("id=%d mod=%5d { ",tk->id,tk->tkr.software_module_id);
      for(i=0;i<tk->nte;i++)
         printf("%d(%d)%d ",tk->te[i]->id,tk->te[i]->labl,tk->te[i]->detcode);
      printf("}\tScore %f Prob %f\n",tk->score,
             (double)FORTRANNAME(fxprob)(&(tk->tkr.chi2),&(tk->tkr.ndf)));
      if(fxdebug>1)
      {  for(j=0;j<MAXEXCL;j++)
         {  printf("                    { ");
            for(i=0;i<tk->nte;i++)
               printf("%d(%d) ",
                      tk->te[i]->nexcl>j?tk->te[i]->excl[j]->id:0,
                      tk->te[i]->nlink>j?tk->te[i]->link[j]->id:0);
            printf("}\n");
         }
      }
   }
   else
   {  printf("(Null TK)\n");
   }

}


/********************************************
 * Create Hash Code from TK excluding the offending TE,
 * if offendingte is >= 0
 ******************************************** */
static
unsigned int
Hash(struct TK *tk,int offendingte[],int noffending)
{  unsigned int hash;
   int i,j;

   hash=0;
   for(i=0,j=0;j<=noffending;j++,i++)
   {  for(;i<offendingte[j];i++)
      {  hash<<=6;
         hash+=(unsigned int)tk->te[i];
         hash%=32771;
      }
   }
   return hash;
}

/********************************************
 *
 ******************************************** */
static
int
DetCode(struct TK *tk)
{  int i;
   int code=0;

   for(i=0;i<tk->nte;i++)
   {  code|=1<<tk->te[i]->detcode;
   }
   return code;
}


/********************************************
 * The scoring function
 ******************************************** */
static double DetScore[]={-1,-1,-1,-1,     /* ME new scoring !!! */
			  10,              /* ITC = 4 */
			  20,              /* TPC = 5 */
			  -1,
			  5,               /* FCH = 7    */
			  5,5,             /* SIT = 8,9  */
			  5,5,5,5,5,       /* VTXi=10-14 */
			  5,5,5,5,5,5,5    /* FTDi=15-21 */
			  };


static
double
Score(struct TK *tk)
{
#define SCORINGNAME "Additive scoring from array +log(Prob)"
   int i;
   double score=0,Prob;

   for(i=0;i<tk->nte;i++)
   {  if(tk->te[i]->detcode<sizeof(DetScore)/sizeof(double)
	 &&tk->te[i]->detcode>0)
      {  if(DetScore[tk->te[i]->detcode]>=0)
	    score+=DetScore[tk->te[i]->detcode];
	 else
	 {   printf("ERROR in FXSCORE: Score(): Unknown detector code: %d (score %f)!\n",
                tk->te[i]->detcode,DetScore[tk->te[i]->detcode]);
	     fflush(stdout);
	     errorbits|=ERR_UNKNOWN_DETCODE;
	 }
      }
      else
      {  printf("ERROR in FXSCORE: Score(): Unknown detector code: %d (out of range)!\n",
                tk->te[i]->detcode);
	 fflush(stdout);
         errorbits|=ERR_UNKNOWN_DETCODE;
      }
   }

   if(abs(score)>1e10)
   {  printf("Score(): score to big before log!!\n");
      fflush(stdout);
   }
   Prob=(double)FORTRANNAME(fxprob)(&(tk->tkr.chi2),&(tk->tkr.ndf));
   score += log((double)Prob)/(LN10);

   nfxscorecalled++;
   return score+TKBONUS;

}

#ifdef MAXSCORECUT
/********************************************
 * An upper bound of the scoring function
 ******************************************** */
static
double
MaxScore(struct TK *tklist[], int ntk)
{  int ite,itk;
   double maxscore;

   if(ntk>6*TEListUsed)
   {  for(ite=0;ite<TEListUsed;ite++)
         TEList[ite].used=FALSE;
   }
   else
   {  for(itk=0;itk<ntk;itk++)
      {  for(ite=0;ite<tklist[itk]->nte;ite++)
         {  tklist[itk]->te[ite]->used=FALSE;
         }
      }
   }
   maxscore=ntk*TKBONUS;

   for(itk=0;itk<ntk;itk++)
   {
      for(ite=0;ite<tklist[itk]->nte;ite++)
      {  if(!tklist[itk]->te[ite]->used)
         {  tklist[itk]->te[ite]->used=TRUE;
            if(sizeof(DetScore)/sizeof(double)<TEList[ite].detcode
               ||DetScore[TEList[ite].detcode]<0)
            {  if(!(errorbits&ERR_UNKNOWN_DETCODE))
		  printf("ERROR in FXSOLV: MaxScore() unknown detector %d\n",
			 TEList[ite].detcode);
               errorbits|=ERR_UNKNOWN_DETCODE;
            }
            else
               maxscore+=DetScore[TEList[ite].detcode];
         }
      }
   }
   /*printf("Maxscore= %f\n",maxscore);*/
   return maxscore;
}
#endif

/********************************************
 *
 ******************************************** */
static
int
TanagraDetUsed(int detcode)
{  int tanagra=0;

   /* ME this gets easy now */
   tanagra = detcode;
   return tanagra;
}

/********************************************
 *
 ******************************************** */
static
int
TEExclByTK(struct TE *te,struct TK *tk)
{  int i,j;
   struct TE *tktei;

   for(i=tk->nte-1;i>=0;i--)
   {  tktei=tk->te[i];
      if(te==tktei) return TRUE;
      for(j=tktei->nexcl-1;j>=0;j--)
         if(te==tktei->excl[j])
            return TRUE;
   }

   return FALSE;
}

/********************************************
 * Check if two TKs are excluded with each other
 * Return the first offending TE of the *second* TK
 * or -1 if not excluded.
 ******************************************** */
static
int
TKExclByTK(struct TK *tk1,struct TK *tk2)
{  int i,i1,i2,TKCon;

   if(tk1==&EmptyTK||tk2==&EmptyTK)
      return -1;


   i1=tk1-TKList;
   i2=tk2-TKList;

   if(i1<MAXTKCON && i2<MAXTKCON)
   {  TKCon=TKConnection[i1][i2];
   }
   else
   {  TKCon=TKC_UNKNOWN;
   }

   if(0<=TKCon)
      return TKCon;
   else if(TKC_UNCONNECTED==TKCon)
      return -1;
   else if(TKC_UNKNOWN==TKCon)
   {  for(i=tk2->nte-1;i>=0;i--)
      {  if(TEExclByTK(tk2->te[i],tk1))
         {  if(i1<MAXTKCON && i2<MAXTKCON)
               TKConnection[i1][i2]=i;

            return i;
         }
      }
      if(i1<MAXTKCON && i2<MAXTKCON)
         TKConnection[i1][i2]=TKC_UNCONNECTED;

      return -1;
   }


   printf("FATALERROR in FXSOLV: TKExclByTK() failed\n");
   exit(EXIT_ERROR);
}


/********************************************
 *
 ******************************************** */
static
int
TEinTK(struct TE *te,struct TK *tk)
{  int i;

   for(i=tk->nte-1;i>=0;i--)
   {  if(te==tk->te[i]) return TRUE;
   }

   return FALSE;
}

/********************************************
 *
 ******************************************** */
#ifdef STANDALONE
static
int
Refit(struct TK *tk)
{  return TRUE;
}
#else
static
int
Refit(struct TK *tk)
{
   INTEGER teidlist[MAXTE];
   INTEGER nte,ier1;
   INTEGER refitdebug;
   struct TKR tsr;
   int i;

   nte=tk->nte;
   for(i=0;i<nte;i++)
      teidlist[i]=tk->te[i]->id;

   tsr=tk->tkr;
   /* A TSR has a different measurement code: */
   if(tsr.measurement_code&1)
      tsr.measurement_code=7;
   else
      tsr.measurement_code=5;

   fflush(stdout);
   refitdebug=fxdebug;		
   FORTRANNAME(fxrefit)(&tsr,&tk->tkr,&nte,teidlist,&ier1,&refitdebug);

   return 0==ier1;
}
#endif

/********************************************
 *
 ******************************************** */
static
int
TKIncomplete(struct TK *tk)
{  int const incompletecode[]=
      {  /* single detector except TPC */
           0
      };
   int j;

   if(!tk->detcode) return TRUE;          /* Empty TK is incomplete */

   if(3>tk->nte)                          /* 2 hits without TPC is incomplete */
     {
       int notpc = TRUE;
       int i;
       for(i=0;i<tk->nte;i++)
	 {
	   if (TPC_CODE==tk->te[i]->detcode) notpc = FALSE;
	 }
       if(notpc)
	 {
	   /* printf("No TPC on SubTK(), reject: NTE ~ %d\n",tk->nte); */
	   return TRUE;
	 }
     }

   for(j=0;j<sizeof(incompletecode)/sizeof(int);j++)
     {  if(incompletecode[j]==tk->detcode)
       {
	 return TRUE;
       }
     }

   return FALSE;
}

/********************************************
 *
 ******************************************** */
static
struct TK *
CreateSubTK(struct TK *tk, int offendingte[], int noffending, struct TK *tkvec,int *ntk)
{  struct TK  *newtk;
   int i,j;

   if(*ntk>=MAXTKINLIST)
   {  printf("TKList overflow.\n");
      return &EmptyTK;    /* KH */
      /*exit(EXIT_ERROR);    KH */
   }

   if(fxdebug>1)
      printf("CreateSubTK(): noffending = %d\n",noffending);

   newtk= tkvec+(*ntk)++;
   newtk->id= --LastTKID;
   newtk->nte= tk->nte-noffending;

   if(*ntk-1<MAXTKCON)
      for(i=0;i<*ntk;i++)
      {  TKConnection[i][*ntk-1]=TKC_UNKNOWN;
         TKConnection[*ntk-1][i]=TKC_UNKNOWN;
      }

   for(i=0,j=0;j<=noffending;j++,i++)
   {
      for(;i<offendingte[j];i++)
      {  newtk->te[i-j]  =tk->te[i];
         newtk->subtk[i-j]=NULL;
      }
   }

   newtk->hash=Hash(newtk,&newtk->nte,0);
   newtk->detcode=DetCode(newtk);
   newtk->tkr=tk->tkr;
   if(TKIncomplete(newtk))
   {  (*ntk)--;
      return &EmptyTK;
   }
   else
   {
      newtk->tkr.software_module_id+=10000;
      newtk->tkr.det_used=TanagraDetUsed(newtk->detcode);

#ifdef REFIT
      if(!Refit(newtk))
      {  (*ntk)--;

         if(fxdebug>1)
	 {  printf("WARNING in FXSOLV: CreateSubTK(): Refit() failed. Returning EmptyTK.\n");
	    PrintTK(newtk);
	    printf("TKIncomplete(newtk)=%d\n",TKIncomplete(newtk));
	 }
         return &EmptyTK;
      }
#endif

      newtk->score=Score(newtk);

      return newtk;
   }
}


/********************************************
 * Find the TK with the offendingte removed in the TKvec.
 ******************************************** */
static
struct TK *
FindSubTK(struct TK *tk, int offendingte[], int noffending, struct TK tkvec[],int ntk)
{  int k,i,j,sametk,hash;

   hash=Hash(tk,offendingte,noffending);

   for(k=ntk-1;k>=0;k--)
   {  if(tkvec[k].hash==hash&&tkvec[k].nte==tk->nte-noffending)
      {  sametk=TRUE;

         for(i=0,j=0;sametk&&j<=noffending;j++,i++)
         {  for(;sametk&&i<offendingte[j];i++)
               if(!TEinTK(tk->te[i],tkvec+k))
               {  sametk=FALSE;
               }
         }
         if(sametk) return tkvec+k;
      }
   }
   return NULL;
}

/********************************************
 * Get the TK with the offendingte removed.
 ******************************************** */
static
struct TK *
SubTK(struct TK *tk, int Offendingte)
{  struct TK *subtk;
   int i,j,noffending,offendingte[MAXLINK+2];

   if(subtk=tk->subtk[Offendingte])
   {  return subtk;
   }
   else
   {
      /* Setup offendingte[] containing the original offendingte
         as well as TEs linked to this TE. */
      noffending=1;
      offendingte[0]=Offendingte;
      if(tk->te[offendingte[0]]->nlink)
      {  for(i=0;i<tk->te[offendingte[0]]->nlink;i++)
         {  struct TE *lte;

            lte=tk->te[offendingte[0]]->link[i];
            for(j=0;j<tk->nte;j++)
            {  if(lte==tk->te[j])
               {  offendingte[noffending++]=j;
                  break;
               }
            }
         }
         iqsort(offendingte,noffending);
      }
      offendingte[noffending]=tk->nte;

      if(subtk=tk->subtk[offendingte[0]])
      {  return subtk;
      }
      if(subtk=FindSubTK(tk,offendingte,noffending,TKList,TKListUsed))
      {  tk->subtk[offendingte[0]]=subtk;
         return subtk;
      }
      else
      {  subtk=CreateSubTK(tk,offendingte,noffending,TKList,&TKListUsed);
         tk->subtk[offendingte[0]]=subtk;
         return subtk;
      }
   }
}

/*****************************************************************
 * Get the TK which has all TEs of first removed from tk.
 *
 * This routine takes care of splits:
 *   A SubTK of TK1 with respect to TK2 is only returned,
 *   if the SubTK of TK2 with respect to TK1 is also valid.
 ***************************************************************** */
static
struct TK *
TKWithout(struct TK *first, struct TK *tk, int flags)
{   int Offendingte;
    struct TK *subtk;


    if(0<=(Offendingte=TKExclByTK(first,tk)))
    {   if(flags&FLAG_NOSUBTK)
           return &EmptyTK;
        else
        {  subtk=tk;
           do
           {  subtk=SubTK(subtk,Offendingte);
           }while(0<=(Offendingte=TKExclByTK(first,subtk)));
#if 0
           return subtk;
#else
           if(!(flags&FLAG_RECURPROT))
           {  if(TKWithout(tk,first,flags|FLAG_RECURPROT)!=&EmptyTK)
                 return subtk;
              else
              {  /* printf("FXSOLV: TKWithout(): Suppressed substring\n");
                 PrintTK(first);
                 PrintTK(tk); */
                 return &EmptyTK;
              }
           }
           else
              return subtk;
#endif
        }
    }
    else
    {   return tk;
    }
}


/********************************************
 *
 ******************************************** */
static
void
CopyWithout(int first,struct TK *tklist[] , struct TK *tmptk[],int ntk, int flags)
{
   int k;
   struct TK *FirstTK;

   FirstTK=tklist[first];
   for(k=0;k<first;k++)
      tmptk[k]  =TKWithout(FirstTK,tklist[k],flags);
   for(k=first+1;k<ntk;k++)
      tmptk[k-1]=TKWithout(FirstTK,tklist[k],flags);

   return;
}

/********************************************
 *  Remove the TKs with indices listed in list[] from tk[].
 *
 *  List must have nlist+1 entries. The last entry must be larger
 *  than any other entry.
 *
 *  !!! list must be sorted !!!
 ******************************************** */
static
void
RemoveList(struct TK *tk[],int *ntk, int list[], int nlist)
{  int i,j,offset;

   offset=0;
   for(i=0;i<nlist;i++)
   {  if(list[i]!=list[i+1])
      {  offset++;
         for(j=list[i]+1;j<list[i+1];j++)
         {  tk[j-offset]=tk[j];
         }
      }
   }
   *ntk-=offset;
}


/********************************************
 * Remove Incomplete TKs from the list tk[].
 * Set flag to 1, if you want warnings about
 * removed TKs
 ******************************************** */
static
void
RemoveIncomplete(struct TK *tk[],int *ntk, int flag)
{  int i;
   int badlist[MAXTK+1],nbad;

   nbad=0;
   for(i=0;i<*ntk;i++)
   {  if(TKIncomplete(tk[i]))
      {  badlist[nbad++]=i;
         /* Print a warning if flag==1 and TK from barrel or if flag==2 */
         if((flag==1&&tk[i]->tkr.software_module_id>99)||(flag==2))
         {  printf("WARNING in FXSOLV: RemoveIncomplete() removed TK\n");
            PrintTK(tk[i]);
         }
      }
   }
   badlist[nbad]=*ntk;
   RemoveList(tk,ntk,badlist,nbad);
}


/********************************************
 *  Check 2 TKs to be substrings or not.
 *
 *  Returns 3 if TKs are equal,
 *          2 if tk2 is substring of tk1,
 *          1 if tk1 is substring of tk2,
 *          else it returns 0.
 ******************************************** */
static
int
IsSubstring(struct TK *tk1, struct TK *tk2)
{  int i;
   int result=3;

   if( (!tk1)|(!tk2) )
   {    printf("ERROR in FXSCORE: IsSubstring() tk1=%d, tk2=%d\n",tk1,tk2);
        errorbits|=ERR_WRONGARGS;
        return 0;
   }
   for(i=0;i<tk1->nte;i++)
   {  if(!TEinTK(tk1->te[i],tk2))
      {  result-=1;
         break;
      }
   }
   for(i=0;i<tk2->nte;i++)
   {  if(!TEinTK(tk2->te[i],tk1))
      {  result-=2;
         break;
      }
   }
   return result;
}

/********************************************
 *
 ******************************************** */
static
void
RemoveSubstrings(struct TK *tk[],int *ntk)
{  int i,j,iadd,substr;
   int badlist[MAXTK+1],nbad;


   for(i=0;i<*ntk;i+=iadd)
   {  nbad=0;
      if(0==tk[i]->nte)  /* Remove EmptyTK */
      {  badlist[nbad++]=i;
         iadd=0;
      }
      else
      {  iadd=1;
         for(j=i+1;j<*ntk;j++)
         {  substr=IsSubstring(tk[i],tk[j]);
            if(1==substr)
            {  badlist[nbad++]=i;
               iadd=0;
            }
            else if(2==substr)
            {  badlist[nbad++]=j;
            }
            else if(3==substr)     /* Strings are equal */
            {  if(tk[j]->score<tk[i]->score)   /* ME remove the one with */
               {  badlist[nbad++]=j;           /* smaller probability    */
               }                               /* ( VD multiplexing )    */
               else
               {  badlist[nbad++]=i;
                  iadd=0;
               }
            }
         }
         /* badlist[] must be sorted!! Sorting is needed in the case,*/
         /* where tk[i] is the substring.                            */
         /* Not sorting will cause wrong results/crashes, when tk[i] */
         /* is entered multiple times into badlist[].                */
         iqsort(badlist,nbad);
      }
      badlist[nbad]=*ntk;
      RemoveList(tk,ntk,badlist,nbad);
   }
}

/*********************************************************
 *  This reorders tk[] such that connected tracks are placed
 *  adjoiningly.
 *
 *  The number of connected sets is returned in nsets, their
 *  positions are returned in start[]. The nth set starts at
 *  position start[n] and ends at position start[n+1]-1.
 ********************************************************* */
static
void
BuildSubsets(struct TK *tk[], int ntk, int start[], int *nsets)
{
#ifdef FULLNFACTORIAL
   *nsets=1; /* All connected unless we have RemoveSubstrings running */
   start[0]=0;
   start[1]=ntk;
#else
   int itk1,itk2,set;

   for(set=0,start[set]=0;start[set]+1<=ntk;set++)
   {
      start[set+1]=start[set]+1;
      for(itk1=start[set];itk1<start[set+1];itk1++)
      {
         for(itk2=start[set+1];itk2<ntk;itk2++)
         {  if(0<=TKExclByTK(tk[itk1],tk[itk2]))
            {
               /* Swap tk[start[set+1]], itk2  and increase start[set+1] */
               if(start[set+1]!=itk2)
               {  struct TK *tmp;

                  tmp=tk[start[set+1]];
                  tk[start[set+1]]=tk[itk2];
                  tk[itk2]=tmp;
               }
               start[set+1]++;
            }
         }
      }
      if(start[set+1]-start[set]>largestsubsetsize)
         largestsubsetsize=start[set+1]-start[set];
   }
   *nsets=set;
#endif
}

/********************************************
 *
 ******************************************** */
static
int
TKIndexById(int id, struct TK tklist[],int n)
{  int index;

   for(index=n-1;index>=0&&tklist[index].id!=id;index--);
   return index;
}

/********************************************
 *
 ******************************************** */
static
int
TKIndexByPtr(struct TK *tk, struct TK *tklist[],int n)
{  int index;

   for(index=n-1;index>=0&&tklist[index]!=tk;index--);
   return index;
}

/********************************************
 *
 ******************************************** */
static
int
TEIndexById(int id, struct TE telist[], int n)
{  int index;

   for(index=n-1;index>=0&&telist[index].id!=id;index--);
   return index;
}


/********************************************
 *
 ******************************************** */
static
void
PrintSolution(struct TK *solution[],int nsolution,double score)
{  int i;
   for(i=0;i<nsolution;i++)
   {  PrintTK(solution[i]);
   }
   printf("Score: %f\n",score);
}


/********************************************
 *
 ******************************************** */
static
void
SolvN2(struct TK *tklist[],int ntk,struct TK *sol[],int *nsol, double *score, int flags, int *error)
{
   int besti,i,tmpntk;
   double bestscore;
   struct TK *tmptk[MAXTK];


   chkabort();

   memcpy(tmptk,tklist,ntk*sizeof(struct TK *));
   tmpntk=ntk;
   *score=0;
   *nsol=0;

   if(fxdebug>1)
      printf("*** SolvN2 starting *** \n");

   while(tmpntk>0)
   {  /* *** Find best TK *** */
      bestscore=-1e10;
      for(i=0;i<tmpntk;i++)
      {  if(tmptk[i]->score>bestscore)
         {  besti=i;
            bestscore=tmptk[besti]->score;
         }
      }
      /* *** Add best to sol[] *** */
      sol[(*nsol)++]=tmptk[besti];
      *score+=bestscore;


      /* *** Remove used TEs from TKs
             and remove incomplete and substrings *** */

      CopyWithout(besti,tmptk,tmptk,tmpntk,flags);
      tmpntk--;

      RemoveSubstrings(tmptk,&tmpntk);

   }

}

/********************************************
 *
 ******************************************** */
static
void
SolvSet(struct TK *tklist[],int ntk,struct TK *sol[],int *nsol, double *score,
        double minscore, int flags, int *error)
{  static level=0;
   struct TK *tmpsol[MAXTK];
   int start[MAXSUB+1],ntmpsol,tmpntk,nsets;
   double tmpscore,subsolscore;
   int nsubsol;
   struct TK *tmptk[MAXTK];
   int first,set;

   chkabort();

   level++;
   if(level>maxlevel)       maxlevel=level;
   if(level>thissetmaxlevel)thissetmaxlevel=level;

   *nsol=0;
   *score=0;

#if 1
#ifdef MAXSCORECUT
   if(MaxScore(tklist,ntk)<minscore)
   {  level--;
      return;
   }
#elif defined(MAXTKCUT)
   /*  printf("%2d (min %f) Total ",level,minscore);
    *  printf("*****************\n");
    */
   if(ntk<minscore)
   {  level--;
      return;
   }

#endif
#endif

   NDone+=ntk;
   for(first=0;first<ntk;first++)
   {
      tmpsol[0]=tklist[first];
      ntmpsol=1;
      tmpscore=tklist[first]->score;

      CopyWithout(first,tklist,tmptk,ntk,flags);
      tmpntk=ntk-1;

      RemoveSubstrings(tmptk,&tmpntk);
      BuildSubsets(tmptk,tmpntk,start,&nsets);

/*       printf("SolvSet(%d tracks)\n",ntk);
 *       printf("SolvSet(%d tracks)\n",tmpntk);
 *       printf("Level %d, First %d (id=%d)\n",level,first,tklist[first]->id);
 *       printf("NSets %d\n",nsets);
 *           printf("TKListUsed=%d\n",TKListUsed);
 *           fflush(stdout);
 *       for(set=0;set<=nsets;set++)
 *           printf(" %3d\n",start[set]);
 *       printf("----\n");
 *       for(i=0;i<tmpntk;i++)
 *       {  printf(" %d ",i);
 *          PrintTK(tmptk[i]);
 *       }
 *
 *       printf("====\n");
 *
 *
 *
 *       maxtmpscore=MaxScore(tmptk,tmpntk)+tmpscore;
 *       printf("",maxtmpscore);
 *       printf("Level %d ",level);
 *       printf("need=%f max=%f score=%f\n",needscore,maxtmpscore,*score);
 */
#ifdef MAXSCORECUT
     if(MaxScore(tmptk,tmpntk)>minscore-tmpscore)
#elif defined(MAXTKCUT)
     if(tmpntk>=minscore-ntmpsol)
#endif
      for(set=0;set<nsets;set++)
      {
#ifndef USENDONE
         if(timeout())
         {  *error |= ERR_TIMEOUT;
            level--;
            return;
         }
#else
         if(NDone>TIMEOUT*NDONEPERSEC)
         {  *error |= ERR_NDONE;
            level--;
            return;
         }
#endif
         if(thissetmaxlevel>=RESIGNLEVEL)
         {  *error |= ERR_RESIGNED;
            level--;
            return;
         }
#if 0
#ifdef MAXSCORECUT
         if(MaxScore(tmptk+start[set],start[nsets]-start[set])<minscore-tmpscore)
         {  break;
         }
#elif defined(MAXTKCUT)
         if(start[nsets]-start[set]<minscore-ntmpsol)
         {  break;
         }
#endif
#endif
         if(start[set+1]-start[set]==0)
         {  ncombtested++;
         }
         else if(start[set+1]-start[set]==1)
         {  ncombtested++;
            tmpsol[ntmpsol++]=tmptk[start[set]];
            tmpscore+=tmptk[start[set]]->score;
         }
         else
         {
            if(level>=N2LEVEL)
            {  nsubsol=0;
               SolvN2(tmptk+start[set],start[set+1]-start[set],tmpsol+ntmpsol,&nsubsol,&subsolscore,flags,error);
            }
            else
            {  nsubsol=0;
               SolvSet(tmptk+start[set],start[set+1]-start[set],tmpsol+ntmpsol,
                       &nsubsol,&subsolscore,
#ifdef MAXSCORECUT
                       minscore-tmpscore,
#elif defined(MAXTKCUT)
                       minscore-tmpntk-start[nsets]+start[set+1],
#else
                       0.0,
#endif
                       flags,error);
            }
            tmpscore+=subsolscore;
            ntmpsol+=nsubsol;
         }
      }


      if(tmpscore>=*score)
      {  memcpy(sol,tmpsol,ntmpsol*sizeof(*tmpsol));
         *nsol=ntmpsol;
         *score=tmpscore;

#ifdef MAXSCORECUT
         if(tmpscore>minscore)
            minscore=tmpscore;
#elif defined(MAXTKCUT)
         if(ntmpsol>minscore)
            minscore=ntmpsol;
#endif
      }
   }
#ifdef MAXSCORE
   /*printf("%2d (min %f) Final %f\n",level,minscore,*score);*/
#endif
   level--;
}


/****************************************************************************************
 *                            Solv unconnected set                                      *
 ****************************************************************************************/
static
void
SolvUnconSet(struct TK *tklist[],int ntk,struct TK *sol[],int *nsol, double *score,
             int flags, int *fxerror )
{  static level=0;
   struct TK *tmpsol[MAXTK];
   int start[MAXSUB+1],ntmpsol,tmpntk,nsets;
   double tmpscore,subsolscore;
   int nsubsol;
   struct TK *tmptk[MAXTK];
   int set,localerror=0;
   int *error=&localerror;

   chkabort();

   if(ntk>MAXTK)
   {  printf("ERROR in FXSOLV: SolvUnconSet(): ntk > MAXTK\n");
      exit(EXIT_ERROR);
      return;
   }
   *score=0;

   memcpy(tmptk,tklist,ntk*sizeof(*tklist));
   tmpntk=ntk;
   ntmpsol=0;
   tmpscore=0;

   RemoveIncomplete(tmptk,&tmpntk,1);

   if(*nsol!=0)
   {  int first;

      for(ntmpsol=0;ntmpsol<*nsol;ntmpsol++)
      {  first=TKIndexByPtr(sol[ntmpsol],tmptk,tmpntk);
         CopyWithout(first,tmptk,tmptk,tmpntk,0);
         *score+=tmptk[first]->score;
         tmpntk--;
         RemoveSubstrings(tmptk,&tmpntk);
      }

   }
   else
   {
      RemoveSubstrings(tmptk,&tmpntk);
   }

   BuildSubsets(tmptk,tmpntk,start,&nsets);

   if(fxdebug>0)
   {  printf("SolvUnconSet(%d tracks)\n",ntk);
      printf("SolvUnconSet(%d tracks)\n",tmpntk);
      printf("NSets %d\n",nsets);
   }

   /*for(set=0;set<=nsets;set++)
    *   printf(" %3d\n",start[set]);
    *printf("----\n");
    *for(i=0;i<tmpntk;i++)
    *{  printf(" %d ",i);
    *   PrintTK(tmptk[i]);
    *}
    *printf("====\n");
    */

   for(set=0;set<nsets;set++)
   {  if(start[set+1]-start[set]==0)
      {  ncombtested++;
      }
      else if(start[set+1]-start[set]==1)
      {  ncombtested++;
         tmpsol[ntmpsol++]=tmptk[start[set]];
         tmpscore+=tmptk[start[set]]->score;
      }
      else
      {
         if(level>=N2LEVEL||flags&FLAG_USEN2)
         {  SolvN2(tmptk+start[set],start[set+1]-start[set],
                   tmpsol+ntmpsol,&nsubsol,&subsolscore,flags,error);
            *fxerror |= *error;
         }
         else
         {  startclock();
            NDone=0;
            *error&=~(ERR_TIMEOUT|ERR_RESIGNED|ERR_NDONE);

#if defined(MAXSCORECUT)|defined(MAXTKCUT)
            qsort(tmptk+start[set],start[set+1]-start[set],sizeof(*tmptk),
                  (cmpfkt)tkcmp);
#endif

            thissetmaxlevel=0;
            nsubsol=0;
            SolvSet(tmptk+start[set],start[set+1]-start[set],
                    tmpsol+ntmpsol,&nsubsol,&subsolscore,0.0,flags,error);
            *fxerror |= *error;

#if 0
            {  REAL CurrentTime;

               FORTRANNAME(timex)(&CurrentTime);
               printf("FXSOLV_TEST: NDone on set %d, using SolvN2()."
                      " Time: %f. NDone= %d. This sets maxlevel: %d\n",
                      set,CurrentTime-StartClockTime,NDone,thissetmaxlevel);
            }
#endif

            if((*error)&(ERR_TIMEOUT|ERR_RESIGNED|ERR_NDONE))
            {  REAL CurrentTime;

               FORTRANNAME(timex)(&CurrentTime);
               if((*error)&ERR_TIMEOUT)
                  printf("FXSOLV: Set %d timed out, using SolvN2()."
                         " Time: %f. NDone= %d. This sets maxlevel: %d\n",
                         set,CurrentTime-StartClockTime,NDone,thissetmaxlevel);
               else if((*error)&ERR_RESIGNED)
                  printf("FXSOLV: Resigned on set %d, using SolvN2()."
                         " Time: %f. NDone= %d. This sets maxlevel: %d\n",
                         set,CurrentTime-StartClockTime,NDone,thissetmaxlevel);
               else if((*error)&ERR_NDONE)
                  printf("FXSOLV: NDone on set %d, using SolvN2()."
                         " Time: %f. NDone= %d. This sets maxlevel: %d\n",
                         set,CurrentTime-StartClockTime,NDone,thissetmaxlevel);
               else
                  printf("Internal ERROR in FXSOLV!!\n");

               nsubsol=0;
               SolvN2(tmptk+start[set],start[set+1]-start[set],
                      tmpsol+ntmpsol,&nsubsol,&subsolscore,flags,error);
               *fxerror |= *error;
            }
         }
         tmpscore+=subsolscore;
         ntmpsol+=nsubsol;
      }
   }
   memcpy(sol,tmpsol,ntmpsol*sizeof(*tmpsol));
   *nsol=ntmpsol;
   *score=tmpscore;
}


/************************************************************
 *  Print banner (version, date and most important parameters)
 ************************************************************ */
void
FORTRANNAME(fxsolvversion)(void)
{
   printf("FXSOLV Version " VERSION " (" VERSIONDATE ")\n");
#ifdef USENDONE
   printf("FXSOLV TIMEOUT = %ds (counter based)\n",TIMEOUT);
#else
   printf("FXSOLV TIMEOUT = %ds\n",TIMEOUT);
#endif
   printf("FXSOLV N2LEVEL = %d\n",N2LEVEL);
   printf("FXSOLV RESIGNLEVEL = %d\n",RESIGNLEVEL);
   printf("FXSOLV %s, +%dperTK.\n",SCORINGNAME,TKBONUS);
   printf("FXSOLV Cutoff Algorithm: %s",CUTOFFNAME);

#ifndef STANDALONE
   if(sizeof(INTEGER) != sizeof(REAL))
   {  printf("\nERROR in FXSOLVVERSION: REAL and INTEGER are of different size!!\n");
      printf(  "                        Check #define's in fxsolv.c !\n");
      exit(EXIT_ERROR);
   }
#endif
   printf("\n");
   fflush(stdout);
   return;
}

/***************************************************
 * Add Exclusion or Link to TE structure
 * Take care not to add things twice.
 *************************************************** */
#define FLAG_EXCL 1
#define FLAG_LINK 2

static
void
AddExclorLink(struct TE *te, int addid, int flag)
{  int k, *nlist, *list;
   char *name;

   if(addid)
   {
      if(FLAG_EXCL==flag)
      {  nlist=&te->nexcl;
         list=te->exclid;
         name="exclusions";
      }
      else if(FLAG_LINK==flag)
      {  nlist=&te->nlink;
         list=te->linkid;
         name="links";
      }
      else
      {  printf("INTERNAL ERROR in FXSOLV: AddExclorLink() flag=%d.\n",flag);
         exit(EXIT_ERROR);
      }

      for(k=*nlist-1;k>=0;k--)
      {  if(addid==list[k]) break;
      }
      if(k<0)
      {  if(fxdebug>0)
	/* printf("Adding %s to %d in %d.\n",name,list,te); */
         if(*nlist<MAXEXCL)
            list[(*nlist)++]=addid;
         else
         {  printf("ERROR in FXSOLV: AddExclorLink(): To many %s.\n",name);
            errorbits|=ERR_TOMANYEXLI;
         }
      }
   }
}

/************************************************************
 *  FXSOLV  DELANA entry point for this ambiguity processor
 *  (valid since version 2.43)
 *
 *  ntk         - Number of TKs to solv.
 *                (Will be increase, when subTKs are used in the solution)
 *  tkid        - The Tanagra IDs of the TKs to solv
 *  tknte       - The number of TEs in each TK.
 *  tkteid      - The Tanagra IDs for each TE in each TK. (2dim)
 *  tktedet     - The detector identifier for each TE. (2dim)
 *  tkteexli    - The Tanagra IDs of an excluded TE for each TE. (3dim)
 *              - 3rd dim is sized 4 integers containing
 *                (nexcl, exclidx, nlink, linkidx), where nexcl/nlink are
 *                the number of excl/links for the corrosponding TE and
 *                exclidx/linkidx are indices into tkteexliids pointing
 *                to the corrospondig excluded/linked TE-ids.
 *  teexliids   - The Tanagra IDs of excluded or linked TEs. (1dim)
 *  tktelabl    - The Tanagra lables fro each TE. (2dim)
 *  tkr         - TKR structure (2 dim).
 *  ntemax      - Max number of TEs per TK. In fact this is the size of the
 *                second (last but one in fortran) dimension of the
 *                multidim arrays.
 *  ntkmax      - Max number of TKs allowed in these arrays.
 *  tkrsize     - Size of TKR structure. In fact this is the size of the
 *                first (last in fortran) dimension in the 2dim array tkr.
 *  solutionindexlist
 *              - Indices into the above arrays describing the solution
 *  nsolution   - Number of TKs in the solution
 *  fxflags     - Integer indicating the debug level
 *                Bit 0: 0=no debug,1=some debug
 *                Bit 1: 1=full debug
 *                Bit 2: 1=1VD removal after solving the event.
 *                Bit 3: 1=Use SolvN2() for solving Subsets.
 *                Bit 4: 1=Do not create subTKs.
 *                All other bits must be set to 0.
 *  fxerror     - Integer indicating warnings (bit 0 to 15)
 *                or error conditions (bit 16 to 31).
 *                (Bit  0 set: Timeout occured,
 *                 bit  1 set: Resign condition occured,
 *                 bit  2 set: NDone-resign condition occured.
 *                 bit 16 set: Unknown detector code.
 *                 bit 17 set: Wrong Arguments.
 *                 bit 18 set: To many exclusions/links.
 ************************************************************ */
void
FORTRANNAME(fxsolv)(int *ntk, int tkid[], int tknte[],
        int tkteid[],  int tktedet[], int tkteexli[], int teexliids[],
        int tktelabl[], int tkr[],
        const int *ntemax, const int *ntkmax, const int *tkrsize,
        int solutionindexlist[],int *nsolution, int *fxflags, int *fxerror)
{

   struct TK *tkptrlist[MAXTK];
   struct TK *solution[MAXTK];

   double score;
   int i,j,k,teid;
   time_t start_time,end_time;;

   errorbits=0;
   fxdebug=*fxflags&FLAG_DEBUGLEVEL;
   start_time=time(NULL);

   if(*tkrsize*sizeof(INTEGER) != sizeof(struct TKR))
   {  printf("ERROR in FXSOLV: Internal size of TKR does not fit the size "
             "received from the steerings.\n");
      printf("Internal Size is %d, Received Size is %d.\n",
             sizeof(struct TKR),*tkrsize*sizeof(INTEGER));
      exit(EXIT_ERROR);
   }
   if(*ntk>MAXTK)
   {  printf("ERROR in fxsolv(): *ntk > MAXTK\n");
      exit(EXIT_ERROR);
      return;
   }


   /* ********** Read the Event ********** */
   TEListUsed=0;
   for(i=0;i<*ntk;i++)
   {
      TKList[i].id =tkid[i];
      TKList[i].nte=tknte[i];
      memcpy(&TKList[i].tkr,tkr+i**tkrsize,*tkrsize*sizeof(INTEGER));

      for(j=0;j<TKList[i].nte;j++)
      {  teid=tkteid[j+i**ntemax];
         k=TEIndexById(teid,TEList,TEListUsed);
         if(k==-1)
         {  TEList[TEListUsed].id=teid;
            TEList[TEListUsed].detcode=0;
            TEList[TEListUsed].nexcl=0;
            TEList[TEListUsed].nlink=0;
            TKList[i].te[j]=TEList+TEListUsed;
            TEListUsed++;
         }
         else
         {  TKList[i].te[j]=TEList+k;
         }
         TKList[i].subtk[j]=NULL;
      }

      /* **** Read DetCodes and Labls of TEs **** */
      for(j=0;j<TKList[i].nte;j++)
      {  TKList[i].te[j]->detcode=tktedet[j+i**ntemax];
         TKList[i].te[j]->labl=tktelabl[j+i**ntemax];

      }

      /* **** Read exclusions/links of TEs **** */
      for(j=0;j<TKList[i].nte;j++)
      {  int nexli,exliidx;
         nexli=tkteexli[0+j*4+i*4**ntemax];
         exliidx=tkteexli[1+j*4+i*4**ntemax];
         for(k=0;k<nexli;k++)
         {  teid=teexliids[exliidx+k];
            AddExclorLink(TKList[i].te[j],teid,FLAG_EXCL);
         }
         nexli=tkteexli[2+j*4+i*4**ntemax];
         exliidx=tkteexli[3+j*4+i*4**ntemax];
         for(k=0;k<nexli;k++)
         {  teid=teexliids[exliidx+k];
            AddExclorLink(TKList[i].te[j],teid,FLAG_LINK);
         }
      }

      TKList[i].hash=Hash(TKList+i,&((TKList+i)->nte),0);
      TKList[i].detcode=DetCode(TKList+i);
      TKList[i].score=Score(TKList+i);


      tkptrlist[i]=&TKList[i];
   }
   TKListUsed=i;

   /* *** Change the exclusion and link list of TEs
          to contain pointers instead Id numbers **** */
   for(i=0;i<TEListUsed;i++)
   {  for(j=0;j<TEList[i].nexcl;j++)
      {  teid=TEList[i].exclid[j];
         if(teid!=0)
         {  k=TEIndexById(teid,TEList,TEListUsed);
            if(k==-1)
            {  /* Create TE entries for TEs not part in a TK */
               TEList[TEListUsed].id=teid;
               TEList[TEListUsed].detcode=0;
               TEList[TEListUsed].nexcl=0;
               TEList[TEListUsed].nlink=0;
               TEList[i].excl[j]=TEList+TEListUsed;
               TEListUsed++;
            }
            else
               TEList[i].excl[j]=TEList+k;
         }
      }
      for(j=0;j<TEList[i].nlink;j++)
      {
         teid=TEList[i].linkid[j];
         if(teid!=0)
         {  k=TEIndexById(teid,TEList,TEListUsed);

            if(k==-1)
            {  /* Create TE entries for TEs not part in a TK */
               TEList[TEListUsed].id=teid;
               TEList[TEListUsed].detcode=0;
               TEList[TEListUsed].nexcl=0;
               TEList[TEListUsed].nlink=0;
               TEList[i].link[j]=TEList+TEListUsed;
               TEListUsed++;
            }
            else
               TEList[i].link[j]=TEList+k;
         }
      }
   }

   /* ********** Clear TKConnected[][] ********** */
   /*    i=TKListUsed<MAXTKCON?TKListUsed:MAXTKCON;
    *    memset(TKConnection,TKC_UNKNOWN,i*i*sizeof(TKConnection[0][0]));
    */
   for(i=0;i<TKListUsed&&i<MAXTKCON;i++)
      for(j=0;j<TKListUsed&&j<MAXTKCON;j++)
         TKConnection[i][j]=TKC_UNKNOWN;


   if(fxdebug>0)
   {  /* ***** Print what we got ***** */
      if(*fxflags&FLAG_USEN2) printf("FXSOLV: Using SolvN2() only!!\n");
      if(*fxflags&FLAG_NOSUBTK) printf("FXSOLV: Creating no SubTKs!!\n");
      printf("Event to solve:\n");
      for(i=0;i<TKListUsed;i++)
      {  PrintTK(TKList+i);
      }
   }

   /* ********** Solv ********** */
   ncombtested=0;
   maxlevel=0;
   nfxscorecalled=0;
   largestsubsetsize=0;
   *nsolution=0;

   SolvUnconSet(tkptrlist,TKListUsed,solution,nsolution,&score,*fxflags,&errorbits);

   if(fxdebug>0)
   {  printf("Solution: %d Tracks, Score %f\n",*nsolution,score);
      PrintSolution(solution,*nsolution,score);
      printf("Maxlevel = %d\n",maxlevel);
      printf("%d combinations tested.\n",ncombtested);
      printf("%d times Score() called.\n",nfxscorecalled);
      printf("%d TKs in TKList.\n",TKListUsed);
      printf("%d TKs in largest connected set.\n",largestsubsetsize);

      for(i=0;i<*nsolution;i++)
         for(j=0;j<*nsolution;j++)
         {  if(i!=j&&0<=TKExclByTK(solution[i],solution[j]))
            {   printf("ERROR in FXSOLV: Problem with %d, %d!!\n",i,j);
                PrintTK(solution[i]);
                PrintTK(solution[j]);
            }
         }
   }

   /* *********** Write result to arrays ************ */
   if(fxdebug>0)
        printf("Solution Index List:\n");

   for(i=0;i<*nsolution;i++)
   {  if(solution[i]->id > 0)    /* These are already in the arrays */
      {  solutionindexlist[i]=solution[i]-TKList+1;
      }
      else                       /* The rest must be created in the arrays */
      {  if(*ntk>=*ntkmax)
         {  printf("ERROR in FXSOLV(): Not enough room to store the complete solution!!\n");
            printf("                   Increase PARAMETER MAXTKR in fxambi.car\n");
            exit(0);
         }
         else
         {  tkid[*ntk]  = solution[i]->id;
            tknte[*ntk] = solution[i]->nte;
            memcpy(tkr+*ntk**tkrsize,&solution[i]->tkr,*tkrsize*sizeof(INTEGER));

            for(j=0;j<solution[i]->nte;j++)
            {  tkteid[j+*ntk**ntemax]=solution[i]->te[j]->id;
               tktedet[j+*ntk**ntemax]=solution[i]->te[j]->detcode;
/*               if(solution[i]->te[j]->nexcl>0)
 *                 tkteexclid[j+*ntk**ntemax]=solution[i]->te[j]->excl[0]->id;
 *              else
 *                 tkteexclid[j+*ntk**ntemax]=0;
 */         }
            solutionindexlist[i]=++*ntk;
         }
      }
      if(fxdebug>0) printf("%d ",solutionindexlist[i]);
   }
   if(fxdebug>0)
   {  printf("\n");
      end_time=time(NULL);
      printf("Used (real) time: %ds\n",end_time-start_time);
   }
   *fxerror=errorbits;
   fflush(stdout);

}



/***************************************************
 *
 * Now follows code for running FXSOLV without DELANA
 *
 *************************************************** */
#ifdef STANDALONE
int
main(int argc,char *argv[])
{
   FILE *event;
   int j,error;
   int itk,ntk;
   int tkid[MAXTK];
   int tknte[MAXTK];
   int tkteid[MAXTK*MAXTE];
   int tktedet[MAXTK*MAXTE];
   int tkteexli[4*MAXTK*MAXTE];
   int tktelabl[MAXTK*MAXTE];
   int teexliids[MAXEXLI];
   int teexliused=0;
   struct TKR tkr[MAXTK];
   const struct TKR tkrnull={0};
   const int ntemax=MAXTE;
   const int ntkmax=MAXTK;
   const int tkrsize=sizeof(struct TKR)/sizeof(int);
   int solutionindexlist[MAXTK];
   int nsolution;
   int fxflags;
   int fxerror;
   time_t start_time,end_time;;

   if(argc<2||argc>3)
   {  printf("Wrong number of parameters\n");
      printf("Usage: %s Event.asc [Debug_level]\n",argv[0]);
      exit(EXIT_ERROR);
   }

   if(argc>2)
      fxflags=FLAG_DEBUGLEVEL&atoi(argv[2]);


   FORTRANNAME(fxsolvversion)();

   TEListUsed=0;
   event=fopen(argv[1],"r");
   for(itk=0;!feof(event)&&itk<MAXTK;itk++)
   {  error=fscanf(event,"%d %d",&tkid[itk],&tknte[itk]);
      if(error!=2) break;
      if(tknte[itk]>MAXTE)
      {  printf("FXSOLV: main(): Too many TEs in TK no %d. Found %d, max. allowed are %d.\n",
                 itk,tknte[itk],MAXTE);
         exit(EXIT_ERROR);
      }
      for(j=0;j<tknte[itk];j++)
      {  error=fscanf(event,"%d",&tkteid[itk*MAXTE+j]);
      }
      for(j=0;j<tknte[itk];j++)
      {  error=fscanf(event,"%d",&tktedet[itk*MAXTE+j]);
      }
      for(j=0;j<tknte[itk];j++)  /* At most 1 exclusion in the STANDALONE version */
      {  int tmpexcl;
         error=fscanf(event,"%d",&tmpexcl);
         if(tmpexcl)
         {  if(teexliused<MAXEXLI)
            {  teexliids[teexliused]=tmpexcl;
               tkteexli[itk*MAXTE*4+j*4]=1;
               tkteexli[itk*MAXTE*4+j*4+1]=teexliused;
               teexliused++;
            }
         }
         else
         {  tkteexli[itk*MAXTE*4+j*4]=0;
            tkteexli[itk*MAXTE*4+j*4+1]=0;
         }
      }
      for(j=0;j<tknte[itk];j++)  /* No links in the STANDALONE version */
      {  tkteexli[itk*MAXTE*4+j*4+2]=0;
         tkteexli[itk*MAXTE*4+j*4+3]=0;
      }
      for(j=0;j<tknte[itk];j++)  /* No labels in the STANDALONE version */
      {  tktelabl[itk*MAXTE+j]=0;
      }
      tkr[itk]=tkrnull;
   }

   ntk=itk;
   fxflags|=FLAG_1VDREMOVAL;
   /*fxflags|=FLAG_USEN2;*/

   start_time=time(NULL);

   FORTRANNAME(fxsolv)(&ntk, tkid, tknte,
        tkteid, tktedet, tkteexli, teexliids,
        tktelabl,(int *)tkr,
        &ntemax, &ntkmax,  &tkrsize,
        solutionindexlist, &nsolution, &fxflags, &fxerror);


   end_time=time(NULL);
   printf("Used (real) time: %ds\n",end_time-start_time);
   printf("IFail %d\nIWarn %d\n",fxerror>>16,fxerror&0x00FF);
   printf("TEListUsed %d\n",TEListUsed);
   printf("TEExLiUsed %d\n",teexliused);
   return 0;
}
#endif
