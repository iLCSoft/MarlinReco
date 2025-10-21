#define NOVARS 16
#define VARIABLES                                                                                                      \
  {"EtoN_ecal", "EtoN_hcal", "ex",    "dmean", "L1", "EecalToEtot", "Eecal", "Ehcal",                                  \
   "Edmean",    "Necal",     "Nhcal", "L2",    "L3", "EL1",         "EL2",   "EL3"}
// #define  RANGES {100, 0.1, 0.0, 100, 0.1, 0.0, 100, 0.001, 0.0, 100, 100., 0., 40, 40., 0., 100, 1.0, 0.0, 100, 5.,
// 0., 100, 5., 0., 100, 20., 3., 100, 100., 20., 100, 100., 20., 40, 40., 0., 40, 40., 0., 30, 30., 0., 30, 30., 0.,
// 30, 30., 0.}

/********************************************************************/
/*                                                                  */
/*  charged particles                                               */
/*                                                                  */
/* CPDFNAME    : name of pdf                                        */
/* CNOCATS     : number of categories                               */
/*                  (default 3 : electron, muon, pion)              */
/* CCATS       : names of categories                                */
/* CNOVARS     : number of availabe variables                       */
/* CVARIABLES  : names of these variables                           */
/*               which variables are taken must be specified in     */
/*               the code the CreatePDFs.cc file                    */
/* availabe variables are (in order, se above):                     */
/*              E/N in Ecal, E/N in hcal (E=energy, N=#hits)        */
/*              excentricity of cluster                             */
/*              mean distance of cluster hits to track              */
/*              last layer of hcal with hit                         */
/*              E in ecal / total cluster energy                    */
/*              E in ecal, E in hcal                                */
/*              energy weighted mean distance of cluster hits to track */
/*              N in ecal, N in hcal                                */
/*              next to last hcal layer with hit                    */
/*              next to next to last hcal layer with hit            */
/*              1st ecal layer with hit                             */
/*              2nd ecal layer with hit                             */
/*              3rd ecal layer with hit                             */
/* CRANGES     : Ranges of the variables                            */
/*               #bins for var 1, max for var1, min for var1,       */
/*               #bins for var 2, ...                               */
/*  !! # bins strongly depends on the number of events per category */
/* CNOHISTS    : number of histograms used                          */
/*               how histograms are implemented -> CreatePDFs.cc    */
/* CHISTDIM    : dimension of the histograms !<= CVARIABLES !!!     */
/*                                                                  */
/********************************************************************/

#define CPDFNAME "charged"
#define CNOCATS 3
#define CCATS {"electron", "muon", "pion"}
#define CNOVARS 6
#define CVARIABLES {"EtoN_ecal", "EtoN_hcal", "ex", "dmean", "L1", "EecalToEtot"}
#define CRANGES {100, 0.1, 0.0, 100, 0.1, 0.0, 100, 0.001, 0.0, 100, 100., 0., 40, 40., 0., 100, 1.0, 0.0}
#define CNOHISTS 6
#define CHISTDIM 1

// neutral particles
#define NPDFNAME "neutral"
#define NNOCATS 2
#define NCATS {"photon", "kaon/neutron"}
#define NNOVARS 7
#define NVARIABLES {"EtoN_ecal", "EtoN_hcal", "ex", "dmean", "L1", "EecalToEtot", "EL1"}
#define NRANGES {100, 0.1, 0.0, 100, 0.1, 0.0, 100, 0.001, 0.0, 100, 100., 0., 40, 40., 0., 100, 1.0, 0.0, 30, 30., 0.}
#define NNOHISTS 7
#define NHISTDIM 1

// Later maybe with real interface
