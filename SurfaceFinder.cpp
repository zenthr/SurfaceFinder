/*ingw*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <float.h>
#include <vector>
#include <time.h>

#include <stdio.h>

# define PIE 3.141592653589793
# define GEVFM 0.1975


#define NumStep 50 // GridPoints in x,y

# define TAUSTART 0.1    // in fm; initial value for tau, cannot be 0 
# define TAUSTEPFM 0.2   // in fm

# define TAUCOUNTMAX  96   //(TAUSTOP-TAUSTART)/TAUSTEPFM; Number of tau steps

# define XMAX 10
# define XSTEP 0.2

# define YMAX 10
# define YSTEP 0.2

# define ETAMAX 2.5
# define ZSTEP 0.2


/*******************************************************************/
/*******************************************************************/
/*   3+1D hypersurface finder for JETSCAPE                         */
/*   v 1.0.0 (Mar 12 2019)                                         */
/*   Author: Steven Rose       ( zenth @ physics.tamu.edu )        */
/*   With support of:                                              */
/*   Dr. Rainer Fries          ( rjfries @ comp.tamu.edu )         */
/*   Dr. Michael Kordell       ( M.Kordell @ wayne.edu )           */
/*   at Texas A&M University Cyclotron Institute                   */
/*******************************************************************/
/*   Current parameters of are located in:                         */
/*   maindef.h                                                     */
/*   There, edit step sizes in tau, x, y, eta                      */
/*   As well as intial tau value (should be > 0 )                  */
/*******************************************************************/
/*   Reads in ????????????????                                     */                   
/*   Expecting format:                                             */
/*   ???????????????????????????????????????????????               */
/*   for each grid point from hydro simulation                     */
/*******************************************************************/
/*   Writes out to ~/surf.dat                                      */
/*   Output format:                                                */
/*   Position in t, x, y, z (averaged position of verteces)        */
/*   dsigma t, x, y, z                                             */
/*   Vx, Vy, Vz                                                    */
/*******************************************************************/

using namespace std;

    /* Global Counters */
    static int GridCount;
    double GridPlot[4000000][8];

    /* Grid-level parameters */
    int hyperbolic;                                // 0 for Cartesian Input, 1 for Hyperbolic
    /*NOT UPDATED*/double SStep = 0.2;
    /*NOT UPDATED*/int SCoarseGrain = 1;                           // Coarse graining in space coords
    /*NOT UPDATED*/int TCoarseGrain = 1;                           // Coarse graining in time coord

    double Crit;                                   // Critical Value
    double Area[9000000][8];                       // List of area vectors and grid coordinates already computed
    double U[4];                                   //four-velocity vector
    int Link[32][2];                               // Which vertecies bound which links
    int GridPos[4];                                // Retains where on the grid is being evaluated (lowest corner's (t,x,y,z))
    typedef double (*ARR) [NumStep][NumStep][int ((2*ETAMAX/ZSTEP)+1+0.5)][1]; //energydensity/temp; Yes that is an array of size 1. Deal With it. 
    ARR Grid;
    typedef double (*ARRV) [NumStep][NumStep][3];//vx,vy,veta
    ARRV V;
    typedef double (*ARR2) [NumStep][NumStep][6]; // pi's (five)
    ARR2 PiGrd;
    typedef double (*ARR3) [NumStep][NumStep][2]; // pi's (five)
    ARR3 EnPr;
    
    /* Cell-level parameters */
    double Vert[16][4];                             // For Each Vertex [Energy, vx, vy, veta]
    bool BVert[16];                                 // Contains High/Low (True/False) information about each vertex
    double v_LE[32][4];                            // List of low energy direction on each link
    double V_LE[4];                                // Total low energy direction for the cell
    double Point[32][4];                           // Lists t,x,y coordinates of intersection points
    double Patch[1][4];                         // List of area patches in the current cell
    double SpanV[3][4];                            // Two vectors that span the current triangle
    int Tetra[3000][6];                            // List of triangles in cube containing points (numeric ordering), opposing point, and "degeneracy"
    double TetraVec[3000][8];                        // Contains vector information (Area vector and Center) of triangles. Updated alongside Tri!
    /*NOT UPDATED*/int NPoint;                                    // Number of points found
    int NTetra;                                      // Number of triangular patches constructed
    /*NOT UPDATED*/int NArea;                                     // Number of viable patches calculated
    static int flatness;

    int GPoints = 50;
    double GWeight[50];   //Gaussian Weights
    double GAbs[50];      //Gaussian Abscissa
    ifstream TFile;
    ifstream VXFile;
    ifstream VYFile;
//    /*NOT UPDATED*/ofstream SpecFile[1001];
    /*NOT UPDATED*/static int GridT;
    /*NOT UPDATED*/double CritTemp;
    /*NOT UPDATED*/static double CellList[1000][2];
    /*NOT UPDATED*/static double ShellEnergy;
    /*NOT UPDATED*/static double ShellPX;
    /*NOT UPDATED*/static double ShellPY;
    /*NOT UPDATED*/static double ShellEMT[3];
    /*NOT UPDATED*/static double GroundEnergy;
    /*NOT UPDATED*/static int brake;
    
    /*NOT UPDATED*/double rem1; double rem2; double rem3; double rem4;
    /*NOT UPDATED*/static double TensorEnergy;
    /*NOT UPDATED*/static double TensorPX;
    /*NOT UPDATED*/static double TensorPY;
    
    /*NOT UPDATED*/static double TensorKinetic;
    /*NOT UPDATED*/static double ShellKinetic;
    
    /*NOT UPDATED*/static int AbsorbtionEvents;
    /*NOT UPDATED*/int rollmax = 100;
    
    /* Outside Functions */
    /*NOT UPDATED*/void BoseTest(double Temp);
    
    /* Levi Civita */
    static double LC[4][4][4][4];
    
    int saveme;

    //TestVars
    int FailSafe;
    static int ITry;
    static int IFale;
    static int ITry2;
    static int IFale2;
    static int spacey;
    static double holddot;
    static int deltcheck;
    static int deltchange;


    /*Array Definitions*/
    
    void DefGrid(ARR* buf)
{
    *buf = new double[2][NumStep][NumStep][int ((2*ETAMAX/ZSTEP)+1+0.5)][1] ;
    
    for(int i1=0; i1<2; i1++)
    {
    for(int i2=0; i2<NumStep; i2++)
    {
    for(int i3=0; i3<NumStep; i3++)
    {
    for(int i4=0; i4<(2*ETAMAX/ZSTEP)+1+0.5; i4++)
    {
   for(int i5=0; i5<1; i5++)
    {
        (*buf)[i1][i2][i3][i4][i5]=0;
    }
    }
    }
    }
    }
}

    void DefV(ARRV* buf)
{
    *buf = new double[2][NumStep][NumStep][3] ;
    
    for(int i1=0; i1<2; i1++)
    {
    for(int i2=0; i2<NumStep; i2++)
    {
    for(int i3=0; i3<NumStep; i3++)
    {
   for(int i5=0; i5<3; i5++)
    {
        (*buf)[i1][i2][i3][i5]=0;
    }
    }
    }
    }
}

    void DefPi(ARR2* buf)
{
    *buf = new double[2][NumStep][NumStep][6] ;
    
    for(int i1=0; i1<2; i1++)
    {
    for(int i2=0; i2<NumStep; i2++)
    {
    for(int i3=0; i3<NumStep; i3++)
    {
   for(int i5=0; i5<5; i5++)
    {
        (*buf)[i1][i2][i3][i5]=0;
    }
    }
    }
    }
}

    void DefEnPr(ARR3* buf)
{
    *buf = new double[2][NumStep][NumStep][2] ;
    
    for(int i1=0; i1<2; i1++)
    {
    for(int i2=0; i2<NumStep; i2++)
    {
    for(int i3=0; i3<NumStep; i3++)
    {
   for(int i5=0; i5<2; i5++)
    {
        (*buf)[i1][i2][i3][i5]=0;
    }
    }
    }
    }
}

    void RelGrid(ARR* buf)
    {
        delete[](*buf);
    }

    /* Cell Functions */

    /* Get v_LE in t direction */
    void TLE (int C1, int C2, double tau, double x, double y, double eta)
    {
       if(hyperbolic==1)
       {
          double tau = Point[NPoint][0];
          double eta = Point[NPoint][3];
          Point[NPoint][0] = tau*cosh(eta);
          Point[NPoint][3] = tau*sinh(eta);
       }
       
       v_LE[NPoint][0] = tau*cosh(eta) - Point[NPoint][0];
       v_LE[NPoint][1] = x - Point[NPoint][1];
       v_LE[NPoint][2] = y - Point[NPoint][2];
       v_LE[NPoint][3] = tau*sinh(eta) - Point[NPoint][3];
    }
    
    /* Get v_LE in x direction */
    void XLE (int C1, int C2, double tau, double x, double y, double eta)
    {
       if(hyperbolic==1)
       {
          double tau = Point[NPoint][0];
          double eta = Point[NPoint][3];
          Point[NPoint][0] = tau*cosh(eta);
          Point[NPoint][3] = tau*sinh(eta);
       }
                          
       v_LE[NPoint][0] = tau*cosh(eta) - Point[NPoint][0];
       v_LE[NPoint][1] = x - Point[NPoint][1];
       v_LE[NPoint][2] = y - Point[NPoint][2];
       v_LE[NPoint][3] = tau*sinh(eta) - Point[NPoint][3];
    }
    
    /* Get v_LE in y direction */
    void YLE (int C1, int C2, double tau, double x, double y, double eta)
    {
       if(hyperbolic==1)
       {
          double tau = Point[NPoint][0];
          double eta = Point[NPoint][3];
          Point[NPoint][0] = tau*cosh(eta);
          Point[NPoint][3] = tau*sinh(eta);
       }
       
       v_LE[NPoint][0] = tau*cosh(eta) - Point[NPoint][0];
       v_LE[NPoint][1] = x - Point[NPoint][1];
       v_LE[NPoint][2] = y - Point[NPoint][2];
       v_LE[NPoint][3] = tau*sinh(eta) - Point[NPoint][3];
    }

    /* Get v_LE in y direction */
    void ZLE (int C1, int C2, double tau, double x, double y, double eta)
    {
       if(hyperbolic==1)
       {
          double tau = Point[NPoint][0];
          double eta = Point[NPoint][3];
          Point[NPoint][0] = tau*cosh(eta);
          Point[NPoint][3] = tau*sinh(eta);
       }
       
       v_LE[NPoint][0] = tau*cosh(eta) - Point[NPoint][0];
       v_LE[NPoint][1] = x - Point[NPoint][1];
       v_LE[NPoint][2] = y - Point[NPoint][2];
       v_LE[NPoint][3] = tau*sinh(eta) - Point[NPoint][3];
    }
    
    /* Analyze cube, find intersections and low energy direction */
    /* May not need absolute position of points- left in as example if needed later */

/******************************/
/*********NEEDS UPDATE*********/
/******************************/

    void CubeAnalyzer (double t, double x, double y, double eta)
    {
         int i;
         int j;
         
         for(i=0; i<32; i++)
         {
              int C1 = Link[i][0];
              int C2 = Link[i][1];
              if (BVert[C1] ^ BVert[C2])
              {
                   /* 0-11 lower t cube */
                   /* If on link 0... */
                   if ( i==0 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = eta;
                        
                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t, x, y, eta);
                        }
                        else
                        {
                           XLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta); 
                        }
                   }
                   
                   /* If on link 1... */                   
                   else if ( i==1 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]))*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t, x, y, eta);
                        }
                        else
                        {
                           YLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }

                   /* If on link 2... */   
                   else if ( i==2 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta);
                        }
                        else
                        {
                           YLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta); 
                        }
                        
                   }

                   /* If on link 3... */   
                   else if ( i==3 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        
                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           XLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }

                   /* If on link 4... */   
                   else if ( i==4 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t, x, y, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t, x, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 5... */   
                   else if ( i==5 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*ZSTEP*SCoarseGrain + eta;

                        
                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 6... */   
                   else if ( i==6 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 7... */   
                   else if ( i==7 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 8... */   
                   else if ( i==8 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t, x, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           XLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 9... */   
                   else if ( i==9 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;
                        
                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t, x, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           YLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 10... */   
                   else if ( i==10 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;
                        
                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           YLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 11... */   
                   else if ( i==11 )
                   {
                        Point[NPoint][0] = t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ))*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;
                        
                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           XLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /************************/
                   /* 12-19 temporal edges */
                   /************************/

                   /* If on link 12... */                   
                   else if ( i==12 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x, y, eta);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta); 
                        }
                   }

                   /* If on link 13... */                   
                   else if ( i==13 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta); 
                        }
                   }

                   /* If on link 14... */                   
                   else if ( i==14 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }

                   /* If on link 15... */                   
                   else if ( i==15 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }

                   /* If on link 16... */                   
                   else if ( i==16 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 17... */                   
                   else if ( i==17 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 18... */                   
                   else if ( i==18 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /* If on link 19... */                   
                   else if ( i==19 )
                   {
                        Point[NPoint][0] = ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) )*TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           TLE (C1, C2, t, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           TLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }

                   /************************/
                   /* 20-31 "upper t" cube */
                   /************************/

                   /* If on link 20... */                   
                   else if ( i==20 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta);
                        }
                        else
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta); 
                        }
                   }
                   
                   /* If on link 21... */                   
                   else if ( i==21 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta);
                        }
                        else
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }
                   
                   /* If on link 22... */                   
                   else if ( i==22 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] =XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta);
                        }
                        else
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }
                   
                   /* If on link 23... */                   
                   else if ( i==23 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = eta;

                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta); 
                        }
                   }
                   
                   /* If on link 24... */                   
                   else if ( i==24 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 25... */                   
                   else if ( i==25 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 26... */                   
                   else if ( i==26 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 27... */                   
                   else if ( i==27 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta);
                        }
                        else
                        {
                           ZLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 28... */                   
                   else if ( i==28 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 29... */                   
                   else if ( i==29 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = x;
                        Point[NPoint][2] =( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 30... */                   
                   else if ( i==30 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           YLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* If on link 31... */                   
                   else if ( i==31 )
                   {
                        Point[NPoint][0] = TAUSTEPFM*TCoarseGrain + t;
                        Point[NPoint][1] = ( ( (Crit - Vert[C1][0])/(Vert[C2][0]-Vert[C1][0]) ) )*XSTEP*SCoarseGrain + x;
                        Point[NPoint][2] = YSTEP*SCoarseGrain + y;
                        Point[NPoint][3] = ZSTEP*SCoarseGrain + eta;

                        if(Crit > Vert[C1][0])
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain);
                        }
                        else
                        {
                           XLE (C1, C2, t + TAUSTEPFM*TCoarseGrain, x + XSTEP*SCoarseGrain, y + YSTEP*SCoarseGrain, eta + ZSTEP*SCoarseGrain); 
                        }
                   }
                   
                   /* Potential errors? */   
                   else
                   {
                        cout << "The cube analyzer has been confused!";
                   }
                   
                   bool Duplicity = false;
                   
                   for(j=0;j<NPoint;j++)
                   {
                       if( (Point[j][0] == Point[NPoint][0]) && (Point[j][1] == Point[NPoint][1]) && (Point[j][2] == Point[NPoint][2]) && (Point[j][3] == Point[NPoint][3]) )
                       {
                           Duplicity = true;
                       }
                   }
                   
                   if(!Duplicity)
                   {
                       NPoint++;
                   }
              }
         }
         
         V_LE[0] = v_LE[0][0];
         V_LE[1] = v_LE[0][1];
         V_LE[2] = v_LE[0][2];
         V_LE[3] = v_LE[0][3];
         
         for(i=1; i<NPoint; i++)
         {
              V_LE[0] = V_LE[0] + v_LE[i][0];
              V_LE[1] = V_LE[1] + v_LE[i][1];
              V_LE[2] = V_LE[2] + v_LE[i][2];
              V_LE[3] = V_LE[3] + v_LE[i][3];
         }
    }


    /* Test two products with dot product */

/******************************/
/*********NEEDS UPDATE*********/
/******************************/

    bool Dot (double at, double ax, double ay, double az, double bt, double bx, double by, double bz, int chcksum)
    {
         double a = sqrt(at*at + ax*ax + ay*ay + az*az);
         double b = sqrt(bt*bt + bx*bx + by*by + bz*bz);
         at = at/a; ax = ax/a; ay = ay/a; az = az/a;
         bt = bt/b; bx = bx/b; by = by/b; bz = bz/b;
         double dot = (at*bt + ax*bx + ay*by + az*bz); // Uses Euclidean metric.
         

         if (dot > (DBL_EPSILON*1000.))
         {
              return(true);
         }
                        
         else
         {

             return(false);
         }
    }
    
    /* Triangle Checker- manages appending and "degeneracy" */
    
/******************************/
/*********NEEDS UPDATE*********/
/******************************/

    void AppTetra (int k, int a, int b, int c, int d, int i, int j)  // k is some tracker. Tetra abci is opposit point d NOTE: j is tetra connected to needed for testing only
    {

         if(k<NTetra)
         {
             if (Tetra[k][3] == i && Tetra[k][0] ==  a && Tetra[k][1] == b && Tetra[k][2] == c) // Does Kth tetra match "new" tetra?
             {
                 /* If this is triggered, we need to force in a degeneracy check as well */
                 if (!(Tetra[k][5] == 0))
                 {
//                    cout << "\nDEGENERACY NEEDS TO BE CHECKED!\n";
//                    cout << "WAS GIVEN DEGENERACY OF: " << Tetra[k][5] << " FOR TETRAHEDRON: " << k <<"\n";
                 }
                
                 Tetra[k][5] = 1;   // Kth tetra degenerates, search ends in dismal failure.
             }
             
             else
             { AppTetra (k+1, a, b, c, d, i, j); }   // Increase k and try again
         }
         
         else if (b>a && c>b)
         {

double TetrajSize = sqrt(TetraVec[j][0]*TetraVec[j][0] + TetraVec[j][1]*TetraVec[j][1] + TetraVec[j][2]*TetraVec[j][2] + TetraVec[j][3]*TetraVec[j][3]);
double DistiLength = sqrt( (TetraVec[j][4]-Point[i][0])*(TetraVec[j][4]-Point[i][0]) + (TetraVec[j][5]-Point[i][1])*(TetraVec[j][5]-Point[i][1]) +
                     (TetraVec[j][6]-Point[i][2])*(TetraVec[j][6]-Point[i][2]) + (TetraVec[j][7]-Point[i][3])*(TetraVec[j][7]-Point[i][3]) );

double prod = (TetraVec[j][0]*(TetraVec[j][4]-Point[i][0]) + TetraVec[j][1]*(TetraVec[j][5]-Point[i][1]) + 
               TetraVec[j][2]*(TetraVec[j][6]-Point[i][2]) + TetraVec[j][3]*(TetraVec[j][7]-Point[i][3]))/(TetrajSize*DistiLength);

if( fabs( prod ) < DBL_EPSILON)
{
   int RandCoord = rand() %4;
   Point[i][0]= Point[i][0] + 10e-6;
   
   ITry2++;
   
if( fabs( (TetraVec[j][0]*(TetraVec[j][4]-Point[i][0]) + TetraVec[j][1]*(TetraVec[j][5]-Point[i][1]) + TetraVec[j][2]*(TetraVec[j][6]-Point[i][2]) + TetraVec[j][3]*(TetraVec[j][7]-Point[i][3]))/(TetrajSize*DistiLength) ) > DBL_EPSILON)
   {
      IFale2++;
   }
}


            Tetra[NTetra][0] = a;
            Tetra[NTetra][1] = b;
            Tetra[NTetra][2] = c; 
            Tetra[NTetra][3] = i; //Tetra ABCI opposite point D
            Tetra[NTetra][4] = d;
            Tetra[NTetra][5] = 0;
            
            for(int l=0; l<4; l++)
            {
                SpanV[0][l] = Point[a][l] - Point[b][l];
                SpanV[1][l] = Point[a][l] - Point[c][l];
                SpanV[2][l] = Point[a][l] - Point[i][l];
            }

            TetraVec[NTetra][0] = 0; TetraVec[NTetra][1] = 0;
            TetraVec[NTetra][2] = 0; TetraVec[NTetra][3] = 0;
            for(int mu=0; mu<4; mu++) // general cross product of vectors
            {
               for(int nu=0; nu<4; nu++)
               {
                  for(int rho=0; rho<4; rho++)
                  {
                     for(int sig=0; sig<4; sig++)
                     {
                        TetraVec[NTetra][sig] = TetraVec[NTetra][sig] + (1./6.)*( LC[mu][nu][rho][sig]*SpanV[0][mu]*SpanV[1][nu]*SpanV[2][rho]  );
                     }
                  }
               }
            }
            
            TetraVec[NTetra][4] = (1./4.)*(Point[a][0] + Point[b][0] + Point[c][0] + Point[i][0]);
            TetraVec[NTetra][5] = (1./4.)*(Point[a][1] + Point[b][1] + Point[c][1] + Point[i][1]);
            TetraVec[NTetra][6] = (1./4.)*(Point[a][2] + Point[b][2] + Point[c][2] + Point[i][2]);
            TetraVec[NTetra][7] = (1./4.)*(Point[a][3] + Point[b][3] + Point[c][3] + Point[i][3]);
            
            double TArea = sqrt( TetraVec[NTetra][0]*TetraVec[NTetra][0] + TetraVec[NTetra][1]*TetraVec[NTetra][1] + TetraVec[NTetra][2]*TetraVec[NTetra][2] + TetraVec[NTetra][3]*TetraVec[NTetra][3] );
            double CellL = sqrt( (GridPos[0]-49)*(GridPos[0]-49) + (GridPos[1]-49)*(GridPos[1]-49) + (GridPos[2]-49)*(GridPos[2]-49) + (GridPos[3]-49)*(GridPos[3]-49) );
                    
            if(!Dot(TetraVec[NTetra][0], TetraVec[NTetra][1], TetraVec[NTetra][2], TetraVec[NTetra][3], (TetraVec[NTetra][4]-Point[d][0]), ((TetraVec[NTetra][5]-Point[d][1])), ((TetraVec[NTetra][6]-Point[d][2])), ((TetraVec[NTetra][7]-Point[d][3])), 10 ))
            {
                TetraVec[NTetra][0] = - TetraVec[NTetra][0];
                TetraVec[NTetra][1] = - TetraVec[NTetra][1];
                TetraVec[NTetra][2] = - TetraVec[NTetra][2];
                TetraVec[NTetra][3] = - TetraVec[NTetra][3];
            }

            NTetra++;    // If not found, add new nondeg tetra!

         }
         
         else
         {
             cout << "\nERROR! AppTetra given incongruous arguments!/n";
         }
    }

    /* Tetrahedron Finder */

/******************************/
/*********NEEDS UPDATE*********/
/******************************/

    void FindTetras ()
    {
//         cout << "Called it." << endl;
        int i; int j; int k;
        /*ro    ????*/
        //Initial pentachoron:
        Tetra[0][0] = 0; // Point of surface vol
        Tetra[0][1] = 1; // Point of surface vol
        Tetra[0][2] = 2; // Point of surface vol
        Tetra[0][3] = 3; // Point of surface vol
        Tetra[0][4] = 4; // Point opposite surface element (set makes pentachoron)
        Tetra[0][5] = 0; // Is initially exterior
        
        Tetra[1][0] = 0;
        Tetra[1][1] = 1;
        Tetra[1][2] = 2;
        Tetra[1][3] = 4;
        Tetra[1][4] = 3;
        Tetra[1][5] = 0;
        
        Tetra[2][0] = 0;
        Tetra[2][1] = 1;
        Tetra[2][2] = 3;
        Tetra[2][3] = 4;
        Tetra[2][4] = 2;
        Tetra[2][5] = 0;
            
        Tetra[3][0] = 0;
        Tetra[3][1] = 2;
        Tetra[3][2] = 3;
        Tetra[3][3] = 4;
        Tetra[3][4] = 1;
        Tetra[3][5] = 0;

        Tetra[4][0] = 1;
        Tetra[4][1] = 2;
        Tetra[4][2] = 3;
        Tetra[4][3] = 4;
        Tetra[4][4] = 0;
        Tetra[4][5] = 0;
            
        NTetra = 5;
        flatness = 0;

        for(i=0; i<NTetra; i++)
        {          
            for(j=0; j<4; j++)
            {
                SpanV[0][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][1] ][j];
                SpanV[1][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][2] ][j];
                SpanV[2][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][3] ][j];
            }
            

            TetraVec[i][0] = 0; TetraVec[i][1] = 0;
            TetraVec[i][2] = 0; TetraVec[i][3] = 0;
            for(int mu=0; mu<4; mu++) // general cross product of vectors
            {
               for(int nu=0; nu<4; nu++)
               {
                  for(int rho=0; rho<4; rho++)
                  {
                     for(int sig=0; sig<4; sig++)
                     {
                        TetraVec[i][sig] = TetraVec[i][sig] + (1./6.)*( LC[mu][nu][rho][sig]*SpanV[0][mu]*SpanV[1][nu]*SpanV[2][rho]  );
                     }
                  }
               }
            }

            TetraVec[i][4] = (1./4.)*( Point[ Tetra[i][0] ][0] + Point[ Tetra[i][1] ][0] + Point[ Tetra[i][2] ][0] + Point[ Tetra[i][3] ][0] );
            TetraVec[i][5] = (1./4.)*( Point[ Tetra[i][0] ][1] + Point[ Tetra[i][1] ][1] + Point[ Tetra[i][2] ][1] + Point[ Tetra[i][3] ][1] );
            TetraVec[i][6] = (1./4.)*( Point[ Tetra[i][0] ][2] + Point[ Tetra[i][1] ][2] + Point[ Tetra[i][2] ][2] + Point[ Tetra[i][3] ][2] );
            TetraVec[i][7] = (1./4.)*( Point[ Tetra[i][0] ][3] + Point[ Tetra[i][1] ][3] + Point[ Tetra[i][2] ][3] + Point[ Tetra[i][3] ][3] );

            Dot(TetraVec[i][0], TetraVec[i][1], TetraVec[i][2], TetraVec[i][3], (TetraVec[i][4]-Point[ Tetra[i][4] ][0]), ((TetraVec[i][5]-Point[ Tetra[i][4] ][1])), ((TetraVec[i][6]-Point[ Tetra[i][4] ][2])), ((TetraVec[i][7]-Point[ Tetra[i][4] ][3])), 10 );


            
            if(!Dot(TetraVec[i][0], TetraVec[i][1], TetraVec[i][2], TetraVec[i][3], (TetraVec[i][4]-Point[ Tetra[i][4] ][0]), ((TetraVec[i][5]-Point[ Tetra[i][4] ][1])), ((TetraVec[i][6]-Point[ Tetra[i][4] ][2])), ((TetraVec[i][7]-Point[ Tetra[i][4] ][3])), 1 ))
            {
                TetraVec[i][0] = - TetraVec[i][0];
                TetraVec[i][1] = - TetraVec[i][1];
                TetraVec[i][2] = - TetraVec[i][2];
                TetraVec[i][3] = - TetraVec[i][3];
            }            
            
            if(TetraVec[i][0] == 0 && TetraVec[i][1] == 0 && TetraVec[i][2] == 0 && TetraVec[i][3] == 0)
            {
               Point[ Tetra[i][0] ][0] = Point[ Tetra[i][0] ][0] + 10e-6*TAUSTEPFM;
               i--;
            }

        }

//Spacial/Planar Checks
        if( flatness > 0)
        {
           flatness =0;
           int RandPoint = rand() %5;
           int RandCoord = rand() %4;
           if(RandCoord ==0)
           {
              Point[RandPoint][RandCoord] = Point[RandPoint][RandCoord] + TAUSTEPFM*10e-6;
           }
           else if(RandCoord ==1)
           {
              Point[RandPoint][RandCoord] = Point[RandPoint][RandCoord] + XSTEP*10e-6;
           }
           else if(RandCoord ==2)
           {
              Point[RandPoint][RandCoord] = Point[RandPoint][RandCoord] + YSTEP*10e-6;
           }
           else if(RandCoord ==3)
           {
              Point[RandPoint][RandCoord] = Point[RandPoint][RandCoord] + XSTEP*10e-6; //uses XSTEP because 
           }
//           Point[0][0] = Point[0][0] + 10e-6;
           
           for(i=0; i<NTetra; i++)
              {

                 for(j=0; j<4; j++)
                 {
                    SpanV[0][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][1] ][j];
                    SpanV[1][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][2] ][j];
                    SpanV[2][j] = Point[ Tetra[i][0] ][j] - Point[ Tetra[i][3] ][j];
                 }
                
                TetraVec[i][0] = 0; TetraVec[i][1] = 0;
                TetraVec[i][2] = 0; TetraVec[i][3] = 0;
                for(int mu=0; mu<4; mu++) // general cross product of vectors
                {
                   for(int nu=0; nu<4; nu++)
                   {
                      for(int rho=0; rho<4; rho++)
                      {
                         for(int sig=0; sig<4; sig++)
                         {
                            TetraVec[i][sig] = TetraVec[i][sig] + (1./6.)*( LC[mu][nu][rho][sig]*SpanV[0][mu]*SpanV[1][nu]*SpanV[2][rho]  );
                         }
                      }
                   }
                }
            
                TetraVec[i][4] = (1./4.)*( Point[ Tetra[i][0] ][0] + Point[ Tetra[i][1] ][0] + Point[ Tetra[i][2] ][0] + Point[ Tetra[i][3] ][0] );
                TetraVec[i][5] = (1./4.)*( Point[ Tetra[i][0] ][1] + Point[ Tetra[i][1] ][1] + Point[ Tetra[i][2] ][1] + Point[ Tetra[i][3] ][1] );
                TetraVec[i][6] = (1./4.)*( Point[ Tetra[i][0] ][2] + Point[ Tetra[i][1] ][2] + Point[ Tetra[i][2] ][2] + Point[ Tetra[i][3] ][2] );
                TetraVec[i][7] = (1./4.)*( Point[ Tetra[i][0] ][3] + Point[ Tetra[i][1] ][3] + Point[ Tetra[i][2] ][3] + Point[ Tetra[i][3] ][3] );

                //If N dot OUT < 0
                
                /*TEST*/Dot(TetraVec[i][0], TetraVec[i][1], TetraVec[i][2], TetraVec[i][3], (TetraVec[i][4]-Point[ Tetra[i][4] ][0]), ((TetraVec[i][5]-Point[ Tetra[i][4] ][1])), ((TetraVec[i][6]-Point[ Tetra[i][4] ][2])), ((TetraVec[i][7]-Point[ Tetra[i][4] ][3])), 10 );


                if(!Dot(TetraVec[i][0], TetraVec[i][1], TetraVec[i][2], TetraVec[i][3], (TetraVec[i][4]-Point[ Tetra[i][4] ][0]), ((TetraVec[i][5]-Point[ Tetra[i][4] ][1])), ((TetraVec[i][6]-Point[ Tetra[i][4] ][2])), ((TetraVec[i][7]-Point[ Tetra[i][4] ][3])), 1 ))
                {
                   TetraVec[i][0] = - TetraVec[i][0];
                   TetraVec[i][1] = - TetraVec[i][1];
                   TetraVec[i][2] = - TetraVec[i][2];
                   TetraVec[i][3] = - TetraVec[i][3];
                }
              }           
              ITry++;
   
              if( flatness > 0)
              {
                 IFale++;
              }
}
                  

        for (i=5; i<NPoint; i++) //For each additional point...
        {      
            for (j=(NTetra-1); j>-1; j--) //Try to connect to every previous exterior teatrahedron
            {
                if( Dot(TetraVec[j][0], TetraVec[j][1], TetraVec[j][2], TetraVec[j][3], (Point[i][0]- TetraVec[j][4]), (Point[i][1]- TetraVec[j][5]), (Point[i][2]- TetraVec[j][6]), (Point[i][3]- TetraVec[j][7]),1010) && Tetra[j][5] == 0 ) // NonDeg and outside
                {
                    Tetra[j][5] = 1;    // This Tetrahedron immeditately degenerates.
                    
                    AppTetra (0, Tetra[j][0], Tetra[j][1], Tetra[j][2], Tetra[j][3], i, j); // Because we haven't nested enough complex statements to do SCIENCE here yet!
                    AppTetra (0, Tetra[j][0], Tetra[j][1], Tetra[j][3], Tetra[j][2], i, j);
                    AppTetra (0, Tetra[j][0], Tetra[j][2], Tetra[j][3], Tetra[j][1], i, j);
                    AppTetra (0, Tetra[j][1], Tetra[j][2], Tetra[j][3], Tetra[j][0], i, j);
                }
            }
        }
    }

    
    /* Read data */

    void ReadGen (int T)
    {
        int i; int j; int k;
        
        double tau = T*TAUSTEPFM+TAUSTART;

        for (i=0; i<NumStep; i++)
        {
            for (j=0; j<NumStep; j++)
            {
                for (k=0; k<((2*ETAMAX/ZSTEP)+1+0.5); k++)
                {
                    Grid[0][i][j][k][0] = Grid[1][i][j][k][0];
                }
                
                V[0][i][j][0] = V[1][i][j][0];
                V[0][i][j][1] = V[1][i][j][1];
                V[0][i][j][2] = V[1][i][j][2];
                EnPr[0][i][j][0] =EnPr[1][i][j][0];
                EnPr[0][i][j][1] =EnPr[1][i][j][1];
                PiGrd[0][i][j][0] = PiGrd[1][i][j][0];
                PiGrd[0][i][j][1] = PiGrd[1][i][j][1];
                PiGrd[0][i][j][2] = PiGrd[1][i][j][2];
                PiGrd[0][i][j][3] = PiGrd[1][i][j][3];
                PiGrd[0][i][j][4] = PiGrd[1][i][j][4];
                PiGrd[0][i][j][5] = PiGrd[1][i][j][5];

            }
        }

/*Sid Reading Direction*/

    int len = snprintf(0,0,"input/Hydro_Sid_Def/tau%2.2ffm.bin",tau);
    char *buf = new char[len+3];
    sprintf(buf,"input/Hydro_Sid_Def/tau%2.2ffm.bin",tau);

    ifstream init;
    init.open(buf,ios::in | ios::binary);

	if(!init.is_open())
	   {
		printf ("\nCould not open file\n");
		cout << buf << endl;
       }


	if(!init.good())
		printf ("\n file stream not in good shape \n");
		

	for( i=0; i< NumStep; i++)
	{
	for( j=0; j< NumStep; j++)
	{
         double skippi;
         init.read((char*)&EnPr[1][i][j][0],sizeof(double));
         init.read((char*)&Grid[1][i][j][0][0],sizeof(double));
         init.read((char*)&EnPr[1][i][j][1],sizeof(double));
         //FOR READING WITH ETA INVARIANCE
         for(int Ent=0; Ent<3; Ent++)
         {
            init.read((char*)&V[1][i][j][Ent],sizeof(double));
         }

         init.read((char*)&PiGrd[1][i][j][0],sizeof(double));
         init.read((char*)&PiGrd[1][i][j][1],sizeof(double));
         init.read((char*)&PiGrd[1][i][j][2],sizeof(double));
         init.read((char*)&PiGrd[1][i][j][3],sizeof(double));
         init.read((char*)&PiGrd[1][i][j][4],sizeof(double));
         init.read((char*)&PiGrd[1][i][j][5],sizeof(double));//bulk
         


	for( k=1; k< ((2*ETAMAX/ZSTEP)+1+0.5); k++) // Propogate T in eta
	{
         for(int Ent=0; Ent<1; Ent++)
         {
            Grid[1][i][j][k][Ent]=Grid[1][i][j][0][Ent];
         }

      }
      }
      }
      init.close();

    }

    double TEMP (double e)
    {
        double Free = 47.5;  //considering all favors or quarks and gluons
        double factor = 3.0*Free*PIE*PIE/90.0;
        cout << e << endl;
        return ( pow(e/factor,0.25));
    }

    void GridGen (int t)
    {
       for (int i=0; i<NumStep; i++)
       {
          for (int j=0; j<NumStep; j++)
          {
             for (int k=0; k<NumStep; k++)
             {
                Grid[0][i][j][k][0] = Grid[1][i][j][k][0];
                Grid[1][i][j][k][0] = sqrt((i-49)*(i-49) + (j-49)*(j-49) + (k-49)*(k-49) + (t-49)*(t-49))*SStep;
                //SStep = 0.01
             }
          }
       }
    }

    void GridGen2 (int t)
    {
       for (int i=0; i<NumStep; i++)
       {
          for (int j=0; j<NumStep; j++)
          {
             for (int k=0; k<((2*ETAMAX/ZSTEP)+1+0.5-SCoarseGrain); k++)
             {
                 Grid[0][i][j][k][0] = Grid[1][i][j][k][0];
                 
                 if(t<=4) //Sets tau = 0.9= 0.1 + 4*TAUSTEPFM
                 {
                    Grid[1][i][j][k][0] = 100;
                 }
                 else
                 {
                    Grid[1][i][j][k][0] = 200;
                 }
             }
          }
       }
    }

    void GridGen3 (int t)
    {
       for (int i=0; i<NumStep; i++)
       {
          for (int j=0; j<NumStep; j++)
          {
             for (int k=0; k<((2*ETAMAX/ZSTEP)+0.5); k++)
             {
                 Grid[0][i][j][k][0] = Grid[1][i][j][k][0];
                 
                 if(t<=4) //Sets tau = 0.9= 0.1 + 4*TAUSTEPFM
                 {
                    Grid[1][i][j][k][0] = 200;
                 }
                 else
                 {
                    Grid[1][i][j][k][0] = 100;
                 }
             }
          }
       }
    }

   void GridSqua(int t)
    {
       for (int i=0; i<NumStep; i++)
       {
          for (int j=0; j<NumStep; j++)
          {
             for (int k=0; k<NumStep; k++)
             {
                Grid[0][i][j][k][0] = Grid[1][i][j][k][0];
                if ( t<88 )
                {
                   Grid[1][i][j][k][0] = 1+(((double)rand() / RAND_MAX)-0.5)*0.01;
                }
                
                else
                {
                   Grid[1][i][j][k][0] = 0;
                }
             }
          }
       }
    }

    /********************/
    /******* MAIN *******/
    /****** STARTS ******/
    /********************/



int main()
{
    srand(time(0));
    int l;
    DefGrid(&Grid);
    DefPi(&PiGrd);
    DefV(&V);
    DefEnPr(&EnPr);
    hyperbolic = 1; //0=cartesian input; 1=hyperbolic input
    
    /* Definition of values */
    /* "Lower t"cube */
    Link[0][0] = 0;
    Link[0][1] = 1; 
    Link[1][0] = 0;
    Link[1][1] = 2;
    Link[2][0] = 1;
    Link[2][1] = 3;
    Link[3][0] = 2;
    Link[3][1] = 3;
    Link[4][0] = 0;
    Link[4][1] = 4;
    Link[5][0] = 1;
    Link[5][1] = 5;
    Link[6][0] = 2;
    Link[6][1] = 6;
    Link[7][0] = 3;
    Link[7][1] = 7;
    Link[8][0] = 4;
    Link[8][1] = 5;
    Link[9][0] = 4;
    Link[9][1] = 6;
    Link[10][0] = 5;
    Link[10][1] = 7;
    Link[11][0] = 6;
    Link[11][1] = 7;

    /* Temporal Edges */
    Link[12][0] = 0;
    Link[12][1] = 8;
    Link[13][0] = 1;
    Link[13][1] = 9;
    Link[14][0] = 2;
    Link[14][1] = 10;
    Link[15][0] = 3;
    Link[15][1] = 11;
    Link[16][0] = 4;
    Link[16][1] = 12;
    Link[17][0] = 5;
    Link[17][1] = 13;
    Link[18][0] = 6;
    Link[18][1] = 14;
    Link[19][0] = 7;
    Link[19][1] = 15;
    
    /* "Upper t" cube */
    Link[20][0] = 8;
    Link[20][1] = 9; 
    Link[21][0] = 8;
    Link[21][1] = 10;
    Link[22][0] = 9;
    Link[22][1] = 11;
    Link[23][0] = 10;
    Link[23][1] = 11;
    Link[24][0] = 8;
    Link[24][1] = 12;
    Link[25][0] = 9;
    Link[25][1] = 13;
    Link[26][0] = 10;
    Link[26][1] = 14;
    Link[27][0] = 11;
    Link[27][1] = 15;
    Link[28][0] = 12;
    Link[28][1] = 13;
    Link[29][0] = 12;
    Link[29][1] = 14;
    Link[30][0] = 13;
    Link[30][1] = 15;
    Link[31][0] = 14;
    Link[31][1] = 15;
    
    /* Levi Civita */

    LC[0][1][2][3]=1.; LC[0][2][3][1]=1.; LC[0][3][1][2]=1.; LC[1][0][3][2]=1.;
    LC[1][2][0][3]=1.; LC[1][3][2][0]=1.; LC[2][0][1][3]=1.; LC[2][1][3][0]=1.; 
    LC[2][3][0][1]=1.; LC[3][0][2][1]=1.; LC[3][1][0][2]=1.; LC[3][2][1][0]=1.; 
    
    LC[0][1][3][2]=-1.; LC[0][2][1][3]=-1.; LC[0][3][2][1]=-1.; LC[1][0][2][3]=-1.;
    LC[1][2][3][0]=-1.; LC[1][3][0][2]=-1.; LC[2][0][3][1]=-1.; LC[2][1][0][3]=-1.; 
    LC[2][3][1][0]=-1.; LC[3][0][1][2]=-1.; LC[3][1][2][0]=-1.; LC[3][2][0][1]=-1.; 



    Crit = 150.;

 //   CritTemp = 0.12; //in inverse gev

    
    /*********/
    /*PRELOAD*/
    /**STEP 0*/
    /*********/

    GridGen3(0);

    /* Find Ground energy level */

    for(int i=0; i<(NumStep); i++)
    {
        for(int j=0; j<(NumStep); j++)
        {
            for(int k=0; k<(21); k++)
            {
                if(Grid[1][i][j][k][0] >= Crit)
                {
                    GroundEnergy = GroundEnergy + Grid[1][i][j][k][0]*SStep*SStep;
                }
            }
        }
    }

    /********************************************************************************************/
    /************************************* PROCESS **********************************************/
    /*************************************  LOOP   **********************************************/
    /************************************* BEGINS! **********************************************/
    /********************************************************************************************/
ofstream PrintArea;
PrintArea.open("surf.dat");

    double totes=0;
    int count = 0;

   for (int GridT=0; GridT<(TAUCOUNTMAX-TCoarseGrain); GridT=GridT+TCoarseGrain)
    {
        /* Initialize values */
        GridPos[0] = GridT;
        GridPos[1] = 0;
        GridPos[2] = 0;
        GridPos[3] = 0;
        NArea = 0;
        
        //define tau so we can use it
        double tau = (GridT*TAUSTEPFM+TAUSTART);
        cout << "tau: " << tau << endl; 
        cout << "Grid Time: " << GridT << "\n";

      /* Read Data for NEW FUTURE timestep*/

      GridGen3(GridT + TCoarseGrain);


      for(GridPos[1]=0; GridPos[1]<(NumStep-SCoarseGrain); GridPos[1]=GridPos[1]+SCoarseGrain)
            {
             double x = GridPos[1]*XSTEP;
             for(GridPos[2]=0; GridPos[2]<(NumStep-SCoarseGrain); GridPos[2]=GridPos[2]+SCoarseGrain)
                {
                  double y = GridPos[2]*YSTEP;
                  for(GridPos[3]=0; GridPos[3]<((2*ETAMAX/ZSTEP)+0.5-SCoarseGrain); GridPos[3]=GridPos[3]+SCoarseGrain)
                    {
                       double eta = (GridPos[3]*ZSTEP-ETAMAX); // def eta
FailSafe=0;

    NPoint = 0;
    NTetra = 0;
    
    /* Get current Vertex Values */
    
    for(int i=0;i<4;i++)
    {
        if(i==0)
        {
        Vert[0][i] = Grid[ 0 ] [ (GridPos[1]) ] [ GridPos[2] ] [ GridPos[3] ] [i];
        Vert[1][i] = Grid[ 0 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [ GridPos[3] ] [i];
        Vert[2][i] = Grid[ 0 ] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain ] [ GridPos[3] ] [i];
        Vert[3][i] = Grid[ 0 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [ GridPos[3] ] [i];
        Vert[4][i] = Grid[ 0] [ (GridPos[1]) ] [ GridPos[2] ] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[5][i] = Grid[ 0] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[6][i] = Grid[ 0] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[7][i] = Grid[ 0] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [ GridPos[3] + SCoarseGrain ] [i];

        Vert[8][i] = Grid[ 1 ] [ (GridPos[1]) ] [ GridPos[2] ] [ GridPos[3] ] [i];
        Vert[9][i] = Grid[ 1 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [ GridPos[3] ] [i];
        Vert[10][i] = Grid[ 1 ] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain ] [ GridPos[3] ] [i];
        Vert[11][i] = Grid[ 1 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [ GridPos[3] ] [i];
        Vert[12][i] = Grid[ 1] [ (GridPos[1]) ] [ GridPos[2] ] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[13][i] = Grid[ 1] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[14][i] = Grid[ 1] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain] [ GridPos[3] + SCoarseGrain ] [i];
        Vert[15][i] = Grid[ 1] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [ GridPos[3] + SCoarseGrain ] [i];
        }
        else
        {
        Vert[0][i] = V[ 0 ] [ (GridPos[1]) ] [ GridPos[2] ] [i-1];
        Vert[1][i] = V[ 0 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [i-1];
        Vert[2][i] = V[ 0 ] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain ] [i-1];
        Vert[3][i] = V[ 0 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [i-1];
        Vert[4][i] = V[ 0] [ (GridPos[1]) ] [ GridPos[2] ] [i-1];
        Vert[5][i] = V[ 0] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [i-1];
        Vert[6][i] = V[ 0] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain] [i-1];
        Vert[7][i] = V[ 0] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [i-1];

        Vert[8][i] = V[ 1 ] [ (GridPos[1]) ] [ GridPos[2] ] [i-1];
        Vert[9][i] = V[ 1 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [i-1];
        Vert[10][i] = V[ 1 ] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain ] [i-1];
        Vert[11][i] = V[ 1 ] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [i-1];
        Vert[12][i] = V[ 1] [ (GridPos[1]) ] [ GridPos[2] ] [i-1];
        Vert[13][i] = V[ 1] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] ] [i-1];
        Vert[14][i] = V[ 1] [ (GridPos[1]) ] [ GridPos[2] + SCoarseGrain] [i-1];
        Vert[15][i] = V[ 1] [ (GridPos[1]) + SCoarseGrain] [ GridPos[2] + SCoarseGrain] [i-1];
        }
    }

    /*Fill out BVert*/
    for(l=0; l<16; l++)
    {     
        BVert[l] = (Vert[l][0] >= Crit);
    }

    
    /**********************/
    /******RELEASE ME!*****/
    /**********************/
    
    CubeAnalyzer (tau,x,y,eta);


    if (NPoint==4)
       {
            for(l=0; l<4; l++) //define spanning vectors
            {
                 SpanV[0][l] = Point[0][l] - Point[1][l];
                 SpanV[1][l] = Point[0][l] - Point[2][l];
                 SpanV[2][l] = Point[0][l] - Point[3][l];
            }
            
            Patch[0][0]=0; Patch[0][1]=0; Patch[0][2]=0; Patch[0][3]=0;
            
            for(int mu=0; mu<4; mu++) // general cross product of vectors
            {
               for(int nu=0; nu<4; nu++)
               {
                  for(int rho=0; rho<4; rho++)
                  {
                     for(int sig=0; sig<4; sig++)
                     {
                        Patch[0][sig] = Patch[0][sig] + (1./6.)*( LC[mu][nu][rho][sig]*SpanV[0][mu]*SpanV[1][nu]*SpanV[2][rho]  );
                     }
                  }
               }
            }
            
            // dot product test
            if(!Dot(Patch[0][0],Patch[0][1],Patch[0][2],Patch[0][3],V_LE[0],V_LE[1],V_LE[2],V_LE[3],0) )
            {
                 Patch[0][0] = - Patch[0][0];
                 Patch[0][1] = - Patch[0][1];
                 Patch[0][2] = - Patch[0][2];
                 Patch[0][3] = - Patch[0][3];
            }



            /* Add final area vector to the list */



            for(l=0; l<4; l++) // record surface vector
            {
                 Area[NArea][l] = Patch[0][l];
            }
            
            for(l=4; l<8; l++) //record surface element location 
            {
                 Area[NArea][l] = (Point[0][l-4] + Point[1][l-4] + Point[2][l-4] + Point[3][l-4])/4.; //Average position
            }
            
            NArea++;
       }



    /* The case of more than four points */
    else if (NPoint > 4)
    {
         
       FindTetras (); //Find tall tetrahedra

       Area[NArea][0] = 0;
       Area[NArea][1] = 0;
       Area[NArea][2] = 0;
       Area[NArea][3] = 0;

       for (l=0; l<NTetra; l++) // Check each tetra
       {
          //if the tetra hedra is exterior AND area points toward low energy, it contributes to the cell's dsigma
          if( (Tetra[l][5] == 0) && Dot(TetraVec[l][0], TetraVec[l][1], TetraVec[l][2], TetraVec[l][3], V_LE[0], V_LE[1], V_LE[2], V_LE[3],0) )
          {
             Area[NArea][0] = Area[NArea][0] + TetraVec[l][0];
             Area[NArea][1] = Area[NArea][1] + TetraVec[l][1];
             Area[NArea][2] = Area[NArea][2] + TetraVec[l][2];
             Area[NArea][3] = Area[NArea][3] + TetraVec[l][3];                 
          }
            
          //if not, mark as pointing in wrong direction
          else if (Tetra[l][5] == 0)
          {
             Tetra[l][5] = 2;
          }
       }
         
       //grid position of the cell
       for(l=4; l<8; l++) //record surface element location 
       {
          Area[NArea][l] = (Point[0][l-4] + Point[1][l-4] + Point[2][l-4] + Point[3][l-4])/4.; //Average position
       }

       NArea++;
    }
    
    if(NPoint > 3) //If there is a surface element only
    {
       GridPlot[GridCount][0]=GridT;
       GridPlot[GridCount][1]=GridPos[1];
       GridPlot[GridCount][2]=GridPos[2];
       GridPlot[GridCount][3]=GridPos[3];
       GridPlot[GridCount][4]=Area[NArea-1][0];
       GridPlot[GridCount][5]=Area[NArea-1][1];
       GridPlot[GridCount][6]=Area[NArea-1][2];
       GridPlot[GridCount][7]=Area[NArea-1][3];

       GridCount++;
    }
   }   //Closes Z Loop
   }   //Closes Y Loop
   }   //Closes X Loop

       count = count + NArea;
       for(int i=0; i<NArea; i++)
       {
          /* Below is only valid for cartesian systems*/
          totes = totes + sqrt( Area[i][0]*Area[i][0] +
                                Area[i][1]*Area[i][1] +
                                Area[i][2]*Area[i][2] +
                                Area[i][3]*Area[i][3]
                              );

          PrintArea << Area[i][4] << " " << Area[i][5] << " " << Area[i][6] << " " << Area[i][7] << " ";
          PrintArea << Area[i][0] << " " << Area[i][1] << " " << Area[i][2] << " " << Area[i][3] << " ";
          PrintArea << U[1] << " " << U[2] << " " << U[3] << endl;
       }
      
       cout << "Total (Cartesian) area:" << "   " << totes << endl;  
    }  //Closes Time Loop

//GridPlot Temp

   ofstream GridPrint;
   ofstream AreaPrint;
   GridPrint.open ("GridTestInfo.txt");
   AreaPrint.open ("AreaTestInfo.txt");

   for(int i=0; i<GridCount; i++)
   {
      GridPrint << GridPlot[i][0] << "   " << GridPlot[i][1] << "   " << GridPlot[i][2] << "   " << GridPlot[i][3] <<  endl;
      AreaPrint << GridPlot[i][4] << "   " << GridPlot[i][5] << "   " << GridPlot[i][6] << "   " << GridPlot[i][7] <<  endl;
   }





PrintArea.close();

cout << "Count: " << count <<endl;

cout << "I tried. " << ITry << endl;
cout << "I faled. " << IFale << endl;
cout << "Round 2" << endl;
cout << "I tried. " << ITry2 << endl;
cout << "I faled. " << IFale2 << endl;
cout << "Head in the clouds: " << spacey << endl;

RelGrid(&Grid);
}
