///////////////////////////Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include <sys/time.h>
#include <pthread.h>
#include <sched.h>
#define NUMT     2

//////////////////////////Directives Definition

#define VER 0.35
#define STEREO 0

//#define VTGRAPE // use Virtualized GRAPE library
#define VTGRAPE
/*#define SOCK_ON*/
#define GL_ON
#define LAP_TIME
#define C_MASS
/*#define TELOP*/
#define SUBWIN
#define CROSS
#define INFO
/*#define SWAP_ENDIAN*/
#if defined(MDGRAPE3) || defined(VTGRAPE)
#define MDM 2      /* 0:host 2:m2 */
#endif
#define SPC 0
#define ST2 0
#define TIP5P 1
#define SYS 0 /* 0:NaCl 1:water(fcc) 2:water(ice) 3:water(ice2) 4:NaCl-water */
#define S_NUM_MAX 10*10*10*8
#define W_NUM_MAX 10*10*10*8
#define ZERO_P 1
#define V_SCALE 0
#define T_CONST 1
#define P_CONST 0
#define KNUM 5                    /* number of particle type */
#define VMAX 462 /*1535*/        /* max value of wave nubmer vector */
#define EFT 12000
#define my_min(x,y) ((x)<(y) ? (x):(y))
#define my_max(x,y) ((x)>(y) ? (x):(y))

#define PI M_PI              /* pi */
#define PIT M_PI*2.0         /* 2 * pi */
#define PI2 M_PI*M_PI        /* pi*pi */
#define IPI M_1_PI           /* 1/pi */
#define ISPI M_2_SQRTPI*0.5  /* 1 / sqrt(pi) */
#define X_PIX_SIZE 1024
#define Y_PIX_SIZE 786
#define DATA_NUM 100
#define DATA_NUM 100
//////////////////////////Global Variables
//////////Used in main
int np;
float sub_x,sub_y,sub_off;
int temp_max = 0,temp_ymax = 10;
/////////Used in display
float clear_color = 0.0;
float eye_width = 0.8;
float eye_len;
float trans[3] = {0.0, 0.0, 0.0};
float angle[3] = {0.0, 0.0, 0.0};
float m_matrix[16];
float i_matrix[16];
/////////Used in Init()
GLfloat color_table[10][4];
float r_table[5];
float trans0[3];
float matrix0[16];
double radius = 0.45;
/////////Used in mouse()
int mpos[2];
int mouse_l = 0;
int mouse_m = 0;
int mouse_r = 0;
////////////////////OpenGL Variables
GLubyte *pix;
double clip[6][4];
int ditail = 15;
int drow_flg[5] = {1,1,1,1,1};
/////////////////////Flags and others
int grape_flg = 0;//Flag for 0=CPU or GPU=1
int bond_flg = 0;
int sc_flg = 0;    /* 0:non  1:server 2:client */
int auto_flg = 0;
#ifdef LAP_TIME
struct timeval time_v;
double md_time,md_time0;
double disp_time,disp_time0;
double sock_time,sock_time0;
#endif
int ini_flg = 1;
static int temp_data[DATA_NUM];
////////////////////////For Init_MD()
/* for MD */

int sys_num = SYS;

int run_flg = 1;
int c_flg = 0;
int c_num = 0;
int velp_flg = 0;
double start_vl = -1;
double t_cd[3];
int w_add,s_add;
#define C_STEP 100

#ifdef LAP_TIME
int vflg = 3;
#endif
int kflg = 0;
int tflg = 0;

char k_file[50];
FILE *fp;

int md_step = 10;
int md_stepf = 0;
int m_clock = 0;
int b_clock = 1;
int timemx = -1;

double avo  = 6.0221367e+23;    /* avogdro's number (mol^-1) */
double kb   = 8.617080363e-5;   /* Boltzmann's number (eV K^-1) */
double e    = 1.60217733e-19;   /* unit charge */

double delt = .5e-15;          /* dt sec */
//double delt = 0.125e-15;          /* dt sec */
double sigma = 1.0e-10;         /* unit of length (m) */
double mass  = 3.8175e-26;      /* unit of mass (Kg) */
double epsv  = 14.39;           /* unit of energy (eV) */
double epsj;

double a_massi[KNUM];
double a_mass[4] = {
  22.989768,   /* Atomic weight of Na */
  35.4527,     /* Atomic weight of Cl */
  15.9994,     /* Atomic weight of O */
  1.00794};    /* Atomic weight of H */

double bond[3] = {.9572, 0.15}; /* distance of O-C and O-M */
double hoh_deg = 104.52;
double m_cdx[4];
double m_cdy[4];
double m_cdz[4];
double moi[3];                  /* moment of inertia */
double temp  = 1800;             /* temperature (K) */
double nden = -1;               /* density \AA^-3 */
double pres;
double ini_temp;
double  *cd;         /* position */
double  *vl;         /* velocity */
double  *fc;         /* force */
double  *fcc;
double *iphi;
double *ang;             /* angle */
double *agv;             /* angular velocity */
double *agvp;            /* angular velocity */
double *angh;            /* angle */
double *agvh;            /* angular velocity */
double *agvph;           /* angular velocity */
double *trq;             /* trque */
int *w_index;
int *w_rindex;
int *w_info;
int w_site;
int w_num,w_num3;
int s_num,s_num3;
int ws_num,ws_num3;
long *nig,*nli;
int *nig_data,*nig_num;
int *atype;          /* particle type */
                     /* 0:Na 1:Cl 2:O 3:H1 4:H2 5:M 6:L1 7:L2 8:C */
int atype_mat[20];
int atype_num[KNUM+4];  /* particle number of each type */
double tmrdp,jrdp;
double crdp,vclrdp;
double erdp;
double side0;
float side[3],sideh[3],iside[3];
double side_s[3],side_e[3];
double h,hsq,hsq2;
double tscale,sc;
double mtemp;
double rtemp;
double ekin,ekin1,ekin2;
double r,rd,rr,inr;
double vir;

double mpres,rpres;
double vol;
double lp=0;
double pist = 0.001;

double xs = 1.0;
double lq = .1;

double center_mass;
GLfloat moji_c[2][4]  = { 0.8, 0.8, 0.8, 1.0, 0.0, 0.0, 0.0, 1.0 };
int np = 2;
int npx,npy,npz;
int n1;
int n2;
int n3;

int nn = 0;
int nw = 0;

double pb;
double pc[2][2],pd[2][2],ipotro[2][2];
double pol[2][2];
double sigm[2][2];
//////////////* local *//////////////////

double neighbor_radius = 3.1;
double min_angle = 15.0;
double max_angle = 75.0;

char keiname[256];
double z[KNUM+4],zz[KNUM+4][KNUM+4];
double wpa,wpc;
double as_s[KNUM][KNUM];
double as_e[KNUM][KNUM];
double as_a[KNUM][KNUM];
double as_c[KNUM][KNUM];
int vmax;
double oalpha = 6, alpha , alpha2, ial2si2;
float *erfct;
int *vecn[VMAX];
int knum=KNUM;
#if MDM != 0
  double gscale[(KNUM+4)*(KNUM+4)];
  double rscale[(KNUM+4)*(KNUM+4)];
  double gscale2[(KNUM+4)*(KNUM+4)];
  double rscale2[(KNUM+4)*(KNUM+4)];

  double charge[(KNUM+4)*(KNUM+4)];
  double roffset[(KNUM+4)*(KNUM+4)];

  double cellsize[3];
  double vecr;
#endif
#if MDM == 2
  double side_min,side_max;
  char f_table_name[50];
  char p_table_name[50];
#endif
double phir_corr;
double phi[3],phir;
int pcun = 1;
int p_count = 0;
////////////////////Structures
#define TIMETABLE_MAX 10000
typedef struct{
  int mouse[3];
  double move[3];
  double rot[3];
  char command;
  double temp;
  double matrix[16];
} TIMETABLE;
TIMETABLE *tt;


struct thread_data{
	int thread_id;

	int range,n1,*nig_num, n3, *atype_mat,*atype;
	double *cd,r, rd, rr, inr,phir,pb,*pol,*sigm;
	double *pc, *pd, *ipotro,*zz;

	double vir,*fc;
};

struct thread_data thread_data_array[NUMT];

////////////////////////Other Functions
//////////////////////////////////////////////
void *CPU_force(void *threadarg){
	int i, j, k, c;
	int i0, i1, i2, i3, i4, i5;
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12;
	double dphir;


	int id,range,n1,*nig_num, n3,*atype_mat,*atype;
	double *cd,r, rd,inr,phir,pb,*pol,*sigm;
	double *pc, *pd, *ipotro,*zz;

	double vir,*fc;

	struct thread_data *my_data;

	my_data = (struct thread_data *) threadarg;
	range = my_data->range;
	id = my_data->thread_id;

	n1 = my_data->n1;
	nig_num = my_data->nig_num;
	n3 = my_data->n3;
	atype_mat = my_data->atype_mat;
	atype = my_data->atype;
	cd = my_data->cd;
	r = my_data->r;
	rd = my_data->rd;

	inr = my_data->inr;
	phir = my_data->phir;
	pb = my_data-> pb;
	pol = my_data-> pol;
	sigm =  my_data->sigm;
	pc = my_data->pc;
	pd = my_data->pd;
	ipotro = my_data->ipotro;
	zz = my_data->zz;
	vir = my_data->vir;
	fc = my_data->fc;



	for (i = 0; i < n1; i++) nig_num[i] = 0;

		for (i = 0; i < n3; i += 3) {
			i0 = atype_mat[atype[i / 3]];
			for (j = i + 3; j < n3; j += 3) {
				d0 = cd[i] - cd[j];
				d1 = cd[i + 1] - cd[j + 1];
				d2 = cd[i + 2] - cd[j + 2];

				rd = d0 * d0 + d1 * d1 + d2 * d2;
				r = sqrt(rd);
				inr = 1. / r;

				i1 = atype_mat[atype[j / 3]];
				d7 = phir;

				if (i0 < 2 && i1 < 2) {
					 	d3 = pb * pol[i0+i1*2] *  exp((sigm[i0+i1*2] - r)  *  ipotro[i0+i1*2]);
						dphir = (d3 * ipotro[i0+i1*2] * inr
								- 6 * pc[i0+i1*2] * pow(inr, 8)
								- 8 * pd[i0+i1*2] * pow(inr, 10)
								+ inr * inr * inr * zz[i0+i1*(KNUM + 4)]);

				}

				vir -= rd * dphir;

				d3 = d0 * dphir;
				d4 = d1 * dphir;
				d5 = d2 * dphir;

				fc[i] += d3;
				fc[i + 1] += d4;
				fc[i + 2] += d5;
				fc[j] -= d3;
				fc[j + 1] -= d4;
				fc[j + 2] -= d5;

				}
		}

	my_data->fc = fc;
	pthread_exit(NULL);
}


void ice_set(float *side)
{
  int i,j,k;
  int i0,i1,i2,i3;
  int c;
  double a;
  double s3,s6;
  double ang0,ang1,ang2,ang3;
  double d0,d1,d2,d3,d4,d5;

  s3 = sqrt(3.0);
  s6 = sqrt(6.0);
  a = 4.52;

  d0 = a;
  d1 = a*s3;
  d2 = a*s6/3.0*2.0;

  if(nden > 0)
    a /= pow(nden/(8.0/(d0*d1*d2)),1./3.);
  else
    printf("You must set number of density");

  for(i = 0; i < n3; i++)
    cd[i] = 0;
  for(i = 0; i < n1*4; i++)
    ang[i] = 0;

  i = 0;
  cd[i] = a*.5; cd[i+1] = s3/6.*a; cd[i+2] = s6/12.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] = s3/6.*a; cd[i+2] = s6/3.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] =-s3/2.*a; cd[i+2] = 0;
  atype[i/3] = 2;
  i += 3;
  cd[i] = a*.5; cd[i+1] =-s3/2.*a; cd[i+2] = (s6/3.+s6/12.)*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] = 0; cd[i+2] = 0;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] = 0; cd[i+2] = (s6/3.+s6/12.)*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] =-s3/3.*a; cd[i+2] = s6/12.*a;
  atype[i/3] = 2;
  i += 3;
  cd[i] = 0; cd[i+1] =-s3/3.*a; cd[i+2] = s6/3.*a;
  atype[i/3] = 2;

  c = 8*3;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 24; i3 += 3){
            cd[c]   = cd[i3]  +i2*a;
            cd[c+1] = cd[i3+1]+i1*a*s3;
            cd[c+2] = cd[i3+2]+i0*a*s6/3*2;
            atype[c/3] = 2;
            /*
               printf("%d 0 % f % f % f % f % f % f %d %d %d\n"
               ,c/3,cd[c],cd[c+1],cd[c+2]
               ,cd[i3],cd[i3+1],cd[i3+2],i0,i1,i2);
            */
            c += 3;
          }

  side[0] = a*npx;
  side[1] = a*s3*npy;
  side[2] = a*s6/3.0*2.0*npz;

  for(i = 0; i < 8*4; i += 4){
    ang0 = 120./180.*PI;
    ang1 = (90.-54.74)/180.*PI;
    ang2 = 0;
    ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
    ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
    ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
    ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  }

  i = 0*4;
  ang0 = 120./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 1*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = -109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 2*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 3*4;
  ang0 = 60./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 4*4;
  ang0 = 30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 5*4;
  ang0 = -60./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 6*4;
  ang0 = 120./180.*PI;
  ang1 = (90.-54.74)/180.*PI;
  ang2 = 0;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
  i = 7*4;
  ang0 = -30./180.*PI;
  ang1 = 90./180.*PI;
  ang2 = 109.47/2/180*PI;
  ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
  ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
  ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
  ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);

  c = 32;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 32; i3++){
            ang[c++] = ang[i3];
          }

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0];
      j = i0*4;
      c = w_info[i/3]*3;
      ang0 = ang[j  ];
      ang1 = ang[j+1];
      ang2 = ang[j+2];
      ang3 = ang[j+3];
      for(k = 0; k < w_site-1; k++){
        d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
            +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
        d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
            +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
        d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
            +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
            +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
        cd[k*3+c  ] = cd[i  ] + d0;
        cd[k*3+c+1] = cd[i+1] + d1;
        cd[k*3+c+2] = cd[i+2] + d2;
      }
    }
    /*
  for(i = 0; i < n3; i += 3)
    printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
  exit(0);
*/
  d0 = 100;  d1 = 0;
  d2 = 100;  d3 = 0;
  d4 = 100;  d5 = 0;
  for(i = 0; i <  n3; i+=3){
    /*    printf("%d %d % f % f % f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);*/
    if(cd[i]   < d0) d0 = cd[i];
    if(cd[i]   > d1) d1 = cd[i];
    if(cd[i+1] < d2) d2 = cd[i+1];
    if(cd[i+1] > d3) d3 = cd[i+1];
    if(cd[i+2] < d4) d4 = cd[i+2];
    if(cd[i+2] > d5) d5 = cd[i+2];
  }
  /*
     printf("%f %f %f  %f %f %f  %f %f %f\n",d0,d1,d1-d0,d2,d3,d3-d2,d4,d5,d5-d4);
     printf("%f %f %f\n",side[0]-(d1-d0),side[1]-(d3-d2),side[2]-(d5-d4));
     */
  for(i = 0; i <  n3; i+=3){
    cd[i]   += (side[0]-(d1-d0))/2.-d0;
    cd[i+1] += (side[1]-(d3-d2))/2.-d2;
    cd[i+2] += (side[2]-(d5-d4))/2.-d4;
    /*    printf("%d %d % f % f % f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);*/
  }

}
void ice_set2(float* side)
{
  int i,j,k;
  int i0,i1,i2,i3;
  int c;
  double ang0,ang1,ang2,ang3;
  double d0,d1,d2,d3,d4,d5;

  double l;

  l = side[0] / (npx * 2);

  cd[0]   = 0;
  cd[1]   = 0;
  cd[2]   = 0;
  cd[3  ] = l;
  cd[3+1] = l;
  cd[3+2] = 0;
  cd[6  ] = l;
  cd[6+1] = 0;
  cd[6+2] = l;
  cd[9  ] = 0;
  cd[9+1] = l;
  cd[9+2] = l;

  for(i = 4*3; i < 8*3; i++)
    cd[i] = cd[i-4*3] + l/2;

  for(i = 0; i < 8*3; i++)
    cd[i] += l/4;

  c = 8*3;
  for(i0 = 0; i0 < npz; i0++)
    for(i1 = 0; i1 < npy; i1++)
      for(i2 = 0; i2 < npx; i2++)
        if(i0 != 0 || i1 != 0 || i2 != 0)
          for(i3 = 0; i3 < 24; i3 += 3){
            cd[c]   = cd[i3]  +(double)i2*side[0]/npx;
            cd[c+1] = cd[i3+1]+(double)i1*side[1]/npy;
            cd[c+2] = cd[i3+2]+(double)i0*side[2]/npz;

            c += 3;
          }
   c = 0;
  for(i = 0; i < n1; i++){
    if(atype[i] == 2){
      if((i % 8) < 4){
        ang0 =-45*PI/180;
        ang1 = 90*PI/180;
        ang2 = 0*PI/180;
      } else {
        ang0 = 45*PI/180;
        ang1 = 90*PI/180;
        ang2 = 0*PI/180;
      }
      ang[c  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[c+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[c+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[c+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[c  ] = ang[c  ];
      angh[c+1] = ang[c+1];
      angh[c+2] = ang[c+2];
      angh[c+3] = ang[c+3];
      c += 4;
    }
  }

}
void fccset2(int lnp,double lside,double cod[])
{
  int     i,j,k,c;
  double  px,py,pz;
  double  l;

  l = lside / (lnp * 2);

  for(i = 0;i < lnp; i++)
    for(j = 0;j < lnp; j++)
      for(k = 0;k < lnp; k++){
        px = k * 2 * l;
        py = j * 2 * l;
        pz = i * 2 * l;
        c  = 4*3*(i*lnp*lnp + j*lnp + k);
        cod[c    ] = px;
        cod[c  +1] = py;
        cod[c  +2] = pz;
        cod[c+3  ] = px + l;
        cod[c+3+1] = py + l;
        cod[c+3+2] = pz;
        cod[c+6  ] = px + l;
        cod[c+6+1] = py;
        cod[c+6+2] = pz + l;
        cod[c+9  ] = px;
        cod[c+9+1] = py + l;
        cod[c+9+2] = pz + l;
      }
  for(i = 0;i < n3; i++)
    cod[i] += l / 2.;
}
void fccset_w(float* side)
{
  int     i,j,k,c,i0;
  double  px,py,pz;
  double  l;

  l = side[0] / (npx * 2);

  for(i = 0;i < npz; i++)
    for(j = 0;j < npy; j++)
      for(k = 0;k < npx; k++){
        px = k * 2 * l;
        py = j * 2 * l;
        pz = i * 2 * l;
        c  = 4*3*(i*npx*npy + j*npx + k);
        cd[c    ] = px;
        cd[c  +1] = py;
        cd[c  +2] = pz;
        cd[c+3  ] = px + l;
        cd[c+3+1] = py + l;
        cd[c+3+2] = pz;
        cd[c+6  ] = px + l;
        cd[c+6+1] = py;
        cd[c+6+2] = pz + l;
        cd[c+9  ] = px;
        cd[c+9+1] = py + l;
        cd[c+9+2] = pz + l;
      }

  for(i0 = 0; i0 < w_num; i0++){
    i = w_index[i0];
    cd[i]   += l / 2.;
    cd[i+1] += l / 2.;
    cd[i+2] += l / 2.;
  }

}
void velset6(double tref,double dh,double tscale,int knum,int num)
{
  int i,j,k,c;
  double u1,u2,v1,v2,s,d;
  double spx=0,spy=0,spz=0;
  double ekin,ts,sc;

  c = 0;
  i = 0;
  while(c < num){
    u1 = (double)rand() / (double)RAND_MAX;
    u2 = (double)rand() / (double)RAND_MAX;
    v1 = 2.*u1-1.;
    v2 = 2.*u2-1.;
    s  = v1*v1 + v2*v2;
    if( s < 1 ){
      while(atype[i/3] > 2) i++;
      vl[i++] = v1*(double)sqrt((-2*log(s))/s);
      c++;
      while(atype[i/3] > 2) i++;
      if(i < n3){
        vl[i++] = v2*(double)sqrt((-2*log(s))/s);
        c++;
      }
    }
  }

  for(k = 0; k < knum; k++){
    j = 0;
    spx = spy = spz = 0;
    for(i = 0; i < n3; i += 3){
      if(atype[i/3] == k){
        spx += vl[i  ];
        spy += vl[i+1];
        spz += vl[i+2];
        j++;
      }
    }
    if(j != 0){
      spx /= j;
      spy /= j;
      spz /= j;
    }
    for(i = 0;i < n3; i += 3){
      if(atype[i/3] == k){
        vl[i  ] -= spx;
        vl[i+1] -= spy;
        vl[i+2] -= spz;
      }
    }
  }

  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i]*vl[i] + vl[i+1]*vl[i+1] + vl[i+2]*vl[i+2])
      *a_mass[atype_mat[atype[i/3]]];
  }


  ts = tscale * ekin;
  sc = sqrt( tref / ts );
  sc *= dh;
  for(i = 0;i < n3; i += 3){
    vl[i  ] *= sc;
    vl[i+1] *= sc;
    vl[i+2] *= sc;
  }
}


double mass_den3(int xp, int xp2, int xm, int xm2, double comp, double temp)
{
  double nden1, nden2, nden;

  double p_atomnum[4] = {6.941, 22.989768, 39.0983, 85.4678}; /* atomic weight + */
  double m_atomnum[4] = {18.9984032, 35.4527, 79.904, 126.90447};

  /*  F       Cl      Br      I   */
  double a[4][4] = {2.3768, 1.8842, 3.0658, 3.7902, /* Li */ /* density */
                      2.655,  2.1393, 3.1748, 3.6274, /* Na */
                      2.6464, 2.1359, 2.9583, 3.3594, /* K  */
                      0.0000, 3.1210, 3.7390, 3.9449}; /* Rb */
  double b[4][4] = {0.4902, 0.4328, 0.6520, 0.9176,
                      0.560,  0.5430, 0.9169, 0.9491,
                      0.6515, 0.5831, 0.8253, 0.9557,
                      0.0000, 0.8832, 1.0718, 1.1435};

  double a_s[4][4] = {0,    0,     0,    0, /* density of solid */
                        0,    2.168, 0,    0,
                        0,    1.985, 0,    0,
                        0,    0,     0,    0};
  double b_s[4][4] = {0,    0,     0,    0,
                        0,    1.267, 0,    0,
                        0,    .5459, 0,    0,
                        0,    0,     0,    0};
  double c_s[4][4] = {0,    0,     0,    0,
                        0,    1.754, 0,    0,
                        0,    1.836, 0,    0,
                        0,    0,     0,    0};

  if(xp2 == -1){
    return((a[xp][xm] - b[xp][xm]*1e-3*temp)/((a_mass[0]+a_mass[1])/2*mass)*1e-27);
  } else {
    nden1 = (a[xp][xm] - b[xp][xm]*1e-3*temp)/((a_mass[0]+a_mass[1])/2*mass)*1e-27;
    nden2 = (a[xp2][xm2] - b[xp2][xm2]*1e-3*temp)/((a_mass[2]+a_mass[3])/2*mass)*1e-27;
    return( nden1*(100-comp)/100 + nden2*comp/100);
  }
}
double nden_set(double tmp)
{
  int i;
  double den[] = {   0.99984, 0.99990, 0.99994, 0.99996, 0.99997,
                       0.99996, 0.99994, 0.99990, 0.99985, 0.99978,
                       0.99970, 0.99961, 0.99949, 0.99938, 0.99924,
                       0.99910, 0.99894, 0.99877, 0.99860, 0.99841,
                       0.99820, 0.99799, 0.99777, 0.99754, 0.99730,
                       0.99704, 0.99678, 0.99651, 0.99623, 0.99594,
                       0.99565, 0.99534, 0.99503, 0.99470, 0.99437,
                       0.99403, 0.99368, 0.99333, 0.99297, 0.99259,
                       0.99222, 0.99183, 0.99144, 0.99104, 0.99063,
                       0.99021, 0.98979, 0.98936, 0.98893, 0.98849,
                       0.98804, 0.98758, 0.98715, 0.98665, 0.98618,
                       0.98570, 0.98521, 0.98471, 0.98422, 0.98371,
                       0.98320, 0.98268, 0.98216, 0.98163, 0.98110,
                       0.98055, 0.98001, 0.97946, 0.97890, 0.97834,
                       0.97777, 0.97720, 0.97662, 0.97603, 0.97544,
                       0.97485, 0.97425, 0.97364, 0.97303, 0.97242,
                       0.97180, 0.97117, 0.97054, 0.96991, 0.96927,
                       0.96862, 0.96797, 0.96731, 0.96665, 0.96600,
                       0.96532, 0.96465, 0.96379, 0.96328, 0.96259,
                       0.96190, 0.96120, 0.96050, 0.95979, 0.95906};

  if(temp >= 273 && temp <= 373)
    return(den[(int)(tmp+.5)]);
  else
    return 0.917;
}

void potpar5(int xp,int xp2,int xm,int xm2, char keiname[]){

  char gpname[4][3] = {"Li", "Na", "K",  "Rb"};
  char gmname[4][3] = {"F",  "Cl", "Br", "I"};

  int nip[4] = {2, 8, 8, 8};
  int nim[4] = {8, 8, 8, 8};

  double sigmp[4] = {.816, 1.17, 1.463, 1.587}; /* 0:Li 1:Na 2:K  3:Rb */
  double sigmm[4] = {1.179, 1.585, 1.716, 1.907}; /* 0:F  1:Cl 2:Br 3:I  */
  /*  F     Cl    Br    I  */
  double rho[4][4] = {.299, .342, .353, .430, /* Li */
                        .330, .317, .340, .386, /* Na */
                        .338, .337, .335, .355, /* K  */
                        .328, .318, .335, .337}; /* Rb */
  double cpp[4] = {0.073, 1.68, 24.3, 59.4};
  double cmm[4][4] = {14.5, 111.0, 185.1, 378.0,
                        16.5, 116.0, 196.0, 392.0,
                        18.6, 124.5, 206.0, 403.0,
                        18.9, 130.0, 215.0, 428.0};
  double cpm[4][4] = { 0.8,  2.0,  2.5,   3.3,
                         4.5, 11.2, 13.0,  19.1,
                         19.5, 48.0, 60.0,  82.0,
                         31.0, 79.0, 99.0, 135.0};

  double dpp[4] = { 0.03, 0.8, 24.0, 82.0};
  double dmm[4][4] = {17, 223, 423, 1060,
                        20, 233, 450, 1100,
                        22, 250, 470, 1130,
                        23, 260, 490, 1200};
  double dpm[4][4] = { 0.6,   2.4,   3.3,   5.3,
                         3.8,  13.9,  19.0,  31.0,
                         21.0,  73.0,  99.0, 156.0,
                         40.0, 134.0, 180.0, 280.0};

  if(xp > 3 || xm > 3){
    printf("error\n");
    exit(1);
  }

  strcpy(keiname,gpname[xp]);
  strcat(keiname,gmname[xm]);
  if(xp2 != -1){
    strcat(keiname,gpname[xp2]);
    strcat(keiname,gmname[xm2]);
  }

  pc[1][1] = cmm[xp][xm];
  pc[1][0] = cpm[xp][xm];
  pc[0][1] = cpm[xp][xm];
  pc[0][0] = cpp[xp];
  if(xp2 != -1){
    pc[1][3] = sqrt(cmm[xp][xm]*cmm[xp2][xm2]);
    pc[3][1] = sqrt(cmm[xp][xm]*cmm[xp2][xm2]);
    pc[3][3] = cmm[xp2][xm2];
    pc[0][3] = cpm[xp][xm2];
    pc[3][0] = cpm[xp][xm2];
    pc[1][2] = cpm[xp2][xm];
    pc[2][1] = cpm[xp2][xm];
    pc[2][3] = cpm[xp2][xm2];
    pc[3][2] = cpm[xp2][xm2];
    pc[0][2] = sqrt(cpp[xp]*cpp[xp2]);
    pc[2][0] = sqrt(cpp[xp]*cpp[xp2]);
    pc[2][2] = cpp[xp2];
  }

  pd[1][1] = dmm[xp][xm];
  pd[0][1] = dpm[xp][xm];
  pd[1][0] = dpm[xp][xm];
  pd[0][0] = dpp[xp];
  if(xp2 != -1){
    pd[1][3] = sqrt(dmm[xp][xm]*dmm[xp2][xm2]);
    pd[3][1] = sqrt(dmm[xp][xm]*dmm[xp2][xm2]);
    pd[3][3] = dmm[xp2][xm2];
    pd[0][3] = dpm[xp][xm2];
    pd[3][0] = dpm[xp][xm2];
    pd[2][1] = dpm[xp2][xm];
    pd[1][2] = dpm[xp2][xm];
    pd[2][3] = dpm[xp2][xm2];
    pd[3][2] = dpm[xp2][xm2];
    pd[0][2] = sqrt(dpp[xp]*dpp[xp2]);
    pd[2][0] = sqrt(dpp[xp]*dpp[xp2]);
    pd[2][2] = dpp[xp2];
  }

  ipotro[1][1] = 1./rho[xp][xm];
  ipotro[0][1] = 1./rho[xp][xm];
  ipotro[1][0] = 1./rho[xp][xm];
  ipotro[0][0] = 1./rho[xp][xm];
  if(xp2 != -1){
    ipotro[1][3] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[3][1] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[3][3] = 1./rho[xp2][xm2];
    ipotro[0][3] = 1./rho[xp][xm2];
    ipotro[3][0] = 1./rho[xp][xm2];
    ipotro[2][1] = 1./rho[xp2][xm];
    ipotro[1][2] = 1./rho[xp2][xm];
    ipotro[2][3] = 1./rho[xp2][xm2];
    ipotro[3][2] = 1./rho[xp2][xm2];
    ipotro[0][2] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[2][0] = 1./((rho[xp][xm]+rho[xp2][xm2])/2);
    ipotro[2][2] = 1./rho[xp2][xm2];
  }

  sigm[1][1] = sigmm[xm]*2;
  sigm[0][1] = sigmp[xp] + sigmm[xm];
  sigm[1][0] = sigmp[xp] + sigmm[xm];
  sigm[0][0] = sigmp[xp]*2;
  if(xp2 != -1){
    sigm[1][3] = sigmm[xm]+sigmm[xm2];
    sigm[3][1] = sigmm[xm]+sigmm[xm2];
    sigm[3][3] = sigmm[xm2]*2;
    sigm[0][3] = sigmp[xp] + sigmm[xm2];
    sigm[3][0] = sigmp[xp] + sigmm[xm2];
    sigm[1][2] = sigmp[xp2] + sigmm[xm];
    sigm[2][1] = sigmp[xp2] + sigmm[xm];
    sigm[2][3] = sigmp[xp2] + sigmm[xm2];
    sigm[3][2] = sigmp[xp2] + sigmm[xm2];
    sigm[0][2] = sigmp[xp]+sigmp[xp2];
    sigm[2][0] = sigmp[xp]+sigmp[xp2];
    sigm[2][2] = sigmp[xp2]*2;
  }

  pol[1][1] = -1./nim[xm] - 1./nim[xm] + 1;
  pol[0][1] =  1./nip[xp] - 1./nim[xm] + 1;
  pol[1][0] =  1./nip[xp] - 1./nim[xm] + 1;
  pol[0][0] =  1./nip[xp] + 1./nip[xp] + 1;
  if(xp2 != -1){
    pol[1][3] =  -1./nim[xm] - 1./nim[xm2] + 1;
    pol[3][1] =  -1./nim[xm] - 1./nim[xm2] + 1;
    pol[3][3] =  -1./nim[xm2] - 1./nim[xm2] + 1;
    pol[0][3] =  1./nip[xp] - 1./nim[xm2] + 1;
    pol[3][0] =  1./nip[xp] - 1./nim[xm2] + 1;
    pol[1][2] =  1./nip[xp2] - 1./nim[xm] + 1;
    pol[2][1] =  1./nip[xp2] - 1./nim[xm] + 1;
    pol[2][3] =  1./nip[xp2] - 1./nim[xm2] + 1;
    pol[3][2] =  1./nip[xp2] - 1./nim[xm2] + 1;
    pol[0][2] =  1./nip[xp] + 1./nip[xp2] + 1;
    pol[2][0] =  1./nip[xp] + 1./nip[xp2] + 1;
    pol[2][2] =  1./nip[xp2] + 1./nip[xp2] + 1;
  }
}
void CircleTable (float **sint,float **cost,const int n){
    int i;
    /* Table size, the sign of n flips the circle direction */
    const int size = abs(n);
    /* Determine the angle between samples */
    const float angle = 2*M_PI/(float)( ( n == 0 ) ? 1 : n );
    /* Allocate memory for n samples, plus duplicate of first entry at the end */
    *sint = (float *) calloc(sizeof(float), size+1);
    *cost = (float *) calloc(sizeof(float), size+1);
    /* Bail out if memory allocation fails, fgError never returns */
    if (!(*sint) || !(*cost))
    {
        free(*sint);
        free(*cost);
        printf("Failed to allocate memory in fghCircleTable");
		exit(0);
    }
    /* Compute cos and sin around the circle */
    (*sint)[0] = 0.0;
    (*cost)[0] = 1.0;
    for (i=1; i<size; i++)
    {
        (*sint)[i] = sin(angle*i);
        (*cost)[i] = cos(angle*i);
    }
    /* Last sample is duplicate of the first */
    (*sint)[size] = (*sint)[0];
    (*cost)[size] = (*cost)[0];
}


void mygluPerspective( GLfloat fovy, GLfloat aspect,
		GLfloat zNear, GLfloat zFar )
{
   GLfloat xmin, xmax, ymin, ymax;
   ymax = zNear * tan( fovy * M_PI / 360.0 );
   ymin = -ymax;
   xmin = ymin * aspect;
   xmax = ymax * aspect;
   glFrustum( xmin, xmax, ymin, ymax, zNear, zFar );

}

void mygluLookAt( GLdouble eyex, GLdouble eyey, GLdouble eyez,
	            GLdouble centerx, GLdouble centery, GLdouble centerz,
	            GLdouble upx, GLdouble upy, GLdouble upz )
	{
	   GLdouble m[16];
	   GLdouble x[3], y[3], z[3];
	   GLdouble mag;
	   /* Make rotation matrix */
	   /* Z vector */
	   z[0] = eyex - centerx;
	   z[1] = eyey - centery;
	   z[2] = eyez - centerz;
	   mag = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] );
	   if (mag) {  /* mpichler, 19950515 */
	      z[0] /= mag;
	      z[1] /= mag;
	      z[2] /= mag;
	   }
	   /* Y vector */
	   y[0] = upx;
	   y[1] = upy;
	   y[2] = upz;
	   /* X vector = Y cross Z */
	   x[0] =  y[1]*z[2] - y[2]*z[1];
	   x[1] = -y[0]*z[2] + y[2]*z[0];
	   x[2] =  y[0]*z[1] - y[1]*z[0];
	   /* Recompute Y = Z cross X */
	   y[0] =  z[1]*x[2] - z[2]*x[1];
	   y[1] = -z[0]*x[2] + z[2]*x[0];
	   y[2] =  z[0]*x[1] - z[1]*x[0];
	   /* mpichler, 19950515 */
	   /* cross product gives area of parallelogram, which is < 1.0 for
	    * non-perpendicular unit-length vectors; so normalize x, y here
	    */
	   mag = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
	   if (mag) {
	      x[0] /= mag;
	      x[1] /= mag;
	      x[2] /= mag;
	   }
	   mag = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] );
	   if (mag) {
	      y[0] /= mag;
	      y[1] /= mag;
	      y[2] /= mag;
	   }
	#define M(row,col)  m[col*4+row]
	   M(0,0) = x[0];  M(0,1) = x[1];  M(0,2) = x[2];  M(0,3) = 0.0;
	   M(1,0) = y[0];  M(1,1) = y[1];  M(1,2) = y[2];  M(1,3) = 0.0;
	   M(2,0) = z[0];  M(2,1) = z[1];  M(2,2) = z[2];  M(2,3) = 0.0;
	   M(3,0) = 0.0;   M(3,1) = 0.0;   M(3,2) = 0.0;   M(3,3) = 1.0;
	#undef M
	   glMultMatrixd( m );
	   /* Translate Eye to Origin */
	   glTranslated( -eyex, -eyey, -eyez );
}
void mat_inv(float a[4][4])
{
  int i,j,k;
  float t, u, det;
  int n = 3;

  det = 1;
  for(k = 0; k < n; k++){
    t = a[k][k]; det *= t;
    for(i = 0; i < n; i++) a[k][i] /= t;
    a[k][k] = 1 / t;
    for(j = 0; j < n; j++)
      if(j != k){
        u = a[j][k];
        for(i = 0; i < n; i++)
          if(i != k) a[j][i] -= a[k][i] * u;
          else       a[j][i] = -u/t;
      }
  }
}
//////////////////////////////// OpenGL Vital Functions
///////////////////////////////////////////////////////


void init(void)
{
  int i,j;

  GLfloat mat_specular[] = {0.2, 0.2, 0.2, 1.0};
  GLfloat mat_ambient[] = {0.1, 0.1, 0.1, 1.0};
  GLfloat mat_shininess[] = {64.0};
  GLfloat light_position[] = {1.0, 1.1, 1.2, 0.0};

  float fm_matrix[16];
  float fi_matrix[16];

#if 0
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glShadeModel(GL_FLAT);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
#endif

#if 0
  glShadeModel(GL_SMOOTH);
/*  glShadeModel(GL_FLAT);*/
  glLightfv(GL_LIGHT0, GL_SPECULAR, mat_specular);
  glLightfv(GL_LIGHT0, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);


//  glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
//  glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#endif
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX,m_matrix);
  glGetFloatv(GL_MODELVIEW_MATRIX,i_matrix);






#if 1
  color_table[0][0] = 0.7;
  color_table[0][1] = 0.38;
  color_table[0][2] = 0.38;
  color_table[0][3] = 1;

  color_table[1][0] = 0.38;
  color_table[1][1] = 0.55;
  color_table[1][2] = 0.38;
  color_table[1][3] = 1;

  for(i = 0; i < 3; i++){
    color_table[0][i] /= 2.0;
    color_table[1][i] /= 2.0;
  }

  color_table[2][0] = 1;
  color_table[2][1] = .4;
  color_table[2][2] = 1;
  color_table[2][3] = 1;

  color_table[3][0] = 0;
  color_table[3][1] = 0.8;
  color_table[3][2] = 1;
  color_table[3][3] = 1;

  color_table[4][0] = 1;
  color_table[4][1] = 1;
  color_table[4][2] = 1;
  color_table[4][3] = 1;

  r_table[0] = 2.443/2;
  r_table[1] = 3.487/2;
  r_table[2] = 3.156/2;
  r_table[3] = .7;
  r_table[4] = .7;

#endif

#if 1
//#ifdef INFO
  for(i = 0; i < 3; i++)
    trans0[i] = 0;
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      if(i == j){
	matrix0[i*4+j] = 1;
      } else {
	matrix0[i*4+j] = 0;
      }
    }
  }
#endif
}
void keep_mem(int num_s, int num_w)
{
  int i,j;
  int add = 100;

  if((nli = (long*)malloc(20000 * sizeof(long))) == NULL){
    printf("memory error\n");
    exit(1);
  }

  if((nig_num = (int*)malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((nig_data = (int*)malloc((num_s+num_w+add)*6 * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }

  if((atype = malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((cd = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((vl = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((fc = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((fcc = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((iphi = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((ang = malloc((num_w+add) * 4 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agv = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvp = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((angh = malloc((num_w+add) * 4 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvh = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((agvph = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((trq = malloc((num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_index = malloc((num_w+add)/w_site * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_rindex = malloc((num_s+num_w+add) * 3 * sizeof(double))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  if((w_info = malloc((num_s+num_w+add) * sizeof(int))) == NULL){
    printf("memory error\n");
    exit(1);
  }

  if((erfct = malloc((EFT+1) * sizeof(float))) == NULL){
    printf("memory error\n");
    exit(1);
  }
  for(i = 0;i < VMAX; i++)
    if((vecn[i] = malloc(4 * sizeof(int))) == NULL){
      printf("memory error\n");
      exit(1);
    }
#ifdef GL_ON
  if((pix=(GLubyte*)malloc((X_PIX_SIZE+3)*(Y_PIX_SIZE)*3*sizeof(GLubyte)))==NULL){
    printf("memory error\n");
    exit(1);
  }
#endif

}
void init_MD(void)
{
  double d0,d1,d2;
  int i,j;

  srand( 1 );
  /*  srand( ( unsigned )time( NULL ));*/

  mass = a_mass[0]/avo*1e-3;

  for(i = 1; i < 4; i++){
    a_mass[i] /= a_mass[0];
  }
  a_mass[0] = 1.0;

  epsj  = epsv*1.60219e-19;
  crdp = sigma * 1e+10;
  tmrdp = sqrt(mass / epsj) * sigma;
  erdp = epsv * 2.30492e+1;       /* for calculate energy(kcal/mol) */

  keiname[0] = 0;

  if(sys_num == 0){
    if(tflg == 0)
      //temp  = 3000;
    delt = 2.0e-15;
  } else {
    if(tflg == 0)
      //temp  = 293;
    delt = .5e-15;
  }

  rtemp = temp / epsv * kb;
  h     = delt / tmrdp;
  hsq   = h * h;

  atype_mat[0] = 0; /* Na */
  atype_mat[1] = 1; /* Cl */
  z[0] = 1.0;
  z[1] =-1.0;


#if TIP5P  == 1
  d0 = 0.16;
  d1 = 3.12;
  wpa = d0*4*pow(d1,12)/2.30492e+1/epsv;
  wpc = d0*4*pow(d1, 6)/2.30492e+1/epsv;
  z[2] = 0; z[3] = 0.241; z[4] =-0.241;
  bond[0] = 0.9572;  bond[1] = 0.7;

  hoh_deg = 104.52;
  w_site = 5;
  atype_mat[2] = 2; /* O  */
  atype_mat[3] = 3; /* H1 */
  atype_mat[4] = 3; /* H2 */
  atype_mat[5] = 4; /* L1 */
  atype_mat[6] = 4; /* L2 */
  atype_mat[8] = 2; /* O  */
#endif

  for(i = 0; i < KNUM+4; i++)
    for(j = 0; j < KNUM+4; j++){
      zz[i][j] = z[i]*z[j];
    }

  pb = 0.338e-19/epsj;
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      pc[i][j] = 0;
      pd[i][j] = 0;
      pol[i][j] = 0;
      sigm[i][j] = 0;
      ipotro[i][j] = 0;
    }
  }
  potpar5(1,-1,1,-1,keiname);

  for(i = 0;i < 2; i++){
    for(j = 0;j < 2; j++){
      pc[i][j] *= 1e-79/epsj/pow(sigma,6);
      pd[i][j] *= 1e-99/epsj/pow(sigma,8);
    }
  }

  /* 0:Na 1:Cl 2:O 3:H */
  as_s[0][0] = 2.443; as_s[0][1] = 2.796; as_s[0][2] = 2.72; as_s[0][3] =1.310;
  as_s[1][0] = 2.796; as_s[1][1] = 3.487; as_s[1][2] = 3.55; as_s[1][3] =2.140;
  as_s[2][0] = 2.72;  as_s[2][1] = 3.55;  as_s[2][2] = 3.156;as_s[2][3] =0.0;
  as_s[3][0] = 1.310; as_s[3][1] = 2.140; as_s[3][2] = 0.0;  as_s[3][3] =0.0;

  as_e[0][0]=0.11913; as_e[0][1]= 0.3526;as_e[0][2]=0.56014;as_e[0][3]=0.56014;
  as_e[1][0]=0.3526;  as_e[1][1]=0.97906;as_e[1][2]=1.50575;as_e[1][3]=1.50575;
  as_e[2][0]=0.56014; as_e[2][1]=1.50575;as_e[2][2]=0.65020;as_e[2][3]=0.0;
  as_e[3][0]=0.56014; as_e[3][1]=1.50575;as_e[3][2] = 0.0;  as_e[3][3]=0.0;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++){
      as_e[i][j] *= 4.*1000. * 1.0364272e-5 /epsv;
      as_a[i][j] = as_e[i][j]*pow(as_s[i][j],12);
      as_c[i][j] = as_e[i][j]*pow(as_s[i][j], 6);
    }

  if(sys_num == 0){
    n1 = np*np*np*8;
  } else if(sys_num == 1 || sys_num == 10){
    n1 = np*np*np*4*w_site;
  } else if(sys_num == 2){
#if ZERO_P == 1
    n1 = (npx*npy*npz*8+npy*npz*4+(npy-1)*npz*8)*w_site;
#else

#endif
  } else if(sys_num == 3){
    n1 = np*np*np*8*w_site;
  } else if(sys_num == 4){
    if(np*np*np*8-nn*2 >= 0)
      n1 = (np*np*np*8-nn*2)*w_site+nn*2;
    else {
      printf("nn is too large!\n");
      exit(0);
    }
  } else if(sys_num == 5){
    if(np*np*np*8-nw >= 0)
      n1 = (np*np*np*8-nw)+nw*w_site;
    else {
      printf("nw is too large!\n");
      exit(0);
    }
  } else if(sys_num == 6){
    if(np-nw >= 0){
      if(nw > 0)
        n1 = np*np*np*8+(np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8)*(w_site-1);
      else
        n1 = np*np*np*8;
    } else {
      printf("nw is too large!\n");
      exit(0);
    }
  }

  n2 = n1 * 2;
  n3 = n1 * 3;

  m_cdx[0] =  bond[0]*sin(hoh_deg/2/180*PI);
  m_cdy[0] = -bond[0]*cos(hoh_deg/2/180*PI);
  m_cdz[0] = 0;
  m_cdx[1] = -bond[0]*sin(hoh_deg/2/180*PI);
  m_cdy[1] = -bond[0]*cos(hoh_deg/2/180*PI);
  m_cdz[1] = 0;

#if TIP5P == 1
  d0 = 109.47;
  m_cdx[2] =  0;
  m_cdy[2] =  bond[1]*cos(d0/2/180*PI);
  m_cdz[2] =  bond[1]*sin(d0/2/180*PI);
  m_cdx[3] =  0;
  m_cdy[3] =  bond[1]*cos(d0/2/180*PI);
  m_cdz[3] = -bond[1]*sin(d0/2/180*PI);
#endif

  d0 = d1 = d2= 0;
  for(i = 0; i < 2; i++){
    d0 += m_cdx[i]*a_mass[3];
    d1 += m_cdy[i]*a_mass[3];
    d2 += m_cdz[i]*a_mass[3];
  }
  center_mass = -d1/(a_mass[2]+a_mass[3]*2.0);

#ifdef C_MASS
  moi[0] = (a_mass[2]*center_mass*center_mass +
	    a_mass[3]*(m_cdy[0]+center_mass)*(m_cdy[0]+center_mass)*2.0);
  moi[1] = a_mass[3]*m_cdx[0]*m_cdx[0]*2.0;
  d0 = sqrt(bond[0]*bond[0]+center_mass*center_mass
	    -2.0*bond[0]*center_mass*cos(hoh_deg/2/180*PI));
  moi[2] = a_mass[2]*center_mass*center_mass+a_mass[3]*d0*d0*2.0;
#else

#endif

  if(sys_num == 0){
    w_site = 1;
    nden = mass_den3(1,-1,1,-1,0,temp);
    tscale = 1. / 3. /((double)n1 - 1);
    w_num = 0;
    w_num3 = 0;
    s_num = n1;
    s_num3 = n3;

    vmax = VMAX;
    oalpha = 6;
    if(np >= 8){
      vmax = 462;
      oalpha = 8.6;
    }
    if(np == 7){
      vmax = 462;
      oalpha = 8.2;
    }
    if(np == 6){
      vmax = 40;
      oalpha = 4.0;
      vmax = 462;
      oalpha = 7.6;
    }
    if(np == 5){
      vmax = 46;
      oalpha = 4.4;
      vmax = 462;
      oalpha = 6.9;
    }
    if(np == 4){
      vmax = 101;
      oalpha = 5.2;
    }
    if(np == 3){
      vmax = 309;
      oalpha = 5.6;
    }

  } else if(sys_num >= 1){

    vmax = VMAX;
    oalpha = 1.5;
    nden = nden_set(temp-273)/((a_mass[2]+a_mass[3]*2)*mass)*1e-27;
    tscale = 1. / 3. /((double)n1/w_site*2 - 1);

    w_num = n1/w_site;
    w_num3= w_num*3;
    s_num = 0;
    s_num3= 0;

    if(sys_num == 4){
      w_num = np*np*np*8-nn*2;
      w_num3= w_num*3;
      s_num = nn*2;
      s_num3= s_num*3;
      printf("np %d nn %d w_num %d s_num %d n1 %d\n",np,nn,w_num,s_num,n1);
    }
    if(sys_num == 5){
      w_num = nw;
      w_num3= w_num*3;
      s_num = np*np*np*8-nw;
      s_num3= s_num*3;
      nden = mass_den3(1,-1,1,-1,0,temp);
      printf("np %d nn %d w_num %d s_num %d n1 %d\n",np,nn,w_num,s_num,n1);
    }
    if(sys_num == 6){
      if(nw > 0)
        w_num = np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8;
      else
        w_num = 0;
      w_num3= w_num*3;
      s_num = np*np*np*8-w_num;
      s_num3= s_num*3;
      nden = mass_den3(1,-1,1,-1,0,temp);
      printf("np %d nw %d w_num %d s_num %d n1 %d\n",np,nw,w_num,s_num,n1);
    }
  }
  ws_num = w_num+s_num;
  ws_num3= ws_num*3;

  tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);
}
void keyboard(unsigned char key, int x, int y)
{
  int i,j,k;
  int i0,i1,i2;
  double d0,d1,d2,d3,d4,d5;
  double ang0,ang1,ang2,ang3;

  int c;
  int cf[6];
  double cfr[6];
  double cp[18];

  double l = 2.285852;

  l /= 2;

  if(key == '?'){
    printf("!   : initilize\n");
    printf("q,Q : quit program\n");
    printf("i,I : (print information of postion and angle)\n");
    printf("p,P : (set cliping area)\n");
    printf("a,A : auto mode on/off\n");
    printf("k,K : simulatin cell display on/off\n");
    printf("W   : start output bmp file\n");
    printf("c,C : change backgrand color\n");
    printf("r   : make radius small\n");
    printf("R   : make radius large\n");
    printf("d   : make ditail donw\n");
    printf("D   : make ditail up\n");
    printf("v,V : varbose on/off\n");
    printf("s   : md_step--\n");
    printf("S   : md_step++\n");
    printf("t,T : temp += 100\n");
    printf("g,G : temp -= 100\n");
    printf("y,Y : temp += 10\n");
    printf("h,H : temp -= 10\n");
    printf("z,Z : stop/restart\n");
    printf("0-9 : chage particle number\n");
    printf("n   : create a new positive ion\n");
    printf("m   : create a new negative ion\n");
    printf("N   : create 4 new ions\n");
    printf("M   : create 27 new ions\n");
    printf("SP  : shoot new ion(s)\n");
  }

  if(key == ' ' && c_flg != C_STEP && start_vl <= 0){
    grape_flg = (grape_flg == 0 ? 1:0);
  }

  if(key == '!'){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

  if((key >= '1' && key <= '9') && c_flg == 0){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    np = key-'0';
    npx = np;
    npy = np;
    npz = np;
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

  if(key == 'q' || key == 'Q'){
    if(kflg == 1)
      fclose(fp);
    exit(0);
  }
#ifdef INFO
  if(key == 'i' || key == 'I'){
    printf("(%f,%f,%f)\n"
	   ,trans[0]-trans0[0],trans[1]-trans0[1],trans[2]-trans0[2]);
    printf("(");
    for(i = 0; i < 4; i++){
      for(j = 0; j < 4; j++){
	if(i == 0 && j == 0)
	  printf("%f",m_matrix[i*4+j]-matrix0[i*4+j]);
	else
	  printf(",%f",m_matrix[i*4+j]-matrix0[i*4+j]);
      }
    }
    printf(")\n");
    for(i = 0; i < 3; i++){
      trans0[i] = trans[i];
    }
    for(i = 0; i < 16; i++)
      matrix0[i] = m_matrix[i];
  }
#endif
  if(key == 'f' || key == 'F') bond_flg = bond_flg == 0 ? 1:0;
#ifdef GL_ON
  if(key == 'p' || key == 'P'){
    //if(clip_flg == 4) clip_flg = 0; else clip_flg++;
  }
  if(key == 'a' || key == 'A') auto_flg = auto_flg == 0 ? 1:0;
  //if(key == 'k' || key == 'K') kabe_flg = kabe_flg == 0 ? 1:0;
  //if(key == 'W') save_flg = save_flg == 0 ? 1:0;
  if(key == 'c' || key == 'C') clear_color = clear_color == 0 ? 1:0;
  if(key == 'r' && (int)(radius*10+.5) > 1) {radius -= .1;}
  if(key == 'R') {radius += .1; }
  if(key == 'd' && ditail > 5)  {ditail -= 1; }
  if(key == 'D' && ditail < 20) {ditail += 1; }
  if(key == 'v' || key == 'V'){
    if(vflg != 0)
      vflg--;
    else
#ifdef LAP_TIME
      vflg = 3;
#else

#endif
  }
#endif
#if defined(MDGRAPE3) || defined(VTGRAPE)
  if(key == 's' && md_step > 10){ md_step -= 10; md_stepf = 10;}
  if(key == 'S'){ md_step += 10; md_stepf = 10;}
#else


#endif
  if(key == 't' || key == 'T'){
    temp += 100;
    rtemp = temp / epsv * kb;
  }
  if(key == 'g' || key == 'G'){
    if(temp > 100){
      temp -= 100;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'h' || key == 'H'){
    if(temp > 10){
      temp -= 10;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'y' || key == 'Y'){
    temp += 10;
    rtemp = temp / epsv * kb;
  }
  if(key == 'z' || key == 'Z'){
    run_flg *= -1;
#ifdef GL_ON
    if(sc_flg == 0){
      if(run_flg == 1);
        //glutIdleFunc(md_run);
      else
        glutIdleFunc(NULL);
    } else {
      if(run_flg == 1)
        glutIdleFunc(NULL);
      else
	glutIdleFunc(NULL);
    }
#endif
  }


  if(key >= '0' && key <= '9' && c_flg == C_STEP){
    start_vl = .3*(key-'0'+1)/10*delt/2e-15;
    velp_flg = 1;
  }


  if((key == 'N' || key == 'M' || key == 'B' ||
      key == 'n' || key == 'm' || key == 'b') && c_flg == 0){

    c_flg = 1;
    w_add = s_add = 0;

    if(key == 'b' || key == 'B'){
      c_num = w_site;
      w_add = 1;
    } else if(key == 'N'){
      c_num = 4;
      s_add = 4;
      r = 3;
    } else if(key == 'M'){
      c_num = 27;
      s_add = 27;
      r = 9;
    } else {
      c_num = 1;
      s_add = 1;
      r = 1;
    }

    d0 = (i_matrix[0]*(-trans[0])+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d1 = (i_matrix[1]*(-trans[0])+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d2 = (i_matrix[2]*(-trans[0])+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d0 += sideh[0];
    d1 += sideh[1];
    d2 += sideh[2];

    d3 = (i_matrix[0]*(-trans[0]+eye_len-10)+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d4 = (i_matrix[1]*(-trans[0]+eye_len-10)+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d5 = (i_matrix[2]*(-trans[0]+eye_len-10)+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d3 += sideh[0];
    d4 += sideh[1];
    d5 += sideh[2];
    t_cd[0] = d3-d0;
    t_cd[1] = d4-d1;
    t_cd[2] = d5-d2;

    if(d3 < side[0]-r && d3 > r &&
       d4 < side[1]-r && d4 > r &&
       d5 < side[2]-r && d5 > r){
    } else {

      for(i = 0; i < 6; i++){
        cf[i] = 0;
        cfr[i] = -1;
      }

      i0 = 0;      /* x > */
      cp[i0]   = side[0]-r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < x */
      cp[i0]   = r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[1]-r) cf[i0/3] = 1;

      i0 += 3;    /* y > */
      cp[i0+1]   = side[1]-r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < y */
      cp[i0+1]   = r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;

      i0 += 3;   /* z > */
      cp[i0+2]   = side[2]-r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;
      i0 += 3;
      cp[i0+2]   = r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;

      for(i = 0; i < 6; i++){
        if(cf[i] == 1){
          cfr[i] = sqrt((cp[i*3]  -d3)*(cp[i*3]  -d3)+
                        (cp[i*3+1]-d4)*(cp[i*3+1]-d4)+
                        (cp[i*3+2]-d5)*(cp[i*3+2]-d5));
        }
      }
      d0 = 10000;
      c = -1;
      for(i = 0; i < 6; i++){
        if(cf[i] == 1 && d0 > cfr[i]){
          d0 = cfr[i];
          c = i;
        }
      }
      if(c == -1) c_num = 0;
      c *= 3;
      d3 = cp[c];
      d4 = cp[c+1];
      d5 = cp[c+2];
    }

    if(key == 'b' || key == 'B'){
      for(i = 0; i < c_num*3/w_site; i+=3*w_site){
        cd[n3+i]   = d3;
        cd[n3+i+1] = d4;
        cd[n3+i+2] = d5;
        atype[n1+i/3] = 2;

        ang0 = PI/2 * (key == 'b' ? 1:-1);
        ang1 = 0;
        ang2 = 0;
        ang[w_num*4  ] = sin(ang1/2)*sin((ang2-ang0)/2);
        ang[w_num*4+1] = sin(ang1/2)*cos((ang2-ang0)/2);
        ang[w_num*4+2] = cos(ang1/2)*sin((ang2+ang0)/2);
        ang[w_num*4+3] = cos(ang1/2)*cos((ang2+ang0)/2);

        w_index[w_num] = n3+i;
        w_info[n1] =  (n3+i)/3+1;
        for(k = 1; k < w_site; k++){
          w_info[n1+k] =  (n3+i)/3;
          atype[n1+i/3+k] = 2+k;
        }
        ang0 = ang[w_num*4  ];
        ang1 = ang[w_num*4+1];
        ang2 = ang[w_num*4+2];
        ang3 = ang[w_num*4+3];
        for(k = 0; k < w_site-1; k++){
          d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
              +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
              +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
          d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
              +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
              +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
          d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
              +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
              +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
          cd[(k+1)*3  +n3+i] = cd[n3+i  ] + d0;
          cd[(k+1)*3+1+n3+i] = cd[n3+i+1] + d1;
          cd[(k+1)*3+2+n3+i] = cd[n3+i+2] + d2;
        }
      }
    } else {


      if(c_num == 1){
        cd[n3]   = d3;
        cd[n3+1] = d4;
        cd[n3+2] = d5;
        atype[n1] = ((key == 'n') ? 0:1);
      } else if(c_num == 4){
        i = 0;
        l *= 1.1;
        for(i0 = 0; i0 < 2; i0++){
          for(i1 = 0; i1 < 2; i1++){
            d0 = i_matrix[4]*(l*(i0*2-1))
                +i_matrix[8]*(l*(i1*2-1));
            d1 = i_matrix[5]*(l*(i0*2-1))
                +i_matrix[9]*(l*(i1*2-1));
            d2 = i_matrix[6]*(l*(i0*2-1))
                +i_matrix[10]*(l*(i1*2-1));
            cd[n3+i]   = d0+d3;
            cd[n3+i+1] = d1+d4;
            cd[n3+i+2] = d2+d5;
            atype[n1+i/3] = (i0+i1) %2;
            i += 3;
          }
        }
      } else if(c_num == 27){
        i = 0;
        l *= 1.2;
        for(i0 = 0; i0 < 3; i0++){
          for(i1 = 0; i1 < 3; i1++){
            for(i2 = 0; i2 < 3; i2++){
              d0 = i_matrix[0]*(l*2*(i0-1))
                  +i_matrix[4]*(l*2*(i1-1))
                  +i_matrix[8]*(l*2*(i2-1));
              d1 = i_matrix[1]*(l*2*(i0-1))
                  +i_matrix[5]*(l*2*(i1-1))
                  +i_matrix[9]*(l*2*(i2-1));
              d2 = i_matrix[2]*(l*2*(i0-1))
                  +i_matrix[6]*(l*2*(i1-1))
                  +i_matrix[10]*(l*2*(i2-1));
              cd[n3+i]   = d0+d3;
              cd[n3+i+1] = d1+d4;
              cd[n3+i+2] = d2+d5;
              atype[n1+i/3] = (i0+i1+i2) %2;
              i += 3;
            }
          }
        }
      }
    }

  }
  if(key == ' ' && c_flg == C_STEP && start_vl > 0){

    w_num += w_add;
    w_num3 = w_num*3;
    s_num += s_add;
    s_num3 = s_num*3;
    ws_num = w_num + s_num;
    ws_num3= ws_num*3;

    n1 = s_num + w_num*w_site;
    n3 = n1*3;
    tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);

    r = sqrt(t_cd[0]*t_cd[0] + t_cd[1]*t_cd[1] + t_cd[2]*t_cd[2]);
    d3 = 0;
    for(i = 0; i < c_num*3; i+=3){
      vl[n3-3-i]   = -t_cd[0]/r *start_vl;
      vl[n3-3-i+1] = -t_cd[1]/r *start_vl;
      vl[n3-3-i+2] = -t_cd[2]/r *start_vl;
      d3 += (vl[n3-3-i]*vl[n3-3-i]+vl[n3-2-i]*vl[n3-2-i]+vl[n3-1-i]*vl[n3-1-i])
        *a_mass[atype_mat[atype[(n3-3-i)/3]]];
    }
    d3 *= tscale/hsq;
    rtemp += d3;
    temp = rtemp * epsv /kb;

    c_flg = 0;
    c_num = 0;
    velp_flg = 0;
    start_vl = -1;

   }

#ifdef GL_ON
  if(sc_flg != 1)
    glutPostRedisplay();
#endif
}

void md_run()
{
  int i,j,k,c;
  int i0,i1,i2,i3,i4,i5;
  double d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12;

  double agv0,agv1,agv2;
  double ang0,ang1,ang2,ang3;
  int md_loop;

  double dphir;


  for(md_loop = 0; md_loop < md_step; md_loop++){

    m_clock++;

#ifdef C_MASS
    for(i = 0; i < n3; i += 3){
      if(atype[i/3] == 2){
	ang0 = ang[w_rindex[i/3]*4];
	ang1 = ang[w_rindex[i/3]*4+1];
	ang2 = ang[w_rindex[i/3]*4+2];
	ang3 = ang[w_rindex[i/3]*4+3];

	d0 = center_mass*(-2)*(ang0*ang1+ang2*ang3);
	d1 = center_mass*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3);
	d2 = center_mass*  2 *(ang1*ang3-ang0*ang2);

	cd[i]   -= d0;
	cd[i+1] -= d1;
	cd[i+2] -= d2;

      }
    }
#endif

      for(i = 0; i < n3; i++){ /* update coordinations */
        if(atype[i/3] <= 2 && atype[i/3] != 8){
#if T_CONST == 1
          vl[i] = (vl[i]*(1-xs)+fc[i])/(1+xs);
#else

#endif
          cd[i]   += vl[i];
        }
      }
      for(i = 0; i < w_num3; i++){ /* half update anglusr velocity  */
        agv[i]   = agvh[i]   + trq[i]  *0.5;
      }
      for(i=0, j=0; i < w_num3; i+=3, j+=4){ /* coordination transformation */
        agv0 = agv[i  ];
        agv1 = agv[i+1];
        agv2 = agv[i+2];

        ang0 = ang[j  ];
        ang1 = ang[j+1];
        ang2 = ang[j+2];
        ang3 = ang[j+3];
        agvp[i  ] =(agv0     *(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
                    +agv1*  2 *( ang2*ang3-ang0*ang1)
                    +agv2*  2 *( ang1*ang2+ang0*ang3)) /moi[0];
        agvp[i+1] =(agv0*(-2)*( ang0*ang1+ang2*ang3)
                    +agv1     *( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
                    +agv2*  2 *( ang1*ang3-ang0*ang2)) /moi[1];
        agvp[i+2] =(agv0*  2 *( ang1*ang2-ang0*ang3)
                    +agv1*(-2)*( ang0*ang2+ang1*ang3)
                    +agv2     *(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3))/moi[2];
      }
      for(i=0, j=0; i < w_num3; i+=3, j+=4){ /* half update of angle */
        ang0 = ang[j  ];
        ang1 = ang[j+1];
        ang2 = ang[j+2];
        ang3 = ang[j+3];
        angh[j]   = ang0+(-ang2*agvp[i]-ang3*agvp[i+1]+ang1*agvp[i+2])*0.25;
        angh[j+1] = ang1+( ang3*agvp[i]-ang2*agvp[i+1]-ang0*agvp[i+2])*0.25;
        angh[j+2] = ang2+( ang0*agvp[i]+ang1*agvp[i+1]+ang3*agvp[i+2])*0.25;
        angh[j+3] = ang3+(-ang1*agvp[i]+ang0*agvp[i+1]-ang2*agvp[i+2])*0.25;
      }

      for(i = 0; i < w_num*4; i += 4){ /* correction of angle */
        ang0 = angh[i  ];
        ang1 = angh[i+1];
        ang2 = angh[i+2];
        ang3 = angh[i+3];
        angh[i]  /= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        angh[i+1]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        angh[i+2]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        angh[i+3]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
      }

      for(i = 0; i < w_num3; i++){ /* update of anglular velocity */
#if T_CONST == 1
        agvh[i] = (agvh[i]*(1-xs)+trq[i])/(1+xs);
#else

#endif
      }

      for(i=0, j=0; i < w_num3; i+=3, j+=4){ /* coordination transformation */
        agv0 = agvh[i  ];
        agv1 = agvh[i+1];
        agv2 = agvh[i+2];
        ang0 = angh[j  ];
        ang1 = angh[j+1];
        ang2 = angh[j+2];
        ang3 = angh[j+3];
        agvph[i  ]=(agv0     *(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
                    +agv1*  2 *( ang2*ang3-ang0*ang1)
                    +agv2*  2 *( ang1*ang2+ang0*ang3)) /moi[0];
        agvph[i+1]=(agv0*(-2)*( ang0*ang1+ang2*ang3)
                    +agv1     *( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
                    +agv2*  2 *( ang1*ang3-ang0*ang2)) /moi[1];
        agvph[i+2]=(agv0*  2 *( ang1*ang2-ang0*ang3)
                    +agv1*(-2)*( ang0*ang2+ang1*ang3)
                    +agv2     *(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3))/moi[2];
      }
      for(i=0, j=0; i < w_num3; i+=3, j+=4){ /* update of angle */
        ang0 = angh[j  ];
        ang1 = angh[j+1];
        ang2 = angh[j+2];
        ang3 = angh[j+3];
        ang[j]  += (-ang2*agvph[i] - ang3*agvph[i+1] + ang1*agvph[i+2])*0.5;
        ang[j+1]+= ( ang3*agvph[i] - ang2*agvph[i+1] - ang0*agvph[i+2])*0.5;
        ang[j+2]+= ( ang0*agvph[i] + ang1*agvph[i+1] + ang3*agvph[i+2])*0.5;
        ang[j+3]+= (-ang1*agvph[i] + ang0*agvph[i+1] - ang2*agvph[i+2])*0.5;
      }
      for(i = 0; i < w_num*4; i += 4){ /* correction of angle */
        ang0 = ang[i  ];
        ang1 = ang[i+1];
        ang2 = ang[i+2];
        ang3 = ang[i+3];
        ang[i]  /= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        ang[i+1]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        ang[i+2]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
        ang[i+3]/= sqrt(ang0*ang0+ang1*ang1+ang2*ang2+ang3*ang3);
      }

#ifdef C_MASS
    for(i = 0; i < n3; i += 3){
      if(atype[i/3] == 2){
	ang0 = ang[w_rindex[i/3]*4];
	ang1 = ang[w_rindex[i/3]*4+1];
	ang2 = ang[w_rindex[i/3]*4+2];
	ang3 = ang[w_rindex[i/3]*4+3];

	d0 = center_mass*(-2)*(ang0*ang1+ang2*ang3);
	d1 = center_mass*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3);
	d2 = center_mass*  2 *(ang1*ang3-ang0*ang2);

	cd[i]   += d0;
	cd[i+1] += d1;
	cd[i+2] += d2;

      }
    }
#endif

      for(i0 = 0; i0 < w_num; i0++){
        i = w_index[i0];
        j = i0*4;
        c = w_info[i/3]*3;
        ang0 = ang[j  ];
        ang1 = ang[j+1];
        ang2 = ang[j+2];
        ang3 = ang[j+3];
        /*      printf("%d % f % f % f % f\n",j/4,ang0,ang1,ang2,ang3);*/
        for(k = 0; k < w_site-1; k++){
          d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
              +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
              +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
          d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
              +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
              +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
          d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
              +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
              +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
          cd[k*3+c  ] = cd[i  ] + d0;
          cd[k*3+c+1] = cd[i+1] + d1;
          cd[k*3+c+2] = cd[i+2] + d2;
        }
      }
#if ZERO_P == 0

#else
      for(i = 0; i < n3; i += 3){
        if(atype[i/3] <= 2){
#if 1
          if (cd[i]   < 0 || cd[i]   > side[0]) vl[i] *= -1;
          if (cd[i+1] < 0 || cd[i+1] > side[1]) vl[i+1] *= -1;
          if (cd[i+2] < 0 || cd[i+2] > side[2]) vl[i+2] *= -1;
#else

#endif
        }
	else{
	  printf("atye[%d] is %d\n",i/3,atype[i/3]);
	}
      }
#endif
      /*****************  calculation of force  *********************/

#ifdef LAP_TIME
#if defined(_WIN32) && !defined(__CYGWIN__)
    md_time0 = (double)timeGetTime()/1000.;
#elif defined(MAC)
    md_time0 = (double)clock()/60.;
#else
    gettimeofday(&time_v,NULL);
    md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif
#endif

      for(i = 0;i < 2; i++){
        phi[i] = 0;
      }
      phir = 0;
      for(i = 0;i < n3; i++){
        fc[i] = 0;
        iphi[i] = 0;
      }
      vir = 0;

      /*#if MDM == 0*/
      for(i = 0;i < s_num3; i++){
        fc[i] = 0;
      }

      if(grape_flg == 1){

    	  int rc;
    	  pthread_t threads[NUMT];
    	  pthread_attr_t attr[NUMT];
    	  size_t stacksize;
    	  void *status;

    	  for(i=0;i<NUMT;i++){
			  pthread_attr_init(&attr[i]);
			  pthread_attr_getstacksize(&attr[i],&stacksize);
			  pthread_attr_setdetachstate(&attr[i], PTHREAD_CREATE_JOINABLE);

			  //printf("Default stack size = %li\n", stacksize);
			  stacksize = sizeof(thread_data_array) * NUMT;
			  //printf("Amount of stack needed per thread = %li\n",stacksize);
			  pthread_attr_setstacksize (&attr[i], stacksize);
			  //printf("Creating threads with stack size = %li bytes\n",stacksize);

			  thread_data_array[i].thread_id = i;
			  thread_data_array[i].n1 = n1;
			  thread_data_array[i].range = n1/NUMT ;
			  thread_data_array[i].nig_num = nig_num;
			  thread_data_array[i].n3 = n3;
			  thread_data_array[i].atype_mat = atype_mat;
			  thread_data_array[i].atype=atype;
			  thread_data_array[i].cd=cd;
			  thread_data_array[i].r=r;
			  thread_data_array[i].rd=rd;
			  thread_data_array[i].inr=inr;
			  thread_data_array[i].phir=phir;
			  thread_data_array[i].pb=pb;
			  thread_data_array[i].pol=(double *)pol;
			  thread_data_array[i].sigm=(double *)sigm;
			  thread_data_array[i].pc=(double *)pc;
			  thread_data_array[i].pd=(double *)pd;
			  thread_data_array[i].ipotro=(double *)ipotro;
			  thread_data_array[i].zz=(double *)zz;
			  thread_data_array[i].vir=vir;
			  thread_data_array[i].fc=fc;

			  rc = pthread_create(&threads[i], &attr[i], CPU_force, (void *) &thread_data_array[i]);
					if (rc){
					   printf("ERROR; return code from pthread_create() is %d\n", rc);
					   exit(0);
					}

			  pthread_attr_destroy(&attr[i]);
			  rc = pthread_join(threads[i],&status);
			  if(rc){
					printf("ERROR; return code from pthread_join()is %d\n", rc);
					exit(0);
			  }
    	  }


    	  fc = thread_data_array[1].fc;

      }
      else {

    	for (i = 0; i < n1; i++)nig_num[i] = 0;

		for (i = 0; i < n3; i += 3) {
			i0 = atype_mat[atype[i / 3]];
			for (j = i + 3; j < n3; j += 3) {
				d0 = cd[i] - cd[j];
				d1 = cd[i + 1] - cd[j + 1];
				d2 = cd[i + 2] - cd[j + 2];

				rd = d0 * d0 + d1 * d1 + d2 * d2;
				r = sqrt(rd);
				inr = 1. / r;

				i1 = atype_mat[atype[j / 3]];
				d7 = phir;

				if (i0 < 2 && i1 < 2) {
						d3 = pb * pol[i0][i1]
								* exp((sigm[i0][i1] - r) * ipotro[i0][i1]);

						dphir = (d3 * ipotro[i0][i1] * inr
								- 6 * pc[i0][i1] * pow(inr, 8)
								- 8 * pd[i0][i1] * pow(inr, 10)
								+ inr * inr * inr * zz[i0][i1]);
					}

				vir -= rd * dphir;

				d3 = d0 * dphir;
				d4 = d1 * dphir;
				d5 = d2 * dphir;

				fc[i] += d3;
				fc[i + 1] += d4;
				fc[i + 2] += d5;
				fc[j] -= d3;
				fc[j + 1] -= d4;
				fc[j + 2] -= d5;

			}
		}

      }

#if 1
      /*#define CHECK*/

#if defined(MDGRAPE3) || defined(VTGRAPE)

      if(bond_flg == 1 && grape_flg==0){
#else

#endif

	for(i = 0; i < n1; i++){
	  i0 = 0;
	  for(j = 0; j < nig_num[i]; j++){
	    d0 = cd[nig_data[i*6+j]*3]  -cd[i*3];
	    d1 = cd[nig_data[i*6+j]*3+1]-cd[i*3+1];
	    d2 = cd[nig_data[i*6+j]*3+2]-cd[i*3+2];
	    for(k = j+1; k < nig_num[i]; k++){
	      d3 = cd[nig_data[i*6+k]*3]  -cd[i*3];
	      d4 = cd[nig_data[i*6+k]*3+1]-cd[i*3+1];
	      d5 = cd[nig_data[i*6+k]*3+2]-cd[i*3+2];
	      d6 = acos(fabs((d0*d3 + d1*d4 + d2*d5)
			     /sqrt(d0*d0+d1*d1+d2*d2)
			     /sqrt(d3*d3+d4*d4+d5*d5)))/PI*180;
	      if(d6 < min_angle || d6 > max_angle) i0++;
	    }
	  }
	  if(i0 < 4 || i0 != nig_num[i]*(nig_num[i]-1)/2){

	    nig_num[i] = 0;
	  } else {

	  }
	}

	for(i = 0; i < n1; i++){
	  if(nig_num[i] < 6){
	    i0 = 0;
	    for(j = 0; j < nig_num[i]; j++){
	      if(nig_num[nig_data[i*6+j]] >= 4) i0++;
	    }
	    if(i0 == 0) nig_num[i] = 0;

	  }
	}

      }

#endif




#ifdef LAP_TIME
#if defined(_WIN32) && !defined(__CYGWIN__)
    md_time = (double)timeGetTime()/1000.;
#elif defined(MAC)
    md_time = (double)clock()/60.;
#else

    //for (i=90;i<100;i++)printf("fc[%d]=%f\n",i,fc[i]);
    gettimeofday(&time_v,NULL);
    md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif
#endif

#ifdef C_MASS
    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0];
      j = i0*4;
      c = w_info[i/3]*3;
      ang0 = ang[j  ];
      ang1 = ang[j+1];
      ang2 = ang[j+2];
      ang3 = ang[j+3];

      trq[i0*3] = trq[i0*3+1] = trq[i0*3+2] = 0;

      d10 = (center_mass*(-2)*(ang0*ang1+ang2*ang3));
      d11 = (center_mass*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3));
      d12 = (center_mass*  2 *(ang1*ang3-ang0*ang2));
      trq[i0*3  ] += d11*fc[i+2] - d12*fc[i+1];
      trq[i0*3+1] += d12*fc[i  ] - d10*fc[i+2];
      trq[i0*3+2] += d10*fc[i+1] - d11*fc[i  ];
      for(k = 0; k < w_site-1; k++){
	d0 = (m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)+
	      (m_cdy[k]+center_mass)*(-2)*(ang0*ang1+ang2*ang3)+
	      m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3));
	d1 = (m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)+
	      (m_cdy[k]+center_mass)*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)+
	      m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3));
	d2 = (m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)+
	      (m_cdy[k]+center_mass)*  2 *(ang1*ang3-ang0*ang2)+
	      m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3));
	trq[i0*3  ] += d1*fc[k*3+c+2] - d2*fc[k*3+c+1];
	trq[i0*3+1] += d2*fc[k*3+c  ] - d0*fc[k*3+c+2];
	trq[i0*3+2] += d0*fc[k*3+c+1] - d1*fc[k*3+c  ];
      }
    }
#else

#endif

      /*
      for(i = 0; i < w_num3; i += 3){
        printf("%3d % e % e % e\n",i/3,trq[i],trq[i+1],trq[i+2]);
      }
      */

      for(i = 0; i < n3; i++){
        if(atype[i/3] == 2)
          fc[i] *= hsq/(a_mass[2]+2*a_mass[3]);
        else if(atype[i/3] == 0 || atype[i/3] == 1)
          fc[i] *= hsq/a_mass[atype_mat[atype[i/3]]];
      }

      for(i = 0; i < w_num3; i++)
        trq[i] *= hsq;

      ekin1 = 0;
      ekin2 = 0;
      for(i = 0; i < n3; i += 3){
        ekin1 += (vl[i  ]*vl[i  ] +
                  vl[i+1]*vl[i+1] +
                  vl[i+2]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
      }
      for(i = 0; i < w_num3; i += 3){
        ekin2 += (moi[0]*agvph[i  ]*agvph[i  ]+
                  moi[1]*agvph[i+1]*agvph[i+1]+
                  moi[2]*agvph[i+2]*agvph[i+2]);
      }

      ekin1 /= hsq;
      ekin2 /= hsq;

      ekin = ekin1 + ekin2;


      mtemp = tscale * ekin;
#if 0
      printf("mtemp=%e temp=%e\n",mtemp,temp);
#endif
      mpres  = nden / 3. * (ekin - vir) / (s_num+w_num);
#if T_CONST == 1
      xs += (mtemp - rtemp)/ lq * hsq *.5;
#endif
/*
      printf("%4d %f %f %f %f %f %f\n",m_clock,mtemp*epsv/kb
             ,(ekin/2.+phir)*erdp/(double)(s_num/2+w_num)
             ,ekin/2.*erdp/(double)(s_num/2+w_num)
             ,phir*erdp/(double)(s_num/2+w_num)
             ,ekin1/(ekin1+ekin2), ekin2/(ekin1+ekin2));
*/
  }

#if defined(SUBWIN) && defined(GL_ON)
  if((b_clock % 10) == 0){
    if(p_count < DATA_NUM){
      temp_data[p_count] = (int)(mtemp*epsv/kb);
      p_count++;
    } else {
      for(i = 0; i < DATA_NUM-1; i++){
	temp_data[i] = temp_data[i+1];
      }
      temp_data[DATA_NUM-1] = (int)(mtemp*epsv/kb);
    }

    temp_max = 0;
    for(i = 0; i < DATA_NUM; i++){
      if(temp_data[i] > temp_max) temp_max=temp_data[i];
    }
    if(temp_ymax < temp_max){ temp_ymax = temp_max*1.5;}
    if(temp_ymax > temp_max*2 && temp_ymax/2 > 2000){ temp_ymax /= 2;}

  }
#endif

#ifdef GL_ON
  if(sc_flg != 1)
    glutPostRedisplay();
#endif

  if(auto_flg == 1){
#ifdef GL_ON
    mouse_l = tt[b_clock].mouse[0];
    mouse_m = tt[b_clock].mouse[1];
    mouse_r = tt[b_clock].mouse[2];
    for(i = 0; i < 3; i++){
      trans[i] += tt[b_clock].move[i];
      angle[i] = tt[b_clock].rot[i];
    }
#endif
    if(sc_flg != 1){
      if(tt[b_clock].command != 0){
    	  //keyboard(tt[b_clock].command,0,0);
      }

    }
    temp += tt[b_clock].temp;
    rtemp = temp / epsv * kb;
    for(i = 0; i < 16; i++)
      m_matrix[i] += tt[b_clock].matrix[i];

  }

  b_clock++;
}


void set_cd(int ini_m2)
{
  int i,j,k,c;
  int i0,i1,i2,i3,i4,i5,i10;
  double d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12;
  double ang0,ang1,ang2,ang3;

  if(ini_m2 == 0){
    temp = ini_temp;
    rtemp = temp / epsv * kb;
    if(sys_num == 0){
      n1 = np*np*np*8;
    } else if(sys_num == 1 || sys_num == 10){
      n1 = np*np*np*4*w_site;
    } else if(sys_num == 2){
#if ZERO_P == 1
      n1 = (npx*npy*npz*8+npy*npz*4+(npy-1)*npz*8)*w_site;
#else
      n1 = np*np*np*8*w_site;
#endif
    } else if(sys_num == 3){
      n1 = np*np*np*8*w_site;
    } else if(sys_num == 4){
      if(np*np*np*8-nn*2 >= 0)
        n1 = (np*np*np*8-nn*2)*w_site+nn*2;
      else {
        printf("nn is too large!\n");
        exit(0);
      }
    } else if(sys_num == 5){
      if(np*np*np*8-nw >= 0)
        n1 = (np*np*np*8-nw)+nw*w_site;
      else {
        printf("nw is too large!\n");
        exit(0);
      }
    } else if(sys_num == 6){
      if(np-nw >= 0){
        if(nw > 0)
          n1 = np*np*np*8+(np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8)*(w_site-1);
        else
          n1 = np*np*np*8;
      } else {
        printf("nw is too large!\n");
        exit(0);
      }
    }
    n2 = n1 * 2;
    n3 = n1 * 3;

    if(sys_num == 0){
      w_num = 0;
      w_num3 = 0;
      s_num = n1;
      s_num3 = n3;
    } else if(sys_num >= 1){
      w_num = n1/w_site;
      w_num3= w_num*3;
      s_num = 0;
      s_num3= 0;
      if(sys_num == 4){
        w_num = np*np*np*8-nn*2;
        w_num3= w_num*3;
        s_num = nn*2;
        s_num3= s_num*3;
      }
      if(sys_num == 5){
        w_num = nw;
        w_num3= w_num*3;
        s_num = np*np*np*8-nw;
        s_num3= s_num*3;
      }
      if(sys_num == 6){
        if(nw > 0)
          w_num = np*np*np*8-(np-nw)*(np-nw)*(np-nw)*8;
        else
          w_num = 0;
        w_num3= w_num*3;
        s_num = np*np*np*8-w_num;
        s_num3= s_num*3;
      }
    }
    ws_num = w_num+s_num;
    ws_num3= ws_num*3;

    tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);
  }

  if(sys_num == 0){
    side[0] = pow(8 / nden, 1./3.) * np;
    side0 = side[0];
    side[1] = side[0];
    side[2] = side[0];
    fccset2(np,side[0],cd);      /* set fcc */
    for(i = 0; i < s_num3/2; i++)
      cd[i+s_num3/2] = cd[i];
    for(i = 0;i < s_num3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = 0; i < s_num/2; i++)
      atype[i] = 0;
    for(i = s_num/2; i < s_num; i++)
      atype[i] = 1;

  } else if (sys_num == 1){
    strcpy(keiname,"water");
    side[0] = pow(4./nden,1./3.)*npx;
    side[1] = pow(4./nden,1./3.)*npy;
    side[2] = pow(4./nden,1./3.)*npz;

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
      w_rindex[i] = i;
    }
    for(i = 0; i < w_num; i++){
      w_info[i] = i*(w_site-1)+w_num;
    }
    for(i = w_num; i < n1; i++){
      w_info[i] = (i-w_num)/(w_site-1);
    }
    /*
  for(i = 0; i < w_num3; i++)
    printf("%d %d\n",i,w_info[i]);
*/
    fccset_w(side);

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }
    /*
    for(i = 0; i < n3; i += 3)
      printf("%d %d %f %f %f\n",i/3,atype[i/3],cd[i],cd[i+1],cd[i+2]);
    */
    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }
  } else if (sys_num == 10){
    strcpy(keiname,"water");
    side[0] = pow(4./nden,1./3.)*npx;
    side[1] = pow(4./nden,1./3.)*npy;
    side[2] = pow(4./nden,1./3.)*npz;

#if 0

#else

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
    }
    fccset_w(side);
    s_num = 0;
    for(i = 0; i < w_num3; i += 3){
      if((cd[i]   > side[0]/np/2*(int)((np-2)/3*2) && cd[i]   < side[0]-side[0]/np/2*(int)((np-2)/3*2) &&
	  cd[i+1] > side[0]/np/2*(int)((np-2)/3*2) && cd[i+1] < side[0]-side[0]/np/2*(int)((np-2)/3*2)
	  ) &&
	 ((cd[i+2] < side[0]/2+side[0]/np/2*3 &&
	   cd[i+2] > side[0]/2+side[0]/np-side[0]/np/2) ||
	  (cd[i+2] > side[0]/2-side[0]/np/2*3 &&
	   cd[i+2] < side[0]/2-side[0]/np+side[0]/np/2))
	 ){
	atype[i/3] = 8;
	s_num++;
      } else {
	atype[i/3] = 2;
      }
    }
    s_num3 = s_num*3;
    w_num -= s_num;
    w_num3 = w_num*3;
    ws_num = w_num + s_num;
    ws_num3 = w_num3 + s_num3;

#ifdef GL_ON
    clip[0][0] = -1.0;
    clip[0][1] =  0.0;
    clip[0][2] =  0.0;
    clip[0][3] =  side[0]-side[0]/np/2*(int)((np-2)/3*2)-side[0]/2;
    clip[1][0] =  1.0;
    clip[1][1] =  0.0;
    clip[1][2] =  0.0;
    clip[1][3] =  clip[0][3];

    clip[2][0] =  0.0;
    clip[2][1] = -1.0;
    clip[2][2] =  0.0;
    clip[2][3] =  side[0]-side[0]/np/2*(int)((np-2)/3*2)-side[0]/2;
    clip[3][0] =  0.0;
    clip[3][1] =  1.0;
    clip[3][2] =  0.0;
    clip[3][3] =  clip[2][3];

    clip[4][0] =  0.0;
    clip[4][1] =  0.0;
    clip[4][2] = -1.0;
    clip[4][3] =  side[0]/2+side[0]/np-side[0]/np/2-side[0]/2;
    clip[5][0] =  0.0;
    clip[5][1] =  0.0;
    clip[5][2] =  1.0;
    clip[5][3] =  side[0]/2+side[0]/np/2*3-side[0]/2;
#endif

    for(i = s_num+w_num; i < n1; i += w_site-1){
      for(j = 0; j < w_site-1; j++)
        atype[i+j] = 3+j;
    }


    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }

    n1 = s_num + w_num*w_site;
    n2 = n1*2;
    n3 = n1*3;

#endif

    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }
  } else if(sys_num == 2){

    strcpy(keiname,"water");
    for(i = 0; i < w_num; i++)
      w_index[i] = i*3;

    for(i = 0; i < w_num; i++)
      w_info[i] = i*(w_site-1)+w_num;

    for(i = w_num; i < n1; i++)
      w_info[i] = (i-w_num)/(w_site-1);

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }

    ice_set(side);

#if ZERO_P == 1

    c = npx*npy*npz*8*3;
    for(i = 0; i < npx*npy*npz*8*3; i += 3){
      if(cd[i] < side[0]/npx/2){
        cd[c]   = cd[i]  +side[0];
        cd[c+1] = cd[i+1];
        cd[c+2] = cd[i+2];
        i0 = w_index[i/3]/3*4;
        ang[c/3*4]   = ang[i0];
        ang[c/3*4+1] = ang[i0+1];
        ang[c/3*4+2] = ang[i0+2];
        ang[c/3*4+3] = ang[i0+3];
        c += 3;
      }
    }
    c = (npx*npy*npz*8+npy*npz*4)*3;
    for(k = 0; k < npz ; k++){
      for(i = 0; i < npy-1 ; i++){
        for(j = 0; j < 4*3; j += 3){
          cd[c]   = cd[j]   - side[0]/npx;
          cd[c+1] = cd[j+1] + (side[1]/npy)*(j >= 6 ? 1:0) + side[1]/npy*i;
          cd[c+2] = cd[j+2] + side[2]/npz*k;
          i0 = w_index[j/3]/3*4;
          ang[c/3*4]   = ang[i0];
          ang[c/3*4+1] = ang[i0+1];
          ang[c/3*4+2] = ang[i0+2];
          ang[c/3*4+3] = ang[i0+3];
          c += 3;
        }
        for(j = 0; j < 4*3; j += 3){
          cd[c]   = cd[j]   + side[0];
          cd[c+1] = cd[j+1] + (side[1]/npy)*(j >= 6 ? 1:0) + side[1]/npy*i;
          cd[c+2] = cd[j+2] + side[2]/npz*k;
          i0 = w_index[j/3]/3*4;
          ang[c/3*4]   = ang[i0];
          ang[c/3*4+1] = ang[i0+1];
          ang[c/3*4+2] = ang[i0+2];
          ang[c/3*4+3] = ang[i0+3];
          c += 3;
        }
      }
    }

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0];
      j = i0*4;
      c = w_info[i/3]*3;
      ang0 = ang[j  ];
      ang1 = ang[j+1];
      ang2 = ang[j+2];
      ang3 = ang[j+3];
      for(k = 0; k < w_site-1; k++){
        d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
            +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
        d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
            +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
            +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
        d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
            +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
            +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
        cd[k*3+c  ] = cd[i  ] + d0;
        cd[k*3+c+1] = cd[i+1] + d1;
        cd[k*3+c+2] = cd[i+2] + d2;
      }
    }

    side[0] += side[0]/npx/2;

#ifdef GL_ON


    ini_flg = 1;
    mouse_l = 1;
    angle[1] = 90;

#endif

#endif

  } else if(sys_num == 3){
    strcpy(keiname,"water");
    side[0] = pow(8./nden,1./3.)*npx;
    side[1] = pow(8./nden,1./3.)*npy;
    side[2] = pow(8./nden,1./3.)*npz;

    for(i = 0; i < w_num; i++){
      w_index[i] = i*3;
    }
    for(i = 0; i < w_num; i++){
      w_info[i] = i*(w_site-1)+w_num;
    }
    for(i = w_num; i < n1; i++){
      w_info[i] = (i-w_num)/(w_site-1);
    }

    for(i0 = 0; i0 < w_num; i0++){
      i = w_index[i0]/3;
      atype[i] = 2;
      for(j = 0; j < w_site-1; j++)
        atype[w_info[i]+j] = j+3;
    }

    ice_set2(side);


  } else if(sys_num == 4){
    strcat(keiname,"-water");
    side[0] = pow(8./nden,1./3.)*npx;
    side[1] = pow(8./nden,1./3.)*npy;
    side[2] = pow(8./nden,1./3.)*npz;

    for(i = 0; i < s_num+w_num; i++)
      atype[i] = 2;
    for(i = s_num+w_num; i < n1; i += w_site-1){
      for(j = 0; j < w_site-1; j++)
        atype[i+j] = 3+j;
    }

    c = 0;
    while(c < s_num){
      i0 = (int)((double)rand()/(double)RAND_MAX*(s_num+w_num));
      if(atype[i0] == 2){
        if(c % 2 == 0) atype[i0] = 0; else atype[i0] = 1;
        c++;
      }
    }

    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }



  }else if (sys_num == 5){

    side[0] = pow(8 / nden, 1./3.) * np;
    side[1] = side[0];
    side[2] = side[0];

    for(i = 0; i < (n1-w_num*(w_site-1))/2; i++)
      atype[i] = 0;
    for(i = (n1-w_num*(w_site-1))/2; i < n1-w_num*(w_site-1); i++)
      atype[i] = 1;

    fccset2(np,side[0],cd);     /* set fcc */
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i++)
      cd[i+(n1-w_num*(w_site-1))*3/2] = cd[i];
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = (n1-w_num*(w_site-1))*3; i < n3; i++)
      cd[i] = 0;

    c = 0;
    while(c < w_num){
      i0 = (int)((double)rand()/(double)RAND_MAX*(s_num+w_num));
      if(atype[i0] == (c%2)){
        atype[i0] = 2;
        for(i = 0; i < w_site-1; i++)
          atype[s_num+w_num+c*(w_site-1)+i] = 3+i;
        c++;
      }
    }


    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++)
          w_info[c+w_num+s_num+j] = i;
        c += w_site-1;
      }

    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }


  }else if (sys_num == 6){

    side[0] = pow(8 / nden, 1./3.) * np;
    side[1] = side[0];
    side[2] = side[0];

    for(i = 0; i < (n1-w_num*(w_site-1))/2; i++)
      atype[i] = 0;
    for(i = (n1-w_num*(w_site-1))/2; i < n1-w_num*(w_site-1); i++)
      atype[i] = 1;

    fccset2(np,side[0],cd);     /* set fcc */
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i++)
      cd[i+(n1-w_num*(w_site-1))*3/2] = cd[i];
    for(i = 0; i < (n1-w_num*(w_site-1))*3/2; i += 3){
      cd[i] += side[0] / np / 2.;
      if(cd[i] < 0)       cd[i] += side[0];
      if(cd[i] > side[0]) cd[i] -= side[0];
    }
    for(i = (n1-w_num*(w_site-1))*3; i < n3; i++)
      cd[i] = 0;

    d0 = side[0]/(npx*2)*nw;
    d1 = side[1]/(npy*2)*nw;
    d2 = side[2]/(npz*2)*nw;

    for(i = 0; i < n3; i += 3){
      if(cd[i]   < d0 || cd[i]   > side[0]-d0) atype[i/3] = 2;
      if(cd[i+1] < d1 || cd[i+1] > side[1]-d1) atype[i/3] = 2;
      if(cd[i+2] < d2 || cd[i+2] > side[2]-d2) atype[i/3] = 2;

    }

    for(i = 0; i < n1; i++)
      w_info[i] = -1;
    c = 0;
    for(i = 0; i < s_num+w_num; i++)
      if(atype[i] == 2){
        w_index[c/(w_site-1)] = i*3;
        w_info[i] = c+w_num+s_num;
        for(j = 0; j < w_site-1; j++){
          atype[c+w_num+s_num+j] = 3+j;
          w_info[c+w_num+s_num+j] = i;
        }
        c += w_site-1;
      }

    for(i = 0; i < w_num*4; i += 4){
      ang0 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang1 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang2 = ((double)rand()/(double)RAND_MAX)*360*PI/180;
      ang[i  ] = sin(ang1/2)*sin((ang2-ang0)/2);
      ang[i+1] = sin(ang1/2)*cos((ang2-ang0)/2);
      ang[i+2] = cos(ang1/2)*sin((ang2+ang0)/2);
      ang[i+3] = cos(ang1/2)*cos((ang2+ang0)/2);
      angh[i  ] = ang[i  ];
      angh[i+1] = ang[i+1];
      angh[i+2] = ang[i+2];
      angh[i+3] = ang[i+3];
    }

  }

#if ZERO_P == 1

  for(i = 0; i < n3; i += 3){
    cd[i]   += 2*side[0];
    cd[i+1] += 2*side[1];
    cd[i+2] += 2*side[2];
  }
  for(i = 0; i < 3; i++)
    side[i] *= 5;
#endif

  alpha = oalpha / side[0];
  alpha2 = alpha*alpha;
  ial2si2 = 1. / (alpha*alpha*side[0]*side[0]);

  for(i = 0; i < 3; i++){
    sideh[i] = side[i] *.5;
    iside[i] = 1./side[i];
  }

  for(i = 0; i < KNUM+2; i++)
    atype_num[i] = 0;
  for(i = 0; i < n1; i++){
    atype_num[atype_mat[atype[i]]]++;
  }
  for(i = 0; i < n3; i++)
    vl[i] = 0;
  velset6(rtemp,h,tscale,knum,s_num*3+w_num*3);

  d6 = d7 = d8 = 0;
  for(i = 0; i < ws_num3; i += 3){
    d6 += cd[i];
    d7 += cd[i+1];
    d8 += cd[i+2];
  }
  d6 /= ws_num;
  d7 /= ws_num;
  d8 /= ws_num;
  for(i = 0; i < ws_num3; i += 3){
    cd[i]   -= d6;
    cd[i+1] -= d7;
    cd[i+2] -= d8;
  }

  d3 = d4 = d5 = 0;
  for(i = 0; i < ws_num3; i += 3){ /* calculae moment of inertia */
    d3 += (cd[i+1]*cd[i+1]+cd[i+2]*cd[i+2])*a_mass[atype_mat[atype[i/3]]];
    d4 += (cd[i]  *cd[i]  +cd[i+2]*cd[i+2])*a_mass[atype_mat[atype[i/3]]];
    d5 += (cd[i]  *cd[i]  +cd[i+1]*cd[i+1])*a_mass[atype_mat[atype[i/3]]];
  }
  d0 = d1 = d2 = 0;
  for(i = 0; i < ws_num3; i += 3){ /* calculate angular velocity */
    d0 +=(cd[i+1]*vl[i+2]-cd[i+2]*vl[i+1])*a_mass[atype_mat[atype[i/3]]]/d3;
    d1 +=(cd[i+2]*vl[i  ]-cd[i  ]*vl[i+2])*a_mass[atype_mat[atype[i/3]]]/d4;
    d2 +=(cd[i  ]*vl[i+1]-cd[i+1]*vl[i  ])*a_mass[atype_mat[atype[i/3]]]/d5;
  }
  for(i = 0; i < ws_num3; i += 3){
    vl[i]   -= d1*cd[i+2] - d2*cd[i+1];
    vl[i+1] -= d2*cd[i  ] - d0*cd[i+2];
    vl[i+2] -= d0*cd[i+1] - d1*cd[i  ];
  }
  for(i = 0; i < ws_num3; i += 3){
    cd[i]   += d6;
    cd[i+1] += d7;
    cd[i+2] += d8;
  }



  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i  ]*vl[i  ] +
              vl[i+1]*vl[i+1] +
              vl[i+2]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
  }
  ekin /= hsq;
  mtemp = tscale * ekin;

  d0 = sqrt(rtemp / mtemp);
  for(i = 0; i < n3; i++){
    vl[i] *= d0;
  }
  ekin = 0;
  for(i = 0; i < n3; i += 3){
    ekin += (vl[i  ]*vl[i  ] +
              vl[i+1]*vl[i+1] +
              vl[i+2]*vl[i+2])*a_mass[atype_mat[atype[i/3]]];
  }
  ekin /= hsq;
  mtemp = tscale * ekin;
  /*  printf("%f %f\n",mtemp*epsv/kb,ekin);*/


#if MDM == 2

  if(sys_num == 0){
#if ZERO_P  == 1
    strcpy(f_table_name,"table/fncl_af.table");
    strcpy(p_table_name,"table/fncl_ap.table");
    side_min = 3.5*side[0];
    side_max = 6.5*side[0];
    i0 = KNUM+4;
    for(i = 0; i < i0*i0; i++){
      gscale[i] = 0;
      rscale[i] = 0;
    }
    gscale[0*i0+0] = 1;
    gscale[0*i0+1] = 1;
    gscale[1*i0+0] = 1;
    gscale[1*i0+1] = 1;
    rscale[0*i0+0] = 1;
    rscale[0*i0+1] = pow(2,21);
    rscale[1*i0+0] = pow(2,21);
    rscale[1*i0+1] = pow(2,42);
#endif
  } else if(sys_num == 1 || sys_num == 2 || sys_num == 3 || sys_num == 10){
#if SPC == 1

#elif TIP5P == 1
    strcpy(f_table_name,"table/fx_ljclf.table");
    strcpy(p_table_name,"table/fx_ljclp.table");

    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));

    i10 = KNUM+4;

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        i0 = atype_mat[i];
        i1 = atype_mat[j];
        charge[i*i10+j] = z[i0]*z[i1];
        roffset[i*i10+j] = pow(2,32);
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale[i*i10+j] = d1*d1/d0;
          rscale[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale[i*i10+j] = 0;
          rscale[i*i10+j] = 0;
        }
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        d0 = wpa;
        d1 = wpc;
        if(d0 != 0){
          gscale2[i*i10+j] = 6.0*d1*pow(d1/d0,4.0/3.0);
          rscale2[i*i10+j] = pow(d1/d0,1.0/3.0);
        } else {
          gscale2[i*i10+j] = 0;
          rscale2[i*i10+j] = 0;
        }
      }
#endif

  } else if(sys_num >= 4){

    strcpy(f_table_name,"table/fx_ljclf.table");
    strcpy(p_table_name,"table/fx_ljclp.table");

    side_min = 0;
    side_max = my_max(side[0],my_max(side[1],side[2]));

    i10 = KNUM+4;

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
        i0 = (i == 4 ? 3:i);
        i1 = (j == 4 ? 3:j);
        charge[i*i10+j] = z[i0]*z[i1];
        roffset[i*i10+j] = pow(2,32);
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
	if(i <= 4 && j <= 4){
	  i0 = (i == 4 ? 3:i);
	  i1 = (j == 4 ? 3:j);
	  d0 = as_a[i0][i1];
	  d1 = as_c[i0][i1];
	  if(d0 != 0){
	    gscale[i*i10+j] = d1*d1/d0;
	    rscale[i*i10+j] = pow(d1/d0,1.0/3.0);
	  } else {
	    gscale[i*i10+j] = 0;
	    rscale[i*i10+j] = 0;
	  }
        } else {
          gscale[i*i10+j] = 0;
          rscale[i*i10+j] = 0;
        }
      }

    for(i = 0; i < i10; i++)
      for(j = 0; j < i10; j++){
	if(i <= 4 && j <= 4){
	  i0 = (i == 4 ? 3:i);
	  i1 = (j == 4 ? 3:j);
	  d0 = as_a[i0][i1];
	  d1 = as_c[i0][i1];
	  if(d0 != 0){
	    gscale2[i*i10+j] = 6.0*d1*pow(d1/d0,4.0/3.0);
	    rscale2[i*i10+j] = pow(d1/d0,1.0/3.0);
	  } else {
	    gscale2[i*i10+j] = 0;
	    rscale2[i*i10+j] = 0;
	  }
        } else {
          gscale2[i*i10+j] = 0;
          rscale2[i*i10+j] = 0;
        }
      }

  }
#endif

#if T_CONST == 1
  xs = 0;
  /*  lq = rtemp*(n3-3)*1e-2;*/
#endif

  phir_corr = 0;
  i0 = 2; i1 = 3; r = bond[0];
/*    phir_corr += (2.*zz[i0][i1]*erfc(alpha * r)/ r)*w_num;*/
  phir_corr -= (2.*zz[i0][i1]/ r)*w_num;

  i0 = 3; i1 = 3; r = bond[0]*sin(hoh_deg/2./180.*PI)*2;
/*    phir_corr += (zz[i0][i1]*erfc(alpha * r)/ r)*w_num;*/
  phir_corr -= (zz[i0][i1]/ r)*w_num;

  for(i0 = 0; i0 < w_num; i0++){
    i = w_index[i0];
    j = i0*4;
    c = w_info[i/3]*3;
    ang0 = ang[j  ];
    ang1 = ang[j+1];
    ang2 = ang[j+2];
    ang3 = ang[j+3];
    for(k = 0; k < w_site-1; k++){
      d0 = m_cdx[k]*(-ang0*ang0+ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdy[k]*(-2)*(ang0*ang1+ang2*ang3)
          +m_cdz[k]*( 2)*(ang1*ang2-ang0*ang3);
      d1 = m_cdx[k]*  2 *(ang2*ang3-ang0*ang1)
          +m_cdy[k]*( ang0*ang0-ang1*ang1-ang2*ang2+ang3*ang3)
          +m_cdz[k]*(-2)*(ang0*ang2+ang1*ang3);
      d2 = m_cdx[k]*  2 *(ang1*ang2+ang0*ang3)
          +m_cdy[k]*  2 *(ang1*ang3-ang0*ang2)
          +m_cdz[k]*(-ang0*ang0-ang1*ang1+ang2*ang2+ang3*ang3);
      cd[k*3+c  ] = cd[i  ] + d0;
      cd[k*3+c+1] = cd[i+1] + d1;
      cd[k*3+c+2] = cd[i+2] + d2;
    }
  }

  d0 = my_max(sideh[0],my_max(sideh[1],sideh[2]))/10;
  eye_len = d0/tan(15./180.*PI)*1.6*2;
  ini_temp = temp;

  for(i = 0; i < w_num3; i++){
    agv[i] = 0;
    agvh[i] = 0;
    trq[i] = 0;
  }
  for(i = 0; i <n3; i++){
    fc[i] = 0;
  }

}




void reshape(int w, int h){
  //glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glViewport(0, 0, w, h);
  glClearColor(clear_color, clear_color, clear_color, 0.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //mygluPerspective(30.0, (float)w / (float)h, 1.0, 800.0);
  mygluPerspective(45.0, (float)w / (float)h, 0.5f, 150.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

GLuint LoadTextureRAW( const char * filename, int wrap )

{
  GLuint texture;

  int width, height;
  unsigned char* data;
  FILE * file;
  // open texture data
  file = fopen( filename, "rb" );
  if ( file == NULL ) return 0;
  // allocate buffer
  width = 256;
  height = 256;
  data = malloc( width * height * 3 );
  // read texture data
  fread( data, width * height * 3, 1, file );
  fclose( file );
  // allocate a texture name
  glGenTextures( 1, &texture );
 // select our current texture
  glBindTexture( GL_TEXTURE_2D, texture );
 // select modulate to mix texture with color for shading
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
  // when texture area is small, bilinear filter the closest MIP map
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST );
  // when texture area is large, bilinear filter the first MIP map

  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  // if wrap is true, the texture wraps over at the edges (repeat)
  glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
  glTexImage2D(GL_TEXTURE_2D,0,3,width,height,0,GL_RGB,GL_UNSIGNED_BYTE,data);
  // free buffer
  free( data );
  return texture;
}
void small_font(double px, double py, double pz, char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, moji[i]);
  }
}
void medium_font(double px, double py, double pz, char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, moji[i]);
  }
}

void large_font(double px, double py, double pz, char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, moji[i]);
  }
}


void display(){
	int which=0;
	float d0,d1,d2,d3,d4;
	int i,j;
	int w=800,h=800;
	char str_buf[256];
	char str_buf2[256];

	GLfloat particle_color[4];
	GLfloat color[4];
	///////For MOuse Movement;


	//////////////Reshape//////////////////////
#if 1
	glViewport(0, 0, w, h);
	glClearColor(clear_color, clear_color, clear_color, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//mygluPerspective(30.0, (float)w / (float)h, 1.0, 800.0);
	mygluPerspective(45.0, (float)w / (float)h, 0.5f, 800.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
#endif
	///////////////////////////////////////////

	//glClearColor(clear_color, clear_color, clear_color, 0.0);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#if 1 /////for color
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glCullFace(GL_BACK);

#endif


	glPushMatrix();
	d3 = atan((eye_width*which)/eye_len);
	d1 = sin(d3)*eye_len;
	d0 = cos(d3)*eye_len;
	//mygluLookAt(d0, d1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	mygluLookAt(d0*3, d1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

#if 1
	//  glTranslated(trans[0], trans[1], trans[2]);
	glTranslatef(trans[0], trans[1], trans[2]);


	glPushMatrix();
	glLoadIdentity();
	glRotatef( angle[0],1.0,0.0,0.0);
	glRotatef( angle[1],0.0,1.0,0.0);
	glRotatef( angle[2],0.0,0.0,1.0);
	//glMultMatrixd(m_matrix);
	glMultMatrixf(m_matrix);

	//glGetDoublev(GL_MODELVIEW_MATRIX, m_matrix);
	glGetFloatv(GL_MODELVIEW_MATRIX, m_matrix);
	glPopMatrix();


	for(i = 0; i < 16; i++)
	i_matrix[i] = m_matrix[i];
	mat_inv((float(*)[4])i_matrix);

	//glMultMatrixd(m_matrix);
	glMultMatrixf(m_matrix);


	//////Mouse Detection

	angle[0] = 0;
	if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
	angle[1] = 0;
	angle[2] = 0;
	}
	if(ini_flg == 1){
	mouse_l = 0;
	ini_flg = 0;
	}

#endif



	d2 = (i_matrix[0]*(1.0)+i_matrix[4]*(1.0)+i_matrix[8]*(1.0));
	d3 = (i_matrix[1]*(1.0)+i_matrix[5]*(1.0)+i_matrix[9]*(1.0));
	d4 = (i_matrix[2]*(1.0)+i_matrix[6]*(1.0)+i_matrix[10]*(1.0));

	d0 = side0/2;
	d1 = 0.3;
	color[0] = d1; color[1] = d1; color[2] = d1; color[3] = 1.0;
	glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,color);
	glLineWidth(2.0);
	glPointSize(9.0);


#if 1
	///////For Drawing with draw arrays
	GLfloat line1[]={0,0,-d0,0,0,d0};
	GLfloat Nline[]={d2,d3,d4};
	GLfloat line2[]={0,-d0,0,0,d0,0};
	GLfloat line3[]={-d0,0,0,d0,0,0};
	GLfloat white[]={1.0,1.0,1.0,1.0};
	float side_s[3],side_e[3];

	for(i = 0; i < 3; i++){
	   side_s[i] = -sideh[i];
	   side_e[i] = side[i]-sideh[i];
	}

	//////Coordinates for Cube
	float cubevert1[]={  side_s[0],side_s[1],side_s[2],
					   side_s[0],side_e[1],side_s[2],
					   side_s[0],side_e[1],side_e[2],
					   side_s[0],side_s[1],side_e[2]
	};

	float cubevert2[]={ side_e[0],side_s[1],side_s[2],
					   side_e[0],side_e[1],side_s[2],
					   side_e[0],side_e[1],side_e[2],
					   side_e[0],side_s[1],side_e[2]
	};

	float cubebody[]={side_e[0],side_s[1],side_s[2],
		   side_s[0],side_s[1],side_s[2],
		   side_e[0],side_e[1],side_s[2],
		   side_s[0],side_e[1],side_s[2],
		   side_e[0],side_e[1],side_e[2],
		   side_s[0],side_e[1],side_e[2],
		   side_e[0],side_s[1],side_e[2],
		   side_s[0],side_s[1],side_e[2]
	};


	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glEnable(GL_COLOR_MATERIAL);

	/////Drawing the Cross 3D Axis X,Y,Z
	glVertexPointer(3,GL_FLOAT,0,line1);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINES,0,2);

	glVertexPointer(3,GL_FLOAT,0,line2);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINES,0,2);

	glVertexPointer(3,GL_FLOAT,0,line3);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINES,0,2);

	///////Drawing the Cube lines
	glVertexPointer(3,GL_FLOAT,0,cubevert1);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINE_LOOP,0,4);

	glVertexPointer(3,GL_FLOAT,0,cubevert2);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINE_LOOP,0,4);

	glVertexPointer(3,GL_FLOAT,0,cubebody);
	glNormalPointer(GL_FLOAT,0,Nline);
	glColorPointer(4,GL_FLOAT,0,white);
	glDrawArrays(GL_LINES,0,8);

	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);
	//glDisableClientState(GL_COLOR_ARRAY);

#endif


	///////////////////////////////////Drawing Spheres with DrawElements
#if 1
	GLuint buf[3];
	glGenBuffers(3, buf);

	int n = n3/3;
	int q,m,l=3,k,p,t=0,r=0,slices =ditail,stacks=ditail/2;
	float z0,z1,r0,r1,radios;

	float *f_pol,*f_pol_n,*f_clr_a,*f_clr_b;
	float *sint1,*cost1;
	float *sint2,*cost2;

	int cuad_mem    =(((ditail/2)-2)*ditail*6*3);//+1; por si ditail=5  =>0
	int pol_mem 	= cuad_mem + (ditail*3*3*2);
	int pol_size 	=(ditail*3*2)+(cuad_mem/3);

	unsigned int size 		= (n*pol_mem)*sizeof(float);
	unsigned int size_color = (n*pol_size*4)*sizeof(float);


	f_pol    = (float*)malloc(n*pol_mem*sizeof(float));
	f_pol_n  = (float*)malloc(n*pol_mem*sizeof(float));
	f_clr_a	 = (float*)malloc(n*pol_size*4*sizeof(float));
	f_clr_b	 = (float*)malloc(n*pol_size*4*sizeof(float));

	/* Pre-computed circle */


	CircleTable(&sint1,&cost1,-slices);
	CircleTable(&sint2,&cost2,stacks*2);

//////Mapping The circle////////////////////////////////////////////////////
#if 1

	for(i=0; i<n3;i+=3){
	   if(drow_flg[atype_mat[atype[i/3]]] == 1){
		   /////////////////////Compute Color
		   for(q=0;q<pol_size;q++){
			d0 = (vl[i]*vl[i]+vl[i+1]*vl[i+1]+vl[i+2]*vl[i+2])*500;
			*(f_clr_a+0+r+q*4) = color_table[atype_mat[atype[i/3]]][0]+d0;
			*(f_clr_a+1+r+q*4) = color_table[atype_mat[atype[i/3]]][1]+d0/3;
			*(f_clr_a+2+r+q*4) = color_table[atype_mat[atype[i/3]]][2]+d0/3;
			*(f_clr_a+3+r+q*4) = color_table[atype_mat[atype[i/3]]][3];
/*
			f_clr_b[0+r+q] = color_table[atype_mat[atype[i/3]]][0]+d0/4;
			f_clr_b[1+r+q] = color_table[atype_mat[atype[i/3]]][1]+d0/12;
			f_clr_b[2+r+q] = color_table[atype_mat[atype[i/3]]][2]+d0/12;
			f_clr_b[3+r+q] = color_table[atype_mat[atype[i/3]]][3];
			*/
		   }

			radios=radius*r_table[atype_mat[atype[i/3]]];
			///////Compute one Circle/////////////////////////
					z0 = 1.0f;
			        z1 = cost2[(stacks>0)?1:0];
			        r0 = 0.0f;
			        r1 = sint2[(stacks>0)?1:0];

			        for(m=0;m<slices;m++){
			                *(f_pol+t+(9*m))     =cd[i]-sideh[0];
			                *(f_pol+t+1+(9*m))   =cd[i+1]-sideh[1];
			                *(f_pol+t+2+(9*m))   =cd[i+2]-sideh[2]+radios;
			                *(f_pol_n+t+(9*m))    =0;
			                *(f_pol_n+t+1+(9*m))  =0;
			                *(f_pol_n+t+2+(9*m))  =1;
			            }

			            l=3;

			            for (j=slices; j>0; j--){
			                *(f_pol+t+l)     =(cd[i]-sideh[0])+(cost1[j]*r1*radios);
			                *(f_pol+t+l+1)   =(cd[i+1]-sideh[1])+(sint1[j]*r1*radios);
			                *(f_pol+t+l+2)   =(cd[i+2]-sideh[2])+ (z1*radios);
			                *(f_pol_n+t+l)   =cost1[j]*r1;
			                *(f_pol_n+t+l+1) =sint1[j]*r1;
			                *(f_pol_n+t+l+2) =z1;

			                *(f_pol+t+l+3)   =(cd[i]-sideh[0])  + (cost1[j-1]*r1*radios);
			                *(f_pol+t+l+4)   =(cd[i+1]-sideh[1])+ (sint1[j-1]*r1*radios);
			                *(f_pol+t+l+5)   =(cd[i+2]-sideh[2])+ (z1*radios);
			                *(f_pol_n+t+l+3) =cost1[j-1]*r1;
			                *(f_pol_n+t+l+4) =sint1[j-1]*r1;
			                *(f_pol_n+t+l+5) =z1;

			                l+=9;
			            }
			/////////////////
			           l-=3;
			            for( k=1; k<stacks-1; k++ ){
			                z0 = z1; z1 = cost2[k+1];
			                r0 = r1; r1 = sint2[k+1];

			                p=0;
			                for(j=0; j<slices; j++){
			                	 //////////////////First Triangle////////////////////////////////
			                    *(f_pol+t+l+p)     = (cd[i]-sideh[0] ) +(cost1[j]*r1*radios);
			                    *(f_pol+t+l+p+1)   = (cd[i+1]-sideh[1])+(sint1[j]*r1*radios);
			                    *(f_pol+t+l+p+2)   = (cd[i+2]-sideh[2])+(z1*radios);
			                    *(f_pol_n+t+l+p)   = cost1[j]*r1;
			                    *(f_pol_n+t+l+p+1) = sint1[j]*r1;
			                    *(f_pol_n+t+l+p+2) = z1;

			                    *(f_pol+t+l+p+3)   = (cd[i]-sideh[0])  +(cost1[j]*r0*radios);
			                    *(f_pol+t+l+p+4)   = (cd[i+1]-sideh[1])+(sint1[j]*r0*radios);
			                    *(f_pol+t+l+p+5)   = (cd[i+2]-sideh[2])+(z0*radios);
			                    *(f_pol_n+t+l+p+3) = cost1[j]*r0;
			                    *(f_pol_n+t+l+p+4) = sint1[j]*r0;
			                    *(f_pol_n+t+l+p+5) = z0;

			                    *(f_pol+t+l+p+6)   = (cd[i]-sideh[0])  +(cost1[j+1]*r1*radios);
			                    *(f_pol+t+l+p+7)   = (cd[i+1]-sideh[1])+(sint1[j+1]*r1*radios);
			                    *(f_pol+t+l+p+8)   = (cd[i+2]-sideh[2])+(z1*radios);
			                    *(f_pol_n+t+l+p+6) = cost1[j+1]*r1;
			                    *(f_pol_n+t+l+p+7) = sint1[j+1]*r1;
			                    *(f_pol_n+t+l+p+8) = z1;

			                    /////////////////////Second Triangle//////////////////////////////

			                    *(f_pol+t+l+p+9)   = *(f_pol+t+l+p+6);
			                    *(f_pol+t+l+p+10)  = *(f_pol+t+l+p+7);
			                    *(f_pol+t+l+p+11)  = *(f_pol+t+l+p+8);
			                    *(f_pol_n+t+l+p+9) = *(f_pol_n+t+l+p+6);
			                    *(f_pol_n+t+l+p+10)= *(f_pol_n+t+l+p+7);
			                    *(f_pol_n+t+l+p+11)= *(f_pol_n+t+l+p+8);

			                    *(f_pol+t+l+p+12)   = *(f_pol+t+l+p+3);
			                    *(f_pol+t+l+p+13)   = *(f_pol+t+l+p+4);
			                    *(f_pol+t+l+p+14)   = *(f_pol+t+l+p+5);
			                    *(f_pol_n+t+l+p+12) = *(f_pol_n+t+l+p+3);
			                    *(f_pol_n+t+l+p+13) = *(f_pol_n+t+l+p+4);
			                    *(f_pol_n+t+l+p+14) = *(f_pol_n+t+l+p+5);

			                    *(f_pol+t+l+p+15)  =(cd[i]-sideh[0] ) +(cost1[j+1]*r0*radios);
			                    *(f_pol+t+l+p+16)  =(cd[i+1]-sideh[1])+(sint1[j+1]*r0*radios);
			                    *(f_pol+t+l+p+17)  =(cd[i+2]-sideh[2])+(z0*radios);
			                    *(f_pol_n+t+l+p+15)=cost1[j+1]*r0;
			                    *(f_pol_n+t+l+p+16)=sint1[j+1]*r0;
			                    *(f_pol_n+t+l+p+17)=z0;

			                    p+=18;
			                }
			            l+=(slices)*6*3;
			            }
			////////////////////
			            z0 = z1;
			            r0 = r1;

			            for(m=0;m<slices;m++){
							*(f_pol+t+l+(9*m))     = cd[i]-sideh[0];
							*(f_pol+t+l+1+(9*m))   = cd[i+1]-sideh[1];
							*(f_pol+t+l+2+(9*m))   = cd[i+2]-sideh[2]-radios;
							*(f_pol_n+t+l+(9*m))   = 0;
							*(f_pol_n+t+l+1+(9*m)) = 0;
							*(f_pol_n+t+l+2+(9*m)) = -1;
			            }

			            p=3;
			            for (j=0; j<slices; j++){
			                *(f_pol+t+l+p)     = (cd[i]-sideh[0])  + (cost1[j]*r0*radios);
			                *(f_pol+t+l+p+1)   = (cd[i+1]-sideh[1])+ (sint1[j]*r0*radios);
			                *(f_pol+t+l+p+2)   = (cd[i+2]-sideh[2])+ (z0*radios);
			                *(f_pol_n+t+l+p)   = cost1[j]*r0;
			                *(f_pol_n+t+l+p+1) = sint1[j]*r0;
			                *(f_pol_n+t+l+p+2) = z0;

			                *(f_pol+t+l+p+3)   = (cd[i]-sideh[0] ) +(cost1[j+1]*r0*radios);
			                *(f_pol+t+l+p+4)   = (cd[i+1]-sideh[1])+(sint1[j+1]*r0*radios);
			                *(f_pol+t+l+p+5)   = (cd[i+2]-sideh[2])+(z0*radios);
			                *(f_pol_n+t+l+p+3) = cost1[j+1]*r0;
			                *(f_pol_n+t+l+p+4) = sint1[j+1]*r0;
			                *(f_pol_n+t+l+p+5) = z0;

			                p+=9;
			            }

	   	   }
	   	   t+=pol_mem;
	   	   r+=pol_size*4;
	}

//////////////////////////Prepare Circle mapped buffers
	//glEnableClientState(GL_VERTEX_ARRAY);
	//glEnableClientState(GL_NORMAL_ARRAY);
	//glEnableClientState(GL_COLOR_ARRAY);
	//glEnable(GL_COLOR_MATERIAL);

	glBindBuffer(GL_ARRAY_BUFFER, buf[0]);
	glBufferData(GL_ARRAY_BUFFER,size,f_pol, GL_DYNAMIC_DRAW);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, buf[1]);
	glBufferData(GL_ARRAY_BUFFER, size,f_pol_n, GL_DYNAMIC_DRAW);
	glNormalPointer(GL_FLOAT,0,0);

	glBindBuffer(GL_ARRAY_BUFFER, buf[2]);
	glBufferData(GL_ARRAY_BUFFER,size_color,f_clr_a, GL_DYNAMIC_DRAW);
	glColorPointer(4,GL_FLOAT,0,0);




	glDrawArrays(GL_TRIANGLES,0,pol_size*n);

#endif
	///////////////////////END of Drawing circle//////////////////////////



#endif


	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisable(GL_COLOR_MATERIAL);

	glPopMatrix();


	gettimeofday(&time_v,NULL);
	disp_time0 = disp_time;
	disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
/////////	////////////Drawing the other part
/////////	//////////////////////////////////
/////////	//////////////////////////////////
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glColor4fv(moji_c[(int)(clear_color+.5)]);
	d0 = -2.4;
	d1 = 2.2;
	d2 = -10;
	d1 = -2.2;
	sprintf(str_buf,"%.3fs/step %.1fGflops",md_time-md_time0
			 ,(double)n1*(double)n1*78/(md_time-md_time0)*1e-9);
	//glRasterPos3d(d0, d1, d2);
	//glCallLists(strlen(str_buf), GL_BYTE, str_buf);
	large_font(d0, d1, d2,str_buf);
	d1 -= .3;
	sprintf(str_buf,"%.3fs/frm %.1ffrm/s",disp_time-disp_time0
	 ,1./(disp_time-disp_time0));
	//glRasterPos3d(d0, d1, d2);
	//glCallLists(strlen(str_buf), GL_BYTE, str_buf);
	large_font(d0, d1, d2,str_buf);

///////////////
	glEnable(GL_LIGHTING);
    d0 = 1.8;
    d1 = -2.4;
    d2 = -10;
    glPushMatrix();
    glTranslated(d0,d1,d2);
    d3 = temp*0.00016;
    particle_color[0] = color_table[atype_mat[0]][0]+d3;
    particle_color[1] = color_table[atype_mat[0]][1]+d3/3;
    particle_color[2] = color_table[atype_mat[0]][2]+d3/3;
    particle_color[3] = color_table[atype_mat[0]][3];
    glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
    particle_color[0] = color_table[atype_mat[0]][0]+d3/4;
    particle_color[1] = color_table[atype_mat[0]][1]+d3/12;
    particle_color[2] = color_table[atype_mat[0]][2]+d3/12;
    particle_color[3] = color_table[atype_mat[0]][3];
    glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
    glutSolidSphere(radius*r_table[atype_mat[0]]/3.5, ditail, ditail/2);
    glPopMatrix();
    glDisable(GL_LIGHTING);
    medium_font(d0-0.08,d1-.04,d2,"Na");
    small_font(d0+0.08,d1+.04,d2,"+");

    glEnable(GL_LIGHTING);
    d0 += 0.5;
    glPushMatrix();
    glTranslated(d0,d1,d2);
    particle_color[0] = color_table[atype_mat[1]][0]+d3;
    particle_color[1] = color_table[atype_mat[1]][1]+d3/3;
    particle_color[2] = color_table[atype_mat[1]][2]+d3/3;
    particle_color[3] = color_table[atype_mat[1]][3];
    glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
    particle_color[0] = color_table[atype_mat[1]][0]+d3/4;
    particle_color[1] = color_table[atype_mat[1]][1]+d3/12;
    particle_color[2] = color_table[atype_mat[1]][2]+d3/12;
    particle_color[3] = color_table[atype_mat[1]][3];
    glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
    glutSolidSphere(radius*r_table[atype_mat[1]]/3.5, ditail, ditail/2);
    glPopMatrix();
    glDisable(GL_LIGHTING);
    medium_font(d0-0.06,d1-.04,d2,"Cl");
    small_font(d0+0.07,d1+.04,d2,"-");



    glDeleteBuffers(3, buf);

	free(f_clr_a);
	free(f_clr_b);
	free(f_pol_n);
	free(f_pol);
	free(sint1);
	free(cost1);
	free(sint2);
	free(cost2);

	glutSwapBuffers();
}
void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_l = 1;
    }
    if (state == GLUT_UP) {
      mouse_l = 0;
    }
    break;
  case GLUT_RIGHT_BUTTON:
     if (state == GLUT_DOWN) {
       mpos[0] = x;
       mpos[1] = y;
       mouse_r = 1;
     }
     if (state == GLUT_UP) {
       mouse_r = 0;
     }
     break;

  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_m = 1;
    }
    if (state == GLUT_UP) {
      mouse_m = 0;
    }
    break;


  default:
    break;
  }
}
void motion(int x, int y)
{
	float d0;
	float len = 10;
	len = eye_len;

	//printf("Mouse X=%i, Y=%i\n",mpos[0],mpos[1]);
	//printf("MOTION X=%i, Y=%i\n",x,y);
	if(mouse_l == 1 && mouse_m == 1){
		trans[0] += (float)(y-mpos[1])*len/150;
		angle[0] = -(float)(x-mpos[0])*0.2;
	} else  if(mouse_m == 1 || (mouse_l == 1 && mouse_r == 1)){
		trans[1] += (float)(x-mpos[0])*len*.001;
		trans[2] -= (float)(y-mpos[1])*len*.001;
	} else if(mouse_r == 1){
		trans[0] -= (float)(y-mpos[1])*len/150;
		angle[0] =  (float)(x-mpos[0])*0.2;
	} else if(mouse_l == 1){
		d0 = len/50;
		if(d0 > 1.0) d0 = 1.0;
			angle[1] = (float)(y-mpos[1])*d0;
			angle[2] = (float)(x-mpos[0])*d0;
		}
	if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
		mpos[0] = x;
		mpos[1] = y;
		glutPostRedisplay();
	}
}

int main(int argc, char **argv){
  //Number of Particles according to the keyboard numbers
	np = 3;
	temp = 300;
	sub_x = 1.5;
	sub_y = 1.5;
	temp_ymax = 2000;


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(100, 0);
	glutCreateWindow(argv[0]);

	init();
	init_MD();
	keep_mem(S_NUM_MAX,W_NUM_MAX*w_site);
	set_cd(1);

	gettimeofday(&time_v,NULL);
	disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(md_run);
	glutMainLoop();



	return 0;
}

