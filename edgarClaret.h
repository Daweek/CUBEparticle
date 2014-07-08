/*
 * edgarClaret.h
 *
 *  Created on: Jul 8, 2014
 *      Author: edgar
 */

///////////////////////////Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include <sys/time.h>
#include <pthread.h>
#include <sched.h>
#define NUMT     1

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
	int sstart,sstop;

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

	sstart = id*range*3;
	sstop  = (id+1)*(range*3);


		for (i = sstart; i < sstop ; i += 3) {
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

