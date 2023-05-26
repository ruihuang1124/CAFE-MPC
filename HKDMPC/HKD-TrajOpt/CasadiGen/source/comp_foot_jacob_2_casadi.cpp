/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) comp_foot_jacob_2_casadi_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const casadi_int casadi_s0[7] = {3, 1, 0, 3, 0, 1, 2};
static const casadi_int casadi_s1[75] = {3, 18, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};

/* comp_foot_jacob_2:(i0[3],i1[3],i2[3])->(o0[3x18]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a6, a7, a8, a9;
  a0=arg[2]? arg[2][2] : 0;
  a1=sin(a0);
  a2=arg[2]? arg[2][1] : 0;
  a3=sin(a2);
  a4=(a1*a3);
  a5=cos(a0);
  a6=cos(a2);
  a7=(a5*a6);
  a4=(a4-a7);
  a7=arg[1]? arg[1][1] : 0;
  a8=cos(a7);
  a9=(a4*a8);
  a10=1.2246467991473532e-16;
  a11=(a10*a6);
  a12=(a5*a11);
  a13=(a10*a3);
  a14=(a1*a13);
  a12=(a12-a14);
  a14=arg[2]? arg[2][0] : 0;
  a15=cos(a14);
  a16=(a12*a15);
  a17=(a5*a3);
  a18=(a1*a6);
  a17=(a17+a18);
  a18=sin(a14);
  a19=(a17*a18);
  a16=(a16+a19);
  a19=arg[1]? arg[1][2] : 0;
  a20=sin(a19);
  a21=(a16*a20);
  a12=(a12*a18);
  a17=(a17*a15);
  a12=(a12-a17);
  a17=cos(a19);
  a22=(a12*a17);
  a21=(a21+a22);
  a22=sin(a7);
  a23=(a21*a22);
  a9=(a9+a23);
  a23=arg[1]? arg[1][0] : 0;
  a24=cos(a23);
  a25=(a9*a24);
  a16=(a16*a17);
  a12=(a12*a20);
  a16=(a16-a12);
  a12=sin(a23);
  a26=(a16*a12);
  a25=(a25-a26);
  a26=cos(a0);
  a27=cos(a2);
  a28=(a26*a27);
  a0=sin(a0);
  a2=sin(a2);
  a29=(a0*a2);
  a28=(a28-a29);
  a29=(a10*a28);
  a30=cos(a14);
  a31=(a29*a30);
  a32=(a26*a2);
  a33=(a0*a27);
  a32=(a32+a33);
  a14=sin(a14);
  a33=(a32*a14);
  a31=(a31+a33);
  a33=sin(a19);
  a34=(a31*a33);
  a29=(a29*a14);
  a35=(a32*a30);
  a29=(a29-a35);
  a19=cos(a19);
  a35=(a29*a19);
  a34=(a34+a35);
  a35=sin(a7);
  a36=(a34*a35);
  a7=cos(a7);
  a37=(a28*a7);
  a36=(a36-a37);
  a37=cos(a23);
  a38=(a36*a37);
  a39=(a31*a19);
  a40=(a29*a33);
  a39=(a39-a40);
  a23=sin(a23);
  a40=(a39*a23);
  a38=(a38-a40);
  a40=(a25*a38);
  a41=-1.2246467991473532e-16;
  a42=(a41*a8);
  a43=(a15*a20);
  a44=(a18*a17);
  a43=(a43+a44);
  a44=(a43*a22);
  a42=(a42-a44);
  a44=(a42*a24);
  a45=(a18*a20);
  a46=(a15*a17);
  a45=(a45-a46);
  a46=(a45*a12);
  a44=(a44-a46);
  a46=(a41*a7);
  a47=(a30*a33);
  a48=(a14*a19);
  a47=(a47+a48);
  a48=(a47*a35);
  a46=(a46-a48);
  a48=(a46*a37);
  a49=(a14*a33);
  a50=(a30*a19);
  a49=(a49-a50);
  a50=(a49*a23);
  a48=(a48-a50);
  a50=(a44*a48);
  a40=(a40+a50);
  a11=(a1*a11);
  a13=(a5*a13);
  a11=(a11+a13);
  a13=(a11*a15);
  a50=(a5*a6);
  a51=(a1*a3);
  a50=(a50-a51);
  a51=(a50*a18);
  a13=(a13-a51);
  a51=(a13*a20);
  a11=(a11*a18);
  a50=(a50*a15);
  a11=(a11+a50);
  a50=(a11*a17);
  a51=(a51+a50);
  a50=(a51*a22);
  a1=(a1*a6);
  a5=(a5*a3);
  a1=(a1+a5);
  a5=(a1*a8);
  a50=(a50-a5);
  a5=(a50*a24);
  a13=(a13*a17);
  a11=(a11*a20);
  a13=(a13-a11);
  a11=(a13*a12);
  a5=(a5-a11);
  a11=(a0*a27);
  a20=(a26*a2);
  a11=(a11+a20);
  a20=(a10*a11);
  a17=(a20*a30);
  a3=(a26*a27);
  a6=(a0*a2);
  a3=(a3-a6);
  a6=(a3*a14);
  a17=(a17-a6);
  a6=(a17*a33);
  a20=(a20*a14);
  a15=(a3*a30);
  a20=(a20+a15);
  a15=(a20*a19);
  a6=(a6+a15);
  a15=(a6*a35);
  a18=(a11*a7);
  a15=(a15-a18);
  a18=(a15*a37);
  a52=(a17*a19);
  a53=(a20*a33);
  a52=(a52-a53);
  a53=(a52*a23);
  a18=(a18-a53);
  a53=(a5*a18);
  a40=(a40+a53);
  if (res[0]!=0) res[0][0]=a40;
  a9=(a9*a12);
  a16=(a16*a24);
  a9=(a9+a16);
  a16=(a9*a38);
  a42=(a42*a12);
  a45=(a45*a24);
  a42=(a42+a45);
  a45=(a42*a48);
  a16=(a16+a45);
  a50=(a50*a12);
  a13=(a13*a24);
  a50=(a50+a13);
  a13=(a50*a18);
  a16=(a16+a13);
  if (res[0]!=0) res[0][1]=a16;
  a21=(a21*a8);
  a4=(a4*a22);
  a21=(a21-a4);
  a38=(a21*a38);
  a4=(a41*a22);
  a43=(a43*a8);
  a4=(a4+a43);
  a48=(a4*a48);
  a38=(a38-a48);
  a1=(a1*a22);
  a51=(a51*a8);
  a1=(a1+a51);
  a18=(a1*a18);
  a38=(a38+a18);
  if (res[0]!=0) res[0][2]=a38;
  a36=(a36*a23);
  a39=(a39*a37);
  a36=(a36+a39);
  a39=(a25*a36);
  a46=(a46*a23);
  a49=(a49*a37);
  a46=(a46+a49);
  a49=(a44*a46);
  a39=(a39+a49);
  a15=(a15*a23);
  a52=(a52*a37);
  a15=(a15+a52);
  a52=(a5*a15);
  a39=(a39+a52);
  if (res[0]!=0) res[0][3]=a39;
  a39=(a9*a36);
  a52=(a42*a46);
  a39=(a39+a52);
  a52=(a50*a15);
  a39=(a39+a52);
  if (res[0]!=0) res[0][4]=a39;
  a36=(a21*a36);
  a46=(a4*a46);
  a36=(a36-a46);
  a15=(a1*a15);
  a36=(a36+a15);
  if (res[0]!=0) res[0][5]=a36;
  a36=(a28*a35);
  a34=(a34*a7);
  a36=(a36+a34);
  a34=(a25*a36);
  a15=(a41*a35);
  a47=(a47*a7);
  a15=(a15+a47);
  a47=(a44*a15);
  a34=(a34-a47);
  a47=(a11*a35);
  a6=(a6*a7);
  a47=(a47+a6);
  a6=(a5*a47);
  a34=(a34+a6);
  if (res[0]!=0) res[0][6]=a34;
  a34=(a9*a36);
  a6=(a42*a15);
  a34=(a34-a6);
  a6=(a50*a47);
  a34=(a34+a6);
  if (res[0]!=0) res[0][7]=a34;
  a36=(a21*a36);
  a15=(a4*a15);
  a36=(a36+a15);
  a47=(a1*a47);
  a36=(a36+a47);
  if (res[0]!=0) res[0][8]=a36;
  a36=-1.9000000000000000e-01;
  a47=(a36*a29);
  a15=-1.9500000000000001e-01;
  a34=-2.0899999999999999e-01;
  a6=(a34*a26);
  a6=(a15+a6);
  a46=(a6*a30);
  a39=6.2000000000000000e-02;
  a52=(a39*a28);
  a37=(a52*a14);
  a46=(a46+a37);
  a47=(a47-a46);
  a46=(a47*a33);
  a52=(a52*a30);
  a37=(a6*a14);
  a52=(a52-a37);
  a37=-4.9000000000000002e-02;
  a28=(a37*a28);
  a52=(a52-a28);
  a28=1.9000000000000000e-01;
  a31=(a28*a31);
  a52=(a52+a31);
  a31=(a52*a19);
  a46=(a46+a31);
  a46=(a46*a7);
  a31=(a41*a6);
  a32=(a39*a32);
  a31=(a31-a32);
  a32=4.9000000000000002e-02;
  a29=(a32*a29);
  a29=(a31+a29);
  a23=(a29*a35);
  a46=(a46-a23);
  a23=(a25*a46);
  a49=1.9500000000000001e-01;
  a26=(a49*a26);
  a38=2.0899999999999999e-01;
  a26=(a26+a38);
  a38=(a26*a27);
  a49=(a49*a0);
  a18=(a49*a2);
  a38=(a38-a18);
  a18=(a32*a14);
  a18=(a38+a18);
  a51=(a18*a35);
  a10=(a10*a38);
  a8=(a10*a30);
  a22=7.5928101547135900e-18;
  a26=(a26*a2);
  a49=(a49*a27);
  a26=(a26+a49);
  a22=(a22-a26);
  a26=(a22*a14);
  a8=(a8-a26);
  a26=(a36*a14);
  a8=(a8-a26);
  a26=(a8*a33);
  a10=(a10*a14);
  a22=(a22*a30);
  a10=(a10+a22);
  a22=6.0007693158220313e-18;
  a10=(a10+a22);
  a22=(a28*a30);
  a10=(a10-a22);
  a22=(a10*a19);
  a26=(a26+a22);
  a26=(a26*a7);
  a51=(a51+a26);
  a26=(a44*a51);
  a23=(a23+a26);
  a36=(a36*a20);
  a34=(a34*a0);
  a0=(a34*a30);
  a26=(a39*a11);
  a22=(a26*a14);
  a0=(a0+a22);
  a36=(a36-a0);
  a0=(a36*a33);
  a26=(a26*a30);
  a14=(a34*a14);
  a26=(a26-a14);
  a37=(a37*a11);
  a26=(a26-a37);
  a28=(a28*a17);
  a26=(a26+a28);
  a28=(a26*a19);
  a0=(a0+a28);
  a0=(a0*a7);
  a41=(a41*a34);
  a39=(a39*a3);
  a41=(a41+a39);
  a32=(a32*a20);
  a32=(a41+a32);
  a35=(a32*a35);
  a0=(a0-a35);
  a35=(a5*a0);
  a23=(a23+a35);
  if (res[0]!=0) res[0][9]=a23;
  a23=(a9*a46);
  a35=(a42*a51);
  a23=(a23+a35);
  a35=(a50*a0);
  a23=(a23+a35);
  if (res[0]!=0) res[0][10]=a23;
  a46=(a21*a46);
  a51=(a4*a51);
  a46=(a46-a51);
  a0=(a1*a0);
  a46=(a46+a0);
  if (res[0]!=0) res[0][11]=a46;
  a47=(a47*a19);
  a52=(a52*a33);
  a47=(a47-a52);
  a52=(a25*a47);
  a8=(a8*a19);
  a10=(a10*a33);
  a8=(a8-a10);
  a10=(a44*a8);
  a52=(a52+a10);
  a36=(a36*a19);
  a26=(a26*a33);
  a36=(a36-a26);
  a26=(a5*a36);
  a52=(a52+a26);
  if (res[0]!=0) res[0][12]=a52;
  a52=(a9*a47);
  a26=(a42*a8);
  a52=(a52+a26);
  a26=(a50*a36);
  a52=(a52+a26);
  if (res[0]!=0) res[0][13]=a52;
  a47=(a21*a47);
  a8=(a4*a8);
  a47=(a47-a8);
  a36=(a1*a36);
  a47=(a47+a36);
  if (res[0]!=0) res[0][14]=a47;
  a47=(a25*a29);
  a36=(a44*a18);
  a47=(a47-a36);
  a36=(a5*a32);
  a47=(a47+a36);
  if (res[0]!=0) res[0][15]=a47;
  a47=(a9*a29);
  a36=(a42*a18);
  a47=(a47-a36);
  a36=(a50*a32);
  a47=(a47+a36);
  if (res[0]!=0) res[0][16]=a47;
  a29=(a21*a29);
  a18=(a4*a18);
  a29=(a29+a18);
  a32=(a1*a32);
  a29=(a29+a32);
  if (res[0]!=0) res[0][17]=a29;
  a29=0.;
  if (res[0]!=0) res[0][18]=a29;
  if (res[0]!=0) res[0][19]=a29;
  if (res[0]!=0) res[0][20]=a29;
  if (res[0]!=0) res[0][21]=a29;
  if (res[0]!=0) res[0][22]=a29;
  if (res[0]!=0) res[0][23]=a29;
  if (res[0]!=0) res[0][24]=a29;
  if (res[0]!=0) res[0][25]=a29;
  if (res[0]!=0) res[0][26]=a29;
  a32=(a25*a31);
  a44=(a44*a38);
  a32=(a32-a44);
  a44=(a5*a41);
  a32=(a32+a44);
  if (res[0]!=0) res[0][27]=a32;
  a32=(a9*a31);
  a42=(a42*a38);
  a32=(a32-a42);
  a42=(a50*a41);
  a32=(a32+a42);
  if (res[0]!=0) res[0][28]=a32;
  a31=(a21*a31);
  a4=(a4*a38);
  a31=(a31+a4);
  a41=(a1*a41);
  a31=(a31+a41);
  if (res[0]!=0) res[0][29]=a31;
  a31=(a25*a6);
  a5=(a5*a34);
  a31=(a31+a5);
  if (res[0]!=0) res[0][30]=a31;
  a31=(a9*a6);
  a50=(a50*a34);
  a31=(a31+a50);
  if (res[0]!=0) res[0][31]=a31;
  a6=(a21*a6);
  a1=(a1*a34);
  a6=(a6+a1);
  if (res[0]!=0) res[0][32]=a6;
  a25=(a15*a25);
  if (res[0]!=0) res[0][33]=a25;
  a9=(a15*a9);
  if (res[0]!=0) res[0][34]=a9;
  a15=(a15*a21);
  if (res[0]!=0) res[0][35]=a15;
  if (res[0]!=0) res[0][36]=a29;
  if (res[0]!=0) res[0][37]=a29;
  if (res[0]!=0) res[0][38]=a29;
  if (res[0]!=0) res[0][39]=a29;
  if (res[0]!=0) res[0][40]=a29;
  if (res[0]!=0) res[0][41]=a29;
  if (res[0]!=0) res[0][42]=a29;
  if (res[0]!=0) res[0][43]=a29;
  if (res[0]!=0) res[0][44]=a29;
  if (res[0]!=0) res[0][45]=a29;
  if (res[0]!=0) res[0][46]=a29;
  if (res[0]!=0) res[0][47]=a29;
  if (res[0]!=0) res[0][48]=a29;
  if (res[0]!=0) res[0][49]=a29;
  if (res[0]!=0) res[0][50]=a29;
  if (res[0]!=0) res[0][51]=a29;
  if (res[0]!=0) res[0][52]=a29;
  if (res[0]!=0) res[0][53]=a29;
  return 0;
}

extern "C" CASADI_SYMBOL_EXPORT int comp_foot_jacob_2(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

extern "C" CASADI_SYMBOL_EXPORT int comp_foot_jacob_2_alloc_mem(void) {
  return 0;
}

extern "C" CASADI_SYMBOL_EXPORT int comp_foot_jacob_2_init_mem(int mem) {
  return 0;
}

extern "C" CASADI_SYMBOL_EXPORT void comp_foot_jacob_2_free_mem(int mem) {
}

extern "C" CASADI_SYMBOL_EXPORT int comp_foot_jacob_2_checkout(void) {
  return 0;
}

extern "C" CASADI_SYMBOL_EXPORT void comp_foot_jacob_2_release(int mem) {
}

extern "C" CASADI_SYMBOL_EXPORT void comp_foot_jacob_2_incref(void) {
}

extern "C" CASADI_SYMBOL_EXPORT void comp_foot_jacob_2_decref(void) {
}

extern "C" CASADI_SYMBOL_EXPORT casadi_int comp_foot_jacob_2_n_in(void) { return 3;}

extern "C" CASADI_SYMBOL_EXPORT casadi_int comp_foot_jacob_2_n_out(void) { return 1;}

extern "C" CASADI_SYMBOL_EXPORT casadi_real comp_foot_jacob_2_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

extern "C" CASADI_SYMBOL_EXPORT const char* comp_foot_jacob_2_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    default: return 0;
  }
}

extern "C" CASADI_SYMBOL_EXPORT const char* comp_foot_jacob_2_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

extern "C" CASADI_SYMBOL_EXPORT const casadi_int* comp_foot_jacob_2_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s0;
    default: return 0;
  }
}

extern "C" CASADI_SYMBOL_EXPORT const casadi_int* comp_foot_jacob_2_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

extern "C" CASADI_SYMBOL_EXPORT int comp_foot_jacob_2_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


