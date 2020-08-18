/*Copyright 2019 DScribe developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <map>
#include <set>
#include "soapGTO.h"
#include "celllist.h"

#define PI2 9.86960440108936
#define PI 3.14159265359
#define PIHalf 1.57079632679490
//===========================================================
inline int getCrosNum(int n)
{
  return n*(n+1)/2;
}
//===========================================================
inline void getReIm2(double* x, double* y, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*x[i]-y[i]*y[i];
    c3[2*i+1] = 2*y[i]*x[i];
  }
}
//===========================================================
inline void getReIm3(double* x, double* y, double* c2, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = x[i]*c2[2*i] - y[i]*c2[2*i + 1];
    c3[2*i+1] = x[i]*c2[2*i+1] + y[i]*c2[2*i  ];
  }
}
//===========================================================
inline void getMulReIm(double* c1, double* c2, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = c1[2*i  ]*c2[2*i  ] - c1[2*i+1]*c2[2*i+1];
    c3[2*i+1] = c1[2*i  ]*c2[2*i+1] + c1[2*i+1]*c2[2*i  ];
  }
}
//===========================================================
inline void getMulDouble(double* c1, double* c3, int Asize){
  for(int i = 0; i < Asize; i++){
    c3[2*i  ] = c1[2*i]*c1[2*i] - c1[2*i+1]*c1[2*i+1];
    c3[2*i+1] = 2*c1[2*i]*c1[2*i+1];
  }
}
//================================================================
inline void getDeltas(double* x, double* y, double* z, const py::array_t<double> &positions, const double ix, const double iy, const double iz, const vector<int> &indices){

    int count = 0;
    auto pos = positions.unchecked<2>();
    for (const int &idx : indices) {
        x[count] = pos(idx, 0) - ix;
        y[count] = pos(idx, 1) - iy;
        z[count] = pos(idx, 2) - iz;
        count++;
    };
}

//================================================================
inline void getRsZs(double* x, double* y, double* z,double* r2,double* r4,double* r6,double* r8,double* z2,double* z4,double* z6,double* z8, int size){
  for(int i = 0; i < size; i++){
    r2[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
    r4[i] = r2[i]*r2[i]; r6[i] = r2[i]*r4[i]; r8[i] = r4[i]*r4[i];
    z2[i] = z[i]*z[i]; z4[i] = z2[i]*z2[i]; z6[i] = z2[i]*z4[i]; z8[i] = z4[i]*z4[i];
  }
}
//================================================================
void getAlphaBeta(double* aOa, double* bOa, double* alphas, double* betas, int Ns,int lMax, double oOeta, double oOeta3O2){

  int  NsNs = Ns*Ns;
  double  oneO1alpha;  double  oneO1alpha2; double  oneO1alpha3;
  double  oneO1alpha4; double  oneO1alpha5; double  oneO1alpha6;
  double  oneO1alpha7; double  oneO1alpha8; double  oneO1alpha9;
  double  oneO1alpha10;
  double  oneO1alphaSqrt;// = (double*) malloc(Ns*sizeof(double));
  double  oneO1alphaSqrtX;

  for(int k = 0; k < Ns; k++){
    oneO1alpha = 1.0/(1.0 + oOeta*alphas[k]);
    oneO1alphaSqrt = sqrt(oneO1alpha);
    aOa[k] = -alphas[k]*oneO1alpha; //got alpha_0k
    oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha;
    for(int n = 0; n < Ns; n++){ bOa[n*Ns + k] = oOeta3O2*betas[n*Ns + k]*oneO1alphaSqrtX;} // got beta_0nk
  }
  if(lMax > 0){
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[Ns + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[Ns + k] = -alphas[Ns + k]*oneO1alpha; //got alpha_1k
      oneO1alpha2 = oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha2;
      for(int n = 0; n < Ns; n++){ bOa[NsNs + n*Ns + k] = oOeta3O2*betas[NsNs + n*Ns + k]*oneO1alphaSqrtX;} // got beta_1nk
    }
  } if(lMax > 1){
    int shift1 = 2*Ns; int shift2 = 2*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_2k
      oneO1alpha3 = oneO1alpha*oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha3;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_2nk
    }
  } if(lMax > 2){
    int shift1 = 3*Ns; int shift2 = 3*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_3k
      oneO1alpha4 = oneO1alpha*oneO1alpha*oneO1alpha*oneO1alpha; oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha4;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_3nk
    }
  } if(lMax > 3){
    int shift1 = 4*Ns; int shift2 = 4*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_4k
      oneO1alpha5 = pow(oneO1alpha,5); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha5;
      for(int n = 0; n < Ns; n++){ bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_4nk
    }
  } if(lMax > 4){
    int shift1 = 5*Ns; int shift2 = 5*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_5k
      oneO1alpha6 = pow(oneO1alpha,6); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha6;
      for(int n = 0; n < Ns; n++){ bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_5nk
    }
  } if(lMax > 5){
    int shift1 = 6*Ns; int shift2 = 6*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_6k
      oneO1alpha7 = pow(oneO1alpha,7); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha7;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_6nk
    }
  } if(lMax > 6){
    int shift1 = 7*Ns; int shift2 = 7*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_7k
      oneO1alpha8 = pow(oneO1alpha,8); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha8;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_7nk
    }
  } if(lMax > 7){
    int shift1 = 8*Ns; int shift2 = 8*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_8k
      oneO1alpha9 = pow(oneO1alpha,9); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha9;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_8nk
    }
  }if(lMax > 8){
    int shift1 = 9*Ns; int shift2 = 9*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_9k
      oneO1alpha10 = pow(oneO1alpha,10); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha10;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_9nk
    }
  }
}
//================================================================
void getCfactors(double* preCoef, int Asize, double* x,double* x2,double* x4,double* x6,double* x8,double* x10,double* x12,double* x14,double* x16,double* x18, double* y,double* y2,double* y4,double* y6,double* y8,double* y10,double* y12,double* y14,double* y16,double* y18, double* z, double* z2, double* z4, double* z6, double* z8,double* z10, double* z12, double* z14, double* z16, double* z18, double* r2,double* r3, double* r4, double* r5, double* r6, double* r7, double* r8, double* r9, double* ReIm2, double* ReIm3, double* ReIm4, double* ReIm5, double* ReIm6, double* ReIm7, double* ReIm8, double* ReIm9,int totalAN, int lMax, int t2, int t3, int t4, int t5, int t6, int t7, int t8, int t9, int t10, int t11, int t12, int t13, int t14, int t15, int t16, int t17, int t18, int t19, int t20, int t21, int t22, int t23, int t24, int t25, int t26, int t27, int t28, int t29, int t30, int t31, int t32, int t33, int t34, int t35, int t36, int t37, int t38, int t39, int t40, int t41, int t42, int t43, int t44, int t45, int t46, int t47, int t48, int t49, int t50, int t51, int t52, int t53, int t54, int t55, int t56, int t57, int t58, int t59, int t60, int t61, int t62, int t63, int t64, int t65, int t66, int t67, int t68, int t69, int t70, int t71, int t72, int t73, int t74, int t75, int t76, int t77, int t78, int t79, int t80, int t81, int t82, int t83, int t84, int t85, int t86, int t87, int t88, int t89, int t90, int t91, int t92, int t93, int t94, int t95, int t96, int t97, int t98, int t99){
  double c20c;double c30c;double c31c;double c40c;double c41c;double c42c;
  double c50c;double c51c;double c52c;double c53c;double c60c;double c61c;
  double c62c;double c63c;double c64c;double c70c;double c71c;double c72c;
  double c73c;double c74c;double c75c;double c80c;double c81c;double c82c;
  double c83c;double c84c;double c85c;double c86c;double c90c;double c91c;
  double c92c;double c93c;double c94c;double c95c;double c96c;double c97c;

  getReIm2(x, y, ReIm2, Asize);
  getReIm3(x, y, ReIm2, ReIm3, Asize);
  getMulDouble(ReIm2, ReIm4, Asize);
  getMulReIm(ReIm2, ReIm3, ReIm5, Asize);
  getMulDouble(ReIm3, ReIm6, Asize);
  getMulReIm(ReIm3, ReIm4, ReIm7, Asize);
  getMulDouble(ReIm4, ReIm8, Asize);
  getMulReIm(ReIm4, ReIm5, ReIm9, Asize);
  int i2;

  for (int i = 0; i < Asize; i++) {
    i2 = 2*i;
    c20c=3*z2[i]-r2[i];
    if (lMax > 2){
      c30c=5*z2[i]-3*r2[i];
      c31c=5*z2[i]-r2[i];
    }
    if (lMax > 3){
      c40c=35*z4[i]-30*z2[i]*r2[i]+3*r4[i];
      c41c=7*z2[i]-3*r2[i];
      c42c=7*z2[i]-r2[i];
    }
    if (lMax > 4){
      c50c=63*z4[i]-70*z2[i]*r2[i]+15*r4[i];
      c51c=21*z4[i]-14*z2[i]*r2[i]+r4[i];
      c52c=3*z2[i]-r2[i];
      c53c=9*z2[i]-r2[i];
    }
    if (lMax > 5){
      c60c=231*z6[i] - 315*z4[i]*r2[i] + 105*z2[i]*r4[i] - 5*r6[i];
      c61c=33*z4[i] - 30*z2[i]*r2[i] + 5*r4[i];
      c62c=33*z4[i] - 18*z2[i]*r2[i] + r4[i];
      c63c=11*z2[i] - 3*r2[i];
      c64c=11*z2[i] - r2[i];
    }
    if (lMax > 6){
      c70c=429*z6[i]-693*z4[i]*r2[i]+315*z2[i]*r4[i]-35*r6[i];
      c71c=429*z6[i]-495*z4[i]*r2[i]+135*z2[i]*r4[i]-5*r6[i];
      c72c=143*z4[i]-110*z2[i]*r2[i]+15*r4[i];
      c73c=143*z4[i]-66*z2[i]*r2[i]+3*r4[i];
      c74c=13*z2[i]-3*r2[i];
      c75c=13*z2[i]-r2[i];
    }
    if (lMax > 7){
      c80c=6435*z8[i]-12012*z6[i]*r2[i]+6930*z4[i]*r4[i]-1260*z2[i]*r6[i]+35*r8[i];
      c81c=715*z6[i]-1001*z4[i]*r2[i]+385*z2[i]*r4[i]-35*r6[i];
      c82c=143*z6[i]-143*z4[i]*r2[i]+33*z2[i]*r4[i]-r6[i];
      c83c=39*z4[i]-26*z2[i]*r2[i]+3*r4[i];
      c84c=65*z4[i]-26*z2[i]*r2[i]+r4[i];
      c85c=5*z2[i]-r2[i];
      c86c=15*z2[i]-r2[i];
    }
    if (lMax > 8){
      c90c=12155*z8[i]-25740*z6[i]*r2[i]+18018*z4[i]*r4[i] -4620*z2[i]*r6[i]+315*r8[i];
      c91c=2431*z8[i]-4004*z6[i]*r2[i]+2002*z4[i]*r4[i]-308*z2[i]*r6[i] + 7*r8[i];
      c92c=221*z6[i]-273*z4[i]*r2[i]+91*z2[i]*r4[i]-7*r6[i];
      c93c=221*z6[i]-195*z4[i]*r2[i]+39*z2[i]*r4[i]-r6[i];
      c94c=17*z4[i]-10*z2[i]*r2[i]+r4[i];
      c95c=85*z4[i]-30*z2[i]*r2[i]+r4[i];
      c96c=17*z2[i]-3*r2[i];
      c97c=17*z2[i]-r2[i];
    }
    /*c20  */  preCoef[        +i] = c20c;
    /*c21Re*/  preCoef[totalAN+i] = z[i]*x[i];
    /*c21Im*/  preCoef[t2+i] = z[i]*y[i];
    /*c22Re*/  preCoef[t3+i] =      ReIm2[2*i];
    /*c22Im*/  preCoef[t4+i] =      ReIm2[i2+1];
    if (lMax > 2){
      /*c30  */  preCoef[t5+i] = c30c*z[i];
      /*c31Re*/  preCoef[t6+i] =       x[i]*c31c;
      /*c31Im*/  preCoef[t7+i] =       y[i]*c31c;
      /*c32Re*/  preCoef[t8+i] = z[i]*ReIm2[i2];
      /*c32Im*/  preCoef[t9+i] = z[i]*ReIm2[i2+1];
      /*c33Re*/  preCoef[t10+i] =      ReIm3[i2  ];
      /*c33Im*/  preCoef[t11+i] =      ReIm3[i2+1];
    }
    if (lMax > 3){
      /*c40  */  preCoef[t12+i] = c40c;
      /*c41Re*/  preCoef[t13+i] = z[i]*x[i]*c41c;
      /*c41Im*/  preCoef[t14+i] = z[i]*y[i]*c41c;
      /*c42Re*/  preCoef[t15+i] =      ReIm2[i2  ]*c42c;
      /*c42Im*/  preCoef[t16+i] =      ReIm2[i2+1]*c42c;
      /*c43Re*/  preCoef[t17+i] = z[i]*ReIm3[i2  ];
      /*c43Im*/  preCoef[t18+i] = z[i]*ReIm3[i2+1];
      /*c44Re*/  preCoef[t19+i] =      ReIm4[i2  ];
      /*c44Im*/  preCoef[t20+i] =      ReIm4[i2+1];
    }
    if (lMax > 4){
      /*c50  */  preCoef[t21+i] = c50c*z[i];
      /*c51Re*/  preCoef[t22+i] =      x[i]*c51c;
      /*c51Im*/  preCoef[t23+i] =      y[i]*c51c;
      /*c52Re*/  preCoef[t24+i] = z[i]*ReIm2[i2  ]*c52c;
      /*c52Im*/  preCoef[t25+i] = z[i]*ReIm2[i2+1]*c52c;
      /*c53Re*/  preCoef[t26+i] =      ReIm3[i2  ]*c53c;
      /*c53Im*/  preCoef[t27+i] =      ReIm3[i2+1]*c53c;
      /*c54Re*/  preCoef[t28+i] = z[i]*ReIm4[i2  ];
      /*c54Im*/  preCoef[t29+i] = z[i]*ReIm4[i2+1];
      /*c55Re*/  preCoef[t30+i] =      ReIm5[i2  ];
      /*c55Im*/  preCoef[t31+i] =      ReIm5[i2+1];
    }
    if (lMax > 5){
      /*c60  */  preCoef[t32+i] = c60c;
      /*c61Re*/  preCoef[t33+i] = z[i]*x[i]*c61c;
      /*c61Im*/  preCoef[t34+i] = z[i]*y[i]*c61c;
      /*c62Re*/  preCoef[t35+i] =      ReIm2[i2  ]*c62c;
      /*c62Im*/  preCoef[t36+i] =      ReIm2[i2+1]*c62c;
      /*c63Re*/  preCoef[t37+i] = z[i]*ReIm3[i2  ]*c63c;
      /*c63Im*/  preCoef[t38+i] = z[i]*ReIm3[i2+1]*c63c;
      /*c64Re*/  preCoef[t39+i] =      ReIm4[i2  ]*c64c;
      /*c64Im*/  preCoef[t40+i] =      ReIm4[i2+1]*c64c;
      /*c65Re*/  preCoef[t41+i] = z[i]*ReIm5[i2  ];
      /*c65Im*/  preCoef[t42+i] = z[i]*ReIm5[i2+1];
      /*c66Re*/  preCoef[t43+i] =      ReIm6[i2  ];
      /*c66Im*/  preCoef[t44+i] =      ReIm6[i2+1];
    }
    if (lMax > 6){
      /*c70  */  preCoef[t45+i] = c70c*z[i];
      /*c71Re*/  preCoef[t46+i] = x[i]*c71c;
      /*c71Im*/  preCoef[t47+i] = y[i]*c71c;
      /*c72Re*/  preCoef[t48+i] = z[i]*ReIm2[i2  ]*c72c;
      /*c72Im*/  preCoef[t49+i] = z[i]*ReIm2[i2+1]*c72c;
      /*c73Re*/  preCoef[t50+i] =      ReIm3[i2  ]*c73c;
      /*c73Im*/  preCoef[t51+i] =      ReIm3[i2+1]*c73c;
      /*c74Re*/  preCoef[t52+i] = z[i]*ReIm4[i2  ]*c74c;
      /*c74Im*/  preCoef[t53+i] = z[i]*ReIm4[i2+1]*c74c;
      /*c75Re*/  preCoef[t54+i] =      ReIm5[i2  ]*c75c;
      /*c75Im*/  preCoef[t55+i] =      ReIm5[i2+1]*c75c;
      /*c76Re*/  preCoef[t56+i] = z[i]*ReIm6[i2  ];
      /*c76Im*/  preCoef[t57+i] = z[i]*ReIm6[i2+1];
      /*c77Re*/  preCoef[t58+i] =      ReIm7[i2  ];
      /*c77Im*/  preCoef[t59+i] =      ReIm7[i2+1];
    }
    if (lMax > 7){
      /*c80  */  preCoef[t60+i] = c80c;
      /*c81Re*/  preCoef[t61+i] = z[i]*x[i]*c81c;
      /*c81Im*/  preCoef[t62+i] = z[i]*y[i]*c81c;
      /*c82Re*/  preCoef[t63+i] =      ReIm2[i2  ]*c82c;
      /*c82Im*/  preCoef[t64+i] =      ReIm2[i2+1]*c82c;
      /*c83Re*/  preCoef[t65+i] = z[i]*ReIm3[i2  ]*c83c;
      /*c83Im*/  preCoef[t66+i] = z[i]*ReIm3[i2+1]*c83c;
      /*c84Re*/  preCoef[t67+i] =      ReIm4[i2  ]*c84c;
      /*c84Im*/  preCoef[t68+i] =      ReIm4[i2+1]*c84c;
      /*c85Re*/  preCoef[t69+i] = z[i]*ReIm5[i2  ]*c85c;
      /*c85Im*/  preCoef[t70+i] = z[i]*ReIm5[i2+1]*c85c;
      /*c86Re*/  preCoef[t71+i] =      ReIm6[i2  ]*c86c;
      /*c86Im*/  preCoef[t72+i] =      ReIm6[i2+1]*c86c;
      /*c87Re*/  preCoef[t73+i] = z[i]*ReIm7[i2  ];
      /*c87Im*/  preCoef[t74+i] = z[i]*ReIm7[i2+1];
      /*c88Re*/  preCoef[t75+i] =      ReIm8[i2  ];
      /*c88Im*/  preCoef[t76+i] =      ReIm8[i2+1];
    }
    if (lMax > 8){ // SOAPGTO to L=<9
      /*c90  */  preCoef[t77+i] = c90c*z[i];
      /*c91Re*/  preCoef[t78+i] = x[i]*c91c;
      /*c91Im*/  preCoef[t79+i] = y[i]*c91c;
      /*c92Re*/  preCoef[t80+i] = z[i]*ReIm2[i2  ]*c92c;
      /*c92Im*/  preCoef[t81+i] = z[i]*ReIm2[i2+1]*c92c;
      /*c93Re*/  preCoef[t82+i] =      ReIm3[i2  ]*c93c;
      /*c93Im*/  preCoef[t83+i] =      ReIm3[i2+1]*c93c;
      /*c94Re*/  preCoef[t84+i] = z[i]*ReIm4[i2  ]*c94c;
      /*c94Im*/  preCoef[t85+i] = z[i]*ReIm4[i2+1]*c94c;
      /*c95Re*/  preCoef[t86+i] =      ReIm5[i2  ]*c95c;
      /*c95Im*/  preCoef[t87+i] =      ReIm5[i2+1]*c95c;
      /*c96Re*/  preCoef[t88+i] = z[i]*ReIm6[i2  ]*c96c;
      /*c96Im*/  preCoef[t89+i] = z[i]*ReIm6[i2+1]*c96c;
      /*c97Re*/  preCoef[t90+i] =      ReIm7[i2  ]*c97c;
      /*c97Im*/  preCoef[t91+i] =      ReIm7[i2+1]*c97c;
      /*c98Re*/  preCoef[t92+i] = z[i]*ReIm8[i2  ];
      /*c98Im*/  preCoef[t93+i] = z[i]*ReIm8[i2+1];
      /*c99Re*/  preCoef[t94+i] =      ReIm9[i2  ];
      /*c99Im*/  preCoef[t95+i] =      ReIm9[i2+1];
    }
    if (lMax > 9){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t96+i] =  1.53479023644398*x[i]*y[i]*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t97+i] = 3.43189529989171*y[i]*z[i]*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t98+i] = -4.45381546176335*x[i]*y[i]*(x2[i] + y2[i] - 18.0*z2[i])*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t99+i] = 1.36369691122981*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 16.0*z2[i])*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t100+i] = 0.330745082725238*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(323.0*z4[i] - 102.0*z2[i]*r[i] + 3.0*r2[i]);
                 preCoef[t101+i] = 0.295827395278969*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(323.0*z4[i] - 170.0*z2[i]*r[i] + 15.0*r2[i]);
                 preCoef[t102+i] = 1.87097672671297*x[i]*y[i]*(x2[i] - y2[i])*(323.0*z6[i] - 255.0*z4[i]*r[i] + 45.0*z2[i]*r2[i] - r3[i]);
                 preCoef[t103+i] = 0.661490165450475*y[i]*z[i]*(3.0*x2[i] - y2[i])*(323.0*z6[i] - 357.0*z4[i]*r[i] + 105.0*z2[i]*r2[i] - 7.0*r3[i]);
                 preCoef[t104+i] = 0.129728894680065*x[i]*y[i]*(4199.0*z8[i] - 6188.0*z6[i]*r[i] + 2730.0*z4[i]*r2[i] - 364.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t105+i] = 0.0748990122652082*y[i]*z[i]*(4199.0*z8[i] - 7956.0*z6[i]*r[i] + 4914.0*z4[i]*r2[i] - 1092.0*z2[i]*r3[i] + 63.0*r4[i]);
                 preCoef[t106+i] = 233.240148813258*z10[i] - 552.410878768242*z8[i]*r[i] + 454.926606044435*z6[i]*r2[i] - 151.642202014812*z4[i]*r3[i] + 17.4971771555552*z2[i]*r4[i] - 0.318130493737367*r5[i];
                 preCoef[t107+i] = 0.0748990122652082*x[i]*z[i]*(4199.0*z8[i] - 7956.0*z6[i]*r[i] + 4914.0*z4[i]*r2[i] - 1092.0*z2[i]*r3[i] + 63.0*r4[i]);
                 preCoef[t108+i] = 0.0648644473400325*(x2[i] - y2[i])*(4199.0*z8[i] - 6188.0*z6[i]*r[i] + 2730.0*z4[i]*r2[i] - 364.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t109+i] = 0.661490165450475*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(323.0*z6[i] - 357.0*z4[i]*r[i] + 105.0*z2[i]*r2[i] - 7.0*r3[i]);
                 preCoef[t110+i] = 0.467744181678242*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(323.0*z6[i] - 255.0*z4[i]*r[i] + 45.0*z2[i]*r2[i] - r3[i]);
                 preCoef[t111+i] = 0.295827395278969*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(323.0*z4[i] - 170.0*z2[i]*r[i] + 15.0*r2[i]);
                 preCoef[t112+i] = 0.165372541362619*(323.0*z4[i] - 102.0*z2[i]*r[i] + 3.0*r2[i])*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t113+i] = 1.36369691122981*x[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 16.0*z2[i])*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i]);
                 preCoef[t114+i] = -0.556726932720418*(x2[i] + y2[i] - 18.0*z2[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t115+i] = 3.43189529989171*x[i]*z[i]*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t116+i] = 0.76739511822199*x10[i] - 34.5327803199895*x8[i]*y2[i] + 161.152974826618*x6[i]*y4[i] - 161.152974826618*x4[i]*y6[i] + 34.5327803199895*x2[i]*y8[i] - 0.76739511822199*y10[i];
    }
    if (lMax > 10){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t117+i] = 0.784642105787197*y[i]*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t118+i] = 7.36059539761062*x[i]*y[i]*z[i]*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t119+i] = -0.567882263783437*y[i]*(x2[i] + y2[i] - 20.0*z2[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t120+i] = 35.1903768038371*x[i]*y[i]*z[i]*(-x2[i] - y2[i] + 6.0*z2[i])*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t121+i] = 0.504576632477118*y[i]*(133.0*z4[i] - 38.0*z2[i]*r[i] + r2[i])*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t122+i] = 0.638244565090152*x[i]*y[i]*z[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(399.0*z4[i] - 190.0*z2[i]*r[i] + 15.0*r2[i]);
                 preCoef[t123+i] = 0.0947934431913346*y[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(2261.0*z6[i] - 1615.0*z4[i]*r[i] + 255.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t124+i] = 4.0127980256608*x[i]*y[i]*z[i]*(x2[i] - y2[i])*(323.0*z6[i] - 323.0*z4[i]*r[i] + 85.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t125+i] = 0.457895832784712*y[i]*(3.0*x2[i] - y2[i])*(969.0*z8[i] - 1292.0*z6[i]*r[i] + 510.0*z4[i]*r2[i] - 60.0*z2[i]*r3[i] + r4[i]);
                 preCoef[t126+i] = 0.489511235746264*x[i]*y[i]*z[i]*(2261.0*z8[i] - 3876.0*z6[i]*r[i] + 2142.0*z4[i]*r2[i] - 420.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t127+i] = 0.0214664877426077*y[i]*(29393.0*z10[i] - 62985.0*z8[i]*r[i] + 46410.0*z6[i]*r2[i] - 13650.0*z4[i]*r3[i] + 1365.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t128+i] = 0.00528468396465431*z[i]*(88179.0*z10[i] - 230945.0*z8[i]*r[i] + 218790.0*z6[i]*r2[i] - 90090.0*z4[i]*r3[i] + 15015.0*z2[i]*r4[i] - 693.0*r5[i]);
                 preCoef[t129+i] = 0.0214664877426077*x[i]*(29393.0*z10[i] - 62985.0*z8[i]*r[i] + 46410.0*z6[i]*r2[i] - 13650.0*z4[i]*r3[i] + 1365.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t130+i] = 0.244755617873132*z[i]*(x2[i] - y2[i])*(2261.0*z8[i] - 3876.0*z6[i]*r[i] + 2142.0*z4[i]*r2[i] - 420.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t131+i] = 0.457895832784712*x[i]*(x2[i] - 3.0*y2[i])*(969.0*z8[i] - 1292.0*z6[i]*r[i] + 510.0*z4[i]*r2[i] - 60.0*z2[i]*r3[i] + r4[i]);
                 preCoef[t132+i] = 1.0031995064152*z[i]*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(323.0*z6[i] - 323.0*z4[i]*r[i] + 85.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t133+i] = 0.0947934431913346*x[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(2261.0*z6[i] - 1615.0*z4[i]*r[i] + 255.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t134+i] = 0.319122282545076*z[i]*(399.0*z4[i] - 190.0*z2[i]*r[i] + 15.0*r2[i])*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t135+i] = 0.504576632477118*x[i]*(133.0*z4[i] - 38.0*z2[i]*r[i] + r2[i])*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i]);
                 preCoef[t136+i] = 4.39879710047964*z[i]*(-x2[i] - y2[i] + 6.0*z2[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t137+i] = -0.567882263783437*x[i]*(x2[i] + y2[i] - 20.0*z2[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t138+i] = 3.68029769880531*z[i]*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t139+i] = 0.784642105787197*x[i]*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
    }
    if (lMax > 11){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t140+i] = 3.20328798313589*x[i]*y[i]*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t141+i] = 3.92321052893598*y[i]*z[i]*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t142+i] = -1.15689166958762*x[i]*y[i]*(x2[i] + y2[i] - 22.0*z2[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t143+i] = 1.56643872562221*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 20.0*z2[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t144+i] = 4.10189944667082*x[i]*y[i]*(161.0*z4[i] - 42.0*z2[i]*r[i] + r2[i])*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t145+i] = 1.0254748616677*y[i]*z[i]*(161.0*z4[i] - 70.0*z2[i]*r[i] + 5.0*r2[i])*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t146+i] = 0.192089041114334*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(3059.0*z6[i] - 1995.0*z4[i]*r[i] + 285.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t147+i] = 1.07809706940566*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(437.0*z6[i] - 399.0*z4[i]*r[i] + 95.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t148+i] = 0.369784244098711*x[i]*y[i]*(x2[i] - y2[i])*(7429.0*z8[i] - 9044.0*z6[i]*r[i] + 3230.0*z4[i]*r2[i] - 340.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t149+i] = 0.12326141469957*y[i]*z[i]*(3.0*x2[i] - y2[i])*(7429.0*z8[i] - 11628.0*z6[i]*r[i] + 5814.0*z4[i]*r2[i] - 1020.0*z2[i]*r3[i] + 45.0*r4[i]);
                 preCoef[t150+i] = 0.301927570987541*x[i]*y[i]*(7429.0*z10[i] - 14535.0*z8[i]*r[i] + 9690.0*z6[i]*r2[i] - 2550.0*z4[i]*r3[i] + 225.0*z2[i]*r4[i] - 3.0*r5[i]);
                 preCoef[t151+i] = 0.0243300170174164*y[i]*z[i]*(52003.0*z10[i] - 124355.0*z8[i]*r[i] + 106590.0*z6[i]*r2[i] - 39270.0*z4[i]*r3[i] + 5775.0*z2[i]*r4[i] - 231.0*r5[i]);
                 preCoef[t152+i] = 931.186918632914*z12[i] - 2672.1015925988*z10[i]*r[i] + 2862.96599207014*z8[i]*r2[i] - 1406.36925926252*z6[i]*r3[i] + 310.228513072616*z4[i]*r4[i] - 24.8182810458093*z2[i]*r5[i] + 0.318183090330888*r6[i];
                 preCoef[t153+i] = 0.0243300170174164*x[i]*z[i]*(52003.0*z10[i] - 124355.0*z8[i]*r[i] + 106590.0*z6[i]*r2[i] - 39270.0*z4[i]*r3[i] + 5775.0*z2[i]*r4[i] - 231.0*r5[i]);
                 preCoef[t154+i] = 0.150963785493771*(x2[i] - y2[i])*(7429.0*z10[i] - 14535.0*z8[i]*r[i] + 9690.0*z6[i]*r2[i] - 2550.0*z4[i]*r3[i] + 225.0*z2[i]*r4[i] - 3.0*r5[i]);
                 preCoef[t155+i] = 0.12326141469957*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(7429.0*z8[i] - 11628.0*z6[i]*r[i] + 5814.0*z4[i]*r2[i] - 1020.0*z2[i]*r3[i] + 45.0*r4[i]);
                 preCoef[t156+i] = 0.0924460610246778*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(7429.0*z8[i] - 9044.0*z6[i]*r[i] + 3230.0*z4[i]*r2[i] - 340.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t157+i] = 1.07809706940566*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(437.0*z6[i] - 399.0*z4[i]*r[i] + 95.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t158+i] = 0.0960445205571672*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(3059.0*z6[i] - 1995.0*z4[i]*r[i] + 285.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t159+i] = 1.0254748616677*x[i]*z[i]*(161.0*z4[i] - 70.0*z2[i]*r[i] + 5.0*r2[i])*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i]);
                 preCoef[t160+i] = 0.512737430833852*(161.0*z4[i] - 42.0*z2[i]*r[i] + r2[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t161+i] = 1.56643872562221*x[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 20.0*z2[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t162+i] = -0.57844583479381*(x2[i] + y2[i] - 22.0*z2[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t163+i] = 3.92321052893598*x[i]*z[i]*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t164+i] = 0.800821995783972*x12[i] - 52.8542517217421*x10[i]*y2[i] + 396.406887913066*x8[i]*y4[i] - 739.95952410439*x6[i]*y6[i] + 396.406887913066*x4[i]*y8[i] - 52.8542517217421*x2[i]*y10[i] + 0.800821995783972*y12[i];
    }
    if (lMax > 12){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t165+i] = 0.816077118837628*y[i]*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t166+i] = 16.6447726141986*x[i]*y[i]*z[i]*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t167+i] = -0.588481579340398*y[i]*(x2[i] + y2[i] - 24.0*z2[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t168+i] = 3.32895452283972*x[i]*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 22.0*z2[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t169+i] = 0.173533750438143*y[i]*(575.0*z4[i] - 138.0*z2[i]*r[i] + 3.0*r2[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t170+i] = 14.560298633254*x[i]*y[i]*z[i]*(115.0*z4[i] - 46.0*z2[i]*r[i] + 3.0*r2[i])*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i]);
                 preCoef[t171+i] = 0.486425436917406*y[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(805.0*z6[i] - 483.0*z4[i]*r[i] + 63.0*z2[i]*r2[i] - r3[i]);
                 preCoef[t172+i] = 2.30218535466597*x[i]*y[i]*z[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(575.0*z6[i] - 483.0*z4[i]*r[i] + 105.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t173+i] = 0.0933659449850856*y[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(10925.0*z8[i] - 12236.0*z6[i]*r[i] + 3990.0*z4[i]*r2[i] - 380.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t174+i] = 0.528157542646753*x[i]*y[i]*z[i]*(x2[i] - y2[i])*(10925.0*z8[i] - 15732.0*z6[i]*r[i] + 7182.0*z4[i]*r2[i] - 1140.0*z2[i]*r3[i] + 45.0*r4[i]);
                 preCoef[t175+i] = 0.0506347929757152*y[i]*(3.0*x2[i] - y2[i])*(37145.0*z10[i] - 66861.0*z8[i]*r[i] + 40698.0*z6[i]*r2[i] - 9690.0*z4[i]*r3[i] + 765.0*z2[i]*r4[i] - 9.0*r5[i]);
                 preCoef[t176+i] = 0.122135716100197*x[i]*y[i]*z[i]*(37145.0*z10[i] - 81719.0*z8[i]*r[i] + 63954.0*z6[i]*r2[i] - 21318.0*z4[i]*r3[i] + 2805.0*z2[i]*r4[i] - 99.0*r5[i]);
                 preCoef[t177+i] = 0.0136551881840328*y[i]*(185725.0*z12[i] - 490314.0*z10[i]*r[i] + 479655.0*z8[i]*r2[i] - 213180.0*z6[i]*r3[i] + 42075.0*z4[i]*r4[i] - 2970.0*z2[i]*r5[i] + 33.0*r6[i]);
                 preCoef[t178+i] = 0.00143145267159059*z[i]*(1300075.0*z12[i] - 4056234.0*z10[i]*r[i] + 4849845.0*z8[i]*r2[i] - 2771340.0*z6[i]*r3[i] + 765765.0*z4[i]*r4[i] - 90090.0*z2[i]*r5[i] + 3003.0*r6[i]);
                 preCoef[t179+i] = 0.0136551881840328*x[i]*(185725.0*z12[i] - 490314.0*z10[i]*r[i] + 479655.0*z8[i]*r2[i] - 213180.0*z6[i]*r3[i] + 42075.0*z4[i]*r4[i] - 2970.0*z2[i]*r5[i] + 33.0*r6[i]);
                 preCoef[t180+i] = 0.0610678580500984*z[i]*(x2[i] - y2[i])*(37145.0*z10[i] - 81719.0*z8[i]*r[i] + 63954.0*z6[i]*r2[i] - 21318.0*z4[i]*r3[i] + 2805.0*z2[i]*r4[i] - 99.0*r5[i]);
                 preCoef[t181+i] = 0.0506347929757152*x[i]*(x2[i] - 3.0*y2[i])*(37145.0*z10[i] - 66861.0*z8[i]*r[i] + 40698.0*z6[i]*r2[i] - 9690.0*z4[i]*r3[i] + 765.0*z2[i]*r4[i] - 9.0*r5[i]);
                 preCoef[t182+i] = 0.132039385661688*z[i]*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(10925.0*z8[i] - 15732.0*z6[i]*r[i] + 7182.0*z4[i]*r2[i] - 1140.0*z2[i]*r3[i] + 45.0*r4[i]);
                 preCoef[t183+i] = 0.0933659449850856*x[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(10925.0*z8[i] - 12236.0*z6[i]*r[i] + 3990.0*z4[i]*r2[i] - 380.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t184+i] = 1.15109267733299*z[i]*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(575.0*z6[i] - 483.0*z4[i]*r[i] + 105.0*z2[i]*r2[i] - 5.0*r3[i]);
                 preCoef[t185+i] = 0.486425436917406*x[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(805.0*z6[i] - 483.0*z4[i]*r[i] + 63.0*z2[i]*r2[i] - r3[i]);
                 preCoef[t186+i] = 1.82003732915675*z[i]*(115.0*z4[i] - 46.0*z2[i]*r[i] + 3.0*r2[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t187+i] = 0.173533750438143*x[i]*(575.0*z4[i] - 138.0*z2[i]*r[i] + 3.0*r2[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t188+i] = 1.66447726141986*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 22.0*z2[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t189+i] = -0.588481579340398*x[i]*(x2[i] + y2[i] - 24.0*z2[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t190+i] = 4.16119315354964*z[i]*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t191+i] = 0.816077118837628*x[i]*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
    }
    if (lMax > 13){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t192+i] = 1.66104416612905*x[i]*y[i]*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t193+i] = 4.39470978027212*y[i]*z[i]*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t194+i] = -2.39217700650788*x[i]*y[i]*(x2[i] + y2[i] - 26.0*z2[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t195+i] = 5.2817838178514*y[i]*z[i]*(-x2[i] - y2[i] + 8.0*z2[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t196+i] = 1.05635676357028*x[i]*y[i]*(225.0*z4[i] - 50.0*z2[i]*r[i] + r2[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t197+i] = 1.92863476060198*y[i]*z[i]*(135.0*z4[i] - 50.0*z2[i]*r[i] + 3.0*r2[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t198+i] = 3.94023104498585*x[i]*y[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(1035.0*z6[i] - 575.0*z4[i]*r[i] + 69.0*z2[i]*r2[i] - r3[i]);
                 preCoef[t199+i] = 0.873160381394213*y[i]*z[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(1035.0*z6[i] - 805.0*z4[i]*r[i] + 161.0*z2[i]*r2[i] - 7.0*r3[i]);
                 preCoef[t200+i] = 0.943121003323134*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(3105.0*z8[i] - 3220.0*z6[i]*r[i] + 966.0*z4[i]*r2[i] - 84.0*z2[i]*r3[i] + r4[i]);
                 preCoef[t201+i] = 1.265329604663*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(1725.0*z8[i] - 2300.0*z6[i]*r[i] + 966.0*z4[i]*r2[i] - 140.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t202+i] = 1.83593315348819*x[i]*y[i]*(x2[i] - y2[i])*(6555.0*z10[i] - 10925.0*z8[i]*r[i] + 6118.0*z6[i]*r2[i] - 1330.0*z4[i]*r3[i] + 95.0*z2[i]*r4[i] - r5[i]);
                 preCoef[t203+i] = 0.0652370439174496*y[i]*z[i]*(3.0*x2[i] - y2[i])*(58995.0*z10[i] - 120175.0*z8[i]*r[i] + 86526.0*z6[i]*r2[i] - 26334.0*z4[i]*r3[i] + 3135.0*z2[i]*r4[i] - 99.0*r5[i]);
                 preCoef[t204+i] = 0.0274050400015396*x[i]*y[i]*(334305.0*z12[i] - 817190.0*z10[i]*r[i] + 735471.0*z8[i]*r2[i] - 298452.0*z6[i]*r3[i] + 53295.0*z4[i]*r4[i] - 3366.0*z2[i]*r5[i] + 33.0*r6[i]);
                 preCoef[t205+i] = 0.0152015810664143*y[i]*z[i]*(334305.0*z12[i] - 965770.0*z10[i]*r[i] + 1062347.0*z8[i]*r2[i] - 554268.0*z6[i]*r3[i] + 138567.0*z4[i]*r4[i] - 14586.0*z2[i]*r5[i] + 429.0*r6[i]);
                 preCoef[t206+i] = 3719.61718745389*z14[i] - 12536.487557715*z12[i]*r[i] + 16548.1635761838*z10[i]*r2[i] - 10792.2805931633*z8[i]*r3[i] + 3597.42686438778*z6[i]*r4[i] - 568.014768061228*z4[i]*r5[i] + 33.4126334153663*z2[i]*r6[i] - 0.318215556336822*r7[i];
                 preCoef[t207+i] = 0.0152015810664143*x[i]*z[i]*(334305.0*z12[i] - 965770.0*z10[i]*r[i] + 1062347.0*z8[i]*r2[i] - 554268.0*z6[i]*r3[i] + 138567.0*z4[i]*r4[i] - 14586.0*z2[i]*45[i] + 429.0*r6[i]);
                 preCoef[t208+i] = 0.0137025200007698*(x2[i] - y2[i])*(334305.0*z12[i] - 817190.0*z10[i]*r[i] + 735471.0*z8[i]*r2[i] - 298452.0*z6[i]*r3[i] + 53295.0*z4[i]*r4[i] - 3366.0*z2[i]*r5[i] + 33.0*r6[i]);
                 preCoef[t209+i] = 0.0652370439174496*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(58995.0*z10[i] - 120175.0*z8[i]*r[i] + 86526.0*z6[i]*r2[i] - 26334.0*z4[i]*r3[i] + 3135.0*z2[i]*r4[i] - 99.0*r5[i]);
                 preCoef[t210+i] = 0.458983288372048*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(6555.0*z10[i] - 10925.0*z8[i]*r[i] + 6118.0*z6[i]*r2[i] - 1330.0*z4[i]*r3[i] + 95.0*z2[i]*r4[i] - r5[i]);
                 preCoef[t211+i] = 1.265329604663*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(1725.0*z8[i] - 2300.0*z6[i]*r[i] + 966.0*z4[i]*r2[i] - 140.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t212+i] = 0.471560501661567*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(3105.0*z8[i] - 3220.0*z6[i]*r[i] + 966.0*z4[i]*r2[i] - 84.0*z2[i]*r3[i] + r4[i]);
                 preCoef[t213+i] = 0.873160381394213*x[i]*z[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(1035.0*z6[i] - 805.0*z4[i]*r[i] + 161.0*z2[i]*r2[i] - 7.0*r3[i]);
                 preCoef[t214+i] = 0.492528880623231*(1035.0*z6[i] - 575.0*z4[i]*r[i] + 69.0*z2[i]*r2[i] - r3[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t215+i] = 1.92863476060198*x[i]*z[i]*(135.0*z4[i] - 50.0*z2[i]*r[i] + 3.0*r2[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t216+i] = 0.52817838178514*(225.0*z4[i] - 50.0*z2[i]*r[i] + r2[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t217+i] = 5.2817838178514*x[i]*z[i]*(-x2[i] - y2[i] + 8.0*z2[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t218+i] = -0.598044251626971*(x2[i] + y2[i] - 26.0*z2[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t219+i] = 4.39470978027212*x[i]*z[i]*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t220+i] = 0.830522083064524*x14[i] - 75.5775095588717*x12[i]*y2[i] + 831.352605147589*x10[i]*y4[i] - 2494.05781544277*x8[i]*y6[i] + 2494.05781544277*x6[i]*y8[i] - 831.352605147589*x4[i]*y10[i] + 75.5775095588717*x2[i]*y12[i] - 0.830522083064524*y14[i];
    }
    if (lMax > 14){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t221+i] = 0.844250650857373*y[i]*(15.0*x14[i] - 455.0*x12[i]*y2[i] + 3003.0*x10[i]*y4[i] - 6435.0*x8[i]*y6[i] + 5005.0*x6[i]*y8[i] - 1365.0*x4[i]*y10[i] + 105.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t222+i] = 9.24830251326002*x[i]*y[i]*z[i]*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t223+i] = -0.60718080651189*y[i]*(x2[i] + y2[i] - 28.0*z2[i])*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t224+i] = 7.41987201697353*x[i]*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 26.0*z2[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t225+i] = 0.53548313829403*y[i]*(261.0*z4[i] - 54.0*z2[i]*r[i] + r2[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t226+i] = 2.44217885935126*x[i]*y[i]*z[i]*(261.0*z4[i] - 90.0*z2[i]*r[i] + 5.0*r2[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t227+i] = 0.49850767216857*y[i]*(1305.0*z6[i] - 675.0*z4[i]*r[i] + 75.0*z2[i]*r2[i] - r3[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t228+i] = 7.38445476455181*x[i]*y[i]*z[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(1305.0*z6[i] - 945.0*z4[i]*r[i] + 175.0*z2[i]*r2[i] - 7.0*r3[i]);
                 preCoef[t229+i] = 0.0680486534764293*y[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(30015.0*z8[i] - 28980.0*z6[i]*r[i] + 8050.0*z4[i]*r2[i] - 644.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t230+i] = 0.638352953401215*x[i]*y[i]*z[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(10005.0*z8[i] - 12420.0*z6[i]*r[i] + 4830.0*z4[i]*r2[i] - 644.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t231+i] = 0.462530657238986*y[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(10005.0*z10[i] - 15525.0*z8[i]*r[i] + 8050.0*z6[i]*r2[i] - 1610.0*z4[i]*r3[i] + 105.0*z2[i]*r4[i] - r5[i]);
                 preCoef[t232+i] = 2.49470484396446*x[i]*y[i]*z[i]*(x2[i] - y2[i])*(10005.0*z10[i] - 18975.0*z8[i]*r[i] + 12650.0*z6[i]*r2[i] - 3542.0*z4[i]*r3[i] + 385.0*z2[i]*r4[i] - 11.0*r5[i]);
                 preCoef[t233+i] = 0.0413039660894728*y[i]*(3.0*x2[i] - y2[i])*(190095.0*z12[i] - 432630.0*z10[i]*r[i] + 360525.0*z8[i]*r2[i] - 134596.0*z6[i]*r3[i] + 21945.0*z4[i]*r4[i] - 1254.0*z2[i]*r5[i] + 11.0*r6[i]);
                 preCoef[t234+i] = 0.0324014967813841*x[i]*y[i]*z[i]*(570285.0*z12[i] - 1533870.0*z10[i]*r[i] + 1562275.0*z8[i]*r2[i] - 749892.0*z6[i]*r3[i] + 171171.0*z4[i]*r4[i] - 16302.0*z2[i]*r5[i] + 429.0*r6[i]);
                 preCoef[t235+i] = 0.00105013854311783*y[i]*(9694845.0*z14[i] - 30421755.0*z12[i]*r[i] + 37182145.0*z10[i]*r2[i] - 22309287.0*z8[i]*r3[i] + 6789783.0*z6[i]*r4[i] - 969969.0*z4[i]*r5[i] + 51051.0*z2[i]*r6[i] - 429.0*r7[i]);
                 preCoef[t236+i] = 0.000766912758094997*z[i]*(9694845.0*z14[i] - 35102025.0*z12[i]*r[i] + 50702925.0*z10[i]*r2[i] - 37182145.0*z8[i]*r3[i] + 14549535.0*z6[i]*r4[i] - 2909907.0*z4[i]*r5[i] + 255255.0*z2[i]*r6[i] - 6435.0*r7[i]);
                 preCoef[t237+i] = 0.00105013854311783*x[i]*(9694845.0*z14[i] - 30421755.0*z12[i]*r[i] + 37182145.0*z10[i]*r2[i] - 22309287.0*z8[i]*r3[i] + 6789783.0*z6[i]*r4[i] - 969969.0*z4[i]*r5[i] + 51051.0*z2[i]*r6[i] - 429.0*r7[i]);
                 preCoef[t238+i] = 0.016200748390692*z[i]*(x2[i] - y2[i])*(570285.0*z12[i] - 1533870.0*z10[i]*r[i] + 1562275.0*z8[i]*r2[i] - 749892.0*z6[i]*r3[i] + 171171.0*z4[i]*r4[i] - 16302.0*z2[i]*r5[i] + 429.0*r6[i]);
                 preCoef[t239+i] = 0.0413039660894728*x[i]*(x2[i] - 3.0*y2[i])*(190095.0*z12[i] - 432630.0*z10[i]*r[i] + 360525.0*z8[i]*r2[i] - 134596.0*z6[i]*r3[i] + 21945.0*z4[i]*r4[i] - 1254.0*z2[i]*r5[i] + 11.0*r6[i]);
                 preCoef[t240+i] = 0.623676210991114*z[i]*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(10005.0*z10[i] - 18975.0*z8[i]*r[i] + 12650.0*z6[i]*r2[i] - 3542.0*z4[i]*r3[i] + 385.0*z2[i]*r4[i] - 11.0*r5[i]);
                 preCoef[t241+i] = 0.462530657238986*x[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(10005.0*z10[i] - 15525.0*z8[i]*r[i] + 8050.0*z6[i]*r2[i] - 1610.0*z4[i]*r3[i] + 105.0*z2[i]*r4[i] - r5[i]);
                 preCoef[t242+i] = 0.319176476700607*z[i]*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(10005.0*z8[i] - 12420.0*z6[i]*r[i] + 4830.0*z4[i]*r2[i] - 644.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t243+i] = 0.0680486534764293*x[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(30015.0*z8[i] - 28980.0*z6[i]*r[i] + 8050.0*z4[i]*r2[i] - 644.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t244+i] = 0.923056845568976*z[i]*(1305.0*z6[i] - 945.0*z4[i]*r[i] + 175.0*z2[i]*r2[i] - 7.0*r3[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t245+i] = 0.49850767216857*x[i]*(1305.0*z6[i] - 675.0*z4[i]*r[i] + 75.0*z2[i]*r2[i] - r3[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t246+i] = 1.22108942967563*z[i]*(261.0*z4[i] - 90.0*z2[i]*r[i] + 5.0*r2[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t247+i] = 0.53548313829403*x[i]*(261.0*z4[i] - 54.0*z2[i]*r[i] + r2[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t248+i] = 1.85496800424338*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 26.0*z2[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t249+i] = -0.60718080651189*x[i]*(x2[i] + y2[i] - 28.0*z2[i])*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t250+i] = 4.62415125663001*z[i]*(x14[i] - 91.0*x12[i]*y2[i] + 1001.0*x10[i]*y4[i] - 3003.0*x8[i]*y6[i] + 3003.0*x6[i]*y8[i] - 1001.0*x4[i]*y10[i] + 91.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t251+i] = 0.844250650857373*x[i]*(x14[i] - 105.0*x12[i]*y2[i] + 1365.0*x10[i]*y4[i] - 5005.0*x8[i]*y6[i] + 6435.0*x6[i]*y8[i] - 3003.0*x4[i]*y10[i] + 455.0*x2[i]*y12[i] - 15.0*y14[i]);
    }
    if (lMax > 15){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t252+i] = 13.7174494214084*x[i]*y[i]*(x14[i] - 35.0*x12[i]*y2[i] + 273.0*x10[i]*y4[i] - 715.0*x8[i]*y6[i] + 715.0*x6[i]*y8[i] - 273.0*x4[i]*y10[i] + 35.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t253+i] = 4.84985075323068*y[i]*z[i]*(15.0*x14[i] - 455.0*x12[i]*y2[i] + 3003.0*x10[i]*y4[i] - 6435.0*x8[i]*y6[i] + 5005.0*x6[i]*y8[i] - 1365.0*x4[i]*y10[i] + 105.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t254+i] = -1.23186332318453*x[i]*y[i]*(x2[i] + y2[i] - 30.0*z2[i])*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t255+i] = 1.94774693364361*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 28.0*z2[i])*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t256+i] = 0.723375051052533*x[i]*y[i]*(899.0*z4[i] - 174.0*z2[i]*r[i] + 3.0*r2[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t257+i] = 0.427954451513054*y[i]*z[i]*(899.0*z4[i] - 290.0*z2[i]*r[i] + 15.0*r2[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t258+i] = 0.2017396631359*x[i]*y[i]*(8091.0*z6[i] - 3915.0*z4[i]*r[i] + 405.0*z2[i]*r2[i] - 5.0*r3[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t259+i] = 0.194401203675805*y[i]*z[i]*(8091.0*z6[i] - 5481.0*z4[i]*r[i] + 945.0*z2[i]*r2[i] - 35.0*r3[i])*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
                 preCoef[t260+i] = 0.549849637559955*x[i]*y[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(40455.0*z8[i] - 36540.0*z6[i]*r[i] + 9450.0*z4[i]*r2[i] - 700.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t261+i] = 0.336712761819039*y[i]*z[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(13485.0*z8[i] - 15660.0*z6[i]*r[i] + 5670.0*z4[i]*r2[i] - 700.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t262+i] = 0.0444043640573282*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(310155.0*z10[i] - 450225.0*z8[i]*r[i] + 217350.0*z6[i]*r2[i] - 40250.0*z4[i]*r3[i] + 2415.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t263+i] = 0.031398626939213*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(310155.0*z10[i] - 550275.0*z8[i]*r[i] + 341550.0*z6[i]*r2[i] - 88550.0*z4[i]*r3[i] + 8855.0*z2[i]*r4[i] - 231.0*r5[i]);
                 preCoef[t264+i] = 0.166145916780101*x[i]*y[i]*(x2[i] - y2[i])*(310155.0*z12[i] - 660330.0*z10[i]*r[i] + 512325.0*z8[i]*r2[i] - 177100.0*z6[i]*r3[i] + 26565.0*z4[i]*r4[i] - 1386.0*z2[i]*r5[i] + 11.0*r6[i]);
                 preCoef[t265+i] = 0.0515196617272515*y[i]*z[i]*(3.0*x2[i] - y2[i])*(310155.0*z12[i] - 780390.0*z10[i]*r[i] + 740025.0*z8[i]*r2[i] - 328900.0*z6[i]*r3[i] + 69069.0*z4[i]*r4[i] - 6006.0*z2[i]*r5[i] + 143.0*r6[i]);
                 preCoef[t266+i] = 0.00631774627238717*x[i]*y[i]*(5892945.0*z14[i] - 17298645.0*z12[i]*r[i] + 19684665.0*z10[i]*r2[i] - 10935925.0*z8[i]*r3[i] + 3062059.0*z6[i]*r4[i] - 399399.0*z4[i]*r5[i] + 19019.0*z2[i]*r6[i] - 143.0*r7[i]);
                 preCoef[t267+i] = 0.00115345738199354*y[i]*z[i]*(17678835.0*z14[i] - 59879925.0*z12[i]*r[i] + 80528175.0*z10[i]*r2[i] - 54679625.0*z8[i]*r3[i] + 19684665.0*z6[i]*r4[i] - 3594591.0*z4[i]*r5[i] + 285285.0*z2[i]*r6[i] - 6435.0*r7[i]);
                 preCoef[t268+i] = 14862.9380228203*z16[i] - 57533.9536367237*z14[i]*r[i] + 90268.7893265838*z12[i]*r2[i] - 73552.3468586979*z10[i]*r3[i] + 33098.5560864141*z8[i]*r4[i] - 8058.77887321386*z6[i]*r5[i] + 959.378437287364*z4[i]*r6[i] - 43.2802302535653*z2[i]*r7[i] + 0.318236987158568*r8[i];
                 preCoef[t269+i] = 0.00115345738199354*x[i]*z[i]*(17678835.0*z14[i] - 59879925.0*z12[i]*r[i] + 80528175.0*z10[i]*r2[i] - 54679625.0*z8[i]*r3[i] + 19684665.0*z6[i]*r4[i] - 3594591.0*z4[i]*r5[i] + 285285.0*z2[i]*r6[i] - 6435.0*r7[i]);
                 preCoef[t270+i] = 0.00315887313619359*(x2[i] - y2[i])*(5892945.0*z14[i] - 17298645.0*z12[i]*r[i] + 19684665.0*z10[i]*r2[i] - 10935925.0*z8[i]*r3[i] + 3062059.0*z6[i]*r4[i] - 399399.0*z4[i]*r5[i] + 19019.0*z2[i]*r6[i] - 143.0*r7[i]);
                 preCoef[t271+i] = 0.0515196617272515*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(310155.0*z12[i] - 780390.0*z10[i]*r[i] + 740025.0*z8[i]*r2[i] - 328900.0*z6[i]*r3[i] + 69069.0*z4[i]*r4[i] - 6006.0*z2[i]*r5[i] + 143.0*r6[i]);
                 preCoef[t272+i] = 0.0415364791950254*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(310155.0*z12[i] - 660330.0*z10[i]*r[i] + 512325.0*z8[i]*r2[i] - 177100.0*z6[i]*r3[i] + 26565.0*z4[i]*r4[i] - 1386.0*z2[i]*r5[i] + 11.0*r6[i]);
                 preCoef[t273+i] = 0.031398626939213*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(310155.0*z10[i] - 550275.0*z8[i]*r[i] + 341550.0*z6[i]*r2[i] - 88550.0*z4[i]*r3[i] + 8855.0*z2[i]*r4[i] - 231.0*r5[i]);
                 preCoef[t274+i] = 0.0222021820286641*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(310155.0*z10[i] - 450225.0*z8[i]*r[i] + 217350.0*z6[i]*r2[i] - 40250.0*z4[i]*r3[i] + 2415.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t275+i] = 0.336712761819039*x[i]*z[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(13485.0*z8[i] - 15660.0*z6[i]*r[i] + 5670.0*z4[i]*r2[i] - 700.0*z2[i]*r3[i] + 21.0*r4[i]);
                 preCoef[t276+i] = 0.0687312046949944*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i])*(40455.0*z8[i] - 36540.0*z6[i]*r[i] + 9450.0*z4[i]*r2[i] - 700.0*z2[i]*r3[i] + 7.0*r4[i]);
                 preCoef[t277+i] = 0.194401203675805*x[i]*z[i]*(8091.0*z6[i] - 5481.0*z4[i]*r[i] + 945.0*z2[i]*r2[i] - 35.0*r3[i])*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
                 preCoef[t278+i] = 0.10086983156795*(8091.0*z6[i] - 3915.0*z4[i]*r[i] + 405.0*z2[i]*r2[i] - 5.0*r3[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t279+i] = 0.427954451513054*x[i]*z[i]*(899.0*z4[i] - 290.0*z2[i]*r[i] + 15.0*r2[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t280+i] = 0.180843762763133*(899.0*z4[i] - 174.0*z2[i]*r[i] + 3.0*r2[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t281+i] = 1.94774693364361*x[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 28.0*z2[i])*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t282+i] = -0.615931661592266*(x2[i] + y2[i] - 30.0*z2[i])*(x14[i] - 91.0*x12[i]*y2[i] + 1001.0*x10[i]*y4[i] - 3003.0*x8[i]*y6[i] + 3003.0*x6[i]*y8[i] - 1001.0*x4[i]*y10[i] + 91.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t283+i] = 4.84985075323068*x[i]*z[i]*(x14[i] - 105.0*x12[i]*y2[i] + 1365.0*x10[i]*y4[i] - 5005.0*x8[i]*y6[i] + 6435.0*x6[i]*y8[i] - 3003.0*x4[i]*y10[i] + 455.0*x2[i]*y12[i] - 15.0*y14[i]);
                 preCoef[t284+i] = 0.857340588838025*x16[i] - 102.880870660563*x14[i]*y2[i] + 1560.35987168521*x12[i]*y4[i] - 6865.5834354149*x10[i]*y6[i] + 11033.9733783454*x8[i]*y8[i] - 6865.5834354149*x6[i]*y10[i] + 1560.35987168521*x4[i]*y12[i] - 102.880870660563*x2[i]*y14[i] + 0.857340588838025*y16[i];
    }
    if (lMax > 16){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t285+i] = 0.869857171920628*y[i]*(17.0*x16[i] - 680.0*x14[i]*y2[i] + 6188.0*x12[i]*y4[i] - 19448.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 12376.0*x6[i]*y10[i] + 2380.0*x4[i]*y12[i] - 136.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t286+i] = 81.1535251976858*x[i]*y[i]*z[i]*(x14[i] - 35.0*x12[i]*y2[i] + 273.0*x10[i]*y4[i] - 715.0*x8[i]*y6[i] + 715.0*x6[i]*y8[i] - 273.0*x4[i]*y10[i] + 35.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t287+i] = -0.624331775925749*y[i]*(x2[i] + y2[i] - 32.0*z2[i])*(15.0*x14[i] - 455.0*x12[i]*y2[i] + 3003.0*x10[i]*y4[i] - 6435.0*x8[i]*y6[i] + 5005.0*x6[i]*y8[i] - 1365.0*x4[i]*y10[i] + 105.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t288+i] = 12.2343542497898*x[i]*y[i]*z[i]*(-x2[i] - y2[i] + 10.0*z2[i])*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t289+i] = 0.549338722534015*y[i]*(341.0*z4[i] - 62.0*z2[i]*r[i] + r2[i])*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t290+i] = 1.79413275488091*x[i]*y[i]*z[i]*(1023.0*z4[i] - 310.0*z2[i]*r[i] + 15.0*r2[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t291+i] = 0.102009639854704*y[i]*(9889.0*z6[i] - 4495.0*z4[i]*r[i] + 435.0*z2[i]*r2[i] - 5.0*r3[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t292+i] = 0.408038559418814*x[i]*y[i]*z[i]*(9889.0*z6[i] - 6293.0*z4[i]*r[i] + 1015.0*z2[i]*r2[i] - 35.0*r3[i])*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
                 preCoef[t293+i] = 0.013881753693839*y[i]*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i])*(267003.0*z8[i] - 226548.0*z6[i]*r[i] + 54810.0*z4[i]*r2[i] - 3780.0*z2[i]*r3[i] + 35.0*r4[i]);
                 preCoef[t294+i] = 1.69879999122657*x[i]*y[i]*z[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(29667.0*z8[i] - 32364.0*z6[i]*r[i] + 10962.0*z4[i]*r2[i] - 1260.0*z2[i]*r3[i] + 35.0*r4[i]);
                 preCoef[t295+i] = 0.0671509657668754*y[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(148335.0*z10[i] - 202275.0*z8[i]*r[i] + 91350.0*z6[i]*r2[i] - 15750.0*z4[i]*r3[i] + 875.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t296+i] = 0.727382699731321*x[i]*y[i]*z[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(40455.0*z10[i] - 67425.0*z8[i]*r[i] + 39150.0*z6[i]*r2[i] - 9450.0*z4[i]*r3[i] + 875.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t297+i] = 0.0656749401202387*y[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(310155.0*z12[i] - 620310.0*z10[i]*r[i] + 450225.0*z8[i]*r2[i] - 144900.0*z6[i]*r3[i] + 20125.0*z4[i]*r4[i] - 966.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t298+i] = 0.0310675249591503*x[i]*y[i]*z[i]*(x2[i] - y2[i])*(3411705.0*z12[i] - 8064030.0*z10[i]*r[i] + 7153575.0*z8[i]*r2[i] - 2960100.0*z6[i]*r3[i] + 575575.0*z4[i]*r4[i] - 46046.0*z2[i]*r5[i] + 1001.0*r6[i]);
                 preCoef[t299+i] = 0.00317081598837913*y[i]*(3.0*x2[i] - y2[i])*(10235115.0*z14[i] - 28224105.0*z12[i]*r[i] + 30045015.0*z10[i]*r2[i] - 15540525.0*z8[i]*r3[i] + 4029025.0*z6[i]*r4[i] - 483483.0*z4[i]*r5[i] + 21021.0*z2[i]*r6[i] - 143.0*r7[i]);
                 preCoef[t300+i] = 0.0219680575732975*x[i]*y[i]*z[i]*(3411705.0*z14[i] - 10855425.0*z12[i]*r[i] + 13656825.0*z10[i]*r2[i] - 8633625.0*z8[i]*r3[i] + 2877875.0*z6[i]*r4[i] - 483483.0*z4[i]*r5[i] + 35035.0*z2[i]*r6[i] - 715.0*r7[i]);
                 preCoef[t301+i] = 0.0006299772562361*y[i]*(64822395.0*z16[i] - 235717800.0*z14[i]*r[i] + 345972900.0*z12[i]*r2[i] - 262462200.0*z10[i]*r3[i] + 109359250.0*z8[i]*r4[i] - 24496472.0*z6[i]*r5[i] + 2662660.0*z4[i]*r6[i] - 108680.0*z2[i]*r7[i] + 715.0*r8[i]);
                 preCoef[t302+i] = 5.09306425332988e-5*z[i]*(583401555.0*z16[i] - 2404321560.0*z14[i]*r[i] + 4071834900.0*z12[i]*r2[i] - 3650610600.0*z10[i]*r3[i] + 1859107250.0*z8[i]*r4[i] - 535422888.0*z6[i]*r5[i] + 81477396.0*z4[i]*r6[i] - 5542680.0*z2[i]*r7[i] + 109395.0*r8[i]);
                 preCoef[t303+i] = 0.0006299772562361*x[i]*(64822395.0*z16[i] - 235717800.0*z14[i]*r[i] + 345972900.0*z12[i]*r2[i] - 262462200.0*z10[i]*r3[i] + 109359250.0*z8[i]*r4[i] - 24496472.0*z6[i]*r5[i] + 2662660.0*z4[i]*r6[i] - 108680.0*z2[i]*r7[i] + 715.0*r8[i]);
                 preCoef[t304+i] = 0.0109840287866487*z[i]*(x2[i] - y2[i])*(3411705.0*z14[i] - 10855425.0*z12[i]*r[i] + 13656825.0*z10[i]*r2[i] - 8633625.0*z8[i]*r3[i] + 2877875.0*z6[i]*r4[i] - 483483.0*z4[i]*r5[i] + 35035.0*z2[i]*r6[i] - 715.0*r7[i]);
                 preCoef[t305+i] = 0.00317081598837913*x[i]*(x2[i] - 3.0*y2[i])*(10235115.0*z14[i] - 28224105.0*z12[i]*r[i] + 30045015.0*z10[i]*r2[i] - 15540525.0*z8[i]*r3[i] + 4029025.0*z6[i]*r4[i] - 483483.0*z4[i]*r5[i] + 21021.0*z2[i]*r6[i] - 143.0*r7[i]);
                 preCoef[t306+i] = 0.00776688123978757*z[i]*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(3411705.0*z12[i] - 8064030.0*z10[i]*r[i] + 7153575.0*z8[i]*r2[i] - 2960100.0*z6[i]*r3[i] + 575575.0*z4[i]*r4[i] - 46046.0*z2[i]*r5[i] + 1001.0*r6[i]);
                 preCoef[t307+i] = 0.0656749401202387*x[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(310155.0*z12[i] - 620310.0*z10[i]*r[i] + 450225.0*z8[i]*r2[i] - 144900.0*z6[i]*r3[i] + 20125.0*z4[i]*r4[i] - 966.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t308+i] = 0.36369134986566*z[i]*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(40455.0*z10[i] - 67425.0*z8[i]*r[i] + 39150.0*z6[i]*r2[i] - 9450.0*z4[i]*r3[i] + 875.0*z2[i]*r4[i] - 21.0*r5[i]);
                 preCoef[t309+i] = 0.0671509657668754*x[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(148335.0*z10[i] - 202275.0*z8[i]*r[i] + 91350.0*z6[i]*r2[i] - 15750.0*z4[i]*r3[i] + 875.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t310+i] = 0.212349998903322*z[i]*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i])*(29667.0*z8[i] - 32364.0*z6[i]*r[i] + 10962.0*z4[i]*r2[i] - 1260.0*z2[i]*r3[i] + 35.0*r4[i]);
                 preCoef[t311+i] = 0.013881753693839*x[i]*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i])*(267003.0*z8[i] - 226548.0*z6[i]*r[i] + 54810.0*z4[i]*r2[i] - 3780.0*z2[i]*r3[i] + 35.0*r4[i]);
                 preCoef[t312+i] = 0.204019279709407*z[i]*(9889.0*z6[i] - 6293.0*z4[i]*r[i] + 1015.0*z2[i]*r2[i] - 35.0*r3[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t313+i] = 0.102009639854704*x[i]*(9889.0*z6[i] - 4495.0*z4[i]*r[i] + 435.0*z2[i]*r2[i] - 5.0*r3[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t314+i] = 0.448533188720228*z[i]*(1023.0*z4[i] - 310.0*z2[i]*r[i] + 15.0*r2[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t315+i] = 0.549338722534015*x[i]*(341.0*z4[i] - 62.0*z2[i]*r[i] + r2[i])*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t316+i] = 6.11717712489491*z[i]*(-x2[i] - y2[i] + 10.0*z2[i])*(x14[i] - 91.0*x12[i]*y2[i] + 1001.0*x10[i]*y4[i] - 3003.0*x8[i]*y6[i] + 3003.0*x6[i]*y8[i] - 1001.0*x4[i]*y10[i] + 91.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t317+i] = -0.624331775925749*x[i]*(x2[i] + y2[i] - 32.0*z2[i])*(x14[i] - 105.0*x12[i]*y2[i] + 1365.0*x10[i]*y4[i] - 5005.0*x8[i]*y6[i] + 6435.0*x6[i]*y8[i] - 3003.0*x4[i]*y10[i] + 455.0*x2[i]*y12[i] - 15.0*y14[i]);
                 preCoef[t318+i] = 5.07209532485536*z[i]*(x16[i] - 120.0*x14[i]*y2[i] + 1820.0*x12[i]*y4[i] - 8008.0*x10[i]*y6[i] + 12870.0*x8[i]*y8[i] - 8008.0*x6[i]*y10[i] + 1820.0*x4[i]*y12[i] - 120.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t319+i] = 0.869857171920628*x[i]*(x16[i] - 136.0*x14[i]*y2[i] + 2380.0*x12[i]*y4[i] - 12376.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 19448.0*x6[i]*y10[i] + 6188.0*x4[i]*y12[i] - 680.0*x2[i]*y14[i] + 17.0*y16[i]);
    }
    if (lMax > 17){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t320+i] = 1.76371153735666*x[i]*y[i]*(9.0*x16[i] - 408.0*x14[i]*y2[i] + 4284.0*x12[i]*y4[i] - 15912.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 15912.0*x6[i]*y10[i] + 4284.0*x4[i]*y12[i] - 408.0*x2[i]*y14[i] + 9.0*y16[i]);
                 preCoef[t321+i] = 5.29113461206997*y[i]*z[i]*(17.0*x16[i] - 680.0*x14[i]*y2[i] + 6188.0*x12[i]*y4[i] - 19448.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 12376.0*x6[i]*y10[i] + 2380.0*x4[i]*y12[i] - 136.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t322+i] = -10.1185847426968*x[i]*y[i]*(x2[i] + y2[i] - 34.0*z2[i])*(x14[i] - 35.0*x12[i]*y2[i] + 273.0*x10[i]*y4[i] - 715.0*x8[i]*y6[i] + 715.0*x6[i]*y8[i] - 273.0*x4[i]*y10[i] + 35.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t323+i] = 2.12901451204377*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 32.0*z2[i])*(15.0*x14[i] - 455.0*x12[i]*y2[i] + 3003.0*x10[i]*y4[i] - 6435.0*x8[i]*y6[i] + 5005.0*x6[i]*y8[i] - 1365.0*x4[i]*y10[i] + 105.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t324+i] = 1.11184156725673*x[i]*y[i]*(385.0*z4[i] - 66.0*z2[i]*r[i] + r2[i])*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t325+i] = 7.03190349956512*y[i]*z[i]*(77.0*z4[i] - 22.0*z2[i]*r[i] + r2[i])*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t326+i] = 2.06241672263956*x[i]*y[i]*(2387.0*z6[i] - 1023.0*z4[i]*r[i] + 93.0*z2[i]*r2[i] - r3[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t327+i] = 1.49436288677056*y[i]*z[i]*(1705.0*z6[i] - 1023.0*z4[i]*r[i] + 155.0*z2[i]*r2[i] - 5.0*r3[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t328+i] = 0.1962194600598*x[i]*y[i]*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i])*(49445.0*z8[i] - 39556.0*z6[i]*r[i] + 8990.0*z4[i]*r2[i] - 580.0*z2[i]*r3[i] + 5.0*r4[i]);
                 preCoef[t329+i] = 0.0247213282718858*y[i]*z[i]*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i])*(346115.0*z8[i] - 356004.0*z6[i]*r[i] + 113274.0*z4[i]*r2[i] - 12180.0*z2[i]*r3[i] + 315.0*r4[i]);
                 preCoef[t330+i] = 0.541617165840082*x[i]*y[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(207669.0*z10[i] - 267003.0*z8[i]*r[i] + 113274.0*z6[i]*r2[i] - 18270.0*z4[i]*r3[i] + 945.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t331+i] = 1.14494717494913*y[i]*z[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(18879.0*z10[i] - 29667.0*z8[i]*r[i] + 16182.0*z6[i]*r2[i] - 3654.0*z4[i]*r3[i] + 315.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t332+i] = 0.132207111932957*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(471975.0*z12[i] - 890010.0*z10[i]*r[i] + 606825.0*z8[i]*r2[i] - 182700.0*z6[i]*r3[i] + 23625.0*z4[i]*r4[i] - 1050.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t333+i] = 0.0898170459553624*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(471975.0*z12[i] - 1051830.0*z10[i]*r[i] + 876525.0*z8[i]*r2[i] - 339300.0*z6[i]*r3[i] + 61425.0*z4[i]*r4[i] - 4550.0*z2[i]*r5[i] + 91.0*r6[i]);
                 preCoef[t334+i] = 0.140148631920647*x[i]*y[i]*(x2[i] - y2[i])*(1550775.0*z14[i] - 4032015.0*z12[i]*r[i] + 4032015.0*z10[i]*r2[i] - 1950975.0*z8[i]*r3[i] + 470925.0*z6[i]*r4[i] - 52325.0*z4[i]*r5[i] + 2093.0*z2[i]*r6[i] - 13.0*r7[i]);
                 preCoef[t335+i] = 0.0192873206845831*y[i]*z[i]*(3.0*x2[i] - y2[i])*(3411705.0*z14[i] - 10235115.0*z12[i]*r[i] + 12096045.0*z10[i]*r2[i] - 7153575.0*z8[i]*r3[i] + 2220075.0*z6[i]*r4[i] - 345345.0*z4[i]*r5[i] + 23023.0*z2[i]*r6[i] - 429.0*r7[i]);
                 preCoef[t336+i] = 0.0063132576421421*x[i]*y[i]*(23881935.0*z16[i] - 81880920.0*z14[i]*r[i] + 112896420.0*z12[i]*r2[i] - 80120040.0*z10[i]*r3[i] + 31081050.0*z8[i]*r4[i] - 6446440.0*z6[i]*r5[i] + 644644.0*z4[i]*r6[i] - 24024.0*z2[i]*r7[i] + 143.0*r8[i]);
                 preCoef[t337+i] = 0.000684768935318508*y[i]*z[i]*(119409675.0*z16[i] - 463991880.0*z14[i]*r[i] + 738168900.0*z12[i]*r2[i] - 619109400.0*z10[i]*r3[i] + 293543250.0*z8[i]*r4[i] - 78278200.0*z6[i]*r5[i] + 10958948.0*z4[i]*r6[i] - 680680.0*z2[i]*r7[i] + 12155.0*r8[i]);
                 preCoef[t338+i] = 59403.1009679377*z18[i] - 259676.412802699*z16[i]*r[i] + 472138.932368544*z14[i]*r2[i] - 461985.406941263*z12[i]*r3[i] + 262853.766018305*z10[i]*r4[i] - 87617.9220061016*z8[i]*r5[i] + 16355.345441139*z6[i]*r6[i] - 1523.7899479322*z4[i]*r7[i] + 54.4210695690072*z2[i]*r8[i] - 0.318251868824604*r9[i];
                 preCoef[t339+i] = 0.000684768935318508*x[i]*z[i]*(119409675.0*z16[i] - 463991880.0*z14[i]*r[i] + 738168900.0*z12[i]*r2[i] - 619109400.0*z10[i]*r3[i] + 293543250.0*z8[i]*r4[i] - 78278200.0*z6[i]*r5[i] + 10958948.0*z4[i]*r6[i] - 680680.0*z2[i]*r7[i] + 12155.0*r8[i]);
                 preCoef[t340+i] = 0.00315662882107105*(x2[i] - y2[i])*(23881935.0*z16[i] - 81880920.0*z14[i]*r[i] + 112896420.0*z12[i]*r2[i] - 80120040.0*z10[i]*r3[i] + 31081050.0*z8[i]*r4[i] - 6446440.0*z6[i]*r5[i] + 644644.0*z4[i]*r6[i] - 24024.0*z2[i]*r7[i] + 143.0*r8[i]);
                 preCoef[t341+i] = 0.0192873206845831*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(3411705.0*z14[i] - 10235115.0*z12[i]*r[i] + 12096045.0*z10[i]*r2[i] - 7153575.0*z8[i]*r3[i] + 2220075.0*z6[i]*r4[i] - 345345.0*z4[i]*r5[i] + 23023.0*z2[i]*r6[i] - 429.0*r7[i]);
                 preCoef[t342+i] = 0.0350371579801619*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(1550775.0*z14[i] - 4032015.0*z12[i]*r[i] + 4032015.0*z10[i]*r2[i] - 1950975.0*z8[i]*r3[i] + 470925.0*z6[i]*r4[i] - 52325.0*z4[i]*r5[i] + 2093.0*z2[i]*r6[i] - 13.0*r7[i]);
                 preCoef[t343+i] = 0.0898170459553624*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(471975.0*z12[i] - 1051830.0*z10[i]*r[i] + 876525.0*z8[i]*r2[i] - 339300.0*z6[i]*r3[i] + 61425.0*z4[i]*r4[i] - 4550.0*z2[i]*r5[i] + 91.0*r6[i]);
                 preCoef[t344+i] = 0.0661035559664783*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(471975.0*z12[i] - 890010.0*z10[i]*r[i] + 606825.0*z8[i]*r2[i] - 182700.0*z6[i]*r3[i] + 23625.0*z4[i]*r4[i] - 1050.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t345+i] = 1.14494717494913*x[i]*z[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(18879.0*z10[i] - 29667.0*z8[i]*r[i] + 16182.0*z6[i]*r2[i] - 3654.0*z4[i]*r3[i] + 315.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t346+i] = 0.0677021457300103*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i])*(207669.0*z10[i] - 267003.0*z8[i]*r[i] + 113274.0*z6[i]*r2[i] - 18270.0*z4[i]*r3[i] + 945.0*z2[i]*r4[i] - 7.0*r5[i]);
                 preCoef[t347+i] = 0.0247213282718858*x[i]*z[i]*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i])*(346115.0*z8[i] - 356004.0*z6[i]*r[i] + 113274.0*z4[i]*r2[i] - 12180.0*z2[i]*r3[i] + 315.0*r4[i]);
                 preCoef[t348+i] = 0.0981097300298999*(49445.0*z8[i] - 39556.0*z6[i]*r[i] + 8990.0*z4[i]*r2[i] - 580.0*z2[i]*r3[i] + 5.0*r4[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t349+i] = 1.49436288677056*x[i]*z[i]*(1705.0*z6[i] - 1023.0*z4[i]*r[i] + 155.0*z2[i]*r2[i] - 5.0*r3[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t350+i] = 0.515604180659891*(2387.0*z6[i] - 1023.0*z4[i]*r[i] + 93.0*z2[i]*r2[i] - r3[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t351+i] = 7.03190349956512*x[i]*z[i]*(77.0*z4[i] - 22.0*z2[i]*r[i] + r2[i])*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t352+i] = 0.555920783628365*(385.0*z4[i] - 66.0*z2[i]*r[i] + r2[i])*(x14[i] - 91.0*x12[i]*y2[i] + 1001.0*x10[i]*y4[i] - 3003.0*x8[i]*y6[i] + 3003.0*x6[i]*y8[i] - 1001.0*x4[i]*y10[i] + 91.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t353+i] = 2.12901451204377*x[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 32.0*z2[i])*(x14[i] - 105.0*x12[i]*y2[i] + 1365.0*x10[i]*y4[i] - 5005.0*x8[i]*y6[i] + 6435.0*x6[i]*y8[i] - 3003.0*x4[i]*y10[i] + 455.0*x2[i]*y12[i] - 15.0*y14[i]);
                 preCoef[t354+i] = -0.632411546418547*(x2[i] + y2[i] - 34.0*z2[i])*(x16[i] - 120.0*x14[i]*y2[i] + 1820.0*x12[i]*y4[i] - 8008.0*x10[i]*y6[i] + 12870.0*x8[i]*y8[i] - 8008.0*x6[i]*y10[i] + 1820.0*x4[i]*y12[i] - 120.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t355+i] = 5.29113461206997*x[i]*z[i]*(x16[i] - 136.0*x14[i]*y2[i] + 2380.0*x12[i]*y4[i] - 12376.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 19448.0*x6[i]*y10[i] + 6188.0*x4[i]*y12[i] - 680.0*x2[i]*y14[i] + 17.0*y16[i]);
                 preCoef[t356+i] = 0.881855768678329*x18[i] - 134.923932607784*x16[i]*y2[i] + 2698.47865215569*x14[i]*y4[i] - 16370.7704897445*x12[i]*y6[i] + 38588.2447258263*x10[i]*y8[i] - 38588.2447258263*x8[i]*y10[i] + 16370.7704897445*x6[i]*y12[i] - 2698.47865215569*x4[i]*y14[i] + 134.923932607784*x2[i]*y16[i] - 0.881855768678329*y18[i];
    }
    if (lMax > 18){ // OBS!!! SOAPGTO L > 9, from here on, switches to tesseral from normal spherical harmonics
                 preCoef[t357+i] = 0.893383784349949*y[i]*(19.0*x18[i] - 969.0*x16[i]*y2[i] + 11628.0*x14[i]*y4[i] - 50388.0*x12[i]*y6[i] + 92378.0*x10[i]*y8[i] - 75582.0*x8[i]*y10[i] + 27132.0*x6[i]*y12[i] - 3876.0*x4[i]*y14[i] + 171.0*x2[i]*y16[i] - y18[i]);
                 preCoef[t358+i] = 11.0143750205445*x[i]*y[i]*z[i]*(9.0*x16[i] - 408.0*x14[i]*y2[i] + 4284.0*x12[i]*y4[i] - 15912.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 15912.0*x6[i]*y10[i] + 4284.0*x4[i]*y12[i] - 408.0*x2[i]*y14[i] + 9.0*y16[i]);
                 preCoef[t359+i] = -0.640197544188601*y[i]*(x2[i] + y2[i] - 36.0*z2[i])*(17.0*x16[i] - 680.0*x14[i]*y2[i] + 6188.0*x12[i]*y4[i] - 19448.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 12376.0*x6[i]*y10[i] + 2380.0*x4[i]*y12[i] - 136.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t360+i] = 35.4833495492953*x[i]*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 34.0*z2[i])*(x14[i] - 35.0*x12[i]*y2[i] + 273.0*x10[i]*y4[i] - 715.0*x8[i]*y6[i] + 715.0*x6[i]*y8[i] - 273.0*x4[i]*y10[i] + 35.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t361+i] = 0.187430649022538*y[i]*(1295.0*z4[i] - 210.0*z2[i]*r[i] + 3.0*r2[i])*(15.0*x14[i] - 455.0*x12[i]*y2[i] + 3003.0*x10[i]*y4[i] - 6435.0*x8[i]*y6[i] + 5005.0*x6[i]*y8[i] - 1365.0*x4[i]*y10[i] + 105.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t362+i] = 4.8875933516657*x[i]*y[i]*z[i]*(259.0*z4[i] - 70.0*z2[i]*r[i] + 3.0*r2[i])*(7.0*x12[i] - 182.0*x10[i]*y2[i] + 1001.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1001.0*x4[i]*y8[i] - 182.0*x2[i]*y10[i] + 7.0*y12[i]);
                 preCoef[t363+i] = 0.521019201915023*y[i]*(2849.0*z6[i] - 1155.0*z4[i]*r[i] + 99.0*z2[i]*r2[i] - r3[i])*(13.0*x12[i] - 286.0*x10[i]*y2[i] + 1287.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 715.0*x4[i]*y8[i] - 78.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t364+i] = 31.1916055279426*x[i]*y[i]*z[i]*(407.0*z6[i] - 231.0*z4[i]*r[i] + 33.0*z2[i]*r2[i] - r3[i])*(3.0*x10[i] - 55.0*x8[i]*y2[i] + 198.0*x6[i]*y4[i] - 198.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - 3.0*y10[i]);
                 preCoef[t365+i] = 0.495167232923569*y[i]*(12617.0*z8[i] - 9548.0*z6[i]*r[i] + 2046.0*z4[i]*r2[i] - 124.0*z2[i]*r3[i] + r4[i])*(11.0*x10[i] - 165.0*x8[i]*y2[i] + 462.0*x6[i]*y4[i] - 330.0*x4[i]*y6[i] + 55.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t366+i] = 0.361619017612871*x[i]*y[i]*z[i]*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i])*(63085.0*z8[i] - 61380.0*z6[i]*r[i] + 18414.0*z4[i]*r2[i] - 1860.0*z2[i]*r3[i] + 45.0*r4[i]);
                 preCoef[t367+i] = 0.0530874997258304*y[i]*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i])*(365893.0*z10[i] - 445005.0*z8[i]*r[i] + 178002.0*z6[i]*r2[i] - 26970.0*z4[i]*r3[i] + 1305.0*z2[i]*r4[i] - 9.0*r5[i]);
                 preCoef[t368+i] = 1.06477924459395*x[i]*y[i]*z[i]*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i])*(232841.0*z10[i] - 346115.0*z8[i]*r[i] + 178002.0*z6[i]*r2[i] - 37758.0*z4[i]*r3[i] + 3045.0*z2[i]*r4[i] - 63.0*r5[i]);
                 preCoef[t369+i] = 0.0665487027871216*y[i]*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i])*(698523.0*z12[i] - 1246014.0*z10[i]*r[i] + 801009.0*z8[i]*r2[i] - 226548.0*z6[i]*r3[i] + 27405.0*z4[i]*r4[i] - 1134.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t370+i] = 0.188228156079767*x[i]*y[i]*z[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(698523.0*z12[i] - 1472562.0*z10[i]*r[i] + 1157013.0*z8[i]*r2[i] - 420732.0*z6[i]*r3[i] + 71253.0*z4[i]*r4[i] - 4914.0*z2[i]*r5[i] + 91.0*r6[i]);
                 preCoef[t371+i] = 0.0352142635297349*y[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(2494725.0*z14[i] - 6135675.0*z12[i]*r[i] + 5785065.0*z10[i]*r2[i] - 2629575.0*z8[i]*r3[i] + 593775.0*z6[i]*r4[i] - 61425.0*z4[i]*r5[i] + 2275.0*z2[i]*r6[i] - 13.0*r7[i]);
                 preCoef[t372+i] = 0.890858231034903*x[i]*y[i]*z[i]*(x2[i] - y2[i])*(498945.0*z14[i] - 1415925.0*z12[i]*r[i] + 1577745.0*z10[i]*r2[i] - 876525.0*z8[i]*r3[i] + 254475.0*z6[i]*r4[i] - 36855.0*z4[i]*r5[i] + 2275.0*z2[i]*r6[i] - 39.0*r7[i]);
                 preCoef[t373+i] = 0.0116097988782603*y[i]*(3.0*x2[i] - y2[i])*(11475735.0*z16[i] - 37218600.0*z14[i]*r[i] + 48384180.0*z12[i]*r2[i] - 32256120.0*z10[i]*r3[i] + 11705850.0*z8[i]*r4[i] - 2260440.0*z6[i]*r5[i] + 209300.0*z4[i]*r6[i] - 7176.0*z2[i]*r7[i] + 39.0*r8[i]);
                 preCoef[t374+i] = 0.00240131363330655*x[i]*y[i]*z[i]*(126233085.0*z16[i] - 463991880.0*z14[i]*r[i] + 695987820.0*z12[i]*r2[i] - 548354040.0*z10[i]*r3[i] + 243221550.0*z8[i]*r4[i] - 60386040.0*z6[i]*r5[i] + 7827820.0*z4[i]*r6[i] - 447304.0*z2[i]*r7[i] + 7293.0*r8[i]);
                 preCoef[t375+i] = 0.000185265368964422*y[i]*(883631595.0*z18[i] - 3653936055.0*z16[i]*r[i] + 6263890380.0*z14[i]*r2[i] - 5757717420.0*z12[i]*r3[i] + 3064591530.0*z10[i]*r4[i] - 951080130.0*z8[i]*r5[i] + 164384220.0*z6[i]*r6[i] - 14090076.0*z4[i]*r7[i] + 459459.0*z2[i]*r8[i] - 2431.0*r9[i]);
                 preCoef[t376+i] = 2.68811250303113e-5*z[i]*(4418157975.0*z18[i] - 20419054425.0*z16[i]*r[i] + 39671305740.0*z14[i]*r2[i] - 42075627300.0*z12[i]*r3[i] + 26466926850.0*z10[i]*r4[i] - 10039179150.0*z8[i]*r5[i] + 2230928700.0*z6[i]*r6[i] - 267711444.0*z4[i]*r7[i] + 14549535.0*z2[i]*r8[i] - 230945.0*r9[i]);
                 preCoef[t377+i] = 0.000185265368964422*x[i]*(883631595.0*z18[i] - 3653936055.0*z16[i]*r[i] + 6263890380.0*z14[i]*r2[i] - 5757717420.0*z12[i]*r3[i] + 3064591530.0*z10[i]*r4[i] - 951080130.0*z8[i]*r5[i] + 164384220.0*z6[i]*r6[i] - 14090076.0*z4[i]*r7[i] + 459459.0*z2[i]*r8[i] - 2431.0*r9[i]);
                 preCoef[t378+i] = 0.00120065681665328*z[i]*(x2[i] - y2[i])*(126233085.0*z16[i] - 463991880.0*z14[i]*r[i] + 695987820.0*z12[i]*r2[i] - 548354040.0*z10[i]*r3[i] + 243221550.0*z8[i]*r4[i] - 60386040.0*z6[i]*r5[i] + 7827820.0*z4[i]*r6[i] - 447304.0*z2[i]*r7[i] + 7293.0*r8[i]);
                 preCoef[t379+i] = 0.0116097988782603*x[i]*(x2[i] - 3.0*y2[i])*(11475735.0*z16[i] - 37218600.0*z14[i]*r[i] + 48384180.0*z12[i]*r2[i] - 32256120.0*z10[i]*r3[i] + 11705850.0*z8[i]*r4[i] - 2260440.0*z6[i]*r5[i] + 209300.0*z4[i]*r6[i] - 7176.0*z2[i]*r7[i] + 39.0*r8[i]);
                 preCoef[t380+i] = 0.222714557758726*z[i]*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(498945.0*z14[i] - 1415925.0*z12[i]*r[i] + 1577745.0*z10[i]*r2[i] - 876525.0*z8[i]*r3[i] + 254475.0*z6[i]*r4[i] - 36855.0*z4[i]*r5[i] + 2275.0*z2[i]*r6[i] - 39.0*r7[i]);
                 preCoef[t381+i] = 0.0352142635297349*x[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(2494725.0*z14[i] - 6135675.0*z12[i]*r[i] + 5785065.0*z10[i]*r2[i] - 2629575.0*z8[i]*r3[i] + 593775.0*z6[i]*r4[i] - 61425.0*z4[i]*r5[i] + 2275.0*z2[i]*r6[i] - 13.0*r7[i]);
                 preCoef[t382+i] = 0.0941140780398835*z[i]*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i])*(698523.0*z12[i] - 1472562.0*z10[i]*r[i] + 1157013.0*z8[i]*r2[i] - 420732.0*z6[i]*r3[i] + 71253.0*z4[i]*r4[i] - 4914.0*z2[i]*r5[i] + 91.0*r6[i]);
                 preCoef[t383+i] = 0.0665487027871216*x[i]*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i])*(698523.0*z12[i] - 1246014.0*z10[i]*r[i] + 801009.0*z8[i]*r2[i] - 226548.0*z6[i]*r3[i] + 27405.0*z4[i]*r4[i] - 1134.0*z2[i]*r5[i] + 7.0*r6[i]);
                 preCoef[t384+i] = 0.133097405574243*z[i]*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i])*(232841.0*z10[i] - 346115.0*z8[i]*r[i] + 178002.0*z6[i]*r2[i] - 37758.0*z4[i]*r3[i] + 3045.0*z2[i]*r4[i] - 63.0*r5[i]);
                 preCoef[t385+i] = 0.0530874997258304*x[i]*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i])*(365893.0*z10[i] - 445005.0*z8[i]*r[i] + 178002.0*z6[i]*r2[i] - 26970.0*z4[i]*r3[i] + 1305.0*z2[i]*r4[i] - 9.0*r5[i]);
                 preCoef[t386+i] = 0.180809508806436*z[i]*(63085.0*z8[i] - 61380.0*z6[i]*r[i] + 18414.0*z4[i]*r2[i] - 1860.0*z2[i]*r3[i] + 45.0*r4[i])*(x10[i] - 45.0*x8[i]*y2[i] + 210.0*x6[i]*y4[i] - 210.0*x4[i]*y6[i] + 45.0*x2[i]*y8[i] - y10[i]);
                 preCoef[t387+i] = 0.495167232923569*x[i]*(12617.0*z8[i] - 9548.0*z6[i]*r[i] + 2046.0*z4[i]*r2[i] - 124.0*z2[i]*r3[i] + r4[i])*(x10[i] - 55.0*x8[i]*y2[i] + 330.0*x6[i]*y4[i] - 462.0*x4[i]*y6[i] + 165.0*x2[i]*y8[i] - 11.0*y10[i]);
                 preCoef[t388+i] = 7.79790138198564*z[i]*(407.0*z6[i] - 231.0*z4[i]*r[i] + 33.0*z2[i]*r2[i] - r3[i])*(x12[i] - 66.0*x10[i]*y2[i] + 495.0*x8[i]*y4[i] - 924.0*x6[i]*y6[i] + 495.0*x4[i]*y8[i] - 66.0*x2[i]*y10[i] + y12[i]);
                 preCoef[t389+i] = 0.521019201915023*x[i]*(2849.0*z6[i] - 1155.0*z4[i]*r[i] + 99.0*z2[i]*r2[i] - r3[i])*(x12[i] - 78.0*x10[i]*y2[i] + 715.0*x8[i]*y4[i] - 1716.0*x6[i]*y6[i] + 1287.0*x4[i]*y8[i] - 286.0*x2[i]*y10[i] + 13.0*y12[i]);
                 preCoef[t390+i] = 2.44379667583285*z[i]*(259.0*z4[i] - 70.0*z2[i]*r[i] + 3.0*r2[i])*(x14[i] - 91.0*x12[i]*y2[i] + 1001.0*x10[i]*y4[i] - 3003.0*x8[i]*y6[i] + 3003.0*x6[i]*y8[i] - 1001.0*x4[i]*y10[i] + 91.0*x2[i]*y12[i] - y14[i]);
                 preCoef[t391+i] = 0.187430649022538*x[i]*(1295.0*z4[i] - 210.0*z2[i]*r[i] + 3.0*r2[i])*(x14[i] - 105.0*x12[i]*y2[i] + 1365.0*x10[i]*y4[i] - 5005.0*x8[i]*y6[i] + 6435.0*x6[i]*y8[i] - 3003.0*x4[i]*y10[i] + 455.0*x2[i]*y12[i] - 15.0*y14[i]);
                 preCoef[t392+i] = 2.21770934683096*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 34.0*z2[i])*(x16[i] - 120.0*x14[i]*y2[i] + 1820.0*x12[i]*y4[i] - 8008.0*x10[i]*y6[i] + 12870.0*x8[i]*y8[i] - 8008.0*x6[i]*y10[i] + 1820.0*x4[i]*y12[i] - 120.0*x2[i]*y14[i] + y16[i]);
                 preCoef[t393+i] = -0.640197544188601*x[i]*(x2[i] + y2[i] - 36.0*z2[i])*(x16[i] - 136.0*x14[i]*y2[i] + 2380.0*x12[i]*y4[i] - 12376.0*x10[i]*y6[i] + 24310.0*x8[i]*y8[i] - 19448.0*x6[i]*y10[i] + 6188.0*x4[i]*y12[i] - 680.0*x2[i]*y14[i] + 17.0*y16[i]);
                 preCoef[t394+i] = 5.50718751027224*z[i]*(x18[i] - 153.0*x16[i]*y2[i] + 3060.0*x14[i]*y4[i] - 18564.0*x12[i]*y6[i] + 43758.0*x10[i]*y8[i] - 43758.0*x8[i]*y10[i] + 18564.0*x6[i]*y12[i] - 3060.0*x4[i]*y14[i] + 153.0*x2[i]*y16[i] - y18[i]);
                 preCoef[t395+i] = 0.893383784349949*x[i]*(x18[i] - 171.0*x16[i]*y2[i] + 3876.0*x14[i]*y4[i] - 27132.0*x12[i]*y6[i] + 75582.0*x10[i]*y8[i] - 92378.0*x8[i]*y10[i] + 50388.0*x6[i]*y12[i] - 11628.0*x4[i]*y14[i] + 969.0*x2[i]*y16[i] - 19.0*y18[i]);
    }
  }
}
//================================================================
void getC(double* C, double* preCoef, double* x, double* y, double* z,double* r2, double* bOa, double* aOa, double* exes,  int totalAN, int Asize, int Ns, int Ntypes, int lMax, int posI, int typeJ, int Nx2, int Nx3, int Nx4, int Nx5, int Nx6, int Nx7, int Nx8, int Nx9, int Nx10, int Nx11, int Nx12, int Nx13, int Nx14, int Nx15, int Nx16, int Nx17, int Nx18, int Nx19, int Nx20, int Nx21, int Nx22, int Nx23, int Nx24, int Nx25, int Nx26, int Nx27, int Nx28, int Nx29, int Nx30, int Nx31, int Nx32, int Nx33, int Nx34, int Nx35, int Nx36, int Nx37, int Nx38, int Nx39, int Nx40, int Nx41, int Nx42, int Nx43, int Nx44, int Nx45, int Nx46, int Nx47, int Nx48, int Nx49, int Nx50, int Nx51, int Nx52, int Nx53, int Nx54, int Nx55, int Nx56, int Nx57, int Nx58, int Nx59, int Nx60, int Nx61, int Nx62, int Nx63, int Nx64, int Nx65, int Nx66, int Nx67, int Nx68, int Nx69, int Nx70, int Nx71, int Nx72, int Nx73, int Nx74, int Nx75, int Nx76, int Nx77, int Nx78, int Nx79, int Nx80, int Nx81, int Nx82, int Nx83, int Nx84, int Nx85, int Nx86, int Nx87, int Nx88, int Nx89, int Nx90, int Nx91, int Nx92, int Nx93, int Nx94, int Nx95, int Nx96, int Nx97, int Nx98, int Nx99, int t2, int t3, int t4, int t5, int t6, int t7, int t8, int t9, int t10, int t11, int t12, int t13, int t14, int t15, int t16, int t17, int t18, int t19, int t20, int t21, int t22, int t23, int t24, int t25, int t26, int t27, int t28, int t29, int t30, int t31, int t32, int t33, int t34, int t35, int t36, int t37, int t38, int t39, int t40, int t41, int t42, int t43, int t44, int t45, int t46, int t47, int t48, int t49, int t50, int t51, int t52, int t53, int t54, int t55, int t56, int t57, int t58, int t59, int t60, int t61, int t62, int t63, int t64, int t65, int t66, int t67, int t68, int t69, int t70, int t71, int t72, int t73, int t74, int t75, int t76, int t77, int t78, int t79, int t80, int t81, int t82, int t83, int t84, int t85, int t86, int t87, int t88, int t89, int t90, int t91, int t92, int t93, int t94, int t95, int t96, int t97, int t98, int t99){

  if(Asize == 0){return;}
  double sumMe = 0; int NsNs = Ns*Ns;  int NsJ = 100*Ns*typeJ; int LNsNs;
  int LNs; int NsTsI = 100*Ns*Ntypes*posI;
  for(int k = 0; k < Ns; k++){
    sumMe = 0; for(int i = 0; i < Asize; i++){ sumMe += exp(aOa[k]*r2[i]);}
    for(int n = 0; n < Ns; n++){ C[NsTsI + NsJ + n] += bOa[n*Ns + k]*sumMe; }
  } if(lMax > 0) { LNsNs=NsNs; LNs=Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c10*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*z[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c11Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*x[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx2 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c11Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*y[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx3 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 1) { LNsNs=2*NsNs; LNs=2*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c20*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx4 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c21Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN + i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx5 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c21Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[t2+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx6 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c22Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[t3+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx7 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c22Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[t4+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx8 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
  }} if(lMax > 2) { LNsNs=3*NsNs; LNs=3*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c30*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t5+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx9 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c31Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t6+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx10 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c31Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t7+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx11 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c32Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t8+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx12 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c32Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t9+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx13 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c33Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t10+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx14 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c33Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t11+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx15 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
  }} if(lMax > 3) { LNsNs=4*NsNs; LNs=4*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c40*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t12+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx16 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c41Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t13+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx17 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c41Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t14+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx18 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c42Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t15+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx19 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c42Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t16+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx20 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c43Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t17+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx21 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c43Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t18+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx22 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c44Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t19+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx23 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c44Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t20+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx24 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }

  }} if(lMax > 4) { LNsNs=5*NsNs; LNs=5*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c50*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t21+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx25 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c51Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t22+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx26 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c51Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t23+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx27 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c52Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t24+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx28 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c52Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t25+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx29 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c53Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t26+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx30 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c53Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t27+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx31 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c54Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t28+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx32 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c54Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t29+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx33 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c55Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t30+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx34 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c55Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t31+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx35 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 5) { LNsNs=6*NsNs; LNs=6*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c60*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t32+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx36 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c61Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t33+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx37 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c61Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t34+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx38 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c62Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t35+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx39 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c62Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t36+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx40 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c63Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t37+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx41 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c63Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t38+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx42 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c64Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t39+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx43 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c64Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t40+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx44 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c65Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t41+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx45 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c65Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t42+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx46 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c66Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t43+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx47 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c66Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t44+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx48 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 6) { LNsNs=7*NsNs; LNs=7*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c70*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t45+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx49 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c71Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t46+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx50 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c71Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t47+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx51 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c72Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t48+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx52 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c72Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t49+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx53 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c73Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t50+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx54 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c73Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t51+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx55 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c74Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t52+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx56 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c74Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t53+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx57 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c75Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t54+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx58 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c75Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t55+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx59 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c76Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t56+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx60 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c76Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t57+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx61 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c77Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t58+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx62 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c77Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t59+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx63 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 7) { LNsNs=8*NsNs; LNs=8*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c80*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t60+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx64 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c81Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t61+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx65 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c81Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t62+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx66 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c82Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t63+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx67 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c82Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t64+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx68 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c83Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t65+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx69 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c83Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t66+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx70 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c84Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t67+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx71 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c84Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t68+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx72 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c85Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t69+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx73 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c85Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t70+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx74 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c86Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t71+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx75 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c86Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t72+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx76 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c87Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t73+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx77 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c87Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t74+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx78 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c88Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t75+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx79 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c88Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t76+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx80 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 8) { LNsNs=9*NsNs; LNs=9*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c90*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t77+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx81 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c91Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t78+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx82 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c91Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t79+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx83 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c92Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t80+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx84 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c92Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t81+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx85 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c93Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t82+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx86 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c93Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t83+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx87 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c94Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t84+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx88 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c94Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t85+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx89 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c95Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t86+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx90 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c95Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t87+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx91 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c96Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t88+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx92 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c96Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t89+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx93 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c97Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t90+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx94 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c97Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t91+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx95 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c98Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t92+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx96 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c98Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t93+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx97 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c99Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t94+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx98 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c99Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[t95+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Nx99 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }}
}
//=======================================================================
/**
 * Used to calculate the partial power spectrum without crossover.
 */
void getPNoCross(double* soapMat, double* Cnnd, int Ns, int Ts, int Hs, int lMax){
  int NsTs100 = Ns*Ts*100;
  int Ns100 = Ns*100;
  int NsNs = (Ns*(Ns+1))/2;
  int NsNsLmax = NsNs*(lMax+1);
  int NsNsLmaxTs = NsNsLmax*Ts;
  int shiftN = 0;

  //  double   cs0  = pow(PIHalf,2);
  //  double   cs1  = pow(2.7206990464,2);
  //  double cs2  = 2*pow(1.9238247452,2); double   cs3  = pow(1.7562036828,2); double cs4  = 2*pow(4.3018029072,2);
  //  double cs5  = 2*pow(2.1509014536,2); double   cs6  = pow(2.0779682205,2); double cs7  = 2*pow(1.7995732672,2);
  //  double cs8  = 2*pow(5.6907503408,2); double cs9  = 2*pow(2.3232390981,2); double   cs10 = pow(0.5890486225,2);
  //  double cs11 = 2*pow(2.6343055241,2); double cs12 = 2*pow(1.8627352998,2); double cs13 = 2*pow(6.9697172942,2);
  //  double cs14 = 2*pow(2.4641671809,2); double   cs15 = pow(0.6512177548,2); double cs16 = 2*pow(1.7834332706,2);
  //  double cs17 = 2*pow(9.4370418280,2); double cs18 = 2*pow(1.9263280966,2); double cs19 = 2*pow(8.1727179596,2);
  //  double cs20 = 2*pow(2.5844403427,2); double   cs21 = pow(0.3539741687,2); double cs22 = 2*pow(2.2940148014,2);
  //  double cs23 = 2*pow(1.8135779397,2); double cs24 = 2*pow(3.6271558793,2); double cs25 = 2*pow(1.9866750947,2);
  //  double cs26 = 2*pow(9.3183321738,2); double cs27 = 2*pow(2.6899707945,2); double   cs28 = pow(0.3802292509,2);
  //  double cs29 = 2*pow(0.3556718963,2); double cs30 = 2*pow(0.8712146618,2); double cs31 = 2*pow(0.6160417952,2);
  //  double cs32 = 2*pow(4.0863589798,2); double cs33 = 2*pow(2.0431794899,2); double cs34 = 2*pow(10.418212089,2);
  //  double cs35 = 2*pow(2.7843843014,2); double   cs36 = pow(0.0505981185,2); double cs37 = 2*pow(0.4293392727,2);
  //  double cs38 = 2*pow(1.7960550366,2); double cs39 = 2*pow(4.8637400313,2); double cs40 = 2*pow(1.8837184141,2);
  //  double cs41 = 2*pow(13.583686661,2); double cs42 = 2*pow(2.0960083567,2); double cs43 = 2*pow(11.480310577,2);
  //  double cs44 = 2*pow(2.8700776442,2); double   cs45 = pow(0.0534917379,2); double cs46 = 2*pow(0.2537335916,2);
  //  double cs47 = 2*pow(2.3802320735,2); double cs48 = 2*pow(1.8179322747,2); double cs49 = 2*pow(16.055543121,2);
  //  double cs50 = 2*pow(1.9190044477,2); double cs51 = 2*pow(4.9548481782,2); double cs52 = 2*pow(2.1455121971,2);
  //  double cs53 = 2*pow(12.510378411,2); double cs54 = 2*pow(2.9487244699,2);
  double cs0=2.4674011003; double cs1=7.4022033011; double cs2=7.4022033005;
  double cs3=3.0842513755; double cs4=37.0110165048; double cs5=9.2527541262;
  double cs6=4.3179519254; double cs7=6.4769278880; double cs8=64.7692788826;
  double cs9=10.7948798139; double cs10=0.3469782797; double cs11=13.8791311886;
  double cs12=6.9395655942; double cs13=97.1539183221; double cs14=12.1442397908;
  double cs15=0.4240845642; double cs16=6.3612684614; double cs17=178.1155169268;
  double cs18=7.4214798715; double cs19=133.5866376943; double cs20=13.3586637700;
  double cs21=0.1252977121; double cs22=10.5250078181; double cs23=6.5781298867;
  double cs24=26.3125195455; double cs25=7.8937558638; double cs26=173.6626290026;
  double cs27=14.4718857505; double cs28=0.1445742832; double cs29=0.2530049956;
  double cs30=1.5180299739; double cs31=0.7590149869; double cs32=33.3966594236;
  double cs33=8.3491648559; double cs34=217.0782862628; double cs35=15.5055918758;
  double cs36=0.0025601696; double cs37=0.3686644222; double cs38=6.4516273890;
  double cs39=47.3119341841; double cs40=7.0967901272; double cs41=369.0330866085;
  double cs42=8.7865020627; double cs43=263.5950618888; double cs44=16.4746913675;
  double cs45=0.0028613660; double cs46=0.1287614710; double cs47=11.3310094474;
  double cs48=6.6097555108; double cs49=515.5609298206; double cs50=7.3651561406;
  double cs51=49.1010409380; double cs52=9.2064451758; double cs53=313.0191359728;
  double cs54=17.3899519988;

  // The power spectrum is multiplied by an l-dependent prefactor that comes
  // from the normalization of the Wigner D matrices. This prefactor is
  // mentioned in the arrata of the original SOAP paper: On representing
  // chemical environments, Phys. Rev. B 87, 184115 (2013). Here the square
  // root of the prefactor in the dot-product kernel is used, so that after a
  // possible dot-product the full prefactor is recovered.

  // SUM M's UP!
  double prel0 = PI*sqrt(8.0/(1.0));
  for(int i = 0; i < Hs; i++){
    for(int j = 0; j < Ts; j++){
      shiftN = 0;
      for(int k = 0; k < Ns; k++){
        for(int kd = k; kd < Ns; kd++){
          soapMat[NsNsLmaxTs*i+ NsNsLmax*j+ 0 +shiftN] = prel0*(
            cs0*Cnnd[NsTs100*i + Ns100*j + 0 + k]*Cnnd[NsTs100*i + Ns100*j + 0 + kd]);
          shiftN++;
        }
      }
    }
  } if (lMax > 0) {
    double prel1 = PI*sqrt(8.0/(2.0*1.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ NsNs + shiftN] = prel1*(
              cs1*Cnnd[NsTs100*i + Ns100*j + 1*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 1*Ns + kd]
             +cs2*Cnnd[NsTs100*i + Ns100*j + 2*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 2*Ns + kd]
             +cs2*Cnnd[NsTs100*i + Ns100*j + 3*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 3*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 1) {
    double prel2 = PI*sqrt(8.0/(2.0*2.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 2*NsNs + shiftN] = prel2*(
              cs3*Cnnd[NsTs100*i + Ns100*j + 4*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 4*Ns + kd]
             +cs4*Cnnd[NsTs100*i + Ns100*j + 5*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 5*Ns + kd]
             +cs4*Cnnd[NsTs100*i + Ns100*j + 6*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 6*Ns + kd]
             +cs5*Cnnd[NsTs100*i + Ns100*j + 7*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 7*Ns + kd]
             +cs5*Cnnd[NsTs100*i + Ns100*j + 8*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 8*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 2) {
    double prel3 = PI*sqrt(8.0/(2.0*3.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 3*NsNs + shiftN] = prel3*(
              cs6*Cnnd[NsTs100*i + Ns100*j + 9*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 9*Ns + kd]
             +cs7*Cnnd[NsTs100*i + Ns100*j + 10*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 10*Ns + kd]
             +cs7*Cnnd[NsTs100*i + Ns100*j + 11*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 11*Ns + kd]
             +cs8*Cnnd[NsTs100*i + Ns100*j + 12*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 12*Ns + kd]
             +cs8*Cnnd[NsTs100*i + Ns100*j + 13*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 13*Ns + kd]
             +cs9*Cnnd[NsTs100*i + Ns100*j + 14*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 14*Ns + kd]
             +cs9*Cnnd[NsTs100*i + Ns100*j + 15*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 15*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 3) {
    double prel4 = PI*sqrt(8.0/(2.0*4.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 4*NsNs + shiftN] = prel4*(
              cs10*Cnnd[NsTs100*i + Ns100*j + 16*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 16*Ns + kd]
             +cs11*Cnnd[NsTs100*i + Ns100*j + 17*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 17*Ns + kd]
             +cs11*Cnnd[NsTs100*i + Ns100*j + 18*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 18*Ns + kd]
             +cs12*Cnnd[NsTs100*i + Ns100*j + 19*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 19*Ns + kd]
             +cs12*Cnnd[NsTs100*i + Ns100*j + 20*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 20*Ns + kd]
             +cs13*Cnnd[NsTs100*i + Ns100*j + 21*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 21*Ns + kd]
             +cs13*Cnnd[NsTs100*i + Ns100*j + 22*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 22*Ns + kd]
             +cs14*Cnnd[NsTs100*i + Ns100*j + 23*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 23*Ns + kd]
             +cs14*Cnnd[NsTs100*i + Ns100*j + 24*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 24*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if(lMax > 4) {
    double prel5 = PI*sqrt(8.0/(2.0*5.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 5*NsNs + shiftN] = prel5*(
              cs15*Cnnd[NsTs100*i + Ns100*j + 25*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 25*Ns + kd]
             +cs16*Cnnd[NsTs100*i + Ns100*j + 26*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 26*Ns + kd]
             +cs16*Cnnd[NsTs100*i + Ns100*j + 27*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 27*Ns + kd]
             +cs17*Cnnd[NsTs100*i + Ns100*j + 28*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 28*Ns + kd]
             +cs17*Cnnd[NsTs100*i + Ns100*j + 29*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 29*Ns + kd]
             +cs18*Cnnd[NsTs100*i + Ns100*j + 30*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 30*Ns + kd]
             +cs18*Cnnd[NsTs100*i + Ns100*j + 31*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 31*Ns + kd]
             +cs19*Cnnd[NsTs100*i + Ns100*j + 32*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 32*Ns + kd]
             +cs19*Cnnd[NsTs100*i + Ns100*j + 33*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 33*Ns + kd]
             +cs20*Cnnd[NsTs100*i + Ns100*j + 34*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 34*Ns + kd]
             +cs20*Cnnd[NsTs100*i + Ns100*j + 35*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 35*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 5) {
    double prel6 = PI*sqrt(8.0/(2.0*6.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 6*NsNs + shiftN] = prel6*(
              cs21*Cnnd[NsTs100*i + Ns100*j + 36*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 36*Ns + kd]
             +cs22*Cnnd[NsTs100*i + Ns100*j + 37*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 37*Ns + kd]
             +cs22*Cnnd[NsTs100*i + Ns100*j + 38*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 38*Ns + kd]
             +cs23*Cnnd[NsTs100*i + Ns100*j + 39*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 39*Ns + kd]
             +cs23*Cnnd[NsTs100*i + Ns100*j + 40*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 40*Ns + kd]
             +cs24*Cnnd[NsTs100*i + Ns100*j + 41*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 41*Ns + kd]
             +cs24*Cnnd[NsTs100*i + Ns100*j + 42*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 42*Ns + kd]
             +cs25*Cnnd[NsTs100*i + Ns100*j + 43*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 43*Ns + kd]
             +cs25*Cnnd[NsTs100*i + Ns100*j + 44*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 44*Ns + kd]
             +cs26*Cnnd[NsTs100*i + Ns100*j + 45*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 45*Ns + kd]
             +cs26*Cnnd[NsTs100*i + Ns100*j + 46*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 46*Ns + kd]
             +cs27*Cnnd[NsTs100*i + Ns100*j + 47*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 47*Ns + kd]
             +cs27*Cnnd[NsTs100*i + Ns100*j + 48*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 48*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 6) {
    double prel7 = PI*sqrt(8.0/(2.0*7.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 7*NsNs + shiftN] = prel7*(
              cs28*Cnnd[NsTs100*i + Ns100*j + 49*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 49*Ns + kd]
             +cs29*Cnnd[NsTs100*i + Ns100*j + 50*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 50*Ns + kd]
             +cs29*Cnnd[NsTs100*i + Ns100*j + 51*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 51*Ns + kd]
             +cs30*Cnnd[NsTs100*i + Ns100*j + 52*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 52*Ns + kd]
             +cs30*Cnnd[NsTs100*i + Ns100*j + 53*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 53*Ns + kd]
             +cs31*Cnnd[NsTs100*i + Ns100*j + 54*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 54*Ns + kd]
             +cs31*Cnnd[NsTs100*i + Ns100*j + 55*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 55*Ns + kd]
             +cs32*Cnnd[NsTs100*i + Ns100*j + 56*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 56*Ns + kd]
             +cs32*Cnnd[NsTs100*i + Ns100*j + 57*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 57*Ns + kd]
             +cs33*Cnnd[NsTs100*i + Ns100*j + 58*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 58*Ns + kd]
             +cs33*Cnnd[NsTs100*i + Ns100*j + 59*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 59*Ns + kd]
             +cs34*Cnnd[NsTs100*i + Ns100*j + 60*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 60*Ns + kd]
             +cs34*Cnnd[NsTs100*i + Ns100*j + 61*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 61*Ns + kd]
             +cs35*Cnnd[NsTs100*i + Ns100*j + 62*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 62*Ns + kd]
             +cs35*Cnnd[NsTs100*i + Ns100*j + 63*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 63*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 7) {
    double prel8 = PI*sqrt(8.0/(2.0*8.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 8*NsNs + shiftN] = prel8*(
              cs36*Cnnd[NsTs100*i + Ns100*j + 64*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 64*Ns + kd]
             +cs37*Cnnd[NsTs100*i + Ns100*j + 65*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 65*Ns + kd]
             +cs37*Cnnd[NsTs100*i + Ns100*j + 66*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 66*Ns + kd]
             +cs38*Cnnd[NsTs100*i + Ns100*j + 67*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 67*Ns + kd]
             +cs38*Cnnd[NsTs100*i + Ns100*j + 68*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 68*Ns + kd]
             +cs39*Cnnd[NsTs100*i + Ns100*j + 69*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 69*Ns + kd]
             +cs39*Cnnd[NsTs100*i + Ns100*j + 70*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 70*Ns + kd]
             +cs40*Cnnd[NsTs100*i + Ns100*j + 71*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 71*Ns + kd]
             +cs40*Cnnd[NsTs100*i + Ns100*j + 72*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 72*Ns + kd]
             +cs41*Cnnd[NsTs100*i + Ns100*j + 73*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 73*Ns + kd]
             +cs41*Cnnd[NsTs100*i + Ns100*j + 74*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 74*Ns + kd]
             +cs42*Cnnd[NsTs100*i + Ns100*j + 75*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 75*Ns + kd]
             +cs42*Cnnd[NsTs100*i + Ns100*j + 76*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 76*Ns + kd]
             +cs43*Cnnd[NsTs100*i + Ns100*j + 77*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 77*Ns + kd]
             +cs43*Cnnd[NsTs100*i + Ns100*j + 78*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 78*Ns + kd]
             +cs44*Cnnd[NsTs100*i + Ns100*j + 79*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 79*Ns + kd]
             +cs44*Cnnd[NsTs100*i + Ns100*j + 80*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 80*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  } if (lMax > 8) {
    double prel9 = PI*sqrt(8.0/(2.0*9.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 9*NsNs + shiftN] = prel9*(
              cs45*Cnnd[NsTs100*i + Ns100*j + 81*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 81*Ns + kd]
             +cs46*Cnnd[NsTs100*i + Ns100*j + 82*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 82*Ns + kd]
             +cs46*Cnnd[NsTs100*i + Ns100*j + 83*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 83*Ns + kd]
             +cs47*Cnnd[NsTs100*i + Ns100*j + 84*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 84*Ns + kd]
             +cs47*Cnnd[NsTs100*i + Ns100*j + 85*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 85*Ns + kd]
             +cs48*Cnnd[NsTs100*i + Ns100*j + 86*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 86*Ns + kd]
             +cs48*Cnnd[NsTs100*i + Ns100*j + 87*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 87*Ns + kd]
             +cs49*Cnnd[NsTs100*i + Ns100*j + 88*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 88*Ns + kd]
             +cs49*Cnnd[NsTs100*i + Ns100*j + 89*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 89*Ns + kd]
             +cs50*Cnnd[NsTs100*i + Ns100*j + 90*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 90*Ns + kd]
             +cs50*Cnnd[NsTs100*i + Ns100*j + 91*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 91*Ns + kd]
             +cs51*Cnnd[NsTs100*i + Ns100*j + 92*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 92*Ns + kd]
             +cs51*Cnnd[NsTs100*i + Ns100*j + 93*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 93*Ns + kd]
             +cs52*Cnnd[NsTs100*i + Ns100*j + 94*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 94*Ns + kd]
             +cs52*Cnnd[NsTs100*i + Ns100*j + 95*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 95*Ns + kd]
             +cs53*Cnnd[NsTs100*i + Ns100*j + 96*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 96*Ns + kd]
             +cs53*Cnnd[NsTs100*i + Ns100*j + 97*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 97*Ns + kd]
             +cs54*Cnnd[NsTs100*i + Ns100*j + 98*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 98*Ns + kd]
             +cs54*Cnnd[NsTs100*i + Ns100*j + 99*Ns + k]*Cnnd[NsTs100*i + Ns100*j + 99*Ns + kd]);
            shiftN++;
          }
        }
      }
    }
  }
}
//=======================================================================
/**
 * Used to calculate the partial power spectrum.
 */
void getPCrossOver(double* soapMat, double* Cnnd, int Ns, int Ts, int Hs, int lMax){
  int NsTs100 = Ns*Ts*100;
  int Ns100 = Ns*100;
  int NsNs = (Ns*(Ns+1))/2;
  int NsNsLmax = NsNs*(lMax+1) ;
  int NsNsLmaxTs = NsNsLmax*getCrosNum(Ts);
  int shiftN = 0;
  int shiftT = 0;

  //  double   cs0  = pow(PIHalf,2);
  //  double   cs1  = pow(2.7206990464,2);
  //  double cs2  = 2*pow(1.9238247452,2); double   cs3  = pow(1.7562036828,2); double cs4  = 2*pow(4.3018029072,2);
  //  double cs5  = 2*pow(2.1509014536,2); double   cs6  = pow(2.0779682205,2); double cs7  = 2*pow(1.7995732672,2);
  //  double cs8  = 2*pow(5.6907503408,2); double cs9  = 2*pow(2.3232390981,2); double   cs10 = pow(0.5890486225,2);
  //  double cs11 = 2*pow(2.6343055241,2); double cs12 = 2*pow(1.8627352998,2); double cs13 = 2*pow(6.9697172942,2);
  //  double cs14 = 2*pow(2.4641671809,2); double   cs15 = pow(0.6512177548,2); double cs16 = 2*pow(1.7834332706,2);
  //  double cs17 = 2*pow(9.4370418280,2); double cs18 = 2*pow(1.9263280966,2); double cs19 = 2*pow(8.1727179596,2);
  //  double cs20 = 2*pow(2.5844403427,2); double   cs21 = pow(0.3539741687,2); double cs22 = 2*pow(2.2940148014,2);
  //  double cs23 = 2*pow(1.8135779397,2); double cs24 = 2*pow(3.6271558793,2); double cs25 = 2*pow(1.9866750947,2);
  //  double cs26 = 2*pow(9.3183321738,2); double cs27 = 2*pow(2.6899707945,2); double   cs28 = pow(0.3802292509,2);
  //  double cs29 = 2*pow(0.3556718963,2); double cs30 = 2*pow(0.8712146618,2); double cs31 = 2*pow(0.6160417952,2);
  //  double cs32 = 2*pow(4.0863589798,2); double cs33 = 2*pow(2.0431794899,2); double cs34 = 2*pow(10.418212089,2);
  //  double cs35 = 2*pow(2.7843843014,2); double   cs36 = pow(0.0505981185,2); double cs37 = 2*pow(0.4293392727,2);
  //  double cs38 = 2*pow(1.7960550366,2); double cs39 = 2*pow(4.8637400313,2); double cs40 = 2*pow(1.8837184141,2);
  //  double cs41 = 2*pow(13.583686661,2); double cs42 = 2*pow(2.0960083567,2); double cs43 = 2*pow(11.480310577,2);
  //  double cs44 = 2*pow(2.8700776442,2); double   cs45 = pow(0.0534917379,2); double cs46 = 2*pow(0.2537335916,2);
  //  double cs47 = 2*pow(2.3802320735,2); double cs48 = 2*pow(1.8179322747,2); double cs49 = 2*pow(16.055543121,2);
  //  double cs50 = 2*pow(1.9190044477,2); double cs51 = 2*pow(4.9548481782,2); double cs52 = 2*pow(2.1455121971,2);
  //  double cs53 = 2*pow(12.510378411,2); double cs54 = 2*pow(2.9487244699,2);
  double cs0=2.4674011003; double cs1=7.4022033011; double cs2=7.4022033005;
  double cs3=3.0842513755; double cs4=37.0110165048; double cs5=9.2527541262;
  double cs6=4.3179519254; double cs7=6.4769278880; double cs8=64.7692788826;
  double cs9=10.7948798139; double cs10=0.3469782797; double cs11=13.8791311886;
  double cs12=6.9395655942; double cs13=97.1539183221; double cs14=12.1442397908;
  double cs15=0.4240845642; double cs16=6.3612684614; double cs17=178.1155169268;
  double cs18=7.4214798715; double cs19=133.5866376943; double cs20=13.3586637700;
  double cs21=0.1252977121; double cs22=10.5250078181; double cs23=6.5781298867;
  double cs24=26.3125195455; double cs25=7.8937558638; double cs26=173.6626290026;
  double cs27=14.4718857505; double cs28=0.1445742832; double cs29=0.2530049956;
  double cs30=1.5180299739; double cs31=0.7590149869; double cs32=33.3966594236;
  double cs33=8.3491648559; double cs34=217.0782862628; double cs35=15.5055918758;
  double cs36=0.0025601696; double cs37=0.3686644222; double cs38=6.4516273890;
  double cs39=47.3119341841; double cs40=7.0967901272; double cs41=369.0330866085;
  double cs42=8.7865020627; double cs43=263.5950618888; double cs44=16.4746913675;
  double cs45=0.0028613660; double cs46=0.1287614710; double cs47=11.3310094474;
  double cs48=6.6097555108; double cs49=515.5609298206; double cs50=7.3651561406;
  double cs51=49.1010409380; double cs52=9.2064451758; double cs53=313.0191359728;
  double cs54=17.3899519988;

  // The power spectrum is multiplied by an l-dependent prefactor that comes
  // from the normalization of the Wigner D matrices. This prefactor is
  // mentioned in the arrata of the original SOAP paper: On representing
  // chemical environments, Phys. Rev. B 87, 184115 (2013). Here the square
  // root of the prefactor in the dot-product kernel is used, so that after a
  // possible dot-product the full prefactor is recovered.

  // SUM M's UP!
  double prel0 = PI*sqrt(8.0/(1.0));
  for(int i = 0; i < Hs; i++){
    shiftT = 0;
    for(int j = 0; j < Ts; j++){
      for(int jd = j; jd < Ts; jd++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            soapMat[NsNsLmaxTs*i + NsNsLmax*shiftT + 0 + shiftN] = prel0*(
                cs0*Cnnd[NsTs100*i + Ns100*j + 0 + k]*Cnnd[NsTs100*i + Ns100*jd + 0 + kd]);
            shiftN++;
          }
        }
        shiftT++;
      }
    }
  } if(lMax > 0){
    double prel1 = PI*sqrt(8.0/(2.0*1.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ NsNs + shiftN] = prel1*(
                  cs1*Cnnd[NsTs100*i + Ns100*j + 1*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 1*Ns + kd]
                  +cs2*Cnnd[NsTs100*i + Ns100*j + 2*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 2*Ns + kd]
                  +cs2*Cnnd[NsTs100*i + Ns100*j + 3*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 3*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 1){
    double prel2 = PI*sqrt(8.0/(2.0*2.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 2*NsNs + shiftN] = prel2*(
                  cs3*Cnnd[NsTs100*i + Ns100*j + 4*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 4*Ns + kd]
                  +cs4*Cnnd[NsTs100*i + Ns100*j + 5*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 5*Ns + kd]
                  +cs4*Cnnd[NsTs100*i + Ns100*j + 6*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 6*Ns + kd]
                  +cs5*Cnnd[NsTs100*i + Ns100*j + 7*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 7*Ns + kd]
                  +cs5*Cnnd[NsTs100*i + Ns100*j + 8*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 8*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 2){
    double prel3 = PI*sqrt(8.0/(2.0*3.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 3*NsNs + shiftN] = prel3*(
                  cs6*Cnnd[NsTs100*i + Ns100*j + 9*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 9*Ns + kd]
                  +cs7*Cnnd[NsTs100*i + Ns100*j + 10*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 10*Ns + kd]
                  +cs7*Cnnd[NsTs100*i + Ns100*j + 11*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 11*Ns + kd]
                  +cs8*Cnnd[NsTs100*i + Ns100*j + 12*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 12*Ns + kd]
                  +cs8*Cnnd[NsTs100*i + Ns100*j + 13*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 13*Ns + kd]
                  +cs9*Cnnd[NsTs100*i + Ns100*j + 14*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 14*Ns + kd]
                  +cs9*Cnnd[NsTs100*i + Ns100*j + 15*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 15*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 3){
    double prel4 = PI*sqrt(8.0/(2.0*4.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 4*NsNs + shiftN] = prel4*(
                  cs10*Cnnd[NsTs100*i + Ns100*j + 16*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 16*Ns + kd]
                  +cs11*Cnnd[NsTs100*i + Ns100*j + 17*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 17*Ns + kd]
                  +cs11*Cnnd[NsTs100*i + Ns100*j + 18*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 18*Ns + kd]
                  +cs12*Cnnd[NsTs100*i + Ns100*j + 19*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 19*Ns + kd]
                  +cs12*Cnnd[NsTs100*i + Ns100*j + 20*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 20*Ns + kd]
                  +cs13*Cnnd[NsTs100*i + Ns100*j + 21*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 21*Ns + kd]
                  +cs13*Cnnd[NsTs100*i + Ns100*j + 22*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 22*Ns + kd]
                  +cs14*Cnnd[NsTs100*i + Ns100*j + 23*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 23*Ns + kd]
                  +cs14*Cnnd[NsTs100*i + Ns100*j + 24*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 24*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 4) {
    double prel5 = PI*sqrt(8.0/(2.0*5.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 5*NsNs + shiftN] = prel5*(
                  cs15*Cnnd[NsTs100*i + Ns100*j + 25*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 25*Ns + kd]
                  +cs16*Cnnd[NsTs100*i + Ns100*j + 26*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 26*Ns + kd]
                  +cs16*Cnnd[NsTs100*i + Ns100*j + 27*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 27*Ns + kd]
                  +cs17*Cnnd[NsTs100*i + Ns100*j + 28*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 28*Ns + kd]
                  +cs17*Cnnd[NsTs100*i + Ns100*j + 29*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 29*Ns + kd]
                  +cs18*Cnnd[NsTs100*i + Ns100*j + 30*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 30*Ns + kd]
                  +cs18*Cnnd[NsTs100*i + Ns100*j + 31*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 31*Ns + kd]
                  +cs19*Cnnd[NsTs100*i + Ns100*j + 32*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 32*Ns + kd]
                  +cs19*Cnnd[NsTs100*i + Ns100*j + 33*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 33*Ns + kd]
                  +cs20*Cnnd[NsTs100*i + Ns100*j + 34*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 34*Ns + kd]
                  +cs20*Cnnd[NsTs100*i + Ns100*j + 35*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 35*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 5){
    double prel6 = PI*sqrt(8.0/(2.0*6.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 6*NsNs + shiftN] = prel6*(
                  cs21*Cnnd[NsTs100*i + Ns100*j + 36*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 36*Ns + kd]
                  +cs22*Cnnd[NsTs100*i + Ns100*j + 37*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 37*Ns + kd]
                  +cs22*Cnnd[NsTs100*i + Ns100*j + 38*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 38*Ns + kd]
                  +cs23*Cnnd[NsTs100*i + Ns100*j + 39*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 39*Ns + kd]
                  +cs23*Cnnd[NsTs100*i + Ns100*j + 40*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 40*Ns + kd]
                  +cs24*Cnnd[NsTs100*i + Ns100*j + 41*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 41*Ns + kd]
                  +cs24*Cnnd[NsTs100*i + Ns100*j + 42*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 42*Ns + kd]
                  +cs25*Cnnd[NsTs100*i + Ns100*j + 43*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 43*Ns + kd]
                  +cs25*Cnnd[NsTs100*i + Ns100*j + 44*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 44*Ns + kd]
                  +cs26*Cnnd[NsTs100*i + Ns100*j + 45*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 45*Ns + kd]
                  +cs26*Cnnd[NsTs100*i + Ns100*j + 46*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 46*Ns + kd]
                  +cs27*Cnnd[NsTs100*i + Ns100*j + 47*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 47*Ns + kd]
                  +cs27*Cnnd[NsTs100*i + Ns100*j + 48*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 48*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 6){
    double prel7 = PI*sqrt(8.0/(2.0*7.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 7*NsNs + shiftN] = prel7*(
                  cs28*Cnnd[NsTs100*i + Ns100*j + 49*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 49*Ns + kd]
                  +cs29*Cnnd[NsTs100*i + Ns100*j + 50*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 50*Ns + kd]
                  +cs29*Cnnd[NsTs100*i + Ns100*j + 51*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 51*Ns + kd]
                  +cs30*Cnnd[NsTs100*i + Ns100*j + 52*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 52*Ns + kd]
                  +cs30*Cnnd[NsTs100*i + Ns100*j + 53*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 53*Ns + kd]
                  +cs31*Cnnd[NsTs100*i + Ns100*j + 54*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 54*Ns + kd]
                  +cs31*Cnnd[NsTs100*i + Ns100*j + 55*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 55*Ns + kd]
                  +cs32*Cnnd[NsTs100*i + Ns100*j + 56*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 56*Ns + kd]
                  +cs32*Cnnd[NsTs100*i + Ns100*j + 57*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 57*Ns + kd]
                  +cs33*Cnnd[NsTs100*i + Ns100*j + 58*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 58*Ns + kd]
                  +cs33*Cnnd[NsTs100*i + Ns100*j + 59*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 59*Ns + kd]
                  +cs34*Cnnd[NsTs100*i + Ns100*j + 60*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 60*Ns + kd]
                  +cs34*Cnnd[NsTs100*i + Ns100*j + 61*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 61*Ns + kd]
                  +cs35*Cnnd[NsTs100*i + Ns100*j + 62*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 62*Ns + kd]
                  +cs35*Cnnd[NsTs100*i + Ns100*j + 63*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 63*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 7){
    double prel8 = PI*sqrt(8.0/(2.0*8.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 8*NsNs + shiftN] = prel8*(
                  cs36*Cnnd[NsTs100*i + Ns100*j + 64*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 64*Ns + kd]
                  +cs37*Cnnd[NsTs100*i + Ns100*j + 65*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 65*Ns + kd]
                  +cs37*Cnnd[NsTs100*i + Ns100*j + 66*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 66*Ns + kd]
                  +cs38*Cnnd[NsTs100*i + Ns100*j + 67*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 67*Ns + kd]
                  +cs38*Cnnd[NsTs100*i + Ns100*j + 68*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 68*Ns + kd]
                  +cs39*Cnnd[NsTs100*i + Ns100*j + 69*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 69*Ns + kd]
                  +cs39*Cnnd[NsTs100*i + Ns100*j + 70*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 70*Ns + kd]
                  +cs40*Cnnd[NsTs100*i + Ns100*j + 71*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 71*Ns + kd]
                  +cs40*Cnnd[NsTs100*i + Ns100*j + 72*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 72*Ns + kd]
                  +cs41*Cnnd[NsTs100*i + Ns100*j + 73*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 73*Ns + kd]
                  +cs41*Cnnd[NsTs100*i + Ns100*j + 74*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 74*Ns + kd]
                  +cs42*Cnnd[NsTs100*i + Ns100*j + 75*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 75*Ns + kd]
                  +cs42*Cnnd[NsTs100*i + Ns100*j + 76*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 76*Ns + kd]
                  +cs43*Cnnd[NsTs100*i + Ns100*j + 77*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 77*Ns + kd]
                  +cs43*Cnnd[NsTs100*i + Ns100*j + 78*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 78*Ns + kd]
                  +cs44*Cnnd[NsTs100*i + Ns100*j + 79*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 79*Ns + kd]
                  +cs44*Cnnd[NsTs100*i + Ns100*j + 80*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 80*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }  if(lMax > 8){
    double prel9 = PI*sqrt(8.0/(2.0*9.0+1.0));
    for(int i = 0; i < Hs; i++){
      shiftT = 0;
      for(int j = 0; j < Ts; j++){
        for(int jd = j; jd < Ts; jd++){
          shiftN = 0;
          for(int k = 0; k < Ns; k++){
            for(int kd = k; kd < Ns; kd++){
              soapMat[NsNsLmaxTs*i+NsNsLmax*shiftT+ 9*NsNs + shiftN] = prel9*(
                  cs45*Cnnd[NsTs100*i + Ns100*j + 81*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 81*Ns + kd]
                  +cs46*Cnnd[NsTs100*i + Ns100*j + 82*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 82*Ns + kd]
                  +cs46*Cnnd[NsTs100*i + Ns100*j + 83*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 83*Ns + kd]
                  +cs47*Cnnd[NsTs100*i + Ns100*j + 84*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 84*Ns + kd]
                  +cs47*Cnnd[NsTs100*i + Ns100*j + 85*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 85*Ns + kd]
                  +cs48*Cnnd[NsTs100*i + Ns100*j + 86*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 86*Ns + kd]
                  +cs48*Cnnd[NsTs100*i + Ns100*j + 87*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 87*Ns + kd]
                  +cs49*Cnnd[NsTs100*i + Ns100*j + 88*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 88*Ns + kd]
                  +cs49*Cnnd[NsTs100*i + Ns100*j + 89*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 89*Ns + kd]
                  +cs50*Cnnd[NsTs100*i + Ns100*j + 90*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 90*Ns + kd]
                  +cs50*Cnnd[NsTs100*i + Ns100*j + 91*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 91*Ns + kd]
                  +cs51*Cnnd[NsTs100*i + Ns100*j + 92*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 92*Ns + kd]
                  +cs51*Cnnd[NsTs100*i + Ns100*j + 93*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 93*Ns + kd]
                  +cs52*Cnnd[NsTs100*i + Ns100*j + 94*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 94*Ns + kd]
                  +cs52*Cnnd[NsTs100*i + Ns100*j + 95*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 95*Ns + kd]
                  +cs53*Cnnd[NsTs100*i + Ns100*j + 96*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 96*Ns + kd]
                  +cs53*Cnnd[NsTs100*i + Ns100*j + 97*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 97*Ns + kd]
                  +cs54*Cnnd[NsTs100*i + Ns100*j + 98*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 98*Ns + kd]
                  +cs54*Cnnd[NsTs100*i + Ns100*j + 99*Ns + k]*Cnnd[NsTs100*i + Ns100*jd + 99*Ns + kd]);
              shiftN++;
            }
          }
          shiftT++;
        }
      }
    }
  }
}

//===========================================================================================
void soapGTO(py::array_t<double> cArr, py::array_t<double> positions, py::array_t<double> HposArr, py::array_t<double> alphasArr, py::array_t<double> betasArr, py::array_t<int> atomicNumbersArr, double rCut, double cutoffPadding, int totalAN, int Nt, int Ns, int lMax, int Hs, double eta, bool crossover) {

  auto atomicNumbers = atomicNumbersArr.unchecked<1>();
  double *c = (double*)cArr.request().ptr;
  double *Hpos = (double*)HposArr.request().ptr;
  double *alphas = (double*)alphasArr.request().ptr;
  double *betas = (double*)betasArr.request().ptr;

  double oOeta = 1.0/eta;
  double oOeta3O2 = sqrt(oOeta*oOeta*oOeta);

  double NsNs = Ns*Ns;
  double* dx  = (double*) malloc(sizeof(double)*totalAN);
  double* dy  = (double*) malloc(sizeof(double)*totalAN);
  double* dz  = (double*) malloc(sizeof(double)*totalAN);
  double* z2 = (double*) malloc(sizeof(double)*totalAN);
  double* z4 = (double*) malloc(sizeof(double)*totalAN);
  double* z6 = (double*) malloc(sizeof(double)*totalAN);
  double* z8 = (double*) malloc(sizeof(double)*totalAN);
  double* r2 = (double*) malloc(sizeof(double)*totalAN);
  double* r4 = (double*) malloc(sizeof(double)*totalAN);
  double* r6 = (double*) malloc(sizeof(double)*totalAN);
  double* r8 = (double*) malloc(sizeof(double)*totalAN);
  double* ReIm2 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm3 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm4 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm5 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm6 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm7 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm8 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm9 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* exes = (double*) malloc (sizeof(double)*totalAN);
  double* preCoef = (double*) malloc(96*sizeof(double)*totalAN);
  double* bOa = (double*) malloc((lMax+1)*NsNs*sizeof(double));
  double* aOa = (double*) malloc((lMax+1)*Ns*sizeof(double));

  int Nx2 = 2*Ns; int Nx3 = 3*Ns; int Nx4 = 4*Ns; int Nx5 = 5*Ns;
  int Nx6 = 6*Ns; int Nx7 = 7*Ns; int Nx8 = 8*Ns; int Nx9 = 9*Ns;
  int Nx10 = 10*Ns; int Nx11 = 11*Ns; int Nx12 = 12*Ns; int Nx13 = 13*Ns;
  int Nx14 = 14*Ns; int Nx15 = 15*Ns; int Nx16 = 16*Ns; int Nx17 = 17*Ns;
  int Nx18 = 18*Ns; int Nx19 = 19*Ns; int Nx20 = 20*Ns; int Nx21 = 21*Ns;
  int Nx22 = 22*Ns; int Nx23 = 23*Ns; int Nx24 = 24*Ns; int Nx25 = 25*Ns;
  int Nx26 = 26*Ns; int Nx27 = 27*Ns; int Nx28 = 28*Ns; int Nx29 = 29*Ns;
  int Nx30 = 30*Ns; int Nx31 = 31*Ns; int Nx32 = 32*Ns; int Nx33 = 33*Ns;
  int Nx34 = 34*Ns; int Nx35 = 35*Ns; int Nx36 = 36*Ns; int Nx37 = 37*Ns;
  int Nx38 = 38*Ns; int Nx39 = 39*Ns; int Nx40 = 40*Ns; int Nx41 = 41*Ns;
  int Nx42 = 42*Ns; int Nx43 = 43*Ns; int Nx44 = 44*Ns; int Nx45 = 45*Ns;
  int Nx46 = 46*Ns; int Nx47 = 47*Ns; int Nx48 = 48*Ns; int Nx49 = 49*Ns;
  int Nx50 = 50*Ns; int Nx51 = 51*Ns; int Nx52 = 52*Ns; int Nx53 = 53*Ns;
  int Nx54 = 54*Ns; int Nx55 = 55*Ns; int Nx56 = 56*Ns; int Nx57 = 57*Ns;
  int Nx58 = 58*Ns; int Nx59 = 59*Ns; int Nx60 = 60*Ns; int Nx61 = 61*Ns;
  int Nx62 = 62*Ns; int Nx63 = 63*Ns; int Nx64 = 64*Ns; int Nx65 = 65*Ns;
  int Nx66 = 66*Ns; int Nx67 = 67*Ns; int Nx68 = 68*Ns; int Nx69 = 69*Ns;
  int Nx70 = 70*Ns; int Nx71 = 71*Ns; int Nx72 = 72*Ns; int Nx73 = 73*Ns;
  int Nx74 = 74*Ns; int Nx75 = 75*Ns; int Nx76 = 76*Ns; int Nx77 = 77*Ns;
  int Nx78 = 78*Ns; int Nx79 = 79*Ns; int Nx80 = 80*Ns; int Nx81 = 81*Ns;
  int Nx82 = 82*Ns; int Nx83 = 83*Ns; int Nx84 = 84*Ns; int Nx85 = 85*Ns;
  int Nx86 = 86*Ns; int Nx87 = 87*Ns; int Nx88 = 88*Ns; int Nx89 = 89*Ns;
  int Nx90 = 90*Ns; int Nx91 = 91*Ns; int Nx92 = 92*Ns; int Nx93 = 93*Ns;
  int Nx94 = 94*Ns; int Nx95 = 95*Ns; int Nx96 = 96*Ns; int Nx97 = 97*Ns;
  int Nx98 = 98*Ns; int Nx99 = 99*Ns;
  int t2 = 2*totalAN;  int t3 = 3*totalAN;  int t4 = 4*totalAN;
  int t5 = 5*totalAN;  int t6 = 6*totalAN;  int t7 = 7*totalAN;
  int t8 = 8*totalAN;  int t9 = 9*totalAN;  int t10 = 10*totalAN;
  int t11 = 11*totalAN;  int t12 = 12*totalAN;  int t13 = 13*totalAN;
  int t14 = 14*totalAN;  int t15 = 15*totalAN;  int t16 = 16*totalAN;
  int t17 = 17*totalAN;  int t18 = 18*totalAN;  int t19 = 19*totalAN;
  int t20 = 20*totalAN;  int t21 = 21*totalAN;  int t22 = 22*totalAN;
  int t23 = 23*totalAN;  int t24 = 24*totalAN;  int t25 = 25*totalAN;
  int t26 = 26*totalAN;  int t27 = 27*totalAN;  int t28 = 28*totalAN;
  int t29 = 29*totalAN;  int t30 = 30*totalAN;  int t31 = 31*totalAN;
  int t32 = 32*totalAN;  int t33 = 33*totalAN;  int t34 = 34*totalAN;
  int t35 = 35*totalAN;  int t36 = 36*totalAN;  int t37 = 37*totalAN;
  int t38 = 38*totalAN;  int t39 = 39*totalAN;  int t40 = 40*totalAN;
  int t41 = 41*totalAN;  int t42 = 42*totalAN;  int t43 = 43*totalAN;
  int t44 = 44*totalAN;  int t45 = 45*totalAN;  int t46 = 46*totalAN;
  int t47 = 47*totalAN;  int t48 = 48*totalAN;  int t49 = 49*totalAN;
  int t50 = 50*totalAN;  int t51 = 51*totalAN;  int t52 = 52*totalAN;
  int t53 = 53*totalAN;  int t54 = 54*totalAN;  int t55 = 55*totalAN;
  int t56 = 56*totalAN;  int t57 = 57*totalAN;  int t58 = 58*totalAN;
  int t59 = 59*totalAN;  int t60 = 60*totalAN;  int t61 = 61*totalAN;
  int t62 = 62*totalAN;  int t63 = 63*totalAN;  int t64 = 64*totalAN;
  int t65 = 65*totalAN;  int t66 = 66*totalAN;  int t67 = 67*totalAN;
  int t68 = 68*totalAN;  int t69 = 69*totalAN;  int t70 = 70*totalAN;
  int t71 = 71*totalAN;  int t72 = 72*totalAN;  int t73 = 73*totalAN;
  int t74 = 74*totalAN;  int t75 = 75*totalAN;  int t76 = 76*totalAN;
  int t77 = 77*totalAN;  int t78 = 78*totalAN;  int t79 = 79*totalAN;
  int t80 = 80*totalAN;  int t81 = 81*totalAN;  int t82 = 82*totalAN;
  int t83 = 83*totalAN;  int t84 = 84*totalAN;  int t85 = 85*totalAN;
  int t86 = 86*totalAN;  int t87 = 87*totalAN;  int t88 = 88*totalAN;
  int t89 = 89*totalAN;  int t90 = 90*totalAN;  int t91 = 91*totalAN;
  int t92 = 92*totalAN;  int t93 = 93*totalAN;  int t94 = 94*totalAN;
  int t95 = 95*totalAN;  int t96 = 96*totalAN;  int t97 = 97*totalAN;
  int t98 = 98*totalAN;  int t99 = 99*totalAN;

  double* cnnd = (double*) malloc(100*Nt*Ns*Hs*sizeof(double));
  for(int i = 0; i < 100*Nt*Ns*Hs; i++){cnnd[i] = 0.0;}

  // Initialize binning
  CellList cellList(positions, rCut+cutoffPadding);

  // Create a mapping between an atomic index and its internal index in the
  // output
  map<int, int> ZIndexMap;
  set<int> atomicNumberSet;
  for (int i = 0; i < totalAN; ++i) {
      atomicNumberSet.insert(atomicNumbers(i));
  };
  int i = 0;
  for (auto it=atomicNumberSet.begin(); it!=atomicNumberSet.end(); ++it) {
      ZIndexMap[*it] = i;
      ++i;
  };

  getAlphaBeta(aOa,bOa,alphas,betas,Ns,lMax,oOeta, oOeta3O2);

  // Loop through the centers
  for (int i = 0; i < Hs; i++) {

    // Get all neighbours for the central atom i
    double ix = Hpos[3*i];
    double iy = Hpos[3*i+1];
    double iz = Hpos[3*i+2];
    CellListResult result = cellList.getNeighboursForPosition(ix, iy, iz);

    // Sort the neighbours by type
    map<int, vector<int>> atomicTypeMap;
    for (const int &idx : result.indices) {
        int Z = atomicNumbers(idx);
        atomicTypeMap[Z].push_back(idx);
    };

    // Loop through neighbours sorted by type
    for (const auto &ZIndexPair : atomicTypeMap) {

      // j is the internal index for this atomic number
      int j = ZIndexMap[ZIndexPair.first];
      int n_neighbours = ZIndexPair.second.size();

      // Save the neighbour distances into the arrays dx, dy and dz
      getDeltas(dx, dy, dz, positions, ix, iy, iz, ZIndexPair.second);

      getRsZs(dx, dy, dz, r2, r4, r6, r8, z2, z4, z6, z8, n_neighbours);
      getCfactors(preCoef, n_neighbours, dx, dy, dz, z2, z4, z6, z8, r2, r4, r6, r8, ReIm2, ReIm3, ReIm4, ReIm5, ReIm6, ReIm7, ReIm8, ReIm9, totalAN, lMax, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99);
      getC(cnnd, preCoef, dx, dy, dz, r2, bOa, aOa, exes, totalAN, n_neighbours, Ns, Nt, lMax, i, j, Nx2, Nx3, Nx4, Nx5, Nx6, Nx7, Nx8, Nx9, Nx10, Nx11, Nx12, Nx13, Nx14, Nx15, Nx16, Nx17, Nx18, Nx19, Nx20, Nx21, Nx22, Nx23, Nx24, Nx25, Nx26, Nx27, Nx28, Nx29, Nx30, Nx31, Nx32, Nx33, Nx34, Nx35, Nx36, Nx37, Nx38, Nx39, Nx40, Nx41, Nx42, Nx43, Nx44, Nx45, Nx46, Nx47, Nx48, Nx49, Nx50, Nx51, Nx52, Nx53, Nx54, Nx55, Nx56, Nx57, Nx58, Nx59, Nx60, Nx61, Nx62, Nx63, Nx64, Nx65, Nx66, Nx67, Nx68, Nx69, Nx70, Nx71, Nx72, Nx73, Nx74, Nx75, Nx76, Nx77, Nx78, Nx79, Nx80, Nx81, Nx82, Nx83, Nx84, Nx85, Nx86, Nx87, Nx88, Nx89, Nx90, Nx91, Nx92, Nx93, Nx94, Nx95, Nx96, Nx97, Nx98, Nx99, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99);
    }
  }

  free(dx);
  free(dy);
  free(dz);
  free(z2);
  free(z4);
  free(z6);
  free(z8);
  free(r2);
  free(r4);
  free(r6);
  free(r8);
  free(ReIm2);
  free(ReIm3);
  free(ReIm4);
  free(ReIm5);
  free(ReIm6);
  free(ReIm7);
  free(ReIm8);
  free(ReIm9);
  free(exes);
  free(preCoef);
  free(bOa);
  free(aOa);

  if (crossover) {
    getPCrossOver(c, cnnd, Ns, Nt, Hs, lMax);
  } else {
    getPNoCross(c, cnnd, Ns, Nt, Hs, lMax);
  };
  free(cnnd);

  return;
}
