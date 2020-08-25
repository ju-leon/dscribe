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
inline void getRsZs(double* x,double* x2,double* x4,double* x6,double* x8,double* x10, double* y,double* y2,double* y4,double* y6,double* y8,double* y10, double* z,double* r2,double* r4,double* r6,double* r8,double* r10,double* z2,double* z4,double* z6,double* z8,double* z10, int size){
  for(int i = 0; i < size; i++){
    r2[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
    r4[i] = r2[i]*r2[i]; r6[i] = r2[i]*r4[i]; r8[i] = r4[i]*r4[i]; r10[i] = r6[i]*r4[i];
    x2[i] = x[i]*x[i]; x4[i] = x2[i]*x2[i]; x6[i] = x2[i]*x4[i]; x8[i] = x4[i]*x4[i];x10[i] = x6[i]*x4[i];
    y2[i] = y[i]*y[i]; y4[i] = y2[i]*y2[i]; y6[i] = y2[i]*y4[i]; y8[i] = y4[i]*y4[i];y10[i] = y6[i]*y4[i];
    z2[i] = z[i]*z[i]; z4[i] = z2[i]*z2[i]; z6[i] = z2[i]*z4[i]; z8[i] = z4[i]*z4[i];z10[i] = z6[i]*z4[i];
  }
}
//================================================================
void getAlphaBeta(double* aOa, double* bOa, double* alphas, double* betas, int Ns,int lMax, double oOeta, double oOeta3O2){

  int  NsNs = Ns*Ns;
  double  oneO1alpha;  double  oneO1alpha2; double  oneO1alpha3;
  double  oneO1alpha4; double  oneO1alpha5; double  oneO1alpha6;
  double  oneO1alpha7; double  oneO1alpha8; double  oneO1alpha9;
  double  oneO1alpha10;
  double  oneO1alpha11;
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
  }
  if(lMax > 8){
    int shift1 = 9*Ns; int shift2 = 9*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_9k
      oneO1alpha10 = pow(oneO1alpha,10); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha10;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_9nk
    }
  }
  if(lMax > 9){ //OBS!!!!! lMax > 9
    int shift1 = 10*Ns; int shift2 = 10*NsNs;
    for(int k = 0; k < Ns; k++){
      oneO1alpha = 1.0/(1.0 + oOeta*alphas[shift1 + k]); oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[shift1 + k] = -alphas[shift1 + k]*oneO1alpha; //got alpha_9k
      oneO1alpha11 = pow(oneO1alpha,11); oneO1alphaSqrtX = oneO1alphaSqrt*oneO1alpha11;
      for(int n = 0; n < Ns; n++){bOa[shift2 + n*Ns + k] = oOeta3O2*betas[shift2 + n*Ns + k]*oneO1alphaSqrtX;} // got beta_9nk
    }
  }
}
//================================================================
void getCfactors(double* preCoef, int Asize, double* x,double* x2, double* x4, double* x6, double* x8, double* x10, double* y,double* y2, double* y4, double* y6, double* y8, double* y10, double* z, double* z2, double* z4, double* z6, double* z8, double* z10, double* r2, double* r4, double* r6, double* r8,double* r10, double* ReIm2, double* ReIm3, double* ReIm4, double* ReIm5, double* ReIm6, double* ReIm7, double* ReIm8, double* ReIm9,int totalAN, int lMax){
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
    /*c21Im*/  preCoef[totalAN*2+i] = z[i]*y[i];
    /*c22Re*/  preCoef[totalAN*3+i] =      ReIm2[2*i];
    /*c22Im*/  preCoef[totalAN*4+i] =      ReIm2[i2+1];
    if (lMax > 2){
      /*c30  */  preCoef[totalAN*5+i] = c30c*z[i];
      /*c31Re*/  preCoef[totalAN*6+i] =       x[i]*c31c;
      /*c31Im*/  preCoef[totalAN*7+i] =       y[i]*c31c;
      /*c32Re*/  preCoef[totalAN*8+i] = z[i]*ReIm2[i2];
      /*c32Im*/  preCoef[totalAN*9+i] = z[i]*ReIm2[i2+1];
      /*c33Re*/  preCoef[totalAN*10+i] =      ReIm3[i2  ];
      /*c33Im*/  preCoef[totalAN*11+i] =      ReIm3[i2+1];
    }
    if (lMax > 3){
      /*c40  */  preCoef[totalAN*12+i] = c40c;
      /*c41Re*/  preCoef[totalAN*13+i] = z[i]*x[i]*c41c;
      /*c41Im*/  preCoef[totalAN*14+i] = z[i]*y[i]*c41c;
      /*c42Re*/  preCoef[totalAN*15+i] =      ReIm2[i2  ]*c42c;
      /*c42Im*/  preCoef[totalAN*16+i] =      ReIm2[i2+1]*c42c;
      /*c43Re*/  preCoef[totalAN*17+i] = z[i]*ReIm3[i2  ];
      /*c43Im*/  preCoef[totalAN*18+i] = z[i]*ReIm3[i2+1];
      /*c44Re*/  preCoef[totalAN*19+i] =      ReIm4[i2  ];
      /*c44Im*/  preCoef[totalAN*20+i] =      ReIm4[i2+1];
    }
    if (lMax > 4){
      /*c50  */  preCoef[totalAN*21+i] = c50c*z[i];
      /*c51Re*/  preCoef[totalAN*22+i] =      x[i]*c51c;
      /*c51Im*/  preCoef[totalAN*23+i] =      y[i]*c51c;
      /*c52Re*/  preCoef[totalAN*24+i] = z[i]*ReIm2[i2  ]*c52c;
      /*c52Im*/  preCoef[totalAN*25+i] = z[i]*ReIm2[i2+1]*c52c;
      /*c53Re*/  preCoef[totalAN*26+i] =      ReIm3[i2  ]*c53c;
      /*c53Im*/  preCoef[totalAN*27+i] =      ReIm3[i2+1]*c53c;
      /*c54Re*/  preCoef[totalAN*28+i] = z[i]*ReIm4[i2  ];
      /*c54Im*/  preCoef[totalAN*29+i] = z[i]*ReIm4[i2+1];
      /*c55Re*/  preCoef[totalAN*30+i] =      ReIm5[i2  ];
      /*c55Im*/  preCoef[totalAN*31+i] =      ReIm5[i2+1];
    }
    if (lMax > 5){
      /*c60  */  preCoef[totalAN*32+i] = c60c;
      /*c61Re*/  preCoef[totalAN*33+i] = z[i]*x[i]*c61c;
      /*c61Im*/  preCoef[totalAN*34+i] = z[i]*y[i]*c61c;
      /*c62Re*/  preCoef[totalAN*35+i] =      ReIm2[i2  ]*c62c;
      /*c62Im*/  preCoef[totalAN*36+i] =      ReIm2[i2+1]*c62c;
      /*c63Re*/  preCoef[totalAN*37+i] = z[i]*ReIm3[i2  ]*c63c;
      /*c63Im*/  preCoef[totalAN*38+i] = z[i]*ReIm3[i2+1]*c63c;
      /*c64Re*/  preCoef[totalAN*39+i] =      ReIm4[i2  ]*c64c;
      /*c64Im*/  preCoef[totalAN*40+i] =      ReIm4[i2+1]*c64c;
      /*c65Re*/  preCoef[totalAN*41+i] = z[i]*ReIm5[i2  ];
      /*c65Im*/  preCoef[totalAN*42+i] = z[i]*ReIm5[i2+1];
      /*c66Re*/  preCoef[totalAN*43+i] =      ReIm6[i2  ];
      /*c66Im*/  preCoef[totalAN*44+i] =      ReIm6[i2+1];
    }
    if (lMax > 6){
      /*c70  */  preCoef[totalAN*45+i] = c70c*z[i];
      /*c71Re*/  preCoef[totalAN*46+i] = x[i]*c71c;
      /*c71Im*/  preCoef[totalAN*47+i] = y[i]*c71c;
      /*c72Re*/  preCoef[totalAN*48+i] = z[i]*ReIm2[i2  ]*c72c;
      /*c72Im*/  preCoef[totalAN*49+i] = z[i]*ReIm2[i2+1]*c72c;
      /*c73Re*/  preCoef[totalAN*50+i] =      ReIm3[i2  ]*c73c;
      /*c73Im*/  preCoef[totalAN*51+i] =      ReIm3[i2+1]*c73c;
      /*c74Re*/  preCoef[totalAN*52+i] = z[i]*ReIm4[i2  ]*c74c;
      /*c74Im*/  preCoef[totalAN*53+i] = z[i]*ReIm4[i2+1]*c74c;
      /*c75Re*/  preCoef[totalAN*54+i] =      ReIm5[i2  ]*c75c;
      /*c75Im*/  preCoef[totalAN*55+i] =      ReIm5[i2+1]*c75c;
      /*c76Re*/  preCoef[totalAN*56+i] = z[i]*ReIm6[i2  ];
      /*c76Im*/  preCoef[totalAN*57+i] = z[i]*ReIm6[i2+1];
      /*c77Re*/  preCoef[totalAN*58+i] =      ReIm7[i2  ];
      /*c77Im*/  preCoef[totalAN*59+i] =      ReIm7[i2+1];
    }
    if (lMax > 7){
      /*c80  */  preCoef[totalAN*60+i] = c80c;
      /*c81Re*/  preCoef[totalAN*61+i] = z[i]*x[i]*c81c;
      /*c81Im*/  preCoef[totalAN*62+i] = z[i]*y[i]*c81c;
      /*c82Re*/  preCoef[totalAN*63+i] =      ReIm2[i2  ]*c82c;
      /*c82Im*/  preCoef[totalAN*64+i] =      ReIm2[i2+1]*c82c;
      /*c83Re*/  preCoef[totalAN*65+i] = z[i]*ReIm3[i2  ]*c83c;
      /*c83Im*/  preCoef[totalAN*66+i] = z[i]*ReIm3[i2+1]*c83c;
      /*c84Re*/  preCoef[totalAN*67+i] =      ReIm4[i2  ]*c84c;
      /*c84Im*/  preCoef[totalAN*68+i] =      ReIm4[i2+1]*c84c;
      /*c85Re*/  preCoef[totalAN*69+i] = z[i]*ReIm5[i2  ]*c85c;
      /*c85Im*/  preCoef[totalAN*70+i] = z[i]*ReIm5[i2+1]*c85c;
      /*c86Re*/  preCoef[totalAN*71+i] =      ReIm6[i2  ]*c86c;
      /*c86Im*/  preCoef[totalAN*72+i] =      ReIm6[i2+1]*c86c;
      /*c87Re*/  preCoef[totalAN*73+i] = z[i]*ReIm7[i2  ];
      /*c87Im*/  preCoef[totalAN*74+i] = z[i]*ReIm7[i2+1];
      /*c88Re*/  preCoef[totalAN*75+i] =      ReIm8[i2  ];
      /*c88Im*/  preCoef[totalAN*76+i] =      ReIm8[i2+1];
    }
    if (lMax > 8){
      /*c90  */  preCoef[totalAN*77+i] = c90c*z[i];
      /*c91Re*/  preCoef[totalAN*78+i] = x[i]*c91c;
      /*c91Im*/  preCoef[totalAN*79+i] = y[i]*c91c;
      /*c92Re*/  preCoef[totalAN*80+i] = z[i]*ReIm2[i2  ]*c92c;
      /*c92Im*/  preCoef[totalAN*81+i] = z[i]*ReIm2[i2+1]*c92c;
      /*c93Re*/  preCoef[totalAN*82+i] =      ReIm3[i2  ]*c93c;
      /*c93Im*/  preCoef[totalAN*83+i] =      ReIm3[i2+1]*c93c;
      /*c94Re*/  preCoef[totalAN*84+i] = z[i]*ReIm4[i2  ]*c94c;
      /*c94Im*/  preCoef[totalAN*85+i] = z[i]*ReIm4[i2+1]*c94c;
      /*c95Re*/  preCoef[totalAN*86+i] =      ReIm5[i2  ]*c95c;
      /*c95Im*/  preCoef[totalAN*87+i] =      ReIm5[i2+1]*c95c;
      /*c96Re*/  preCoef[totalAN*88+i] = z[i]*ReIm6[i2  ]*c96c;
      /*c96Im*/  preCoef[totalAN*89+i] = z[i]*ReIm6[i2+1]*c96c;
      /*c97Re*/  preCoef[totalAN*90+i] =      ReIm7[i2  ]*c97c;
      /*c97Im*/  preCoef[totalAN*91+i] =      ReIm7[i2+1]*c97c;
      /*c98Re*/  preCoef[totalAN*92+i] = z[i]*ReIm8[i2  ];
      /*c98Im*/  preCoef[totalAN*93+i] = z[i]*ReIm8[i2+1];
      /*c99Re*/  preCoef[totalAN*94+i] =      ReIm9[i2  ];
      /*c99Im*/  preCoef[totalAN*95+i] =      ReIm9[i2+1];
    }
    if (lMax > 9){
      /**/  preCoef[totalAN*96+i] = 1.53479023644398*x[i]*y[i]*(5.0*x8[i] - 60.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 60.0*x2[i]*y6[i] + 5.0*y8[i]);
      /**/  preCoef[totalAN*97+i] = 3.43189529989171*y[i]*z[i]*(9.0*x8[i] - 84.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 36.0*x2[i]*y6[i] + y8[i]);
      /**/  preCoef[totalAN*98+i] = -4.45381546176335*x[i]*y[i]*(x2[i] + y2[i] - 18.0*z2[i])*(x6[i] - 7.0*x4[i]*y2[i] + 7.0*x2[i]*y4[i] - y6[i]);
      /**/  preCoef[totalAN*99+i] = 1.36369691122981*y[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 16.0*z2[i])*(7.0*x6[i] - 35.0*x4[i]*y2[i] + 21.0*x2[i]*y4[i] - y6[i]);
      /**/  preCoef[totalAN*100+i] = 0.330745082725238*x[i]*y[i]*(3.0*x4[i] - 10.0*x2[i]*y2[i] + 3.0*y4[i])*(323.0*z4[i] - 102.0*z2[i]*r2[i] + 3.0*r4[i]);
      /**/  preCoef[totalAN*101+i] = 0.295827395278969*y[i]*z[i]*(5.0*x4[i] - 10.0*x2[i]*y2[i] + y4[i])*(323.0*z4[i] - 170.0*z2[i]*r2[i] + 15.0*r4[i]);
      /**/  preCoef[totalAN*102+i] = 1.87097672671297*x[i]*y[i]*(x2[i] - y2[i])*(323.0*z6[i] - 255.0*z4[i]*r2[i] + 45.0*z2[i]*r4[i] - r6[i]);
      /**/  preCoef[totalAN*103+i] = 0.661490165450475*y[i]*z[i]*(3.0*x2[i] - y2[i])*(323.0*z6[i] - 357.0*z4[i]*r2[i] + 105.0*z2[i]*r4[i] - 7.0*r6[i]);
      /**/  preCoef[totalAN*104+i] = 0.129728894680065*x[i]*y[i]*(4199.0*z8[i] - 6188.0*z6[i]*r2[i] + 2730.0*z4[i]*r4[i] - 364.0*z2[i]*r6[i] + 7.0*r8[i]);
      /**/  preCoef[totalAN*105+i] = 0.0748990122652082*y[i]*z[i]*(4199.0*z8[i] - 7956.0*z6[i]*r2[i] + 4914.0*z4[i]*r4[i] - 1092.0*z2[i]*r6[i] + 63.0*r8[i]);
      /**/  preCoef[totalAN*106+i] = 233.240148813258*z10[i] - 552.410878768242*z8[i]*r2[i] + 454.926606044435*z6[i]*r4[i] - 151.642202014812*z4[i]*r6[i] + 17.4971771555552*z2[i]*r8[i] - 0.318130493737367*r10[i];
      /**/  preCoef[totalAN*107+i] = 0.0748990122652082*x[i]*z[i]*(4199.0*z8[i] - 7956.0*z6[i]*r2[i] + 4914.0*z4[i]*r4[i] - 1092.0*z2[i]*r6[i] + 63.0*r8[i]);
      /**/  preCoef[totalAN*108+i] = 0.0648644473400325*(x2[i] - y2[i])*(4199.0*z8[i] - 6188.0*z6[i]*r2[i] + 2730.0*z4[i]*r4[i] - 364.0*z2[i]*r6[i] + 7.0*r8[i]);
      /**/  preCoef[totalAN*109+i] = 0.661490165450475*x[i]*z[i]*(x2[i] - 3.0*y2[i])*(323.0*z6[i] - 357.0*z4[i]*r2[i] + 105.0*z2[i]*r4[i] - 7.0*r6[i]);
      /**/  preCoef[totalAN*110+i] = 0.467744181678242*(x4[i] - 6.0*x2[i]*y2[i] + y4[i])*(323.0*z6[i] - 255.0*z4[i]*r2[i] + 45.0*z2[i]*r4[i] - r6[i]);
      /**/  preCoef[totalAN*111+i] = 0.295827395278969*x[i]*z[i]*(x4[i] - 10.0*x2[i]*y2[i] + 5.0*y4[i])*(323.0*z4[i] - 170.0*z2[i]*r2[i] + 15.0*r4[i]);
      /**/  preCoef[totalAN*112+i] = 0.165372541362619*(323.0*z4[i] - 102.0*z2[i]*r2[i] + 3.0*r4[i])*(x6[i] - 15.0*x4[i]*y2[i] + 15.0*x2[i]*y4[i] - y6[i]);
      /**/  preCoef[totalAN*113+i] = 1.36369691122981*x[i]*z[i]*(-3.0*x2[i] - 3.0*y2[i] + 16.0*z2[i])*(x6[i] - 21.0*x4[i]*y2[i] + 35.0*x2[i]*y4[i] - 7.0*y6[i]);
      /**/  preCoef[totalAN*114+i] = -0.556726932720418*(x2[i] + y2[i] - 18.0*z2[i])*(x8[i] - 28.0*x6[i]*y2[i] + 70.0*x4[i]*y4[i] - 28.0*x2[i]*y6[i] + y8[i]);
      /**/  preCoef[totalAN*115+i] = 3.43189529989171*x[i]*z[i]*(x8[i] - 36.0*x6[i]*y2[i] + 126.0*x4[i]*y4[i] - 84.0*x2[i]*y6[i] + 9.0*y8[i]);
      /**/  preCoef[totalAN*116+i] = 0.76739511822199*x10[i] - 34.5327803199895*x8[i]*y2[i] + 161.152974826618*x6[i]*y4[i] - 161.152974826618*x4[i]*y6[i] + 34.5327803199895*x2[i]*y8[i] - 0.76739511822199*y10[i];

    }
  }
}
//================================================================
void getC(double* C, double* preCoef, double* x, double* y, double* z,double* r2, double* bOa, double* aOa, double* exes,  int totalAN, int Asize, int Ns, int Ntypes, int lMax, int posI, int typeJ){

  if(Asize == 0){return;}
  double sumMe = 0; int NsNs = Ns*Ns;  int NsJ = ((lMax+1)*(lMax+1))*Ns*typeJ; int LNsNs;
  int LNs; int NsTsI = ((lMax+1)*(lMax+1))*Ns*Ntypes*posI;
  for(int k = 0; k < Ns; k++){
    sumMe = 0; for(int i = 0; i < Asize; i++){ sumMe += exp(aOa[k]*r2[i]);}
    for(int n = 0; n < Ns; n++){ C[NsTsI + NsJ + n] += bOa[n*Ns + k]*sumMe; }
  } if(lMax > 0) { LNsNs=NsNs; LNs=Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c10*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*z[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c11Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*x[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*2 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c11Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*y[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*3 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 1) { LNsNs=2*NsNs; LNs=2*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c20*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*4 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c21Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN + i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*5 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c21Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*2+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*6 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c22Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*3+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*7 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c22Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*preCoef[totalAN*4+ i];}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*8 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
  }} if(lMax > 2) { LNsNs=3*NsNs; LNs=3*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c30*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*5+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*9 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c31Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*6+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*10 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c31Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*7+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*11 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c32Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*8+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*12 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c32Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*9+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*13 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c33Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*10+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*14 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c33Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*11+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*15 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
  }} if(lMax > 3) { LNsNs=4*NsNs; LNs=4*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c40*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*12+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*16 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c41Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*13+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*17 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c41Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*14+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*18 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c42Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*15+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*19 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c42Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*16+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*20 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c43Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*17+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*21 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c43Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*18+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*22 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c44Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*19+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*23 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }
    sumMe = 0;/*c44Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*20+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*24 + n] += bOa[LNsNs + n*Ns + k]*sumMe; }

  }} if(lMax > 4) { LNsNs=5*NsNs; LNs=5*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c50*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*21+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*25 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c51Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*22+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*26 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c51Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*23+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*27 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c52Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*24+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*28 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c52Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*25+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*29 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c53Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*26+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*30 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c53Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*27+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*31 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c54Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*28+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*32 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c54Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*29+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*33 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c55Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*30+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*34 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c55Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*31+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*35 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 5) { LNsNs=6*NsNs; LNs=6*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c60*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*32+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*36 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c61Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*33+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*37 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c61Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*34+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*38 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c62Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*35+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*39 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c62Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*36+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*40 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c63Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*37+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*41 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c63Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*38+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*42 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c64Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*39+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*43 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c64Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*40+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*44 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c65Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*41+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*45 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c65Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*42+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*46 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c66Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*43+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*47 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c66Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*44+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*48 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 6) { LNsNs=7*NsNs; LNs=7*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c70*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*45+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*49 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c71Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*46+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*50 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c71Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*47+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*51 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c72Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*48+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*52 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c72Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*49+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*53 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c73Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*50+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*54 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c73Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*51+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*55 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c74Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*52+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*56 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c74Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*53+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*57 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c75Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*54+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*58 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c75Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*55+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*59 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c76Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*56+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*60 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c76Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*57+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*61 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c77Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*58+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*62 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c77Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*59+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*63 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }} if(lMax > 7) { LNsNs=8*NsNs; LNs=8*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c80*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*60+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*64 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c81Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*61+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*65 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c81Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*62+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*66 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c82Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*63+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*67 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c82Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*64+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*68 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c83Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*65+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*69 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c83Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*66+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*70 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c84Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*67+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*71 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c84Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*68+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*72 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c85Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*69+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*73 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c85Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*70+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*74 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c86Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*71+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*75 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c86Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*72+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*76 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c87Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*73+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*77 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c87Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*74+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*78 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c88Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*75+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*79 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c88Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*76+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*80 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }}
  if(lMax > 8) { LNsNs=9*NsNs; LNs=9*Ns;
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
    sumMe = 0;/*c90*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*77+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*81 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c91Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*78+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*82 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c91Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*79+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*83 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c92Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*80+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*84 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c92Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*81+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*85 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c93Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*82+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*86 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c93Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*83+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*87 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c94Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*84+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*88 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c94Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*85+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*89 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c95Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*86+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*90 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c95Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*87+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*91 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c96Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*88+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*92 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c96Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*89+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*93 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c97Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*90+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*94 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c97Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*91+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*95 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c98Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*92+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*96 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c98Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*93+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*97 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c99Re*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*94+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*98 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
    sumMe = 0;/*c99Im*/ for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*95+i]);}
    for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*99 + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
  }}
//  double shiftBuffer = 96;
  if(lMax > 9) { LNsNs=10*NsNs; LNs=10*Ns; // OBS!!!!!! lMax > 9 Case!
  for(int k = 0; k < Ns; k++){
    for(int i = 0; i < Asize; i++){exes[i] = exp(aOa[LNs + k]*r2[i]);}//exponents
      for(int sumems = 100; sumems < 121; sumems++){
        sumMe = 0; for(int i = 0; i < Asize; i++){sumMe += exes[i]*(preCoef[totalAN*(sumems - 4)+i]);}
        for(int n = 0; n < Ns; n++){C[NsTsI + NsJ + Ns*sumems + n] += bOa[LNsNs + n*Ns + k]*sumMe;}
	//shiftBuffer++; WRONG LOGIC, but not in use anyway ( not considering k++)
      }
  }}
}
//=======================================================================
/**
 * Used to calculate the partial power spectrum without crossover.
 */
void getPNoCross(double* soapMat, double* Cnnd, int Ns, int Ts, int Hs, int lMax){
  int NsTs100 = Ns*Ts*((lMax+1)*(lMax+1)); // Used to be NsTs100 = Ns*Ts*100, but 100 is a waste of memory if not lMax = 9, and can't do over that, so changed.
  int Ns100 = Ns*((lMax+1)*(lMax+1));
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
  double cs54=17.3899519988; // constants only up to lMax = 9. After that, it is switched to tesseral, so no need.

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
            cs0*Cnnd[NsTs100*i + Ns100*j + 0 + k]*Cnnd[NsTs100*i + Ns100*j + 0*Ns + kd]);
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
   if (lMax > 9) { // OBS!!!! LMAX > 9 ------
    double prel9 = PI*sqrt(8.0/(2.0*10.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            double buffDouble = 0;
            for(int buffShift = 100; buffShift < 121; buffShift++){
              buffDouble += Cnnd[NsTs100*i + Ns100*j + buffShift*Ns + k]*Cnnd[NsTs100*i + Ns100*j + buffShift*Ns + kd];
	    }
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 10*NsNs + shiftN] = prel9*buffDouble;
            shiftN++;
          }
        }
      }
    }
  } // Correct logic? buffDouble?
}
//=======================================================================
/**
 * Used to calculate the partial power spectrum.
 */
void getPCrossOver(double* soapMat, double* Cnnd, int Ns, int Ts, int Hs, int lMax){
  int NsTs100 = Ns*Ts*((lMax+1)*(lMax+1));
  int Ns100 = Ns*((lMax+1)*(lMax+1));
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

   if (lMax > 9) { // OBS!!!! LMAX > 9 ------
    double prel9 = PI*sqrt(8.0/(2.0*10.0+1.0));
    for(int i = 0; i < Hs; i++){
      for(int j = 0; j < Ts; j++){
        shiftN = 0;
        for(int k = 0; k < Ns; k++){
          for(int kd = k; kd < Ns; kd++){
            double buffDouble = 0;
            for(int buffShift = 100; buffShift < 121; buffShift++){
              buffDouble += Cnnd[NsTs100*i + Ns100*j + buffShift*Ns + k]*Cnnd[NsTs100*i + Ns100*j + buffShift*Ns + kd];
	    }
            soapMat[NsNsLmaxTs*i+NsNsLmax*j+ 10*NsNs + shiftN] = prel9*buffDouble;
            shiftN++;
          }
        }
      }
    }
  } // Correct logic? buffDouble?
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
  double* x2 = (double*) malloc(sizeof(double)*totalAN);
  double* x4 = (double*) malloc(sizeof(double)*totalAN);
  double* x6 = (double*) malloc(sizeof(double)*totalAN);
  double* x8 = (double*) malloc(sizeof(double)*totalAN);
  double* x10 = (double*) malloc(sizeof(double)*totalAN);
  double* y2 = (double*) malloc(sizeof(double)*totalAN);
  double* y4 = (double*) malloc(sizeof(double)*totalAN);
  double* y6 = (double*) malloc(sizeof(double)*totalAN);
  double* y8 = (double*) malloc(sizeof(double)*totalAN);
  double* y10 = (double*) malloc(sizeof(double)*totalAN);
  double* z2 = (double*) malloc(sizeof(double)*totalAN);
  double* z4 = (double*) malloc(sizeof(double)*totalAN);
  double* z6 = (double*) malloc(sizeof(double)*totalAN);
  double* z8 = (double*) malloc(sizeof(double)*totalAN);
  double* z10 = (double*) malloc(sizeof(double)*totalAN);
  double* r2 = (double*) malloc(sizeof(double)*totalAN);
  double* r4 = (double*) malloc(sizeof(double)*totalAN);
  double* r6 = (double*) malloc(sizeof(double)*totalAN);
  double* r8 = (double*) malloc(sizeof(double)*totalAN);
  double* r10 = (double*) malloc(sizeof(double)*totalAN);
  double* ReIm2 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm3 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm4 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm5 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm6 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm7 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm8 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* ReIm9 = (double*) malloc(2*sizeof(double)*totalAN);// 2 -> Re + ixIm
  double* exes = (double*) malloc (sizeof(double)*totalAN);
  double* preCoef = (double*) malloc(((lMax + 1)*(lMax + 1) - 4)*sizeof(double)*totalAN);
  double* bOa = (double*) malloc((lMax+1)*NsNs*sizeof(double));
  double* aOa = (double*) malloc((lMax+1)*Ns*sizeof(double));


  double* cnnd = (double*) malloc(((lMax+1)*(lMax+1))*Nt*Ns*Hs*sizeof(double));
  for(int i = 0; i < ((lMax+1)*(lMax+1))*Nt*Ns*Hs; i++){cnnd[i] = 0.0;}

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

      getRsZs(dx,x2,x4,x6,z8,x10, dy,y2,y4,y6,y8,y10, dz, r2, r4, r6, r8,r10, z2, z4, z6, z8,z10, n_neighbours);
      getCfactors(preCoef, n_neighbours, dx,x2, x4, x6, x8,x10, dy,y2, y4, y6, y8,y10, dz, z2, z4, z6, z8,z10, r2, r4, r6, r8,r10, ReIm2, ReIm3, ReIm4, ReIm5, ReIm6, ReIm7, ReIm8, ReIm9, totalAN, lMax); // Erased tn
      getC(cnnd, preCoef, dx, dy, dz, r2, bOa, aOa, exes, totalAN, n_neighbours, Ns, Nt, lMax, i, j); //erased tn and Nx
    }
  }

  free(dx);
  free(x2);
  free(x4);
  free(x6);
  free(x8);
  free(x10);
  free(dy);
  free(y2);
  free(y4);
  free(y6);
  free(y8);
  free(y10);
  free(dz);
  free(z2);
  free(z4);
  free(z6);
  free(z8);
  free(z10);
  free(r2);
  free(r4);
  free(r6);
  free(r8);
  free(r10);
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
