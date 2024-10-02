import logging
import numpy as np

from math import pi
import matplotlib.pyplot as plt
from textwrap import indent


logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

CONSTANTS = {
    'toRAD': 180/pi,   
}

toRAD = 180/pi

J= [
    [100, 0, 0],
    [0, 200, 0],
    [0, 0, 150]
]
Jinv = [
    [1/100, 0, 0],
    [0, 1/200, 0],
    [0, 0, 1/150]
]

Ms = [0,0,0]
Mv = [0.2,-0.5,0.3]

def cross_prod_vec(vec_1, vec_2):
    """ a = [a[0], a[1], a[2]]  b = [b[0], b[1], b[2]] """
    
    return [
        vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1],
        vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2],
        vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0]
    ]

def mult_matrix_vec(matr, vec):
    #Входные параметры не проверяются. Память под вектор R должна быть выделена. Вектор R должен иметь число элементов равное числу строк матрицы A. 
    #Вектор B должен иметь число элементов равное числу столбцов матрицы A.
    if len(matr) != len(vec):
        logging.error(f'def mult_matrix_vec(): некорректные входные данные (matr = {matr}, vec = {vec})')
        raise ValueError(f'Некорректные входные данные: вектор должен иметь число элементов равное числу столбцов матрицы')
    res = []
    for row in range(len(matr)):
        tmp = 0
        for col in range(len(vec)):
            tmp += matr[row][col] * vec[col]
        res.append(tmp)
    
    return res

def div_vec(vec_1, vec_2):
    vec = []
    
    for i in range(len(vec_1)):
        vec.append(vec_1[i] - vec_2[i])
    
    return vec

def RightPart(t, w):
    Msum = []
    
    jw = mult_matrix_vec(J, w[:3])
    tmp = cross_prod_vec(w, jw)
    
    for i in range(3):
        Msum.append(Ms[i] + Mv[i])
        
    tmp2 = div_vec(Msum, tmp)
    dw = mult_matrix_vec(Jinv, tmp2)
    
    for i in range(3):
        dw.append(w[i])
    
#    for(int i = 0; i <3; i++) dw[i+3] = w[i];
    
    return dw

def RK4(t, w, h):
    POR = 6
    k1t = []
    k2t = []
    k3t = []
    wk = []
    
    
    k1 = RightPart(t,w)
    
    for i in range(POR):
        k1t.append(w[i] + h * 0.5 * k1[i])
    k2 = RightPart(t + h * 0.5, k1t)
    
    for i in range(POR):
        k2t.append(w[i] + h * 0.5 * k2[i])
    k3 = RightPart(t + h * 0.5, k2t)
    
    for i in range(POR):
        k3t.append(w[i] + h * 0.5 * k3[i])
    k4 = RightPart(t + h * 0.5, k3t)
    
    for i in range(POR):
        wk.append(w[i] + (h/6.0)*(k1[i] + 2.0*(k2[i] + k3[i]) + k4[i]))
        
    return wk
    

w_plot = [[],[],[],[],[],[]]
t_plot = []

w = [1.0/toRAD, -1.0/toRAD, 0, 2.5/toRAD, 1.0/toRAD, -3.5/toRAD]
h = 1 #0.01
t = 0
tk = 150
Kfi = 20
Kw = 40
Ki = 1.0

fi_tr = [0/toRAD, -0/toRAD, 0/toRAD]
In = [0, 0, 0]

while(t < tk):
    wk = RK4(t, w, h)
    for i in range(6):
        w[i] = wk[i]
         
    for i in range(3):
        In[i] = In[i] +(w[i+3] - fi_tr[i])*h
        Ms[i] = -(Kfi*(w[i+3] - fi_tr[i]) + Kw*w[i] + Ki*In[i])

    
    print(t, w[0]*toRAD, w[1]*toRAD, w[2]*toRAD, w[3]*toRAD, w[4]*toRAD, w[5]*toRAD)
    
    t_plot.append(t)
    
    for i in range(6):
        w_plot[i].append(w[i]*toRAD)
    

    #//Application->ProcessMessages();
    
    t = t + h
    
plt.plot(t_plot, w_plot[0], t_plot, w_plot[1], t_plot, w_plot[2])
plt.show()
"""
double w[6] = {1.0/toRAD,-1.0/toRAD,0, 2.5/toRAD, 1.0/toRAD, -3.5/toRAD},
                  wk[6], h = 0.01, t = 0, tk = 150;

   double Kfi = 20, Kw = 40, Ki = 1.0;

   double fi_tr[3] = {0/toRAD, -0/toRAD, 0/toRAD};
   double in[3] = {0,0,0};

   FILE *f = fopen("res.csv","w");

   while(t < tk)
   {
     RK4(t, w, wk, h);
     for(int i =0; i<6; i++)    w[i] = wk[i];

     for(int i =0; i<3; i++)
     {
        in[i] = in[i] +(w[i+3] - fi_tr[i])*h;
        Ms[i] = -(Kfi*(w[i+3] - fi_tr[i]) + Kw*w[i] + Ki*in[i]);
     }

     fprintf(f,"%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
             t, w[0]*toRAD, w[1]*toRAD, w[2]*toRAD,
             w[3]*toRAD, w[4]*toRAD, w[5]*toRAD);
     //Application->ProcessMessages();

     t = t + h;

"""


"""
a = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
]
b = [1,2,3]

print(np.dot(a,b))
try:
    print(mult_matrix_vec(a, b))
except Exception as err:
    logging.critical(f'def mult_matrix_vec() не отработала: "{err}"')
"""




"""
//---------------------------------------------------------------------------
#include <vcl.h>
#include <stdio.h>
#include <math.h>
#pragma hdrstop

#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
const double toRAD = 180.0/M_PI;
//---------------------------------------------------------------------------
static double J[9]=
{
  100, 0,   0,
  0,   200, 0,
  0,   0,   150
};
//---------------------------------------------------------------------------
static double Jinv[9]=
{
  1.0/100, 0,           0,
  0,       1.0/200,     0,
  0,       0,           1.0/150
};
//---------------------------------------------------------------------------
double Ms[3] = {0,0,0};
double Mv[3] = {0.2,-0.5,0.3};
//---------------------------------------------------------------------------


void mvec(double *a, double *b, double *c)
{
  // a = a[0], a[1], a[2]
  // b = b[0], b[1], b[2]
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[2]*b[0];
}
//---------------------------------------------------------------------------
void mprd3x1(double *a, double *b, double *c)
{
  // a -         3 3, b -        3 1
   c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
   c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
}
//---------------------------------------------------------------------------
void msub(double *a, double *b, double *c)
{
  for(int i=0; i<3; i++)  c[i] = a[i] - b[i];
}
//---------------------------------------------------------------------------
void RightPart(double t, double *w, double *dw)
{
   double jw[3], tmp[3], tmp2[3], Msum[3];
   mprd3x1(J,w,jw);      // Jw

   mvec(w,jw,tmp);       // w x Jw
   for(int i =0; i<3; i++) Msum[i] = Ms[i] + Mv[i];
   msub(Msum,tmp,tmp2);    // Ms - (w x Jw)
   mprd3x1(Jinv,tmp2,dw);   // dw = w`
   for(int i = 0; i <3; i++) dw[i+3] = w[i];

}
//---------------------------------------------------------------------------
void RK4(double t, double *w, double *wk, double h)
{
  double k1[6];          RightPart(t,w,k1);
  double k2[6],k1t[6];
  for(int i=0; i<6; i++) k1t[i] = w[i] + h*0.5*k1[i];
                         RightPart(t+h*0.5,k1t,k2);
  double k3[6],k2t[6];
  for(int i=0; i<6; i++) k2t[i] = w[i] + h*0.5*k2[i];
                         RightPart(t+h*0.5,k2t,k3);
  double k4[6],k3t[6];
  for(int i=0; i<6; i++) k3t[i] = w[i] + h*k3[i];
                         RightPart(t+h,k3t,k4);

  for(int i =0; i<6; i++)
        wk[i] = w[i] + (h/6.0)*(k1[i] + 2.0*(k2[i] + k3[i]) + k4[i]);  
}
//---------------------------------------------------------------------------
void __fastcall TForm1::GO1Click(TObject *Sender)
{
   double w[6] = {1.0/toRAD,-1.0/toRAD,0, 2.5/toRAD, 1.0/toRAD, -3.5/toRAD},
                  wk[6], h = 0.01, t = 0, tk = 150;

   double Kfi = 20, Kw = 40, Ki = 1.0;

   double fi_tr[3] = {0/toRAD, -0/toRAD, 0/toRAD};
   double in[3] = {0,0,0};

   FILE *f = fopen("res.csv","w");

   while(t < tk)
   {
     RK4(t, w, wk, h);
     for(int i =0; i<6; i++)    w[i] = wk[i];

     for(int i =0; i<3; i++)
     {
        in[i] = in[i] +(w[i+3] - fi_tr[i])*h;
        Ms[i] = -(Kfi*(w[i+3] - fi_tr[i]) + Kw*w[i] + Ki*in[i]);
     }

     fprintf(f,"%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
             t, w[0]*toRAD, w[1]*toRAD, w[2]*toRAD,
             w[3]*toRAD, w[4]*toRAD, w[5]*toRAD);
     //Application->ProcessMessages();

     t = t + h;
   }

   fclose(f);

}


//---------------------------------------------------------------------------
"""
