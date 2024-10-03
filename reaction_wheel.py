import logging
import numpy as np

from math import pi, cos, sin, asin, sqrt
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

Idm = 0.00016
Om = [0, 0, 0]
Ms_max = 0.05
Om_max = 600

Ms = [0,0,0]
Mv = [0, 0, 0]  #[0.2,-0.5,0.3]

def sign(x):
    if x == 0:
        return 0
    elif x > 0:
        return 1
    return -1


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

def quatmul(quat1, quat2):
    return[
        quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] - quat1[3] * quat2[3],
        quat1[0] * quat2[1] + quat1[1] * quat2[0] + quat1[3] * quat2[2] - quat1[2] * quat2[3],
        quat1[0] * quat2[2] + quat1[2] * quat2[0] + quat1[1] * quat2[3] - quat1[3] * quat2[1],
        quat1[0] * quat2[3] + quat1[3] * quat2[0] + quat1[2] * quat2[1] - quat1[1] * quat2[2]
    ]

def right_part(t, w):
    Msum = []
    I_Om = []
    tmp4 = []
    
    jw = mult_matrix_vec(J, w[:3])
    for i in range(3):
        Om[i] += Ms[i]/Idm
        I_Om.append(Idm * Om[i])
        tmp4.append(jw[i] + I_Om[i])

    tmp = cross_prod_vec(w, tmp4)
    
    for i in range(3):
        Msum.append(Mv[i] + Ms[i])
        
        
    tmp2 = div_vec(Msum, tmp)
    dw = mult_matrix_vec(Jinv, tmp2)
    
#    for i in range(3):
#        dw.append(w[i])            # Это если fi_tck = omega
    
#    for(int i = 0; i <3; i++) dw[i+3] = w[i];

    Lw = [0, w[0], w[1], w[2]]
    
    tmp3 = quatmul(w[3:], Lw)
    
    for i in range(len(tmp3)):
        dw.append(tmp3[i])
        
    for i in range(4):
        dw[3 + i] *= 0.5 
    
    #logging.info(dw)
    
    return dw

def RK4(t, w, h):
    count_equations = 7
    k1t = []
    k2t = []
    k3t = []
    wk = []
    
    
    k1 = right_part(t,w)
    
    for i in range(count_equations):
        k1t.append(w[i] + h * 0.5 * k1[i])
    k2 = right_part(t + h * 0.5, k1t)
    
    for i in range(count_equations):
        k2t.append(w[i] + h * 0.5 * k2[i])
    k3 = right_part(t + h * 0.5, k2t)
    
    for i in range(count_equations):
        k3t.append(w[i] + h * 0.5 * k3[i])
    k4 = right_part(t + h * 0.5, k3t)
    
    for i in range(count_equations):
        wk.append(w[i] + (h/6.0)*(k1[i] + 2.0*(k2[i] + k3[i]) + k4[i]))
        
    return wk
    

w_plot = {
    'dwx': [],
    'dwy': [],
    'dwz': [],
    'q0': [],
    'qx': [],
    'qy': [],
    'qz': [],
}

t_plot = []

w = [0.15/toRAD, 0, 0, 1, 0, 0, 0]
h = 0.01
t = 0
tk = 7200
Kfi = 20
Kw = 40
Ki = 1.0

#fi_tr = [0/toRAD, -0/toRAD, 0/toRAD]
In = [0, 0, 0]

#e = [0.7071, 0, -0.7071] #орт вектора вокруг которого хотим повернуться
#angle = 20/toRAD
#q = [cos(angle), e[0] * sin(angle), e[1] * sin(angle), e[2] * sin(angle)] # Кватернион который хотим занять


e = [0, 0, 0] #орт вектора вокруг которого хотим повернуться
angle = 20/toRAD
q = [1, 0, 0, 0]

#logging.info(w[0])


while(t < tk):
    wk = RK4(t, w, h)
    for i in range(7):
        w[i] = wk[i]
        
    
    dl_q = sqrt(w[3]**2 + w[4]**2 + w[5]**2 + w[6]**2)  #
    for i in range(len(w[3:])):                         # Нормировка кватерниона
        w[3 + i] /= dl_q                                #
    
    qs = [q[0], -q[1], -q[2], -q[3]]
    dL = quatmul(qs, w[3:])
    
    #logging.info(f't = {t} - dwx = {w[0]}  qz = {w[6]}')
    
    
    for i in range(3):
        #In[i] = In[i] +(w[i+3] - fi_tr[i])*h
        #Ms[i] = -(Kfi*(w[i+3] - fi_tr[i]) + Kw*w[i] + Ki*In[i])
        
        #In[i] = In[i] + 2 * sign(dL[0]) * dL[1 + i] * h
        #Ms[i] = -(Kfi * 2 * sign(dL[0]) * dL[1 + i] + Kw * w[i] + Ki*In[i])
        
        Ms[i] = -(Kfi * 2 * sign(dL[0]) * dL[1 + i] + Kw * w[i])
        
        for i in range(len(Ms)):
            if abs(Ms[i]) >= Ms_max:
                Ms[i] = Ms_max * sign(Ms[i])  #предельное ускорение ДМ
            if abs(Om[i]) >= Om_max:
                Ms[i] = 0
                #logging.info(f't = {t} уперлись')
    
    
    
    #print(t, w[0]*toRAD, w[1]*toRAD, w[2]*toRAD, w[3]*toRAD, w[4]*toRAD, w[5]*toRAD)
    
    t_plot.append(t)
    
    w_plot['dwx'].append(w[0])
    w_plot['dwy'].append(w[1])
    w_plot['dwz'].append(w[2])
    w_plot['q0'].append(w[3])
    w_plot['qx'].append(w[4])
    w_plot['qy'].append(w[5])
    w_plot['qz'].append(w[6])
    
    

    #//Application->ProcessMessages();
    
    t = t + h

output = f'{t_plot[-1]} '

#print(t, end=' ')
for key in w_plot.keys():
    #print(w_plot[key][-1], end=' ')
    output += f'{w_plot[key][-1]}'
    
logging.info(output)

output = f'{t_plot[0]} '

for key in w_plot.keys():
    output += f'{w_plot[key][0]}'

logging.info(output)

        
plt.plot(
    t_plot, w_plot['q0'], 'red',
    t_plot, w_plot['qx'], 'green',
    t_plot, w_plot['qy'], 'blue',
    t_plot, w_plot['qz'], 'black'
)
#plt.plot(t_plot, w_plot['dwx'])
plt.grid()
plt.show()
