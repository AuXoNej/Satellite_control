import matplotlib.pyplot as plt
import logging

from math import cos, pi, sin, sqrt


logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    level=logging.ERROR
)

CONSTANTS = {
    'TO_RAD': 180/pi,
    'J':[
        [100, 0, 0],
        [0, 200, 0],
        [0, 0, 150]
    ],
    'Jinv':[
        [1/100, 0, 0],
        [0, 1/200, 0],
        [0, 0, 1/150]
    ],
    'Idm': 0.00016,
    'MS_MAX': 0.05,
    'OM_MAX': 600,
}

TO_RAD = 180/pi

J = [
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
MS_MAX = 0.05
OM_MAX = 600

Ms = [0, 0, 0]
Mv = [0, 0, 0]  #[0.2,-0.5,0.3]

class Math_Calc():
    def sign(x):
        try:
            if x == 0:
                return 0
            elif x > 0:
                return 1
            return -1
        except Exception as error:
            logging.error(f'Ошибка выполнения sign(): {error}')
            raise error
    
    def mult_prod_vec3(vec_1, vec_2):
        """ Векторное произведение трехмерных векторов """
        
        if len(vec_1) != len(vec_2) != 3:
            logging.error(
                f'Ошибка выполнения mult_prod_vec3(): ' 
                f'векторы должны быть трехмерные'
            )
            raise ValueError(
                f'Некорректные входные данные: '
                f'векторы должны быть трехмерные'
            )
        try:
            return [
                vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1],
                vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2],
                vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0]
            ]
        except Exception as error:
            logging.error(f'Ошибка выполнения mult_prod_vec3(): {error}')
            
            raise error
    
    def mult_matrix_vec(matr, vec):
        """ Произведение матрицы на вектор """
    
        if len(matr) != len(vec):
            logging.error(
                f'Ошибка выполнения mult_matrix_vec(): '
                f'некорректные входные данные '
                f'вектор должен иметь число элементов '
                f'равное числу столбцов матрицы\n'
                f'(matr = {matr}, vec = {vec})'
            )
    
            raise ValueError(
                f'Некорректные входные данные: '
                f'вектор должен иметь число элементов '
                f'равное числу столбцов матрицы'
            )
        
        try:
            res = []
            for row in range(len(matr)):
                tmp = 0
                for col in range(len(vec)):
                    tmp += matr[row][col] * vec[col]
                res.append(tmp)
        
            return res
        except Exception as error:
            logging.error(
                f'Ошибка выполнения mult_matrix_vec(): '
                f'{error}'
            )
            raise error
    
    def div_vec(vec_1, vec_2):
        """ Разность векторов """
        
        if len(vec_1) != len(vec_2):
            logging.error(
                f'Ошибка выполнения div_vec(): '
                f'разная размерость векторов')
    
            raise ValueError(
                f'Некорректные входные данные: '
                f'разная размерность векторов'
            )
    
        try:
            vec = []
            
            for i in range(len(vec_1)):
                vec.append(vec_1[i] - vec_2[i])
            
            return vec
        except Exception as error:
            logging.error(
                f'Ошибка выполнения div_vec(): '
                f'{error}'
            )
            raise error
    
    def quatmul(quat1, quat2):
        try:
            return[
                quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] - quat1[3] * quat2[3],
                quat1[0] * quat2[1] + quat1[1] * quat2[0] + quat1[3] * quat2[2] - quat1[2] * quat2[3],
                quat1[0] * quat2[2] + quat1[2] * quat2[0] + quat1[1] * quat2[3] - quat1[3] * quat2[1],
                quat1[0] * quat2[3] + quat1[3] * quat2[0] + quat1[2] * quat2[1] - quat1[1] * quat2[2]
            ]
        except Exception as error:
            logging.error(
                f'Ошибка выполнения quatmul(): '
                f'{error}'
            )
            raise error

def right_part(t, w):
    Msum = []
    I_Om = []
    tmp4 = []
    
    jw = Math_Calc.mult_matrix_vec(CONSTANTS['J'], w[:3])

    for i in range(3):
        I_Om.append(CONSTANTS['Idm'] * Om[i])
        tmp4.append(jw[i] + I_Om[i])

    tmp = Math_Calc.mult_prod_vec3(w, tmp4)
    
    for i in range(3):
        Msum.append(Mv[i] + Ms[i])
        
        
    tmp2 = Math_Calc.div_vec(Msum, tmp)
    dw = Math_Calc.mult_matrix_vec(CONSTANTS['Jinv'], tmp2)

    Lw = [0, w[0], w[1], w[2]]
        
    dw += Math_Calc.quatmul(w[3:], Lw)
        
    for i in range(4):
        dw[3 + i] *= 0.5 
    
    for i in range(3):
        Om[i] = -Ms[i]/Idm
    
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

def main():
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
    
    w = [
        1.0/CONSTANTS['TO_RAD'], 0, 0,
        1, 0, 0, 0
    ]                   # Начальные условия (dwx, dwy, dwz, quat[4])
    h = 0.01            # Шаг
    t = 0               # Начальное время
    tk = 400            # Конечное время
    Kfi = 20
    Kw = 40
    Ki = 1.0
    
    In = [0, 0, 0]
    
    #e = [0.7071, 0, -0.7071]             # Орт вектора вокруг которого хотим повернуться
    #angle = 30/TO_RAD                    # Угол на который хотим повернуть / 2
    #q = [cos(angle), e[0]*sin(angle), e[1]*sin(angle), e[2]*sin(angle)]   # Кватернион который хотим занять
    
    e = [1, 0, 0]             
    q = [1, 0, 0, 0]
    
    logging.info(w[0])
    
    while(t < tk):
        wk = RK4(t, w, h)
        for i in range(7):
            w[i] = wk[i]
            
        
        dl_q = sqrt(w[3]**2 + w[4]**2 + w[5]**2 + w[6]**2)  #
        for i in range(len(w[3:])):                         # Нормировка кватерниона
            w[3 + i] /= dl_q                                #
        
        qs = [q[0], -q[1], -q[2], -q[3]]
        dL = Math_Calc.quatmul(qs, w[3:])
        
        #logging.info(f't = {t} - dwx = {w[0]}  qz = {w[6]}')
        
        
        for i in range(3):
            
            Ms[i] = -(Kfi * 2 * Math_Calc.sign(dL[0]) * dL[1 + i] + Kw * w[i])
            
            for i in range(len(Ms)):
                if abs(Ms[i]) >= MS_MAX:
                    Ms[i] = MS_MAX * Math_Calc.sign(Ms[i])  #предельное ускорение ДМ
                if abs(Om[i]) >= OM_MAX:
                    Ms[i] = 0
            
        t_plot.append(t)
        
        w_plot['dwx'].append(w[0])
        w_plot['dwy'].append(w[1])
        w_plot['dwz'].append(w[2])
        w_plot['q0'].append(w[3])
        w_plot['qx'].append(w[4])
        w_plot['qy'].append(w[5])
        w_plot['qz'].append(w[6])
        
        t = t + h
    
    output = f't = {t_plot[-1]} '
    for key in w_plot.keys():
        output += f'{key} = {w_plot[key][-1]} '
    
    logging.info(f'{output}')
    
    plt.plot(
        t_plot, w_plot['q0'], 'red',
        t_plot, w_plot['qx'], 'green',
        t_plot, w_plot['qy'], 'blue',
        t_plot, w_plot['qz'], 'black'
    )
    
    plt.grid()
    plt.show()


if __name__ == '__main__':
    
    main()
