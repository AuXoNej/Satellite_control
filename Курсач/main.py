file = open('результаты.txt')

t1_mass = []
rx_mass = []
ry_mass = []
rz_mass = []
Vx_mass = []
Vy_mass = []
Vz_mass = []
Lam1_mass = []
Lam2_mass = []
Lam3_mass = []
omx_mass = []
omy_mass = []
omz_mass = []
Omx_mass = []
Omy_mass = []
Omz_mass = []
Lam_tr_mas1 = []
Lam_tr_mas2 = []
Lam_tr_mas3 = []
Fi_mas = []
Fi_mas_tr = []

file.readline()

for line in file:
    line = line.split(',')
    t1_mass.append(float(line[0]))
    rx_mass.append(float(line[1]))
    ry_mass.append(float(line[2]))
    rz_mass.append(float(line[3]))
    Vx_mass.append(float(line[4]))
    Vy_mass.append(float(line[5]))
    Vz_mass.append(float(line[6]))
    Lam1_mass.append(float(line[7]))
    Lam2_mass.append(float(line[8]))
    Lam3_mass.append(float(line[9]))
    omx_mass.append(float(line[10]))
    omy_mass.append(float(line[11]))
    omz_mass.append(float(line[12]))
    Omx_mass.append(float(line[13]))
    Omy_mass.append(float(line[14]))
    Omz_mass.append(float(line[15]))


import matplotlib.pyplot as plt
plt.plot(t1_mass, rx_mass, t1_mass, ry_mass, t1_mass, rz_mass)
plt.ylabel('some numbers')
plt.show()

plt.plot(t1_mass, omx_mass, t1_mass, omy_mass, t1_mass, omz_mass)
plt.show()
