#############################
#                           #
#    POINT 5 IN PROBLEM 8   #
#                           #
#############################

# borrowed to test
import numpy as np
import matplotlib.pyplot as plt

dt, t, r_ana, r_num, rel_err = np.loadtxt('data_r_err_RK4.txt', unpack=True)

total_time = 50

dt1 = dt[0]
N1 = int(total_time / dt1)
t1 = t[:N1]
r_ana1 = r_ana[:N1]
r_num1 = r_num[:N1]
rel_err1 = rel_err[:N1]

dt2 = dt[N1]
N2 = int(total_time / dt2)
t2 = t[N1:N1 + N2]
r_ana2 = r_ana[N1:N1 + N2]
r_num2 = r_num[N1:N1 + N2]
rel_err2 = rel_err[N1:N1 + N2]

dt3 = dt[N1 + N2]
N3 = int(total_time / dt3)
t3 = t[N1 + N2:N1 + N2 + N3]
r_ana3 = r_ana[N1 + N2:N1 + N2 + N3]
r_num3 = r_num[N1 + N2:N1 + N2 + N3]
rel_err3 = rel_err[N1 + N2:N1 + N2 + N3]

dt4 = dt[N1 + N2 + N3]
N4 = int(total_time / dt4)
t4 = t[N1 + N2 + N3:N1 + N2 + N3 + N4]
r_ana4 = r_ana[N1 + N2 + N3:N1 + N2 + N3 + N4]
r_num4 = r_num[N1 + N2 + N3:N1 + N2 + N3 + N4]
rel_err4 = rel_err[N1 + N2 + N3:N1 + N2 + N3 + N4]

dt5 = dt[N1 + N2 + N3 + N4]
N5 = int(total_time / dt5)
t5 = t[N1 + N2 + N3 + N4:N1 + N2 + N3 + N4 + N5]
r_ana5 = r_ana[N1 + N2 + N3 + N4:N1 + N2 + N3 + N4 + N5]
r_num5 = r_num[N1 + N2 + N3 + N4:N1 + N2 + N3 + N4 + N5]
rel_err5 = rel_err[N1 + N2 + N3 + N4:N1 + N2 + N3 + N4 + N5]

plt.figure(figsize=(8, 4.5))
# plt.title('Relative error in $\mathbf{r}_i$ using FE')
plt.title('Relative error in $\mathbf{r}_i$ using RK4')
plt.plot(t1, rel_err1, color='black', lw=1, label=f'dt = {dt1:.1e}')
plt.plot(t2, rel_err2, color='red', lw=1, label=f'dt = {dt2:.1e}')
plt.plot(t3, rel_err3, color='royalblue', lw=1, label=f'dt = {dt3:.1e}')
plt.plot(t4, rel_err4, color='green', lw=1, label=f'dt = {dt4:.1e}')
plt.xlabel('Time [$\mu s$]')
# plt.ylim([0, 1000])
plt.ylabel('$\\frac{\|\|\mathbf{r}_{ap}-\mathbf{r}_{ex}\|\|}{\|\|\mathbf{r}_{ex}\|\|}$', fontsize=16)

plt.legend(loc='upper left')
# plt.savefig('relerror_FE.pdf')
# plt.savefig('relerror_RK4.pdf')
plt.show()