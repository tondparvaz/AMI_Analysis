import numpy as np
import matplotlib.pyplot as plt
from AMI_oo_v1 import *

# > AMI analysis inputs
ell_ind = 2 # Choice of AMI length scale: 1>dlta1, 1>dlta2
dlta_per = 0.99 # Choice of BL thickness, e.g., delta99 > 99%
window = 1.5 # Windowing the field: 0>off, 1.5>windowing 1.5*delta_99

# # > Test Dataset inputs
# dat_ind = 20 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -10 # Last x-index for the AMI analysis
# x1 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# x1.initialize()
# x1.ReadData()
# x1.AMI_Budget()
# x1.Beta()
# # # x1.SkinFriction()
# # # x1.ProfPlot('asd', 2341.0)
# # > Test Dataset inputs
# dat_ind = 3 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -10 # Last x-index for the AMI analysis
# x2 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# x2.initialize()
# x2.ReadData()
# x2.AMI_Budget()
# x2.Beta()
# # plt.show()

# # > Dataset inputs
# dat_ind = 10 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -200 # Last x-index for the AMI analysis
# zpg = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# zpg.initialize()
# zpg.ReadData()
# zpg.AMI_Budget()
# # t0.Beta()
# # t0.AMI_Turb()
# # t0.Plotting2()

# # > Dataset inputs
# dat_ind = 20 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -20 # Last x-index for the AMI analysis
# w = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# w.initialize()
# w.ReadData()
# w.AMI_Budget()
# w.Plotting1()
# w.Beta()
# t6.AMI_Turb()
# t6.Plotting2()

# > Dataset inputs
dat_ind = 3 # Numerical Dataset index
xs = 1000 # First x-index for the AMI analysis
xe = -10 # Last x-index for the AMI analysis
b = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
b.initialize()
b.ReadData()
b.AMI_Budget()
# b.Plotting1()
b.AMI_Turb()

# w.ProfPlot('B-W-1/', 6260, 80, -15, 60, '')
# b.ProfPlot('B-W-1/', 6260, 80, -15, 60, '') 
# b.Beta()
# t7.AMI_Turb()
# t7.Plotting2()
# plt.show()

# # > Dataset inputs
# dat_ind = 11 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -400 # Last x-index for the AMI analysis
# b1 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# b1.initialize()
# b1.ReadData()
# b1.AMI_Budget()
# # b1.Beta()
# # t1.AMI_Turb()
# # t1.Plotting2()

# # > Dataset inputs
# dat_ind = 12 # Numerical Dataset index
# xs = 200 # First x-index for the AMI analysis
# xe = -600 # Last x-index for the AMI analysis
# m13 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# m13.initialize()
# m13.ReadData()
# m13.AMI_Budget()
# # m13.Beta()
# # t2.AMI_Turb()
# # t2.Plotting2()

# # > Dataset inputs
# dat_ind = 13 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -400 # Last x-index for the AMI analysis
# b2 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# b2.initialize()
# b2.ReadData()
# b2.AMI_Budget() 
# # b2.Beta()
# # t3.AMI_Turb()
# # t3.Plotting2()

# # > Dataset inputs
# dat_ind = 14 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -400 # Last x-index for the AMI analysis
# m16 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# m16.initialize()
# m16.ReadData()
# m16.AMI_Budget()
# # m16.Beta()
# # t4.AMI_Turb()
# # t4.Plotting2()

# # > Dataset inputs
# dat_ind = 15 # Numerical Dataset index
# xs = 0 # First x-index for the AMI analysis
# xe = -400 # Last x-index for the AMI analysis
# m18 = AMI_Analysis(dat_ind, xs, xe, ell_ind, dlta_per, window)
# m18.initialize()
# m18.ReadData()
# m18.AMI_Budget()
# m18.Beta()
# t5.AMI_Turb()
# t5.Plotting2()

plt.show()

# FileName = t1.dir + "LC.eps"
# ind = 232
# lbl = r"$x/L=$" + "{:.2f}".format(t1.xi[ind])
# # lbl = r"$x/c=$" + "{:.2f}".format(t1.xi[ind])
# plt.figure(31, figsize=(3.2,3.1))
# plt.semilogx(t1.y/t1.Lscl, t1.U[ind, :]/t1.Uscl,linewidth=0.75, color=t1.col
#         , marker='^',markevery=500, label=lbl)
# plt.semilogx(t1.y/t1.Lscl, t1.Ui[ind, :]/t1.Uscl, '--',linewidth=0.75, color=t1.col
#         , marker='^',markevery=500)
# plt.semilogx([t1.dlta[ind]/t1.Lscl, t1.dlta[ind]/t1.Lscl], [0.2, 1], 'k', linewidth=0.5, marker='^')

# ind = 358
# lbl = r"$x/L=$" + "{:.2f}".format(t1.xi[ind])
# # lbl = r"$x/c=$" + "{:.2f}".format(t1.xi[ind])
# plt.figure(31, figsize=(3.2,3.1))
# plt.semilogx(t1.y/t1.Lscl, t1.U[ind, :]/t1.Uscl,linewidth=0.75, color=t1.col
#         , marker='v',markevery=500, label=lbl)
# plt.semilogx(t1.y/t1.Lscl, t1.Ui[ind, :]/t1.Uscl, '--',linewidth=0.75, color=t1.col
#         , marker='v',markevery=500)
# plt.semilogx([t1.dlta[ind]/t1.Lscl, t1.dlta[ind]/t1.Lscl], [0.2, 1], 'k', linewidth=0.5, marker='v')

# ind = 466
# lbl = r"$x/L=$" + "{:.2f}".format(t1.xi[ind])
# # lbl = r"$x/c=$" + "{:.2f}".format(t1.xi[ind])
# plt.figure(31, figsize=(3.2,3.1))
# plt.semilogx(t1.y/t1.Lscl, t1.U[ind, :]/t1.Uscl,linewidth=0.75, color=t1.col
#         , marker='D',markevery=500, label=lbl)
# plt.semilogx(t1.y/t1.Lscl, t1.Ui[ind, :]/t1.Uscl, '--',linewidth=0.75, color=t1.col
#         , marker='D',markevery=500)
# plt.semilogx([t1.dlta[ind]/t1.Lscl, t1.dlta[ind]/t1.Lscl], [0.2, 1], 'k', linewidth=0.5, marker='D')
# plt.ylabel(r"$\overline{U}/U_\infty$", fontsize=10)
# plt.xlabel(r'$y/L$', fontsize=10)
# plt.xlim(1e-4, 1.5e-1)

# legend = plt.legend(fontsize=8)
# legend.get_frame().set_alpha(1)  # Set alpha to 1 (fully opaque)
# plt.savefig(FileName, format='eps', bbox_inches='tight')
# plt.close(31)
# plt.show()
