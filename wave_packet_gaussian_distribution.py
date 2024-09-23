# Wave packet - gaussian distribution

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from mpl_toolkits.mplot3d import proj3d

'''
def get_gaussian(line_space, sigma, mu):
    gauss = 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (line_space - mu) ** 2 / (2 * sigma ** 2))
    return gauss
'''


def get_gaussian(point, sigma, mu):
    gauss = 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (point - mu) ** 2 / (2 * sigma ** 2))
    return gauss


def update_diagrams():
    global x, k, xx, kk, yy, superposed, yy_cos_k0
    global plt_array, plt_superposed, plt_cos_k0
    x = np.linspace(range_x_min, range_x_max, 1000)
    k = np.linspace(- k0, k0, 1000)
    xx, kk = np.meshgrid(x, k)
    yy = np.cos((kk - k0) * xx) * get_gaussian(kk, sigma0, 0)

    plt_array.remove()
    plt_array = ax1.plot_wireframe(xx, kk, yy, linewidth=1, rstride=50, cstride=50)
    ax1.set_zlim(yy.min(), yy.max())
    ax1.set_ylim(- k0, k0)

    superposed = np.sum(yy, axis=0) / 1000
    plt_superposed.set_data(x, superposed)

    yy_cos_k0 = superposed.max() * np.cos(k0 * x)
    plt_cos_k0.set_data(x, yy_cos_k0)


# Setter at tkinter
def set_k0(value):
    global k0
    k0 = float(value)
    update_diagrams()


'''
def set_k1(value):
    global k1
    k1 = float(value)
    update_diagrams()
'''


def set_sigma0(value):
    global sigma0
    sigma0 = float(value)
    update_diagrams()


'''
def set_sigma1(value):
    global sigma1
    sigma1 = float(value)
    update_diagrams()
'''


def set_mu0(value):
    global mu0
    mu0 = float(value)
    update_diagrams()


'''
def set_mu1(value):
    global mu1
    mu1 = float(value)
    update_diagrams()
'''

'''
# Animation control
def step():
    global cnt_step
    pass
    cnt_step += 1


def reset():
    global is_play, cnt
    global cnt_step
    is_play = False
    # cnt = 0
    cnt_step = 1


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True
'''


def update(f):
    pass
    # global cnt
    # global txt_step
    # txt_step.set_text(str(cnt))
    # if is_play:
    #   cnt += 1
    #   pass


# Global variables
'''
# Animation control
cnt = 0
is_play = False
cnt_step = 1
'''

# Parameters
range_x_min = - 20.
range_x_max = 20.
range_y_min = -.2
range_y_max = .2
range_z_min = -.2
range_z_max = .2

range_k_min = -10
range_k_max = 10.

k0 = 5.
omega0 = 0.1
sigma0 = 0.5
mu0 = 0.

'''
k1 = 18.
omega1 = 0.1
sigma1 = 2.
mu1 = 0.
'''

# Generate figure and axes
title_ax0 = "Superposed wave"
title_ax1 = "Element waves"
title_tk = "Wave packet - gaussian distribution"
x_min0 = range_x_min
x_max0 = range_x_max
x_width0 = x_max0 - x_min0
y_min0 = range_y_min
y_max0 = range_y_max
y_width0 = y_max0 - y_min0
'''
z_min0 = range_z_min
z_max0 = range_z_max
z_width0 = z_max0 - z_min0
'''

x_min1 = range_x_min
x_max1 = range_x_max
x_width1 = x_max0 - x_min0
y_min1 = range_y_min
y_max1 = range_y_max
y_width1 = y_max1 - y_min1
z_min1 = range_z_min
z_max1 = range_z_max
z_width1 = z_max1 - z_min1

k_min = range_k_min
k_max = range_k_max

fig = Figure()
ax0 = fig.add_subplot(121)
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel('x')
ax0.set_ylabel('y')
ax0.set_xlim(x_min0, x_max0)
ax0.set_ylim(y_min0, y_max0)

ax1 = fig.add_subplot(122, projection='3d')
ax1.set_box_aspect((4, 4, 4))
ax1.grid()
ax1.set_title(title_ax1)
ax1.set_xlabel('x')
ax1.set_ylabel('delta k')
ax1.set_zlabel('y')
ax1.set_xlim(x_min1, x_max1)
ax1.set_ylim(k_min, k_max)
ax1.set_zlim(z_min1, z_max1)

# Text items
# txt_step = ax0.text(x_min0, y_max0, str(cnt_step))

# Plot items
x = np.linspace(range_x_min, range_x_max, 1000)
k = np.linspace(- k0, k0, 1000)
xx, kk = np.meshgrid(x, k)
yy = np.cos((kk - k0) * xx) * get_gaussian(kk, sigma0, 0)
plt_array = ax1.plot_wireframe(xx, kk, yy, linewidth=1, rstride=50, cstride=50)
ax1.set_zlim(yy.min(), yy.max())
ax1.set_ylim(- k0, k0)

superposed = np.sum(yy, axis=0) / 1000
plt_superposed, = ax0.plot(x, superposed, linestyle='-', label='Superposed')

yy_cos_k0 = superposed.max() * np.cos(k0 * x)
plt_cos_k0, = ax0.plot(x, yy_cos_k0, linestyle=':', linewidth=1, label='Guide (peak * cos(k0))')

# Legend
ax0.legend(loc='lower left')

# Embed in Tkinter
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Animation
'''
frm_anim = ttk.Labelframe(root, relief='ridge', text='Animation', labelanchor='n')
frm_anim.pack(side='left', fill=tk.Y)
# btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
# btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text="Step", command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
btn_reset.pack(side='left')
'''

# Parameters
frm_wave0 = ttk.Labelframe(root, relief="ridge", text="Parameter", labelanchor="n", width=100)
frm_wave0.pack(side='left')

label_k0 = tk.Label(frm_wave0, text="k0(base wave number)")
label_k0.pack(side='left')
var_k0 = tk.StringVar(root)
var_k0.set(str(k0))
s_k0 = tk.Spinbox(frm_wave0, textvariable=var_k0, format="%.2f", from_=0.1, to=20., increment=0.1,
                  command=lambda: set_k0(float(var_k0.get())), width=4)
s_k0.pack(side='left')

label_sigma0 = tk.Label(frm_wave0, text="sigma0(gaussian distribution of delta k)")
label_sigma0.pack(side='left')
var_sigma0 = tk.StringVar(root)
var_sigma0.set(str(sigma0))
s_sigma0 = tk.Spinbox(frm_wave0, textvariable=var_sigma0, format="%.2f", from_=0.1, to=20., increment=0.1,
                      command=lambda: set_sigma0(float(var_sigma0.get())), width=4)
s_sigma0.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
root.mainloop()

