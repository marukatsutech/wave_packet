# Wave packet
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk


def update_curves(t):
    global curves, z_superposed, superposed_curve
    z_superposed = x * 0.
    for i in range(num_of_waves):
        if function_number == 0:
            zc = amp[i] * np.sin(k[i] * x * np.pi - omega[i] * t)
        else:
            zc = amp[i] * np.cos(k[i] * x * np.pi - omega[i] * t)
        curves[i].set_xdata(x)
        curves[i].set_ydata(x * 0. + i)
        curves[i].set_3d_properties(zc)
        if superposition_method == 1:
            z_superposed = np.copy(z_superposed) + zc / num_of_waves
        else:
            z_superposed = np.copy(z_superposed) + zc
    superposed_curve.set_xdata(x)
    superposed_curve.set_ydata(x * 0. + y_max)
    superposed_curve.set_3d_properties(z_superposed)


def generate_amp():
    global amp
    amp.clear()
    for i in range(num_of_waves):
        if superposition_method == 2:
            amp.append(1. / (np.sqrt(2. * np.pi) * sigma) * np.exp(- (i - mu) ** 2. / (2. * sigma ** 2.)))
        else:
            amp.append(1.)


def update_gaussian():
    global gaussian, zz, amp
    gaussian.remove()
    zz = xx * 0. + 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (yy - mu) ** 2 / (2 * sigma ** 2))
    gaussian = ax1.plot_wireframe(xx, yy, zz, rcount=0, ccount=1, linewidth=1, linestyle='--')


def generate_omega():
    global omega
    omega.clear()
    if condition_of_omega == 0:
        for i in range(num_of_waves):
            omega.append(omega0)
            tx_omega.set_text("Omega=" + str(omega0))
    else:
        for i in range(num_of_waves):
            omega.append(c * k[i])
            tx_omega.set_text("Omega=c*k")


def change_sigma(sg):
    global sigma
    sigma = sg
    update_gaussian()
    generate_amp()
    update_curves(cnt)


def change_mu(mm):
    global mu
    mu = mm
    update_gaussian()
    generate_amp()
    update_curves(cnt)


def change_c(cc):
    global c
    c = cc
    generate_omega()
    update_curves(cnt)


def change_omega0(o0):
    global omega0
    omega0 = o0
    generate_omega()
    update_curves(cnt)


def generate_k():
    global k, k_round
    k.clear()
    k_round.clear()
    for i in range(num_of_waves):
        if progression_of_k == 0:
            k.append(k0 * k_diff_ratio ** i)
            k_round.append(round((k0 * k_diff_ratio ** i), 1))
        else:
            k.append(k0 + i * k_diff_ratio)
            k_round.append(round(k0 + i * k_diff_ratio, 1))
    tx_k.set_text("k=" + str(k_round))


def change_omega_condition():
    global condition_of_omega
    condition_of_omega = var_co.get()
    generate_omega()
    update_curves(cnt)


def change_progression():
    global progression_of_k, k, k0
    progression_of_k = var_p.get()
    generate_k()
    generate_omega()
    update_curves(cnt)


def change_k_diff_ratio(df):
    global k_diff_ratio
    k_diff_ratio = df
    generate_k()
    generate_omega()
    update_curves(cnt)


def change_k0(kk0):
    global k0
    k0 = kk0
    generate_k()
    update_curves(cnt)


def change_superposition():
    global superposition_method, gaussian, zz
    superposition_method = var_sp.get()
    if superposition_method == 2:
        update_gaussian()
    else:
        gaussian.remove()
        zz = xx * 0. + yy * 0. + 1.
        gaussian = ax1.plot_wireframe(xx, yy, zz, rcount=0, ccount=1, linewidth=1, linestyle='--')
    generate_amp()
    update_curves(cnt)


def change_function():
    global function_number
    function_number = var_function.get()
    if function_number == 0:
        ax1.set_title('Wave packet (sine wave)')
    else:
        ax1.set_title('Wave packet (cosine wave)')
    update_curves(cnt)


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global tx_step, cnt
    tx_step.set_text("Step=" + str(cnt))
    if is_play:
        update_curves(f)
        cnt += 1


# Global variables
num_of_points = 2000

x_min = 0.
x_max = 6.
y_min = 0.
y_max = 10.
z_min = -2.
z_max = 2.

is_play = False
cnt = 0

function_number = 0
progression_of_k = 1
condition_of_omega = 0

k0_init = 1.
k0 = k0_init
k_diff_ratio_init = 2.
k_diff_ratio = k_diff_ratio_init

omega0_init = 1.
omega0 = omega0_init
c_init = 1.
c = c_init

num_of_waves = 8
k = []
k_round = []
omega = []
amp = []

sigma_init = 0.5
sigma = sigma_init
mu_init = 4.
mu = mu_init

# Generate figure and axes
fig = Figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.set_box_aspect((4, 4, 2))
ax1.grid()
ax1.set_title('Wave packet (sine wave)')
ax1.set_xlabel('x * pi')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_zlim(z_min, z_max)

# Generate items
tx_step = ax1.text(x_min, y_max, z_max * 1.6, "Step=" + str(0))
tx_k = ax1.text(x_min, y_max, z_max * 1.3, "k=" + str(0))
tx_omega = ax1.text(x_min, y_max, z_max * 1., "Omega=" + str(omega0))

x = np.linspace(x_min, x_max, num_of_points)
y = x * 0.
z = x * 0.
z_superposed = x * 0.

superposition_method = 1

curves = []
for i_g in range(num_of_waves):
    y = x * 0. + i_g
    z = x * 0.
    curve, = ax1.plot(x, y, z, linewidth=1, linestyle='-')
    curves.append(curve)
    z_superposed = np.copy(z_superposed) + z / num_of_waves

superposed_curve, = ax1.plot(x, x * 0. + y_max, z_superposed, linestyle='-', label='Superposed')

ax1.legend(loc='upper left')

generate_k()
generate_omega()
generate_amp()
update_curves(0)

# Gaussian
xg = np.linspace(x_min, x_max, 10)
yg = np.linspace(y_min, y_max, 200)
xx, yy = np.meshgrid(xg, yg)
# zz = xx * 0. + 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (yy - mu) ** 2 / (2 * sigma ** 2))
zz = xx * 0. + yy * 0. + 1.
gaussian = ax1.plot_wireframe(xx, yy, zz, rcount=0, ccount=1, linewidth=1, linestyle='--')

# Embed in Tkinter
root = tk.Tk()
root.title("Wave packet")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

frm_f = ttk.Labelframe(root, relief="ridge", text="Function", labelanchor="n", width=100)
frm_f.pack(side='left')
var_function = tk.IntVar(value=function_number)
rdb0_function = tk.Radiobutton(frm_f, text="sine", command=change_function, variable=var_function, value=0)
rdb1_function = tk.Radiobutton(frm_f, text="cosine", command=change_function, variable=var_function, value=1)
rdb0_function.pack(anchor=tk.W)
rdb1_function.pack(anchor=tk.W)

frm_k0 = ttk.Labelframe(root, relief="ridge", text="k[0]", labelanchor="n", width=100)
frm_k0.pack(side='left')
label_k0 = tk.Label(frm_k0, text="k[0]:")
label_k0.pack(anchor=tk.W)
var_k0 = tk.StringVar(root)  # variable for spinbox-value
var_k0.set(k0_init)  # Initial value
s_k0 = tk.Spinbox(frm_k0, textvariable=var_k0, format="%.2f", from_=0.1, to=6., increment=0.1,
                  command=lambda: change_k0(float(var_k0.get())), width=4)
s_k0.pack(anchor=tk.W)
label_diff_ratio = tk.Label(frm_k0, text="Ratio or diff of k:")
label_diff_ratio.pack(anchor=tk.W)
var_diff_ratio = tk.StringVar(root)  # variable for spinbox-value
var_diff_ratio.set(k_diff_ratio_init)  # Initial value
s_diff_ratio = tk.Spinbox(frm_k0, textvariable=var_diff_ratio, format="%.2f", from_=0.1, to=6., increment=0.1,
                          command=lambda: change_k_diff_ratio(float(var_diff_ratio.get())), width=4)
s_diff_ratio.pack(anchor=tk.W)

frm_p = ttk.Labelframe(root, relief="ridge", text="Progression of k", labelanchor="n", width=100)
frm_p.pack(side='left')
var_p = tk.IntVar(value=progression_of_k)
rdb0_p = tk.Radiobutton(frm_p, text="Geometric", command=change_progression, variable=var_p, value=0)
rdb1_p = tk.Radiobutton(frm_p, text="Arithmetic", command=change_progression, variable=var_p, value=1)
rdb0_p.pack(anchor=tk.W)
rdb1_p.pack(anchor=tk.W)

frm_omega0 = ttk.Labelframe(root, relief="ridge", text="Omega", labelanchor="n", width=100)
frm_omega0.pack(side='left')
label_omega0 = tk.Label(frm_omega0, text="Omega:")
label_omega0.pack(anchor=tk.W)
var_omega0 = tk.StringVar(root)  # variable for spinbox-value
var_omega0.set(omega0_init)  # Initial value
s_omega0 = tk.Spinbox(frm_omega0, textvariable=var_omega0, format="%.2f", from_=0.1, to=6., increment=0.1,
                      command=lambda: change_omega0(float(var_omega0.get())), width=4)
s_omega0.pack(anchor=tk.W)
label_c = tk.Label(frm_omega0, text="c:")
label_c.pack(anchor=tk.W)
var_c = tk.StringVar(root)  # variable for spinbox-value
var_c.set(c_init)  # Initial value
s_c = tk.Spinbox(frm_omega0, textvariable=var_c, format="%.2f", from_=0.1, to=6., increment=0.1,
                 command=lambda: change_c(float(var_c.get())), width=4)
s_c.pack(anchor=tk.W)

frm_co = ttk.Labelframe(root, relief="ridge", text="Condition of omega", labelanchor="n", width=100)
frm_co.pack(side='left')
var_co = tk.IntVar(value=condition_of_omega)
rdb0_co = tk.Radiobutton(frm_co, text="Constant", command=change_omega_condition, variable=var_co, value=0)
rdb1_co = tk.Radiobutton(frm_co, text="Omega = c * k", command=change_omega_condition, variable=var_co, value=1)
rdb0_co.pack(anchor=tk.W)
rdb1_co.pack(anchor=tk.W)

frm_gaussian = ttk.Labelframe(root, relief="ridge", text="Gaussian", labelanchor="n", width=100)
frm_gaussian.pack(side='left')
label_sigma = tk.Label(frm_gaussian, text="Sigma:")
label_sigma.pack(anchor=tk.W)
var_sigma = tk.StringVar(root)  # variable for spinbox-value
var_sigma.set(sigma_init)  # Initial value
s_sigma = tk.Spinbox(frm_gaussian, textvariable=var_sigma, format="%.2f", from_=0.1, to=6., increment=0.1,
                     command=lambda: change_sigma(float(var_sigma.get())), width=4)
s_sigma.pack(anchor=tk.W)
label_mu = tk.Label(frm_gaussian, text="Mu:")
label_mu.pack(anchor=tk.W)
var_mu = tk.StringVar(root)  # variable for spinbox-value
var_mu.set(mu_init)  # Initial value
s_mu = tk.Spinbox(frm_gaussian, textvariable=var_mu, format="%.2f", from_=0., to=8., increment=0.1,
                  command=lambda: change_mu(float(var_mu.get())), width=4)
s_mu.pack(anchor=tk.W)

frm_sp = ttk.Labelframe(root, relief="ridge", text="Superposition", labelanchor="n", width=100)
frm_sp.pack(side='left')
var_sp = tk.IntVar(value=superposition_method)
rdb0_sp = tk.Radiobutton(frm_sp, text="Sum", command=change_superposition, variable=var_sp, value=0)
rdb1_sp = tk.Radiobutton(frm_sp, text="Average", command=change_superposition, variable=var_sp, value=1)
rdb2_sp = tk.Radiobutton(frm_sp, text="Gaussian", command=change_superposition, variable=var_sp, value=2)
rdb0_sp.pack(anchor=tk.W)
rdb1_sp.pack(anchor=tk.W)
rdb2_sp.pack(anchor=tk.W)

btn_play = tk.Button(root, text="Play/Pause", command=switch)
btn_play.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=200)
root.mainloop()
