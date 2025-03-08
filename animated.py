import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# Funcție helper pentru actualizarea unui Entry din slider
def update_entry(entry, value):
    entry.delete(0, tk.END)
    entry.insert(0, f"{float(value):.2f}")

# Widget personalizat pentru switch-ul de polaritate
class PolaritySwitch(tk.Canvas):
    def __init__(self, parent, initial="south", width=40, height=20, **kwargs):
        super().__init__(parent, width=width, height=height, highlightthickness=0, **kwargs)
        self.width = width
        self.height = height
        self.value = initial  # "north" sau "south"
        self.draw_switch()
        self.bind("<Button-1>", self.toggle)

    def draw_switch(self):
        self.delete("all")
        if self.value == "north":
            # Partea de sus albastră, partea de jos gri
            self.create_rectangle(0, 0, self.width, self.height/2, fill="blue", outline="")
            self.create_rectangle(0, self.height/2, self.width, self.height, fill="grey", outline="")
            self.create_text(self.width/2, self.height/4, text="NORD", fill="white", font=("Arial", 8))
            self.create_text(self.width/2, 3*self.height/4, text="SUD", fill="black", font=("Arial", 8))
        else:
            # Partea de sus gri, partea de jos roșie
            self.create_rectangle(0, 0, self.width, self.height/2, fill="grey", outline="")
            self.create_rectangle(0, self.height/2, self.width, self.height, fill="red", outline="")
            self.create_text(self.width/2, self.height/4, text="NORD", fill="black", font=("Arial", 8))
            self.create_text(self.width/2, 3*self.height/4, text="SUD", fill="white", font=("Arial", 8))

    def toggle(self, event=None):
        self.value = "north" if self.value == "south" else "south"
        self.draw_switch()
        if hasattr(self, "command") and callable(self.command):
            self.command(self.value)

# Clasa pentru fiecare magnet
class Magnet:
    def __init__(self, distance=10.0, k=2e4, phi=0.0, polarity='south'):
        self.distance = distance  # dacă este 10, magnetul este "ascuns"
        self.k = k
        self.phi = phi
        self.polarity = polarity  # "south" sau "north"
        self.show_vector = True   # Controlează vizibilitatea vectorului

# Clasa aplicației
class App:
    def __init__(self, master):
        self.master = master
        self.running = False

        # Parametrii de simulare
        self.dt = 0.05   # pasul de timp [s]
        self.mass = 0.1
        self.length = 0.5
        self.damping = 0.25
        self.g = 9.81
        self.mu0 = 4*np.pi*1e-7

        # Condiții inițiale
        self.initial_theta = 90 * np.pi/180
        self.initial_omega = 0.0
        self.theta = self.initial_theta
        self.omega = self.initial_omega

        # Variabile pentru parametrii (folosite în UI)
        self.mass_var = tk.DoubleVar(value=self.mass)
        self.length_var = tk.DoubleVar(value=self.length)
        self.damping_var = tk.DoubleVar(value=self.damping)
        self.theta_var = tk.DoubleVar(value=self.initial_theta * 180/np.pi)
        self.omega_var = tk.DoubleVar(value=self.initial_omega)
        self.dt_var = tk.DoubleVar(value=self.dt)  # slider pentru viteza animației

        # Lista de magneți (10 predefiniți, inițial inactivi)
        self.magnets = [Magnet(distance=10.0, k=2e4, phi=0.0, polarity='south') for _ in range(10)]
        self.active_magnets = 0  # 0 magneți la start

        # Variabile pentru vizibilitatea vectorilor
        self.show_grav = tk.BooleanVar(value=True)
        self.show_damp = tk.BooleanVar(value=True)
        self.show_iner = tk.BooleanVar(value=True)
        self.show_magnetic = tk.BooleanVar(value=True)

        self.setup_menu()
        self.setup_ui()
        self.setup_plots()
        self.magnet_arrows = []

    def setup_menu(self):
        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Reset", command=self.reset_simulation)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.master.quit)

    def setup_ui(self):
        style = ttk.Style(self.master)
        style.theme_use("clam")

        self.master.configure(bg="#f0f0f0")
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=0)
        self.master.rowconfigure(1, weight=1)

        # Cadru superior: controale, legendă, grafic energii
        top_frame = ttk.Frame(self.master, padding=10)
        top_frame.grid(row=0, column=0, sticky="ew")
        top_frame.columnconfigure(0, weight=1)
        top_frame.columnconfigure(1, weight=1)
        top_frame.columnconfigure(2, weight=1)

        # Panou de control (stânga)
        self.controls_frame = ttk.Labelframe(top_frame, text="Control Pendul", padding=10)
        self.controls_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        # Parametrii: M, L, damping, unghi, viteză
        self.create_param_control(self.controls_frame, "Masa [kg]:", self.mass_var, 0, 0.01, 5.0)
        self.create_param_control(self.controls_frame, "Lungime [m]:", self.length_var, 1, 0.1, 5.0)
        self.create_param_control(self.controls_frame, "Coef. amortizare:", self.damping_var, 2, 0.0, 1.0)
        self.create_param_control(self.controls_frame, "Unghi start [°]:", self.theta_var, 3, -180, 180)
        self.create_param_control(self.controls_frame, "Viteza inițială [rad/s]:", self.omega_var, 4, -10, 10)

        # Parametru suplimentar: dt (viteza animației)
        self.create_param_control(self.controls_frame, "Viteză animație (dt):", self.dt_var, 5, 0.01, 0.2)

        # Butoane Start, Pauză, Stop, Save
        btn_frame = ttk.Frame(self.controls_frame)
        btn_frame.grid(row=6, column=0, columnspan=3, pady=5)
        self.start_button = ttk.Button(btn_frame, text="Start", command=self.start_sim)
        self.start_button.grid(row=0, column=0, padx=2)
        self.pause_button = ttk.Button(btn_frame, text="Pauză", command=self.pause_sim)
        self.pause_button.grid(row=0, column=1, padx=2)
        self.stop_button = ttk.Button(btn_frame, text="Stop", command=self.stop_sim)
        self.stop_button.grid(row=0, column=2, padx=2)
        self.save_button = ttk.Button(btn_frame, text="Save Snapshot", command=self.save_snapshot)
        self.save_button.grid(row=0, column=3, padx=2)

        # Controale pentru magneți
        ttk.Label(self.controls_frame, text="Magneți:").grid(row=7, column=0, columnspan=3, pady=(10, 0))
        magnet_buttons = ttk.Frame(self.controls_frame)
        magnet_buttons.grid(row=8, column=0, columnspan=3, pady=2)
        self.mag_plus_button = ttk.Button(magnet_buttons, text="+", command=self.add_magnet)
        self.mag_plus_button.grid(row=0, column=0, padx=2)
        self.mag_minus_button = ttk.Button(magnet_buttons, text="-", command=self.remove_magnet)
        self.mag_minus_button.grid(row=0, column=1, padx=2)

        # Parametrii individuali ai magneților
        self.magnet_params_frame = ttk.Frame(self.controls_frame)
        self.magnet_params_frame.grid(row=9, column=0, columnspan=3, sticky="ew", pady=5)
        self.update_magnet_controls()

        # Panou legendă (centru)
        self.legend_frame = ttk.Labelframe(top_frame, text="Legenda vectorilor", padding=10)
        self.legend_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        # Forța gravitațională
        grav_canvas = tk.Canvas(self.legend_frame, width=20, height=20)
        grav_canvas.create_rectangle(0, 0, 20, 20, fill="yellow", outline="")
        grav_canvas.grid(row=0, column=0, padx=2)
        ttk.Checkbutton(self.legend_frame, text="Forța gravitațională", variable=self.show_grav,
                        command=self.update_display).grid(row=0, column=1, sticky="w")
        # Forța de amortizare
        damp_canvas = tk.Canvas(self.legend_frame, width=20, height=20)
        damp_canvas.create_rectangle(0, 0, 20, 20, fill="green", outline="")
        damp_canvas.grid(row=1, column=0, padx=2)
        ttk.Checkbutton(self.legend_frame, text="Forța de amortizare", variable=self.show_damp,
                        command=self.update_display).grid(row=1, column=1, sticky="w")
        # Forța inerțială
        iner_canvas = tk.Canvas(self.legend_frame, width=20, height=20)
        iner_canvas.create_rectangle(0, 0, 20, 20, fill="magenta", outline="")
        iner_canvas.grid(row=2, column=0, padx=2)
        ttk.Checkbutton(self.legend_frame, text="Forța inerțială", variable=self.show_iner,
                        command=self.update_display).grid(row=2, column=1, sticky="w")
        # Forțele magnetice
        mag_canvas = tk.Canvas(self.legend_frame, width=20, height=20)
        mag_canvas.create_rectangle(0, 0, 10, 20, fill="red", outline="")
        mag_canvas.create_rectangle(10, 0, 20, 20, fill="blue", outline="")
        mag_canvas.grid(row=3, column=0, padx=2)
        ttk.Checkbutton(self.legend_frame, text="Forțele magnetice", variable=self.show_magnetic,
                        command=self.update_display).grid(row=3, column=1, sticky="w")

        # Panou pentru graficul energiilor (dreapta)
        self.energy_frame = ttk.Labelframe(top_frame, text="Grafic Energetic", padding=10)
        self.energy_frame.grid(row=0, column=2, sticky="nsew", padx=5, pady=5)

        # Panou pentru animație (jos)
        self.anim_frame = ttk.Frame(self.master, padding=10)
        self.anim_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def create_param_control(self, parent, label_text, var, row, min_val, max_val):
        """Creează un control (label, entry, slider cu min și max) pentru un parametru."""
        frame = ttk.Frame(parent)
        frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=2)
        ttk.Label(frame, text=label_text).grid(row=0, column=0, sticky="w")
        entry = ttk.Entry(frame, textvariable=var, width=8)
        entry.grid(row=0, column=1, sticky="w")
        slider_frame = ttk.Frame(frame)
        slider_frame.grid(row=0, column=2, sticky="ew", padx=5)
        # Eticheta stânga (min)
        min_label = ttk.Label(slider_frame, text=str(min_val))
        min_label.grid(row=0, column=0, sticky="w")
        # Slider
        slider = ttk.Scale(slider_frame, from_=min_val, to=max_val, orient="horizontal", variable=var,
                           command=lambda v, e=entry: update_entry(e, v))
        slider.grid(row=0, column=1, sticky="ew")
        slider_frame.columnconfigure(1, weight=1)
        # Eticheta dreapta (max)
        max_label = ttk.Label(slider_frame, text=str(max_val))
        max_label.grid(row=0, column=2, sticky="e")

    def update_magnet_controls(self):
        for widget in self.magnet_params_frame.winfo_children():
            widget.destroy()
        self.magnet_controls = []
        for i in range(self.active_magnets):
            frame = ttk.Frame(self.magnet_params_frame, relief="sunken", padding=5)
            frame.pack(fill=tk.X, pady=2)
            ttk.Label(frame, text=f"Magnet {i+1}").grid(row=0, column=0, columnspan=4)

            # Dist
            dist_var = tk.DoubleVar(value=self.magnets[i].distance)
            ttk.Label(frame, text="Dist [m]:").grid(row=1, column=0, sticky="w")
            dist_entry = ttk.Entry(frame, textvariable=dist_var, width=6)
            dist_entry.grid(row=1, column=1, sticky="w")
            slider_frame = ttk.Frame(frame)
            slider_frame.grid(row=1, column=2, sticky="ew", padx=2)
            ttk.Label(slider_frame, text="0.1").grid(row=0, column=0, sticky="w")
            dist_slider = ttk.Scale(slider_frame, from_=0.1, to=5.0, orient="horizontal", variable=dist_var,
                                    command=lambda v, e=dist_entry: update_entry(e, v))
            dist_slider.grid(row=0, column=1, sticky="ew")
            slider_frame.columnconfigure(1, weight=1)
            ttk.Label(slider_frame, text="5.0").grid(row=0, column=2, sticky="e")

            # k
            k_var = tk.DoubleVar(value=self.magnets[i].k)
            ttk.Label(frame, text="K:").grid(row=2, column=0, sticky="w")
            k_entry = ttk.Entry(frame, textvariable=k_var, width=6)
            k_entry.grid(row=2, column=1, sticky="w")
            slider_frame2 = ttk.Frame(frame)
            slider_frame2.grid(row=2, column=2, sticky="ew", padx=2)
            ttk.Label(slider_frame2, text="1e3").grid(row=0, column=0, sticky="w")
            k_slider = ttk.Scale(slider_frame2, from_=1e3, to=1e5, orient="horizontal", variable=k_var,
                                 command=lambda v, e=k_entry: update_entry(e, v))
            k_slider.grid(row=0, column=1, sticky="ew")
            slider_frame2.columnconfigure(1, weight=1)
            ttk.Label(slider_frame2, text="1e5").grid(row=0, column=2, sticky="e")

            # phi
            phi_var = tk.DoubleVar(value=self.magnets[i].phi)
            ttk.Label(frame, text="Phi [°]:").grid(row=3, column=0, sticky="w")
            phi_entry = ttk.Entry(frame, textvariable=phi_var, width=6)
            phi_entry.grid(row=3, column=1, sticky="w")
            slider_frame3 = ttk.Frame(frame)
            slider_frame3.grid(row=3, column=2, sticky="ew", padx=2)
            ttk.Label(slider_frame3, text="-180").grid(row=0, column=0, sticky="w")
            phi_slider = ttk.Scale(slider_frame3, from_=-180, to=180, orient="horizontal", variable=phi_var,
                                   command=lambda v, e=phi_entry: update_entry(e, v))
            phi_slider.grid(row=0, column=1, sticky="ew")
            slider_frame3.columnconfigure(1, weight=1)
            ttk.Label(slider_frame3, text="180").grid(row=0, column=2, sticky="e")

            # Switch polaritate + checkbutton vector
            row_switch = 4
            pol_label = ttk.Label(frame, text="Polaritate:")
            pol_label.grid(row=row_switch, column=0, sticky="w")
            pol_switch = PolaritySwitch(frame, initial=self.magnets[i].polarity)
            pol_switch.grid(row=row_switch, column=1, padx=5)
            pol_switch.command = lambda val, idx=i: self.set_magnet_polarity(idx, val)

            show_vec = tk.BooleanVar(value=self.magnets[i].show_vector)
            vector_cb = ttk.Checkbutton(frame, text="Afișează vector", variable=show_vec,
                                        command=lambda var=show_vec, idx=i: self.set_magnet_vector(idx, var.get()))
            vector_cb.grid(row=row_switch, column=2, sticky="w")

            self.magnet_controls.append({
                'dist_var': dist_var,
                'k_var': k_var,
                'phi_var': phi_var,
                'polarity_switch': pol_switch,
                'show_vec': show_vec
            })

    def set_magnet_polarity(self, idx, value):
        self.magnets[idx].polarity = value
        self.update_display()

    def set_magnet_vector(self, idx, show):
        self.magnets[idx].show_vector = show
        self.update_display()

    def add_magnet(self):
        if self.active_magnets < len(self.magnets):
            self.magnets[self.active_magnets].distance = 0.75
            self.magnets[self.active_magnets].k = 2e4
            self.magnets[self.active_magnets].phi = 0.0
            self.magnets[self.active_magnets].polarity = 'south'
            self.magnets[self.active_magnets].show_vector = True
            self.active_magnets += 1
            self.update_magnet_controls()
            self.update_display()

    def remove_magnet(self):
        if self.active_magnets > 0:
            self.active_magnets -= 1
            self.magnets[self.active_magnets].distance = 10.0
            self.update_magnet_controls()
            self.update_display()

    def setup_plots(self):
        self.fig_anim, self.ax_anim = plt.subplots(figsize=(5, 5))
        self.ax_anim.set_aspect('equal')
        self.ax_anim.grid(True)
        self.ax_anim.set_xlim([-1.5, 1.5])
        self.ax_anim.set_ylim([-1.5, 1.5])
        self.ax_anim.set_xlabel('x [m]')
        self.ax_anim.set_ylabel('y [m]')
        self.ax_anim.set_title('Pendul')
        self.canvas_anim = FigureCanvasTkAgg(self.fig_anim, master=self.anim_frame)
        self.canvas_anim.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.fig_energy, self.ax_energy = plt.subplots(figsize=(4, 3))
        self.ax_energy.set_xlabel('Unghi [°]')
        self.ax_energy.set_ylabel('Energie [J]')
        self.ax_energy.set_title('Energii')
        self.ax_energy.grid(True)
        self.canvas_energy = FigureCanvasTkAgg(self.fig_energy, master=self.energy_frame)
        self.canvas_energy.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        x_bob = self.length * np.sin(self.theta)
        y_bob = -self.length * np.cos(self.theta)
        self.pendulum_line, = self.ax_anim.plot([0, x_bob], [0, y_bob], 'b-', linewidth=2)
        self.bob_marker, = self.ax_anim.plot(x_bob, y_bob, 'bo', markersize=8)
        self.quiv_grav = self.ax_anim.quiver(x_bob, y_bob, 0, 0, color='yellow', angles='xy', scale_units='xy', scale=1)
        self.quiv_damp = self.ax_anim.quiver(x_bob, y_bob, 0, 0, color='green', angles='xy', scale_units='xy', scale=1)
        self.quiv_iner = self.ax_anim.quiver(x_bob, y_bob, 0, 0, color='magenta', angles='xy', scale_units='xy', scale=1)
        self.magnet_arrows = []
        for i in range(self.active_magnets):
            mag = self.magnets[i]
            color = 'red' if mag.polarity=='south' else 'blue'
            q = self.ax_anim.quiver(x_bob, y_bob, 0, 0, color=color, angles='xy', scale_units='xy', scale=1)
            self.magnet_arrows.append(q)
        self.pivot_marker, = self.ax_anim.plot(0, 0, 'ko', markersize=10)
        self.magnet_markers = []
        for i in range(len(self.magnets)):
            phi_rad = np.deg2rad(self.magnets[i].phi)
            x_m = self.magnets[i].distance * np.sin(phi_rad)
            y_m = -self.magnets[i].distance * np.cos(phi_rad)
            marker, = self.ax_anim.plot(x_m, y_m, 'ro', markersize=12)
            marker.set_visible(i < self.active_magnets)
            marker.set_color('red' if self.magnets[i].polarity=='south' else 'blue')
            self.magnet_markers.append(marker)

        self.theta_array = np.linspace(-np.pi/2, np.pi/2, 200)
        self.energy_line_total, = self.ax_energy.plot([], [], 'k-', linewidth=2, label='E_total')
        self.energy_line_grav, = self.ax_energy.plot([], [], 'b--', linewidth=1.5, label='E_grav')
        self.energy_line_mag, = self.ax_energy.plot([], [], 'r--', linewidth=1.5, label='E_mag')
        self.energy_marker, = self.ax_energy.plot([], [], 'mo', markersize=6)
        self.ax_energy.legend()

    def compute_energy(self, theta):
        m = self.mass
        g = self.g
        l = self.length
        U_grav = -m*g*l*np.cos(theta)
        U_mag = 0
        for i in range(self.active_magnets):
            mag = self.magnets[i]
            d = mag.distance
            k = mag.k
            phi = np.deg2rad(mag.phi)
            r = np.sqrt(l**2 + d**2 - 2*l*d*np.cos(theta - phi))
            U_i = - (self.mu0/(4*np.pi))*(k/(r**3))*(3*np.cos(theta-phi)**2 - 1)
            if mag.polarity=='north':
                U_i *= -1
            U_mag += U_i
        return U_grav, U_mag, U_grav + U_mag

    def compute_dU_dtheta(self, theta):
        eps = 1e-5
        _, _, U_plus = self.compute_energy(theta+eps)
        _, _, U_minus = self.compute_energy(theta-eps)
        return (U_plus - U_minus)/(2*eps)

    def update_display(self):
        # Actualizăm valorile parametri
        self.mass = float(self.mass_var.get())
        self.length = float(self.length_var.get())
        self.damping = float(self.damping_var.get())
        self.dt = float(self.dt_var.get())

        x_bob = self.length*np.sin(self.theta)
        y_bob = -self.length*np.cos(self.theta)
        self.pendulum_line.set_data([0,x_bob],[0,y_bob])
        self.bob_marker.set_data(x_bob,y_bob)

        # Forțe
        scale = 0.5
        F_grav = np.array([0, -self.mass*self.g])
        if self.show_grav.get():
            self.quiv_grav.set_visible(True)
            self.quiv_grav.set_offsets([x_bob, y_bob])
            self.quiv_grav.set_UVC(F_grav[0]*scale, F_grav[1]*scale)
        else:
            self.quiv_grav.set_visible(False)

        tangent = np.array([np.cos(self.theta), np.sin(self.theta)])
        F_damp = - self.damping*(self.length*self.omega)*tangent
        if self.show_damp.get():
            self.quiv_damp.set_visible(True)
            self.quiv_damp.set_offsets([x_bob, y_bob])
            self.quiv_damp.set_UVC(F_damp[0]*scale, F_damp[1]*scale)
        else:
            self.quiv_damp.set_visible(False)

        def f(state):
            theta, omega = state
            dU_dtheta = self.compute_dU_dtheta(theta)
            return np.array([omega, - dU_dtheta/(self.mass*self.length**2) - self.damping*omega])
        state = np.array([self.theta, self.omega])
        angular_acc = f(state)[1]
        F_iner = self.mass*self.length*angular_acc*tangent
        if self.show_iner.get():
            self.quiv_iner.set_visible(True)
            self.quiv_iner.set_offsets([x_bob, y_bob])
            self.quiv_iner.set_UVC(F_iner[0]*scale, F_iner[1]*scale)
        else:
            self.quiv_iner.set_visible(False)

        # Magneți
        for i in range(len(self.magnets)):
            mag = self.magnets[i]
            if i < self.active_magnets:
                phi_rad = np.deg2rad(mag.phi)
                x_m = mag.distance*np.sin(phi_rad)
                y_m = -mag.distance*np.cos(phi_rad)
                self.magnet_markers[i].set_data(x_m,y_m)
                self.magnet_markers[i].set_visible(True)
                self.magnet_markers[i].set_color('red' if mag.polarity=='south' else 'blue')
            else:
                self.magnet_markers[i].set_visible(False)

            if i < self.active_magnets:
                # Citim valorile actuale din slider/entry
                try:
                    mag.distance = float(self.magnet_controls[i]['dist_var'].get())
                    mag.k = float(self.magnet_controls[i]['k_var'].get())
                    mag.phi = float(self.magnet_controls[i]['phi_var'].get())
                except:
                    pass
                phi_rad = np.deg2rad(mag.phi)
                x_m = mag.distance*np.sin(phi_rad)
                y_m = -mag.distance*np.cos(phi_rad)
                r = np.sqrt(self.length**2 + mag.distance**2 - 2*self.length*mag.distance*np.cos(self.theta - phi_rad))
                F_mag = - (self.mu0/(4*np.pi))*(mag.k/(r**4))*(3*np.cos(self.theta-phi_rad)**2 - 1)
                bob_pos = np.array([x_bob,y_bob])
                magnet_pos = np.array([x_m,y_m])
                dir_vec = magnet_pos - bob_pos
                norm_dir = np.linalg.norm(dir_vec)
                unit_dir = dir_vec/norm_dir if norm_dir!=0 else np.array([0,0])
                if mag.polarity=='north':
                    arrow = - unit_dir*abs(F_mag)
                else:
                    arrow = unit_dir*abs(F_mag)
                if i < len(self.magnet_arrows):
                    q = self.magnet_arrows[i]
                    q.set_color('red' if mag.polarity=='south' else 'blue')
                    q.set_offsets(bob_pos)
                    q.set_UVC(arrow[0]*scale, arrow[1]*scale)
                    q.set_visible(self.show_magnetic.get() and mag.show_vector)
                else:
                    q = self.ax_anim.quiver(x_bob, y_bob, arrow[0]*scale, arrow[1]*scale,
                                            color='red' if mag.polarity=='south' else 'blue',
                                            angles='xy', scale_units='xy', scale=1)
                    q.set_visible(self.show_magnetic.get() and mag.show_vector)
                    self.magnet_arrows.append(q)

        self.canvas_anim.draw()

        # Grafic energii
        theta_arr = np.linspace(-np.pi/2, np.pi/2, 200)
        U_grav_arr = - self.mass*self.g*self.length*np.cos(theta_arr)
        U_mag_arr = np.zeros_like(theta_arr)
        for i in range(self.active_magnets):
            mag = self.magnets[i]
            phi_rad = np.deg2rad(mag.phi)
            r_arr = np.sqrt(self.length**2 + mag.distance**2 - 2*self.length*mag.distance*np.cos(theta_arr - phi_rad))
            U_i_arr = - (self.mu0/(4*np.pi))*(mag.k/(r_arr**3))*(3*np.cos(theta_arr - phi_rad)**2 - 1)
            if mag.polarity=='north':
                U_i_arr *= -1
            U_mag_arr += U_i_arr
        U_total_arr = U_grav_arr + U_mag_arr
        self.energy_line_total.set_data(theta_arr*180/np.pi, U_total_arr)
        self.energy_line_grav.set_data(theta_arr*180/np.pi, U_grav_arr)
        self.energy_line_mag.set_data(theta_arr*180/np.pi, U_mag_arr)
        _, _, U_total_curr = self.compute_energy(self.theta)
        self.energy_marker.set_data([self.theta*180/np.pi],[U_total_curr])
        self.ax_energy.relim()
        self.ax_energy.autoscale_view()
        self.canvas_energy.draw()

    def sim_step(self):
        if not self.running:
            return
        # Citim dt la fiecare pas, pentru a putea fi modificat din slider
        self.dt = float(self.dt_var.get())

        def f(state):
            theta, omega = state
            dU_dtheta = self.compute_dU_dtheta(theta)
            return np.array([omega, - dU_dtheta/(self.mass*self.length**2) - self.damping*omega])
        state = np.array([self.theta,self.omega])
        k1 = f(state)
        k2 = f(state + 0.5*self.dt*k1)
        k3 = f(state + 0.5*self.dt*k2)
        k4 = f(state + self.dt*k3)
        state_new = state + (self.dt/6)*(k1+2*k2+2*k3+k4)
        self.theta, self.omega = state_new
        self.update_display()
        self.master.after(int(self.dt*1000), self.sim_step)

    def start_sim(self):
        if not self.running:
            self.running = True
            self.sim_step()

    def pause_sim(self):
        self.running = False

    def stop_sim(self):
        self.running = False
        try:
            self.theta = float(self.theta_var.get())*np.pi/180
            self.omega = float(self.omega_var.get())
        except:
            pass
        self.update_display()

    def save_snapshot(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                 filetypes=[("PNG files", "*.png"),("All Files","*.*")],
                                                 title="Save Snapshot")
        if file_path:
            self.fig_anim.savefig(file_path)
            messagebox.showinfo("Save Snapshot", f"Snapshot saved to {file_path}")

    def reset_simulation(self):
        self.theta = self.initial_theta
        self.omega = self.initial_omega
        self.mass = float(self.mass_var.get())
        self.length = float(self.length_var.get())
        self.damping = float(self.damping_var.get())
        self.dt = float(self.dt_var.get())
        self.update_display()

if __name__ == '__main__':
    root = tk.Tk()
    root.title("Pendul cu Magneți")
    app = App(root)
    root.mainloop()