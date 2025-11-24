import tkinter as tk
from tkinter import ttk
import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# =========================================================
# Physical constants
# =========================================================
c = 299_792_458                      # m/s
h = 6.62607015e-34                  # J·s
J_per_eV = 1.602176634e-19

# Atomic units
AU_FREQ = 4.134137e16
AU_TIME = 2.418884e-17
AU_LENGTH = 5.291772e-11
AU_ENERGY = 4.3597447222071e-18

def safe_float(x):
    try:
        return float(x)
    except:
        return None

active_field = None
last_frequency = 3e14  # default frequency


# =========================================================
# Functions to update all fields
# =========================================================
def update_all_from_frequency(f):
    global active_field

    lam = c / f                           # meters
    lam_nm = lam * 1e9                    # nm
    T = 1 / f                             # seconds
    T_fs = T * 1e15                       # fs
    E_J = h * f
    E_eV = E_J / J_per_eV

    # SI-like units (nm, fs) with 5 significant digits
    if active_field != "freq_si": freq_si.set(f"{f:.5g}")
    if active_field != "wavelength_si": wavelength_si.set(f"{lam_nm:.5g}")
    if active_field != "period_si": period_si.set(f"{T_fs:.5g}")
    if active_field != "energy_si": energy_si.set(f"{E_J:.5g}")

    # atomic units
    if active_field != "freq_au": freq_au.set(f"{f / AU_FREQ:.5g}")
    if active_field != "wavelength_au": wavelength_au.set(f"{lam / AU_LENGTH:.5g}")
    if active_field != "period_au": period_au.set(f"{T / AU_TIME:.5g}")
    if active_field != "energy_au": energy_au.set(f"{E_J / AU_ENERGY:.5g}")

    # eV
    if active_field != "energy_eV": energy_eV.set(f"{E_eV:.5g}")

    update_plot(f)

# =========================================================
# Callbacks
# =========================================================
def cb_freq_si(event=None):
    global active_field
    active_field = "freq_si"
    f = safe_float(freq_si.get())
    if f: update_all_from_frequency(f)
    active_field = None

def cb_freq_au(event=None):
    global active_field
    active_field = "freq_au"
    x = safe_float(freq_au.get())
    if x: update_all_from_frequency(x * AU_FREQ)
    active_field = None

def cb_wavelength_si(event=None):
    global active_field
    active_field = "wavelength_si"
    x = safe_float(wavelength_si.get())
    if x: update_all_from_frequency(c / (x * 1e-9))
    active_field = None

def cb_wavelength_au(event=None):
    global active_field
    active_field = "wavelength_au"
    x = safe_float(wavelength_au.get())
    if x: update_all_from_frequency(c / (x * AU_LENGTH))
    active_field = None

def cb_period_si(event=None):
    global active_field
    active_field = "period_si"
    x = safe_float(period_si.get())
    if x: update_all_from_frequency(1 / (x * 1e-15))
    active_field = None

def cb_period_au(event=None):
    global active_field
    active_field = "period_au"
    x = safe_float(period_au.get())
    if x: update_all_from_frequency(1 / (x * AU_TIME))
    active_field = None

def cb_energy_si(event=None):
    global active_field
    active_field = "energy_si"
    x = safe_float(energy_si.get())
    if x: update_all_from_frequency(x / h)
    active_field = None

def cb_energy_au(event=None):
    global active_field
    active_field = "energy_au"
    x = safe_float(energy_au.get())
    if x: update_all_from_frequency((x * AU_ENERGY) / h)
    active_field = None

def cb_energy_eV(event=None):
    global active_field
    active_field = "energy_eV"
    x = safe_float(energy_eV.get())
    if x: update_all_from_frequency((x * J_per_eV) / h)
    active_field = None

def cb_cycles(event=None):
    update_plot(last_frequency)

def cb_phase(event=None):
    update_plot(last_frequency)

# =========================================================
# Plotting function with sin² envelope and CEP
# =========================================================
def update_plot(frequency):
    global last_frequency
    last_frequency = frequency

    # cycles
    try:
        N = float(cycles.get())
    except:
        N = 5.0

    # CEP
    try:
        phi = float(phase_pi.get()) * np.pi
    except:
        phi = 0.0

    T = 1 / frequency
    total_duration = N * T
    t = np.linspace(0, total_duration, 2000)

    envelope = np.sin(np.pi * t / (N * T)) ** 2
    envelope[(t <= 0) | (t >= total_duration)] = 0.0

    E = np.sin(2 * np.pi * frequency * t + phi) * envelope

    ax.clear()
    ax.plot(t * 1e15, E)
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Electric Field (arb. units)")
    ax.set_title(f"Laser pulse with sin² envelope, {N} cycles")
    ax.grid(True)
    fig.subplots_adjust(bottom=0.22)
    canvas.draw()

# =========================================================
# GUI setup
# =========================================================
root = tk.Tk()
root.title("Laser Calculator with sin² Envelope and Phase")
root.geometry("1000x850")
root.configure(bg="white")

# =========================================================
# GUI variables
# =========================================================
freq_si = tk.StringVar()
wavelength_si = tk.StringVar()
period_si = tk.StringVar()
energy_si = tk.StringVar()

freq_au = tk.StringVar()
wavelength_au = tk.StringVar()
period_au = tk.StringVar()
energy_au = tk.StringVar()

energy_eV = tk.StringVar()
cycles = tk.StringVar(value="5")
phase_pi = tk.StringVar(value="0")  # CEP in units of π


frame = tk.Frame(root, bg="white")
frame.pack(fill="both", expand=True, padx=10, pady=10)

# prevent grid stretching
for col in range(7):
    frame.grid_columnconfigure(col, weight=0, minsize=5)

tk.Label(frame, text="Laser Calculator", bg="white",
         font=("Arial", 20, "bold")).grid(row=0, column=0, columnspan=7, pady=10)

# =========================================================
# Helper to make a tight row
# =========================================================
def make_row(label, var_si, var_au, var_ev, row,
             cb_si, cb_au, cb_ev=None,
             unit_si="", unit_au="", unit_ev=""):
    tk.Label(frame, text=label, bg="white").grid(row=row, column=0, sticky="w", pady=3, padx=(0,5))

    # SI
    f1 = tk.Frame(frame, bg="white")
    f1.grid(row=row, column=1, sticky="w", padx=(0,1))
    e1 = ttk.Entry(f1, textvariable=var_si, width=12)
    e1.pack(side="left")
    tk.Label(f1, text=unit_si, bg="white").pack(side="left", padx=(2,0))
    e1.bind("<KeyRelease>", cb_si)

    # AU
    f2 = tk.Frame(frame, bg="white")
    f2.grid(row=row, column=2, sticky="w", padx=(0,1))
    e2 = ttk.Entry(f2, textvariable=var_au, width=12)
    e2.pack(side="left")
    tk.Label(f2, text=unit_au, bg="white").pack(side="left", padx=(2,0))
    e2.bind("<KeyRelease>", cb_au)

    # eV
    if var_ev is not None:
        f3 = tk.Frame(frame, bg="white")
        f3.grid(row=row, column=3, sticky="w", padx=(0,1))
        e3 = ttk.Entry(f3, textvariable=var_ev, width=12)
        e3.pack(side="left")
        tk.Label(f3, text=unit_ev, bg="white").pack(side="left", padx=(2,0))
        e3.bind("<KeyRelease>", cb_ev)

# =========================================================
# Create rows
# =========================================================
make_row("Frequency", freq_si, freq_au, None, 2,
         cb_freq_si, cb_freq_au, unit_si="Hz", unit_au="a.u.")

make_row("Wavelength", wavelength_si, wavelength_au, None, 3,
         cb_wavelength_si, cb_wavelength_au, unit_si="nm", unit_au="a.u.")

make_row("Period", period_si, period_au, None, 4,
         cb_period_si, cb_period_au, unit_si="fs", unit_au="a.u.")

make_row("Energy", energy_si, energy_au, energy_eV, 5,
         cb_energy_si, cb_energy_au, cb_energy_eV,
         unit_si="J", unit_au="a.u.", unit_ev="eV")

make_row("Number of cycles", cycles, tk.StringVar(), None, 6,
         cb_cycles, lambda e: None)

# CEP row
tk.Label(frame, text="Phase (π)", bg="white").grid(row=7, column=0, sticky="w", pady=3)
f_phase = tk.Frame(frame, bg="white")
f_phase.grid(row=7, column=1, sticky="w")
e_phase = ttk.Entry(f_phase, textvariable=phase_pi, width=12)
e_phase.pack(side="left")
tk.Label(f_phase, text="π", bg="white").pack(side="left", padx=(2,0))
e_phase.bind("<KeyRelease>", cb_phase)

# =========================================================
# Matplotlib plot
# =========================================================
fig = Figure(figsize=(8, 3.5), dpi=100)
ax = fig.add_subplot(111)

canvas = FigureCanvasTkAgg(fig, master=frame)
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=8, column=0, columnspan=7, pady=10)

toolbar_frame = tk.Frame(frame, bg="white")
toolbar_frame.grid(row=9, column=0, columnspan=7)
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
toolbar.update()

update_plot(last_frequency)

root.mainloop()
