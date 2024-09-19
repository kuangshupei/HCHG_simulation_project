import pynlo
import threading
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
from utils import *
from PIL import Image, ImageTk
from tkinter import ttk, messagebox
from gas_property_data import noble_gases, loadGasParameters
from simulation import GasPropertyBuilder, PeakIntensityBuilder

class SimulationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Highly Cascaded Harmonic Generation (HCHG) Simulation Software")
        
        # set the default window size
        window_width = 1600
        window_height = 900
        center_window(self.root, window_width, window_height)
        # track the progress window
        self.progress_window = None

        self.create_widgets()
        self.add_scroll_binding()

        # override the close window protocol
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def create_widgets(self):

        # Create a canvas and scrollbars for the entire window
        self.canvas = tk.Canvas(self.root)
        self.v_scrollbar = tk.Scrollbar(self.root, orient="vertical", command=self.canvas.yview)
        self.h_scrollbar = tk.Scrollbar(self.root, orient="horizontal", command=self.canvas.xview)
        self.scrollable_frame = tk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(
                scrollregion=self.canvas.bbox("all")
            )
        )

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set, xscrollcommand=self.h_scrollbar.set)

        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.v_scrollbar.grid(row=0, column=1, sticky="ns")
        self.h_scrollbar.grid(row=1, column=0, sticky="ew")

        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        self.entries = {}

        # create a frame for all input sections
        input_frame = tk.Frame(self.scrollable_frame)
        input_frame.grid(row=0, column=0, padx=10, pady=5, sticky="nsew")

        # laser properties section
        laser_frame = tk.LabelFrame(input_frame, text="Laser Properties", padx=10, pady=10)
        laser_frame.grid(row=0, column=0, padx=10, pady=5, sticky="ew")
        
        self.add_entry(laser_frame, 'Pulse Wavelength [nm]', 0)
        self.add_entry(laser_frame, 'FWHM Duration [ps]', 1)
        self.add_entry(laser_frame, 'Repetition Rate [MHz]', 2)
        self.add_entry(laser_frame, 'Energy Per Pulse [mJ]', 3)
        self.add_entry(laser_frame, 'Pulse Shape', 4)
        self.add_entry(laser_frame, 'Group Delay Dispersion [ps^2]', 5)
        self.add_entry(laser_frame, 'Third Order Dispersion [ps^3]', 6)

        # fiber properties section
        fiber_frame = tk.LabelFrame(input_frame, text="Fiber Properties", padx=10, pady=10)
        fiber_frame.grid(row=1, column=0, padx=10, pady=5, sticky="ew")
        
        self.add_entry(fiber_frame, 'Fiber Length [mm]', 0)
        self.add_entry(fiber_frame, 'Fiber Radius [microns]', 1)
        self.add_entry(fiber_frame, 'Fiber Wavelength [nm]', 2)
        self.add_entry(fiber_frame, 'Attentuation Coefficient [dB/cm]', 3)
        
        # gas properties section
        gas_frame = tk.LabelFrame(input_frame, text="Gas Properties", padx=10, pady=10)
        gas_frame.grid(row=2, column=0, padx=10, pady=5, sticky="ew")
        
        self.add_entry(gas_frame, 'Gas Name', 0)

        self.differential_pumping_var = tk.BooleanVar()
        self.differential_pumping_checkbox = tk.Checkbutton(gas_frame, text="Differential Pumping", variable=self.differential_pumping_var, command=self.toggle_pressure_inputs)
        self.differential_pumping_checkbox.grid(row=1, column=0, columnspan=2, pady=10)

        self.pressure_unit_var = tk.StringVar(value="Bar")
        self.pressure_unit_label = tk.Label(gas_frame, text="Pressure Units")
        self.pressure_unit_label.grid(row=2, column=0, padx=10, pady=5)

        self.pressure_unit_menu = ttk.Combobox(gas_frame, textvariable=self.pressure_unit_var, values=["Bar", "Torr"], state="readonly")
        self.pressure_unit_menu.grid(row=2, column=1, padx=10, pady=5)

        self.constant_pressure_label = tk.Label(gas_frame, text="Constant Pressure")
        self.constant_pressure_label.grid(row=3, column=0, padx=10, pady=5)
        self.constant_pressure_entry = tk.Entry(gas_frame)
        self.constant_pressure_entry.grid(row=3, column=1, padx=10, pady=5)

        self.pressure_start_label = tk.Label(gas_frame, text="Pressure At Entrance")
        self.pressure_start_label.grid(row=4, column=0, padx=10, pady=5)
        self.pressure_boundaries_start = tk.Entry(gas_frame)
        self.pressure_boundaries_start.grid(row=4, column=1, padx=10, pady=5)

        self.pressure_end_label = tk.Label(gas_frame, text="Pressure At Exit")
        self.pressure_end_label.grid(row=5, column=0, padx=10, pady=5)
        self.pressure_boundaries_end = tk.Entry(gas_frame)
        self.pressure_boundaries_end.grid(row=5, column=1, padx=10, pady=5)

        self.entries['Constant Pressure'] = self.constant_pressure_entry
        self.entries['Pressure At Entrance'] = self.pressure_boundaries_start
        self.entries['Pressure At Exit'] = self.pressure_boundaries_end

        self.toggle_pressure_inputs()

        # nonlinear optics section
        NLO_frame = tk.LabelFrame(input_frame, text="Nonlinear Optics", padx=10, pady=10)
        NLO_frame.grid(row=3, column=0, padx=10, pady=5, sticky="ew")
        
        self.add_entry(NLO_frame, 'SHG Conversion Efficiency', 0)

        self.raman_var = tk.BooleanVar()
        self.raman_checkbox = tk.Checkbutton(NLO_frame, text="Enable Raman Effect", variable=self.raman_var)
        self.raman_checkbox.grid(row=1, column=0, columnspan=2, pady=10)

        self.steep_var = tk.BooleanVar()
        self.steep_checkbox = tk.Checkbutton(NLO_frame, text="Enable Self-Steepening", variable=self.steep_var)
        self.steep_checkbox.grid(row=2, column=0, columnspan=4, pady=10)

        # allocate EPP manually section
        epp_frame = tk.LabelFrame(input_frame, text="Manual EPP Allocation", padx=10, pady=10)
        epp_frame.grid(row=4, column=0, padx=10, pady=5, sticky="ew")

        self.allocate_epp_var = tk.BooleanVar()
        self.allocate_epp_checkbox = tk.Checkbutton(epp_frame, text="Allocate EPP Manually", variable=self.allocate_epp_var, command=self.toggle_epp_inputs)
        self.allocate_epp_checkbox.grid(row=0, column=0, columnspan=2, pady=10)

        self.fundamental_pulse_energy_label = tk.Label(epp_frame, text="Fundamental Pulse Energy [μJ]")
        self.fundamental_pulse_energy_label.grid(row=1, column=0, padx=10, pady=5)
        self.fundamental_pulse_energy_entry = tk.Entry(epp_frame)
        self.fundamental_pulse_energy_entry.grid(row=1, column=1, padx=10, pady=5)

        self.second_harmonic_pulse_energy_label = tk.Label(epp_frame, text="Second Harmonic Pulse Energy [μJ]")
        self.second_harmonic_pulse_energy_label.grid(row=2, column=0, padx=10, pady=5)
        self.second_harmonic_pulse_energy_entry = tk.Entry(epp_frame)
        self.second_harmonic_pulse_energy_entry.grid(row=2, column=1, padx=10, pady=5)

        self.toggle_epp_inputs()

        # simulation options section
        simulation_frame = tk.LabelFrame(input_frame, text="Simulation Options", padx=10, pady=10)
        simulation_frame.grid(row=5, column=0, padx=10, pady=5, sticky="ew")
        
        self.add_entry(simulation_frame, 'Simulation Steps', 0)
        self.add_entry(simulation_frame, 'Simulation Window [ps]', 1)
        self.add_entry(simulation_frame, 'Simulation Points', 2)

        self.reload_fiber_each_step_var = tk.BooleanVar()
        self.reload_fiber_each_step_checkbox = tk.Checkbutton(simulation_frame, text="Reload Fiber Each Step", variable=self.reload_fiber_each_step_var)
        self.reload_fiber_each_step_checkbox.grid(row=3, column=0, columnspan=2, pady=10)

        self.run_button = tk.Button(simulation_frame, text="Run Simulation", command=self.run_simulation)
        self.run_button.grid(row=4, column=0, columnspan=2, pady=10)
    
        # plotting section
        plots_frame = tk.Frame(self.scrollable_frame, width=600)
        plots_frame.grid(row=0, column=1, rowspan=6, columnspan=2, padx=10, pady=5, sticky="nsew")

        self.tab_control = ttk.Notebook(plots_frame)
        self.tab1 = ttk.Frame(self.tab_control)
        self.tab2 = ttk.Frame(self.tab_control)
        self.tab3 = ttk.Frame(self.tab_control)
        self.tab4 = ttk.Frame(self.tab_control)
        self.tab5 = ttk.Frame(self.tab_control)
        self.tab_control.add(self.tab1, text='Plot 1')
        self.tab_control.add(self.tab2, text='Plot 2')
        self.tab_control.add(self.tab3, text='Plot 3')
        self.tab_control.add(self.tab4, text='Plot 4')
        self.tab_control.add(self.tab5, text='Combined Plot')
        self.tab_control.pack(expand=1, fill="both")

    def add_scroll_binding(self):
        # Bind mouse wheel events for Windows, Mac, and Linux
        self.canvas.bind_all("<MouseWheel>", self.on_mouse_wheel)  # Windows and Linux
        self.canvas.bind_all("<Button-4>", self.on_mouse_wheel)    # Linux
        self.canvas.bind_all("<Button-5>", self.on_mouse_wheel)    # Linux
        self.canvas.bind_all("<Shift-MouseWheel>", self.on_shift_mouse_wheel)  # Windows and Linux
        self.canvas.bind_all("<Shift-Button-4>", self.on_shift_mouse_wheel)    # Linux
        self.canvas.bind_all("<Shift-Button-5>", self.on_shift_mouse_wheel)    # Linux

        # Special bindings for Mac touchpad scrolling
        self.root.bind_all("<Command-MouseWheel>", self.on_mouse_wheel)  # Mac
        self.root.bind_all("<Command-Shift-MouseWheel>", self.on_shift_mouse_wheel)  # Mac

    def on_mouse_wheel(self, event):
        # For Mac touchpad, use event.delta directly
        if event.num == 4 or event.delta > 0:
            self.canvas.yview_scroll(-1, "units")
        elif event.num == 5 or event.delta < 0:
            self.canvas.yview_scroll(1, "units")

    def on_shift_mouse_wheel(self, event):
        # For Mac touchpad, use event.delta directly
        if event.num == 4 or event.delta > 0:
            self.canvas.xview_scroll(-1, "units")
        elif event.num == 5 or event.delta < 0:
            self.canvas.xview_scroll(1, "units")


    def add_entry(self, parent, label_text, row, callback=None):
        label = tk.Label(parent, text=label_text)
        label.grid(row=row, column=0, padx=10, pady=5)

        entry = tk.Entry(parent)
        entry.grid(row=row, column=1, padx=10, pady=5)

        self.entries[label_text] = entry
        if callback:
            entry.bind("<KeyRelease>", callback)
        return entry
    
    def toggle_pressure_inputs(self):
        if self.differential_pumping_var.get():
            self.constant_pressure_label.grid_remove()
            self.constant_pressure_entry.grid_remove()
            self.pressure_start_label.grid()
            self.pressure_boundaries_start.grid()
            self.pressure_end_label.grid()
            self.pressure_boundaries_end.grid()
        else:
            self.constant_pressure_label.grid()
            self.constant_pressure_entry.grid()
            self.pressure_start_label.grid_remove()
            self.pressure_boundaries_start.grid_remove()
            self.pressure_end_label.grid_remove()
            self.pressure_boundaries_end.grid_remove()

    def toggle_epp_inputs(self):
        if self.allocate_epp_var.get():
            self.fundamental_pulse_energy_label.grid()
            self.fundamental_pulse_energy_entry.grid()
            self.second_harmonic_pulse_energy_label.grid()
            self.second_harmonic_pulse_energy_entry.grid()
        else:
            self.fundamental_pulse_energy_label.grid_remove()
            self.fundamental_pulse_energy_entry.grid_remove()
            self.second_harmonic_pulse_energy_label.grid_remove()
            self.second_harmonic_pulse_energy_entry.grid_remove()
    

    def main(self, SHG_efficiency, FWHM, fiber_len, fiber_rad, fiberWL, alpha, peak_intensity_CGS,
             pulseWL, pulse_shape, frep_MHz, GDD, TOD, gas, constant_pressure, pressure_boundaries, differential_pumping, 
             manually_allocate_epp, fundamental_pulse_energy, second_harmonic_pulse_energy,
             Window, Points, Steps, Raman, Steep, pressure_in_torr, reload_fiber_each_step):
        
        if not manually_allocate_epp:
            # build the pulse properties
            intensitybuilder = PeakIntensityBuilder(SHG_efficiency=SHG_efficiency, 
                                            pulse_duration=FWHM*1e-12, 
                                            peak_intensity_CGS=peak_intensity_CGS)
            # allocate pulse energy
            total_EPP, EPP_FD, EPP_SH = intensitybuilder.calcEPP(fiber_radius=fiber_rad*1e-6)
        else:
            # manuallu set pulse energy
            EPP_FD, EPP_SH = fundamental_pulse_energy, second_harmonic_pulse_energy
        
        # load gas parameters
        n2_at_one_bar, betas_at_one_bar = loadGasParameters(gas, noble_gases)
        
        # build the gas properties using static field
        gasbuilder = GasPropertyBuilder(fiber_length=fiber_len*1e-3, fiber_radius=fiber_rad*1e-6, pulse_wavelength=pulseWL, 
                                        betas_at_one_bar=betas_at_one_bar, constant_pressure=constant_pressure, pressure_boundaries=pressure_boundaries, 
                                        n2_at_one_bar=n2_at_one_bar, differential_pumping=differential_pumping, pressure_in_Torr=pressure_in_torr)

        # create the fiber
        fiber = pynlo.media.fibers.fiber.FiberInstance()
        if differential_pumping:
            fiber.generate_fiber(fiber_len*1e-3, center_wl_nm=fiberWL, betas=gasbuilder.dispersionFunction(0), 
                                 gamma_W_m=gasbuilder.gammaDistribution(0), gvd_units='ps^n/m', gain=-alpha)
            fiber.set_gamma_function(gasbuilder.gammaDistribution)
            fiber.set_dispersion_function(gasbuilder.dispersionFunction, dispersion_format='GVD')
        else:
            fiber.generate_fiber(fiber_len*1e-3, center_wl_nm=fiberWL, betas=gasbuilder.dispersionFunction(), 
                                 gamma_W_m=gasbuilder.gammaDistribution(), gvd_units='ps^n/m', gain=-alpha)

        # create the fundamental pulse
        if pulse_shape.lower()=="sech":
            pulse_FD = pynlo.light.DerivedPulses.SechPulse(
                1, FWHM/1.76, pulseWL, time_window_ps=Window,
                GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=frep_MHz, power_is_avg=True
            )
        elif pulse_shape.lower()=="gaussian":
            pulse_FD = pynlo.light.DerivedPulses.GaussianPulse(
                1, FWHM, pulseWL, time_window_ps=Window,
                GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=frep_MHz, power_is_avg=True
            )
        pulse_FD.set_epp(EPP_FD)

        # create the second harmonic pulse
        if pulse_shape.lower()=="sech":
            pulse_SH = pynlo.light.DerivedPulses.SechPulse(
                1, FWHM/1.76, pulseWL/2, time_window_ps=Window,
                GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=frep_MHz, power_is_avg=True
            )
        elif pulse_shape.lower()=="gaussian":
            pulse_SH = pynlo.light.DerivedPulses.GaussianPulse(
                1, FWHM, pulseWL/2, time_window_ps=Window,
                GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=frep_MHz, power_is_avg=True
            )

        pulse_SH = pulse_SH.interpolate_to_new_center_wl(pulseWL)
        pulse_SH.set_epp(EPP_SH)

        # create the combined pulse
        pulse = pynlo.light.PulseBase.Pulse()
        pulse.set_NPTS(Points)
        pulse.set_time_window_ps(Window)
        pulse.set_frep_MHz(frep_MHz)
        pulse.set_center_wavelength_nm(fiberWL)
        pulse.set_AW(pulse_FD.AW + pulse_SH.AW)

        # propagation
        evol = pynlo.interactions.FourWaveMixing.SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True, 
                                                           disable_Raman=np.logical_not(Raman), 
                                                           disable_self_steepening=np.logical_not(Steep))

        y, AW, AT, pulse_out = evol.propagate(pulse_in=pulse, fiber=fiber, n_steps=Steps, reload_fiber_each_step=reload_fiber_each_step)
        return y, AW, AT, pulse, pulse_out, pulse_FD, pulse_SH
    

    def plot_results(self, y, AW, AT, pulse, pulse_out, pulse_FD, pulse_SH, fiber_len):
        y_mm = y * 1e3
        F = (pulse_SH.W_mks / (2 * np.pi)) * 1e-12  # Convert to THz
        W_eV = photonEnergyGrid(pulse_SH)
        
        zW = dB(np.transpose(AW)[:, (W_eV > 0)])
        zT = dB(np.transpose(AT))
        
        # Spectral power density plot (output pulse)
        plt.figure(figsize=(8, 5))
        plt.plot(W_eV[W_eV > 0], np.abs(pulse_out.AW[W_eV > 0])**2, color='darkred', label='Output Pulse')
        plt.xlabel('Photon Energy (eV)', fontsize=14)
        plt.ylabel('Spectral Power Density (W/eV)', fontsize=14)
        plt.yscale('log')
        plt.xlim(0, 20)
        plt.ylim(1e-15, None)
        plt.grid(True, which='both', linestyle='--', linewidth=0.7)
        plt.tight_layout()
        plt.savefig('plot1.png', dpi=300)
        plt.close()
        
        # Temporal power density plot (output pulse)
        plt.figure(figsize=(8, 5))
        plt.plot(pulse_SH.T_ps, np.abs(pulse_out.AT)**2, color='darkred', label='Output Pulse')
        x_ticks = np.arange(np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps) + 1, 1)
        plt.xticks(x_ticks)
        plt.xlim(np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps))
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('Temporal Power Density (W/eV)', fontsize=14)
        plt.yscale('log')
        plt.grid(True, which='both', linestyle='--', linewidth=0.7)
        plt.tight_layout()
        plt.savefig('plot2.png', dpi=300)
        plt.close()
        
        # Spectral broadening along fiber plot
        plt.figure(figsize=(10, 5))
        extent = (np.min(W_eV[W_eV > 0]), np.max(W_eV[W_eV > 0]), 0, fiber_len)
        plt.imshow(zW, extent=extent, vmin=np.max(zW) - 250.0, vmax=np.max(zW), aspect='auto', origin='lower', cmap='RdBu_r')
        plt.xlim(0, 20)
        plt.xlabel('Photon Energy (eV)', fontsize=14)
        plt.ylabel('Propagation Distance (mm)', fontsize=14)
        plt.colorbar(label='Spectral Power Density (dB)')
        plt.tight_layout()
        plt.savefig('plot3.png', dpi=300)
        plt.close()
        
        # Temporal broadening along fiber plot
        plt.figure(figsize=(9, 6))
        extent = (np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps), 0, fiber_len)
        plt.imshow(zT, extent=extent, vmin=np.max(zT) - 60.0, vmax=np.max(zT), aspect='auto', origin='lower', cmap='RdBu_r')
        x_ticks = np.arange(np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps) + 1, 1)
        plt.xticks(x_ticks)
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('Propagation Distance (mm)', fontsize=14)
        plt.colorbar(label='Temporal Power Density (dB)')
        plt.tight_layout()
        plt.savefig('plot4.png', dpi=300)
        plt.close()
        
        # Combined plot
        fig = plt.figure(figsize=(15, 10))
        ax0 = plt.subplot2grid((2, 2), (0, 0))
        ax1 = plt.subplot2grid((2, 2), (0, 1))
        ax2 = plt.subplot2grid((2, 2), (1, 0), sharex=ax0)
        ax3 = plt.subplot2grid((2, 2), (1, 1), sharex=ax1)

        # Spectral power density
        ax0.plot(W_eV[W_eV > 0], np.abs(pulse_out.AW[W_eV > 0])**2, color='darkred')
        ax0.set_xlim(0, 20)
        ax0.set_yscale('log')
        ax0.set_ylabel('Spectral Power Density (W/eV)', fontsize=14)

        # Temporal power density
        ax1.plot(pulse_SH.T_ps, np.abs(pulse_out.AT)**2, color='darkred')
        ax1.set_yscale('log')
        ax1.set_xticks(np.arange(np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps) + 1, 1))
        ax1.set_ylabel('Temporal Power Density (W/eV)', fontsize=14)

        # Spectral broadening along fiber
        extent = (np.min(W_eV[W_eV > 0]), np.max(W_eV[W_eV > 0]), 0, fiber_len)
        ax2.imshow(zW, extent=extent, vmin=np.max(zW) - 250.0, vmax=np.max(zW), aspect='auto', origin='lower', cmap='RdBu_r')
        ax2.set_ylabel('Propagation Distance (mm)', fontsize=14)
        ax2.set_ylim(1e-18,None)

        # Temporal broadening along fiber
        extent = (np.min(pulse_SH.T_ps), np.max(pulse_SH.T_ps), 0, fiber_len)
        ax3.imshow(zT, extent=extent, vmin=np.max(zT) - 60.0, vmax=np.max(zT), aspect='auto', origin='lower', cmap='RdBu_r')
        ax3.set_ylabel('Propagation Distance (mm)', fontsize=14)

        ax2.set_xlabel('Photon Energy (eV)', fontsize=14)
        ax3.set_xlabel('Time (ps)', fontsize=14)

        plt.tight_layout()
        plt.savefig('combined_plot.png', dpi=300)
        plt.close()
        
        # Display the plots in the respective tabs
        self.display_plot_in_tab('plot1.png', self.tab1)
        self.display_plot_in_tab('plot2.png', self.tab2)
        self.display_plot_in_tab('plot3.png', self.tab3)
        self.display_plot_in_tab('plot4.png', self.tab4)
        self.display_plot_in_tab('combined_plot.png', self.tab5)

    def display_plot_in_tab(self, img_path, tab):
        for widget in tab.winfo_children():
            widget.destroy()

        img = Image.open(img_path)
        img = ImageTk.PhotoImage(img)
        panel = tk.Label(tab, image=img)
        panel.image = img
        panel.pack(side="top", fill="both", expand="yes")

    def run_simulation(self):

        # close any existing progress window
        if self.progress_window is not None:
            self.progress_window.destroy()

        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Simulation in Progress")

        center_window(self.progress_window, 300, 100)
        tk.Label(self.progress_window, text="Simulation is running, please wait...").pack(padx=20, pady=20)
        self.progress_window.update()
        self.run_simulation_thread()
    
    def run_simulation_thread(self):
        # get laser properties
        pulseWL = float(self.entries['Pulse Wavelength [nm]'].get())
        FWHM = float(self.entries['FWHM Duration [ps]'].get())
        frep_MHz = float(self.entries['Repetition Rate [MHz]'].get())
        total_EPP = float(self.entries['Energy Per Pulse [mJ]'].get()) * 1e-3
        pulse_shape = self.entries['Pulse Shape'].get()
        GDD = float(self.entries['Group Delay Dispersion [ps^2]'].get())
        TOD = float(self.entries['Third Order Dispersion [ps^3]'].get())

        # get fiber properties
        fiber_len = float(self.entries['Fiber Length [mm]'].get())
        fiber_rad = float(self.entries['Fiber Radius [microns]'].get())
        fiberWL = float(self.entries['Fiber Wavelength [nm]'].get())
        Alpha = float(self.entries['Attentuation Coefficient [dB/cm]'].get())
        alpha = np.log((10**(Alpha * 0.1))) * 100  # convert from dB/cm to 1/m

        # get nonlinear optics parameters
        peak_intensity_CGS = 2e12  # peak intensity for nonlinear processes [W/cm^2]
        SHG_efficiency = float(self.entries['SHG Conversion Efficiency'].get())
        Raman = self.raman_var.get()
        Steep = self.steep_var.get()

        # get gas properties
        gas = self.entries['Gas Name'].get().lower()
        differential_pumping = self.differential_pumping_var.get()
        pressure_unit = self.pressure_unit_var.get()
        pressure_in_torr = (pressure_unit == "Torr")

        if differential_pumping:
            pressure_boundaries = [
                float(self.entries['Pressure At Entrance'].get()),
                float(self.entries['Pressure At Exit'].get())
            ]
            constant_pressure = None
        else:
            constant_pressure = float(self.entries['Constant Pressure'].get())
            pressure_boundaries = None

        # get EPP allocation parameters
        manually_allocate_epp = self.allocate_epp_var.get()
        fundamental_pulse_energy = float(self.fundamental_pulse_energy_entry.get())*1e-6 if manually_allocate_epp else None
        second_harmonic_pulse_energy = float(self.second_harmonic_pulse_energy_entry.get())*1e-6 if manually_allocate_epp else None
        
        # get simulation parameters
        Steps = float(self.entries['Simulation Steps'].get())
        Window = float(self.entries['Simulation Window [ps]'].get())
        Points = float(self.entries['Simulation Points'].get())
        reload_fiber_each_step = self.reload_fiber_each_step_var.get()

        y, AW, AT, pulse, pulse_out, pulse_FD, pulse_SH = self.main(SHG_efficiency, FWHM, fiber_len, fiber_rad, fiberWL, alpha, peak_intensity_CGS,
                                                                    pulseWL, pulse_shape, frep_MHz, GDD, TOD, gas, constant_pressure, pressure_boundaries, differential_pumping, 
                                                                    manually_allocate_epp, fundamental_pulse_energy, second_harmonic_pulse_energy,
                                                                    Window, Points, Steps, Raman, Steep, pressure_in_torr, reload_fiber_each_step)
        # plot results
        self.plot_results(y, AW, AT, pulse, pulse_out, pulse_FD, pulse_SH, fiber_len)
        self.progress_window.destroy()
    
    def on_closing(self):
        self.root.quit()
        self.root.destroy()
        
if __name__ == '__main__':
    root = tk.Tk()
    app = SimulationApp(root)
    root.mainloop()