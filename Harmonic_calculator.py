"""
This program solves the differential equation of harmonic motion then plots displacement and velocity curves over a 100 second time period. 
Complete with GUI for user convenience.
"""

# import libraries
from scipy.integrate import odeint
import numpy as np, matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import font as tkFont


def submit(): # define the main loop function with ODE solver and plotter
    try:
        # User inputs
        g = float(g_entry.get()) # viscous damping coefficient
        k = float(k_entry.get()) # spring constant
        m = float(m_entry.get()) # mass
        initial_displacement = float(id_entry.get()) # initial displacement

        # Set up the 2nd order differential equation of motion as a system of 1st order equations to pass into odeint
        def F(y, t):
            zeta = g/(2*np.sqrt(m*k))                   # Define the damping ratio
            w0 = np.sqrt(k/m)                           # Define the natural frequency
            dy = [y[1], -2*zeta*w0*y[1] - w0**2*y[0]]   # Creates a new list with elements from input vector y
        
            return dy # Sets up the ODE y" = -2*(zeta)*w0*y' - w0**2*y for odeint to solve below

        # Define calculated parameters
        behavior = str()             # Initializes an empty string
        zeta = g/(2*np.sqrt(m*k))    # Defines the damping ratio, zeta
        w0 = np.sqrt(m/k)            # Defines the natural frequency, omega zero

        # Categorize behavior
        if zeta > 1.0: 
                behavior = 'Overdamped'
        elif zeta == 1:
                behavior = 'Critically Damped'
        elif zeta == 0.0:
            behavior = 'Undamped'
        else: behavior = 'Underdamped'

        # # Solve the ODE and plot the solution
        t_min = 0; t_max = 100; dt = 0.01  # Defines the min and max values of time as well as spacing for odeint to run through
        t = np.arange(t_min, t_max, dt)   # Creates an array from 0 to t_max with t_max/dt elements

        # # Initial conditions
        initial_conditions = (initial_displacement, 0.0) # Creates a tuple (y(t), t) as the given set of initial conditions.

        # Solve ODE - odeint takes the function name and initial conditions and solves it numerically for given t values
        y = odeint(F, initial_conditions, t)

        # Format plot
        ax.clear()

        ax.set_title(behavior + ' Harmonic Oscillator for $\zeta =$' + str(np.round(zeta, 6)) +\
                    ', $\omega_0 =$' + str(np.round(w0, 6)), fontsize=28 ) # Creates a plot title
        ax.set_xlabel("Time $(s)$", fontsize=20)       # Generates x-axis label
        ax.set_ylabel("Amplitude $(m)$", fontsize=20)  # Generates y-axis label
        ax.plot(t, y[:,0])                   # Plots solution to ODE
        ax.plot(t, y[:,1])                   # Plots derivative of solution to ODE
        ax.legend(("Displacement $(m)$", "Velocity $(m/s)$"), fontsize=20)   # Adds a legend to the graph after plotting
        ax.grid()                            # Adds a background grid
        canvas.draw()

        # display calculated values
        result_label.config(text=f"Zeta = {zeta:.5f}, Omega = {w0:.5f}, Behavior is {behavior}")
        


    except ValueError:
        result_label.config(text=f"Please enter valid numbers")
 


### GUI 
        
# create GUI window with title
root = tk.Tk()
root.title("Harmonic Motion Calculator")
root.geometry("1600x1200") # window size
font = tkFont.Font(size=30) # font size

# damping coefficient entry
tk.Label(root, text="Viscous damping coefficient:", font=font).grid(row=0, column=0)
g_entry = tk.Entry(root, font=font)
g_entry.grid(row=0, column=1)
g_entry.insert(0, "2.11")  # Default value

# spring constant entry
tk.Label(root, text="Spring constant:", font=font).grid(row=1, column=0)
k_entry = tk.Entry(root, font=font)
k_entry.grid(row=1, column=1)
k_entry.insert(0, "1.5")  # Default value

# mass entry
tk.Label(root, text="Mass (kg):", font=font).grid(row=2, column=0)
m_entry = tk.Entry(root, font=font)
m_entry.grid(row=2, column=1)
m_entry.insert(0, "20.0")  # Default value

# initial displacement entry
tk.Label(root, text="Initial displacement (m):", font=font).grid(row=3, column=0)
id_entry = tk.Entry(root, font=font)
id_entry.grid(row=3, column=1)
id_entry.insert(0, "0.3")  # Default value

# submit button
submit_button = tk.Button(root, text="Submit", command=submit, font=font)
submit_button.grid(row=4, column=1)

# label for result
result_label = tk.Label(root, text="", font=font)
result_label.grid(row=5, column=1)

# matplotlib plot
fig, ax = plt.subplots(figsize=(20,12))
root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)
canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=6, column=0, columnspan=2, sticky='ew')

# start the GUI
root.mainloop()

