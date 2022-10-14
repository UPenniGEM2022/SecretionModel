"""
-list all rate constants for creation, degradation, transport, etc.

-list all initial concentrations

-this must be done for both lER and mER systems

-put in a light administration regimen (matrix?)
======

for loop to step through all compartment equations (ER, cis Golgi, trans Golgi, membrane/cell external)

the PhoCl matrix should have function that can switch protein from membrane to lumen, etc. based on PhoCl transfer function 
"""

#PACKAGES
from __future__ import annotations
import webbrowser
from openpyxl import load_workbook
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import math
import numpy as np
from tkinter import *
from tkinter import ttk
import tkinter.messagebox
from tkinter.tix import *  #used for the indicator window
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo
from PIL import ImageTk, Image, ImageDraw
from scipy import integrate #use solve_ivp()
from tkinter.font import Font, nametofont
import platform

# CONSTANTS

# Determine platform
SYSTEM = platform.system().lower()

# Font sizes
if SYSTEM == "linux":
    REG_FONT_SIZE = 8
    HEADER_FONT_SIZE = 10
else:
    REG_FONT_SIZE = 20
    HEADER_FONT_SIZE = 20

# Fonts
DEF_REG_FONT = "lucida"
HEADER_FONT = "lucida"

# Window size
WX = 2400
WY = 1400

# Other sizes
SLIDER_ENTRY_WIDTH = 6
SLIDER_LABEL_WIDTH = 15

# Website
WEBSITE_ROOT = "https://google.com/"
EQUATIONS_PAGE = ""

# Default values
def_k_prod: float = 0
def_k_cell_to_plasma: float = 0
def_k_plasma_to_cell: float = 0     
def_k_lumen_to_out: float = 0     
def_k_out_to_lumen: float = 0       
def_c_cell: float = 0
def_c_plasma: float = 0
def_k_degradation: float = 0
C_CELL0: float = 0
C_PLASMA0: float = 0
C_LUMEN0: float = 0
C_OUT0: float = 0

#FIRST, DEFINE ALL FUNCTIONS FOR ALL ORGANELLES FOR THE SOLVE_IVP CODE
#stand-in function for now:
def f_shed(t, y, k_prod, k_cell_to_plasma, k_plasma_to_cell, k_lumen_to_out, k_out_to_lumen, c_cell, c_plasma, k_deg, light_regimen: list[tuple[str, float]]):
    #concentrations of proteins in various compartments
    mem_conc = y[0]
    plasma_conc = y[1]
    lumen_conc = y[2]
    out_conc = y[3]

    #light regimen indicator function
    chi = 0
    if t == 0:
        chi = 1
    else:
        for i in range(len(light_regimen)):
            if light_regimen[i][1] > t:
                if light_regimen[i-1][0] == "on":
                    chi = 1
                break

    #final derivative equations for the concentrations
    d_cell = k_prod + k_plasma_to_cell*plasma_conc - k_cell_to_plasma*mem_conc - chi*c_cell*mem_conc
    d_plasma = k_cell_to_plasma*mem_conc - k_plasma_to_cell*plasma_conc - chi*c_plasma*plasma_conc
    d_lumen = k_out_to_lumen*out_conc - k_lumen_to_out*lumen_conc + chi*c_cell*mem_conc
    d_out = k_lumen_to_out*lumen_conc - k_out_to_lumen*out_conc + chi*c_plasma*plasma_conc

    out = np.array([d_cell, d_plasma, d_lumen, d_out])
    out += np.array([mem_conc, plasma_conc, lumen_conc, out_conc])*k_deg
    
    return out
   #https://www.youtube.com/watch?v=Gg--FOdupwY super helpful video to set up the final form of the function :)
 

def rmseCalc(te, y_actual, ts, y_predicted):
    # Do a linear interpolation of predicted data to check against experimental data
    n = len(te)
    ns = len(ts)
    if n != len(y_actual):
        raise ValueError("Something wrong with experimental data in RMSE")

    if ns != len(y_predicted):
        raise ValueError("Something wrong with simulated data in RMSE")

    sidx = 0
    mse = 0
    for i in range(n):
        t = te[i]
        for j in range(sidx, ns):
            if ts[j] > t:
                sidx = j
                break
        
        # Linear interpolation
        t0 = ts[sidx-1]
        t1 = ts[sidx]
        y0 = y_predicted[sidx-1]
        y1 = y_predicted[sidx]
        y = (y0*(t1 - t) + y1*(t - t0))/(t1 - t0)
        
        mse += (y_actual[i] - y)**2
    
    mse /= n
    return math.sqrt(mse)

class VarSlider(Frame):
    inputWindow: InputWindow

    def __init__(self, win, inputWindow: InputWindow, length: int, bg: str, label: str, callback=lambda x: None, from_def: float=0, to_def: float=100):
        super().__init__(win, bg=bg)
        self.inputWindow = inputWindow
        self.callback = callback

        self.label = Label(self, text=label + ":", bg=bg, width=SLIDER_LABEL_WIDTH)
        self.from_input = Entry(self, width=SLIDER_ENTRY_WIDTH)
        self.from_input.insert(0, str(from_def))
        self.slider = Scale(self, length=length, bg=bg, orient=HORIZONTAL, command=self.handle_callback)
        self.to_input = Entry(self, width=SLIDER_ENTRY_WIDTH)
        self.to_input.insert(0, str(to_def))

        # self.columnconfigure(0, weight=20)
        # self.columnconfigure(1, weight=1)
        # self.columnconfigure(2, weight=2)
        # self.columnconfigure(3, weight=1)

        self.label.grid(column=0, row=0)
        self.from_input.grid(column=1, row=0, sticky="e")
        self.slider.grid(column=2, row=0)
        self.to_input.grid(column=3, row=0, sticky="w")
    
    def disable(self):
        self.slider['state'] = 'disabled'
        self.from_input['state'] = 'disabled'
        self.to_input['state'] = 'disabled'
    
    def enable(self):
        self.slider['state'] = 'normal'
        self.to_input['state'] = 'normal'
        self.from_input['state'] = 'normal'
    
    def handle_callback(self, value):
        self.slider["from"] = self.from_input.get()
        self.slider["to"] = self.to_input.get()
        self.callback(value)
    
    def get(self):
        return self.slider.get()

class InputWindow:
    var_sliders: list[VarSlider]
    rates: dict[str, float]

    def __init__(self, win):
        # Instantiate rates
        self.rates = {"k_prod": def_k_prod,
                    "k_cell_to_plasma": def_k_cell_to_plasma,
                    "k_plasma_to_cell": def_k_plasma_to_cell,
                    "k_lumen_to_out": def_k_lumen_to_out,
                    "k_out_to_lumen": def_k_out_to_lumen,
                    "c_cell": def_c_cell,
                    "c_plasma": def_c_plasma,
                    "k_degradation": def_k_degradation}
        
        # Create objects of tkinter ImageTk for use in GUI
        # Load the image
        #info_img=Image.open("./assets/info_img.jpg")
        #info_img = info_img.resize((20, 20))
        #info_img=ImageTk.PhotoImage(info_img)
        # Parameters
        #self.font = Font(size=12)
        


        #LIGHT FRAME 
        m = 0 #row index for the light frame
        #define the subframe
        light_frame = Frame(win)
        #define labels, text entries, buttons
        light_title=Label(light_frame, text = 'Light Administration Regimen      ', font = f'{HEADER_FONT} {HEADER_FONT_SIZE} bold')
        lbl1=Label(light_frame, text='Time on (min, sec):')
        lbl2=Label(light_frame, text='Time off (min, sec):')
        lbl3=Label(light_frame, text='# Cycles:')
        lbl4=Label(light_frame, text = 'Time after last cycle (min, sec):')
        self.t1=Entry(light_frame, bd=7)
        self.t1.insert(0, "0,0")
        self.t2=Entry(light_frame, bd=7)
        self.t2.insert(0, "0,0")
        self.t3=Entry(light_frame, bd=7)
        self.t3.insert(0, "1")
        self.t4=Entry(light_frame, bd=7)
        self.t4.insert(0, "0,0")
        #Place all elements in subframe
        light_title.grid(row = m, column = 0, columnspan = 2, ipady = 20)
        m =+ 1
        lbl1.grid(row = m, column = 0)
        self.t1.grid(row = m, column = 1, padx=20)
        m += 1
        lbl2.grid(row = m, column = 0)
        self.t2.grid(row = m, column = 1)
        m += 1
        lbl3.grid(row = m, column = 0)
        self.t3.grid(row = m, column = 1)
        m += 1
        lbl4.grid(row = m, column = 0)
        self.t4.grid(row = m, column = 1)
        b2=Button(light_frame, bd=7, text='Reset', command = self.reset, cursor = 'hand2')
        m += 1
        b2.grid(row = m, column = 0, pady = 20, columnspan=2)
        light_frame.grid(row = 0, column = 0)
        Grid.columnconfigure(light_frame, 0, weight=1)
        Grid.columnconfigure(light_frame, 1, weight=1)
        # Create a Label Widget to display helpful tip - TODO
        # light_info = Label(light_frame, image = info_img)
        # light_info.image = info_img #keep a reference apparently??? ASK SAM
        # light_info.grid(row = 0, column = 1)
        #info text       
        # light_tip=Balloon(win)
        # light_tip.bind_widget(light_info, balloonmsg="Python is an interpreted, high-level and general-purpose programming language \n when does it switch to next line") #use \n to get to the next line
        

        
        #SELECT EXPERIMENTAL DATA
        experiment_frame = Frame(win)
        exp_title=Label(experiment_frame, text = 'Experimental Data (Excel Sheet)', font = f'{HEADER_FONT} {HEADER_FONT_SIZE} bold')
        exp_title.grid(row = 0, column = 0, columnspan = 2)
        bexp=Button(experiment_frame, bd = 7, text='Import Data', command=self.import_data, cursor = 'hand2')
        bexp.grid(row = 1, column = 0, columnspan = 2, pady= 20)
        experiment_frame.grid(row = 1, column = 0, pady = 20)
        #Grid.columnconfigure(experiment_frame, 0, weight=1)
        #Grid.columnconfigure(experiment_frame, 1, weight=1)


        #CONSTRUCT SELECTION
        construct_frame = Frame(win)
        #title the subframe
        con_title = Label(construct_frame, text = 'Select PhoCl Construct', font = f'{HEADER_FONT} {HEADER_FONT_SIZE} bold')
        con_title.grid(row = 0, column = 0, columnspan = 3)
        #define construct options
        var = IntVar()
        R_secrete = Radiobutton(construct_frame, text="PhoCl Secrete", variable=var, value=1, command=lambda: self.update_sliders_frame("secrete"))
        R_shed = Radiobutton(construct_frame, text="PhoCl Shed", variable=var, value=2, command=lambda: self.update_sliders_frame("shed"))
        R_disp = Radiobutton(construct_frame, text="PhoCl display", variable=var, value=3, command=lambda: self.update_sliders_frame("display"))
        #place construct options in frame
        R_secrete.grid(row = 1, column = 0, sticky = 'w')
        R_shed.grid(row = 1, column = 1, sticky = 'w')
        R_disp.grid(row = 1, column = 2, sticky = 'w')
        #place subframe in main window
        construct_frame.grid(row = 2, column = 0, pady = 20)
        


        #VARIABLE TUNERS
        s_length = 600
        s_color = 'light gray'
        tuner_frame = Frame(win, relief = 'sunken', bd = 10, bg = s_color)
        m = 0 #row index
        #title
        self.tuner_title=Label(tuner_frame, text = 'Rate Constants (PhoCl Secrete)', font = f'{HEADER_FONT} {HEADER_FONT_SIZE} bold', bg = s_color)
        self.tuner_title.grid(row = m, column = 0, columnspan = 2)
        m += 1
        #protein production slider
        self.k_prod_slider = VarSlider(tuner_frame, self, length=s_length, label='k_prod', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_prod", float(x)))
        self.k_prod_slider.grid(row=m, column = 0, columnspan = 2, sticky="ew")
        m += 1 

        self.k_c2p_slider = VarSlider(tuner_frame, self, length=s_length, label='k_mi->mo', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_cell_to_plasma", float(x)))
        self.k_c2p_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.k_p2c_slider = VarSlider(tuner_frame, self, length=s_length, label='k_mo->mi', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_plasma_to_cell", float(x)))
        self.k_p2c_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.k_l2o_slider = VarSlider(tuner_frame, self, length=s_length, label='k_li->lo', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_lumen_to_out", float(x)))
        self.k_l2o_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.k_o2l_slider = VarSlider(tuner_frame, self, length=s_length, label='k_lo->li', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_lumen_to_out", float(x)))
        self.k_o2l_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.c_cell_slider = VarSlider(tuner_frame, self, length=s_length, label='L_i', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("c_cell", float(x)))
        self.c_cell_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.c_plasma_slider = VarSlider(tuner_frame, self, length=s_length, label='L_o', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("c_plasma", float(x)))
        self.c_plasma_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1

        self.k_degrad_slider = VarSlider(tuner_frame, self, length=s_length, label='k_degradation', from_def=0, to_def=100, bg=s_color, callback=lambda x: self.update_rates("k_degradation", float(x)))
        self.k_degrad_slider.grid(row=m, column=0, columnspan=2, sticky="ew")
        m += 1


        #clickable links
        # self.scale_range = Label(tuner_frame, text = 'Change range of sliders', fg = 'blue', bg = s_color, font = "lucida 16 underline", pady = 10, cursor = 'hand2')
        # self.scale_range.grid(row = m, column = 0)
        self.eq_list = Label(tuner_frame, text = 'See construct equations', fg = 'blue', bg = s_color, font = "lucida 16 underline", pady = 10, cursor = 'hand2')
        self.eq_list.grid(row = m, column = 0, columnspan=2)
        # self.scale_range.bind("<Button-1>", self.onTextClick)
        self.eq_list.bind("<Button-1>", lambda _: webbrowser.open(f"{WEBSITE_ROOT}/{EQUATIONS_PAGE}", new=True))
        m += 1
        #final button
        self.bs=Button(tuner_frame, bd=7, text='Run Simulation', command=self.submit, cursor = 'hand2')
        self.bs.grid(row = m, column = 0, columnspan = 2, padx = 60, pady = 20, sticky = 'nsew')
        tuner_frame.grid(row = 3, column = 0, pady = 20)

        # Put sliders into list for ease of use later
        self.var_sliders = [self.k_prod_slider, self.k_c2p_slider, self.k_degrad_slider, self.k_l2o_slider, self.k_o2l_slider, self.k_p2c_slider, self.c_cell_slider, self.c_plasma_slider]

        # Select correct construct radio button
        R_shed.select()
        self.update_sliders_frame("shed")

        #PLOT
        matplotlib_frame = Frame(win)
        # Some example data to display
        self.x_exp = np.linspace(0, 2 * np.pi, 400)
        self.y_exp = self.x_exp * 0
        self.x_sim = np.linspace(0, 2 * np.pi, 400)
        self.y_sim = self.x_sim * 0
        fig, self.ax = plt.subplots()
        self.ax.plot(self.x_exp, self.y_exp, label="experimental")
        self.ax.plot(self.x_sim, self.y_sim, label="simulation")
        self.ax.set_xlabel('Time', fontsize=16)
        self.ax.set_ylabel('Concentration/Intensity', fontsize=16)
        self.ax.legend()
        self.ax.set_title('Experimental and Simulation Plots', fontsize=20)             
        self.canvas = FigureCanvasTkAgg(fig, master=matplotlib_frame)  # A tk.DrawingArea.
        self.canvas.draw()
        # pack_toolbar=False will make it easier to use a layout manager later on.
        toolbar = NavigationToolbar2Tk(self.canvas, matplotlib_frame, pack_toolbar=True)
        toolbar.update()
        matplotlib_frame.grid(row = 0, column = 1, rowspan = 4, sticky='nsew')
        #toolbar.pack(side=BOTTOM, fill=X)
        self.canvas.get_tk_widget().pack(side=TOP, expand=True, fill="both") #fill=BO



        #RMSE FRAME
        rmse_frame = Frame(win)
        self.rmse_label = Label(rmse_frame, text = 'RMSE: ')
        self.rmse_label.grid(row = 0, column = 0)
        rmse_frame.grid(row = 5, column = 1, sticky='ew')
        
     #===============================FUNCTIONS=====================================================   

    def import_data(self):
        #delete current experimental values
        self.x_exp = []
        self.y_exp = []
        #ask user to input excel file
        file_name = fd.askopenfilename(
            title='Open a file',
            initialdir='/')
        showinfo(
            title='Selected File',
            message=file_name
        )
        workbook = load_workbook(filename=file_name)
        sheet = workbook.active #assumption that there is only one sheet in a given excel file
        #take cell values and put into x_exp and y_exp
        columnA = sheet['A']  # Column
        self.x_exp = [columnA[i].value for i in range(len(columnA))]
        columnB = sheet['B']  # Column
        self.y_exp = [columnB[i].value for i in range(len(columnB))]
        
        #clear the plot of old values
        plt.cla()
        #replot everything
        self.ax.set_xlabel('Time', fontsize=16)
        self.ax.set_ylabel('Concentration/Intensity', fontsize=16)
        self.ax.plot(self.x_exp, self.y_exp, marker='o', label="experimental")
        self.ax.plot(self.x_sim, self.y_sim, label="simulation")
        self.ax.legend()
        self.ax.set_title('Experimental and Simulation Plots', fontsize=20) 
        self.canvas.draw()
        #clear RMSE values from display
        self.rmse_label["text"] = "RMSE: "
        


    def submit(self):
        # STUFF FROM OLD submit FUNCTION
        dt = 1 #TEMPORARILY
        #take each text entry to get string
        onTime = str(self.t1.get())
        offTime = str(self.t2.get())
        numCycles = int(self.t3.get())
        afterTime = str(self.t4.get())
        #separate minute and sec values
        onTime = onTime.split(",")
        offTime = offTime.split(",")
        afterTime = afterTime.split(",")
        #convert the str values to int
        onTime = [int(i) for i in onTime] #would prefer to only use the first 2 values
        offTime = [int(i) for i in offTime]
        afterTime = [int(i) for i in afterTime]
        # Convert to seconds
        onTime = onTime[0]*60 + onTime[1]
        offTime = offTime[0]*60 + offTime[1]
        afterTime = afterTime[0]*60 + afterTime[1]

        # Get list of on and off times
        lightCourse: list[tuple[str, float]] = []
        t0 = t = 0
        for _ in range(numCycles):
            # Light starts on
            lightCourse.append(("on", t))
            # Turns off
            t += onTime
            lightCourse.append(("off", t))
            t += offTime
        
        t1 = t + afterTime
        
        #get array of milliseconds that say when the light is on... if the current time point evaluation is between 
        #light on or off... if not between 2 nonzero values, then light is off?

        # STUFF FROM OLD UPDATE_VAR FUNCTION
        #delete current simulation values
        self.x_sim = []
        self.y_sim = []
        #set times
        t_span = np.array([t0, t1])
        #set initial conditions
        y0 = np.array ([C_CELL0, C_LUMEN0, C_PLASMA0, C_OUT0]) #will probably somehow need to come from the system of equations themselves, solved with an external funciton
        #Solve IVP both at intended points and at 
        res = integrate.solve_ivp(f_shed, t_span, y0, args=[self.rates[x] for x in self.rates.keys()] + [lightCourse]) #TO HAVE AN ADAPTIVE SMOOTH TIME.... STUPID???
        ts = res.t
        c_cell, c_plasma, c_lumen, c_out = res.y
        #clear the plot of old values
        plt.cla()
        #replot 
        self.x_sim = ts
        self.y_sim = c_out
        self.ax.set_xlabel('Time', fontsize=16)
        self.ax.set_ylabel('Concentration/Intensity', fontsize=16)
        #replot experimental times
        self.ax.plot(self.x_exp, self.y_exp, marker='o', label="experimental", color = 'blue')
        #replot t_eval sim times
        self.ax.scatter(self.x_sim, self.y_sim, label="simulation", color = 'red')
        self.ax.legend()
        self.ax.set_title('Experimental and Simulation Plots', fontsize=20) 
        self.canvas.draw()
        #recalculate RMSE with updated experimental data
        RMSE_val = rmseCalc(self.x_exp, self.y_exp, ts, c_out)
        RMSE = str(RMSE_val)
        self.rmse_label["text"] = f"RMSE: {RMSE}"
        
    def reset(self):
        #get rid of all text entries
        self.t1.delete(0, 'end')
        self.t2.delete(0, 'end')
        self.t3.delete(0, 'end')
        self.t4.delete(0, 'end')
        #replace default entries
        self.t1.insert(0, "0,0")
        self.t2.insert(0, "0,0")
        self.t3.insert(0, "1")
        self.t4.insert(0, "0,0")
    
    def update_sliders_frame(self, construct: str):
        if construct == "shed":
            for slider in self.var_sliders:
                slider.enable()
            self.tuner_title['text'] = 'Rate Constants (PhoCl Shed)'
            self.bs['state'] = 'normal'
        else:
            for slider in self.var_sliders:
                slider.disable()
            self.tuner_title['text'] = f'Rate Constants (PhoCl {"".join([construct[x].upper() if x == 0 else construct[x] for x in range(len(construct))])}) - UNIMPLEMENTED'
            self.bs['state'] = 'disabled'
    
    def update_slider_lims(self, n: int, slider_from: float|None, slider_to: float|None):
        if n < 0 or n > len(self.var_sliders):
            return
        
        if slider_from:
            self.var_sliders[n]['from'] = slider_from
        if slider_to:
            self.var_sliders[n]['to'] = slider_to
    
    def update_rates(self, rate: str, value: float):
        self.rates[rate] = value

    def onTextClick(self, event):
        tkinter.messagebox.showinfo("Welcome to GFG.",  "Hi I'm your message")
        

    def onConstructChoice(self, event):
        tkinter.messagebox.showinfo("Welcome to GFG.",  "Hi I'm your message")
        #if radiobutton is 1, or 2, or 3, redo the grid placing of the slider frame!

        
       

window=Tk()

#prevent fullscreen
window.resizable(False, False)

#set default font of the GUI
window.option_add( "*font", f"lucida {REG_FONT_SIZE}" )

inputwin=InputWindow(window)
window.title('Secretion Sim')

#Set the geometry of frame
window.geometry(f"{WX}x{WY}")

#Change the default Font that will affect in all the widgets
default_font = nametofont("TkDefaultFont")
default_font.configure(size=16)
text_font = nametofont("TkTextFont")
text_font.configure(size=16)
fixed_font = nametofont("TkFixedFont")
fixed_font.configure(size=16)

#change the row and column weightings 
Grid.rowconfigure(window,0,weight=1)
#Grid.columnconfigure(window,0,weight=1)
Grid.rowconfigure(window,1,weight=1)
Grid.columnconfigure(window, 1, weight=3)

window.mainloop()

#EVENTUALLY, USE THE MAIN WINDOW CLASS TO CREATE TABS FOR BOTH THE SPATIAL AND TEMPORAL DYNAMICS



#VARIABLE CONCENTRATIONS
#ER
P_mER = 4
P_lER = 4
#CG
P_mCG = 4
P_lCG = 4
#TG
P_mTG = 4
P_lTG = 4
#EC (extracellular)
P_mEC = 4
P_lEC = 4

'''CONSTANTS===============================

#PRODUCTION 
alpha_mER = 4 #rate of protein production / transport to ER that is NOT from cis Golgi
alpha_lER = 4 #rate of degradation/transport to lysozome from ER

#DEGRADATION
#ER
delta_mER = 4 
delta_lER = 4
#Cis Golgi
delta_mCG = 4
delta_lCG = 4
#Trans Golgi
delta_mTG = 4
delta_lTG = 4
#Extracellular
delta_mEC = 4
delta_lEC = 4

#TRANSPORT
#ER <=> CG
tau_mCG_mER = 4 #rate of retrograde transport of POI from mCG to mER
tau_mER_mCG = 4 #rate of anterograde transport of POI from mER to mCG
tau_lCG_lER = 4 
tau_lER_lCG = 4 
#CG <=> TG
tau_mCG_mTG = 4 
tau_mTG_mCG = 4 
tau_lCG_lTG = 4 
tau_lTG_lCG = 4 
#TG <=> EC
tau_mTG_mEC = 4 
tau_mEC_mTG = 4 
tau_lTG_lEC = 4 
tau_lEC_lTG = 4  

#PHOCL STRENGTH MULTIPLIER
#ER
C_mER = 4
C_lER = 4
#Cis Golgi
C_mCG = 4
C_lCG = 4
#Trans Golgi
C_mTG = 4
C_lTG = 4
#Extracellular
C_mEC = 4


#COMBINE VARIABLES TO SIMPLIFY FINAL EQUATIONS:

'''
#https://www.codegrepper.com/code-examples/python/slider+python




   

#graph the measured concentrations and the ODE predictions in the same graph
#sliders will rerun the ODE sim and replot.

"""
SLiders should adjust all variables, then rerun the entire simulation with precalculated timestep
https://www.pythontutorial.net/tkinter/tkinter-slider/
measurements should be compared to the theoretical model, and NRMSE value should be calculated and displayed by interpolating for the same time points? 
"""



#calculate RMSE values :)


