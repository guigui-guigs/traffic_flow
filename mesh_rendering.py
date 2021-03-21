from tkinter import *
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import numpy as np
#from LWR_new_schema_perodic import schemas_test

def test(event):
    print("Hello")

def func(alpha):
    fig = Figure(figsize=(5, 4), dpi=100)
    t = np.arange(0, 3, .01)
    fig.add_subplot(111).plot(t, 2 * np.sin(alpha * np.pi * t))
    return fig

def run(widget,canvas):
        
    #schemas_test("fixes","creneau","upwind", 50,0.004, 0.01, 1)
    clear_canvas(canvas)
    alpha = scale_vitesse.get()
    figure = func(alpha)
    canvas = FigureCanvasTkAgg(figure, master=mainapp)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def show_modal_window():
    messagebox.showerror("ERREUR","Une erreur est survenue lors de la compilation")

def clear_canvas(canvas):
    canvas.get_tk_widget().destroy()

mainapp = Tk()
mainapp.title("Trafic routier")
mainapp.minsize(1080,720)

screen_x = int(mainapp.winfo_screenwidth())
screen_y = int(mainapp.winfo_screenheight())
window_x = 1080
window_y = 720
posX = screen_x//2 - window_x//2
posY = screen_y//2 - window_y//2
geo = "{}x{}+{}+{}".format(window_x, window_y, posX, posY)
mainapp.geometry(geo)

# Création des frames 
frame_parameters = LabelFrame(mainapp,text="Paramètres",width=500,height=100,borderwidth=1)
frame_parameters.pack(side="right",fill="y",padx=20, pady=20)

msg="Nombre de mailles"
label = Label(mainapp,text=msg)
label.pack()

entry_nbmailles = Entry(mainapp)
entry_nbmailles.pack()

button_run = Button(mainapp, text="Run", width=15, command=lambda: run(mainapp, canvas))
button_run.pack()

button_pause = Button(mainapp, text="Pause", width=15)
button_pause.pack()

check_test = Checkbutton(frame_parameters, text="test")
check_test.pack()

radio_graph = Radiobutton(frame_parameters, text="Density chart", value=0)
radio_graph.pack()

radio_circuit = Radiobutton(frame_parameters, text="Circuit rendering", value=1)
radio_circuit.pack()

radio_mesh = Radiobutton(frame_parameters, text="Mesh rendering", value=2)
radio_mesh.pack()

scale_vitesse = Scale(frame_parameters,from_=10, to=100)
scale_vitesse.pack()

list_schemas = Listbox(frame_parameters)
list_schemas.insert(1, "Type upwind")
list_schemas.insert(2, "Type Lax-Friedrichs")
list_schemas.pack()

# Ici mettre un canvas vierge au début
figure = func(10)
canvas = FigureCanvasTkAgg(figure, master=mainapp)  
canvas.draw()
canvas.get_tk_widget().pack()

mainapp.mainloop()