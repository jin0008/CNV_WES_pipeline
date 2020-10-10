# coding: utf-8

# Author: Julie bogoin

import  tkinter as tk
from tkinter import filedialog
import os

# Fenetre principale 
racine = tk.Tk()
racine.title("Pipeline CNV")    
racine.geometry("500x220+5+5")

def browse_button():
    # Allow user to select a directory and store it in global var called folder_path
    global folder_path
    filename = filedialog.askdirectory()
    os.chdir(filename)
    folder_path.set(filename)

def pipeline():
    os.system('sudo bash ~/CNV_WES_pipeline/pipeline_cnv.sh')

def resultats():
    os.system('sudo bash ~/CNV_WES_pipeline/pipeline_results.sh')

# Zone  dans  la  fenêtre  principale  où  on  écrit  un  texte
pathlabel = tk.Label(racine, text='\nQuel est le run a analyser?\n')
pathlabel.pack()

# Bouton parcourir
folder_path = tk.StringVar()
lbl1 = tk.Label(racine, textvariable = folder_path)
lbl1.pack()
button2 = tk.Button(text="Parcourir", command=browse_button)
button2.pack()

# Zone  dans  la  fenêtre  principale  où  on  écrit  un  texte
label = tk.Label(racine , text="\nQue souhaitez-vous faire?\n")
label.pack()

# Bouton pipeline
bouton_pipeline = tk.Button(racine , text="Pipeline CVN", command=pipeline)
bouton_pipeline["fg"] = "blue"
bouton_pipeline.pack()

# Bouton resultats
bouton_resultats = tk.Button(racine , text="Resultats CNV", command=resultats)
bouton_resultats["fg"] = "green"
bouton_resultats.pack()

# Zone  dans  la  fenêtre  principale  où  on  écrit  un  texte
quitlabel = tk.Label(racine)
quitlabel.pack()

# Bouton Quitter
bouton_quit = tk.Button(racine , text=" Quitter", command=racine.quit)
bouton_quit["fg"] = "red"
bouton_quit.pack()

racine.mainloop()

print("\nJob done!\n")