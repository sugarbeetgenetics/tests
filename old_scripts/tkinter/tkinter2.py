from tkinter import *
from tkinter import ttk

root = Tk()
root.title("Test GUI")
frame = Frame(root)

Label(frame, text="Buttons").pack()

Button(frame, text="B1").pack(side=LEFT, fill=Y)
Button(frame, text="B2").pack(side=RIGHT, fill=X)
Button(frame, text="B3").pack(side=TOP, fill=X)
Button(frame, text="B4").pack(side=BOTTOM, fill=Y)

frame.pack()

root.mainloop()