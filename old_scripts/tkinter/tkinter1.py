from tkinter import *
from tkinter import ttk
root = Tk()

root.title("Test GUI")
frame = Frame(root)
labelText = StringVar()
label = Label(frame,  textvariable=labelText)
button = Button(frame, text="Click me")

labelText.set("I am a label")
label.pack()
button.pack()
frame.pack()
root.mainloop()