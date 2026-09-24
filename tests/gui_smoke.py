"""Optional desktop smoke test; requires a display and Tk."""
import sys
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tkinter as tk
from tkinter import ttk
import RHEEDsim


def exercise(root):
    buttons = {}
    def visit(widget):
        for child in widget.winfo_children():
            if isinstance(child, ttk.Button):
                buttons[child.cget('text')] = child
            visit(child)
    visit(root)
    buttons['Load model'].invoke()
    buttons['Line profile'].invoke()
    buttons['Streaks'].invoke()
    buttons['LEED'].invoke()
    buttons['Line: new window'].invoke()
    root.update()
    root.destroy()


def fail(title, message):
    raise AssertionError(f'{title}: {message}')


with patch.object(tk.Tk, 'mainloop', exercise), patch('tkinter.filedialog.askopenfilename', return_value=str(RHEEDsim.HERE/'examples/sqrt3HD.txt')), patch('tkinter.messagebox.showerror', side_effect=fail):
    RHEEDsim.launch_gui()
print('GUI smoke passed: load, lattices, line profile, streaks, LEED, new window.')
