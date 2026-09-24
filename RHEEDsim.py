"""Run without arguments for the desktop GUI; --help lists batch options."""
import argparse
from pathlib import Path
import numpy as np
from core import Surface, lattices, simulate, local_LEEDgen

HERE = Path(__file__).resolve().parent


def colormap():
    from matplotlib.colors import ListedColormap
    return ListedColormap(np.load(HERE / 'RHEEDcmap.npy'), name='RHEEDsim')


def plot_data(ax, kind, data):
    ax.clear()
    if kind == 'line':
        ax.plot(data[:, 0], data[:, 1])
        ax.set(xlabel='Reciprocal distance (Å⁻¹)', ylabel='Intensity', title='RHEED line profile')
    else:
        x, y, density = data
        ax.pcolormesh(x, y, density, cmap=colormap(), shading='gouraud')
        ax.set_aspect('equal')
        ax.set(xlabel='Reciprocal distance (Å⁻¹)', ylabel='Reciprocal distance (Å⁻¹)',
               title='RHEED streaks' if kind == 'streaks' else 'LEED')


def launch_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    class App:
        def __init__(self, root):
            self.root = root
            self.surface = None
            self.results = {}
            root.title('RHEEDsim — Python port of Wang & Smith (2011)')
            root.geometry('1250x900')
            controls = ttk.Frame(root, padding=8)
            controls.pack(fill='x')
            self.fields = {}
            self.status = tk.StringVar(value='Load a surface model to begin.')
            ttk.Button(controls, text='Load model', command=self.load).grid(row=0, column=0, padx=4)
            self.info = tk.StringVar()
            ttk.Label(controls, textvariable=self.info).grid(row=0, column=1, columnspan=12, sticky='w')
            for row, entries in enumerate([
                [('X repeats', '5'), ('Y repeats', '5')],
                [('Beam angle (°)', '0'), ('RHEED sigma', '')],
                [('LEED radius', '6'), ('LEED sigma', ''), ('Pixels', '300')],
            ], start=1):
                for i, (label, default) in enumerate(entries):
                    ttk.Label(controls, text=label).grid(row=row, column=2*i, sticky='e')
                    var = tk.StringVar(value=default)
                    self.fields[label] = var
                    ttk.Entry(controls, textvariable=var, width=10).grid(row=row, column=2*i+1, padx=5, pady=3)
            for row, label, callback in [(1, 'Plot lattices', self.plot_lattices), (2, 'Line profile', self.line), (3, 'LEED', self.leed)]:
                ttk.Button(controls, text=label, command=lambda f=callback: self.guard(f)).grid(row=row, column=6, padx=4)
            ttk.Button(controls, text='Streaks', command=lambda: self.guard(self.streaks)).grid(row=2, column=7)
            self.figure = Figure(figsize=(11, 7), constrained_layout=True)
            grid = self.figure.add_gridspec(2, 6)
            self.axes = {'real': self.figure.add_subplot(grid[0, :2]),
                         'reciprocal': self.figure.add_subplot(grid[0, 2:4]),
                         'line': self.figure.add_subplot(grid[0, 4:]),
                         'streaks': self.figure.add_subplot(grid[1, :3]),
                         'leed': self.figure.add_subplot(grid[1, 3:])}
            self.canvas = FigureCanvasTkAgg(self.figure, master=root)
            self.canvas.get_tk_widget().pack(fill='both', expand=True)
            NavigationToolbar2Tk(self.canvas, root)
            exports = ttk.Frame(root)
            exports.pack(fill='x')
            for i, kind in enumerate(('line', 'streaks', 'leed')):
                ttk.Button(exports, text=f'Save {kind} data', command=lambda k=kind: self.guard(lambda: self.save(k))).grid(row=0, column=2*i, padx=4)
                ttk.Button(exports, text=f'{kind.title()}: new window', command=lambda k=kind: self.guard(lambda: self.new_window(k))).grid(row=0, column=2*i+1, padx=4)
            ttk.Label(root, textvariable=self.status, padding=6).pack(fill='x')
            self.progress = ttk.Progressbar(root, maximum=100)
            self.progress.pack(fill='x')

        def guard(self, callback):
            try:
                if self.surface is None:
                    raise ValueError('Load a surface model first.')
                self.root.config(cursor='watch')
                self.root.update_idletasks()
                callback()
                self.canvas.draw_idle()
            except Exception as exc:
                messagebox.showerror('RHEEDsim', str(exc))
            finally:
                self.root.config(cursor='')

        def load(self):
            path = filedialog.askopenfilename(initialdir=HERE / 'examples', filetypes=[('Surface model', '*.txt')])
            if not path:
                return
            try:
                model = Surface.load(path)
                self.surface = model
                self.results.clear()
                self.fields['RHEED sigma'].set('')
                self.fields['LEED sigma'].set('')
                self.info.set(f'{model.title} — {len(model.coordinates)} atoms')
                self.status.set(f'a = {model.cell.tolist()}   b = {np.round(model.reciprocal, 5).tolist()}')
                for ax in self.axes.values():
                    ax.clear()
                self.plot_lattices()
                self.canvas.draw_idle()
            except Exception as exc:
                messagebox.showerror('Cannot load model', str(exc))

        def number(self, name, default=None, integer=False):
            value = self.fields[name].get().strip()
            if not value:
                if default is None:
                    return None
                value = default
            result = float(value)
            if not np.isfinite(result) or (integer and result != int(result)):
                raise ValueError(f'{name} must be a finite {"integer" if integer else "number"}.')
            return int(result) if integer else result

        def plot_lattices(self):
            real, reciprocal, basis = lattices(self.surface, self.number('X repeats', 5, True), self.number('Y repeats', 5, True))
            for key, points, title in [('real', real, 'Real lattice (Å)'), ('reciprocal', reciprocal, 'Reciprocal lattice (Å⁻¹)')]:
                ax = self.axes[key]
                ax.clear()
                if key == 'real':
                    ax.plot(basis[:, 0], basis[:, 1], 'g+')
                ax.plot(points[:, 0], points[:, 1], 'o', color='blue', markerfacecolor='blue' if key == 'real' else 'red', markersize=4)
                ax.set_aspect('equal')
                ax.set_title(title)

        def calculate(self):
            angle = self.number('Beam angle (°)', 0)
            width = self.number('RHEED sigma')
            s, im, fine, streaks = simulate(self.surface, angle, width)
            if width is None:
                self.fields['RHEED sigma'].set(f'{s[2]/5:.16g}')
            self.results['line'] = fine
            self.results['streaks'] = (fine[:, 0], np.linspace(-10, 0, 400), streaks)
            self.status.set(f'Perpendicular vector: {int(s[0])} b1 + {int(s[1])} b2; spacing = {s[2]:.8g} Å⁻¹')
            self.progress['value'] = 100

        def line(self):
            self.calculate()
            plot_data(self.axes['line'], 'line', self.results['line'])

        def streaks(self):
            self.line()
            plot_data(self.axes['streaks'], 'streaks', self.results['streaks'])

        def leed(self):
            radius = self.number('LEED radius', 6)
            width = self.number('LEED sigma', np.linalg.norm(self.surface.reciprocal[0])/5)
            pixels = self.number('Pixels', 300, True)
            self.fields['LEED sigma'].set(f'{width:.16g}')
            def progress(fraction):
                self.progress['value'] = fraction*100
                self.root.update_idletasks()
            _, density = local_LEEDgen(radius, width, pixels, self.surface.reciprocal, self.surface.coordinates, progress)
            axis = np.linspace(-radius, radius, pixels)
            self.results['leed'] = (axis, axis, density)
            plot_data(self.axes['leed'], 'leed', self.results['leed'])
            self.status.set(f'LEED map: {pixels} × {pixels}')

        def result(self, kind):
            if kind not in self.results:
                raise ValueError(f'Compute {kind} first.')
            return self.results[kind]

        def save(self, kind):
            data = self.result(kind)
            path = filedialog.asksaveasfilename(defaultextension='.txt', initialfile=f'{kind}.txt')
            if path:
                np.savetxt(path, data if kind == 'line' else data[2], fmt='%.18e')
                self.status.set(f'Saved {path}')

        def new_window(self, kind):
            data = self.result(kind)
            window = tk.Toplevel(self.root)
            window.title(kind.title())
            figure = Figure(figsize=(8, 6), constrained_layout=True)
            plot_data(figure.add_subplot(), kind, data)
            canvas = FigureCanvasTkAgg(figure, master=window)
            canvas.get_tk_widget().pack(fill='both', expand=True)
            NavigationToolbar2Tk(canvas, window)
            canvas.draw()

    root = tk.Tk()
    App(root)
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, help='Surface file; enables batch mode')
    parser.add_argument('--angle', type=float, default=0)
    parser.add_argument('--width', type=float, help='RHEED Gaussian sigma; default spacing/5')
    parser.add_argument('--leed', action='store_true', help='Also calculate LEED')
    parser.add_argument('--radius', type=float, default=6)
    parser.add_argument('--leed-width', type=float)
    parser.add_argument('--pixels', type=int, default=300)
    parser.add_argument('--output', type=Path, default=HERE / 'output')
    args = parser.parse_args()
    if args.input is None:
        launch_gui()
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    surface = Surface.load(args.input)
    s, peaks, fine, streaks = simulate(surface, args.angle, args.width)
    args.output.mkdir(parents=True, exist_ok=True)
    datasets = {'line': fine, 'streaks': (fine[:, 0], np.linspace(-10, 0, 400), streaks)}
    np.savetxt(args.output / 'peaks.txt', peaks)
    if args.leed:
        width = args.leed_width if args.leed_width is not None else np.linalg.norm(surface.reciprocal[0])/5
        leed_peaks, density = local_LEEDgen(args.radius, width, args.pixels, surface.reciprocal, surface.coordinates)
        np.savetxt(args.output / 'leed_peaks.txt', leed_peaks)
        axis = np.linspace(-args.radius, args.radius, args.pixels)
        datasets['leed'] = (axis, axis, density)
    for kind, data in datasets.items():
        np.savetxt(args.output / f'{kind}.txt', data if kind == 'line' else data[2])
        fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
        plot_data(ax, kind, data)
        fig.savefig(args.output / f'{kind}.png', dpi=160)
        plt.close(fig)
    print(f'{surface.title}\nReciprocal vector: {s}\nSaved results to {args.output.resolve()}')


if __name__ == '__main__':
    main()
