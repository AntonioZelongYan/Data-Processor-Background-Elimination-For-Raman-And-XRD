import tkinter as tk
from tkinter import filedialog
from tkinter import Scale
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from functools import partial
import pandas as pd
from scipy.signal import find_peaks

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        #GUI setting
        self.master.geometry("1500x800")
        self.master.minsize(1500, 800)
        self.my_title = "Data Processer"
        self.master.title(self.my_title)
        #Initial global variables
        self.x_data = []
        self.y_data = []
        self.baseline_x_data = []
        self.baseline_y_data = []
        self.preview_x_data = []
        self.preview_y_data = []
        self.peakx = []
        self.peaky = []
        self.upper = tk.StringVar(value='0')
        self.lower = tk.StringVar(value='0')
        self.upperData = 0
        self.lowerData = 0
        #Initial global parameters setting
        self.p_value = 0.0001
        self.m_value = 100
        self.MA_value = 0
        self.peak_value = 0
        #Figure initialization
        self.dataTitle = ''
        self.fig1 ,self.ax1 = plt.subplots(figsize=(8, 6))
        self.fig2,self.ax2 = plt.subplots(figsize=(8, 6))
        #GUI setting
        self.create_menu() 
        self.create_widget()
    
    #File opening
    def menu_open_clicked(self, event=None):
        file_path = tk.filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
        self.set_data(file_path)
    #App closing
    def menu_quit_clicked(self):
        self.master.destroy() 
    #Menu setting
    def create_menu(self):
        self.menu_bar = tk.Menu(self)
        #Menu setting
        self.file_menu = tk.Menu(self.menu_bar, tearoff = tk.OFF)
        self.menu_bar.add_cascade(label="File", menu=self.file_menu)
        #Menu function
        self.file_menu.add_command(label="Open", command = self.menu_open_clicked, accelerator="Ctrl+O")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command = self.menu_quit_clicked)
        self.menu_bar.bind_all("<Control-o>", self.menu_open_clicked)

        self.master.config(menu=self.menu_bar)
    #The Alternating Least Square algorithm
    def baseline_als_optimized(self,y, p, lam, niter=10):
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
        w = np.ones(L)
        W = sparse.spdiags(w, 0, L, L)
        for i in range(niter):
            W.setdiag(w) # Do not create a new matrix, just update diagonal values
            Z = W + D
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return list(z)
    #Data loading
    def set_data(self, filename):
        self.x_data = []
        self.y_data = []
        self.baseline_x_data = []
        self.baseline_y_data = []
        self.preview_x_data = []
        self.preview_y_data = []
        self.peakx = []
        self.peaky = []
        self.upper.set('0')
        self.lower.set('0')
        self.upperData = 0 
        self.lowerData = 0 
        temp = filename.split('/')
        temp = temp[-1][:-4]
        self.dataTitle = temp
        if not filename:
            return
        with open(filename, 'r') as file:
            for line in file:
                x, y = map(float, line.strip().split())
                self.x_data.append(x)
                self.y_data.append(y)
                self.baseline_x_data.append(x)
                self.baseline_y_data.append(y)
        self.baseline_y_data = self.baseline_als_optimized(self.baseline_y_data,self.p_value,self.m_value)
        self.upper.set(str(max(self.x_data)))
        self.upperData = max(self.x_data)
        self.lower.set(str(min(self.x_data)))
        self.lowerData = min(self.x_data)
        self.update_plot(1)
    #Moving average function
    def moving_average(self,data,parameter):
        output = []
        for i in range(len(data)):
            if i == 0:
                output.append(data[0])
            elif i < parameter:
                temp = data[:i]
                output.append(np.mean(temp))
            elif i > len(data) - parameter - 1:
                temp = data[i:]
                output.append(np.mean(temp))
            else:
                temp = data[i-parameter:i+parameter]
                output.append(np.mean(temp))
        return output
    #Plot updating
    def update_plot(self,val):
        self.p_value,  self.m_value , self.MA_value ,self.peak_value = self.slider_p.get(), self.slider_n.get() , self.slider_MA.get(), self.slider_peak.get()
        upper, lower = float(self.upper_range.get()), float(self.lower_range.get())
        #Range check
        if upper > self.upperData or upper < lower:
            upper = self.upperData
        if lower < self.lowerData or lower > upper:
            lower = self.lowerData
        px_data = []
        py_data = []
        new_base_x_data = []
        temp = []
        for i in range(len(self.x_data)):
            if self.x_data[i] > lower and self.x_data[i] < upper:
                px_data.append(self.x_data[i])
                py_data.append(self.y_data[i])
                new_base_x_data.append(self.baseline_x_data[i])
                temp.append(self.baseline_y_data[i])
        #Origin figure ploting
        new_base_y_data = self.baseline_als_optimized(temp,self.p_value,self.m_value)
        self.ax1.clear()
        self.ax1.plot(px_data, py_data, color='black')
        self.ax1.plot(new_base_x_data, new_base_y_data, color='red')
        self.ax1.set_xlabel('Wavelength')
        self.ax1.set_ylabel('Intensity')
        self.ax1.set_title(self.dataTitle)
        self.ax1.relim()  # Recalculate the data limits
        self.ax1.autoscale_view()  # Auto-adjust the axis limits
        self.canvas.draw()
        #Preview figure ploting
        temp_x_data = []
        temp_y_data = []
        for i in range(len(self.x_data)):
            if self.x_data[i] > lower and self.x_data[i] < upper:
                temp_x_data.append(self.x_data[i])
                temp_y_data.append(self.y_data[i])
        if self.MA_value == 0:
            temp_y_data = list((np.array(temp_y_data) - np.array(new_base_y_data))/max(np.array(temp_y_data) - np.array(new_base_y_data)))
        else:
            temp_y_data = list((np.array(temp_y_data) - np.array(new_base_y_data)))
            temp_y_data = self.moving_average(temp_y_data,self.MA_value)
            temp_y_data = temp_y_data / max(temp_y_data)
        peaks, _ = find_peaks(temp_y_data,prominence=self.peak_value)
        self.peakx = []
        self.peaky = []
        for i in peaks:
            self.peakx.append(temp_x_data[i])
            self.peaky.append(temp_y_data[i])
        self.ax2.clear()
        self.ax2.plot(temp_x_data, temp_y_data, color='black')
        self.ax2.scatter(self.peakx, self.peaky, "x", c = 'red')
        self.ax2.set_xlabel('Wavelength')
        self.ax2.set_ylabel('Intensity')
        self.ax2.set_title('Preview')
        self.ax2.relim()  # Recalculate the data limits
        self.ax2.autoscale_view()  # Auto-adjust the axis limits
        self.canvasPreview.draw()
        return
    #Reset function
    def reset_parameters(self):
        self.slider_p.set(0.00001)
        self.slider_n.set(10000)
    #txt output function
    def save_txt_data(self):
        if self.x_data and self.y_data:
            file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
            if file_path:
                with open(file_path, 'w') as file:
                    base_y_data = self.baseline_als_optimized(self.y_data,self.p_value,self.m_value)
                    if self.MA_value == 0:
                        temp_y_data = list((np.array(self.y_data) - np.array(base_y_data))/max(np.array(self.y_data) - np.array(base_y_data)))
                    else:
                        temp_y_data = list((np.array(self.y_data) - np.array(base_y_data)))
                        temp_y_data = self.moving_average(temp_y_data,self.MA_value)
                        temp_y_data = temp_y_data / max(temp_y_data)
                    for x, y in zip(self.x_data, temp_y_data):
                        new_x = x
                        new_y = y
                        file.write(f"{new_x} {new_y}\n")
    #csv output function
    def save_csv_data(self):
        if self.x_data and self.y_data:
            file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("Data Files", "*.csv")])
            if file_path:
                base_y_data = self.baseline_als_optimized(self.y_data,self.p_value,self.m_value)
                if self.MA_value == 0:
                    temp_y_data = list((np.array(self.y_data) - np.array(base_y_data))/max(np.array(self.y_data) - np.array(base_y_data)))
                else:
                    temp_y_data = list((np.array(self.y_data) - np.array(base_y_data)))
                    temp_y_data = self.moving_average(temp_y_data,self.MA_value)
                    temp_y_data = temp_y_data / max(temp_y_data)
                    
                df = pd.DataFrame({'x':self.x_data,'y':temp_y_data})
                df.to_csv(file_path, index=False)
    #Peak saving function
    def save_peak(self):
        if self.peakx  and self.peaky:
            file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("Data Files", "*.csv")])
            if file_path:
                df = pd.DataFrame({'Wavelength':self.peakx,'Intensity':self.peaky})
                df.to_csv(file_path, index=False)
    #Widget Info
    def create_widget(self):
        self.app_label = tk.Label(self.master, text="The algorithm is based on Alternating Least Squares and the Find-peak function is built in Scipy", font=("Helvetica", 16))
        self.app_label.pack(side=tk.TOP)
        
        self.gap1 = tk.Label(self.master, text="", font=("Helvetica", 10))
        self.gap1.pack(side=tk.TOP)
        
        
        self.frameSlider = tk.Frame(self.master)
        
        self.label_p = tk.Label(self.frameSlider, text="Parameter p:", font=("Helvetica", 16))
        self.label_p.pack(side=tk.LEFT)
        
        self.slider_p = Scale(self.frameSlider, from_=0.00001, to=0.001, resolution=0.00001, orient='horizontal', command=partial(self.update_plot))
        self.slider_p.set(0.00001)  
        self.slider_p.pack(side=tk.LEFT)
        
        self.reset_button = tk.Button(self.frameSlider, text="Reset", command=self.reset_parameters)
        self.reset_button.pack(side=tk.RIGHT)
        
        self.slider_n = Scale(self.frameSlider, from_=10000, to=1000000, resolution=10000, orient='horizontal', command=partial(self.update_plot))
        self.slider_n.set(10000)  
        self.slider_n.pack(side=tk.RIGHT)
        
        self.label_n = tk.Label(self.frameSlider, text="Parameter λ:", font=("Helvetica", 16))
        self.label_n.pack(side=tk.RIGHT)
        
        self.frameSlider.pack(side=tk.TOP)
        
        self.gap2 = tk.Label(self.master, text="", font=("Helvetica", 10))
        self.gap2.pack(side=tk.TOP)
        
        self.frameSave = tk.Frame(self.master)
        
        self.save_txt_button = tk.Button(self.frameSave, text="Save txt Data", command=self.save_txt_data)
        self.save_txt_button.pack(side=tk.LEFT)
        
        self.save_csv_button = tk.Button(self.frameSave, text="Save csv Data", command=self.save_csv_data)
        self.save_csv_button.pack(side=tk.RIGHT)
        
        self.frameSave.pack(side=tk.TOP)
        
        self.frameFunction = tk.Frame(self.master)
        
        self.label_peak = tk.Label(self.frameFunction, text="Find-Peak Prominence:", font=("Helvetica", 16))
        self.label_peak.pack(side=tk.LEFT)
        
        self.slider_peak = Scale(self.frameFunction, from_=0, to=1, resolution=0.01, orient='horizontal', command=partial(self.update_plot))
        self.slider_peak.set(0)  
        self.slider_peak.pack(side=tk.LEFT)
        
        self.save_peak_button = tk.Button(self.frameFunction, text="Save Peaks", command=self.save_peak)
        self.save_peak_button.pack(side=tk.LEFT)
        
        self.gap3 = tk.Label(self.frameFunction, text="       ", font=("Helvetica", 10))
        self.gap3.pack(side=tk.LEFT)
        
        self.slider_MA = Scale(self.frameFunction, from_=0, to=15, resolution=1, orient='horizontal', command=partial(self.update_plot))
        self.slider_MA.set(0)  
        self.slider_MA.pack(side=tk.RIGHT)
        
        self.label_MA = tk.Label(self.frameFunction, text="Moving Average Parameter:", font=("Helvetica", 16))
        self.label_MA.pack(side=tk.RIGHT)
        
        self.gap4 = tk.Label(self.frameFunction, text="       ", font=("Helvetica", 10))
        self.gap4.pack(side=tk.RIGHT)
        
        self.updateRange = tk.Button(self.frameFunction, text="Update", command=partial(self.update_plot,1))
        self.updateRange.pack(side=tk.RIGHT)
        
        
        self.upper_range = tk.Entry(self.frameFunction, width=8, textvariable=self.upper)
 
        self.upper_range.pack(side=tk.RIGHT)
        
        self.label_upper = tk.Label(self.frameFunction, text="Upper range:", font=("Helvetica", 16))
        self.label_upper.pack(side=tk.RIGHT)
        

        self.lower_range = tk.Entry(self.frameFunction, width=8, textvariable=self.lower)
 
        self.lower_range.pack(side=tk.RIGHT)
        
        self.label_lower = tk.Label(self.frameFunction, text="Lower range:", font=("Helvetica", 16))
        self.label_lower.pack(side=tk.RIGHT)
        
        self.frameFunction.pack(side=tk.TOP)
        
        
        self.frameView = tk.Frame(self.master)
        
        self.canvas = FigureCanvasTkAgg(self.fig1, master=self.frameView)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(expand=True,side=tk.LEFT,  fill=tk.BOTH)
        
        self.canvasPreview = FigureCanvasTkAgg(self.fig2, master=self.frameView)
        self.canvasPreview.draw()
        self.canvasPreview.get_tk_widget().pack(expand=True,side=tk.RIGHT,  fill=tk.BOTH)
        
        self.frameView.pack(expand=True,  fill=tk.BOTH)
        
        self.frameView = tk.Frame(self.master)


if __name__ == "__main__":
    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()