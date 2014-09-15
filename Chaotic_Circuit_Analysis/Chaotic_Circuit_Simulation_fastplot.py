#Author: Connor McPartland 
#University of Pittsburgh
#2013


#load libraries
#
import numpy as np
import peakdetect as pkdt
from scipy import *

#
# basic functionality: initialize, start, stop, clear
#
def initialize(shell, **kwargs):
    shell.interact(kwargs.copy())

def clearAll(messages, bifurcationDiagram, phasePortrait, poincarePlot, timedomainWaveform, FFT, **kwargs):
    #bifurcationDiagram.clear()
    phasePortrait.figure.clear()
    timedomainWaveform.clear()
    FFT.clear()
    messages.clear()
    poincarePlot.clear()

def clear_bifurcationDiagram(bifurcationDiagram, messages, **kwargs):
    bifurcationDiagram.clear()
    messages.clear()

def clear_phasePortrait(phasePortrait, messages, **kwargs):
    #phasePortrait.clear()
    phasePortrait.figure.clear()
    messages.clear()

def clear_timedomainWaveform(timedomainWaveform, messages, **kwargs):
    timedomainWaveform.clear()
    messages.clear()

def clear_FFT(FFT, messages, **kwargs):
    FFT.clear()
    messages.clear()

def plot_RV(poincarePlot, plotPhasePortrait, plotWaveform, plotFFT, open_file, params_file, FFTthreshold_box,
        phasePortrait, timedomainWaveform, FFT, messages, stop, **kwargs):
    p_file = open(params_file.value, 'r')
    C1 = float(p_file.readline().split(':')[1])
    C2 = float(p_file.readline().split(':')[1])
    C3 = float(p_file.readline().split(':')[1])
    R1 = float(p_file.readline().split(':')[1])
    R2 = float(p_file.readline().split(':')[1])
    R5 = float(p_file.readline().split(':')[1])
    R6 = float(p_file.readline().split(':')[1])
    R7 = float(p_file.readline().split(':')[1])
    R8 = float(p_file.readline().split(':')[1])
    R9 = float(p_file.readline().split(':')[1])
    R10 = float(p_file.readline().split(':')[1])
    Rin = float(p_file.readline().split(':')[1])
    Vin = float(p_file.readline().split(':')[1])
    RVstart = float(p_file.readline().split(':')[1])
    RVend = float(p_file.readline().split(':')[1])
    RVdelta = float(p_file.readline().split(':')[1])
    #DPC_start = float(p_file.readline().split(':')[1])
    #DPC_end = float(p_file.readline().split(':')[1])
    t_start = float(p_file.readline().split(':')[1])
    t_end = float(p_file.readline().split(':')[1])
    t_delta = float(p_file.readline().split(':')[1])
    Lambda =  1/(C3*R10)
    
    def F(t, x, f_return): 
        f_return[0] = x[1]      #x
        f_return[1] = x[2]      #x''
        f_return[2] = A1*x[2] + A2*x[1] + A3*D(x[0]) + A4      #x'''
    
    def PowerSpectrum(f):
        return (f*f.conjugate()).real/len(f)**2

    #phasePortrait.clear()
    timedomainWaveform.clear()
    FFT.clear()
    messages.clear()
    
    
    data_file = open(open_file.value, 'r')
    messages.write('Starting simulation...\n')
    data_file.readline()
    data = data_file.readlines()
    data.pop(0)
    data_array = np.zeros((len(data),4))
    for i in range(len(data)):
        data_array[i] = tuple(data[i].split('\t '))
    #ts = len(data_array)
    #t_start = data_array[int(ts*80/100),0]
    #t_end = data_array[ts-1,0]
    #t_delta = data_array[1,0]-data_array[0,0]
    waveformSize = np.int(t_end/t_delta)
    # data_array[0] = t, data_array[1] = x, data_array[2] = xdot, data_array[3] = xdotdot
    t_index = np.where(data_array[:,0]>t_start)[0][0]
    
    RV = float(open_file.value.split('_')[-1].strip('.txt'))
    #global data_type
    threshold = float(FFTthreshold_box.value)
        
    #bifurcationDiagram.set_plot_properties(
    #    title='Bifurcation Diagram',
    #    x_label='Rv [kOhms]',
    #    y_label='x_maxima [V]',
    #    x_limits=(RVstart*1e-3-.1, 1e-3*RVend+1),
    #    aspect_ratio='auto')

    # Set up plotting objects
    # Initialize a Poincare Plot
    poincarePlot.set_plot_properties(
        title='Poincare Plot',
        x_label='x [V]',
        y_label='x_dot [V]',
        aspect_ratio='auto')

    # Initialize a time-domain plot
    timedomainWaveform.set_plot_properties(
        title='Time-Domain Waveform',
        x_label='t (s)',
        y_label='x [V]',
        aspect_ratio='auto')  
    
    # Initialize a FFT plot
    FFT.set_plot_properties(
        title='FFT',
        x_label='f [Hz]',
        y_label='Power',
        y_scale='log',
        aspect_ratio='auto')
    
    #bifurcationDiagram.new_curve(key='bif', memory='growable', length=10000, animated=False,
    #            marker_style='.', line_color='black', line_style='')
    poincarePlot.new_curve(key='pp', memory='growable', length=10000, animated=False,
                line_color='red', line_style='', marker_style='.')
    timedomainWaveform.new_curve(key='tdw', memory='growable', length=10000, animated=False,
                line_color='green', marker_style='.', line_style='')
    timedomainWaveform.new_curve(key='tdw_max', memory='growable', length=10000, animated=False,
                marker_style='o', line_color='red', line_style='')
    FFT.new_curve(key='fft', memory='growable', length=10000, animated=False,
                line_color='black', line_width=1., line_style='-')
                
    if plotPhasePortrait.value:
        phasePortrait_axes = phasePortrait.figure.add_subplot(111, projection='3d')
        phasePortrait_axes.set_xlabel('x[v]')
        phasePortrait_axes.set_ylabel('xdot[v]')
        phasePortrait_axes.set_zlabel('xdotdot[v]')
        phasePortrait_axes.set_title('Phase Portrait')
        X = data_array[t_index:,1]
        Y = -1*data_array[t_index:,2]/Lambda
        Z = -1*data_array[t_index:,3]/Lambda**2
        poincare_points = []
        for i in range(len(Z)-1):
            if Z[i+1] > 0 and Z[i] < 0:
                poincare_points.append((X[i+1], Y[i+1]))
        points = np.array(poincare_points)
        poincarePlot.set_data('pp', np.column_stack((points[:,0], points[:,1])))
        poincarePlot.set_plot_properties(x_limits=(np.amin(X)-.1, np.amax(X)+.1), y_limits=(np.amin(Y)-.1, np.amax(Y)+.1))
        #3d plot in phase-space - illustrates the fact that even though the function looks like it returns to its initial value on a 2d-map,
        # in reality it is 'traveling' in a 3rd dimension, and must travel in the negative direction of the 3rd dimension to truly return
        # to initial value
        phasePortrait_axes.plot(xs=X, ys=Y, zs=Z, color='red')
        #3d plot is flattened onto bottom-most x-y plane of the graph
        phasePortrait_axes.plot(xs=X, ys=Y, zs=np.amin(Z)*np.ones(len(Z)), color='blue', linestyle='dotted', linewidth=0.5)
        #phasePortrait_axes.plot(xs=points[:,0], ys=points[:,1], zs=np.amin(Z)*np.ones(len(points)), color='blue', linestyle='dotted', linewidth=0.5)
        phasePortrait.canvas.draw()

    #Compute FFT. 
    if (plotFFT.value):
        temp = np.fft.rfft(data_array[:,1])
        ffts = np.delete(temp, len(temp)-1)
    #Use only positive, non-zero frequencies.
        fs = np.fft.fftfreq(waveformSize, t_delta)[1:waveformSize/2]
        ps = PowerSpectrum(ffts[1:])
    #Cut off plotting frequencies with no significant amplitude in f-space.
        fft_indices = np.where(ps > threshold)[0] 

    
    if plotWaveform.value:
        timedomainWaveform.set_data('tdw', np.column_stack((data_array[t_index:,0],data_array[t_index:,1])))
    if plotFFT.value:
        FFT.set_data('fft', np.column_stack((fs[fft_indices[0]-20:fft_indices[-1]+20], ps[fft_indices[0]-20:fft_indices[-1]+20])))
    window_length = 10
    sensitivity = 0
        
    # Skip the initial waveform point, i.e., (0, 0)
    _max, _min = pkdt.peakdetect(data_array[:,1], data_array[:,0], window_length, sensitivity)
    #messages.write('I executed peakdetect\n')
    #messages.write('I found peaks = %g\n' % size(_max))
    #for p in _max :
    #    messages.write('max point: (%g, ' % p[0])
    #    messages.write('%g)\n' % p[1])
    #Get rid of maxima that occur at endpoints
    if(_max[0][0] == data_array[0,0]):
        _max = _max[1:]
    if(_max[-1][0] == t_end):
        _max = _max[0:-1]
    times = [p[0] for p in _max]
    maxima = [p[1] for p in _max]
    times_index = 0
    
    for i in range(len(times)):
        if times[i] > t_start:
            times_index = i
            break
    if plotWaveform.value:
        timedomainWaveform.set_data('tdw_max', np.column_stack((times[times_index:], maxima[times_index:])))
    messages.write('Data plotting done.\n')


# Bifurcation plot      
def plot_bifurcation(bifurcationDiagram, bifurcation_file, messages, **kwargs):
    
    #bifurcationDiagram_points = np.column_stack((RV*np.ones(size(maxima))/1.0e3, maxima))
    bifurcationDiagram.figure.clear()
    b_file = open(bifurcation_file.value, 'r')
    first_line = b_file.readline().split('\t')
    data_type = first_line[0]
        
    b_data = b_file.readlines()
    b_points = np.zeros((len(b_data), 2))
    for i in range(len(b_data)):
        b_points[i] = tuple(b_data[i].split('\t'))
    if data_type == 'DPC':
        f = open("saved_RV_and_DPC_Thu Apr 11 18-06-59 2013.txt", 'r')
        f.readline()
        f_data = f.readlines()
        dpcs = np.zeros(len(f_data))
        Rvs = np.zeros(len(f_data))
        Rv_values = np.zeros(len(b_points))
        for i in range(len(f_data)):
            temp = f_data[i].split('\t')
            dpcs[i] = temp[0]
            Rvs[i] = temp[1]
        for i in range(len(b_points)):
            Rv_values[i] = Rvs[np.where(b_points[i][0] == dpcs)[0]]
        start_index = np.where(dpcs == b_points[:,0][0])[0]
        end_index = np.where(dpcs == b_points[:,0][-1])[0]
        b_points[:,0] = Rv_values/1000.
    bifurcation_axes = bifurcationDiagram.figure.add_subplot(111)
    bifurcation_axes.set_xlabel('Rv[kOhm]')
    bifurcation_axes.set_xlimits=(b_points[:,0][0], b_points[:,0][-1])
    bifurcation_axes.set_ylabel('x_maxima [v]')
    bifurcation_axes.set_title('Bifurcation Diagram')
    bifurcation_axes.scatter(b_points[:,0], b_points[:,1], color='black', marker='.')
    bifurcationDiagram.canvas.draw()

    #bifurcationDiagram.set_plot_properties(
    #    title='Bifurcation Diagram',
    #    x_label='DPC',
    #    y_label='x_maxima [V]',
    #    x_limits=(b_points[0,0]-.1, b_points[len(b_points)-1,0]+.1),
    #    aspect_ratio='auto')
    #bifurcationDiagram.new_curve(key='bif', memory='growable', length=10000, animated=False,
    #            marker_style='.', line_color='black', line_style='')
    #bifurcationDiagram.set_data('bif', b_points)
    b_file.close()
    messages.write('Bifurcation plot done.\n')
