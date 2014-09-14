# load libraries
#
import numpy as np
import peakdetect as pkdt
from scipy import *

#
# basic functionality: initialize, start, stop, clear
#
def initialize(shell, **kwargs):
    shell.interact(kwargs.copy())

def clearAll(messages, bifurcationDiagram, phasePortrait, timedomainWaveform, FFT, **kwargs):
    bifurcationDiagram.clear()
    phasePortrait.figure.clear()
    timedomainWaveform.clear()
    FFT.clear()
    messages.clear()

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

#
# run: the simulation
#

n_depvars = 3
f_return = np.zeros((n_depvars), dtype=np.float)
x_temp = np.zeros((n_depvars), dtype=np.float)
k1 = np.zeros((n_depvars), dtype=np.float)
k2 = np.zeros((n_depvars), dtype=np.float)
k3 = np.zeros((n_depvars), dtype=np.float)
k4 = np.zeros((n_depvars), dtype=np.float)

def RK4(f, dt, t, x):
    #global f_return, x_temp, k1, k2, k3, k4
    f(t, x, f_return)
    k1[:] = dt*f_return 
    x_temp[:] = x + k1/2. 
    f(t + dt/2., x_temp, f_return) 
    k2[:] = dt*f_return 
    x_temp[:] = x + k2/2. 
    f(t + dt/2., x_temp, f_return)
    k3[:] = dt*f_return  
    x_temp[:] = x + k3 
    f(t + dt, x_temp, f_return) 
    k4[:] = dt*f_return  
    x += (k1+2.*(k2+k3)+k4)/6.

def run(prob_t_start, prob_t_end, prob_t_delta, 
        prob_RVstart, prob_RVend, prob_RVdelta, fft_threshold,
        prob_C1, prob_C2, prob_C3,
        prob_R1, prob_R2, prob_R5, prob_R6, prob_R7, prob_R8, prob_R9, prob_R10, 
        prob_Vin, prob_Rin, 
        saveBifurcation, savePhasePortrait, saveWaveform, saveFFT,
        plotBifurcation, plotPhasePortrait, plotWaveform, plotFFT,
        bifurcationDiagram, phasePortrait, timedomainWaveform, FFT, messages, stop, **kwargs):

                
    C1 = prob_C1.value*1.0e-6     #Convert to Farads
    C2 = prob_C2.value*1.0e-6     #Convert to Farads
    C3 = prob_C3.value*1.0e-6     #Convert to Farads
    R1 = prob_R1.value*1.0e3      #Convert to Ohms
    R2 = prob_R2.value*1.0e3      #Convert to Ohms
    R5 = prob_R5.value*1.0e3      #Convert to Ohms
    R6 = prob_R6.value*1.0e3      #Convert to Ohms
    R7 = prob_R7.value*1.0e3      #Convert to Ohms
    R8 = prob_R8.value*1.0e3      #Convert to Ohms
    R9 = prob_R9.value*1.0e3      #Convert to Ohms
    R10 = prob_R10.value*1.0e3    #Convert to Ohms
    Rin = prob_Rin.value*1.0e3    #Convert to Ohms
    Vin = prob_Vin.value
    
    Lambda = 1/(C3*R10)
    A1 = -1.0    #This gets redefined inside the loop below   
    A2 = -R7/(C1*C2*R5*R6*R8)
    A3 = +R7/(C1*C2*C3*R5*R6*R9*R10)
    A4 = -Vin/(C1*C2*C3*Rin*R6*R10)
    
    def D(x):
        return -(R2/R1)*np.min((x, 0.0))

    def F(t, x, f_return): 
        f_return[0] = x[1]      #x
        f_return[1] = x[2]      #x''
        f_return[2] = A1*x[2] + A2*x[1] + A3*D(x[0]) + A4      #x'''
    
    def PowerSpectrum(f):
        return (f*f.conjugate()).real/len(f)**2

    RVstart = float(prob_RVstart.value)*1.0e3      #Convert to Ohms
    RVend = float(prob_RVend.value)*1.0e3      #Convert to Ohms
    RVdelta = float(prob_RVdelta.value)*1.0e3      #Convert to Ohms


    
    #Output files
    bifurcationDiagram_outfile = open('bifurcationDiagram.txt', 'w')
    phasePortrait_outfile = open('phasePortrait.txt', 'w')
    timedomainWaveform_outfile = open('timedomainWaveform.txt', 'w')
    FFT_outfile = open('fft.txt', 'w')

    
    bifurcationDiagram.clear()
    #phasePortrait.clear()
    timedomainWaveform.clear()
    messages.clear()
    
    #Set up plotting objects
    bifurcationDiagram.set_plot_properties(
        title='Bifurcation Diagram',
        x_label='Rv [kOhms]',
        y_label='x_maxima [V]',
        x_limits = (RVstart*1e-3, RVend*1e-3),
        aspect_ratio='auto')
    # For 2d plot of phase diagram only. 
    #phasePortrait.set_plot_properties(
    #    title='Phase Portrait',
    #    x_label='x [V]',
    #    y_label='x_dot [V]',
    #    aspect_ratio='auto')
    timedomainWaveform.set_plot_properties(
        title='Time-Domain Waveform',
        x_label='t (s)',
        y_label='x [V]',
        aspect_ratio='auto')  
    FFT.set_plot_properties(
        title='FFT',
        x_label='f [Hz]',
        y_label='Power',
        aspect_ratio='auto')
    
    bifurcationDiagram.new_curve(key='bif', memory='growable', length=10000, animated=False,
                marker_style='.', line_color='black', line_style='')
    #The next two lines are for a 2d plot of phasePortrait only.
    #phasePortrait.new_curve(key='phase', memory='growable', length=10000, animated=False,
    #            line_color='red', line_width=1., line_style='-')
    timedomainWaveform.new_curve(key='tdw', memory='growable', length=10000, animated=False,
                line_color='green', marker_style='.', line_style='')
    timedomainWaveform.new_curve(key='tdw_max', memory='growable', length=10000, animated=False,
                marker_style='o', line_color='red', line_style='')
    FFT.new_curve(key='fft', memory='growable', length=10000, animated=False,
                line_color='black', line_width=1., line_style='-')
    
    
    # initial dependent variable values [x(0), x'(0), x''(0)]
    x = np.array([0.0, 0.0, 0.0])
    
    messages.write('Starting simulation...\n')
    for RV in np.arange(RVstart, RVend*(1.001), RVdelta):
        if stop.value: break
        
        messages.write('RV = %g\n' % (RV/1.0e3))

        t_delta = prob_t_delta.value
        t_start = prob_t_start.value
        t_end = prob_t_end.value
        waveformSize = np.int(t_end/t_delta)
        phasePortrait_points = np.zeros((waveformSize, 3), dtype='float')
        timedomainWaveform_points = np.zeros((waveformSize, 2), dtype='float')
        A1 = -1.0/(C1*RV)
        t = 0
        t_index = 0
        threshold = fft_threshold.value
        for i in range(waveformSize): 
            phasePortrait_points[i] = (x[0], -x[1]/Lambda, -x[2]/Lambda**2)
            timedomainWaveform_points[i] = (t, x[0])
            RK4(F, t_delta, t, x)
            t += t_delta
        #phasePortrait.clear()
        #timedomainWaveform.clear()
        #FFT.clear()
        #phasePortrait.figure.clear()
        phasePortrait_axes = phasePortrait.figure.add_subplot(111, projection='3d')
        phasePortrait_axes.set_xlabel('x[v]')
        phasePortrait_axes.set_ylabel('xdot[v]')
        phasePortrait_axes.set_zlabel('xdotdot[v]')
        phasePortrait_axes.set_title('Phase Portrait')
        #Compute FFT. 
        if (plotFFT.value or saveFFT.value):
            ffts = np.delete(np.fft.rfft(timedomainWaveform_points[:,1]), int(waveformSize)/2)
        #Use only positive, non-zero frequencies.
            fs = np.fft.fftfreq(waveformSize, t_delta)[1:waveformSize/2]
            ps = PowerSpectrum(ffts[1:])
        #Cut off plotting frequencies with no significant amplitude in f-space.
            fft_indices = np.where(ps > threshold)[0] 

        t_index = np.where(timedomainWaveform_points[:,0]>t_start)[0][0]
        if plotPhasePortrait.value:
            X = phasePortrait_points[t_index:,0]
            Y = phasePortrait_points[t_index:,1]
            Z = phasePortrait_points[t_index:,2]
            #3d plot in phase-space
            phasePortrait_axes.plot(xs=X, ys=Y, zs=Z, color='red')
            #3d plot is flattened onto bottom-most x-y plane of the graph
            phasePortrait_axes.plot(xs=X, ys=Y, zs=np.amin(Z)*np.ones(len(Z)), color='blue', linestyle='dotted', linewidth=0.5)
            phasePortrait.canvas.draw()
            #the next line is for a 2d plot only
            #phasePortrait.set_data('phase', np.column_stack((phasePortrait_points[t_index:,0], phasePortrait_points[t_index:,1])))
        if plotWaveform.value:
            timedomainWaveform.set_data('tdw', np.column_stack((timedomainWaveform_points[t_index:,0],timedomainWaveform_points[t_index:,1])))
        if plotFFT.value:
            FFT.set_data('fft', np.column_stack((fs[fft_indices[0]-20:fft_indices[-1]+20], ps[fft_indices[0]-20:fft_indices[-1]+20])))
        window_length = 10
        sensitivity = 0
        
        # Skip the initial waveform point, i.e., (0, 0)
        _max, _min = pkdt.peakdetect(timedomainWaveform_points[:,1], timedomainWaveform_points[:,0], window_length, sensitivity)
        #messages.write('I executed peakdetect\n')
        #messages.write('I found peaks = %g\n' % size(_max))
        #for p in _max :
        #    messages.write('max point: (%g, ' % p[0])
        #    messages.write('%g)\n' % p[1])
        #Get rid of maxima that occur at endpoints
        if(_max[0][0] == timedomainWaveform_points[0,0]):
            _max = _max[1:]
        if(_max[-1][0] == t):
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

        if savePhasePortrait.value:
            phasePortrait_outfile.write('RV = {}\n'.format(RV/1.0e3))
            for p in phasePortrait_points:
                phasePortrait_outfile.write('{}\t{}\t{}\n'.format(p[0], p[1], p[2]))

        if saveWaveform.value:
            timedomainWaveform_outfile.write('RV = {}\n'.format(RV/1.0e3))
            for p in timedomainWaveform_points:
                timedomainWaveform_outfile.write('{}\t{}\n'.format(p[0], p[1]))
        
        if saveFFT.value:
            FFT_outfile.write('RV = %f\n' % (RV/1.0e3))
            for i in range(len(fs)):
                FFT_outfile.write('%f\t%f\n' % (fs[i], ps[i]))
        
        # Bifurcation plot
        bifurcationDiagram_points = np.column_stack((RV*np.ones(size(maxima))/1.0e3, maxima))
        if plotBifurcation.value:
            bifurcationDiagram.append_data('bif', bifurcationDiagram_points)
        if saveBifurcation.value:
            for p in bifurcationDiagram_points:
                bifurcationDiagram_outfile.write('{}\t{}\n'.format(p[0], p[1]))
    
    bifurcationDiagram_outfile.close()
    phasePortrait_outfile.close()
    timedomainWaveform_outfile.close()
    FFT_outfile.close()
        
    # reset the stop button in case it was pushed
    stop.value = False
    messages.write('Done.\n')
