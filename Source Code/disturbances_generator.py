# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 12:29:40 2018

@author: Rene
"""
import numpy as np
import math
from random import uniform
from random import randint

class Disturbances(object):
    
    def __init__(self, quant, tam, periodo, freq, snr=0):
        self.num = quant;
        self.tam = tam;
        self.periodo = periodo;
        self.freq = freq;
        self.snr = snr;
        self.t = np.linspace(0.0, self.tam/10000.0, self.tam);

    def pure_sine(self): 
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            #Frquency Value
            freq = self.freq + uniform(-self.freq*0.1,self.freq*0.1);
            #Amplitude Reduction
            mult_amp = uniform(0.999, 1.001);
            seno = mult_amp*np.sin(2*np.pi*freq*self.t);
            matrix_r[i] = seno;
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
        return matrix_r;
    
    def interrupt(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            #Define where begins and ends the interruption
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Amplitude Reduction
            mult_amp = uniform(0, 0.1);
            
            #Interruption Signal
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):            
            #Begin and End of the disturbance
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Amplitude Reduction
            mult_amp = uniform(0.1, 0.9);
            
            #Sag Signal
            sag_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = sag_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            #Begin and End of the disturbance
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Amplitude growth
            mult_amp = uniform(1.1, 1.8);
            
            #Swell Signal
            swell_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = swell_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def oscilation(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
        
            #Parameters
            #Oscilation Amplitude    
            amp_osc = uniform(0.3, 0.5)
            beta = 1
            gama = 1/uniform(0.003, 0.05)
            t1_size = randint(30,100)
            t1 = np.linspace(0.001, 0.002, t1_size)
            #Oscilation Frequency
            freq_o = randint(100, 4999)
            
            #Begin and End of the disturbance
            i_begin = randint(0, self.tam-200);
            
            #Oscilation
            h1 = beta*(np.exp(-gama*t1));
            h2 = np.sin(2*math.pi*t1*freq_o);
            h3 = h1*h2;
            
            osc = np.concatenate((np.zeros(i_begin), h3, np.zeros(self.tam-(i_begin+t1_size))), axis=0);
          
            #Oscilation Signal
            osc_seno = seno + amp_osc*osc; 
            
            matrix_r[i] = osc_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def oscilation_without_sine(self):
        matrix_r = np.empty([self.num, self.tam])
        for i in range(0, self.num):
                        
            #Parameters
            #Oscilation Amplitude    
            amp_osc = uniform(0.3, 0.5)
            beta = 1
            gama = 1/uniform(0.003, 0.05)
            t1_size = randint(30,100)
            t1 = np.linspace(0.001, 0.002, t1_size)
            #Oscilation Frequency
            freq_o = randint(100, 4999)
            
            #Begin and End of the disturbance
            i_begin = randint(0, self.tam-200);
            
            #Oscilation
            h1 = beta*(np.exp(-gama*t1));
            h2 = np.sin(2*math.pi*t1*freq_o);
            h3 = h1*h2;
            
            osc = np.concatenate((np.zeros(i_begin), h3, np.zeros(self.tam-(i_begin+t1_size))), axis=0);
          
            #Oscilation Signal
            osc_seno = amp_osc*osc; 
            
            matrix_r[i] = osc_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def flicker(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Flicker Parameters
            alpha = uniform(0.1, 0.2);
            freq_fck = uniform(0.1, 0.5);
            
            #Flicker
            flck = (1+(alpha*np.sin(2*np.pi*self.freq*self.t*freq_fck)));
            flck = flck*seno;
          
            #Flicker Signal        
            matrix_r[i] = flck;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def harmonic(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Number of harmonics
            num_harm = randint(1, 6);
            harm = seno;
            
            for j in range(0, num_harm):
                #Amplitude
                amp_h = uniform(0.05,0.2);
                #Frequency
                freq_h = randint(2, 40);
                #Harmonic plus signal
                harm = harm + amp_h*np.sin(2*np.pi*self.t*self.freq*freq_h);
            
            #Harmonic Signal        
            matrix_r[i] = harm;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Notching Amplitude
            alpha = uniform(0.01, 0.2);
            
            #Number of Notchings (To define the frequency)
            n_notch = randint(5,20);
            
            #Notching empty
            nch = np.zeros(self.tam);
            
            for j in range(0, n_notch):
                nch[(self.tam/n_notch)*j] = 1;
                    
            #Notching Signal        
            matrix_r[i] = alpha*nch + seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def notching_without_sine(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Amplitude Notching
            alpha = uniform(0.01, 0.2);
            
            #Number of Notchings (to define the frequency)
            n_notch = randint(5,20);
            
            #Notching empty
            nch = np.zeros(self.tam);
            
            for j in range(0, n_notch):
                nch[(self.tam/n_notch)*j] = 1;
                    
            #Notching without sine       
            matrix_r[i] = alpha*nch;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Amplitude Spike
            alpha = uniform(0.1, 0.5);
            
            #Begin and End of the disturbance
            i_low = randint(0,480);
            i_high = min(499,i_low+randint(2,20))
            
            #Spike empty
            spk = np.zeros(self.tam);
            spk[i_low:i_high] = 1;
                    
            #Spike Signal          
            matrix_r[i] = alpha*spk + seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def spike_without_sine(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Amplitude Spike
            alpha = uniform(0.1, 0.5);
            
            ##Begin and End of the disturbance
            i_low = randint(0,480);
            i_high = min(499,i_low+randint(2,20))
            
            #Spike empty
            spk = np.zeros(self.tam);
            spk[i_low:i_high] = 1;
                    
            #Spike without sine         
            matrix_r[i] = alpha*spk
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def dc_level(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #DC Variation
            nivel = uniform(-0.3, 0.3);
            nivel = nivel if nivel != 0 else 0.25;
            #Begin and End of the disturbance
            i_begin = randint(1,50);
            i_end = i_begin + randint(150,self.tam-50);
            
            #DC Level
            dc = np.concatenate((seno[0:i_begin], 
                                 seno[i_begin:i_end]+nivel, 
                                 seno[i_end:self.tam]), axis=0);
                    
            #DC-Level Signal          
            matrix_r[i] = dc;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            #Harmonic
            harmonic = self.harmonic();
                    
            #Sag + Harmonic Signal     
            matrix_r[i] = (sag[0] + harmonic[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell
            self.quant = 1;
            swell = self.swell();
            #Harmonic
            harmonic = self.harmonic();
                    
            #Swell + Harmonic Signal    
            matrix_r[i] = (swell[0] + harmonic[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            self.num = 1;
            seno = self.harmonic();
            seno = seno[0];
            
            #Begin and End of the disturbance
            i_begin = randint(0, 300);
            i_end = min(499,i_begin + randint(84,500));
            
            #Amplitude Reduction
            mult_amp = uniform(0, 0.1);
            
            #Interruption + Harmonic
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            #Interruption + Harmonic Signal
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Interruption
            self.quant = 1
            inte = self.interrupt
            
            #Notching
            notch = self.notching_without_sine
            
            #Interruption + Notching Signal
            matrix_r[i] = (inte[0] + notch[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            
            #Notching
            notching = self.notching_without_sine();
                    
            #Sag + Notching Signal            
            matrix_r[i] = (sag[0] + notching[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell
            self.quant = 1;
            swell = self.swell();
            
            #Notching
            notching = self.notching_without_sine();
                    
            #Swell + Notching Signal       
            matrix_r[i] = (swell[0] + notching[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_osc(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Interruption
            self.quant = 1;
            inte = self.interrupt
            
            #Oscilation
            osc = self.oscilation_without_sine
            
            #Interruption + Oscilation Signal            
            matrix_r[i] = (inte[0] + osc[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_oscilation(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            
            #Oscilation
            oscilation = self.oscilation_without_sine();
                    
            #Sag + Oscilation Signal        
            matrix_r[i] = (sag[0] + oscilation[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_oscilation(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell
            self.quant = 1;
            swell = self.swell();
            
            #Oscilation
            oscilation = self.oscilation_without_sine();
                    
            #Swell + Oscilation Signal    
            matrix_r[i] = (swell[0] + oscilation[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            
            #Spike Whithou Sine
            spike = self.spike_without_sine();
                    
            #Sag + Spike Signal          
            matrix_r[i] = (sag[0] + spike[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell
            self.quant = 1;
            swell = self.swell();
            
            #Spike Without Sine
            spike = self.spike_without_sine();
                    
            #Swell + Spike Signal       
            matrix_r[i] = (swell[0] + spike[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harmonics_oscilation(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag + Harmonics + Oscilation
            self.quant = 1;
            sag = self.sag();
            harmonic = self.harmonic();
            oscilation = self.oscilation_without_sine();
            
                    
            #Join Disturbances        
            matrix_r[i] = (sag[0] + harmonic[0] + oscilation[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harmonics_oscilation(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell + Harmonics + Oscilation
            self.quant = 1;
            swell = self.swell();
            harmonic = self.harmonic();
            oscilation = self.oscilation_without_sine();
                    
            #Join Disturbances  
            matrix_r[i] = (swell[0] + harmonic[0] + oscilation[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harmonics_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag + Harmonics + Notching
            self.quant = 1;
            sag = self.sag();
            harmonic = self.harmonic();
            notching = self.notching_without_sine();
            
                    
            ##Join Disturbances       
            matrix_r[i] = (sag[0] + harmonic[0] + notching[0]*2.0)/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harmonics_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell + Hamonics + Notching
            self.quant = 1;
            swell = self.swell();
            harmonic = self.harmonic()
            notching = self.notching_without_sine();
                    
            ##Join Disturbances    
            matrix_r[i] = (swell[0] + harmonic[0] + notching[0]*2.0)/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def noise(self, Signal, Desired_SNR_dB):
        Npts = Signal.shape[0] # Number of input time samples
        Noise = np.random.randn(Npts); # Generate initial noise; mean zero, variance one
        
        Signal_Power = np.sum(np.abs(Signal)*np.abs(Signal))/Npts;
        Noise_Power = np.sum(np.abs(Noise)*np.abs(Noise))/Npts;
        
        K = (Signal_Power/Noise_Power)*math.pow(10,(-Desired_SNR_dB/10));  #Scale factor
        
        New_Noise = np.sqrt(K)*Noise; #Change Noise level
        
        Noisy_Signal = Signal + New_Noise;
        return Noisy_Signal