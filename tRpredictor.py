# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 11:21:34 2021

@author: Jakob Häggström

The following class stores multiple sklearn models that predicts the retention
factor for a given oligonucleotide. The features of the input is divided into
several groups as per Strum article.

Features:
    
    'count': Frequency of each base
    
    'contact': Frequency of all possible di-nucleotides (AA, AC, CA, AT etc)
    
    'scontact': Frequency of all possible di-nucleotides neglecting
    the position (AT + TA, AC + CA, AA, etc.)
    
    'hairpin': Secondary structure features, number of nucleotides in stem, in
    loop and free.
    
The user can choose between these at the moment. More features might be added.

The models has to be trained before used, and saved in a .pkl file. 

"""

from pickle import load
import sklearn.svm as svm
import numpy as np
import tRpredictor_func as tR_f




class tR_predictor:
    
    def __init__(self, pathvec, featvec, t0 = None, td = None):
        
        """
        ATTRIBUTES:
            
            pathvec: A vector of paths to the saved Sklearn models in .pkl 
            format. type: list or ndarray of str.
            
            featvec: A vector of the chosen features described as above.
            type: list or ndarray of str.
            
        OPTIONAL:
            
            t0: The measured t0 of the system.
            
            td: The dwell time.
        """
        
        self.featvec = featvec
        
        self.modelvec = []
        
        for path in pathvec:
            
            self.modelvec.append(load(open(path, 'rb'))) # Loads all the models. 
        
        self.t0 = t0
        
        self.td = td
    
    def __str__(self):
        pass
    
    def predict_tR(self, olig, ret_factor = True):
        
        """
        INPUT:
            
            olig: str or vector of strings with oligonucleotides.
            
        OPTIONAL:
            
            ret_factor: Choice of unit in the predicted retention time.
            Retention factor if True (default), minutes if False. Though
            td and t0 is required if minutes is desired. (bool) 
        
        RETURNS:
            
            The retention time in retention factor or minutes. (float)
            
        """
        
        olig_input = tR_f.encode_input(olig, self.featvec)
        olig_input = olig_input.values
        
        if type(olig) == str: 
            
            tR = np.zeros((1, len(self.modelvec)))
            
        else:
            
            tR = np.zeros((len(olig), len(self.modelvec)))
        
        for i, model in enumerate(self.modelvec): # Makes prediction for each gradient.
            
            tR[:, i] = model.predict(olig_input)
        
        if ret_factor:
            
            return tR
        
        else:
            
            if self.t0 and self.td != None:
            
                return tR * self.t0 + self.t0 + self.td
            
            else:
                
                raise ValueError('No t0 or td is given.')
        
    
    
    