# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 22:10:37 2021

@author: Jakob Häggström

For more information see:
Journal of Chromatography A, xxx (2022) xxxx–xxxx
doi:xxx

The following program is more or less a user interface that uses the trained 
models to predict retentiontimes for a given oligonucleotide sequence and dataset
(PM (C18 method) or IPM (Ion pair method)) for three different gradients. 
The gradients in the dataset were:
PM: 
[G1, G2, G3] = [2.22, 1.23, 0.81] v%MeCN min^-1
IPM: 
[G1, G2, G3] = [0.32, 0.16, 0.08] mM TEA min^-1
The output from the model will follow the same order.
"""


import tRpredictor as tR_pre
import os



def get_model_path(MAIN_DIR,dataset):
    
    gradient = [f"G{i + 1}" for i in range(3)]
    
    pathvec = []
    
    for G in gradient:
        
        filename = f"SVR_{dataset}_{G}.sav"
        
        pathvec.append(os.path.join(MAIN_DIR, 'SVR_models', filename))
    
    return pathvec


def load_models(MAIN_DIR,dataset):
    
    pathvec = get_model_path(MAIN_DIR, dataset)
    
    featvec = ['count']
    
    models = tR_pre.tR_predictor(pathvec, featvec)

    return models


def predict(dataset,seq):
    """
    

    Parameters
    ----------
    dataset : String, either PM or IPM
    
    seq : String, or list of strings, string of combinations of ATGC.

    Returns
    -------
    1X3 list of floats, where each elements correnspond to each gradient dataset. 
    I.e [[time gradient 1, time gradient 2, time gradient 3]]

    """
    
    MAIN_DIR = os.path.split(os.path.abspath(__file__))[0] 
    
    model = load_models(MAIN_DIR, dataset)
    
    return model.predict_tR(seq)


def main():
    
    # Here you can input arbitrary sequence(s) 
    #The most important is that you only make it a single list
    # of sequences.
    
    
    #Example sequences Sample S12A and S16C (Supplementary Material Table S1)
    seq = ['CCCACACCCAAC','ATTTTTGTGCGCTCTA']
    
    #Expected output
    #[[ 7.21372445  8.74126582  9.93287348]
    # [ 8.44104151 11.15892085 13.67501449]]
    #[[ 6.42627625  7.09895575  7.6196043 ]
    # [ 8.96175395 11.58446108 14.66242965]]
    
    # Here you input the sequence or sequences into the model, and you get 
    # the retention time for three different gradients for the PM or the IPM
    # dataset. 
    out_PM = predict('PM',seq) # Output for PM 
    out_IPM = predict('IPM',seq) # Output for IPM
    
    # Print the output.
    print(f"{out_PM}")
    print(f"{out_IPM}")

if __name__ == '__main__':
    
    main()
    
    
    
    
    
    
    