# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:21:51 2021

@author: Jakob Häggström

The following scripts is a collection of functions that are used to encode
the different features. Note that seqfold is needed in order to use this 
script. Run 'pip install seqfold' before using this. 


"""

from seqfold import fold

import pandas as pd

import numpy as np

from itertools import product, combinations


def prod(seq, repeat):  # Returns all the possible di-nucleotides.

    return ["".join([c for c in tup]) for tup in list(product(seq, repeat=repeat))]


def pair_comb(seq, k, duplicates=False):

    if type(seq) == str:

        seq = list(seq)

    if duplicates:

        tempcomb = [sorted([tup[0], tup[1]])
                    for tup in list(combinations(seq, k))]

        for combi in [[c, c] for c in seq]:

            tempcomb.append(combi)

        return tempcomb

    else:

        return [sorted([tup[0], tup[1]]) for tup in list(combinations(seq, k))]


def count_contact_scontact(seq, featvec, get_feat=False):
    """
    INPUT: 

        seq: string of oligonucleotides.

        featvec: list of the chosen features.

    RETURNS:

        A list of the encoded features.
    """

    N = len(seq)

    seqlist = list(seq)

    dictvec = []  # Dict that stores all dicts of encoded features.

    if 'count' in featvec:  # This section counts the freq of each nucleotides.

        countdict = {'Fa': None, 'Ft': None, 'Fc': None, 'Fg': None}

        Na = seq.count('A')

        Nt = seq.count('T')

        Nc = seq.count('C')

        Ng = seq.count('G')

        count = np.array([Na, Nt, Nc, Ng])

        for i, key in enumerate(countdict.keys()):

            countdict[key] = count[i]

        dictvec.append(countdict)

    if 'contact' in featvec:  # Encodes the contact features.

        combvec = prod('ATGC', repeat=2)

        contactdict = {}

        # Iterating through all combinations. (AA, TC etc.)
        for comb in combvec:

            contactdict[comb] = 0

            # Iterating through the nucleotides two by two.
            for i in range(N-1):

                tempstr = "".join([seqlist[i], seqlist[i+1]])

                if comb == tempstr:  # If the dinucleotide is found, increment
                    # the frequency.
                    contactdict[comb] += 1

        dictvec.append(contactdict)

    if 'scontact' in featvec:

        scontactdict = {}

        if 'contact' in featvec:  # If you have chosen contact, you don't
            # need to create a new dict of the contact frequencys.
            un_combvec = pair_comb('ATGC', 2)

            # Done by simply adding freq(nuc) + freq(nuc_reverse)
            for comb in un_combvec:

                tempstr = f"{comb[0]}{comb[-1]}_{comb[-1]}{comb[0]}"

                scontactdict[tempstr] = (
                    contactdict[f"{comb[0]}{comb[-1]}"] + contactdict[f"{comb[-1]}{comb[0]}"])

            dictvec.append(scontactdict)

        else:  # Does the same but generates a dict of contact, since you need it
            # in order to calculate s_contact.
            combvec = prod('ATGC', repeat=2)

            un_combvec = pair_comb('ATGC', 2)

            contactdict = {}

            scontactdict = {}

            for comb in combvec:

                contactdict[comb] = 0

                for i in range(N-1):

                    tempstr = "".join([seqlist[i], seqlist[i+1]])

                    if comb == tempstr:

                        contactdict[comb] += 1

            for comb in un_combvec:

                tempstr = f"{comb[0]}{comb[-1]}_{comb[-1]}{comb[0]}"

                scontactdict[tempstr] = (
                    contactdict[f"{comb[0]}{comb[-1]}"] + contactdict[f"{comb[-1]}{comb[0]}"])

            for dup_comp in ['AA', 'TT', 'CC', 'GG']:  # Adds the duplicates.

                scontactdict[dup_comp] = contactdict[dup_comp]

            dictvec.append(scontactdict)

    input_list = []
    feat_list = []

    for in_dict in dictvec:  # Sums all features in one list.

        for key in in_dict.keys():

            input_list.append(in_dict[key])
            
            if get_feat:
                
                feat_list.append(key)
                
    if get_feat:
           
        return input_list, feat_list
    
    else:
        
        return input_list


def gen_hairpin(seq, temp=50):
    """

    INPUT: 

        seq: string of oligonucleotides.

    OPTIONAL:

        temp: Folding temperature in Celsius. default is 50. 

    RETURNS:

        A string that represents the secondary structure of the nucleotides.
        This script uses the notation '.' for nucleotides not binding to
        another nucleotides except for the adjacent ones, '(' and ')', for the
        nucletides within a stem. Example '...(....)..'.

    """

    # Uses seqfold.fold in order to calcutate the
    structs = fold(seq, temp=temp)
    # secoundary structure.
    m = ['.' for i in range(len(seq))]

    for struct in structs:

        if 'STACK' in struct.desc:  # Sets a stem if stem is found.

            m[struct.ij[0][0]] = '('
            m[struct.ij[0][1]] = ')'

    return "".join(m)


def count_secstruct(nucmap, get_feat = False):
    """

    INPUT: 

        seq: string of secoundary structure with the '.', '(' and ')' notation.


    RETURNS:

        A list with the encoded hairpin features. The order is:
            [n_free, n_stems, n_loop]

    """

    # The number of bases in stem is simply counted
    n_stems = 2 * nucmap.count('(')
    # by findind number of '(' and multiplying it with to, since '(' always has
    # a counterpart.
    if n_stems == 0:

        n_free = len(nucmap)

        n_loop = 0

    else:

        nuclist = list(nucmap)

        n_free = len(nuclist[:nuclist.index('(')])  # Counting number of free
        # nucleotides left of the first stem.
        # Same for the last
        ind = [i for i, n in enumerate(nuclist) if n == ')'][-1]

        # If there are nucleotides to the right of last stem.
        if ind + 1 < len(nuclist):

            n_free += len(nuclist[ind + 1:])

        # The number of nucleotides within loops.
        n_loop = len(nuclist) - n_free - n_stems
    
    if get_feat:
        
        return [n_free, n_stems, n_loop], ['n_free', 'n_stems', 'n_loop']
    
    else:
        
        return [n_free, n_stems, n_loop]


# This function combines gen_hairpin and count_secstruct.
def count_hairpin(seq, temp=50, get_feat = False):

    nucmap = gen_hairpin(seq, temp=temp)

    return count_secstruct(nucmap, get_feat)


def encode_input(seq, featvec, temp=50):

    count = []

    hair = []
    
    features = []

    check = [(feat in ['count', 'contact', 'scontact']) for feat in featvec]

    if True in check:

        if type(seq) == str:
            
            countvec, feat = count_contact_scontact(seq, featvec, get_feat = True)
            
            count.append(countvec)
            
            features = np.concatenate((features, feat))

        else:

            for i, nuc in enumerate(seq):
                
                if i == 0:
                    
                    countvec, feat = count_contact_scontact(nuc, featvec, get_feat = True)
                    
                    count.append(countvec)
                    
                    features = np.concatenate((features, feat))
                
                else:
                    
                    count.append(count_contact_scontact(nuc, featvec))

    if 'hairpin' in featvec:

        if type(seq) == str:
            
            countvec, feat = count_hairpin(seq, temp=temp, get_feat = True)
            
            hair.append(countvec)
            
            features = np.concatenate((features, feat))

        else:

            for i, nuc in enumerate(seq):
                
                if i == 0:
                    
                    countvec, feat = count_hairpin(nuc, temp=temp, get_feat = True)
                    
                    hair.append(countvec)
            
                    features = np.concatenate((features, feat))
                
                else:
                    
                    hair.append(count_hairpin(nuc, temp=temp))

        if True in check:
            
            data = np.concatenate((count, hair), axis=1)

            return pd.DataFrame(data, columns = features)

        else:

            return pd.DataFrame(np.array(hair), columns = features)

    return pd.DataFrame(np.array(count), columns = features)


if __name__ == "__main__":
    featvec = ['count', 'contact','scontact', 'hairpin']
    l = ['ATGC', 'AAATTTATCTGA']
    f = encode_input(l, featvec)
