# function definition script for jupyter notebook scripts
# provide absolute path in initialization cell

print('Defining functions...')

import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn', turns off some warnings from pandas

# function to return distance maps from chain pair from PDB files
def distmap_calc(chain_ids,pdbstructure):
    from evcouplings.compare.pdb import ClassicPDB
    from evcouplings.compare import (
        pdb, DistanceMap, intra_dists,
        multimer_dists, coupling_scores_compared
    )
    
    chain1 = pdbstructure.get_chain(str(chain_ids[0]),0) #grab the first chain in tuple
    #chain1 = chain1.filter_atoms(atom_name='CB') # uncomment to only get C_beta distances
    
    if len(chain_ids) == 2: # for pairwise chain distances
        chain2 = pdbstructure.get_chain(str(chain_ids[1]),0) #grab the second chain in tuple
        #chain2 = chain2.filter_atoms(atom_name='CB') #uncomment to only get C_beta distances

        distmap_pdb = DistanceMap.from_coords(chain1, chain2)
        
    else: #for intra-chain distances
        distmap_pdb = DistanceMap.from_coords(chain1)
        
    return distmap_pdb

# rescaling of couplings in contact maps based on specified probability in 'size' column of dataframes
def rescale_size(df,flag):
    if len(df) > 5:
        df.loc[:,flag] = (
                            0.25+0.75* #rescale to a coupling point scaling size based on probability
                            (df.loc[:,flag] - min(df.loc[:,flag]))
                            /
                            (max(df.loc[:,flag]) - min(df.loc[:,flag]))
        )
    return df

# functions to grab the boundaries for the plots
def boundary_returner_intra(boundaries, intra_ecs, distance_intra_i):
    from evcouplings.visualize import find_boundaries
    detected_boundaries = list(
        find_boundaries(
            boundaries, intra_ecs, distance_intra_i,
            symmetric=True, multimer=None
        )
    )     
    return detected_boundaries

def boundary_returner_inter(intraA_boundaries, intraB_boundaries):
    from evcouplings.visualize import find_boundaries
    detected_boundaries = [#[(A[0][0],A[0][1]),(B[0][0],B[0][1])]
    (intraA_boundaries[0][0],intraA_boundaries[0][1]),
    (intraB_boundaries[0][0],intraB_boundaries[0][1])
    ]
    return detected_boundaries

# threshold a couplings list by minimum probability for plotting
def Pthresh(couplings,Pcutoff):
    couplings_out = couplings.loc[couplings['probability']>Pcutoff]
    return couplings_out
    
# function for automated annotation of pie chart with percentages + number of couplings per slice
def make_autopct(values): #function def for annotation
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.1f}%\n({v:d})'.format(p=pct,v=val)
    return my_autopct    
    
# generates ChimeraX distance command between residue i and j in two chains
distfun = lambda i,j,ch1,ch2,col:\
"distance #1/"+ch1+":"+str(i)+"@CA #1/"\
+ch2+":"+str(j)+"@CA color "+col+" radius .33 symbol false dashes 0\n"

# generates ChimeraX command to show atoms in all chains
showfun = lambda i,j,chainlist: "show #1/"+','.join(chainlist)+":"+str(i)+","+str(j)+ " atoms\n"

# writer function to write a ChimeraX command that draws couplings in specified chains 
def writerfun(tbl,idx,chain1,chain2,col,arg,flg):
    writer = distfun(tbl['i'][idx],tbl['j'][idx],chain1,chain2,col)
    if flg:
        #write the header only if flg is passed as True
        header = ('# '+arg+' - '+str(idx+1)+', '
        'P: '+str(round(tbl['probability'][idx],4))+
        ', cn: '+str(round(tbl['cn'][idx],4))+
        ', dist: '+str(round(tbl['dist'][idx],4))+
        '\n')
        writer = header+writer
    return writer
print('Functions initalized!')