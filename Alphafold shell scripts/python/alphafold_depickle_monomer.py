# modification of https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py
# saves PAEs and pLDDTs of AF2 multimer run as a json file and makes pretty plots as in ColabFold
# operates on all result_.pkl files in an input folder as long as
#
# INPUT: specify the input in the command line according to python alphafold_depickle -h with --var='value' arguments
#
# OUTPUT in output_folder (default is the same as the specified input folder):
# - a PNG called Coverage_LDDT.png showing MSA coverage and per-residue LDDT predictions
# - a PNG called PAE.png showing pairwise PAE matrices, plotted on a grid with the size models x predictions
# - .json files for each prediction that contain PAE as a list of rows, pLDDT as a list, the
#   pTM and maximal PAE as integers, as well as a list of protein names per chain and their lengths (starting from chain A)
#   to plot the chain boundaries in the complex. The names of the .jsons are the same as the model names to allow easy
#   visualization in ChimeraX via "Tools -> Structure Prediction -> AlphaFold Error Plot" for each model.
#   .json outputs are formatted as in ColabFold, https://github.com/sokrypton/ColabFold, 
#   and can be read in with a custom Jupyter Notebook script to make pretty plots of the PAE data

import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle
import json

def get_pae_plddt(model_names): #grabs values from the pickle files and returns them as a dict
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
        out[name[:-4]] = {'plddt': d['plddt'], 'pae':d['predicted_aligned_error'], 
                               'max_pae':d['max_predicted_aligned_error'],'ptm':d['ptm']}
    return out

def generate_output_images(feature_dict, out_dir, name, count, num_model, pae_plddt_per_model): # plots data
    msa = feature_dict['msa'] #acquire axis limits for msa plot
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    plt.figure(figsize=(14, 4), dpi=200)
	
    ##################################################################
    plt.subplot(1, 2, 1)
    plt.title("Sequence coverage")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    ##################################################################
    plt.subplot(1, 2, 2)
    plt.title("Predicted LDDT per position")
    
    model_indices = np.repeat(list(range(1,num_model+1)),count)
    run_indices = np.tile(list(range(0,count)),num_model)
    
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.plot(value["plddt"],label='model_'+str(model_indices[n])+'_pred_'+str(run_indices[n]))
        
    #plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}{name+('_' if name else '')}coverage_LDDT.png")
    ##################################################################

    ##################################################################
    plt.figure(figsize=(3*num_model*count, 2), dpi=200)
  
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(1, count*num_model, n + 1)
        plt.title('model_'+str(model_indices[n])+
                               '_pred_'+str(run_indices[n]), fontsize = 7)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=35)
        plt.colorbar()
        
        #save information as JSON in ColabFold format
        with open(model_name+'.json','w') as f: 
            print('Saving '+model_name+'.json to output folder')
            
            #round the values from 10 decimals to 2 to reduce filesize
            pae = value["pae"].tolist()
            pae_rounded = np.round(pae,2).tolist()
            
            plddt = value["plddt"].tolist()
            plddt_rounded = np.round(plddt,2).tolist()
            
            #pack into json
            json_data = {'pae': pae_rounded,'plddt':plddt_rounded, 
                         'max_pae':value["max_pae"].tolist(),'ptm':value["ptm"].tolist()}
            json.dump(json_data, f, ensure_ascii=False)

    plt.savefig(f"{out_dir}{name+('_' if name else '')}PAE.png")
    ##################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',dest='input_dir',required=True, 
                    help='(Required input) Input directory with .pkl files, given as a str containing an absolute path.')
parser.add_argument('--name',dest='name', 
                    help='Custom filename extension as a str. Default will be the filename of the .pkl files for output .json.')
parser.set_defaults(name='')
parser.add_argument('--num_model',dest='num_model', type=int,
                    help = 'Number of models used by AF2. Default is 5. This parameter should only change if you want to analyze a subset of the .pkl files.')
parser.set_defaults(num_model=5)
parser.add_argument('--run_count',dest='run_count', type=int,
                    help = "Number of iterations per model, starting from 0 to run-count-1 as an index in the filename.")
parser.set_defaults(run_count=1)
parser.add_argument('--output_dir',dest='output_dir',
                    help = 'Outpath directory path for .json files and plots, given as a str containing an absolute path. Defaults to the input directory if not specified.')
parser.set_defaults(output_dir='')
args = parser.parse_args()

print('> Loading AF2 run information in '+args.input_dir)

feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb')) #get MSA information

model_names = [] #get names of files to extract PAE from
for model_set in range(1,args.num_model+1):
	string = str(args.input_dir)+"/result_model_"+str(model_set)+"_ptm_pred_0.pkl"
	print('Loaded .pkl file '+string)
	model_names.append(string)

pae_plddt_per_model = get_pae_plddt(model_names)
print('> Run info extracted, making output images and saving .json files')
generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, 
                       args.name, args.run_count, args.num_model, pae_plddt_per_model)
