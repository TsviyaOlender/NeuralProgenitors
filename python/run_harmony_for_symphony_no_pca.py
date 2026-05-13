#from tkinter import W
import numpy as np
import pandas as pd
import scanpy as sc
#import shoji
import symphonypy as sp
import os.path
import sys
import argparse

parser = argparse.ArgumentParser(description='Run harmony before symphony')
parser.add_argument('-ar', '--anndata_ref', type=str)
#parser.add_argument('-ws', '--workspace', type=str)
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args()
#workspace_name = args.workspace
output_file = args.output_file

if(args.anndata_ref is not None):
    adata_ref_file = args.anndata_ref
    adata_ref = sc.read(adata_ref_file)

else:
    sys.exit("Can't find workspace ot anndata")

sc.pp.scale(adata_ref, max_value=10)
sc.pp.pca(adata_ref, n_comps=30, zero_center=False,svd_solver='arpack')
sc.pp.neighbors(adata_ref, n_neighbors=15)
sc.tl.umap(adata_ref,min_dist=0.1)
adata_ref.write(adata_ref_file)
# You can skip Harmony if you have only one batch in reference
sp.pp.harmony_integrate(adata_ref, key='Chemistry', max_iter_harmony=30)
#sp.pp.harmony_integrate(adata_ref, key='SampleID', max_iter_harmony=30)
sc.pp.neighbors(adata_ref, use_rep="X_pca_harmony")
sc.tl.umap(adata_ref,min_dist=0.1)
#sc.tl.louvain(adata_ref, resolution=2)
sc.tl.leiden(adata_ref, resolution=0.6)
adata_ref.write(output_file+'.h5ad')
