"""
This module contains implementations of the ML representations for compositions
that were used throughout the hyperparameter tuning and model testing.
"""

import json
import torch
import os
import ast
import numpy as np
from copy import deepcopy
from functools import partial
from pandarallel import pandarallel

pandarallel.initialize(progress_bar=False, nb_workers=32)

# e.g. '~/venv/lib/python3.10/site-packages/ElMD/el_lookup/mod_petti.json'
ELMD_PATH = '<INSERT ElMD INSTALL PATH>/el_lookup/mod_petti.json'

def tryexcept(func, *errors):
    def featurize(value):
        try:
            return func(value)
        except errors:
            return np.NaN
    return featurize

def sort_formula_by_enegs(c: str, enegs: dict):
    # Tokenize composition
    comp = c.replace('\n', ' ')
    comp_list = comp.split(' ')
    # Get elements
    el_list = [el.rstrip('0123456789+-.') for el in comp_list]
    # Get electronegativities
    eneg_list = [enegs[el] for el in el_list]
    # Exception catcher
    if not all(eneg_list):
        return comp
    # Sort composition
    comp_list_sorted = [
        el for _, el in sorted(zip(eneg_list, comp_list))
    ]
    return ' '.join(comp_list_sorted)

def get_embedding_vectors(c, embedding=""):
    from elementembeddings.composition import CompositionalEmbedding
    compemb = CompositionalEmbedding(c, embedding=embedding)
    return np.multiply(compemb.norm_stoich_vector[:,np.newaxis],
                       compemb.el_matrix)

def pad_embedding_vectors(vs, npad=16):
    rep = np.zeros((npad, vs.shape[1]))
    if vs.shape[0] > npad:
        rep[:,:] = vs[:npad,:]
    else:
        rep[:vs.shape[0],:] = vs
    return rep.tolist()


class RepresentationGenerator:
    EMBEDDING_REPS = [
        'cgnf',
        'crystallm',
        'magpie',
        'mat2vec',
        'matscholar',
        'megnet',
        'oliynik',
        'random',
        'skipatom',
        'xenopy'
    ]
    def __init__(self, df, dim=1) -> None:
        self.df = deepcopy(df)
        self.dim = dim
        

    def get_representations(self, **kwargs):
        if self.dim == 1:
            return self.get_1d_rep(**kwargs)
        elif self.dim == 2:
            return self.get_2d_rep(**kwargs)
        elif self.dim == -1:
            return self.get_graph_rep()
        else:
            print(f'Dimensionality of data not implemented!')
            return None, None

    # Ouptut of shape (Ndata, Nrep)
    def get_1d_rep(self, rep_type='', **kwargs):
        if rep_type == 'comp-embedding':
            return self._from_element_embeddings(**kwargs)
        elif rep_type == 'structure':
            return self._from_matminer()
        elif rep_type == 'joint-embedding':
            return self._combined_element_embeddings_matminer(**kwargs)
        else:
            print(f'Representation type {rep_type} unknown!')
            return None, None
        
    def _from_element_embeddings(self, embedding='', **kwargs):
        from elementembeddings.composition import composition_featuriser
        from elementembeddings.core import Embedding
        # Load embedding data
        emb = Embedding.load_data(embedding)
        # Add symbols M, X, L
        null = np.zeros_like(emb.embeddings['H'])
        emb.embeddings['M'] = deepcopy(null)
        emb.embeddings['X'] = deepcopy(null)
        emb.embeddings['L'] = deepcopy(null)
        emb.element_list.extend(['M', 'X', 'L'])
        # Featurize the dataframe 
        df_featurized = composition_featuriser(
            self.df['composition'].tolist(),
            embedding=emb,
            stats=kwargs.pop('stats', ['mean'])
        )
        # Convert dataframe to tensor
        X = torch.tensor(np.array(df_featurized))
        y = torch.tensor(self.df['disordered'].values, dtype=int) # .float()
        return X, y

    def _from_matminer(self):
        from matminer.featurizers.site import CrystalNNFingerprint
        from matminer.featurizers.structure import SiteStatsFingerprint
        from pymatgen.core import Structure
        # Convert structure strings to Structure objects
        self.df['structure'] = self.df['structure'].parallel_apply(json.loads)
        self.df['structure'] = self.df['structure'].parallel_apply(Structure.from_dict)
        # Prepare featurization
        ssf = SiteStatsFingerprint(
            CrystalNNFingerprint.from_preset(
                'ops',
                distance_cutoffs=None,
                x_diff_weight=0,
                porous_adjustment=False
            ),
            stats=('mean', 'std_dev', 'minimum', 'maximum')
        )
        featurize = tryexcept(ssf.featurize, ValueError)
        # Featurization
        self.df['structure'] = self.df['structure'].parallel_apply(featurize)
        # Take care of NaN values
        Nnan = self.df['structure'].isna().sum()
        if Nnan > 0:
            print(f'Masking {Nnan} entries with null vector!')
            idx = self.df['structure'].first_valid_index()
            nullvec = np.zeros_like(self.df['structure'].loc[idx]).tolist()
            self.df['structure'] = self.df['structure'].parallel_apply(
                lambda d: d if isinstance(d, list) else nullvec
            )
        # Convert to tensor
        X = torch.tensor(self.df['structure'].values.tolist())
        y = torch.tensor(self.df['disordered'].values, dtype=int)
        return X, y
    
    def _combined_element_embeddings_matminer(self, **kwargs):
        X_comp, y = self._from_element_embeddings(**kwargs)
        X_struc, _ = self._from_matminer()
        return torch.concat([X_comp, X_struc], dim=-1), y


    # Ouptut of shape (Ndata, Nrep1, Nrep2)
    def get_2d_rep(self, rep_type='', **kwargs):
        if rep_type == 'chemsys':
            return self._chem_sys_matrix()
        elif rep_type == 'chemsys-petti':
            return self._chem_sys_matrix_petti()
        elif rep_type == 'comp-rnn':
            return self._comp_weighted_embedding_sequence(**kwargs)
        else:
            print(f'Representation type {rep_type} unknown!')
            return None, None

    def _chem_sys_matrix_petti(self):
        # Get periodic table
        mod_petti_file = ELMD_PATH
        with open(mod_petti_file, 'r') as jfile:
            mod_petti = json.load(jfile)
        masked_symbols = ['X', 'M', 'L']
        for k in masked_symbols:
            mod_petti[k] = 0
        petti_nums = list(mod_petti.values())
        nperiodic = max(petti_nums) + 1
        # Initialize representation tensor(s)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        y = torch.zeros((size[0], 1)) # .long()
        # Fill tensors from dataframe
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs_petti'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            for elpair in elpairs:
                indices.append([i, elpair[0], elpair[1]])
            y[i] = int(row['disordered'])
        vsize = len(indices)
        indices = torch.Tensor(indices)
        values = torch.ones(vsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        return X.to_dense(), y

    def _chem_sys_matrix(self):
        from ase.data import chemical_symbols
        nperiodic = len(chemical_symbols)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        y = torch.zeros((size[0], 1))
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            for elpair in elpairs:
                indices.append([i, elpair[0], elpair[1]])
            y[i] = int(row['disordered'])
        vsize = len(indices)
        indices = torch.Tensor(indices)
        values = torch.ones(vsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        return X.to_dense(), y

    def _comp_weighted_embedding_sequence(self, npad=16, embedding=""):
        ### Load electronegativities for sorting
        from elementembeddings import __file__ as ELEMB_DIR
        module_directory = os.path.abspath(os.path.dirname(ELEMB_DIR))
        data_directory = os.path.join(module_directory, "data")
        pt_dir = os.path.join(data_directory, "element_data", "periodic-table-lookup-symbols.json")
        with open(pt_dir) as f:
            pt = json.load(f)
        enegs = {el: v['electronegativity_pauling'] for el, v in pt.items()}
        sort_by_eneg = partial(sort_formula_by_enegs, enegs=enegs)
        ### Sort compositions by electronegativities
        self.df['composition'] = self.df['composition'].parallel_apply(
            sort_by_eneg
        )
        ### Tokenize compositions
        get_embeddings = partial(get_embedding_vectors, embedding=embedding)
        self.df['composition'] = self.df['composition'].parallel_apply(
            get_embeddings
        )
        ### Padding
        pad_embeddings = partial(pad_embedding_vectors, npad=npad)
        self.df['composition'] = self.df['composition'].parallel_apply(
            pad_embeddings
        )
        ### Convert dataframes to tensors
        X = torch.tensor(self.df['composition'].values.tolist())
        print(f'Input shape: {X.shape}')
        y = torch.tensor(self.df['disordered'].values, dtype=int)
        print(f'Target shape: {y.shape}')
        return X, y


class BaselineRepresentationGenerator:
    def __init__(self, df):
        self.df = deepcopy(df)

    def get_representations(self, split='', **kwargs):
        if split == 'train':
            return self.get_train_rep(**kwargs)
        elif split == 'test':
            return self.get_test_rep(**kwargs)
        else:
            print(f'Representation not implemented!')
            return None, None
    
    ### Representations for baseline model evaluation
    # Ouptut of shape (Ndata, Nrep1, Nrep2), (Ndata,)
    def get_test_rep(self, rep_type=''):
        if rep_type == 'chemsys':
            return self._chem_sys_matrix()
        elif rep_type == 'chemsys-petti':
            return self._chem_sys_matrix_petti()
        else:
            print(f'Representation type {rep_type} unknown!')
            return None, None
        
    def _chem_sys_matrix_petti(self):
        # Get periodic table
        mod_petti_file = ELMD_PATH
        with open(mod_petti_file, 'r') as jfile:
            mod_petti = json.load(jfile)
        masked_symbols = ['X', 'M', 'L']
        for k in masked_symbols:
            mod_petti[k] = 0
        petti_nums = list(mod_petti.values())
        nperiodic = max(petti_nums) + 1
        # Initialize representation tensor(s)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        y = torch.zeros((size[0], 1))
        # Fill tensors from dataframe
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs_petti'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            for elpair in elpairs:
                indices.append([i, elpair[0], elpair[1]])
            y[i] = int(row['disordered'])
        vsize = len(indices)
        indices = torch.Tensor(indices)
        values = torch.ones(vsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        return X.to_dense(), y
    
    def _chem_sys_matrix(self):
        from ase.data import chemical_symbols
        nperiodic = len(chemical_symbols)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        y = torch.zeros((size[0], 1)) # .long()
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            for elpair in elpairs:
                indices.append([i, elpair[0], elpair[1]])
            y[i] = int(row['disordered'])
        vsize = len(indices)
        indices = torch.Tensor(indices)
        values = torch.ones(vsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        return X.to_dense(), y    


    ### Representations to generate element pair probability matrix
    # Ouptut of shape (Ndata, Nrep1, Nrep2), (Ndata, Nrep1, Nrep2)
    def get_train_rep(self, rep_type=''):
        if rep_type == 'chemsys':
            return self._site_label_chem_sys_matrix()
        elif rep_type == 'chemsys-petti':
            return self._site_label_chem_sys_matrix_petti()
        else:
            print(f'Representation type {rep_type} unknown!')
            return None, None
    
    def _site_label_chem_sys_matrix(self):
        from ase.data import chemical_symbols
        nperiodic = len(chemical_symbols)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        disorder_indices = []
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            elpair_disorder = ast.literal_eval(row['element_pairs_disorder'])
            for elpair, disorder in zip(elpairs, elpair_disorder):
                indices.append([i, elpair[0], elpair[1]])
                if disorder:
                    disorder_indices.append([i, elpair[0], elpair[1]])
        vsize = len(indices)
        dsize = len(disorder_indices)
        indices = torch.Tensor(indices)
        disorder_indices = torch.Tensor(disorder_indices)
        values = torch.ones(vsize)
        disorder_values = torch.ones(dsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        y = torch.sparse_coo_tensor(disorder_indices.t(), disorder_values, size, dtype=torch.int)
        return X, y
    
    def _site_label_chem_sys_matrix_petti(self):
        # Get periodic table
        mod_petti_file = ELMD_PATH
        with open(mod_petti_file, 'r') as jfile:
            mod_petti = json.load(jfile)
        masked_symbols = ['X', 'M', 'L']
        for k in masked_symbols:
            mod_petti[k] = 0
        petti_nums = list(mod_petti.values())
        nperiodic = max(petti_nums) + 1
        # Initialize representation tensor(s)
        size = (self.df.shape[0], nperiodic, nperiodic)
        indices = []
        disorder_indices = []
        # Fill tensors from dataframe
        for i, (idx, row) in enumerate(self.df.iterrows()):
            elpairs_str = row['element_pairs_petti'].replace('(', '[').replace(')', ']')
            elpairs = json.loads(elpairs_str)
            elpair_disorder = ast.literal_eval(row['element_pairs_disorder'])
            for elpair, disorder in zip(elpairs, elpair_disorder):
                indices.append([i, elpair[0], elpair[1]])
                if disorder:
                    disorder_indices.append([i, elpair[0], elpair[1]])
        vsize = len(indices)
        dsize = len(disorder_indices)
        indices = torch.Tensor(indices)
        disorder_indices = torch.Tensor(disorder_indices)
        values = torch.ones(vsize)
        disorder_values = torch.ones(dsize)
        X = torch.sparse_coo_tensor(indices.t(), values, size, dtype=torch.double)
        y = torch.sparse_coo_tensor(disorder_indices.t(), disorder_values, size, dtype=torch.int)
        return X, y


class TorchStandardScaler:
  def __init__(self) -> None:
      self.mean = None
      self.std = None

  def fit(self, x):
    self.mean = x.mean(0, keepdim=True)
    self.std = x.std(0, unbiased=False, keepdim=True)

  def transform(self, x):
    x -= self.mean
    x /= (self.std + 1e-7)
    return x
