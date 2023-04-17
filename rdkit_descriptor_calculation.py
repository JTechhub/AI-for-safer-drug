# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:37:59 2019

@author: Gil

sd: list of rdkit Chem mol instance.
"""

from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys                       #MACCS Keys (166 bit)
from rdkit.Chem.Fingerprints import FingerprintMols    #Topological Fingerprints (RDKFingerprint, 2048 bit)
from rdkit.Avalon import pyAvalonTools                 #512 bit
from rdkit.Chem.AtomPairs import Pairs, Torsions       #Atom pair & torsions fingerprints

from rdkit.Chem.EState.Fingerprinter import FingerprintMol # https://www.rdkit.org/docs/source/rdkit.Chem.EState.Fingerprinter.html (79 bit)
from rdkit.Chem.EState import EState_VSA # https://www.rdkit.org/docs/source/rdkit.Chem.EState.EState_VSA.html
from rdkit.Chem import GraphDescriptors #https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html
from rdkit.Chem import Descriptors #https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
from rdkit.Chem import Crippen #https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
from rdkit.Chem import MolSurf #https://www.rdkit.org/docs/source/rdkit.Chem.MolSurf.html
from rdkit.Chem import rdMolDescriptors #https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
from rdkit.Chem import Lipinski

from rdkit import Chem, DataStructs

import numpy as np
import pandas as pd

def fingerprint_bitgenerator(fps):
    np_fps = []
    
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp,arr)
        np_fps.append(arr)
        
    return np_fps

def atompairfingerprint_keyextractor(fps):
    fpkeys = []
    
    for fp in fps:
        fpkeys += list(fp.GetNonzeroElements().keys())
    
    fpkeys = list(set(fpkeys))
    fpkeys.sort()
    return fpkeys

def atompair_bitgenerator(fps, keys):
    
    keylen = len(keys)
    np_fps = []
    for fp in fps:
        arr = np.zeros(keylen)
        for i in range(keylen):
            arr[i]=fp[keys[i]]
        np_fps.append(arr)
        
    return np_fps

def description2D_calc(sd):
    #Topliss fragments / BCUT2D
    desc_2d = []
    for m in sd:
        each_desc = []
        
        desc = [rdMolDescriptors.CalcNumAliphaticCarbocycles(m),
                rdMolDescriptors.CalcNumAliphaticHeterocycles(m),
                rdMolDescriptors.CalcNumAliphaticRings(m),
                rdMolDescriptors.CalcNumAmideBonds(m),
                rdMolDescriptors.CalcNumAromaticCarbocycles(m),
                rdMolDescriptors.CalcNumAromaticHeterocycles(m),
                rdMolDescriptors.CalcNumAromaticRings(m),
                rdMolDescriptors.CalcNumAtomStereoCenters(m),
                rdMolDescriptors.CalcNumBridgeheadAtoms(m),
                rdMolDescriptors.CalcNumHBA(m),
                rdMolDescriptors.CalcNumHBD(m),
                rdMolDescriptors.CalcNumHeteroatoms(m),
                rdMolDescriptors.CalcNumHeterocycles(m),
                rdMolDescriptors.CalcNumLipinskiHBA(m),
                rdMolDescriptors.CalcNumLipinskiHBD(m),
                rdMolDescriptors.CalcNumRings(m),
                rdMolDescriptors.CalcNumRotatableBonds(m),
                rdMolDescriptors.CalcNumSaturatedCarbocycles(m),
                rdMolDescriptors.CalcNumSaturatedHeterocycles(m),
                rdMolDescriptors.CalcNumSaturatedRings(m),
                rdMolDescriptors.CalcNumSpiroAtoms(m),
                
                Lipinski.FractionCSP3(m),
                Lipinski.HeavyAtomCount(m),
                Lipinski.NHOHCount(m),
                Lipinski.NOCount(m),

                Descriptors.ExactMolWt(m),
                Descriptors.HeavyAtomMolWt(m),
                Descriptors.MaxAbsEStateIndex(m),
                Descriptors.MaxEStateIndex(m),
                Descriptors.MinAbsEStateIndex(m),
                Descriptors.MinEStateIndex(m),
                Descriptors.MaxAbsPartialCharge(m),
                Descriptors.MaxPartialCharge(m),
                Descriptors.MinAbsPartialCharge(m),
                Descriptors.MinPartialCharge(m),
                Descriptors.NumRadicalElectrons(m),
                Descriptors.NumValenceElectrons(m),
                
                Crippen.MolLogP(m),
                Crippen.MolMR(m),
                    
                MolSurf.LabuteASA(m),
                MolSurf.TPSA(m)]
        
        graph_desc = [GraphDescriptors.BalabanJ(m),
                      GraphDescriptors.BertzCT(m),
                      GraphDescriptors.Ipc(m),
                      GraphDescriptors.Chi0(m),
                      GraphDescriptors.Chi0n(m),
                      GraphDescriptors.Chi0v(m),
                      GraphDescriptors.Chi1(m),
                      GraphDescriptors.Chi1n(m),
                      GraphDescriptors.Chi1v(m),
                      GraphDescriptors.Chi2n(m),
                      GraphDescriptors.Chi2v(m),
                      GraphDescriptors.Chi3n(m),
                      GraphDescriptors.Chi3v(m),
                      GraphDescriptors.Chi4n(m),
                      GraphDescriptors.Chi4v(m),
                      GraphDescriptors.HallKierAlpha(m),                        
                      GraphDescriptors.Kappa1(m),
                      GraphDescriptors.Kappa2(m),
                      GraphDescriptors.Kappa3(m)]
        
        estate_vsa = list(EState_VSA.EState_VSA_(m))
        vsa_estate = list(EState_VSA.VSA_EState_(m))
        
        peoe_vsa = list(MolSurf.pyPEOE_VSA_(m))
        smr_vsa = list(MolSurf.pySMR_VSA_(m))
        logp_vsa = list(MolSurf.pySlogP_VSA_(m))

        mqns = rdMolDescriptors.MQNs_(m)

        each_desc = desc + graph_desc + estate_vsa + vsa_estate + peoe_vsa + smr_vsa + logp_vsa + mqns
        
        desc_2d.append(each_desc)
        
    desc_indices = ['num_aliphatic_carbocycles','num_alphatic_heterocycles',
                    'num_alphatic_ring','num_amidebond',
                    'num_aromatic_carbocycle','num_aromatic_heterocycle',
                    'num_aromatic_ring','num_atom_stereo_center',
                    'num_bridge_head_atom','num_HBA','num_HBD',
                    'num_heteroatoms','num_heterocycle',
                    'lipinski_HBA','lipinski_HBD','num_ring',
                    'num_rotatable_bond','num_saturated_carbocycle',
                    'num_saturated_heterocycle','num_saturated_ring',
                    'num_spiro_atoms','fraction_csp3','num_heavy_atom',
                    'num_NHOH','num_NO','mw','heavy_atom_mw',
                    'max_abs_estate_index','max_estate_index',
                    'min_abs_estate_index','min_estate_index',
                    'max_abs_partial_charge','max_partial_charge',
                    'min_abs_partial_charge','min_partial_charge',
                    'num_radical_electron','num_valence_electron',
                    'logP','MR','ASA','TPSA',
                    
                    'balabanJ','bertzCT','Ipc','chi0','chi0n','chi0v',
                    'chi1','chi1n','chi1v','chi2n','chi2v','chi3n','chi3v',
                    'chi4n','chi4v','hallkierAlpha','kapp1','kappa2','kappa3']

    desc_indices_estate_vsa = []
    for i in range(1,12):
        desc_indices_estate_vsa.append('estate_vsa'+str(i))

    desc_indices_vsa_estate = []
    for i in range(1,11):
        desc_indices_vsa_estate.append('vsa_estate'+str(i))
        
    desc_indices_peoe_vsa = []
    for i in range(1,15): desc_indices_peoe_vsa.append('peoe_vsa'+str(i))
    
    desc_indices_smr_vsa = []
    for i in range(1,11): desc_indices_smr_vsa.append('smr_vsa'+str(i))
    
    desc_indices_logp_vsa = []
    for i in range(1,13): desc_indices_logp_vsa.append('logp_vsa'+str(i))

    desc_indices_mqns = ['mqns'+str(i) for i in range(42)]

    total_desc_indices = desc_indices + desc_indices_estate_vsa + desc_indices_vsa_estate + desc_indices_peoe_vsa + desc_indices_smr_vsa + desc_indices_logp_vsa + desc_indices_mqns
    
    return pd.DataFrame(desc_2d, columns=total_desc_indices)

def estateFP_calc(sd):
    return pd.DataFrame([FingerprintMol(m)[-1] for m in sd])

def autocorr2D_calc(sd):
    return pd.DataFrame([rdMolDescriptors.CalcAUTOCORR2D(m) for m in sd])

def maccsFP_calc(sd):
    maccfps = [MACCSkeys.GenMACCSKeys(m) for m in sd]
    maccfps_bit = fingerprint_bitgenerator(maccfps)
    return pd.DataFrame(maccfps_bit)

def avalonFP_calc(sd):
    avalonfps = [pyAvalonTools.GetAvalonFP(m) for m in sd]
    avalonfps_bit = fingerprint_bitgenerator(avalonfps)
    return pd.DataFrame(avalonfps_bit)

def avalonCountFP_calc(sd):
    avaloncountfps = [pyAvalonTools.GetAvalonCountFP(m) for m in sd]
    avaloncountfps_bit = fingerprint_bitgenerator(avaloncountfps)
    return pd.DataFrame(avaloncountfps_bit)
        
def layerFP_calc(sd):
    layerfps = [Chem.rdmolops.LayeredFingerprint(m) for m in sd]
    layerfps_bit = fingerprint_bitgenerator(layerfps)
    return pd.DataFrame(layerfps_bit)

def morganFP_calc(sd, level):
    circularfps = [AllChem.GetMorganFingerprintAsBitVect(m,level) for m in sd]
    circularfps_bit = fingerprint_bitgenerator(circularfps)
    return pd.DataFrame(circularfps_bit)

def atompairFP_calc(sd):
    atompairfp = [Pairs.GetAtomPairFingerprint(m) for m in sd]
    atompairfp_key = atompairfingerprint_keyextractor(atompairfp)
    atompairfp_bit = atompair_bitgenerator(atompairfp, atompairfp_key)
    return pd.DataFrame(atompairfp_bit,columns=list(range(len(atompairfp_bit[0]))))

def torsionFP_calc(sd):
    torsionfp = [Torsions.GetTopologicalTorsionFingerprint(m) for m in sd]
    torsionfp_key = atompairfingerprint_keyextractor(torsionfp)
    torsionfp_bit = atompair_bitgenerator(torsionfp, torsionfp_key)
    return pd.DataFrame(torsionfp_bit,columns=list(range(len(torsionfp_bit[0]))))