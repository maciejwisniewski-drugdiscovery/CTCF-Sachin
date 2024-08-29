import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pytraj as pt

def load_trajectory(trajectory_filepath, topology_filepath=None):
    trajectory_type = trajectory_filepath.split('.')[-1]
    if trajectory_type == 'xtc':
        trajectory = pt.load(trajectory_filepath, topology_filepath)
    elif trajectory_type == 'pdb':
        trajectory = pt.load(trajectory_filepath)
    return trajectory
def rigid_transform_Kabsch_3D_numpy(A, B, include_reflections=True):
    # Ensure A and B are 3xN pointclouds
    if A.shape[0] != 3:
        raise ValueError(f"matrix A is not 3xN, it is {A.shape[0]}x{A.shape[1]}")
    if B.shape[0] != 3:
        raise ValueError(f"matrix B is not 3xN, it is {B.shape[0]}x{B.shape[1]}")

    # Compute centroids of A and B
    centroid_A = np.mean(A, axis=1, keepdims=True)
    centroid_B = np.mean(B, axis=1, keepdims=True)

    # Center the point clouds
    Am = A - centroid_A
    Bm = B - centroid_B

    # Compute the covariance matrix H
    H = Am @ Bm.T

    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation matrix R
    R = Vt.T @ U.T

    # Handle reflection case
    if include_reflections and np.linalg.det(R) < 0:
        SS = np.diag([1, 1, -1])
        R = (Vt.T @ SS) @ U.T

    # Compute translation
    t = -R @ centroid_A + centroid_B

    return R, t.flatten()
def split_trajectory(trajectory,start=None,end=None,strand=None,only_CA=True):

    protein_mask = '1'
    if only_CA == False:
        protein_mask = '1'
    ligand_mask = '2'

    protein_trajectory = trajectory[protein_mask][start:end:strand]
    ligand_trajectory = trajectory[ligand_mask][start:end:strand]

    return protein_trajectory, ligand_trajectory
def calculate_pairwise_RMSD(trajectory,reference_trajectory_snapshot):
    return pairwise_rmsd_data
def calculate_RMSD(trajectory):
    return rmsd_data
def calculate_RMSF(trajectory):
    return rmsf_data
def pairwise_RMSD_heatmap():
    return None
def rmsd_plot():
    return None
def rmsf_plot():
    return None

