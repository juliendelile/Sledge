
import sys
sys.path.append('/data/Work/Projets/Single_cells_GRN/synthetic_data')
import ND_lineage_trajectories_func as ndlt

ndlt.generateDataset(
    n_cells = 1000,
    n_lin_states = [10, 3],
    n_genes_per_lin_state = [800, 1500],
    n_cc_states = 4,
    n_genes_per_cc_phase = 1000,
    n_unexpressed_genes = 2002,
    p_branching=[.5, 0])


