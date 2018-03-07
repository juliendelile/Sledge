
import sys
sys.path.append('/data/Work/Projets/Single_cells_GRN/synthetic_data')
import ND_lineage_trajectories_func as ndlt
import git

# git commit version
repo = git.Repo()

if 'ND_lineage_trajectories_func.py' in [ item.a_path for item in repo.index.diff(None) ]:
	print "commit changes first"
	quit()

commit_hash = repo.head.object.hexsha[0:7]


# Lineagity test: Lin Lon CC
ndlt.generateDataset(
    n_cells = 1000,
    n_lin_states = [10,3],
    n_genes_per_lin_state = [800, 1500],
    n_cc_states = 4,
    n_genes_per_cc_phase = 1000,
    n_unexpressed_genes = 2000,
    p_branching=[0.5, 0],
    common_branch_ratio=[0, 0],
    n_genes_per_common_state=0,
    num_common_state=0,
    noise_intensity=.5,
    commit_hash=commit_hash,
    seed=1234
    )

# # Multiscale test: Common genes setup
# ndlt.generateDataset(
#     n_cells = 1000,
#     n_lin_states = [10],
#     n_genes_per_lin_state = [800],
#     n_cc_states = 0,
#     n_genes_per_cc_phase = 0,
#     n_unexpressed_genes = 2000,
#     p_branching=[.5],
#     common_branch_ratio=[.67],
#     n_genes_per_common_state=6000,
#     num_common_state=2,
#     noise_intensity=.5,
#     commit_hash=commit_hash
#     )

# # Multiscale test: complex lineage
# ndlt.generateDataset(
#     n_cells = 1000,
#     n_lin_states = [30],
#     n_genes_per_lin_state = [300],
#     n_cc_states = 0,
#     n_genes_per_cc_phase = 0,
#     n_unexpressed_genes = 2000,
#     p_branching=[.5],
#     common_branch_ratio=[0],
#     n_genes_per_common_state=0,
#     num_common_state=2,
#     noise_intensity=.75,
#     commit_hash=commit_hash,
#     seed=123
#     )

# # Smoothing test: complex lineage with more noise
# ndlt.generateDataset(
#     n_cells = 1000,
#     n_lin_states = [30],
#     n_genes_per_lin_state = [300],
#     n_cc_states = 0,
#     n_genes_per_cc_phase = 0,
#     n_unexpressed_genes = 2000,
#     p_branching=[.5],
#     common_branch_ratio=[0],
#     n_genes_per_common_state=0,
#     num_common_state=2,
#     noise_intensity=1,
#     commit_hash=commit_hash,
#     seed=123
#     )

