# Test dataset generation

import sledge

if __name__ == "__main__":

    sledge.generate_dataset(
        n_cells                  = 500,
        n_tree_states            = [21, 3],
        p_branching              = [0.3, 0],
        n_genes_per_tree_state   = [20, 50],
        n_cc_states              = 4,
        n_genes_per_cc_state     = 50,
        n_unexpressed_genes      = 1000,
        noise_intensity          = 0.3,
        seed                     = 123,
        plot_size_factor         = 1)

# END