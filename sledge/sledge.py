# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import igraph
import os
import pandas as pd
import json
import copy
import random


def generate_dataset(
        n_cells=1000,
        n_tree_states=[20, 3],
        p_branching=[0.3, 0],
        n_genes_per_tree_state=[200, 500],
        n_cc_states=4,
        n_genes_per_cc_state=500,
        n_unexpressed_genes=1000,
        noise_intensity=0.3,
        seed=123,
        plot_size_factor=1):
    """Generate a synthetic RNA-seq dataset.

    This function generates gene expression profiles according to the
    following assumptions:

    * The genes are grouped into independent sets of modules.
    * A module contains genes upregulated in a specific cell "state".
    * Each set of gene modules is associated with a specific type of trajectories.
    * The trajectories can be either:
        * linear: gene expression is cascading in-between a linear sequence of
        transitions between "cell states" (e.g. A->B, B->C between cell states
        A, B and C).
        * bifurcating: a state can be followed by two states (e.g. A->B, A->C).
        * cycling: the genes are expressed along a cycle of states (A->B, B->C, C->A)
    * The depth of the trajectory associated with the first set of modules
    determines the temporal coordinates of the cells.
    * Each gene profile follow a normal distribution along its trajectory.
    * A gaussian noise term is added to gene expression.
    * A set of unexpressed genes containing only gaussian noise can be integrated.

    Args:
        n_cells (int): the number of cells.
        n_tree_states (list of int): each integer determines the number of states
        composing a specific non-cycling trajectory. 
        n_genes_per_tree_state (list of int): the number of genes associated with
        each state of a non-cycling trajectory.
        n_cc_states (int): the number of cell cycle states.
        n_genes_per_cc_state (int): the number of genes associated with each
        state of the cell cycle.
        n_unexpressed_genes (int): the number of unexpressed genes.
        p_branching (list of float): the probability that a (non-cyclic) state is
        bifurcating. The sequence of transition is linear if equal to 0.
        noise_intensity (float): weight scaling the amplitude of the gaussian
        noise added to the gene expression levels.
        seed (int): seed of the random number generator.
        plot_size_factor (float): scaling factor of the labels in the plots.

    Returns:
        None. A set of csv files and plots are produced in the working
            directory.
    """

    # store arguments for writing json below
    argmts = locals()
    args_clean = copy.deepcopy(argmts)

    lin_names = [str(nls)+'x'+str(n_genes_per_tree_state[id])+'p' +
                 str(p_branching[id])+'lin'+str(id) for id, nls in enumerate(n_tree_states)]

    directory_name = 'synthetic_data_' + "_".join(lin_names) + '_'+str(n_cc_states)+'x'+str(
        n_genes_per_cc_state)+'cc_'+str(n_cells)+'cells_'+str(n_unexpressed_genes)+'unexp_'+str(noise_intensity)+'noise'

    if not os.path.exists(directory_name):
        os.mkdir(directory_name, 0o755)

    with open(directory_name+'/settings.json', 'w') as outfile:
        json.dump(args_clean, outfile, indent=4)

    def r(): return random.randint(0, 255)

    # colors from R colorbrewer + 200 random colors
    # c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"))
    scolors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999', '#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#EBFF00', '#8F00FF', '#14FF00', '#00FFFF', '#24FF00', '#FF00E6', '#FF000F', '#FF3D00', '#3300FF', '#00B2FF', '#1400FF', '#0075FF', '#FF7A00', '#9E00FF', '#FF8A00', '#CCFF00', '#0066FF', '#FF0F00', '#FF9900', '#05FF00', '#FF004C', '#FFF500', '#FF005C', '#0047FF', '#FFD600', '#00A3FF', '#2400FF', '#FF00D6', '#00FF57', '#00FF47', '#FF001F', '#FF00B8', '#ADFF00', '#00FF75', '#70FF00', '#FF0099', '#00FFA3', '#00FF29', '#FFE500', '#0500FF', '#FF6B00', '#00E0FF',
               '#00FF1A', '#0085FF', '#00FFF0', '#FF007A', '#BD00FF', '#0057FF', '#8000FF', '#00FFE0', '#00FF85', '#FF00C7', '#FFA800', '#0094FF', '#FFC700', '#DBFF00', '#EB00FF', '#00FFC2', '#BDFF00', '#33FF00', '#0038FF', '#FF0000', '#FA00FF', '#FF003D', '#00C2FF', '#4200FF', '#52FF00', '#00FFD1', '#9EFF00', '#00FF94', '#8FFF00', '#FF5C00', '#FFB800', '#FF002E', '#000AFF', '#FF00A8', '#CC00FF', '#00FF38', '#00FF66', '#00FFB3', '#00FF0A', '#FF00F5', '#00D1FF', '#00F0FF', '#5200FF', '#0019FF', '#80FF00', '#7000FF', '#0029FF', '#FF1F00', '#FF008A', '#FF4D00', '#AD00FF', '#FAFF00', '#DB00FF', '#42FF00', '#6100FF', '#FF006B', '#FF2E00', '#61FF00'] + ['#%02X%02X%02X' % (r(), r(), r()) for i in range(200)]

    # Generate Lineage data
    #######################

    lin_list = list()

    for id, nls in enumerate(n_tree_states):

        lin_list.append(generate_lineage_data(
            n_cells=n_cells,
            n_states=nls,
            n_genes_per_state=n_genes_per_tree_state[id],
            p_branching=p_branching[id],
            seed=seed,
            basename='l'+'iouea'[id]+'n'))

        igraph.plot(
            lin_list[id]['state_tree'],
            directory_name + "/artificial_lineage_"+str(id)+".pdf",
            layout=lin_list[id]['state_tree'].layout_reingold_tilford(root=[0]),
            vertex_color=[scolors[i] for i in range(nls)],
            vertex_size=plot_size_factor*30,
            vertex_label_size=plot_size_factor*15,
            edge_width=3,
            margin=100)

    lin_state_status = list()

    for id, lin in enumerate(lin_list):

        # state status for each cell
        lin_state_status.append(np.zeros((n_cells)).astype(int))
        for k, v in lin['cells_by_state'].items():
            lin_state_status[id][np.hstack(list(v.values()))] = k

        # order by 1. ordered edge 2. time
        lin_tree_leaves = np.where(np.asarray(lin['state_tree'].degree(
            mode='in')) - np.asarray(lin['state_tree'].degree(mode='out')) == 1)[0]
        all_paths_vertices = lin['state_tree'].get_all_shortest_paths(
            0, to=lin_tree_leaves)
        all_paths_edges = [list(zip(p, p[1:])) for p in all_paths_vertices]

        edges_ordered = list(OrderedDict.fromkeys(sum(all_paths_edges, [])))
        edges = lin['state_tree'].get_edgelist()
        edges_ordered_idx = [edges.index(e) for e in edges_ordered]
        cells_by_edges = np.asarray(edges_ordered_idx).argsort()[
            lin['cells_edgeid']]
        cells_order = np.lexsort((lin['cells_t'], cells_by_edges))

        if id == 0:
            cells_order_global = cells_order

            # Plot lineage distance
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
            plt.matshow(lin_list[id]['ldistance'][cells_order, :][:, cells_order],
                        aspect='auto', origin='upper', cmap=plt.cm.Oranges)
            plt.savefig(directory_name+"/distance_lineage_" +
                        str(id)+".pdf", bbox_inches='tight')
            plt.close(fig)

        # branch status for each cell
        branch_status = np.zeros((n_cells))
        for k1, v1 in lin['cells_by_state'].items():
            for k2, v2 in v1.items():
                branch_status[v2] = k1 + [.0 if k2 ==
                                          'up' else .33 if k2 == 'down1' else .66][0]

        # label each branch at the middle cell position
        labels = list()
        labels_pos = list()
        for k1, v1 in lin['cells_by_state'].items():
            for k2, v2 in v1.items():
                if len(v2) > 0:
                    labels.append(str(k1) + k2)
                    labels_pos.append(cells_order.argsort()
                                      [v2[int(.5*len(v2))]])

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(np.vstack((branch_status/float(n_tree_states[id]), lin['data']))[
                   :, cells_order], aspect='auto', interpolation=None)
        plt.xticks(labels_pos, labels, rotation='vertical')
        plt.savefig(directory_name+'/artificial_data_lin'+str(id)+'.pdf')
        plt.close(fig)

    # Generate Cell Cycle data
    ##########################

    # state status for each cell
    if n_cc_states > 0 and n_genes_per_cc_state > 0:

        cc2 = generate_cell_cycle_data(n_cells=n_cells, n_cc_phases=n_cc_states,
                                       n_genes_per_phase=n_genes_per_cc_state, seed=seed)

        cc_state_status = np.zeros((n_cells)).astype(int)
        for k, v in cc2['cells_by_state'].items():
            cc_state_status[np.hstack(list(v.values()))] = k

        igraph.plot(cc2['cc_graph'], directory_name+"/cc_graph.pdf", vertex_color=[scolors[i] for i in range(n_cc_states)],
                    vertex_size=60,
                    vertex_label_size=30,
                    edge_width=3,
                    margin=100
                    )

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(cc2['data'][:, cc2['cells_t'].argsort()],
                   aspect='auto', interpolation=None)
        # plt.show()
        plt.savefig(directory_name+'/cc_data.pdf')
        plt.close(fig)

    else:
        cc2 = {'gene_names': []}

    # Merge data
    ############

    data = lin_list[0]['data']

    # shuffle other "lin" data
    if len(lin_list) > 1:
        for id in range(1, len(lin_list)):
            shuffling = np.random.permutation(list(range(n_cells)))
            lin_data_perm = lin_list[id]['data'][:, shuffling]
            lin_state_status[id] = lin_state_status[id][shuffling]
            data = np.vstack((data, lin_data_perm))

    # shuffle CC data
    if n_cc_states > 0 and n_genes_per_cc_state > 0:
        shuffling = np.random.permutation(list(range(n_cells)))
        cc_data_perm = cc2['data'][:, shuffling]
        cc_state_status = cc_state_status[shuffling]

        data = np.vstack((data, cc_data_perm))

    # Add unexpressed genes
    unexpressed_names = ['unexp'+str(i).zfill(5)
                         for i in range(n_unexpressed_genes)]
    data = np.vstack((data, np.zeros((n_unexpressed_genes, n_cells))))

    # Add technical noise
    #####################

    data_noise = data + noise_intensity * \
        np.random.randn(data.shape[0], data.shape[1])
    data_noise[data_noise < 0] = 0

    # Cast to int
    #############

    data_int = (10 * data_noise / data_noise.max()).astype(int)

    # Plot final data
    #################

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.imshow(data_noise[:, cells_order_global],
               aspect='auto', interpolation=None)
    # plt.show()
    plt.savefig(directory_name+'/artificial_data.pdf')
    plt.close(fig)

    # Export as ExpressionSet
    #########################

    cell_names = ['cell' + str(i).zfill(5) for i in range(n_cells)]
    gene_order = np.random.permutation(
        list(range(data_int.shape[0])))  # shuffle gene order

    all_gene_names = sum([sum([lin_list[i]['gene_names'] for i in range(
        len(lin_list))], []), cc2['gene_names'], unexpressed_names], [])

    assay_df = pd.DataFrame(
        data_int[gene_order, :],
        index=[all_gene_names[i] for i in gene_order],
        columns=cell_names)

    assay_df.to_csv(directory_name + '/assayData.csv',
                    index=True, header=True, sep='\t')

    tp = lin_list[0]['cells_t'].astype(int) + 1
    pheno_df = pd.DataFrame({
        'timepoint': tp,
        'replicate_id': [1 for i in range(n_cells)],
        'treatment': ['_'.join('N' for i in range(int(tp[j]))) for j in range(n_cells)]
    }, index=cell_names)

    if n_cc_states > 0 and n_genes_per_cc_state > 0:
        pheno_df['cc_state'] = cc_state_status

    for id in range(len(lin_list)):
        pheno_df['l'+'iouea'[id]+'n_state'] = lin_state_status[id]

    pheno_df.to_csv(directory_name + '/phenoData.csv',
                    index=True, header=True, sep='\t')

    # Plot output dataset
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.imshow(data_int[gene_order, :], aspect='auto', interpolation=None)
    # plt.show()
    plt.savefig(directory_name+'/artificial_data_shuffled.pdf')
    plt.close(fig)


def generate_lineage_tree(
        n_states=20,
        p_branching=.5,
        seed=12345,
        filename=None):

    np.random.seed(seed)

    tree = igraph.Graph(directed=True)
    tree.add_vertices(n_states)

    tree.add_edge(0, 1)

    if n_states > 2:

        available_nodes = list(range(2, n_states))

        while(True):

            current_leafs = np.where(np.asarray(tree.degree(
                mode='in')) - np.asarray(tree.degree(mode='out')) == 1)[0]

            for cl in current_leafs:

                branching = (np.random.rand() < p_branching)

                # print("Current node ", cl, ['(branching)' if branching else ''][0])

                if branching and len(available_nodes) < 2:
                    branching = False

                for i in range(1 * branching + 1):
                    tree.add_edge(cl, available_nodes[0])
                    available_nodes.pop(0)

                if len(available_nodes) == 0:
                    break
            else:
                continue
            break

    tree.vs["label"] = list(range(n_states))

    if filename != None:
        igraph.plot(tree, filename,
                    layout=tree.layout_reingold_tilford(root=[0]))

    return tree


def generate_lineage_data(
        n_cells=300,
        n_states=20,
        n_genes_per_state=100,
        p_branching=.5,
        seed=12345,
        filename=None,
        basename='lin'):

    np.random.seed(seed)

    lin_tree = generate_lineage_tree(n_states=n_states, p_branching=p_branching,
                                     seed=seed, filename=filename)

    nodes_depth = np.asarray([len(p) for p
                              in lin_tree.get_shortest_paths(
        v=0,
        to=list(range(n_states)))]) - 1

    cells_edgeid = np.random.choice(
        list(range(n_states-1)), size=n_cells, replace=True)

    edges = lin_tree.get_edgelist()

    cells_t = np.asarray([np.random.uniform(nodes_depth[edges[cells_edgeid[i]][0]],
                                            nodes_depth[edges[cells_edgeid[i]][1]]) for i in range(n_cells)])

    # lineage distance: temporal distance on lineage tree
    # ###################################################
    tree = igraph.Graph(directed=True)
    tree.add_vertices(n_cells)

    # connect temporally ordered cells in each branch
    cells_by_branch = dict()
    for e_id in range(len(edges)):
        c_ids = np.where(cells_edgeid == e_id)[0]
        c_ids_ordered = c_ids[cells_t[c_ids].argsort()]
        for i in range(len(c_ids) - 1):
            tree.add_edge(c_ids_ordered[i], c_ids_ordered[i + 1])
        cells_by_branch[e_id] = c_ids_ordered

    # connect last cell from parent edge to first cell of child edge
    for s in range(n_states):

        for s_in in lin_tree.neighbors(s, mode='Int'):

            parent_edge_id = edges.index((s_in, s))
            last_cell_id = cells_by_branch[parent_edge_id][-1]

            for s_out in lin_tree.neighbors(s, mode='Out'):

                child_edge_id = edges.index((s, s_out))
                first_cell_id = cells_by_branch[child_edge_id][0]

                tree.add_edge(last_cell_id, first_cell_id)

    t_lineage_distance = np.zeros((n_cells, n_cells))
    for v in range(n_cells):

        to = np.arange(v + 1, n_cells)

        for p in tree.get_shortest_paths(v=v, to=to, mode='All'):

            t_lineage_distance[v, p[-1]] = sum([abs(t - s)
                                                for s, t in list(zip(cells_t[p], cells_t[p][1:]))])

    t_lineage_distance = t_lineage_distance + np.transpose(t_lineage_distance)

    # Cells branch by state
    cells_by_state = dict()
    for s in range(n_states):
        cells_by_state[s] = dict()

    for e_id in range(len(edges)):

        branch_time = cells_t[cells_by_branch[e_id]] - \
            np.min(cells_t[cells_by_branch[e_id]]).astype(int)
        down = cells_by_branch[e_id][np.where(branch_time < .5)[0]]
        up = cells_by_branch[e_id][np.where(branch_time >= .5)[0]]

        # upstream part of the branch
        cells_by_state[edges[e_id][1]]['up'] = up

        # downstream part of the branch
        k = ['down2' if 'down1' in cells_by_state[edges[e_id][0]].keys()
             else 'down1'][0]
        cells_by_state[edges[e_id][0]][k] = down

    # Rule I: each state is associated with a unique set of genes
    data = np.zeros((n_states * n_genes_per_state, n_cells))
    sigma = 1
    gene_names = list()
    for s in range(n_states):

        s_cells = np.hstack(list(cells_by_state[s].values()))

        if len(s_cells) > 0:  # happens if many states and few cells

            g_cellmax = np.random.choice(
                s_cells, size=n_genes_per_state, replace=True)

            gene_names = gene_names + [basename + str(s).zfill(len(str(n_states))) + 'g' + str(
                i).zfill(len(str(n_genes_per_state))) for i in range(n_genes_per_state)]

            for g in range(n_genes_per_state):

                data[(s * n_genes_per_state + g), :] = np.exp(-0.5/sigma**2 *
                                                              t_lineage_distance[g_cellmax[g], :] * t_lineage_distance[g_cellmax[g], :])

    return {
        'data': data,
        'state_tree': lin_tree,
        'cells_edgeid': cells_edgeid,
        'cells_t': cells_t,
        'cells_by_state': cells_by_state,
        'gene_names': gene_names,
        'ldistance': t_lineage_distance}


def generate_cell_cycle_data(
        n_cells=100,
        n_cc_phases=7,
        n_genes_per_phase=100,
        seed=12345):

    np.random.seed(seed)

    cc_adj = np.zeros((n_cc_phases, n_cc_phases))

    for cc_p in range(n_cc_phases - 1):
        cc_adj[cc_p, (cc_p + 1)] = 1
    cc_adj[-1, 0] = 1

    cc_graph = igraph.Graph.Adjacency((cc_adj > 0).tolist())
    cc_graph.vs['label'] = list(range(n_cc_phases))

    cells_edgeid = np.random.choice(
        list(range(n_cc_phases)), size=n_cells, replace=True)

    edges = cc_graph.get_edgelist()

    # here is the main difference with tree time
    # cells_t = np.asarray([np.random.uniform(nodes_depth[edges[cells_edgeid[i]][0]], nodes_depth[edges[cells_edgeid[i]][1]]) for i in range(n_cells)])

    cells_t = np.asarray([np.random.uniform(
        edges[cells_edgeid[i]][0], edges[cells_edgeid[i]][0]+1) for i in range(n_cells)])

    # # lineage distance: temporal distance on lineage tree
    # # ###################################################
    cc_full_graph = igraph.Graph(directed=True)
    cc_full_graph.add_vertices(n_cells)

    # connect temporally ordered cells in each branch
    cells_by_branch = dict()
    for e_id in range(len(edges)):
        c_ids = np.where(cells_edgeid == e_id)[0]
        c_ids_ordered = c_ids[cells_t[c_ids].argsort()]
        for i in range(len(c_ids) - 1):
            cc_full_graph.add_edge(c_ids_ordered[i], c_ids_ordered[i + 1])
        cells_by_branch[e_id] = c_ids_ordered

    cc_distance = np.zeros((n_cells, n_cells))
    for i in range(n_cells):
        for j in range(n_cells):
            if i < j:
                d = abs(cells_t[i] - cells_t[j])
                if d >= .5 * n_cc_phases:
                    d = n_cc_phases - d
                cc_distance[i, j] = d

    cc_distance = cc_distance + np.transpose(cc_distance)

    # Cells branch by state (like lineage)
    cells_by_state = dict()
    for s in range(n_cc_phases):
        cells_by_state[s] = dict()

    for e_id in range(len(edges)):

        branch_time = cells_t[cells_by_branch[e_id]] - \
            np.min(cells_t[cells_by_branch[e_id]]).astype(int)
        down = cells_by_branch[e_id][np.where(branch_time < .5)[0]]
        up = cells_by_branch[e_id][np.where(branch_time >= .5)[0]]

        # upstream part of the branch
        cells_by_state[edges[e_id][1]]['up'] = up

        # downstream part of the branch
        k = ['down2' if 'down1' in cells_by_state[edges[e_id][0]].keys()
             else 'down1'][0]
        cells_by_state[edges[e_id][0]][k] = down

    # Hypothesis: each state is associated with a unique set of genes (like lineage)
    data = np.zeros((n_cc_phases * n_genes_per_phase, n_cells))
    sigma = .5
    gene_names = list()
    for s in range(n_cc_phases):

        s_cells = np.hstack(list(cells_by_state[s].values()))

        g_cellmax = np.random.choice(
            s_cells, size=n_genes_per_phase, replace=True)

        gene_names = gene_names + ['cc' + str(s).zfill(len(str(n_cc_phases))) + 'g' + str(
            i).zfill(len(str(n_genes_per_phase))) for i in range(n_genes_per_phase)]

        for g in range(n_genes_per_phase):

            data[(s * n_genes_per_phase + g), :] = np.exp(-0.5/sigma**2 *
                                                          cc_distance[g_cellmax[g], :] * cc_distance[g_cellmax[g], :])

    return {'data': data, 'gene_names': gene_names, 'cells_by_state': cells_by_state, 'cc_graph': cc_graph, 'cells_t': cells_t}
