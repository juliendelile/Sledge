
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from collections import OrderedDict, Counter
import igraph
import os
import pandas as pd


def generateStateLineageTree(n_states=20, p_branching=.5, seed=12345,
                             filename=None):

    np.random.seed(seed)

    tree = igraph.Graph(directed=True)
    tree.add_vertices(n_states)

    tree.add_edge(0, 1)

    if n_states > 2:

        available_nodes = range(2, n_states)

        while(True):

            current_leafs = np.where(np.asarray(tree.degree(mode='in')) - np.asarray(tree.degree(mode='out')) == 1)[0]

            for cl in current_leafs:

                branching = (np.random.rand() < p_branching)

                print "Current node ", cl, ['(branching)' if branching else ''][0]

                if branching and len(available_nodes) < 2:
                    branching = False

                for i in range(1 * branching + 1):
                    tree.add_edge(cl, available_nodes[0])
                    available_nodes.pop(0)

                print len(available_nodes)
                if len(available_nodes) == 0:
                    break
            else:
                continue
            break

    tree.vs["label"] = range(n_states)
    
    if filename != None:
        igraph.plot(tree, filename,
                    layout=tree.layout_reingold_tilford(root=[0]))

    return tree


def generateLineageData2(n_cells=300,
                         n_states=20, n_genes_per_state=100,
                         p_branching=.5, seed=12345,
                         filename=None,
                         basename='lin',
                         branch_ratio = .0,
                         n_genes_per_common_state=0,
                         num_common_state=2):

    np.random.seed(seed)

    lin_tree = generateStateLineageTree(n_states=n_states, p_branching=p_branching,
                                        seed=seed, filename=filename)

    # cell_state_density = 'uniform'

    nodes_depth = np.asarray([len(p) for p
                              in lin_tree.get_shortest_paths(
                                v=0,
                                to=range(n_states))]) - 1

    cells_edgeid = np.random.choice(range(n_states-1), size=n_cells, replace=True)

    edges = lin_tree.get_edgelist()

    cells_t = np.asarray([np.random.uniform(nodes_depth[edges[cells_edgeid[i]][0]], nodes_depth[edges[cells_edgeid[i]][1]]) for i in range(n_cells)])

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

    # igraph.plot(tree, 'artificial_lineage_all_cells.pdf',
    #             layout=tree.layout_reingold_tilford(root=[cells_by_branch[0][0]]),
    #             vertex_size=1, edge_arrow_size=0)


    t_lineage_distance = np.zeros((n_cells, n_cells))
    for v in range(n_cells):

        to = np.arange(v + 1, n_cells)

        for p in tree.get_shortest_paths(v=v, to=to, mode='All'):
            # print v, p[-1]
            t_lineage_distance[v, p[-1]] = sum([abs(t - s) for s, t in zip(cells_t[p], cells_t[p][1:])])    

    t_lineage_distance = t_lineage_distance + np.transpose(t_lineage_distance)

    # gaussKern = np.exp(-0.5/sigma**2 * t_lineage_distance * t_lineage_distance)

    # Cells branch by state
    cells_by_state = dict()
    for s in range(n_states):
        cells_by_state[s] = dict()

    for e_id in range(len(edges)):
        
        branch_time = cells_t[cells_by_branch[e_id]] - np.min(cells_t[cells_by_branch[e_id]]).astype(int)
        down = cells_by_branch[e_id][np.where(branch_time < .5)[0]]
        up = cells_by_branch[e_id][np.where(branch_time >= .5)[0]]

        # upstream part of the branch
        cells_by_state[edges[e_id][1]]['up'] = up

        # downstream part of the branch
        k = ['down2' if 'down1' in cells_by_state[edges[e_id][0]].keys() else 'down1'][0]
        cells_by_state[edges[e_id][0]][k] = down


    # Rule I: each state is associated with a unique set of genes
    data = np.zeros((n_states * n_genes_per_state, n_cells))
    sigma = 1
    gene_names = list()
    for s in range(n_states):

        s_cells = np.hstack(cells_by_state[s].values())
        
        if len(s_cells) > 0: # happens if many states and few cells

            g_cellmax = np.random.choice(s_cells, size=n_genes_per_state, replace=True)
        
            gene_names = gene_names + [basename+ str(s).zfill(len(str(n_states))) + 'g' + str(i).zfill(len(str(n_genes_per_state))) for i in range(n_genes_per_state)]

            for g in range(n_genes_per_state):

                data[(s * n_genes_per_state + g), :] = np.exp(-0.5/sigma**2 * t_lineage_distance[g_cellmax[g], :] * t_lineage_distance[g_cellmax[g], :])

    # Optional rule II: set of genes expressed commonly in different branches

    if branch_ratio > 0:
        # find leaves
        tree_leaves = np.where(np.asarray(lin_tree.degree(mode='out')) == 0)[0]
        # select the k longest branches, k given by ratio
        leaves_branches = lin_tree.get_shortest_paths(v=0, to=tree_leaves)
        leaves_branches_length = [[lin_tree.degree(mode='all')[i] for i in lb][::-1].index(3) for lb in leaves_branches]
        k = int(len(tree_leaves) * branch_ratio)
        used_branches_id = np.argsort(leaves_branches_length)[-k:]
        used_branches = [leaves_branches[i] for i in used_branches_id]
        
        data_com = np.zeros((num_common_state * n_genes_per_common_state, n_cells))

        for d in range(num_common_state):

            s_cells_ids = np.random.choice(range(1000), size=n_genes_per_common_state, replace=True) / 1000.0

            # ss = [ub[-(d+1)] for ub in used_branches]

            if d == 0:
                # s_cells = [cells_by_state[s]['up'] for s in ss]
                s_cells = [cells_by_state[ub[-(d+1)]]['up'] for ub in used_branches]
            else:
                # s_cells = [np.hstack([cells_by_state[s]['up'], cells_by_state[s]['down1']]) for s in ss]
                s_cells = [np.hstack([cells_by_state[ub[-(d+1)]]['up'], cells_by_state[ub[-(d+1)]]['down1']]) for ub in used_branches]
            
            if min([len(sc) for sc in s_cells]) > 0: # happens if many states and few cells

                # g_cellmax = np.random.choice(s_cells, size=n_genes_per_state, replace=True)
                
                g_cellmaxs = [s_c[np.floor(len(s_c) * s_cells_ids).astype(int)] for s_c in s_cells]
            
                gene_names = gene_names + ['com' + str(d).zfill(len(str(n_states))) + 'g' + str(i).zfill(len(str(n_genes_per_common_state))) for i in range(n_genes_per_common_state)]

                for g in range(n_genes_per_common_state):

                    data_com[(d * n_genes_per_common_state + g), :] = reduce(lambda x,y: x+y, [np.exp(-0.5/sigma**2 * t_lineage_distance[g_cellmax[g], :] * t_lineage_distance[g_cellmax[g], :]) for g_cellmax in g_cellmaxs])

        data = np.vstack((data, data_com))

    return {'data': data, 'state_tree': lin_tree, 'cells_edgeid': cells_edgeid, 'cells_t': cells_t, 'cells_by_state': cells_by_state, 'gene_names': gene_names}

# Keeping temporal distance strategy similar to the lineage case, if we generalize later for any graph (circle or tree)
def generateCellCycleData2(n_cells=171, n_cc_phases=7,
                          n_genes_per_phase=100, seed=12345):

    np.random.seed(seed)

    cc_adj = np.zeros((n_cc_phases, n_cc_phases))

    for cc_p in range(n_cc_phases - 1):
        cc_adj[cc_p, (cc_p + 1)] = 1
    cc_adj[-1, 0] = 1

    cc_graph = igraph.Graph.Adjacency((cc_adj > 0).tolist())
    cc_graph.vs['label'] = range(n_cc_phases)

    # igraph.plot(cc_graph, "cc_graph.pdf")


    cells_edgeid = np.random.choice(range(n_cc_phases), size=n_cells, replace=True)

    edges = cc_graph.get_edgelist()

    # here is the main difference with tree time
    # cells_t = np.asarray([np.random.uniform(nodes_depth[edges[cells_edgeid[i]][0]], nodes_depth[edges[cells_edgeid[i]][1]]) for i in range(n_cells)])

    cells_t = np.asarray([np.random.uniform(edges[cells_edgeid[i]][0], edges[cells_edgeid[i]][0]+1) for i in range(n_cells)])

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

    # # connect last cell from parent edge to first cell of child edge
    # for s in range(n_cc_phases):

    #     for s_in in cc_graph.neighbors(s, mode='Int'):

    #         parent_edge_id = edges.index((s_in, s))
    #         last_cell_id = cells_by_branch[parent_edge_id][-1]

    #         for s_out in cc_graph.neighbors(s, mode='Out'):

    #             child_edge_id = edges.index((s, s_out))
    #             first_cell_id = cells_by_branch[child_edge_id][0]

    #             cc_full_graph.add_edge(last_cell_id, first_cell_id)

    # igraph.plot(cc_full_graph, 'artificial_CC_all_cells.pdf',
    #             layout=tree.layout_reingold_tilford(root=[cells_by_branch[0][0]]),
    #             vertex_size=1, edge_arrow_size=0)

    # cc_distance = np.zeros((n_cells, n_cells))
    # for v in range(n_cells):

    #     to = np.arange(v + 1, n_cells)

    #     for p in cc_full_graph.get_shortest_paths(v=v, to=to, mode='All'):
    #         # print v, p[-1]
    #         cc_distance[v, p[-1]] = sum([abs(t - s) for s, t in zip(cells_t[p], cells_t[p][1:])])    

    # cc_distance = cc_distance + np.transpose(cc_distance)

    # simplified for speed sake
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
        
        branch_time = cells_t[cells_by_branch[e_id]] - np.min(cells_t[cells_by_branch[e_id]]).astype(int)
        down = cells_by_branch[e_id][np.where(branch_time < .5)[0]]
        up = cells_by_branch[e_id][np.where(branch_time >= .5)[0]]

        # upstream part of the branch
        cells_by_state[edges[e_id][1]]['up'] = up

        # downstream part of the branch
        k = ['down2' if 'down1' in cells_by_state[edges[e_id][0]].keys() else 'down1'][0]
        cells_by_state[edges[e_id][0]][k] = down


    # Hypothesis: each state is associated with a unique set of genes (like lineage)
    data = np.zeros((n_cc_phases * n_genes_per_phase, n_cells))
    sigma = .5
    gene_names = list()
    for s in range(n_cc_phases):

        s_cells = np.hstack(cells_by_state[s].values())
        
        g_cellmax = np.random.choice(s_cells, size=n_genes_per_phase, replace=True)
        
        gene_names = gene_names + ['cc'+ str(s).zfill(len(str(n_cc_phases))) + 'g' + str(i).zfill(len(str(n_genes_per_phase))) for i in range(n_genes_per_phase)]

        for g in range(n_genes_per_phase):

            data[(s * n_genes_per_phase + g), :] = np.exp(-0.5/sigma**2 * cc_distance[g_cellmax[g], :] * cc_distance[g_cellmax[g], :])

    return {'data':data, 'gene_names': gene_names, 'cells_by_state': cells_by_state, 'cc_graph': cc_graph, 'cells_t': cells_t}



def generateDataset(n_cells, n_lin_states, n_genes_per_lin_state, n_cc_states, n_genes_per_cc_phase, n_unexpressed_genes, p_branching, common_branch_ratio, n_genes_per_common_state, num_common_state):

    # import pdb; pdb.set_trace()

    lin_names = [str(nls)+'x'+str(n_genes_per_lin_state[id])+'p'+str(p_branching[id])+'lin'+str(id) for id, nls in enumerate(n_lin_states)]

    lin_names2 = [l + '+' + str(num_common_state) + 'x' + str(n_genes_per_common_state) + 'r' + str(common_branch_ratio[i]) if common_branch_ratio[i] > 0 else l for i,l in enumerate(lin_names)]

    directory_name = 'synthetic_data_' + "_".join(lin_names2) +  '_'+str(n_cc_states)+'x'+str(n_genes_per_cc_phase)+'cc_'+str(n_cells)+'cells_'+str(n_unexpressed_genes)+'unexp'

    if not os.path.exists(directory_name):
        os.mkdir(directory_name, 0755)
    
    import random
    r = lambda: random.randint(0,255)
    
    scolors = ['#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999','#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3'] + ['#%02X%02X%02X' % (r(),r(),r()) for i in range(200)] #colors from R colorbrewer

        # c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"))

    # Generate Lineage data
    #######################

    lin_list = list()

    for id, nls in enumerate(n_lin_states):

        lin_list.append(generateLineageData2(n_cells=n_cells, n_states=nls, n_genes_per_state=n_genes_per_lin_state[id], p_branching=p_branching[id], branch_ratio=common_branch_ratio[id], n_genes_per_common_state=n_genes_per_common_state, num_common_state=num_common_state, seed=1234, basename='l'+'iouea'[id]+'n'))

        igraph.plot(lin_list[id]['state_tree'], directory_name+"/artificial_lineage_"+str(id)+".pdf",
                    layout=lin_list[id]['state_tree'].layout_reingold_tilford(root=[0]),
                    vertex_color=[scolors[i] for i in range(nls)]
                    )

    lin_state_status = list()

    for id, lin in enumerate(lin_list):

        # state status for each cell
        lin_state_status.append(np.zeros((n_cells)).astype(int))
        for k, v in lin['cells_by_state'].iteritems():
            lin_state_status[id][np.hstack(v.values())] = k

        # order by 1. ordered edge 2. time
        lin_tree_leaves = np.where(np.asarray(lin['state_tree'].degree(mode='in')) - np.asarray(lin['state_tree'].degree(mode='out')) == 1)[0]
        all_paths_vertices = lin['state_tree'].get_all_shortest_paths(0, to=lin_tree_leaves)
        all_paths_edges = [zip(p, p[1:]) for p in all_paths_vertices]

        edges_ordered = list(OrderedDict.fromkeys(sum(all_paths_edges, [])))
        edges = lin['state_tree'].get_edgelist()
        edges_ordered_idx = [edges.index(e) for e in edges_ordered]
        cells_by_edges = np.asarray(edges_ordered_idx).argsort()[lin['cells_edgeid']]
        cells_order = np.lexsort((lin['cells_t'], cells_by_edges))

        if id==0:
            cells_order_global = cells_order
        
        # order by time only
        # cells_order = lin['cells_t'].argsort()

        # branch status for each cell
        branch_status = np.zeros((n_cells))
        for k1, v1 in lin['cells_by_state'].iteritems():
            for k2, v2 in v1.iteritems():
                branch_status[v2] = k1 + [.0 if k2=='up' else .33 if k2=='down1' else .66][0]

        # label each branch at the middle cell position
        labels = list()
        labels_pos = list()
        for k1, v1 in lin['cells_by_state'].iteritems():
            for k2, v2 in v1.iteritems():
                if len(v2) > 0:
                    labels.append(str(k1) + k2)
                    labels_pos.append( cells_order.argsort()[v2[int(.5*len(v2))]] )

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(np.vstack((branch_status/float(n_lin_states[id]), lin['data']))[:, cells_order], aspect='auto', interpolation=None)
        plt.xticks(labels_pos, labels, rotation='vertical')
        plt.savefig(directory_name+'/artificial_data_lin'+str(id)+'.pdf')



    # Generate Cell Cycle data
    ##########################

    # state status for each cell
    if n_cc_states > 0 and n_genes_per_cc_phase > 0:

        cc2 = generateCellCycleData2(n_cells=n_cells, n_cc_phases=n_cc_states,
                                           n_genes_per_phase=n_genes_per_cc_phase, seed=12345)

        cc_state_status = np.zeros((n_cells)).astype(int)
        for k, v in cc2['cells_by_state'].iteritems():
            cc_state_status[np.hstack(v.values())] = k

        igraph.plot(cc2['cc_graph'], directory_name+"/cc_graph.pdf", vertex_color=[scolors[i] for i in range(n_cc_states)])

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(cc2['data'][:, cc2['cells_t'].argsort()], aspect='auto', interpolation=None)
        # plt.show()
        plt.savefig(directory_name+'/cc_data.pdf')
    else:
        cc2 = {'gene_names':[]}
        # {'data':data, 'gene_names': gene_names, 'cells_by_state': cells_by_state, 'cc_graph': cc_graph, 'cells_t': cells_t}

    # Merge data
    ############

    data = lin_list[0]['data']

    # shuffle other "lin" data
    if len(lin_list) > 1:
        for id in range(1, len(lin_list)):
            shuffling = np.random.permutation(range(n_cells))
            lin_data_perm = lin_list[id]['data'][:, shuffling]
            lin_state_status[id] = lin_state_status[id][shuffling]
            data = np.vstack((data, lin_data_perm))

    # shuffle CC data
    if n_cc_states > 0 and n_genes_per_cc_phase > 0:
        shuffling = np.random.permutation(range(n_cells))
        cc_data_perm = cc2['data'][:, shuffling]
        cc_state_status = cc_state_status[shuffling]

        data = np.vstack((data, cc_data_perm))

    # Add unexpressed genes
    unexpressed_names = ['unexp'+str(i).zfill(5) for i in range(n_unexpressed_genes)]
    data = np.vstack((data, np.zeros((n_unexpressed_genes, n_cells))))

    # Add technical noise
    #####################

    noise_intensity = .5
    data_noise = data + noise_intensity * np.random.randn(data.shape[0], data.shape[1])
    data_noise[data_noise < 0] = 0

    # Plot final data
    #################

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.imshow(data_noise[:, cells_order_global], aspect='auto', interpolation=None)
    # plt.show()
    plt.savefig(directory_name+'/artificial_data.pdf')

    # Export as ExpressionSet 
    #########################

    cell_names = ['cell' + str(i).zfill(5) for i in range(n_cells)]
    gene_order = np.random.permutation(range(data_noise.shape[0])) # shuffle gene order

    all_gene_names = sum([sum([lin_list[i]['gene_names'] for i in range(len(lin_list))], []), cc2['gene_names'], unexpressed_names], [])

    assay_df = pd.DataFrame(data_noise[gene_order, :], index=[all_gene_names[i] for i in gene_order], columns=cell_names)

    assay_df.to_csv(directory_name + '/assayData.csv', index=True, header=True, sep='\t')

    tp = lin_list[0]['cells_t'].astype(int) + 1
    pheno_df = pd.DataFrame({
                        'timepoint': tp,
                        'replicate_id': [1 for i in range(n_cells)],
                        'treatment': ['_'.join('N' for i in range(int(tp[j]))) for j in range(n_cells)]
        }, index=cell_names)

    if n_cc_states > 0 and n_genes_per_cc_phase > 0:
        pheno_df['cc_state'] = cc_state_status

    for id in range(len(lin_list)):
        pheno_df['l'+'iouea'[id]+'n_state'] = lin_state_status[id]


    pheno_df.to_csv(directory_name + '/phenoData.csv', index=True, header=True, sep='\t')
    # pheno_df.loc[:, ['timepoint', 'replicate_id', 'treatment', 'lin_state', 'cc_state']].to_csv(directory_name + '/phenoData.csv', index=True, header=True, sep='\t')

    # Plot output dataset
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.imshow(data_noise[gene_order, :], aspect='auto', interpolation=None)
    # plt.show()
    plt.savefig(directory_name+'/artificial_data_shuffled.pdf')
