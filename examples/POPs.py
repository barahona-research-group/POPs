from Bio.PDB import *
from collections import defaultdict
import networkx as nx
import pandas as pd
import numpy as np
import os
import shutil
from pandas import DataFrame
from collections import defaultdict
import math
from operator import itemgetter
import plotly.figure_factory as ff
import plotly.graph_objects as go

class PropOptPaths():
    def __init__(self, pdb_file, source_lig, allo_lig, graph_bond, bond_propensity):
        self.pdb_file = pdb_file
        residue_list = self.pdb_residue_list()
        self.res_list_with_lig = residue_list
        self.source = self.site_format(source_lig)
        self.allo_lig = self.site_format(allo_lig)
        self.res_list_no_lig = []
        for residue in residue_list:
            if residue.split(' ')[0] + residue.split(' ')[1] + ' ' + residue.split(' ')[2] not in self.source and residue.split(' ')[0] + residue.split(' ')[1] + ' ' + residue.split(' ')[2] not in self.allo_lig:
                self.res_list_no_lig.append(residue)

        self.bond = pd.read_csv(graph_bond, dtype = str)
        self.bond_propensity = pd.read_csv(bond_propensity)
        self.all_bonds_qs = self.bond_propensity[['bond_name', 'qs']].values
        neighbour_res_summary = self.neighbour_res_summary()
        self.cov_bond_res_summary = neighbour_res_summary
        self.all_weak_edges = self.weak_edges_list()

    # To extract residues from a pdb file and put them in the format of eg. 'GLU 170 A' etc.
    def pdb_residue_list(self):
        parser = PDBParser()
        structure = parser.get_structure('protein', self.pdb_file)
        residues = structure.get_residues()
        res_list = []
        for res in residues:
            if str(res).split('=')[3][0] == ' ':
                res_list.append(str(res).split('=')[0].split(' ')[-2] + ' ' 
                                + str(res).split('=')[2].split(' ')[0] 
                                + str(res).split('=')[3][0] + str(res.get_parent())[-2])
            else:
                res_list.append(str(res).split('=')[0].split(' ')[-2] + ' ' 
                                + str(res).split('=')[2].split(' ')[0] 
                                + str(res).split('=')[3][0] + ' ' + str(res.get_parent())[-2])
        return(res_list)

    # To change residues from [('500', 'A')] to ['TRP500 A']
    def site_format(self, site):
        index_list = []
        for residue in self.res_list_with_lig:
            for res_site in site:
                if f' {res_site[0]} {res_site[1]}' in residue:
                    index_list.append(self.res_list_with_lig.index(residue))
        site_formatted = []
        for i in index_list:
            site_formatted.append(self.res_list_with_lig[i].split(' ')[0] + self.res_list_with_lig[i].split(' ')[1] + ' ' + self.res_list_with_lig[i].split(' ')[2])
        return(site_formatted)
    
    # find the residues which are connected to the residue of interest through weak bonds, which includes the source residues
    def bonded_residues(self, residues, qs_filter, chain_filter):
        # remove the bonds which involve the source ligand(s) if the source used is not site resdiues
        all_weak_bond = []
        for res in residues:
            all_weak_bond = all_weak_bond + [_ for _ in self.all_bonds_qs if res in _[0]]
            # this is to select the weak bonds and followed by select those bonds with qs of a certain range
        weak_bond_selected = [bond[0] for bond in all_weak_bond if qs_filter[0] <= bond[1] <= qs_filter[1]]
        all_res = []
        for bond in weak_bond_selected:
            all_res = all_res + [bond.split(' ')[0] + ' ' + bond.split(' ')[1]] + [bond.split(' ')[4] + ' ' + bond.split(' ')[5]]
        bonded_res = [res for res in all_res if res not in residues]
        bonded_res_processed = []
        for res in bonded_res:
            if res not in bonded_res_processed and res.split(' ')[1] not in chain_filter:
                bonded_res_processed.append(res)

        # this part is to rank the residues from low to high res num
        res_list_updated = [_.split(' ')[0] + _.split(' ')[1] + ' ' + _.split(' ')[2] for _ in self.res_list_with_lig]
        site_pre_process = []
        for i in bonded_res_processed:
            if i not in site_pre_process:
                site_pre_process = site_pre_process + [(int(self.res_list_with_lig[res_list_updated.index(i)].split(' ')[1]), 
                                                    i.split(' ')[1])]
        site_rank = sorted(site_pre_process, key = itemgetter(0), reverse = False)
        site_rank_process = [(str(_[0]), _[1]) for _ in site_rank]
        final_ranked_res = self.site_format(site_rank_process)
        return(final_ranked_res)
    
    def neighbour_res_summary(self):
        res_1 = self.bond['atom1_res_name'].values + self.bond['atom1_res_num'].values + ' ' + self.bond['atom1_chain'].values
        res_2 = self.bond['atom2_res_name'].values + self.bond['atom2_res_num'].values + ' ' + self.bond['atom2_chain'].values
        self.bond['res_1'] = res_1
        self.bond['res_2'] = res_2
        df_cov = self.bond[self.bond['bond_type'] == 'COVALENT']
        cov_sum = df_cov.res_1.values + ',' + df_cov.res_2.values
        cov_sum_filtered = [_ for _ in cov_sum if _.split(',')[0] != _.split(',')[1]]
        summary = defaultdict(list)
        for connection in cov_sum_filtered:
            summary[connection.split(',')[0]].append(connection.split(',')[1])
            summary[connection.split(',')[1]].append(connection.split(',')[0])
        return(summary)

    def edge_selection(self, edges, weight_cutoff, chain_filter):
        selected_edges = [_ for _ in edges if _[2]['weight'] >= weight_cutoff and _[0].split(' ')[1] not in chain_filter and _[1].split(' ')[1] not in chain_filter]
        return(selected_edges)

    def score_path(self, paths):
        G = nx.Graph()
        G.add_nodes_from([_.split(' ')[0] + _.split(' ')[1] + ' ' + _.split(' ')[2] for _ in self.res_list_no_lig])
        G.add_edges_from(self.all_weak_edges)

        if paths == []:
            max_score = 0
            path_max_score = 'No pathway found.'
        else:
            each_path_score = {}
            for path in paths:
                if len(path) == 2 and path[0] in self.cov_bond_res_summary[path[1]]:
                    each_path_score['_'.join(path)] = 1
                else:
                    qs_path = [G[path[n]][path[n + 1]]['weight'] for n in range(len(path) - 1) if path[n] not in self.cov_bond_res_summary[path[n + 1]]]
                    each_path_score['_'.join(path)]= np.prod(np.array(qs_path) ** (1 / len(qs_path)))
            max_score = each_path_score[sorted(each_path_score, key = each_path_score.get, reverse = True)[0]]
            path_max_score = [sorted(each_path_score, key = each_path_score.get, reverse = True)[0]]
        return(max_score, path_max_score)

    def weak_edges_list(self):
        bonds = []
        for _ in self.all_bonds_qs:
            # ignore any bond containg source residues
            if any(item in _[0] for item in self.source) == False:
                res_1 = _[0].split(' ')[0] + ' ' + _[0].split(' ')[1]
                res_2 = _[0].split(' ')[4] + ' ' + _[0].split(' ')[5]
                if res_1 != res_2:
                    bonds.append(((res_1, res_2), _[1]))
            
        weight_dict = defaultdict(list)
        for i in bonds:
            weight_dict[i[0]].append(i[1])

        # for two residues connected by multiple bonds, only keep the bond with the highest bond_qs
        for j in weight_dict.keys():
            weight_dict[j] = max(weight_dict[j])

        edges = [(bonds[0], bonds[1], {'weight': weight}) for bonds, weight in weight_dict.items() if weight != 0]

        # below is to further selected the highest scoring bond as there are two sequences per pair of residues
        index_list = []
        for _ in edges:
            for i in range(len(edges)):
                if _[0] == edges[i][1] and _[1] == edges[i][0]:
                    if edges.index(_) < i and (edges.index(_), i) not in index_list:
                        index_list.append((edges.index(_), i))
                    elif edges.index(_) > i and (i, edges.index(_)) not in index_list:
                        index_list.append((i, edges.index(_)))
                        
        edges_to_remove = []
        for indices in index_list:
            if edges[indices[0]][2]['weight'] > edges[indices[1]][2]['weight']:
                edges_to_remove.append(edges[indices[1]])
            else:
                edges_to_remove.append(edges[indices[0]])

        for edge in edges_to_remove:
            edges.remove(edge)
        
        return(edges)

    def cov_edges_list(self):
        cov_edges = []
        for res in self.cov_bond_res_summary.keys():
            cov_edge = [(res, _, {'weight': float(1)}) for _ in self.cov_bond_res_summary[res]]
            cov_edges += cov_edge
        return(cov_edges)

    def bond_qs_selection_cutoff(self):
        qs_lower_bound_all = []
        for qs in self.bond_propensity['qs'].values:
            if qs not in qs_lower_bound_all:
                qs_lower_bound_all.append(qs)
        qs_lower_bound_all = sorted(qs_lower_bound_all, reverse = True)
        percentage = np.arange(0, 1.05, 0.05)
        def normal_round(value):
            return(int(value + 0.5))
        qs_lower_bound = [float(1)] + [qs_lower_bound_all[normal_round((len(qs_lower_bound_all) - 1) * _)] 
                                    for _ in percentage if _ != 0 and _ != 1] + [float(0)]
        return(qs_lower_bound)

    # generate a dictionary of all the residues with corresponding residues (bond through weak bonds and the neighbours)
    def simplified_pg(self, qs_lower_bound, chain_filter):
        simplified_graph = nx.Graph()
        simplified_graph.add_nodes_from([_.split(' ')[0] + _.split(' ')[1] + ' ' + _.split(' ')[2] for _ in self.res_list_no_lig])
        simplified_graph.add_edges_from(self.edge_selection(self.all_weak_edges, qs_lower_bound, chain_filter))
        simplified_graph.add_edges_from(self.cov_edges_list())
        removed_edges = []
        for res in self.allo_lig:
            removed_edges = removed_edges + [edge for edge in list(simplified_graph.edges([res]))]
        simplified_graph.remove_edges_from(removed_edges)
        simplified_graph.remove_nodes_from(self.allo_lig)
        return(simplified_graph) 
        
    # find the shortest path from the above dictionary generated network with a specific qs range
    def shortest_paths(self, simplified_graph, start_res, target_res):
        if start_res == target_res:
            selected_paths = ['identical start and target residue']
            no_steps = 0
        else:
            try:
                shortest_paths = list(nx.all_shortest_paths(simplified_graph, start_res, target_res))
            except Exception:
                shortest_paths = []   
                
            selected_paths = []
            for path in shortest_paths:
            # allow only maximum one movement for consecutive residues i.e. move only one step if covalently bonded
                check = []
                for n in range(len(path) - 1):
                    if path[n - 1] in self.cov_bond_res_summary[path[n]] and path[n + 1] in self.cov_bond_res_summary[path[n]] and len(path) != 2:
                        check.append('remove')
                if 'remove' not in check:
                    selected_paths.append(path)

            # the number bonds between the source and target in the path
            if selected_paths == []:
                no_steps = 'none as there is no pathway'
            else:
                no_steps = int(len(selected_paths[0])) - 1
        return(selected_paths, no_steps)      

    # final function, given all the results and the starting residue(s), compute qs_path for targetting residues respectively
    
    def calculate_POP(self, PDB, chain_filter, start_res, target_res):
        qs_lower_bound = self.bond_qs_selection_cutoff()
        
        qs_protein = []
        all_paths = []
        for residue in target_res:
            qs_each_res = []
            path_each_res = []
            if residue in self.source:
                qs_each_res.append(0)
                path_each_res.append('source residue')
            elif start_res == residue:
                qs_each_res.append(1)
                path_each_res.append('identical start and target residue')
            else:
                paths_protein = []
                for i in range(len(qs_lower_bound)):
                    paths_protein = self.shortest_paths(self.simplified_pg(qs_lower_bound[i], chain_filter), start_res, residue)[0]
                    if paths_protein != [] and paths_protein != ['identical start and target residue']:
                        paths_scores = self.score_path(paths_protein)
                        qs_each_res.append(paths_scores[0])
                        path_each_res = path_each_res + paths_scores[1]
                        break
            qs_protein.append(qs_each_res)
            all_paths.append(path_each_res)
        # process the final results, assign residues with no pathways as qs_path = 0
        qs_protein_updated = [[0] if qs == [] else qs for qs in qs_protein]
        paths_protein_updated = ['no pathway found' if path == [] else path[0] for path in all_paths]
        qs_path_scoring = [qs[0] for qs in qs_protein_updated] # this is just to remove the brackets []
        qs_path_scoring = [1 if math.isnan(x) else x for x in qs_path_scoring]
        qs_path_table = {'res': target_res, 'qs_path': qs_path_scoring, 'path': paths_protein_updated}
        qs_path_summary = DataFrame(qs_path_table, columns = ['res', 'qs_path', 'path'])
        print(f'Calculation finished for {start_res} as the starting residue.')
        qs_path_summary.to_csv(f'{PDB}_qs_path_{start_res[:-2]}_{start_res[-1]}.csv')

    # this function to compute characteristic path upon removal of any residue in the protein to quantify mutations
    def calculate_characteristic_path(self, PDB, res_remove, start_res, target_res, chain_filter, results_dir = None):
        qs_lower_bound = self.bond_qs_selection_cutoff()

        os.chdir(results_dir)

        #create a new folder for results
        directory = f'{res_remove}_removed'
        if os.path.exists(directory):
            shutil.rmtree(directory)
            os.mkdir(directory)
        else:
            os.mkdir(directory)
        os.chdir(directory)
        print(f'This is the removed residue {res_remove}.')

        for starting_res in start_res:
            qs_protein = []
            all_paths = []
            no_steps = []
            for residue in target_res:
                qs_each_res = []
                path_each_res = []
                no_steps_each_res = []
                for i in range(len(qs_lower_bound)):
                    ori_graph = self.simplified_pg(qs_lower_bound[i], chain_filter)
                    removed_edges = [edge for edge in list(ori_graph.edges([res_remove])) if edge[0] not in self.cov_bond_res_summary[edge[1]]]
                    ori_graph.remove_edges_from(removed_edges)
                    paths_protein = self.shortest_paths(ori_graph, starting_res, residue)
                    if paths_protein[0] != []:
                        #find the shortest path with the highest qs_path and note down the path
                        paths_scores = self.score_path(paths_protein[0])
                        qs_each_res.append(paths_scores[0])
                        path_each_res = path_each_res + paths_scores[1]
                        no_steps_each_res = no_steps_each_res + [len(paths_scores[1][0].split('_')) - 1]
                        break
                qs_protein.append(qs_each_res)
                all_paths.append(path_each_res[0])
                no_steps.append(no_steps_each_res[0])
            print(f'Calculation finished for {starting_res} as the staring residue.')
            # process the final results, assign residues with no pathways as qs_path = 0
            qs_protein_updated = [[0] if qs == [] else qs for qs in qs_protein]
            qs_path_scoring = [qs[0] for qs in qs_protein_updated] # this is just to remove the brackets []
            qs_path_table = {'res': target_res, 'qs_path': qs_path_scoring, 'path': all_paths, 'path_length': no_steps}
            qs_path_summary = DataFrame(qs_path_table, columns = ['res', 'qs_path', 'path', 'path_length'])
            qs_path_summary.to_csv(f'{PDB}_qs_path_{starting_res[:-2]}_{starting_res[-1]}_with_{res_remove[:-2]}_{res_remove[-1]}_removed.csv')
        os.chdir(results_dir)

    
    def plot_heatmap(self, act_site, allo_site, POP_results_folder = None):
        file_path = POP_results_folder
        results = os.listdir(file_path)
        os.chdir(file_path)
        qs_path_all = []
        for act_res in act_site:
            file = [_ for _ in results if f'{act_res[:-2]}_{act_res[-1]}' in _][0]
            df = pd.read_csv(file)
            qs_path_values = np.array([df[df['res'] == allo_res].qs_path.values[0] for allo_res in allo_site])
            qs_path_all.append(qs_path_values)

        path_score = []
        for qs in qs_path_all:
            path_score.append([1 if math.isnan(x) else x for x in qs])
            
        path_score_round = []
        for score in path_score:
            path_score_round.append([round(qs, 3) for qs in score])

        z = path_score
        x = allo_site
        y = act_site
        z_text = path_score_round

        fig = go.FigureWidget(ff.create_annotated_heatmap(z, x = x, y = y, annotation_text = z_text, 
                                                          colorscale = 'Viridis', showscale = True))

        fig.update_layout(
            title = {'text': f'Active site to allosteric site connection'},
            xaxis = {'title': 'Allosteric site residues',
                    'side': 'bottom'},
            yaxis_title = 'Active site residues',
            width = len(allo_site) * 120, height = len(act_site) * 40,
            font = dict(family = 'Arial', size = 20, color = 'black'))

        fig.data[0].colorbar = dict(title = 'POP', titleside = 'top')

        fig.write_image(f'{POP_results_folder}/POPs.png', scale = 3, width = len(allo_site) * 120, height = len(act_site) * 40)
