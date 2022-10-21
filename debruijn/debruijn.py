#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
from importlib.resources import path
import os
from pickle import FALSE
from re import T
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Adriana Lecourieux"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Adriana Lecourieux"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Adriana Lecourieux"
__email__ = "adriana.lecourieux@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()

#dans fastq plein d'info, on veut lire la seq
#la premier fonction fait yield de tous les reads
#lit une ligne puis on passe les 3 suivantes
#with open est un generateur de fichier, si on fait next donne ligne suivante
#for i du fichier ouvert (générateur)

def read_fastq(fastq_file):
    """Read fastq file and conserve sequence in a generator.

    Parameters
    ----------
    fastq_file : fastq file
        Name of the fastq file

    Returns
    -------
    yield
        A sequences generator 
    """
    with open (fastq_file, "r") as file: # open the file
        for ligne in file: # iteration on each line
            yield next(file).strip() # save generator
            next(file) # skip lines
            next(file)
            

def cut_kmer(read, kmer_size):
    """Create a k-mer generator.

    Parameters
    ----------
    read :
        sequence from a generator 
    kmer_size : int
        length of k-mer

    Returns
    -------
    yield
        A k-mers generator 
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i+kmer_size]
        
        
    

def build_kmer_dict(fastq_file, kmer_size):
    """Create a kmer dictionnary.

    Parameters
    ----------
    fastq_file  :
        fatq file
    kmer_size : int
        length of k-mer

    Returns
    -------
    kmer_dict :
        dictionnary of all the k-mers with their occurence
    """
    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] +=1
    return(kmer_dict)
            

def build_graph(kmer_dict):
    """Create Digraph from k-mers prefixes and suffixes described.

    Parameters
    ----------
    kmer_dict  :
        dictionnary of all the k-mers with their occurence

    Returns
    -------
    digraph
        tree of k-mers prefixes and suffixes described
    """
    digraph = nx.DiGraph()
    for kmer, occurence in kmer_dict.items():
        digraph.add_edge(kmer[:-1], kmer[1:], weight = occurence)
    return(digraph)
        
        

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Create a kmer dictionnary.

    Parameters
    ----------
    graph  : digraph
        tree of k-mers prefixes and suffixes described
    path_list : list
        list of paths
    delete_entry_node : boolean
        true if enter node of each path will be remove
    delete_sink_node : boolean
        true if output node of each path will be remove
    Returns
    -------
    graph : digraph
        tree of k-mers prefixes and suffixes described
    """
    for node in path_list:
        if delete_entry_node == True and delete_sink_node == True:
            graph.remove_nodes_from(node)
        elif delete_entry_node == True:
            graph.remove_nodes_from(node[:-1])
        elif delete_sink_node == True:
            graph.remove_nodes_from(node[1:])
        elif delete_entry_node == False and delete_sink_node == False:
            graph.remove_nodes_from(node[1:-1])
    return(graph)

def std(data):
    """Calcul standard deviation.

    Parameters
    ----------
    data  : list
        list 

    Returns
    -------
    data : float
        Result of standard deviation of data
    """
    return (statistics.stdev(data))

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    weight_std = std(weight_avg_list)
    length_std = std(path_length)
    
    if weight_std != 0:
        best_path = weight_avg_list.index(max(weight_avg_list))
    elif weight_std == 0:
        if length_std > 0:
            best_path = path_length.index(max(path_length))
        elif length_std == 0:
            best_path = randint(0, len(path_list)-1)
    
    for index in range(len(path_list)):
        if path_list[index] != path_list[best_path]:
            graph = remove_paths(graph, [path_list[index]], delete_entry_node, delete_sink_node)
    return(graph)

def path_average_weight(graph, path):
    """Calcul average weight of paths.

    Parameters
    ----------
    graph  : digraph from networkx  

    Returns
    -------
    float
        Result of average weight of paths
    """
    return(statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]))
    
def solve_bubble(graph, ancestor_node, descendant_node):
    """Solve bubbles.

    Parameters
    ----------
    graph  : digraph from networkx  
    ancestor_node :
        the ancestor node
    descendant_node :
        the descendant node
    Returns
    -------
    graph
        digraph without bubble
    """
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = []
    weight_avg_list = []
    
    for path in path_list:
        weight_avg_list.append(path_average_weight(graph, path))
        path_length.append(len(path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    
    return(graph)

def simplify_bubbles(graph):
    """Solve bubbles.

    Parameters
    ----------
    graph  : digraph from networkx  
    
    Returns
    -------
    graph
        digraph without bubble
    """
    bubble = False
    graph_node = graph.nodes()
    for node in graph_node:
        predecessor_list = list(graph.predecessors(node))
        if len(predecessor_list) > 1:
            for i in range(len(predecessor_list)-1):
                node_i = predecessor_list[i]
                for j in range(i+1, len(predecessor_list)):
                    node_j = predecessor_list[j]
                    ancestor_node = nx.lowest_common_ancestor(graph,node_i, node_j)
                    if ancestor_node != None:
                        bubble = True
                        break
            if bubble == True:
                break         
    if bubble == True:                
      graphe = simplify_bubbles(solve_bubble(graph, ancestor_node, node))      
    return(graph)

def solve_entry_tips(graph, starting_nodes):    
    # arriver à générer les combinaisons de points d'entrees
    
    predecessor = []
    path_list = []
    path_length = []
    weight_avg_list = []
    graph_node = graph.nodes()
    
    for node in starting_nodes:
        for desc in nx.descendants(graph, node):
            if len(graph.pred[desc]) > 2 and desc not in predecessor:
                predecessor.append(desc)
    for i in starting_nodes:            
        for pred in predecessor:
            for path in nx.all_simple_paths(graph, i, pred):
                path_list.append(path)
                path_length.append(path)
                if len(path) > 2:
                    weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=True, delete_sink_node=False)         
    return(graph)

def solve_out_tips(graph, ending_nodes):
    descendant = []
    path_list = []
    path_length = []
    weight_avg_list = []
    graph_node = graph.nodes()
    
    for node in ending_nodes:
        for desc in nx.ancestors(graph, node):
            if len(graph.pred[desc]) > 2 and desc not in descendant:
                descendant.append(desc)
    for i in ending_nodes:            
        for pred in descendant:
            for path in nx.all_simple_paths(graph, i, pred):
                path_list.append(path)
                path_length.append(path)
                if len(path) > 2:
                    weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=True)         
    return(graph)

def get_starting_nodes(graph):
    """Create a list of all the enters nodes.

    Parameters
    ----------
    graph  :
        Digraph from networkx 

    Returns
    -------
    starting_node : list
        list of enter nodes
    """
    starting_node = []
    for node in graph.nodes():
        if not list(graph.predecessors(node)):
            starting_node.append(node)
    return(starting_node)
                       

def get_sink_nodes(graph):
    """Create a list of all the output nodes.

    Parameters
    ----------
    graph  :
        Digraph from networkx 

    Returns
    -------
    ending_node : list
        list of output nodes
    """
    ending_node = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            ending_node.append(node)
    return(ending_node)

def get_contigs(graph, starting_nodes, ending_nodes):
    """Create a tuple of contigs with their length.

    Parameters
    ----------
    graph  :
        Digraph from networkx 
    starting_node
        list of enter nodes
    ending_node :
        list of output nodes
    
    Returns
    -------
    contig : tuple
        tuple of contigs with their length
    """
    contig = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node, end_node) == True:
                for path in nx.all_simple_paths(graph, start_node,end_node): #pour le premier chemin
                    seq = path[0] #1er node
                    for node in path[1:]: #pour chaque nouveau node on prend que la pos 2 pour "merge"
                        seq += node[-1]
                    contig.append(tuple((seq, len(seq))))
    return(contig)

def save_contigs(contigs_list, output_file):
    """Create a file with the contigs in fasta format.
    Parameters
    ----------
    contigs_list  : list
        Digraph from networkx 
    output_file : string
        name of output file
    """
    with open(output_file, "w") as file:
        for i, (contig, length) in enumerate(contigs_list):
            file.write(f'>contig_{i} len={length}\n')
            file.write(f'{fill(contigs_list[i][0], width=80)}\n')
            


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, "contigs_list.txt")
    
    path_average_weight(graph, path)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
