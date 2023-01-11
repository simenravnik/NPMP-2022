import networkx as nx
import matplotlib.pyplot as plt
import os
import re
import pandas as pd
import numpy as np
from collections import Counter
import math
import time
import random

def getEdgesGold(timeSeriesFile, goldStandardFile): 
    df = pd.read_csv(timeSeriesFile, sep="\t", decimal=",")  
    df = df.apply(pd.to_numeric)   
    df = df.dropna() 
    df = df.drop(columns=["Time"]) 
	
    columns = df.columns 
    columnNumbers = np.arange(1, len(df.columns)+1) 
    geneNums = {}

    for c, n in zip(columns,columnNumbers): 
        geneNums[c] = n               
    
    edges = []   

    with open(goldStandardFile) as goldStructureFile:   
        lines = goldStructureFile.read().splitlines()  

    #print(timeSeriesFile)
    #print(goldStandardFile) 

    for line in lines:  
        items = re.findall(r'[^\s]+', line)   
        if items[2] == '1':   
            edges.append((geneNums[items[1]], geneNums[items[0]])) #switch positions to target <- regulator                

    return edges 

def getEdges(DNfilePath):
    # print(DNfilePath)

    edges = [] #list of tuples
    with open(DNfilePath) as structureFile:
        lines = structureFile.read().splitlines()

    for line in lines:
        edge = tuple(re.findall(r'[0-9]+', line))
        edges.append(edge)
    return edges

def plotDirectedNetwork(edges, filename, reverse=False):
    G = nx.DiGraph()
    G.add_edges_from(edges)

    PR = nx.pagerank(G)

    d = dict(G.degree)

    if reverse:
        G = G.reverse(copy=True)
    pos = nx.spring_layout(G) # positions for all nodes
    # nx.draw(G, pos, nodelist=list(d.keys()), node_size=[v * 100 for v in d.values()], with_labels=True)
    PR_values = [v * 10000 for v in PR.values()]
    nx.draw(G, pos, nodelist=list(d.keys()), node_size=PR_values, node_color=PR_values, cmap=plt.cm.coolwarm, with_labels=True)

    plt.subplots_adjust(wspace=0.15, hspace=0.1)
    plt.savefig("./img/visualizations/" + filename, facecolor='white', bbox_inches="tight")
    # displaying the title
    plt.show()

def calculate_degrees(G, degree_list):
    degree_list = sorted([i for i in degree_list if i != 0])

    degree_dict = Counter(degree_list)
    k = list(degree_dict.keys())
    Pk = [i / G.number_of_nodes() for i in degree_dict.values()]

    return k, Pk


def calculate_y(G, degrees, k_min = 1):
    n = 0
    sum_k = 0
    for _, k in degrees:
        if k >= k_min:
            sum_k += math.log(k / (k_min - 1/2))
            n += 1
    return 1 + n/sum_k

def plot_distributions(edges, filename):

    G = nx.DiGraph()
    G.add_edges_from(edges)

    # Calculate all degrees
    x, y = calculate_degrees(G, list(dict(G.degree).values()))
    x_in, y_in = calculate_degrees(G, list(dict(G.in_degree).values()))
    x_out, y_out = calculate_degrees(G, list(dict(G.out_degree).values()))

    fig, axs = plt.subplots(2, 2, figsize=(15,8), dpi=300, facecolor="white")

    color = ["#ff71ce", "#01cdfe", "#05ffa1", "#b967ff", "#fffb96"]

    # All degrees
    axs[0][0].loglog(x, y, "o", label='Ecoli', color=color[0], alpha=0.6)
    axs[0][1].loglog(x, y, "o", label='Ecoli', color=color[0], alpha=0.4)

    # Seperated degrees
    axs[0][0].loglog(x_in, y_in, "o", label='Ecoli indegree', color=color[1], alpha=0.6)
    axs[1][0].loglog(x_in, y_in, "o", label='Ecoli indegree', color=color[1], alpha=0.4)
    
    axs[0][0].loglog(x_out, y_out, "o", label='Ecoli outdegree', color=color[2], alpha=0.6)
    axs[1][1].loglog(x_out, y_out, "o", label='Ecoli outdegree', color=color[2], alpha=0.4)

    # Set labels and legend
    # axs = [ax11]
    for i in axs:
        for ax in i:
            ax.set_xlabel("k")
            ax.set_ylabel("Pk")
            ax.legend()

    # Saving and showing fig
    plt.subplots_adjust(wspace=0.15, hspace=0.1)
    plt.savefig("./img/degree-distribuitons/" + filename, facecolor='white', bbox_inches="tight")
    # plt.show()

def distances(G, n = 100):
    D = []
    for i in G.nodes() if len(G) <= n else random.choices(list(G.nodes()), k=n):
        D.extend([d for d in nx.shortest_path_length(G, source = i).values() if d > 0])
    return D

def info(edges, reverse=False):

    G = nx.DiGraph()
    G.add_edges_from(edges)

    if reverse:
        G = G.reverse(copy=True)

    n = G.number_of_nodes()
    m = G.number_of_edges()

    D = distances(G)
        
    info_dict = dict()

    info_dict["nodes"] = n
    info_dict["isolates"] = nx.number_of_isolates(G)
    info_dict["edges"] = m
    info_dict["self_loops"] = nx.number_of_selfloops(G)
    info_dict["degree"] = 2 * m / n if n != 0 else 0
    info_dict["max_degree"] =  max([k for _, k in G.degree()]) if len(G.degree()) > 0 else 0
    info_dict["density"] = 2 * m / n / (n - 1) if n > 2 else 0
    info_dict["distance"] = sum(D) / len(D) if len(D) > 0 else 0
    info_dict["max_distance"] = max(D) if len(D) > 0 else 0
    if isinstance(G, nx.MultiGraph):
        G = nx.Graph(G)
    info_dict["clustering"] = nx.average_clustering(G) if n > 0 else 0
    
    return info_dict


def plot_networks(organism, networkSizes, methods):
    for networkSize in networkSizes:
        for method in methods:
            data_path = os.path.join("data", "original", organism, str(networkSize))
            results_path = os.path.join("data", "results", organism, str(networkSize), method)

            networkNum = 1
            crossIteration = 1

            DNfilePath = os.path.join(results_path, organism + "-" + str(networkNum) + "_" + str(crossIteration) + "_structure.tsv")
            timeSeriesFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_dream4_timeseries.tsv")
            goldStandardFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_goldstandard.tsv")

            edgesGold = getEdgesGold(timeSeriesFilePath, goldStandardFilePath)
            plotDirectedNetwork(edgesGold, filename=str(str(networkSize) + "-" + organism + "-" + str(networkNum) +"_goldstandard.png"))

            edges = getEdges(DNfilePath)
            plotDirectedNetwork(edges, reverse=True)

            # plot_distributions(edges=edgesGold, filename=str(str(networkSize) + "-" + organism + "-" + str(networkNum) +"_goldstandard.png"))

if __name__ == "__main__":

    organism = "Ecoli"
    networkSizes = [16, 32, 64]
    methods = ["GABNI", "MIBNI", "BestFit", "REVEAL", "ATEN"]

    # plot_networks(organism, networkSizes, methods)

    for networkSize in networkSizes:

        # Gold standard
        data_path = os.path.join("data", "original", organism, str(networkSize))
        gold_standard_infos = []
        for networkNum in range(1, 11):

            data_path = os.path.join("data", "original", organism, str(networkSize))
            timeSeriesFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_dream4_timeseries.tsv")
            goldStandardFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_goldstandard.tsv")

            edgesGold = getEdgesGold(timeSeriesFilePath, goldStandardFilePath)
            gold_standard_info = info(edgesGold)
            gold_standard_infos.append(gold_standard_info)

        df_gold = pd.DataFrame.from_dict(gold_standard_infos)
        df_gold.to_csv(os.path.join(data_path, "gold_standard_network_info_" + str(networkNum) + ".csv"), index=False)

        # Results
        for method in methods: 
            results_path = os.path.join("data", "results", organism, str(networkSize), method)
            average_infos = []
            for networkNum in range(1, 11):
                infos = []
                for crossIteration in range(0, 10):

                    DNfilePath = os.path.join(results_path, organism + "-" + str(networkNum) + "_" + str(crossIteration) + "_structure.tsv")
                    # timeSeriesFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_dream4_timeseries.tsv")
                    # goldStandardFilePath = os.path.join(data_path, organism + "-" + str(networkNum) +"_goldstandard.tsv")

                    # edgesGold = getEdgesGold(timeSeriesFilePath, goldStandardFilePath)
                    
                    if os.path.exists(DNfilePath):
                        edges = getEdges(DNfilePath)
                        infos.append(info(edges, reverse=True))

                df = pd.DataFrame.from_dict(infos)
                df.to_csv(os.path.join(results_path, "results_network_info_" + str(networkNum) + ".csv"), index=False)
                df2 = df.mean(axis=0)
                average_infos.append(df2)
            
            df3 = pd.DataFrame.from_dict(average_infos)
            df3.to_csv(os.path.join(results_path, "results_network_info_all" + ".csv"), index=False)
