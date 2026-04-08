import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import cobra
import numpy as np
import seaborn as sb  # for correlation matrix
from scipy.cluster.hierarchy import dendrogram, linkage  # for hierarchical clustering
import matplotlib.pyplot as plt  # for hierarchical clustering
import sklearn.preprocessing as prep  # for PCA
from sklearn.decomposition import PCA  # for PCA
from sklearn.cluster import KMeans, AgglomerativeClustering  # for PCA
from collections import defaultdict  # for summarize_pathways
import networkx as nx  # for flux_coupling_network
import matplotlib.colors as mcolors  # for color_network_by_communities
from networkx import (
    edge_betweenness_centrality as betweenness,
)  # for color_networks_by_communities including edge weights
from plotly.subplots import make_subplots  # for community box and bar plots
from math import ceil  # for community box and bar plots
import hvplot.pandas  # for more hover data in PCA
import contextlib  # for community_pw_summary to suppress printing
import io  # for community_pw_summary to suppress printing

# from IPython.core.display import display, HTML # for printing example


## code by Nina Engelkamp ##


def help():
    """
    Function to provide an overview of all the available functions and how they can be called

    Return: str
    ----------
    information text
    """
    text = "Available functions:\n- build_dataframe(model, solutions, conditions)\n- scatter_plot_fluxes(fba, conditions) \
        \n- correlation_matrix(fba, conditions, [opt: color])\n- filter_dataframe(fba, conditions)\n- bar_plot_flux_variation(fba, [opt: start, end]) \
        \n- hierarchical_clustering_std_dev(fba, [opt: start, end, d_max])\n- pca(fba, base_condition, other_conditions, [opt: clusters, dimensions]) \
        \n- normalize_dataframe_rows(fba, conditions [opt: method])\n- normalize_dataframe_cols(fba, conditions, obj_values) \
        \n- pca_color_pathways(fba, base_condition, other_conditions, pathway, pathway_dict, [opt: dimensions=2]) \
        \n- summarize_pathways(id_list, id_dict)\n- flux_coupling_matrix(fba, conditions, remove_unchanging = 0.01, abs_values = True)\
        \n- flux_coupling_network(flux_coupling_matrix = pd.DataFrame(), pc = 0.8)\n- color_network_by_pathway(nw, pathway, pathway_dict)\
        \n- color_network_by_communities(nw, iterations)\n- visualize_network(nw, [opt: layout])\n- get_community_list(nw, iterations)\
        \n- community_bar_plots(fba, conditions, com_list)\n- community_box_plots(fba, conditions, com_list)\
        \n- community_pw_summary(com_list, id_dict, pw_dict, [opt: lower_bound, upper_bound, cut_off])\
        \n- community_sankey_diagram(com_list, pw_list, fba, id_dict, [opt: diagram_height, pw_hierarchy])\
        \n- color_network(nw, [opt: color])\
        \n\n For an example of how these functions can be used call example()"

    print(text)


def example():
    """
    Function to show an example of how this program can be used

    Return: void
    ----------
    displays information text
    """
    display(
        HTML(
            "<u>Example on how this program can be used:</u><br>Here are some examples how to use the different available functions and some suggestions on how to use the parameters.<br><br>\
    <u>1.) Required data</u><br><div style='margin-left: 20px; margin-top: 0; margin-bottom: 0; padding: 0'>\
    <span style='color:blue;'>model</span> &rarr; can be created with e.g. load_json_model or read_sbml_model with cobrapy\
    <br><span style='color:blue;'>solutions = [sol1, sol2, sol3, ...]</span> &rarr; list of solutions for every condition, can be created with model.optimize() in cobrapy\
    <br><span style='color:blue;'>conditions = [\"name1\", \"name2\", \"name3\", ...]</span> &rarr; list of names for the conditions (in the same order as the solutions), \
    these names will be the names used in the created plots</div><u>\
    <br>\
    2.) Create a dataframe<br></u><div style='margin-left: 20px; margin-top: 0; margin-bottom: 0; padding: 0'><span style='color:blue;'>fba = build_dataframe(model, solutions, conditions)</span></div><u>\
    <br>\
    3.) Preprocess the data with methods of your choice e.g.:</u><br><div style='margin-left: 20px; margin-top: 0; margin-bottom: 0; padding: 0'>\
    <ul style='margin: 0; padding-left: 20px; line-height: 1.0;'><li><span style='color:blue;'>filtered_fba = filter_dataframe(fba, conditions)</span> &rarr; removes all reactions where all conditions have zero flux or very close to zero flux</li>\
    <br><li><span style='color:blue;'>obj_values = [obj_value1, obj_value2, obj_value3, ...]</span> &rarr; list of objective values for each condition, obj_value1 can be created with sol1.objective_value\
    <br><span style='color:blue;'>normalized_fba = normalize_dataframe_cols(filtered_fba, conditions, obj_values)</span>  &rarr; performs normalization of the conditions to make them more comparable to eachother</li>\
    <br><li><span style='color:blue;'>reaction_wise_normalized_fba = normalize_dataframe_rows(normalized_fba, condition, method = \"div_by_max\")</span>\
         &rarr; performs normalization of the reactions, so that they all have the same flux range from -1 to 1. This can be helpful for some plots, e.g. the correlation matrices</li></ul></div><br>\
    <br>\
    <u>4.) Create statistics to analyze and visualize the metabolic changes in the conditions</u><br><div style='margin-left: 20px; margin-top: 0; margin-bottom: 0; padding: 0'>\
    There are many different statistics available to visualize different aspects of the comparison. Here is an overview. <br>Note: Some of the statistics use a pathway_dict or an id_dict. The pathway_dict is \
    a python dictionary mapping biological pathway names to reaction IDs. The id_dict is a pathway dictionary mapping reactions to the pathways they are in. These dicts can for example be created using\
    the pathways from the KEGG database.<br>\
    <br>\
    Statistics to compare the entire conditions against each other:\
    <ul style='margin: 0; padding-left: 20px; line-height: 1.0;'><br><li><span style='color:blue;'>correlation_matrix(reaction_wise_normalized_fba, conditions)</span></li>\
    <br><li><span style='color:blue;'>scatter_plot_fluxes(normalized_fba, ['name1', 'name2']]) </span></li></ul>\
    <br>\
    Statistics to identify most/least variable reactions:<ul style='margin: 0; padding-left: 20px; line-height: 1.0;'>\
    <br><li><span style='color:blue;'>bar_plot_flux_variation(normalized_fba, start=0, end=20)</span> &rarr; shows the 20 most variable reactios</li></ul>\
    <br>\
    Statistics to group reactions by variability:<ul style='margin: 0; padding-left: 20px; line-height: 1.0;'>\
    <br><li><span style='color:blue;'>hierarchical_clustering_std_dev(normalized_fba, start=0, end=50, d_max=0.3]) </span></li>\
    <br><li><span style='color:blue;'>pca(normalized_fba, conditions,  method=\"Differences\", base_condition=\"Mean\") </span></li>\
    <br><li><span style='color:blue;'>pca_color_by_pathway(normalized_fba, conditions, \"Glycolysis / Gluconeogenesis\", pathway_dict, method=\"Differences\", base_condition=\"Mean\") </span></li>\
    <br><li><span style='color:blue;'>pca_color_by_list(normalized_fba, conditions, ['reaction1', 'reaction2', ...], method=\"Differences\", base_condition=\"Mean\") </span></li></ul>\
    <br>\
    Statistics to identify most/least variable pathways:<ul style='margin: 0; padding-left: 20px; line-height: 1.0;'>\
    <br><li><span style='color:blue;'>pathway_variability(normalized_fba, pathway_dict, num=10)</span> &rarr; prints 10 most variable pathways</li>\
    <br><li><span style='color:blue;'>pathway_variability(normalized_fba, pathway_dict, ascending=False, num=10)</span> &rarr; prints 10 least variable pathways</li>\
    <br><li><span style='color:blue;'>summarize_pathways(['reaction1', 'reaction2', ...], id_dict) </span></li></ul>\
    <br>\
    Statistics to identify coupled reactions/pathways:<ul style='margin: 0; padding-left: 20px; line-height: 1.0;'>\
    <br><li><span style='color:blue;'>pca(normalized_fba, conditions,  method=\"Correaltion\") </span></li>\
    <br><li><span style='color:blue;'>pca_color_by_pathway(normalized_fba, conditions, \"Carbon metabolism\", pathway_dict, method=\"Correlation\") </span></li>\
    <br><li><span style='color:blue;'>pca_color_by_list(normalized_fba, conditions,['reaction1', 'reaction2', ...], method=\"Correlation\") </span></li>\
    <br><li><span style='color:blue;'>corr_matrix = flux_coupling_matrix(reaction_wise_normalized_fba, conditions) </span></li>\
    <br><li><span style='color:blue;'>nw = flux_coupling_network(flux_coupling_matrix = corr_matrix) </span></li>\
    <br><li><span style='color:blue;'>color_network(nw, color=\"blue\") </span></li>\
    <br><li><span style='color:blue;'>color_network_by_pathway(nw, 'Citrate cycle (TCA cycle)', pathway_dict) </span></li>\
    <br><li><span style='color:blue;'>color_network_by_communities(nw, 5, include_weights = True) </span></li>\
    <br><li><span style='color:blue;'>visualize_network(nw) </span></li>\
    <br><li><span style='color:blue;'>com_list = get_community_list(nw) </span></li>\
    <br><li><span style='color:blue;'>color_map = get_color_map(com_list) </span></li>\
    <br><li><span style='color:blue;'>community_box_plots(reaction_wise_normalized_fba, conditions, com_list)</span></li>\
    <br><li><span style='color:blue;'>community_sankey_diagram(com_list, ['Carbon metabolism', 'Nucleotide metabolsm', ...], reaction_wise_normalized_fba, id_dict) </span></li>\
    <br><li><span style='color:blue;'>community_pw_summary(com_list, id_dict, pathway_dict) </span></li>\
    <br><li><span style='color:blue;'>summarize_pathways(color_map[3], id_dict) </span></li></ul></div>"
        )
    )


def build_dataframe(model, solutions, conditions):
    """
    Function to build data frame from multiple FBA solutions for different conditions

    Parameters
    ----------
    model : cobra.core.model.Model
        common cobra metabolic model (needs to be the same for all conditions [-> have the same reactions])

    solutions: List[cobra.core.solution.Solution]
        list of FBA solutions for all conditions of interest

    conditions: List[str]
        list of condition names in the order of the matching solutions

    Return: pandas.DataFrame
    ----------
    Data frame which lists for each reaction the index, reaction id, reaction formula, fluxes for every condition, the mean of the fluxes and the std dev of the fluxes

    """
    ## check input
    if not (len(solutions) == len(conditions)):
        raise ValueError(
            f"The number of provided conditions ({len(conditions)}) does not match the number of solutions ({len(solutions)})"
        )

    ## create empty dataframe
    empty_df = pd.DataFrame()
    col_names = ["ID", "Reaction", "Name"]
    col_names = col_names + conditions
    fba = pd.DataFrame(columns=col_names)

    ## fill dataframe
    fba["ID"] = [r.id for r in model.reactions]
    fba["Reaction"] = model.reactions
    fba["Name"] = [r.name for r in model.reactions]
    for i, cond in enumerate(
        conditions
    ):  # there is probaby a more effective way to do this
        if not (len(solutions[i].fluxes) == len(fba["Reaction"])):
            raise ValueError(
                f"The condition {cond} has {len(solutions[i].fluxes)} reactions. This does not match the provided model, \
                             which has {len(fba['Reaction'])} reactions."
            )
        fluxes = []
        for r, (reaction, flux) in zip(fba["ID"], solutions[i].fluxes.items()):
            if not (r == reaction):
                raise ValueError(
                    f"The reation {reaction} in condition {cond} does not match the reaction {r} in list from the model!"
                )
            else:
                fluxes.append(flux)
        fba[cond] = fluxes
    fba["Std_dev"] = fba[conditions].std(
        axis=1
    )  # calculate the std_dev for each reaction
    fba["Mean"] = fba[conditions].mean(axis=1)  # calculate the mean for each reaction
    return fba


def filter_dataframe(fba, conditions, rounding=False):
    """
    Function to filter the data

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be filtered

    rounding: bool
        default: False
        rounds all values that are very close to zero (np.isclose(0) == True) to zero

    Return: pandas.DataFrame
    ----------
    Filtered DataFrame
    """
    if rounding:
        fba[conditions] = fba[conditions].where(~np.isclose(fba[conditions], 0), 0)
    return fba.loc[
        (fba[conditions] != 0).any(axis=1)
    ]  # remove reactions where all fluxes are zero


def normalize_dataframe_cols(fba, conditions, obj_values):
    """
    Function to normalize the data between conditions

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be normalized

    obj_values: List[float]
        list of the objective values of each condition in the order of the conditions

    Return: pandas.DataFrame
    ----------
    Filtered DataFrame
    """
    ## check input
    if not (len(conditions) == len(obj_values)):
        raise ValueError(
            "The number of conditions does not match the number of objective values"
        )

    normalized_fba = fba.copy()
    for c, v in zip(conditions, obj_values):
        normalized_fba[c] = np.divide(normalized_fba[c], v)
    normalized_fba["Std_dev"] = normalized_fba[conditions].std(
        axis=1
    )  # calculate new mean and std_dev
    normalized_fba["Mean"] = normalized_fba[conditions].mean(axis=1)
    return normalized_fba


def normalize_dataframe_rows(fba, conditions, method="div_by_max"):
    """
    Function to normalize the data between reactions
    (important note: this function does not normalize the stde or mean of the reactions)

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be normalized

    method: "standard" / "minmax" / "div_by_sum" / "div_by_max"
        (optional) normalization method, default is "div_by_max"
        - "standard" : substract mean and divide by standard deviation (biased) -> values cetered around 0 with std_dev ~ 1
        - "minmax": substract min and divide by (max - min) -> values in range [0, 1]
        - "div_by_sum": divide by sum and multiply with 100 -> values scaled in range [0], 100]
        - "div_by_max": divide by absolute value of max -> values in [-1, 1], keeps signs

    Return: pandas.DataFrame
    ----------
    FLux DataFrame with normalized reactions
    """

    ## initialize new dataframe and copy information columns
    fba_normalized = pd.DataFrame()
    fba_normalized["ID"] = fba["ID"]
    fba_normalized["Reaction"] = fba["Reaction"]
    fba_normalized["Name"] = fba["Name"]

    ## normalize based on selected method
    if method == "div_by_sum":
        min = fba[conditions].min(axis=None, numeric_only=True)
        shifted = fba[conditions]
        if min < 0:
            shifted -= min
        normalized = (
            shifted[conditions]
            .div(shifted[conditions].sum(axis=1), axis=0)
            .mul(100, axis=0)
        )
        for c in conditions:
            fba_normalized[c] = normalized[c]
        fba_normalized["Std_dev"] = fba["Std_dev"]
        fba_normalized["Mean"] = fba["Mean"]
        return fba_normalized

    if method == "div_by_max":
        normalized = fba[conditions].div(fba[conditions].abs().max(axis=1), axis=0)
        for c in conditions:
            fba_normalized[c] = normalized[c]
        fba_normalized["Std_dev"] = fba["Std_dev"]
        fba_normalized["Mean"] = fba["Mean"]
        return fba_normalized

    if method == "div_by_max_pos":
        min = fba[conditions].min(axis=None, numeric_only=True)
        shifted = fba[conditions]
        if min < 0:
            shifted -= min
        normalized = shifted[conditions].div(shifted[conditions].max(axis=1), axis=0)
        for c in conditions:
            fba_normalized[c] = normalized[c]
        fba_normalized["Std_dev"] = fba["Std_dev"]
        fba_normalized["Mean"] = fba["Mean"]
        return fba_normalized

    scaler = prep.StandardScaler()
    if method == "minmax":
        scaler = prep.MinMaxScaler()
    elif method == "standard":
        pass
    else:
        raise ValueError(
            'Invalid mathod. Valid methods are "minmax", "standard", "div_by_sum" and "div_by_max"'
        )

    normalized = scaler.fit_transform(fba[conditions].values.T)
    for i, c in enumerate(conditions):
        fba_normalized[c] = normalized[i]
    fba_normalized["Std_dev"] = fba["Std_dev"]
    fba_normalized["Mean"] = fba["Mean"]
    return fba_normalized


def scatter_plot_fluxes(fba, conditions, font_size=8):
    """
    Function to create pairwise scatter plot of the reaction fluxes

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be compared

    font_size: int
        optional, default is 8
        adjust font size in scatterplots compating many conditions to avoid overlapping labels

    Return: plotly.graph_objs._figure.Figure
    ----------
    Scatter plot that can be viewed with .show()
    """

    if len(conditions) == (0 or 1):
        raise ValueError(f"The condition list needs to contain at least 2 arguments")
    if len(conditions) == 2:
        diff = [abs(x - y) for (x, y) in zip(fba[conditions[0]], fba[conditions[1]])]
        fig = px.scatter(
            fba,
            x=conditions[0],
            y=conditions[1],
            hover_name="ID",
            hover_data="Name",
            color=diff,
            title="Scatter plot of reaction fluxes, colored by absolute flux difference",
        )
        fig.update_layout(coloraxis_colorbar_title_text="diff")
    else:
        fig = px.scatter_matrix(
            fba,
            dimensions=conditions,
            title="Pairwise scatter plot of reaction fluxes",
            hover_name="ID",
        )
        fig.update_traces(showupperhalf=False)
        fig.update_traces(diagonal_visible=False)
        fig.update_layout(yaxis={"automargin": True})
        fig.update_layout(font=dict(size=font_size))
    return fig


def correlation_matrix(fba, conditions, color="Blues"):
    """
    Function to create a correlation matrix between the fluxes of coditions

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be compared

    color: str
        [optional] color scheme for the plot, valid colors are for example "Blues", "Greens", "Oranges", "Reds"

    Return: matplotlib.axes._axes.Axes
    ----------
    Correlation matrix
    """
    if len(conditions) == 0 or len(conditions) == 1:
        raise ValueError("There must be at least two conditions in the condition list.")
    else:
        r = fba[conditions].corr()
        return sb.heatmap(r, cmap=color, annot=True)


def bar_plot_flux_variation(
    fba, conditions, start=0, end=20, ascending=False, color="#636EFA"
):
    """
    Function to create a bar plot of the reactions sorted by descending std_dev

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of condition names that will be compared

    start: int
        (optional) can be used to display only a window of the bar plot, starting at position start

    end: int
        (optional) can be used to display only a window of the bar plot, ending at position end

    ascending: bool
        (optional) if set to True then the Std dev will be considered ascending (from lowest to highest) instead of descending

    color: string
        (optional, default is blue) color of the bars in the bar plot (e.g 'blue', 'red', 'green', 'orange')

    Return: plotly.graph_objs._figure.Figure
    ----------
    Barplot that can be viewed with .show()
    """
    ## check for invalid start and end position
    if start < 0 or start > len(fba) - 1:
        raise ValueError("The start position is invalid.")
    if end < start or end > len(fba):
        raise ValueError("The end position is invalid.")
    if end == 0:
        end = len(fba)

    ## sort reactions by descending std_dev
    # sorted_flux = fba[["ID", "Name", "Std_dev"]].sort_values(by=['Std_dev'], ascending=ascending)
    sorted_flux = fba.sort_values(by=["Std_dev"], ascending=ascending)

    ## create and return bar plot
    fig = px.bar(
        sorted_flux[start:end],
        x="ID",
        y="Std_dev",
        hover_data=["Name"] + conditions,
        hover_name="ID",
        color_discrete_sequence=[color],
    )
    return fig


def box_plot_flux_variation(fba, conditions, start=0, end=0, ascending=False):
    fba["ratio"] = abs(fba["Std_dev"] / fba["Mean"])
    sorted_flux = fba.sort_values(by=["ratio"], ascending=False)  ## div by zero???

    # Prepare data for the box plot
    box_data = []
    for idx, row in sorted_flux[start:end].iterrows():
        box_data.append(go.Box(y=row[conditions], name=row["ID"]))

    # Create the figure
    fig = go.Figure(box_data)

    # Update layout for better visualization
    fig.update_layout(
        xaxis_title="Reactions",
        yaxis_title="Flux",
        boxmode="group",  # Optional, for grouping if needed
    )

    # Show the plot
    return fig


""" helper function"""


def _fancy_dendrogram(
    *args, **kwargs
):  # adapted from https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    max_d = kwargs.pop("max_d", None)
    if max_d and "color_threshold" not in kwargs:
        kwargs["color_threshold"] = max_d
    annotate_above = kwargs.pop("annotate_above", 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get("no_plot", False):
        for i, d, c in zip(ddata["icoord"], ddata["dcoord"], ddata["color_list"]):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, "o", c=c)
                plt.annotate(
                    "%.3g" % y,
                    (x, y),
                    xytext=(0, -5),
                    textcoords="offset points",
                    va="top",
                    ha="center",
                )
        if max_d:
            plt.axhline(y=max_d, c="k")
    return ddata


def hierarchical_clustering_std_dev(fba, start=0, end=0, d_max=0):
    """
    Function to create a dendogram according to hierarchcal clustering of the std_dev

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    start: int
        (optional) can be used to only select a window of reactions starting st position start (reactions are sorted by descending std_dev)

    end: int
        (optional) can be used to only select a window of reactions ending st position end (reactions are sorted by descending std_dev)

    d_max: int
        (optional) can be used to create a more annotated version of the dendogram with a distance cut-off d_max for the clustering

    Return: module
    ----------
    Matplotlip dendogram
    """
    ## check input
    if start < 0 or start > len(fba) - 1:
        raise ValueError("The start position is invalid.")
    if end < start or end > len(fba):
        raise ValueError("The end position is invalid.")
    if end == 0:
        end = len(fba)
    if d_max < 0:
        raise ValueError("d_max is not allowed to be negative")

    ## sort reactions by descending std_dev (in case start and end are set to custom values)
    sorted_flux = fba[["ID", "Std_dev"]].sort_values(by=["Std_dev"], ascending=False)
    sorted_flux = sorted_flux[start:end]
    sorted_flux.index = sorted_flux["ID"]

    ## hierarchical clustering
    linkage_data = linkage(sorted_flux[["Std_dev"]], method="ward", metric="euclidean")

    ## create and return dendogram
    plt.figure(figsize=(25, 7))
    plt.xlabel("Reactions")
    plt.ylabel("Distance btw std_dev")
    plt.title("Hierarchical Clustering of the Std_dev")
    if d_max == 0:
        dendrogram(linkage_data, labels=sorted_flux.index, leaf_rotation=90.0)
    else:
        _fancy_dendrogram(
            linkage_data,
            labels=sorted_flux.index,
            leaf_rotation=90.0,
            annotate_above=10,
            max_d=d_max,
        )
    return plt


def pca(
    fba,
    conditions,
    dimensions=2,
    clusters=1,
    method="Differences",
    base_condition="",
    remove_unchanging=0.01,
    abs_values=True,
    display="with_size",
):
    """
    Function to perform PCA colored by K-means clustering

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of all the conditions that will be considered in the PCA

    dimensions: int
        (optional) number of principal components, can be 2 or 3, default is 2

    clusters: int
        (optional) number of clusters, default is 1

    method: str
        (optional) The method to use the data, can be 'Correlation' or 'Differences', default is correlation
        - 'Correlation': for each reaction use a vector of pairwise correlations to all other reactions as data
        - 'Differences': for each reaction use a vector of differences in flux to the base_condition as data

    base_condition: str
        (optional, only used if method = correlation) name of the base condition (the difference from the other conditions will be calculated from this condition)
        default is the first condition in the list of conditions

    remove_unchanging: float
        (optional, only used if method="Correlation", default = 0.01)
        Reactions with a smaller Std_dev than this threshold will be removed from the plot

    abs_values: bool
        (optional, only used if method="Correlation", default = True)
        will consider absolute values when calculating the correaltions and remove reactions where the sign swaps across conditions

    display: str
        default: "regular", options: "regular", "with_size", "hvplot"
        the type of plot that will be returned
        "regular" returns normal plotly plot (unfortunately you can not see data points with the same coordinates)
        "with_size" returns plotly plot, the size of a dot indicate sthe number of datapoints at this position (but hovering will only display info for one datapoint)
        "hvplot" returns hvplot with sizing as in "overlap_by_side" and hovering displays info for all nodes, but not as well designed as a plotly plot

    Return: plotly.graph_objs._figure.Figure
    ----------
    Plotly express scatter plot showing the PCA values
    """
    # check input
    if not ((dimensions == 2) or (dimensions == 3)):
        raise ValueError(
            "Dimensions is only allowed to take on values 2 and 3 (for visualization)"
        )
    if clusters < 1:
        raise ValueError("The number of clusters must be at least 1")
    if not ((method == "Correlation") or (method == "Differences")):
        raise ValueError(
            "The parameter method can only be 'Correlation' or 'Differences'"
        )
    if method == "Differences":
        if base_condition == "":
            base_condition = conditions[0]
        if base_condition in conditions:
            conditions = conditions[:]
            conditions.remove(base_condition)
    if not (display == "regular" or display == "with_size" or display == "hvplot"):
        raise ValueError(
            f'The display parameter value {display} is invalid. Valid display values are "regular", "with_size" and "hvplot"'
        )
    if display == "hvplot" and dimensions == 3:
        raise ValueError(
            'For the display type "hvplot" only a 2 dimensional version is available.'
        )

    ## calculate data according to chosen method

    # for correlation method
    # calculate pearson correlation coefficients
    if method == "Correlation":
        data = flux_coupling_matrix(
            fba, conditions, remove_unchanging=remove_unchanging, abs_values=abs_values
        )  # create correlation matrix
        ## remove and print all reactions for which the pearson correlation is NaN
        flag = True
        for index, el in data.isna().all(axis=1).items():
            if el == True:
                if flag:
                    print(
                        "The following reactions were removed beacuse the pearson correlation was Nan (probably because there was no change in flux)"
                    )
                    flag = False
                print(index)
        data.dropna(axis=0, how="all", inplace=True)
        data.dropna(axis=1, how="all", inplace=True)

    # for difference method
    if method == "Differences":
        data = pd.DataFrame(columns=conditions)  # difference matrix
        for c in conditions:
            data[c] = abs(fba[c] - fba[base_condition])
        data.index = fba["ID"]

    ## PCA    (this part is adapted from Aabhas work)
    pca = PCA(n_components=dimensions)
    pca_data = pca.fit_transform(data)

    total_var = pca.explained_variance_ratio_.sum() * 100
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100
    pc3_var = 0
    if dimensions == 3:
        pc3_var = pca.explained_variance_ratio_[2] * 100

    # create a DataFrame with the PCA results
    dim = ["PC1", "PC2"]
    if dimensions == 3:
        dim.append("PC3")
    pca_df = pd.DataFrame(data=pca_data, columns=dim)

    ## KMeans Clustering
    kmeans = KMeans(n_clusters=clusters, random_state=0)
    kmeans.fit(pca_data)
    labels = kmeans.labels_
    # centers = kmeans.cluster_centers_

    ## assemble dataframe wih all relevant data
    all_conditions = conditions[:]
    if not (base_condition == ""):
        all_conditions.append(base_condition)
        base_condition = "(to " + base_condition + ") "
    fba_infos = fba[np.bitwise_or.reduce([fba["ID"] == name for name in data.index])]
    fba_infos = fba_infos.reset_index(drop=True)
    pca_df_with_infos = pd.concat(
        [
            fba_infos[["ID", "Reaction", "Name", "Std_dev"]],
            pca_df,
            pd.Series(labels),
            fba_infos[all_conditions],
        ],
        axis=1,
    )
    if dimensions == 2:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "Cluster",
        ] + all_conditions
    else:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "PC3",
            "Cluster",
        ] + all_conditions

    ## get sizing info for nodes
    if not (display == "regular"):
        pca_df_with_infos["PC1"] = pca_df_with_infos["PC1"].round(6)
        pca_df_with_infos["PC2"] = pca_df_with_infos["PC2"].round(6)
        pca_df_with_infos["Count"] = pca_df_with_infos.groupby(["PC1", "PC2"])[
            "ID"
        ].transform("count")
        pca_df_with_infos["Node_size"] = [
            c if c < 40 else 40 for c in pca_df_with_infos["Count"]
        ]
    if display == "hvplot":
        pca_df_with_infos["Node_size"] = pca_df_with_infos["Node_size"] * 100

    ## create and return plot
    if display == "hvplot":
        return pca_df_with_infos.hvplot.scatter(
            x="PC1",
            y="PC2",
            size="Node_size",
            by="Cluster",
            hover_cols=["ID", "Name", "Count", "Std_dev"] + conditions,
            hover=True,
            alpha=0.6,
            title=f"{method} {base_condition}PCA colored by K-means clustering with {clusters} clusters\n"
            + f"Total explained variance: {total_var:.2f}%\n"
            + f"Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
        )
    if dimensions == 2:
        if display == "regular":
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                color="Cluster",
                title=f"<b>{method} {base_condition}PCA colored by K-means clustering with {clusters} clusters</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                color="Cluster",
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by K-means clustering with {clusters} clusters</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
            # fig.update_traces(marker=dict(sizemin=5))
    else:
        if display == "regular":
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="Cluster",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                title=f"<b>{method} {base_condition}PCA colored by K-means clustering with {clusters} clusters</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="Cluster",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by K-means clustering with {clusters} clusters</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
    return fig


def pca_color_by_pathway(
    fba,
    conditions,
    pathway,
    pathway_dict,
    dimensions=2,
    method="Correlation",
    base_condition="",
    remove_unchanging=0.01,
    abs_values=True,
    display="with_size",
):
    """
    Function to perform PCA colored by presence of the reactions in a pathway

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of all the conditions that will be considered in the PCA

    pathway: str
        Name of the pathway that wil be used to color the plot

    pathway_dict: Dict[str -> List[str]]
        Dictionary that maps pathways to the reactions of the model belonging to the pathway

    dimensions: int
        (optional) number of principal components, can be 2 or 3, default is 2

    method: str
        (optional) The method to use the data, can be 'Correlation' or 'Differences', default is correlation
        - 'Correlation': for each reaction use a vector of pairwise correlations to all other reactions as data
        - 'Differences': for each reaction use a vector of differences in flux to the base_condition as data

    base_condition: str
        (optional, only used if method = correlation) name of the base condition (the difference from the other conditions will be calculated from this condition)
        default is the first condition in the list of conditions

    remove_unchanging: float
        (optional, only used if method="Correlation", default = 0.01)
        Reactions with a smaller Std_dev than this threshold will be removed from the plot

    abs_values: bool
        (optional, only used if method="Correlation", default = True)
        will consider absolute values when calculating the correaltions and remove reactions where the sign swaps across conditions

    display: str
        default: "regular", options: "regular", "with_size", "hvplot"
        the type of plot that will be returned
        "regular" returns normal plotly plot (unfortunately you can not see data points with the same coordinates)
        "with_size" returns plotly plot, the size of a dot indicate sthe number of datapoints at this position (but hovering will only display info for one datapoint)
        "hvplot" returns hvplot with sizing as in "overlap_by_side" and hovering displays info for all nodes, but not as well designed as a plotly plot

    Return: plotly.graph_objs._figure.Figure
    ----------
    Plotly express scatter plot showing the PCA values
    """
    # check input
    if not ((dimensions == 2) or (dimensions == 3)):
        raise ValueError(
            f"Dimensions is only allowed to take on values 2 and 3 (for visualization), current value is {dimensions}"
        )
    if not ((method == "Correlation") or (method == "Differences")):
        raise ValueError(
            "The parameter method can only be 'Correlation' or 'Differences'"
        )
    if method == "Differences":
        if base_condition == "":
            base_condition = conditions[0]
        if base_condition in conditions:
            conditions = conditions[:]
            conditions.remove(base_condition)
    if not (display == "regular" or display == "with_size" or display == "hvplot"):
        raise ValueError(
            f'The display parameter value {display} is invalid. Valid display values are "regular", "with_size" and "hvplot"'
        )
    if display == "hvplot" and dimensions == 3:
        raise ValueError(
            'For the display type "hvplot" only a 2 dimensional version is available.'
        )

    ## calculate data according to chosen method

    # for correlation method
    # calculate pearson correlation coefficients
    if method == "Correlation":
        data = flux_coupling_matrix(
            fba, conditions, remove_unchanging=remove_unchanging, abs_values=abs_values
        )  # create correlation matrix
        ## remove and print all reactions for which the pearson correlation is NaN
        flag = True
        for index, el in data.isna().all(axis=1).items():
            if el == True:
                if flag:
                    print(
                        "The following reactions were removed beacuse the pearson correlation was Nan (probably because there was no change in flux)"
                    )
                    flag = False
                print(index)
        data.dropna(axis=0, how="all", inplace=True)
        data.dropna(axis=1, how="all", inplace=True)

    # for difference method
    if method == "Differences":
        data = pd.DataFrame(columns=conditions)  # difference matrix
        for c in conditions:
            data[c] = abs(fba[c] - fba[base_condition])
        data.index = fba["ID"]

    ## PCA
    pca = PCA(n_components=dimensions)
    pca_data = pca.fit_transform(data)

    total_var = pca.explained_variance_ratio_.sum() * 100
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100
    pc3_var = 0
    if dimensions == 3:
        pc3_var = pca.explained_variance_ratio_[2] * 100

    # create a DataFrame with the PCA results
    dim = ["PC1", "PC2"]
    if dimensions == 3:
        dim.append("PC3")
    pca_df = pd.DataFrame(data=pca_data, columns=dim)

    ## identify reactions that are part of the pathway of interest
    fba_infos = fba[np.bitwise_or.reduce([fba["ID"] == name for name in data.index])]
    fba_infos = fba_infos.reset_index(drop=True)
    mask = [False] * len(fba_infos)
    for el in pathway_dict[pathway]:
        mask = np.bitwise_or.reduce([mask, fba_infos["ID"] == el])

    ## assemble dataframe wih all relevant data
    all_conditions = conditions[:]
    if not (base_condition == ""):
        all_conditions.append(base_condition)
        base_condition = "(to " + base_condition + ") "
    pca_df_with_infos = pd.concat(
        [
            fba_infos[["ID", "Reaction", "Name", "Std_dev"]],
            pca_df,
            pd.Series(mask),
            fba_infos[all_conditions],
        ],
        axis=1,
        ignore_index=True,
    )
    if dimensions == 2:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "In_pathway",
        ] + all_conditions
    else:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "PC3",
            "In_pathway",
        ] + all_conditions

    ## get sizing info for nodes
    if not (display == "regular"):
        pca_df_with_infos["PC1"] = pca_df_with_infos["PC1"].round(6)
        pca_df_with_infos["PC2"] = pca_df_with_infos["PC2"].round(6)
        pca_df_with_infos["Count"] = pca_df_with_infos.groupby(["PC1", "PC2"])[
            "ID"
        ].transform("count")
        pca_df_with_infos["Node_size"] = [
            c if c < 40 else 40 for c in pca_df_with_infos["Count"]
        ]
    if display == "hvplot":
        pca_df_with_infos["Node_size"] = pca_df_with_infos["Node_size"] * 100

    ## create and return plot
    if display == "hvplot":
        return pca_df_with_infos.hvplot.scatter(
            x="PC1",
            y="PC2",
            size="Node_size",
            by="In_pathway",
            hover_cols=["ID", "Name", "Count", "Std_dev"] + conditions,
            hover=True,
            alpha=0.4,
            title=f"{method} {base_condition}PCA colored by presence in pathway {pathway}\n"
            + f"Total explained variance: {total_var:.2f}%\n"
            + f"Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
        )
    if dimensions == 2:
        if display == "regular":
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                color="In_pathway",
                title=f"<b>{method} {base_condition}PCA colored by presence in pathway {pathway}</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                color="In_pathway",
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by presence in pathway {pathway}</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
        fig.update_layout(legend_title=f"In {pathway}")
    else:
        if display == "regular":
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="In_pathway",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                title=f"<b>{method} {base_condition}PCA colored by presence in pathway {pathway}</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="In_pathway",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by presence in pathway {pathway}</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
        fig.update_layout(legend_title=f"In {pathway}")
    return fig


def pca_color_by_list(
    fba,
    conditions,
    reaction_list,
    dimensions=2,
    method="Correlation",
    base_condition="",
    remove_unchanging=0.01,
    abs_values=True,
    display="with_size",
):
    """
    Function to perform PCA colored by presence of the reactions in a list

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        list of all the conditions that will be considered in the PCA

    reaction_list: List[str]
        List of rection IDs that will be colored red in the plot

    dimensions: int
        (optional) number of principal components, can be 2 or 3, default is 2

    method: str
        (optional) The method to use the data, can be 'Correlation' or 'Differences', default is correlation
        - 'Correlation': for each reaction use a vector of pairwise correlations to all other reactions as data
        - 'Differences': for each reaction use a vector of differences in flux to the base_condition as data

    base_condition: str
        (optional, only used if method = correlation) name of the base condition (the difference from the other conditions will be calculated from this condition)
        default is the first condition in the list of conditions

    remove_unchanging: float
        (optional, only used if method="Correlation", default = 0.01)
        Reactions with a smaller Std_dev than this threshold will be removed from the plot

    abs_values: bool
        (optional, only used if method="Correlation", default = True)
        will consider absolute values when calculating the correaltions and remove reactions where the sign swaps across conditions

    display: str
        default: "regular", options: "regular", "with_size", "hvplot"
        the type of plot that will be returned
        "regular" returns normal plotly plot (unfortunately you can not see data points with the same coordinates)
        "with_size" returns plotly plot, the size of a dot indicate sthe number of datapoints at this position (but hovering will only display info for one datapoint)
        "hvplot" returns hvplot with sizing as in "overlap_by_side" and hovering displays info for all nodes, but not as well designed as a plotly plot

    Return: plotly.graph_objs._figure.Figure
    ----------
    Plotly express scatter plot showing the PCA values
    """
    # check input
    if not ((dimensions == 2) or (dimensions == 3)):
        raise ValueError(
            "Dimensions is only allowed to take on values 2 and 3 (for visualization)"
        )
    if not ((method == "Correlation") or (method == "Differences")):
        raise ValueError(
            "The parameter method can only be 'Correlation' or 'Differences'"
        )
    if method == "Differences":
        if base_condition == "":
            base_condition = conditions[0]
        if base_condition in conditions:
            conditions = conditions[:]
            conditions.remove(base_condition)
    if not (display == "regular" or display == "with_size" or display == "hvplot"):
        raise ValueError(
            f'The display parameter value {display} is invalid. Valid display values are "regular", "with_size" and "hvplot"'
        )
    if display == "hvplot" and dimensions == 3:
        raise ValueError(
            'For the display type "hvplot" only a 2 dimensional version is available.'
        )

    ## calculate data according to chosen method

    # for correlation method
    # calculate pearson correlation coefficients
    if method == "Correlation":
        data = flux_coupling_matrix(
            fba, conditions, remove_unchanging=remove_unchanging, abs_values=abs_values
        )  # create correlation matrix
        ## remove and print all reactions for which the pearson correlation is NaN
        flag = True
        for index, el in data.isna().all(axis=1).items():
            if el == True:
                if flag:
                    print(
                        "The following reactions were removed beacuse the pearson correlation was Nan (probably because there was no change in flux)"
                    )
                    flag = False
                print(index)
        data.dropna(axis=0, how="all", inplace=True)
        data.dropna(axis=1, how="all", inplace=True)

    # for difference method
    if method == "Differences":
        data = pd.DataFrame(columns=conditions)  # difference matrix
        for c in conditions:
            data[c] = abs(fba[c] - fba[base_condition])
        data.index = fba["ID"]

    ## PCA
    pca = PCA(n_components=dimensions)
    pca_data = pca.fit_transform(data)

    total_var = pca.explained_variance_ratio_.sum() * 100
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100
    pc3_var = 0
    if dimensions == 3:
        pc3_var = pca.explained_variance_ratio_[2] * 100

    # create a DataFrame with the PCA results
    dim = ["PC1", "PC2"]
    if dimensions == 3:
        dim.append("PC3")
    pca_df = pd.DataFrame(data=pca_data, columns=dim)

    ## identify reactions that are part of the list
    fba_infos = fba[np.bitwise_or.reduce([fba["ID"] == name for name in data.index])]
    fba_infos = fba_infos.reset_index(drop=True)
    mask = [False] * len(fba_infos)
    for el in reaction_list:
        mask = np.bitwise_or.reduce([mask, fba_infos["ID"] == el])

    ## assemble dataframe wih all relevant data
    all_conditions = conditions[:]
    if not (base_condition == ""):
        all_conditions.append(base_condition)
        base_condition = "(to " + base_condition + ") "
    pca_df_with_infos = pd.concat(
        [
            fba_infos[["ID", "Reaction", "Name", "Std_dev"]],
            pca_df,
            pd.Series(mask),
            fba_infos[all_conditions],
        ],
        axis=1,
    )
    if dimensions == 2:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "In_list",
        ] + all_conditions
    else:
        pca_df_with_infos.columns = [
            "ID",
            "Reaction",
            "Name",
            "Std_dev",
            "PC1",
            "PC2",
            "PC3",
            "In_list",
        ] + all_conditions

    ## get sizing info for nodes
    if not (display == "regular"):
        pca_df_with_infos["PC1"] = pca_df_with_infos["PC1"].round(6)
        pca_df_with_infos["PC2"] = pca_df_with_infos["PC2"].round(6)
        pca_df_with_infos["Count"] = pca_df_with_infos.groupby(["PC1", "PC2"])[
            "ID"
        ].transform("count")
        pca_df_with_infos["Node_size"] = [
            c if c < 40 else 40 for c in pca_df_with_infos["Count"]
        ]
    if display == "hvplot":
        pca_df_with_infos["Node_size"] = pca_df_with_infos["Node_size"] * 100

    ## create and return plot
    if display == "hvplot":
        return pca_df_with_infos.hvplot.scatter(
            x="PC1",
            y="PC2",
            size="Node_size",
            by="In_list",
            hover_cols=["ID", "Name", "Count", "Std_dev"] + conditions,
            hover=True,
            alpha=0.4,
            title=f"{method} {base_condition}PCA colored by presence in reaction list\n"
            + f"Total explained variance: {total_var:.2f}%\n"
            + f"Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
        )
    if dimensions == 2:
        if display == "regular":
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                color="In_list",
                title=f"<b>{method} {base_condition}PCA colored by presence in reaction list</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                color="In_list",
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by presence in reaction list</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%",
            )
        fig.update_layout(legend_title=f"in_list")
    else:
        if display == "regular":
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="In_list",
                hover_name="ID",
                hover_data=["Name", "Std_dev"] + all_conditions,
                title=f"<b>{method} {base_condition}PCA colored by presence in reaction list</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
        else:  # display = with_size
            fig = px.scatter_3d(
                pca_df_with_infos,
                x="PC1",
                y="PC2",
                z="PC3",
                color="In_list",
                hover_name="ID",
                hover_data=["Name", "Count", "Std_dev"] + all_conditions,
                opacity=0.5,
                size="Node_size",
                size_max=40,
                title=f"<b>{method} {base_condition}PCA colored by presence in reaction list</b>"
                + f"<br>Total explained variance: {total_var:.2f}%"
                + f"<br>Explained variance by PC: PC1: {pc1_var:.2f}%, PC2: {pc2_var:.2f}%, PC3: {pc3_var:.2f}%",
            )
        fig.update_layout(legend_title=f"in_list")
    return fig


def summarize_pathways(id_list, id_dict):
    """
    Function to summarize pathway info for multiple reactions in a dictionary

    Parameters
    ----------
    id_list : List[str]
        List of reaction IDs for which the pathways will be summarized

    id_dict: Dict: str -> List[str]
        Dictionary that contains for each ID the list of pathways corresponding to that reaction (if available)

    Return: Dict: str -> List[str]
    ----------
    Dictionary that contains for each pathway all reactions ids of the id_list that are part of this pathway

    """
    pathway_dict = defaultdict(list)
    for id in id_list:
        if id in id_dict:
            for pw in id_dict[id]:
                pathway_dict[pw].append(id)
    for key, val in pathway_dict.items():
        print(f"{key}: {len(val)} -> {val}")
    return pathway_dict


def pathway_variability(fba, pathway_dict, ascending=False, num=10):
    """
    Function to calculate the average "Pathway variability" for each pathway

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    pathway_dict: Dict: str -> List[str]
        Dictionary that contains for each Pathway the list of reactions in the model that belong to this pathway

    ascending: bool
        [optional] False, if the list should be sorted from highest to lowest Std_dev
                   True, if the list should be sorted from lowest to highest Std_dev

    num: int
        [optional] The number of pahways that should be displayed

    Return: List[Tuple(str, float, int)]
    ----------
    List of the {num} most/least variable pathways, the "average std_dev" of the pathway and the number of reactions in the pathway
    """
    ## check input
    if num > len(pathway_dict):
        num = len(pathway_dict)

    relevant_ids = set(fba["ID"])
    avg_dev = []
    for p, ids in pathway_dict.items():
        count = 0
        sum = 0
        for i in ids:
            if i in relevant_ids:
                count += 1
                sum += fba["Std_dev"][fba[fba["ID"] == i].index[0]]
        if count == 0:
            continue
        avg_dev.append((p, sum / count, count))

    sorted_avg_dev = sorted(avg_dev, key=lambda x: x[1], reverse=ascending)
    for p, avg, count in sorted_avg_dev[0:num]:
        print(f"{p}\t\t -> {avg:.2f} (contains {count} reactions)")

    return sorted_avg_dev[0:num]


def flux_coupling_matrix(fba, conditions, remove_unchanging=0.01, abs_values=True):
    """
    Function to create a flux coupling matrix

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        List of the names of the different conditions

    remove_unchanging: float
        remove all reactions with a std_dev < remove_unchanging
        (-> because reactions that do not change across conditions are typically not interesting for flux coupling)

    abs_values: bool
        If set to True, absolute flux values will be used and rections that have both positive and negative flux values are removed


    Return: pd.DataFrame
    ----------
    Correlation matrix between the flux values of the reactions
    """
    unchanging_fba = fba[
        fba["Std_dev"] >= remove_unchanging
    ]  # remove reactions that do not change across conditions
    num_unchanging = len(fba) - len(unchanging_fba)
    if not num_unchanging == 0:
        print(
            f"{num_unchanging} reaction were removed, because their Std_dev is lower than {remove_unchanging}"
        )

    if (
        abs_values
    ):  # get absolute values and remove reactions with swapping signs if the option abs_values is selected
        num_swapping = len(unchanging_fba)
        swapping_fba = unchanging_fba[
            ~(
                (unchanging_fba[conditions] >= 0).all(axis=1)
                | (unchanging_fba[conditions] <= 0).all(axis=1)
            )
        ]
        unchanging_fba = unchanging_fba[
            (unchanging_fba[conditions] >= 0).all(axis=1)
            | (unchanging_fba[conditions] <= 0).all(axis=1)
        ]
        unchanging_fba.loc[:, conditions] = unchanging_fba[conditions].abs().to_numpy()
        num_swapping -= len(unchanging_fba)
        if not num_swapping == 0:
            print(
                f"The following {num_swapping} reactions were removed because they swap signs across conditions: {list(swapping_fba["ID"])}"
            )

    corr_matrix = unchanging_fba[conditions].T.corr()  # calculate correlation matrix
    corr_matrix.columns = unchanging_fba["ID"]
    corr_matrix.index = unchanging_fba["ID"]

    return corr_matrix


def flux_coupling_network(flux_coupling_matrix, pc=0.8):
    """
    Function to create a flux coupling network from a flux coupling matrix

    Parameters
    ----------
    flux_coupling_matrix: pandas:DataFrame
        correlation matrix between the flux values of the reactions, cab be created using flux_coupling_matrix(fba, conditions)
        default: will be automatically created with default parameters

    pc: float
        Pearson correlation threahold: 2 reactions in the network will be connected if their pairwise pearson correlation is above this treshold


    Return: networkx.classes.graph.Graph
    ----------
    Flux coupling graph that connects 2 recations if their pairwise pearson corr is above the threshold pc
    """
    if pc > 1 or pc < -1:
        raise ValueError(
            f"Invalid value {pc} for pc, the pc threshold must lie between -1 and 1"
        )

    nw = nx.Graph()
    nw.name = f"Flux coupling graph with pc >= {pc}"
    nw.graph["info"] = ""
    nw.add_nodes_from(flux_coupling_matrix.columns, color="blue")
    rows, cols = np.tril_indices(
        flux_coupling_matrix.shape[0], k=-1
    )  # get lower triangle of correlation matrix
    for i, j in zip(rows, cols):
        value = flux_coupling_matrix.iloc[i, j]
        if value >= pc:
            reaction1 = flux_coupling_matrix.index[i]
            reaction2 = flux_coupling_matrix.columns[j]
            nw.add_edge(reaction1, reaction2, weight=value)

    print(
        f"Created Flux Coupling network with pc threshold {pc}.\nNumber of nodes: {nw.number_of_nodes()}\nNumber of edges: {nw.number_of_edges()}\nNumber of connected components: {nx.number_connected_components(nw)}"
    )
    return nw


def color_network(nw, color="blue"):
    """
    Function to color all nodes in a chosen color (can for example be used to reset all nodes to the color blue after coloring them by community or pathway)

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network to be adapted

    color: string
        name of the color, default: "blue"


    Return: networkx.classes.graph.Graph
    ----------
    network with all nodes stoing the chosen color, can be visualized with visualize_network(nw, layout) or visualize_interactive_network(nw, fba, condiions, layout)
    """
    nw.graph["info"] = f" colored in {color}"
    for node in nw.nodes():
        nw.nodes[node]["color"] = color
    return nw


def color_network_by_pathway(nw, pathway, pathway_dict):
    """
    Function to color nodes of reactions belonging to a particular pathway red

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network to be adapted

    pathway: string
        KEGG name of the pathway

    pathway_dict: Dict[str -> List[str]]
        Dictionary that maps pathways to the reactions of the model belonging to the pathway

    Return: networkx.classes.graph.Graph
    ----------
    adapted network (with coloring information matching the pathway, can ve visualized with visualize_network(nw, layout))
    """
    nw.graph["info"] = f" colored by presence in {pathway}"
    for node in nw.nodes():
        if node in pathway_dict[pathway]:
            nw.nodes[node]["color"] = "red"
    return nw


""" helper function """


def _most_valuable_edge_with_weights(G):
    centrality = nx.edge_betweenness_centrality(G)
    for (
        edge,
        btwn,
    ) in (
        centrality.items()
    ):  # divide each betweeness measure by scaled weight to simulate multigraph
        u, v = edge
        weight = G[u][v].get("weight", None)
        if weight != 0:
            centrality[edge] = btwn / (weight * 10)
    return max(centrality, key=centrality.get)


def color_network_by_communities(nw, iterations, include_weights=False):
    """
    Function to color the network by communities determined by the Girvan Newman algorithm

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network to be adapted

    iterations: int
        number of iterations the Girvan Newman algorithm should be run

    include_Weights: bool
        default: False
        if set to True the weigtts of the edges are considered whe constructing communities

    Return: networkx.classes.graph.Graph
    ----------
    adapted network (with coloring information matching the communities, can ve visualized with visualize_network(nw, layout))
    """
    if iterations < 1:
        raise ValueError(
            f"{iterations} is not a valid value for iterations, the value must be positive"
        )

    # get communities with Girvan_Newman algorithm
    if include_weights:
        communities_generator = nx.community.girvan_newman(
            nw, most_valuable_edge=_most_valuable_edge_with_weights
        )
    else:
        communities_generator = nx.community.girvan_newman(nw)
    communities = next(communities_generator)
    for _ in range(iterations - 1):
        communities = next(communities_generator)

    nw.graph["info"] = f" colored by {len(communities)} communities"
    color_names = px.colors.qualitative.Light24 + list(
        mcolors.CSS4_COLORS.keys()
    )  # find more colors here: https://plotly.com/python/discrete-color/

    for i, community in enumerate(communities):
        if i > 157:
            raise ValueError(
                "There are too many communities to display (more than 157), try using less iterations of the Girvan Newman algorithm"
            )
        for node in community:
            nw.nodes[node]["color"] = color_names[i]

    return nw


def visualize_network(nw, layout="kamada_kawai"):
    """
    Function to visualize a flux coupling network

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network to visualize

    layout: string
        valid options: "kamada_kawai", "spring", "random", "circular"
        default: "kamda_kawai"
        layout of the graph


    Return: plotly.graph_objs._figure.Figure
    ----------
    Flux coupling network as plot
    """
    if layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(nw)
    elif layout == "spring":
        pos = nx.spring_layout(nw)
    elif layout == "random":
        pos = nx.random_layout(nw)
    elif layout == "circular":
        pos = nx.circular_layout(nw)
    else:
        raise ValueError(
            f"{layout} is not a valid layout. Valid layouts are kamada_kawai, spring, random and circular."
        )

    node_colors = [nw.nodes[node]["color"] for node in nw.nodes()]

    plt.figure(figsize=(10, 10))
    nx.draw(
        nw,
        pos,
        node_size=25,
        edge_color="gray",
        width=0.5,
        node_color=node_colors,
        alpha=0.5,
        with_labels=True,
        font_size=7,
    )
    plt.title(nw.name + nw.graph["info"])
    return plt


"""helper function"""


def _flux_info(fba, i, conditions):
    text = ""
    for c in conditions:
        text += f"<br>{c}: {fba[c][i]:.3f}"
    return text


def visualize_interactive_network(
    nw, fba, conditions, layout="kamada_kawai", edge_weights=True
):
    """
    Function to visualize a flux coupling network interactively (e.g. allows zooming, hovering over nodes)

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network to visualize

    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        List of the names of the different conditions

    layout: string
        valid options: "kamada_kawai", "spring", "random", "circular"
        default: "kamda_kawai"
        layout of the graph

    edge_weights: boolean
        default: True
        if set to True, the edge thickness will depend on the edge weight, so reactions with higher pearson correlation will be connected by a thicker line


    Return: plotly.graph_objs._figure.Figure
    ----------
    Flux coupling network as interactive plot
    """
    if layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(nw)
    elif layout == "spring":
        pos = nx.spring_layout(nw)
    elif layout == "random":
        pos = nx.random_layout(nw)
    elif layout == "circular":
        pos = nx.circular_layout(nw)
    else:
        raise ValueError(
            f"{layout} is not a valid layout. Valid layouts are kamada_kawai, spring, random and circular."
        )

    node_colors = [nw.nodes[node]["color"] for node in nw.nodes()]
    index = [fba[fba["ID"] == node].index[0] for node in nw.nodes()]
    node_info = [
        f"{fba["ID"][i]}<br>{fba["Name"][i]}{_flux_info(fba, i, conditions)}<br>Std_dev: {fba["Std_dev"][i]:.3f}"
        for i in index
    ]

    fig = go.Figure()
    edge_width = 0.35

    # add edges (this part is adapted from ChatGPT)
    for edge in nw.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        if edge_weights:
            scaled_weight = (
                edge[2].get("weight", 1) ** 5
            )  # scale edge thickness by the correlation coefficient
            edge_width = max(0.1, float(scaled_weight))  # ensure minimum thickness

        fig.add_trace(
            go.Scatter(
                x=[x0, x1, None],  # separate edges with None
                y=[y0, y1, None],
                mode="lines",
                line=dict(width=edge_width, color="gray"),
            )
        )

    # add nodes (this part is adapted from ChatGPT)
    node_x, node_y = zip(*[pos[node] for node in nw.nodes()])
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            marker=dict(
                size=9, color=node_colors, opacity=0.5
            ),  # or use size 15 for bigger nodes, normal size is 9
            text=list(nw.nodes),
            textposition="middle center",
            hovertext=node_info,
            hoverinfo="text",
        )
    )

    fig.update_layout(
        title=dict(text=nw.name + nw.graph["info"], x=0.5),
        title_font_size=18,
        showlegend=False,
        hovermode="closest",
        width=1000,
        height=1000,
        plot_bgcolor="white",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=100, r=100, t=70, b=0),
        font=dict(size=10),
    )  # or use size 20 (or 12) for a bigger font size, normal size is 10
    return fig


def get_community_list(nw):
    """
    Function to get a dictionary of the communities determined by the Girvan_Newman algorithm and plot a legend for the network based on the communities

    Parameters
    ----------
    nw : networkx.classes.graph.Graph
        network

    Return: dict(str -> List["str"])
    ----------
    Dictionary that maps each community (named by color it would have in the network) to the reactions assigned to this community, also gets displayed as a legend for the network
    """
    community_map = defaultdict(list)  # store color info in node for each community

    for node in nw.nodes():
        community_map[nw.nodes[node]["color"]].append(node)

    legend_text = []  # create legend
    for index, (key, val) in enumerate(community_map.items(), start=1):
        dots = ""
        if len(val) > 10:
            dots = "..."
        legend_text.append(
            f"{index}. ({key}) with {len(val)} elements: {val[0:9]}{dots}"
        )
    legend_color = list(community_map.keys())

    fig = (
        go.Figure()
    )  # create empty dummy figure (bc I don' know how else to create a legend without a figure)
    for text, color in zip(legend_text, legend_color):
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=color),
                name=text,
            )
        )
    fig.update_layout(
        showlegend=True,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        legend=dict(x=0.5, y=0.5, xanchor="center", yanchor="middle", orientation="v"),
        width=1000,
        height=300,
        plot_bgcolor="white",
        margin=dict(l=0, r=0, t=0, b=0),
    )

    fig.show()  # display legend
    return community_map


def get_color_map(com_list):
    """
    Function to get a dictionary that maps the community number to the community color

    Parameters
    ----------
    com_list : dict(str -> List["str"])
        community map, can be created with get_community_list(nw)

    Return: dict(int -> str)
    ----------
    Dictionary that maps each community number to its color
    """
    color = {}
    for i, c in enumerate(com_list.keys(), start=1):
        color[i] = c
    return color


def community_bar_plots(fba, conditions, com_list):
    """
    Function to create bar plots with the flux value means of every condition for every community (idea: recognize patterns in the communities)

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        List of the names of the different conditions

    com_List: dict(str -> List["str"])
        dictionary that maps commuity color to a list of reaction IDs that are in the community, can be created with get_community_list(nw, iterations)

    Return: plotly.graph_objs._figure.Figure
    ----------
    Figure with one subplot for each community showing a bar plot of the mean of the flux values for each condition
    """
    # create data of the flux means for each community
    flux_dist = pd.DataFrame(
        {"Condition": conditions}
    )  # dataframe with just the flux values
    for com in com_list.keys():
        indices = [fba[fba["ID"] == id].index[0] for id in com_list[com]]
        means = []
        for c in conditions:
            c_fluxes = [abs(fba[c][i]) for i in indices]
            means.append(np.average(c_fluxes))
        flux_dist[com] = means

    num_cols = 3  # create plot
    rows = ceil(len(com_list.keys()) / num_cols)
    fig = make_subplots(
        rows=rows,
        cols=num_cols,
        subplot_titles=[
            str(i) + f" with {len(com_list[c])} reactions"
            for i, c in enumerate(com_list.keys(), start=1)
        ],
    )

    for i, com in enumerate(com_list.keys()):
        row = (i // num_cols) + 1  # calculate row index
        col = (i % num_cols) + 1  # calculate column index
        fig.add_trace(
            go.Bar(
                x=flux_dist["Condition"],
                y=flux_dist[com],
                marker=dict(color=com),
                opacity=0.5,
            ),
            row=row,
            col=col,
        )

    fig.update_layout(
        title="Flux mean per condition of each community ",
        height=rows * 200,
        width=1000,
        showlegend=False,
    )
    return fig


def community_box_plots(fba, conditions, com_list):
    """
    Function to create box plots with the flux value means and std_devs of every condition for every community (idea: recognize patterns in the communities)

    Parameters
    ----------
    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    conditions: List[str]
        List of the names of the different conditions

    com_List: dict(str -> List["str"])
        dictionary that maps commuity color to a list of reaction IDs that are in the community, can be created with get_community_list(nw, iterations)

    Return: plotly.graph_objs._figure.Figure
    ----------
    Figure with one subplot for each community showing a box plot of the flux value distribution (with mean and std_dev) for each condition
    """
    num_cols = 3  # create plot
    num_rows = -(-len(com_list.keys()) // num_cols)
    fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        subplot_titles=[
            "community " + str(i) + f" with {len(com_list[c])} reactions"
            for i, c in enumerate(com_list.keys(), start=1)
        ],
    )

    for i, com in enumerate(com_list.keys()):  # create sub plots
        row = i // num_cols + 1
        col = i % num_cols + 1
        indices = [
            fba[fba["ID"] == id].index[0] for id in com_list[com]
        ]  # create df with flux values for each community
        com_df = fba[fba.index.isin(indices)]
        com_df.loc[:, conditions] = com_df[conditions].abs()

        for c in conditions:  # add box for each condition
            fig.add_trace(
                go.Box(
                    y=com_df[c],
                    name=c,
                    boxmean=True,
                    marker=dict(color=com),
                    opacity=0.6,
                ),
                row=row,
                col=col,
            )

    fig.update_layout(
        title="Flux value distribution for the different communities",
        height=200 * num_rows,
        width=1000,
        showlegend=False,
    )
    return fig


def community_pw_summary(
    com_list, id_dict, pw_dict, lower_bound=5, upper_bound=100, cut_off=0.1
):
    """
    Function to get pathway summary info on the communities
    note: if you don't want any constraints set lower_bound = 0, upper_bound = 10000, cut_off = 0

    Parameters
    ----------
    com_list: dict(str -> List["str"])
        dictionary that maps commuity color to a list of reaction IDs that are in the community, can be created with get_community_list(nw, iterations)

    id_dict: Dict: str -> List[str]
        Dictionary that contains for each ID the list of pathways corresponding to that reaction (if available)

    pathway_dict: Dict[str -> List[str]]
        Dictionary that maps pathways to the reactions of the model belonging to the pathway

    lower_bound: float
        num of reactions that are pathway needs to at least have in order to be considered

    upper_bound: float
        num of reactions that are pathway is allowed to have at most in order to be considered

    cut_off: float
        minimum percentage of reactions in the pathway that are in the community so that this pathway will be printed

    Return: void
    ----------
    Prints pathway infos
    """
    count = 0
    for color, reactions in com_list.items():
        count += 1
        print(f"Community {count} ({color}, {len(reactions)} elements):")
        with contextlib.redirect_stdout(io.StringIO()):
            com_dict = summarize_pathways(reactions, id_dict)
        info = []
        for key, val in com_dict.items():
            if len(pw_dict[key]) > lower_bound and len(pw_dict[key]) < upper_bound:
                if (len(val) / len(pw_dict[key])) > cut_off:
                    info.append((key, len(val), len(val) / len(pw_dict[key])))
        info.sort(key=lambda x: -x[2])
        for name, num, per in info:
            print(f"{name}: {num} -> {per:.2f} %")
        print("\n-------------------\n")


def community_sankey_diagram(
    com_list,
    pw_list,
    fba,
    id_dict,
    diagram_height=500,
    pw_hierarchy={},
    unchanging=0.01,
    text_size=14,
):
    """
    Function to create Sankey diagram that shows the distribution of reactions of specific pathways into the communities

    Parameters
    ----------
    com_list: dict(str -> List[str])
        dictionary that maps commuity color to a list of reaction IDs that are in the community, can be created with get_community_list(nw, iterations)

    pw_list: List[str]
        list of KEGG pathway names that will be considered in the diagram

    fba : pandas.DataFrame
        DataFrame of recation fluxes for different conditions, can be created using the build_dataframe(model, solutions, conditions) function

    id_dict: Dict: str -> List[str]
        Dictionary that contains for each ID the list of pathways corresponding to that reaction (if available)

    diagram_height: int (typically in range 100 to 1000)
        default: 500
        pixel height of the resulting sankey diagram

    pw_hierarchy: dict(str -> List[str])
        optional, dictionary mapping a pathway overview to all the subpathways included in the overview pathways
        intended use: put in overview of KEGG pathway hierarchy

    unchanging: float
        default: 0.01
        all reactions with a std_dev lower than this threshold will be conidered as unchanging

    text_size: int
        default: 14
        font size of the text in the plot


    Return: plotly.graph_objs._figure.Figure
    ----------
    Sankey diagram
    """
    if pw_hierarchy:
        for pathway in pw_list:
            if not pathway in pw_hierarchy.keys():
                raise ValueError(
                    f"The pw_list must be a subset of the pw_hierarchy. The entry {pathway} was not found in the pw_hierarchy."
                )

    # create data
    reactions_per_pw = [0] * len(pw_list)  # calculate num reactions per pathway
    for id in fba["ID"]:
        if id in id_dict:
            for pw, pathway in enumerate(pw_list):
                if pw_hierarchy:
                    for sub_pw in pw_hierarchy[pathway]:
                        if sub_pw in id_dict[id]:
                            reactions_per_pw[pw] += 1
                else:
                    if pathway in id_dict[id]:
                        reactions_per_pw[pw] += 1

    opacity = 0.5
    sources, targets, values, colors, labels = [], [], [], [], []

    unchanging_list = fba[fba["Std_dev"] <= unchanging][
        "ID"
    ]  # identify pathways corresponding to unchanging reactions
    unchanging_pws = [[] for _ in range(len(pw_list))]
    for reaction in unchanging_list:
        if reaction in id_dict:
            for pw, pathway in enumerate(pw_list):
                if pw_hierarchy:
                    for sub_pw in pw_hierarchy[pathway]:
                        if sub_pw in id_dict[reaction]:
                            unchanging_pws[pw].append(reaction)
                            break
                else:
                    if pathway in id_dict[reaction]:
                        unchanging_pws[pw].append(reaction)
    for pw, reactions in enumerate(unchanging_pws):
        if len(reactions) > 0:
            sources.append(pw)
            targets.append(len(pw_list))
            values.append(len(reactions))
            colors.append(f"rgba{mcolors.to_rgba("gray", opacity)}")
            r_labels = ""
            for r, reaction in enumerate(reactions):
                if r == 0:
                    r_labels += f"{reaction}"
                elif r > 15:
                    r_labels += ", ..."
                    break
                else:
                    r_labels += f", {reaction}"
            labels.append(f"reactions: {r_labels}")

    for c, (com, reactions) in enumerate(
        com_list.items(), start=len(pw_list) + 1
    ):  # identify pathways correspnding tpo communities
        community_pws = [[] for _ in range(len(pw_list))]
        for reaction in reactions:
            if reaction in id_dict:
                for pw, pathway in enumerate(pw_list):
                    if pw_hierarchy:
                        for sub_pw in pw_hierarchy[pathway]:
                            if sub_pw in id_dict[reaction]:
                                community_pws[pw].append(reaction)
                                break
                    else:
                        if pathway in id_dict[reaction]:
                            community_pws[pw].append(reaction)
        for pw, reactions in enumerate(community_pws):
            if len(reactions) > 0:
                sources.append(pw)
                targets.append(c)
                values.append(len(reactions))
                colors.append(
                    f"rgba{mcolors.to_rgba(com, opacity)}"
                )  # make colors more transparent
                r_labels = ""
                for r, reaction in enumerate(reactions):
                    if r == 0:
                        r_labels += f"{reaction}"
                    elif r > 15:
                        r_labels += ", ..."
                        break
                    else:
                        r_labels += f", {reaction}"
                labels.append(
                    f"reactions: {r_labels}"
                )  # add reaction list to hoverdata

    # create plot
    fig = go.Figure(
        data=[
            go.Sankey(
                valueformat=".0f",
                valuesuffix=" reactions",
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.2),
                    label=pw_list
                    + ["unchanging reactions"]
                    + [str(i) for i in range(1, len(com_list.keys()) + 1)],
                    color=["gray"] * (len(pw_list) + 1) + list(com_list.keys()),
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=colors,
                    label=labels,
                ),
            )
        ]
    )
    fig.update_layout(
        title_text="Sankey Diagram of pathway distribution across communities",
        font_size=text_size,
        height=diagram_height,
    )
    return fig


if __name__ == "__main__":
    print(
        "Currently only intended as import into jupyter notebooks.\nThe following functions are available:\n \
          - build_dataframe(model, solutions, conditions)\n \
          - scatter_plot_fluxes(fba, conditions)\n \
          - correlation_matrix(fba, conditions[opt:, color])\n \
          - filter_dataframe(fba, conditions)\n \
          - bar_plot_flux_variation(fba[opt:, start, end])\n \
          - hierarchical_clustering_std_dev(fba[opt: ,start, end, d_max]\n \
          - pca(fba, base_condition, other_conditions[opt:, clusters, dimensions)\n \
          - normalize_dataframe_rows(fba, conditions[, opt: method])\n \
          - normalize_dataframe_cols(fba, conditions, obj_values)\n \
          - pca_color_pathways(fba, base_condition, other_conditions, pathway, pathway_dict, [opt: dimensions=2])\n \
          - summarize_pathways(id_list, id_dict)\n \
          - flux_coupling_matrix(fba, conditions, remove_unchanging = 0.01, abs_values = True)\n \
          - flux_coupling_network(flux_coupling_matrix = pd.DataFrame(), pc = 0.8)\n \
          - color_network_by_pathway(nw, pathway, pathway_dict)\n \
          - color_network_by_communities(nw, iterations)\n \
          - visualize_network(nw, [opt: layout])\n \
          - get_community_list(nw, iterations)\n \
          - community_bar_plots(fba, conditions, com_list)\n \
          - community_box_plots(fba, conditions, com_list)\n \
          - community_pw_summary(com_list, id_dict, pw_dict, [opt: lower_bound, upper_bound, cut_off])\n \
          - community_sankey_diagram(com_list, pw_list, fba, id_dict, [opt: diagram_height, pw_hierarchy])\n \
          - color_network(nw, [opt: color])"
    )
