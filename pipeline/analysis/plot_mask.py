"""
Script file to create colour mask, colour mask plot and one hot encoding mask 

Args: 
    - input_fn: the input mask file in tiff 
    - type_asn_csv: the csv containing the cell type assignments 
    - pair_csv: the csv containing the ObjectNumber field and ObjectNumber_renamed field
    - output_dir: the directory path of where we want to store the three output files 
    
Flags: 
    --make_dir or -md : to make a directory containing all output file
    --colour_mask or -cm : to output colour mask file
    --colour_plot or -cp : to output colour plot file
    --one_hot_enc_mask or -hm : to output one hot encoding mask file 
    --mask_and_annotation or -ma : to output a plot file which has the mask file and the annotation file side by side 

Output:
    - <core_name>_colour_plot.pdf 
    - <core_name>_colour_mask.tiff 
    - <core_name>_one_hot_encoding_mask.tiff 

Example: 
    $ input_fn="/Users/sunyunlee/Desktop/ltri/share/projects/imc-2020/data-raw/jackson-2020/OMEnMasks/Basel_Zuri_masks/BaselTMA_SP43_25.8kx22ky_10500x6500_8_20170928_114_115_X4Y8_262_a0_full_maks.tiff"
    
    $ type_asn_fn="/Users/sunyunlee/Desktop/ltri/share/projects/imc-2020/output/squirrel/astir_assignments/basel_astir_assignments.csv"
    $ type_asn_fn="/Users/sunyunlee/Desktop/ltri/share/projects/imc-2020/output/squirrel/astir_assignments/zurich1_astir_assignments.csv"

    $ pair_fn="/Users/sunyunlee/Desktop/ltri/share/projects/imc-2020/data-raw/jackson-2020/spatial-locations/Basel_SC_locations.csv"
    $ pair_fn="/Users/sunyunlee/Desktop/ltri/share/projects/imc-2020/data-raw/jackson-2020/spatial-locations/Zurich_SC_locations.csv"

    $ output_dir="/Users/sunyunlee/Desktop"

    $ python3 plot_mask.py $input_fn $type_asn_fn $pair_fn $output_dir -md -cm -cp -hm
    python3 plot_mask.py $input_fn $type_asn_fn $pair_fn $output_dir -cm 
    python3 plot_mask.py $input_fn $type_asn_fn $pair_fn $output_dir -pm
"""

""" 
Example: 
    $ input_fn=/Users/sunyunlee/Desktop/masks/ZTMA208_slide_28.23kx22.4ky_7000x7000_5_20171115_114_71_Ay15x2_397_a0_full_maks.tiff
    $ type_asn_fn=/Users/sunyunlee/Documents/astir_campbell_lab/astir-repos/img-data-scripts/data-raw/zurich1_astir_assignments.csv
    $ pair_fn=/Users/sunyunlee/Documents/astir_campbell_lab/astir-repos/img-data-scripts/data-raw/zurich1_SC_locations.csv
    $ output_dir="/Users/sunyunlee/Desktop"
    python3 plot_mask.py $input_fn $type_asn_fn $pair_fn $output_dir -md -cm -pm
"""

import argparse
import os
import imghdr
import pandas as pd
import numpy as np 
import tifffile as tf
from PIL import ImageColor
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager as font_manager


""" function definitions """
CONSTANTS = {
    "PALETTE": [
        "#8B5B42", "#AF4EA9", "#FFB60A", "#0AC694", "#0024DD", "#6CC1FF",
        "#0496FF", "#E11E00"
    ],
    "CELL_NAMES": [
        "Stromal", "B cells", "T cells", "Macrophage", "Epithelial (basal)", "Epithelial (luminal)", "Epithelial (other)", "Endothelial"
    ]
}

CONSTANTS["CELL_CATEGORIES"] = CONSTANTS["CELL_NAMES"] + ["Other", "Unknown", "Not cell"]

def extract_sample_name_from_filename(df_asn: pd.DataFrame, filename: str):
    """ Takes in a filename and extracts sample name

    :param filename: the name of the file to extract sample name from 
    :param type: basel or zur1

    :return: the sample name 
    """
    fn = filename.split("_")
    core_name = fn[1] + "_" + fn[7] + "_" + fn[8]

    core_asns = df_asn[df_asn.index.str.contains(core_name)]

    # Check if there is at least one core with core_name
    if not core_asns.shape[0] > 0:
        raise Exception("Bad core name: {} is not a core in the given assignment file".format(core_name))
    
    # Return core_name according to whether or not there are several cores with the same name 
    core = core_asns.index.values[0]
    if len(core.split("_")) == 5: 
        return fn[1] + "_" + fn[7] + "_" + fn[8]
    elif len(core.split("_")) == 6: 
        return fn[1] + "_" + fn[7] + "_" + fn[8] + "_" + fn[9]


def most_likely_celltype(row, threshold, cell_types):
    """ Given a row of the assignment matrix return the most likely cell type

    :param row: 
    :param threshold: 
    :param cell_types: 

    :return: 
    """
    row = row.to_numpy()
    max_prob = np.max(row)

    if max_prob < threshold:
        return "Unknown"

    return cell_types[np.argmax(row)]


def get_cell_type_assignments(df_asn: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    """ Get cell type assignments 

    :param df_asn: pandas dataframe containing assignments of all cells for each cell type
    :param sample_name: name of the sample to extract

    :return: pandas dataframe containing assignments of just sample_name sample 
    """
    cell_types = list(df_asn.columns)
    cell_type_assignments = df_asn.apply(
        most_likely_celltype,
        axis=1,
        threshold=0.7,
        cell_types=cell_types
    )
    df_type_asn = pd.DataFrame(cell_type_assignments)
    df_type_asn.columns = ["cell_type"]
    df_type_asn = df_type_asn[df_type_asn.index.str.contains(sample_name)]
    df_type_asn.index = [name.split("_")[-1] for name in df_type_asn.index]

    return df_type_asn


def get_cell_type_probabilities(df_asn: pd.DataFrame, sample_name: str) -> pd.DataFrame: 
    """ Get probabilities of cell type assignments cells in sample_name
    """
    df_core_asn = df_asn[df_asn.index.str.contains(sample_name)]
    df_core_asn = df_core_asn.max(axis=1)
    df_core_asn.index = [name.split("_")[-1] for name in df_core_asn.index]
    
    return df_core_asn


def map_mask_to_cell_ids(df_pairs: pd.DataFrame, df_mask: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    """ Replace cell ids in mask to their actual cell ids 

    :param df_mask: mask dataframe to modify
    :param sample_name: 

    :return: mask file in panda dataframe
    """
    # df_pairs = pd.read_csv(INPUT_FILES["MASK_TO_CELL_ID_PAIRS_CSV_BASEL"])
    # df_pairs = pd.read_csv(INPUT_FILES["MASK_TO_CELL_ID_PAIRS_CSV_ZUR1"])
    df_pairs = df_pairs[df_pairs["core"].str.contains(sample_name)]
    # print(df_pairs)

    if not (np.sort(df_pairs["ObjectNumber"].to_numpy()) == 
    np.sort(np.delete(np.unique(df_mask), np.argwhere(np.unique(df_mask) == 0)))).all():
        raise Exception("Mask files have unexpected cell ids for sample_name {}".format(sample_name))

    from_id = list(map(int, df_pairs["ObjectNumber"]))
    to_id = list(map(int, df_pairs["ObjectNumber_renamed"]))

    df_mask = df_mask.replace(to_replace=from_id, value=to_id)
    return df_mask


def map_mask_id_to_cell_types(df_mask: pd.DataFrame, df_asn: pd.DataFrame) -> pd.DataFrame: 
    """ Replace cell ids in mask dataframe to their cell names 

    :param df_mask: mask dataframe to modify 
    :
    """
    nums_to_names = df_asn.to_dict()["cell_type"]
    nums_to_names["0"] = "Not cell"

    df_mask = df_mask.replace(to_replace=list(map(int, nums_to_names.keys())), value=list(nums_to_names.values()))

    return df_mask 

def map_mask_id_to_probabilities(df_mask: pd.DataFrame, df_prob: pd.DataFrame) -> pd.DataFrame:
    """ Replace cell ids in mask dataframe to their cell type probabilities

    :param df_mask: mask dataframe to modify 
    :param df_prob: cell type probabilities
    """
    nums_to_prob = df_prob.to_dict()
    nums_to_prob["0"] = 0

    df_mask = df_mask.replace(to_replace=list(map(int, nums_to_prob.keys())), value=list(nums_to_prob.values()))

    return df_mask


def map_mask_cell_names_to_colors(df_mask: pd.DataFrame) -> np.array:
    """ Returns a mask pandas dataframe that 

    :param df_mask: mask file containing cell type names

    :return: three color channel numpy array 
    """
    palette_rgba = [ImageColor.getcolor(color, "RGBA") for color in CONSTANTS["PALETTE"]]
    map_type_color = dict(zip(CONSTANTS["CELL_NAMES"], palette_rgba))

    # Add "Other" and "Unknown" categories to cell types 
    map_type_color["Other"] = ImageColor.getcolor("#B3B3B3", "RGBA") # Grey 70
    map_type_color["Unknown"] = ImageColor.getcolor("#E5E5E5", "RGBA") # Grey 90

    # map_type_color["Not cell"] = [255, 255, 255, 255] # White background
    map_type_color["Not cell"] = [0, 0, 0, 255] # Black background

    np_mask = df_mask.to_numpy()
    np_rgba_mask = np.array([np.zeros_like(np_mask) for i in range(4)])

    for type, color in map_type_color.items(): 
        for i, channel in enumerate(color):
            np_rgba_mask[i] = np.where(df_mask == type, channel, np_rgba_mask[i])
    
    return np_rgba_mask.transpose(1, 2, 0).astype(np.uint8)


def save_mask_as_tiff(np_mask: np.array, filename: str) -> None:
    """ Save mask np array to tiff file  

    :param np_mask: numpy array of channel mask to save as tiff
    """
    tf.imsave(filename, np_mask)


def map_mask_cell_names_to_one_hot_encoding(df_mask: pd.DataFrame) -> np.array:
    """ Returns a mask pandas dataframe that 

    :param df_mask: mask file containing cell type names

    :return: three color channel numpy array 
    """
    cell_names = CONSTANTS["CELL_CATEGORIES"]
    # map_type_num = dict(zip(cell_names, range(0, len(cell_names))))
    np_mask = df_mask.to_numpy()
    np_ch_mask = np.array([np.zeros_like(np_mask) for i in range(len(cell_names))])

    for i, type in enumerate(cell_names): 
        np_ch_mask[i] = np.where(df_mask == type, 1, np_ch_mask[i])
    
    return np_ch_mask


def test_hot_enc_correct(np_ch_mask, np_type_name_mask):
    cell_names = CONSTANTS["CELL_CATEGORIES"]

    if not len(np_ch_mask) == len(cell_names):
        return False
    if not (np.unique(np_ch_mask) == np.array([0, 1])).all():
        return False 

    for i in range(len(np_ch_mask)):
        np_curr_ch_mask = np_ch_mask[i]
        if not np_curr_ch_mask.shape == np_type_name_mask.shape:
            return False 
        
        hot_enc_indices = np.where(np_curr_ch_mask == 1)
        mask_indices = np.where(np_type_name_mask == cell_names[i])

        if not (len(hot_enc_indices) == 2 and len(mask_indices) == 2):
            return False 
        if not (hot_enc_indices[0] == mask_indices[0]).all(): 
            return False 
        if not (hot_enc_indices[1] == mask_indices[1]).all():
            return False
        
    return True


def plot_mask(np_rgba_mask: np.array, core_name: str, output_dir: str) -> None:
    palette_rgba = [ImageColor.getcolor(color, "RGBA") for color in CONSTANTS["PALETTE"]]
    map_type_color = dict(zip(CONSTANTS["CELL_NAMES"], palette_rgba))

    # Add "Other" and "Unknown" categories to cell types 
    map_type_color["Other"] = ImageColor.getcolor("#B3B3B3", "RGBA") # Grey 70
    map_type_color["Unknown"] = ImageColor.getcolor("#E5E5E5", "RGBA") # Grey 90

    # map_type_color["Not cell"] = [255, 255, 255, 255] # White background
    map_type_color["Not cell"] = [0, 0, 0, 255] # Black background
    # create patches as legend
    patches =[mpatches.Patch(color=np.array(map_type_color[type]) / 255,label=type) for type in map_type_color]
    f = plt.figure()
    ax = plt.subplot(111)
    plt.axis("off")
    plt.imshow(np_rgba_mask)
    plt.title(core_name)
    # ax.legend(handles=patches, borderaxespad=0., loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=4)
    font_legend = font_manager.FontProperties(family="Arial", size=7)
    plt.legend(handles=patches, loc=10, bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize="x-small", prop=font_legend, frameon=False)
    output_filename = core_name + "_color_plot.pdf"
    f.savefig(os.path.join(output_dir, output_filename), bbox_inches='tight', dpi=5000)
    plt.close()
# def plot_mask_and_annotation(np_mask: np.array, np_ann: np.array, core_name: str, output_dir: str) -> None:
#     palette_rgba = [ImageColor.getcolor(color, "RGBA") for color in CONSTANTS["PALETTE"]]
#     map_type_color = dict(zip(CONSTANTS["CELL_NAMES"], palette_rgba))

#     # Add "Other" and "Unknown" categories to cell types 
#     map_type_color["Other"] = ImageColor.getcolor("#B3B3B3", "RGBA") # Grey 70
#     map_type_color["Unknown"] = ImageColor.getcolor("#E5E5E5", "RGBA") # Grey 90
    
#     ann_color = {}
#     ann_color["duct"] = ImageColor.getcolor("#FFFF00", "RGBA")
#     ann_color["not duct"] = ImageColor.getcolor("#0000FF", "RGBA")

#     # create patches as legend
#     fig, (ax1, ax2) = plt.subplots(1, 2)
    
#     plt1 = ax1.imshow(np_mask)
#     ax1.axis('off')
#     ax1.set_title("Mask", fontdict={'fontsize': 10, 'fontweight': 'medium', 'fontname': 'Arial'})

#     plt2 = ax2.imshow(np_ann)
#     ax2.axis('off')
#     ax2.set_title("Duct Annotation", fontdict={'fontsize': 10, 'fontweight': 'medium', 'fontname': 'Arial'})
    
#     mask_patches =[mpatches.Patch(color=np.array(map_type_color[type]) / 255, label=type) for type in map_type_color]
#     mask_legend = font_manager.FontProperties(family="Arial", size=5)
#     ax1.legend(handles=mask_patches, loc=10, bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize="x-small", prop=mask_legend, frameon=False)

#     ann_patches =[mpatches.Patch(color=np.array(ann_color[type]) / 255, label=type) for type in ann_color]
#     ann_legend = font_manager.FontProperties(family="Arial", size=5)
#     ax2.legend(handles=ann_patches, loc=10, bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize="x-small", prop=ann_legend, frameon=False)
#     output_filename = core_name + "_mask_and_ann.pdf"
#     plt.suptitle(core_name, fontsize=13, y=0.8)
#     fig.savefig(os.path.join(output_dir, output_filename), bbox_inches='tight', dpi=300)
#     plt.close()


""" Script arguments """
parser = argparse.ArgumentParser(description="Create a plot of the mask")
parser.add_argument("input_fn", type=str, help="Input mask file")
parser.add_argument("type_asn_csv", type=str, help="CSV file that contains cell type assignments")
parser.add_argument("pair_csv", type=str, help="CSV file that contains the pairs of cell ids in masks to their actual cell ids")
parser.add_argument("output_dir", type=str, help="Output colour mask filename")
parser.add_argument("--make_dir", "-md", default=False, action="store_true")
parser.add_argument("--colour_mask", "-cm", default=False, action="store_true")
parser.add_argument("--colour_plot", "-cp", default=False, action="store_true")
parser.add_argument("--one_hot_enc_mask", "-hm", default=False, action="store_true")
parser.add_argument("--prob", "-pm", default=False, action="store_true")
# parser.add_argument("--mask_and_annotation", "-ma", default=False, action="store_true")


args = parser.parse_args()

input_fn = args.input_fn
type_asn_csv = args.type_asn_csv
pair_csv = args.pair_csv
output_dir = args.output_dir
is_make_dir = args.make_dir
is_colour_mask = args.colour_mask
is_colour_plot = args.colour_plot
is_one_hot_enc_mask = args.one_hot_enc_mask
is_prob = args.prob
# is_mask_and_annotation = args.mask_and_annotation

""" Error Checking """
# Check if the input file exists
if not os.path.exists(input_fn):
    raise Exception("No such file: {}".format(input_fn))

# Check if the file is a tiff format 
ft = imghdr.what(input_fn)
if not ft == "tiff":
    raise Exception("Expected file type tiff but got {}".format(ft))

# Check if the output directory exists 
if not os.path.isdir(output_dir):
    raise Exception("Output directory does not exist")


""" Set up constants """
df_asn = pd.read_csv(type_asn_csv, index_col=0)
df_pairs = pd.read_csv(pair_csv)
fn = os.path.basename(input_fn)
core_name = extract_sample_name_from_filename(df_asn, fn)

""" Creating output files """
# Create a directory to save all outputs if it does not exist in the output directory 
if is_make_dir: 
    output_dir_core = os.path.join(output_dir, core_name)
    if not os.path.isdir(output_dir_core):
        os.mkdir(output_dir_core)
else: 
    output_dir_core = output_dir

""" Terminate early if the file exists """
if is_colour_mask + is_one_hot_enc_mask + is_colour_plot == 1: 
    if is_colour_mask:
        fn = os.path.join(output_dir_core, core_name + "_colour_mask.tiff")
        if os.path.exists(fn):
            # print("----------", fn, "----------")
            exit()
        else: 
            print("++++++++++", fn, "++++++++++")
    if is_one_hot_enc_mask:
        fn = os.path.join(output_dir_core, core_name + "_one_hot_encoding_mask.tiff")
        if os.path.exists(fn):
            # print("----------", fn, "----------")
            exit()
        else:
            print("++++++++++", fn, "++++++++++")


""" Create data """
mask = tf.imread(input_fn)
df_mask = pd.DataFrame(mask)

# Convert cell ids in mask file to the matching actual cell ids 
df_mask = map_mask_to_cell_ids(df_pairs, df_mask, core_name)


# if not is_prob: 
df_sample_asn = get_cell_type_assignments(df_asn, core_name)

# Convert mask containing cell ids to cell type names 
df_type_name_mask = map_mask_id_to_cell_types(df_mask, df_sample_asn)

# Numpy array of colour mask 
np_rgba_mask = map_mask_cell_names_to_colors(df_type_name_mask)


""" 1. Get colour mask file """
if is_colour_mask: 
    # Save mask 
    save_mask_as_tiff(np_rgba_mask, os.path.join(output_dir_core, core_name + "_colour_mask.tiff"))


""" 2. Create one hot encoding tiff file where each channel is the cell type """
if is_one_hot_enc_mask:
    np_ch_mask = map_mask_cell_names_to_one_hot_encoding(df_type_name_mask)
    if not (test_hot_enc_correct(np_ch_mask, df_type_name_mask.to_numpy())):
        raise Exception("One hot encoding mask is incorrectly defined for mask {}".format(tiff_file))
    np_ch_mask = np_ch_mask.transpose(1, 2, 0).astype(np.uint8)
    filename = os.path.join(output_dir_core, core_name + ".tiff")

    save_mask_as_tiff(np_ch_mask, os.path.join(output_dir_core, core_name + "_one_hot_encoding_mask.tiff"))


""" 3. Create plot """
if is_colour_plot:
    plot_mask(np_rgba_mask, core_name, output_dir_core)

""" 4. Create probabilities plot """ 
if is_prob: 
    print(core_name)
    # Get probabilities 
    df_max_prob = get_cell_type_probabilities(df_asn, core_name)
    # print(df_max_prob)
    # # 
    df_prob_mask = map_mask_id_to_probabilities(df_mask, df_max_prob)

    np_prob_mask = df_prob_mask.to_numpy()
    min_prob = np.min(np_prob_mask[np.nonzero(np_prob_mask)])
    a = np.ma.masked_where(np_prob_mask < 0.3927083311428659, np_prob_mask)
    # import copy
    # cmap = copy.copy(plt.cm.get_cmap("viridis"))
    cmap = plt.cm.viridis
    cmap.set_bad(color='black')

    print(min_prob)
    plt.imsave(os.path.join(output_dir_core, core_name + ".png"), a, cmap="viridis", vmin=0.3927083311428659, vmax=1)
