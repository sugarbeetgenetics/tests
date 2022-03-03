######## improvements #########
# Draw gene structure 

# usage 
#/media/pflanz/Avneesh/Beet_translocation/ProperDataSet/repeats$ python3 /media/pflanz/Hs1-2/Scripts/AdvancedPlots/PIL/bsasicPILdraw_V2.py -l to_be_drawn.txt -m 3 -t 2 -lc blue -mc red -o test -l_c 0.01 -l_s 20

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
import numpy as np
import pandas as pd
import math
import argparse
from matplotlib import colors as cls
import os
import sys
from tqdm import tqdm
import color_nameTOdecimal as nTd

version="V2.1"
cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Draw lines with markings')
parser.add_argument('-l', '--lines', type=str, metavar='', required=True, help='lines to draw')
parser.add_argument('-o', '--out_file', type=str, metavar='', required=True, help='out plot file name')
parser.add_argument('-m', '--multiply', type=int, metavar='', required=False, help='Set size multiplyer (Default : 1)', default=1)
parser.add_argument('-l_c', '--length_correct', type=float, metavar='', required=False, help='Set line length correction (Default : 1)', default=1)
parser.add_argument('-t', '--thickness', type=int, metavar='', required=False, help='Set size mark thickness (Default : 1)', default=1)
parser.add_argument('-l_d', '--line_distance', type=int, metavar='', required=False, help='Set distance between two lines (Default : 20)', default=20)
parser.add_argument('-l_s', '--label_size', type=int, metavar='', required=False, help='Set label size (Default : 16)', default=16)
parser.add_argument('-i_h', '--image_height', type=int, metavar='', required=False, help='Set image height (Default : 3000)', default=3000)
parser.add_argument('-i_w', '--image_width', type=int, metavar='', required=False, help='Set image width (Default : 4000)', default=4000)
#parser.add_argument('-order', '--order_file', type=str, metavar='', required=False, help='Set order', default='None')
parser.add_argument('-lc', '--line_color', type=str, metavar='', required=False, choices=["red", "blue", "black", "coral", "tomato", "lime", "darkgreen", "cyan"], help='Set line color (Default : blue) choices: red, blue, black, coral, tomato, lime, darkgreen, cyan', default="cyan")
parser.add_argument('-mc', '--mark_color', type=str, metavar='', required=False, choices=["red", "blue", "black", "coral", "tomato", "lime", "darkgreen", "cyan"], help='Set marking color (Default : red) choices: red, blue, black, coral, tomato, lime, darkgreen, cyan', default="red")

if len(sys.argv)<1:
    #parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

linesFile =  args.lines
multiply_size =  args.multiply
length_correction = args.length_correct
thickness_mark = args.thickness
line_distance = args.line_distance
label_size_defined = args.label_size
defined_height = args.image_height
defined_width = args.image_width
#order_defined = args.order_file
#order_defined_lst = pd.read_csv(str(order_defined), sep="\t", comment="#", low_memory=False, names=["name_order"])
out_file_name = cwd + "/" + version + args.out_file + str(".png")

line_items_df = pd.read_csv(str(linesFile), sep=" ", comment="#", low_memory=False, names=["chr","prog","annotation","start","stop","score","strand","dot","info"])
#line_items_df["start"] = line_items_df["start"].astype(int) #make marker column as numeric
#line_items_df.stop = line_items_df.stop.astype(int) #make marker column as numeric
#line_items_df.chr = line_items_df.chr.astype(str) #make marker column as string
#line_items_df.annotation = line_items_df.annotation.astype(str) #make marker column as string

#print (order_defined_lst)
#order_list = list(order_defined_lst['name_order'])
#print (order_list)
#if (args.order_file != 'None'):
    # reorder data based on given order
#    line_items_df.set_index('name',drop=False, inplace=True)
#    line_items_df.reindex(order_list)

line_items_df_np = np.array(line_items_df)

############### fix data if it is distributed in several lines for each scaffold ###############
"""
_, idx = np.unique(line_items_df_np[:,0], return_index=True)
col0_uniq = line_items_df_np[np.sort(idx),0] # uniq scaffolds ids
col1_uniq = line_items_df_np[np.sort(idx),1] # line length for each uniq scaffold id

# combine points for each scaffold into one line
new_col2 = []
#print (col1_uniq)
for col1 in col1_uniq[:]:
    #print (col1)
    #print (line_items_df_np[line_items_df_np[:,0]==col1,1])
    joined_col2 = ",".join(line_items_df_np[line_items_df_np[:,1]==col1,2])
    #print (joined_col2)
    new_col2.append(str(joined_col2))

combined_data = pd.DataFrame(list(zip(col0_uniq, col1_uniq, new_col2)), columns=['ids', 'length', 'marker'])

combined_data_np = np.array(combined_data)
with open(str(cwd + '/' + version + "combined_data_plotted" + '.out'), 'w') as out_file:
    combined_data.to_csv(out_file, sep="\t", index=False)
################################################################################################
"""
line_color_name = nameTOdecimal(str(args.line_color))
mark_color_name = nameTOdecimal(str(args.mark_color))

#mark_items_df = pd.read_csv(str(markFile), sep="\t", comment="#", low_memory=False, header=0)
#mark_items_df_np = np.array(mark_items_df).squeeze()
def write_text(image, text, yaxis, xaxis=10, text_size=16):
    font = ImageFont.truetype("arial.ttf", text_size)
    image.text((xaxis, yaxis), text, (0,0,0), font=font)

def straight_line(data, start, end, yaxis, thickness, color="blue"):
    data[yaxis-thickness:yaxis+thickness, start:end] = color

def start_marking(data, start, yaxis):
    data[yaxis-(thickness+6):yaxis+(thickness+6), start] = [0,0,0]

def end_marking(data, end, yaxis):
    data[yaxis-(thickness+6):yaxis+(thickness+6), end] = [0,0,0]

def markings(data, position, yaxis, thickness=1, color="red"):
    data[yaxis-(thickness+2):yaxis+(thickness+2), position] = color

def nested_change(input, func): #define function to change list data type, replace func with desired dtype
    if isinstance(input, list):
        return [nested_change(x, func) for x in input]
    return func(input)

#### define canvas parameters #######
thickness =  thickness_mark
padding_h = 100
effective_space_one_line = 8 + thickness + line_distance
padding_w_l = 250
padding_w_r = 10



'''
#height = (2*padding_h + (combined_data_np.shape[0])*effective_space_one_line + (combined_data_np.shape[0])*gap)
#width = (2*padding_w + (combined_data_np[:,0].max(axis=0)))
'''

height = defined_height
width = defined_width
###############################

######################## draw on canvas #############################
 
canvas = np.zeros((height, width, 3), dtype=np.uint8)

print ("canvas height: " + str(height))
print ("canvas width: " + str(width))
canvas.fill(255) # paint canvas white

# split marking column
label_up = [] # label upper limit of yaxis
label_bottom = [] # label lower limit yaxis
id_name_all = []
i=0
yaxis = 0
previous_line_end = 0
gene_start = 0
gene_end = 0
gene_length = 0
per_base_pixels = 0
for i in range(0,line_items_df_np.shape[0],1):
    print(line_items_df_np[i,0])
    chr = line_items_df_np[i,0]
    start = line_items_df_np[i,4]
    stop = line_items_df_np[i,5]
    annotation = line_items_df_np[i,3]
    strand = line_items_df_np[i,7]

    if(annotation=="gene"):
        yaxis = yaxis + 300
        line_start = padding_w_l
        line_end = width - padding_w_r
        straight_line(canvas, line_start, line_end, yaxis, thickness, line_color_name)
        per_base_pixels = math.floor((stop-start)/(line_end - line_start))
    print(line_items_df_np[i,0],stop, stop-start, per_base_pixels)
    markings(canvas,(gene_start-start)*per_base_pixels, yaxis, thickness, mark_color_name)
    

#print (width, hight)

#white_canvas[10,10] = []

img = Image.fromarray(canvas, 'RGB')
img.save(out_file_name)

# write text on created image
read_img = Image.open(out_file_name)
draw = ImageDraw.Draw(read_img)
for text_idx in range(0,len(label_up[:]),1):
    write_text(draw,str(id_name_all[text_idx]), label_up[text_idx], 2, label_size_defined)
    #print ("id " + str(id_name_all[text_idx]) + "----- yaxis " + str(label_up[text_idx]))
read_img.save(out_file_name)

