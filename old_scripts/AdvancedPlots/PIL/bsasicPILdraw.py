# version 1
# can draw lines and marking on them. input file must be a table with two columns, where first column must tell the line length and next column must tell the marking positions (comma separated)
# its limitation is when the data is too big to be seen in normal size image.
# next version will make lines in continuation
#usage python3 /media/pflanz/Hs1-2/Scripts/AdvancedPlots/PIL/bsasicPILdraw.py -l to_be_drawn.txt -lc darkgreen -t 10 -l_d 4 -m 1 -l_c 0.001

from PIL import Image
import numpy as np
import pandas as pd
import math
import argparse
from matplotlib import colors as cls
import os
from tqdm import tqdm

version="v1"
cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Draw lines with markings')
parser.add_argument('-l', '--lines', type=str, metavar='', required=True, help='lines to draw')
parser.add_argument('-m','--multiply', type=int, metavar='', required=False, help='Set size multiplyer (Default : 1)', default=1)
parser.add_argument('-l_c','--length_correct', type=float, metavar='', required=False, help='Set line length correction (Default : 1)', default=1)
parser.add_argument('-t','--thickness', type=int, metavar='', required=False, help='Set size mark thickness (Default : 2)', default=1)
parser.add_argument('-l_d','--line_dist', type=int, metavar='', required=False, help='Set distance between two lines (Default : 20)', default=20)
parser.add_argument('-lc','--line_color', type=str, metavar='', required=False, help='Set line color (Default : blue) choices: red, black, coral, tomator, lime, darkgreen, cyan', default="cyan")
parser.add_argument('-mc','--mark_color', type=str, metavar='', required=False, help='Set marking color (Default : red) choices: red, black, coral, tomator, lime, darkgreen, cyan', default="red")
args = parser.parse_args()
linesFile =  args.lines
multiply_size =  args.multiply
length_correction = args.length_correct
thickness_mark = args.thickness
line_distance = args.line_dist
line_items_df = pd.read_csv(str(linesFile), sep="\t", comment="#", low_memory=False, names=["length","marker"])
line_items_df.marker = line_items_df.marker.astype(str)
line_items_df_np = np.array(line_items_df)

_, idx = np.unique(line_items_df_np[:,0], return_index=True)
col1_uniq = line_items_df_np[np.sort(idx),0]
new_col2 = []
#print (col1_uniq)
for col1 in col1_uniq[:]:
    #print (col1)
    #print (line_items_df_np[line_items_df_np[:,0]==col1,1])
    joined_col2 = ",".join(line_items_df_np[line_items_df_np[:,0]==col1,1])
    #print (joined_col2)
    new_col2.append(str(joined_col2))

combined_data = pd.DataFrame(list(zip(col1_uniq, new_col2)), columns=['length', 'marker'])

combined_data_np = np.array(combined_data)
with open(str(cwd + '/' + version + "combined_data_plotted" + '.out'), 'w') as out_file:
    combined_data.to_csv(out_file, sep="\t", index=False)


#color name to decimal conversion
def nameTOdecimal(name):
    if (name=="red"):
        return (255,0,0)
    elif (name == "black"):
        return (0,0,0)
    elif (name == "coral"):
        return (255,127,80)
    elif (name == "tomato"):
        return (255,99,71)
    elif (name == "lime"):
        return (0,255,0)
    elif (name == "darkgreen"):
        return (0,100,0)
    elif (name == "cyan"):
        return (0,255,255)

line_color_name = nameTOdecimal(str(args.line_color))
mark_color_name = nameTOdecimal(str(args.mark_color))

#mark_items_df = pd.read_csv(str(markFile), sep="\t", comment="#", low_memory=False, header=0)
#mark_items_df_np = np.array(mark_items_df).squeeze()

def straight_line(data, start, end, yaxis, thickness, color="blue"):
    data[yaxis-thickness:yaxis+thickness, start:end] = color

def markings(data, position, yaxis, thickness, color="red"):
    data[yaxis-thickness-2:yaxis+thickness+2, position-thickness:position+thickness] = color

# define canvas size
thickness =  min(1, 4*multiply_size)
padding_h = 100
effective_space_one_line = 8
gap = 20 * line_distance
padding_w = 50

'''
height = (2*padding_h + (combined_data_np.shape[0])*effective_space_one_line + (combined_data_np.shape[0])*gap)
width = (2*padding_w + (combined_data_np[:,0].max(axis=0)))
'''

height = 3000
width = 10000

# draw canvas
canvas = np.zeros((height, width, 3), dtype=np.uint8)
#create background white
print ("height: " + str(height))
print ("width: " + str(width))
canvas.fill(255)

# split marking column
i=0
yaxis = 0
for line,mark in tqdm(combined_data_np):
    mark_split = mark.split(sep=",")
    #line = math.floor(float(line))
    #print (type(mark_split[2]))
    #print (type(line))
    line = math.floor(line * length_correction)
    if(i==0):
        yaxis = padding_h + yaxis
    else:
        yaxis = gap + yaxis + effective_space_one_line
    
    line_start = padding_w
    line_end = padding_w + line
    #print ( line_color_name)
    #black_canvas[yaxis-8:yaxis+8, line_start:line_end] = [255,0,0]
    straight_line(canvas, line_start, line_end, yaxis, thickness, line_color_name)
    for item in mark_split:
        item = float(item) * length_correction
        #print (item)
        markings(canvas, padding_w+math.floor(float(item)), yaxis, 2*thickness, mark_color_name)
    #print (yaxis)
    i=i+1

#print (width, hight)

#white_canvas[10,10] = []

img = Image.fromarray(canvas, 'RGB')
img.save('my.png')
#img.show()

#print (line_color_name)
#print (mark_color_name)
#print (multiply_size)
