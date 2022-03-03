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

version="v2"
cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Draw lines with markings')
parser.add_argument('-l', '--lines', type=str, metavar='', required=True, help='lines to draw')
parser.add_argument('-o', '--out_file', type=str, metavar='', required=True, help='out plot file name')
parser.add_argument('-m', '--multiply', type=int, metavar='', required=False, help='Set size multiplyer (Default : 1)', default=1)
parser.add_argument('-l_c', '--length_correct', type=float, metavar='', required=False, help='Set line length correction (Default : 1)', default=1)
parser.add_argument('-t', '--thickness', type=int, metavar='', required=False, help='Set size mark thickness (Default : 1)', default=1)
parser.add_argument('-l_d', '--line_distance', type=int, metavar='', required=False, help='Set distance between two lines (Default : 20)', default=20)
parser.add_argument('-l_s', '--label_size', type=int, metavar='', required=False, help='Set label size (Default : 16)', default=16)
parser.add_argument('-lc', '--line_color', type=str, metavar='', required=False, choices=["red", "blue", "black", "coral", "tomato", "lime", "darkgreen", "cyan"], help='Set line color (Default : blue) choices: red, blue, black, coral, tomato, lime, darkgreen, cyan', default="cyan")
parser.add_argument('-mc', '--mark_color', type=str, metavar='', required=False, choices=["red", "blue", "black", "coral", "tomato", "lime", "darkgreen", "cyan"], help='Set marking color (Default : red) choices: red, blue, black, coral, tomato, lime, darkgreen, cyan', default="red")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

linesFile =  args.lines
multiply_size =  args.multiply
length_correction = args.length_correct
thickness_mark = args.thickness
line_distance = args.line_distance
label_size_defined = args.label_size
out_file_name = cwd + "/" + version + args.out_file + str(".png")

line_items_df = pd.read_csv(str(linesFile), sep="\t", comment="#", low_memory=False, names=["name","length","marker"])
line_items_df.marker = line_items_df.marker.astype(str)
line_items_df_np = np.array(line_items_df)

############### fix data if it is distibuted in several lines for each scaffold ###############
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
    elif (name == "blue"):
        return (0,0,255)

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
height = (2*padding_h + (combined_data_np.shape[0])*effective_space_one_line + (combined_data_np.shape[0])*gap)
width = (2*padding_w + (combined_data_np[:,0].max(axis=0)))
'''

height = 3000
width = 4000
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
for id_name,line,marking in tqdm(combined_data_np):
    id_name_all.append(id_name)
    #previous_line_end = 0
    mark_split_0 = marking.split(sep=",")
    #line = math.floor(float(line))
    mark_split = nested_change(mark_split_0, int)
    #print ("unsorted mark list: " + str(mark_split))
    mark_split.sort()
    #print ("sorted mark list: " + str(mark_split))
    line = math.floor(line * length_correction)
    #print ("Line size: " + str(line))
    if(i==0):
        yaxis = padding_h + yaxis
    else:
        yaxis = yaxis

    line_left = line
    available_width = (width-padding_w_l-padding_w_r)
    #print (available_width)
    label_up.append(yaxis)
    item = 0
    while_count = 0
    line_start = 0
    line_end = 0
    mark = 0
    while line_left >= available_width:
        line_start = padding_w_l
        if (line_left==line):

            start_marking_start = line_start
            start_marking(canvas, start_marking_start, yaxis)
               
        line_end = width - padding_w_r
        line_left = line_left - (line_end - line_start)
        straight_line(canvas, line_start, line_end, yaxis, thickness, line_color_name)
        #print ("Drawn length: " + str(line_end-line_start))
        #print("Line left: " + str(line_left))
        #print("Yaxis : " + str(yaxis))

        for item in range(item,len(mark_split)):
            mark = math.floor((float(mark_split[item]) * length_correction) + padding_w_l - while_count*(width-padding_w_l-padding_w_r))
            #print ("WHILE line 153 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
            if (mark <= (width - padding_w_r)) and (mark >= padding_w_l) :
                #print ("WHILE line 157 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
                markings(canvas, mark, yaxis, thickness, mark_color_name)
                item = item + 1
            else:
                break
        
        yaxis = yaxis + effective_space_one_line
        while_count = while_count + 1
    #print ("WHILE line 162 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
    item_count = item
    #print("before If " + str(item))
    if line_left < available_width:
        line_start = padding_w_l   
        if (while_count == 0):
            start_marking(canvas,line_start, yaxis)

        line_end = padding_w_l + line_left
        #print ( line_color_name)
        #black_canvas[yaxis-8:yaxis+8, line_start:line_end] = [255,0,0]
        straight_line(canvas, line_start, line_end, yaxis, thickness, line_color_name)
        end_marking_end = line_end
        end_marking(canvas, end_marking_end, yaxis)
        for item in range(0,len(mark_split)):
            
#            if (while_count > 0):
            mark = math.floor((float(mark_split[item]) * length_correction) + padding_w_l - while_count*(width-padding_w_l-padding_w_r))
            #print("If " + str(item) + " While count: " + str(while_count) + "mark: " + str(mark))
#            else:
#                mark = math.floor((float(mark_split[item]) * length_correction) + padding_w_l)
            #print ("IF line 180 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
            if (mark <= (width - padding_w_r)) and (mark >= padding_w_l) :
                #print("If marking " + str(item))
                #print ("IF line 184 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
                markings(canvas, mark, yaxis, thickness, mark_color_name)
                item = item + 1
            else:
                pass
        line_left = line_left - (line_end - line_start)
        label_bottom.append(yaxis)
        yaxis = yaxis +  effective_space_one_line + 2*line_distance
        #print("final " + str(item))        
        #print("Line left: " + str(line_left))
        #print("Yaxis : " + str(yaxis))
    #print ("IF line 191 - line start: " + str(line_start) + ", mark: " + str(mark) + ", while count: " + str(while_count) + ", Line left: " + str(line_left) + ", item count: " + str(item) + ", i: " + str(i))
    #print (yaxis)
    while_count = 0
    item = 0
    i=i+1

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