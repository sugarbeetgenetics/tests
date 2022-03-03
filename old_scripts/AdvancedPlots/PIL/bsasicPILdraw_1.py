from PIL import Image, ImageDraw
import numpy as np
import pandas as pd
import math
import argparse
from matplotlib import colors as cls

parser = argparse.ArgumentParser(description='Draw lines with markings')
parser.add_argument('-l', '--lines', type=str, metavar='', required=True, help='lines to draw')
parser.add_argument('-m','--multiply', type=int, metavar='', required=False, help='Set size multiplyer (Default : 1)', default="1")
parser.add_argument('-t','--thickness', type=int, metavar='', required=False, help='Set size mark thickness (Default : 2)', default="2")
parser.add_argument('-lc','--line_color', type=str, metavar='', required=False, help='Set line color (Default : blue) choices: red, black, coral, tomator, lime, darkgreen, cyan', default="blue")
parser.add_argument('-mc','--mark_color', type=str, metavar='', required=False, help='Set line color (Default : red) choices: red, black, coral, tomator, lime, darkgreen, cyan', default="red")
args = parser.parse_args()
linesFile =  args.lines
multiply_size =  args.multiply
thickness_mark = args.thickness
line_items_df = pd.read_csv(str(linesFile), sep="\t", comment="#", low_memory=False, header=0)
line_items_df_np = np.array(line_items_df)
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

def straight_line(data, start, end, yaxis, color, thickness):
    data[yaxis-thickness:yaxis+thickness, start:end] = color

def markings(data, position, yaxis, color, thickness):
    data[yaxis-thickness:yaxis+thickness, position-thickness:position+thickness] = color

# define canvas size
thickness =  min(4, 4*multiply_size)
padding_h = 50 * multiply_size
effective_space_one_line = 8  * multiply_size
gap = 10  * multiply_size
padding_w = 25  * multiply_size
if (line_items_df_np[:,0].max(axis=0)>4000):
    height = (2*padding_h + (line_items_df_np.shape[0])*effective_space_one_line + (line_items_df_np.shape[0])*gap)
    width = (2*padding_w + (line_items_df_np[:,0].max(axis=0)))
else:
    height = 2164
    width = 4096

# draw canvas
canvas = np.zeros((height, width, 3), dtype=np.uint8)
#create background white
canvas.fill(255)


img = Image.fromarray(canvas, 'RGB')
img.save('my.png')

# open image to modify
im = Image.open("my.png")
draw = ImageDraw.Draw(im)
draw.line((0, 0) + im.size, fill=128)
draw.line((0, im.size[1], im.size[0], 0), fill=128)

# write to stdout
im.save("new_my.png")


'''
# split marking column
i=0
yaxis = 0
for line,mark in line_items_df_np:
    mark_split = mark.split(sep=",")
    #print (mark_split[2])
    #print (line)
    line2 =  int(math.floor(math.log10(line))) * 2000
    if(i==0):
        yaxis = padding_h + yaxis
    else:
        yaxis = gap + yaxis + effective_space_one_line
    
    line_start = padding_w
    line_end = padding_w + line2

    #black_canvas[yaxis-8:yaxis+8, line_start:line_end] = [255,0,0]
    straight_line(canvas, line_start, line_end, yaxis, line_color_name, thickness)
    for item in mark_split:
        #print (item)
        item = int(math.floor(math.log10(float(item) + float(line)))) 
        markings(canvas, padding_w+math.floor(float(item)), yaxis, mark_color_name, thickness)
    #print (yaxis)
    i=i+1

#print (width, hight)

#white_canvas[10,10] = []

#img.show()

print (line_color_name)
print (mark_color_name)
print (multiply_size)
'''