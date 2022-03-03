
import multidict as multidict

import numpy as np

import os
import re
from PIL import Image
from os import path
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='make tag-frequency plot')
parser.add_argument('-d','--data', type=str, metavar='', required=True, help='data to plot')
parser.add_argument('-o','--out_file', type=str, metavar='', required=True, help='out file name')
args = parser.parse_args()
inputFile =  args.data

version = "v1.1"
cwd = os.getcwd()
out_file_name = cwd + "/" + version + args.out_file + str(".png")

def getFrequencyDictForText(sentence):
    fullTermsDict = multidict.MultiDict()
    tmpDict = {}

    # making dict for counting frequencies
    for text in sentence.split("\n"):
        if re.match("a|the|an|the|to|in|for|of|or|by|with|is|on|that|be", text):
            continue
        val = tmpDict.get(text, 0)
        tmpDict[text.lower()] = val + 1
    for key in tmpDict:
        fullTermsDict.add(key, tmpDict[key])
    return fullTermsDict


def makeImage(text,out):
    beet_mask = np.array(Image.open("/media/pflanz/Hs1-2/Scripts/AdvancedPlots/sugarbeet.png"))

    wc = WordCloud(background_color="white", max_words=1000, mask=beet_mask)
    # generate word cloud
    wc.generate_from_frequencies(text)

    # show
    plt.imshow(wc, interpolation="bilinear")
    plt.axis("off")
    plt.savefig(out, dpi=440,)

text = open(inputFile, encoding='utf-8')
text = text.read()
makeImage(getFrequencyDictForText(text), out_file_name)