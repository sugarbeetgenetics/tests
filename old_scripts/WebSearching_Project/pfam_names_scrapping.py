
from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re

version="v1_3"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Search for pfam domains online and assign them into a category!')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='pfam domain list')
parser.add_argument('-o','--out_file', type=str, metavar='', required=True, help='out file name')

args = parser.parse_args()

list_items = args.list
list_items_df = pd.read_table(list_items, delimiter="\t", comment="#", names=["locus","annotation","pfam_id"], low_memory=False)
list_items_df_np = np.array(list_items_df).squeeze()

version = "v1.1"
cwd = os.getcwd()
out_file_name = cwd + "/" + version + args.out_file + str(".txt")

#html parse
from requests import get
from requests.exceptions import RequestException
from contextlib import closing
from bs4 import BeautifulSoup

def simple_get(url):
    """
    Attempts to get the content at `url` by making an HTTP GET request.
    If the content-type of response is some kind of HTML/XML, return the
    text content, otherwise return None.
    """
    try:
        with closing(get(url, stream=True)) as response:
            if is_good_response(response):
                return response.content
            else:
                return None

    except RequestException as e:
        log_error('Error during requests to {0} : {1}'.format(url, str(e)))
        return None


def is_good_response(response):
    """
    Returns True if the response seems to be HTML, False otherwise.
    """
    content_type = response.headers['Content-Type'].lower()
    return (response.status_code == 200 
            and content_type is not None 
            and content_type.find('html') > -1)


def log_error(e):
    """
    It is always a good idea to log errors. 
    This function just prints them, but you can
    make it do anything.
    """
    print(e)


ids_points_all = [] #in order to save ids and corresponding points

#print (list_items_df_np.shape[0])

for i in tqdm(range(0,list_items_df_np[:,2].shape[0],1)):
    ids = list_items_df_np[i,2]
    #print (ids)
    site = "http://pfam.xfam.org/family/" + str(ids)
    raw_html = simple_get(str(site))
    html = BeautifulSoup(raw_html,'lxml')
	#print (site)
    tags = html.find_all('div', class_ = 'pfamData')
    #print ("id: " + str(ids) + " " + str(tags.h1.text))

    ids_points = ""
    for div in tags:
        text = div.h1.text
        ids_points = str(site) + "\t" + list_items_df_np[i,0] + "\t" + list_items_df_np[i,1] + "\t" + str(text)
    ids_points_all.append(ids_points)

ids_points_pd = pd.DataFrame(ids_points_all)

with open(out_file_name, 'w') as out_file:
    ids_points_pd.to_csv(out_file, sep="\t", index=False)
