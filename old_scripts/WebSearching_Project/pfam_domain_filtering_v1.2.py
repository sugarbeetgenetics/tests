# usage
# python3 pfam_domain_filtering.py -l /media/pflanz/Avneesh/Beet_translocation/ProperDataSet/TRregion520/Tr520region_withedges.ORFfinder.PfamScan.bed.pfamIDs


from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re

version="v1_2"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Search for pfam domains online and assign them into a category!')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='pfam domain list')

args = parser.parse_args()

list_items = args.list
list_items_df = pd.read_table(list_items, delimiter="\t", comment="#", names=["locus","annotation","pfam_id"], low_memory=False)
list_items_df_np = np.array(list_items_df).squeeze()



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


ids_points = np.empty((0,5)) #in order to save ids and corresponding points

#print (list_items_df_np.shape[0])

for i in tqdm(range(0,list_items_df_np[:,1].shape[0],1)):
    ids = list_items_df_np[i,1]
    #print (ids)
    site = "http://pfam.xfam.org/family/" + str(ids)
    raw_html = simple_get(str(site))
    html = BeautifulSoup(raw_html,'lxml')
    #print (site)
    tags = html.find_all('div', attrs={'class':'pfamData'})
    for div in tags:
        full_text = ""
        text = div.get_text()
        full_text = full_text + " BLOCK END START " +  text
        #print (full_text)
    score = ""
    points = 0
    if len(re.findall('[Dd]efense',full_text)) > 0:
        points = points + 2
        score = score + "Defense, "
    if len(re.findall('[Rr]esistance',full_text)) > 0:
        points = points + 2
        score = score + "Resistance, "
    if len(re.findall('[Aa]poptosis',full_text)) > 0:
        points = points + 1
        score = score + "Apoptosis, "
    if len(re.findall('[Cc]ell death',full_text)) > 0:
        points = points + 1
        score = score + "Cell death, "
    if len(re.findall('[Nn]ematode',full_text)) > 0:
        points = points + 2
        score = score + "Nematode, "
    if len(re.findall('[Cc]yst nematode',full_text)) > 0:
        points = points + 3
        score = score + "Cyst nematode, "
    if len(re.findall('[Ss]ignal',full_text)) > 0:
        points = points + 1
        score = score + "Signal, "
    if len(re.findall('WRKY',full_text)) > 0:
        points = points + 2
        score = score + "WRKY, "        
    if len(re.findall('[Kk]inase',full_text)) > 0:
        points = points + 1
        score = score + "Kinase, "
    if len(re.findall('[Pp]rotease',full_text)) > 0:
        points = points + 1
        score = score + "Protease, "
    if len(re.findall('[Ll]ethal',full_text)) > 0:
        points = points + 1
        score = score + "Lethal, "
    if score == "":
        score = "None"

    ids_points = np.vstack((ids_points,np.array([list_items_df_np[i,0],list_items_df_np[i,1] ,str(ids), points, score])))

ids_points_pd = pd.DataFrame(ids_points)

with open(str(cwd + '_' + 'ids_points.out'), 'w') as out_file:
    ids_points_pd.to_csv(out_file, sep="\t", index=False)
