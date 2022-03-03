# usage
# python3 pfam_domain_filtering.py -l /media/pflanz/Avneesh/Beet_translocation/ProperDataSet/TRregion520/Tr520region_withedges.ORFfinder.PfamScan.bed.pfamIDs


from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re


parser = argparse.ArgumentParser(description='Search for SRR ids for respective R ids on ncbi!')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='IDs to be scanned')
parser.add_argument('-o','--out_file', type=str, metavar='', required=True, help='out file name')

args = parser.parse_args()

version="v1_3"
cwd = os.getcwd() + "/" + version
out_file_name = cwd + "/" + version + args.out_file + str(".txt")

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

id_file = pd.read_csv("/mnt/HDD-01/sarah/991accessions_mapping/DataFromChina/spring_minDepth9.R_ids.lst", sep="\t", comment="#", low_memory=False, header=0)


#print (list_items_df_np.shape[0])
for id in id_file.iloc[:,0]:
    print(id)
    print("\n")
    site = "https://www.ncbi.nlm.nih.gov/sra/?term="+id
    raw_html = simple_get(str(site))
    html = BeautifulSoup(raw_html,'html')
    # get rslt tags
    rslt_tags = html.find_all('div', attrs={'class':"rslt"})
    # process each tag to get title tag
    count = 0
    for ptag in rslt_tags:
        count += 1
        title = ptag.find('p', class_='title')
        
        # check if title tag matches id
        if (id == str(title.get_text())):
            # get rprtid tag to get SRX number
            rprtid = ptag.find('dl', class_='rprtid')
            link = "https://www.ncbi.nlm.nih.gov/sra/"+ str(rprtid.get_text()).split(" ")[2] + str("[accn]")
            print(str(rprtid.get_text()).split(" ")[2])
            print(link)
            second_raw_html = simple_get(str(link))
            second_html = BeautifulSoup(second_raw_html,'html')
            second_rprt_tags = second_html.find_all('div', attrs={'class':"rprt"})
            
            for a_tag in second_rprt_tags:
                SRR_id = a_tag.find('td', align='left')
                SRR = str(SRR_id.get_text())
                print (SRR)
    if count <= 1:
        rprt_tags = html.find_all('div', attrs={'class':"rprt"})
        for a_tag in rprt_tags:
            SRR_id = a_tag.find('td', align='left')
            SRR = str(SRR_id.get_text())
            print (SRR)
        
ids_points_pd = pd.DataFrame(ids_points)

name_outfile = "root_search"

with open(out_file_name, 'w') as out_file:
    ids_points_pd.to_csv(out_file, sep="\t", index=False)

