# usage
# python3 pfam_domain_filtering.py -l /media/pflanz/Avneesh/Beet_translocation/ProperDataSet/TRregion520/Tr520region_withedges.ORFfinder.PfamScan.bed.pfamIDs


from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re

version="v1_3"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Search for InterPro domains online and assign them into a category!')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='InterPro domain list')

args = parser.parse_args()

list_items = args.list
list_items_df = pd.read_table(list_items, delimiter="\t", comment="#", names=["locus","annotation","IPR_id"], low_memory=False)
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


ids_points = np.empty((0,6)) #in order to save ids and corresponding points

#print (list_items_df_np.shape[0])

for i in tqdm(range(0,list_items_df_np[:,2].shape[0],1)):
	ids = list_items_df_np[i,2]
	#print (ids)
	site = "https://www.ebi.ac.uk/interpro/entry/" + str(ids)
	raw_html = simple_get(str(site))
	html = BeautifulSoup(raw_html,'lxml')
	#print (site)
	tags = html.find_all('div', attrs={'class':'goSections'})
	tagsWP = html.find_all('div', attrs={'class':'wpData'})
	#print (tagsPfam)
	#tags.append(tagsWP)
	
	for div in tagsWP:
		full_text = ""
		text = div.get_text()
		#print (full_text)
		full_text = text
		score = ""
		points = 0
		print (div)
		print (text)
