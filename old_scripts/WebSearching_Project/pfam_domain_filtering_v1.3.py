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

parser = argparse.ArgumentParser(description='Search for pfam domains online and assign them into a category!')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='pfam domain list')
parser.add_argument('-o','--out_file', type=str, metavar='', required=True, help='out file name')

args = parser.parse_args()
version = "v1.1"
cwd = os.getcwd()
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


ids_points = np.empty((0,6)) #in order to save ids and corresponding points

#print (list_items_df_np.shape[0])

for i in tqdm(range(0,list_items_df_np[:,2].shape[0],1)):
	ids = list_items_df_np[i,2]
	#print (ids)
	site = "http://pfam.xfam.org/family/" + str(ids)
	raw_html = simple_get(str(site))
	html = BeautifulSoup(raw_html,'lxml')
	#print (site)
	tags = html.find_all('div', attrs={'class':'pfamData'})
	tagsWP = html.find_all('div', attrs={'class':'wpData'})
	#print (tagsPfam)
	
	#tags.append(tagsWP)
	
	for div in tagsWP:
		full_text = ""
		text = div.get_text()
		full_text = full_text + " BLOCK END START " +  text
		#print (full_text)
		full_text = text
		score = ""
		points = 0
		if len(re.findall('[Dd]efense\D',full_text)) > 0:
			points = points + 2
			score = score + "Defense, "
		elif len(re.findall('[Rr]esistance\D',full_text)) > 0:
			points = points + 2
			score = score + "Resistance, "
		elif len(re.findall('[Aa]poptosis\D',full_text)) > 0:
			points = points + 1
			score = score + "Apoptosis, "
		elif len(re.findall('[Cc]ell death\D',full_text)) > 0:
			points = points + 1
			score = score + "Cell death, "
		elif len(re.findall('[Nn]ematode\D',full_text)) > 0:
			points = points + 2
			score = score + "Nematode, "
		elif len(re.findall('[Cc]yst nematode\D',full_text)) > 0:
			points = points + 3
			score = score + "Cyst nematode, "
		elif len(re.findall('[Ss]ignal\D',full_text)) > 0:
			points = points + 1
			score = score + "Signal, "
		elif len(re.findall('WRKY\D',full_text)) > 0:
			points = points + 2
			score = score + "WRKY, "        
		elif len(re.findall('[Kk]inase\D',full_text)) > 0:
			points = points + 1
			score = score + "Kinase, "
		elif len(re.findall('[Pp]rotease\D',full_text)) > 0:
			points = points + 1
			score = score + "Protease, "
		elif len(re.findall('[Ee]mbryo\D',full_text)) > 0:
			points = points + 1
			score = score + "Embryo, "
		elif len(re.findall('[Ll]ethal\D[Ee]mbryo',full_text)) > 0:
			points = points + 1
			score = score + "Lethal, "
		elif len(re.findall('EMB\W',full_text)) > 0:
			points = points + 10
			score = score + "EMB, "
		elif len(re.findall('[Rr]oot\D',full_text)) > 0:
			points = points + 2
			score = score + "Root, "
		elif len(re.findall('[Hh]istone\D',full_text)) > 0:
			points = points + 2
			score = score + "Histone modification, "
		elif len(re.findall('[Hh]omozygous recessive\D',full_text)) > 0:
			points = points + 2
			score = score + "Homozygous recessive, "
		elif len(re.findall('[Hh]omozygous\D',full_text)) > 0:
			points = points + 2
			score = score + "Homozygous, "
		elif len(re.findall('[Rr]ecessive\D',full_text)) > 0:
			points = points + 2
			score = score + "Recessive, "
		elif len(re.findall('[Dd]ie\D',full_text)) > 0:
			points = points + 2
			score = score + "Die, "
		elif score == "":
			score = "None"
		ids_points = np.vstack((ids_points,np.array([str(site),list_items_df_np[i,0],list_items_df_np[i,1] ,str(ids), points, score])))

ids_points_pd = pd.DataFrame(ids_points)

name_outfile = "root_search"

with open(out_file_name, 'w') as out_file:
    ids_points_pd.to_csv(out_file, sep="\t", index=False)

