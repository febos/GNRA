import os,requests,json,glob


date   = "2025-01-13"
folder = "./PDB"

def DownloadPDB():

    empty = False
    start = 0
    rows  = 25

    while not empty:

        empty = True
        print(start)

        query = "%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A"+\
                "%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22"+\
                "attribute%22%3A%22entity_poly.rcsb_entity_polymer_type%22%2C%22operator%22%3A%22exact_match%22%2C%22negation%22%3Afalse"+\
                "%2C%22value%22%3A%22DNA%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%"+\
                "22attribute%22%3A%22entity_poly.rcsb_entity_polymer_type%22%2C%22operator%22%3A%22exact_match%22%2C%22negation%22%3Afalse"+\
                "%2C%22value%22%3A%22RNA%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22"+\
                "attribute%22%3A%22entity_poly.rcsb_entity_polymer_type%22%2C%22operator%22%3A%22exact_match%22%2C%22negation%22%3Afalse%2C%22"+\
                "value%22%3A%22NA-hybrid%22%7D%7D%5D%2C%22logical_operator%22%3A%22or%22%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22"+\
                "%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%2C%22operator%22%3A%22"+\
                "less_or_equal%22%2C%22negation%22%3Afalse%2C%22value%22%3A%22{}%22%7D%7D%5D%2C%22label%22%3A%22text%22%7D%2C%22".format(date)+\
                "return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A{}%2C%22rows%22%3A{}%7D%2C%22".format(start,
                                                                                                                                                   rows)+\
                "results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22"+\
                "desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%7D"

        x = requests.get("https://search.rcsb.org/rcsbsearch/v2/query?json="+query)

        if "result_set" in x.json():
            
            for row in x.json()["result_set"]:
                print("downloading {}...".format(row['identifier']),end=' ')

                #skip if exists
                if os.path.exists(os.path.join(folder,"{}.cif").format(row['identifier']))\
                   and os.path.getsize(os.path.join(folder,"{}.cif".format(row['identifier']))) > 0:
                    print('SKIPPED')
                    empty = False
                    continue
                    
                r = requests.get("https://files.rcsb.org/download/{}.cif.gz".format(row['identifier']))
                filename = os.path.join(folder,"{}.cif.gz".format(row['identifier']))
                open(filename, 'wb').write(r.content)
                os.system('gunzip ' + filename)
                print('DONE')
                empty = False

            start += rows
        


if __name__ == "__main__":

    os.makedirs(folder, exist_ok=True)
    DownloadPDB()
