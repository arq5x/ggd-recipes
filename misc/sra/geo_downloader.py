#!/usr/bin/env python
"""Script for automated downloading of
raw data from GEO database
Example:
    python geo_downloader.py --dfile dataset.md [From Methbase-tracker]
    python geo_downloader.py
"""
import argparse
from Bio import Entrez
from bs4 import BeautifulSoup
import csv
from ftplib import FTP
import futures
import hashlib
import json
import math
import os
import sys
import tempfile
import yaml

from hurry.filesize import size

__GEO_LINK_NAME__ = 'GSE Link'
__NCBI_FTP__ = 'ftp-trace.ncbi.nlm.nih.gov'
__GEO_LINK_COLUMN__ = 3
"""
Checks to run on NCBI output for a Geo project
Example: Ensure the results were generated from a Methylation experiment
{'gdsType': 'Methylation profiling by high throughput sequencing'}"""
__DATA_CHECKS__ = {}
# Absolute path where files are downloaded
__ROOT_DOWNLOAD_LOCATION__ = None
__RETMAX__ = 10**9


def _set_root_download_path(path):
    global __ROOT_DOWNLOAD_LOCATION__
    __ROOT_DOWNLOAD_LOCATION__ = path
class GEOQuery:
    def __init__(self, search_term=None, email="all@smithlabresearch.org"):
        Entrez.email = email
        # Geo Datasets database
        self.database = 'gds'
        # Refer http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html
        # on how to search
        # For e.g. "Methylation profiling by high throughput sequencing"[DataSet Type] AND "Homo Sapiens"[Organism]
        # or "Methylation profiling by high throughput sequencing"[DataSet Type] AND "Mus Musculus"[Organism]
        self.search_term = search_term
        self.max_records = None

    def submit_query(self, retstart=0):
        """
        Submit Query to GEODatasets
        Returns records
        """
        query = Entrez.esearch(db=self.database, term=self.search_term, retmax=__RETMAX__, retstart=retstart)
        record = Entrez.read(query)
        if self.max_records is None:
            self.max_records = record['Count']
        return record

    def get_titles(self, uid_list=None):
        """
        Returns titles
        """
        query = Entrez.esummary(db=self.database, id=(",").join(uid_list))
        records = Entrez.read(query)
        return records

    def get_authors(self, pubmed_id_list):
        """Get authors form pubmed id list
        Data can exist without the publication having a pubmed id.
        So a more concrete way to store the data should be 'Contributor' field.
        """
        query = Entrez.esummary(db="pubmed", id=(",").join(pubmed_id_list))
        records = Entrez.read(query)
        return records

    def get_gsm_from_gse(self, geosample_id):
        # GDS = Geo DataSets
        #         ^
        #         |
        # GPL = Geo PLatform
        #         ^
        #         |
        # GSE = Geo Sample sEries
        #         ^
        #         |
        # GSM = Geo SaMple
        #         ^
        #         |
        #       Raw SRA

        # Get all GSM from each GSE
        query = Entrez.esearch(db='gds', term=geosample_id + " AND \"gsm\"", retmax=__RETMAX__)
        records = Entrez.read(query)
        return records

    def get_sra_from_gsm(self, gsm_id_list):
        """
        Get SRA links from GSM
        """
        query = Entrez.esummary(db='gds', id=(",").join(gsm_id_list))
        records = Entrez.read(query)
        return records

    def get_sra_metadata(self, sra_id):
        """Get SRA metadata"""
        query = Entrez.esearch('sra', term=sra_id)
        record = Entrez.read(query)
        query = Entrez.esummary(db='sra', id=record['IdList'])
        records = Entrez.read(query)
        return records


def stop_err(message):
     sys.stderr.write('ERROR: ' + message + '\n')
     sys.exit(1)


class MarkdownParser:

    """
    This class implements a very raw parser for parsing
    our markdown fikes. The current format used:

        |Project name|Content|GSE link|paper link|Method|Who|Finished|
        |:---:|:---:|:---:|:---:|:---:|:---:|:---:|

    The implementation assumes that the first occurence of '|' is a good enough
    catch to estimate the number of columns, then reads this line as its header
    and looks for the column number of 'GSE link'. The GEO links are then stored
    as a {'project_name': 'GEO accession key'}
    """

    def __init__(self, file_location):
        self.file_location = file_location
        assert(os.path.isfile(file_location))

    def is_row_junk(self, row):
        """
        Return True if a row has  at least one occurence of the  separator: ':---:'
        """
        if len(row) <= 7:
            return True
        for element in row:
            if element == ':---:':
                return True
        return False

    def get_geo_column_number(self, header_row):
        """
        Returns the column number where the __GEO_LINK_NAME__
        appears in the header row
        """
        for id, element in enumerate(header_row):
            if element == __GEO_LINK_NAME__:
                return id
        return None

    def get_header(self):
        """
        Moves the file pointer to next line till it encounters a header column as described in __init__
        """
        row = self.reader.next()
        while len(row) < 9:
            row = self.reader.next()
            if len(row) == 0:
                row = self.reader.next()
        return row

    def get_geo_links(self):
        """
        Returns a dictiornary with keys as the project name
        and values as the geo links
        """
        self.reader = csv.reader(open(self.file_location, 'r'), delimiter='|')
        header = self.get_header()
        assert(header is not None)
        column_number = __GEO_LINK_COLUMN__
        assert(column_number is not None)
        geo_links = {}
        for row in self.reader:
            if not self.is_row_junk(row):
                project_name = row[1]
                geo_link = row[column_number]
                geo_link = geo_link.split('=')
                try:
                    geo_link = geo_link[1]
                    geo_link = geo_link.replace(')', '')
                    geo_links[project_name] = geo_link
                except IndexError:
                    sys.stderr.write(
                        'ERROR parsing record for {} \n'.format(project_name))
        return geo_links


def is_already_downloaded(download_to):
    """
    Checks if the SRA being downloaded is already present,
    """
    if os.path.isfile(download_to):
        return True
    return False

def safe_to_download(temp_file, download_to):
    """
    Verifies if the checksums are different[To ensure reproducubility]
    """
    if not os.path.isfile(download_to):
        return True
    existing_hash = hashlib.sha256(open(download_to, 'rb').read()).digest()
    moving_hash = hashlib.sha256(open(temp_file, 'rb').read()).digest()
    if existing_hash == moving_hash:
        return True
    else:
        return False

class QueryProcessor:
    def __init__(self, ids_to_download):
        self.ids_to_download = ids_to_download
        self.download_queue = []
        self.geoQ = GEOQuery()
        QueryProcessor.ncbi_ftp = FTP(__NCBI_FTP__)
        # Set passive mode to keep Luigi et al. happy
        QueryProcessor.ncbi_ftp.set_pasv(True)
        QueryProcessor.ncbi_ftp.login()
        QueryProcessor.ncbi_ftp.voidcmd('TYPE I')

    def process(self):
        for gse_id in self.ids_to_download:
            gse_location = os.path.join(__ROOT_DOWNLOAD_LOCATION__, gse_id)
            if not os.path.exists(gse_location):
                os.makedirs(gse_location)
            else:
                sys.stderr.write("WARNING: {} location already exists\n".format(gse_location))
            gsm_records = self.geoQ.get_gsm_from_gse(gse_id)
            id_list = gsm_records['IdList']
            sra_records = geoQ.get_sra_from_gsm(id_list)
            for sra_record in sra_records:
                ext_relations = sra_record['ExtRelations']
                sra_title = sra_record['title'].replace(" ", "__") + "___" + sra_record['Accession']
                gsm_location = os.path.join(gse_location, sra_title)
                if not os.path.exists(gsm_location):
                    os.makedirs(gsm_location)
                sra_links = filter(lambda x: x['RelationType']=='SRA', ext_relations)
                assert(len(sra_links) == 1)
                target_ftp_link = sra_links[0]['TargetFTPLink']
                QueryProcessor.ncbi_ftp.cwd("/")
                srp = target_ftp_link.split("/")[-2]
                QueryProcessor.ncbi_ftp.cwd(target_ftp_link[32:])
                sra_ids = QueryProcessor.ncbi_ftp.nlst()
                for index, sra_id in enumerate(sra_ids):
                    title = sra_title
                    location = gsm_location
                    technical_replicate = 0
                    if len(sra_ids) > 1:
                        print "Found technical replicates at", target_ftp_link
                        title += "___R" + str(index+1)
                        location = os.path.join(gsm_location, "___R" + str(index+1))
                        technical_replicate = str(index+1)
                    QueryProcessor.ncbi_ftp.cwd(sra_id)
                    metadata = geoQ.get_sra_metadata(sra_id)
                    list_sras = QueryProcessor.ncbi_ftp.nlst()
                    assert(len(metadata)==1)
                    assert(len(list_sras)==1)

                    QueryProcessor.ncbi_ftp.voidcmd('TYPE I')
                    upstream_file_size = QueryProcessor.ncbi_ftp.size(list_sras[0])
                    strategy = str(BeautifulSoup(metadata[0]['ExpXml']).find('library_strategy').string)
                    download_dict = {'download_location': location, 'title': title, 'strategy': strategy, 'metadata': metadata,
                                     'target_ftp_link': target_ftp_link[32:] + sra_id + '/' + list_sras[0],
                                     'sra': list_sras[0],
                                     'upstream_file_size': upstream_file_size,
                                     'technical_replicate': technical_replicate,
                                     'sra_id': srp
                                     }
                    self.download_queue.append(download_dict)
                    QueryProcessor.ncbi_ftp.cwd("../")
        walkthrough(self.download_queue)
        with futures.ThreadPoolExecutor(max_workers=8) as executor:
            jobs =  [executor.submit(DownloadManager().download(record)) for record in self.download_queue]
        with open('result.yml', 'w') as yaml_file:
            yaml_file.write( yaml.dump(self.download_queue, default_flow_style=False))

def walkthrough(queue):
    sys.stdout.write('The following SRAs are to be downloaded.\nSRA\t\t\tFTPSize\tPartialDownloadSize\tSample ID\tTechnical Replicate(0=no technical replicate)\n')
    for record in queue:
        download_location = os.path.join(record['download_location'], record['sra'])
        upstream_size = record['upstream_file_size']
        downstream_size = 0
        if os.path.isfile(download_location):
            downstream_size = os.path.getsize(download_location)
        sys.stdout.write("{0}\t\t{1}\t{2}\t{3}\t{4}\n".format(record['sra'], size(upstream_size), size(downstream_size), record['sra_id'], record['technical_replicate']))


class DownloadManager(QueryProcessor):

    def __init__(self):
        pass

    def download(self, record):#download_location, upstream_file, sra_filename):
        """
        Pulls SRA from the FTP
        """
        self.download_location = record['download_location']
        self.ncbi_ftp = QueryProcessor.ncbi_ftp
        self.upstream_file = record['target_ftp_link']#upstream_file
        self.sra = record['sra']
        self.downstream_file = None
        self.upstream_file_size = record['upstream_file_size']
        self.download_to = os.path.join(self.download_location, self.sra)
        if not os.path.exists(self.download_location):
            os.makedirs(self.download_location)
        downstream_file_size = 0
        if is_already_downloaded(self.download_to):
            downstream_file_size = os.path.getsize(self.download_to)
        with open(self.download_to+'.metadata', 'wb') as f:
            json.dump(record, f)
        with open(self.download_to, 'w+b') as downstream_file:
            self.downstream_file = downstream_file
            self.downstream_file.seek(downstream_file_size)
            if downstream_file_size!=0:
                self.ncbi_ftp.retrbinary('RETR {}'.format(self.upstream_file), callback=self.download_block, blocksize=1024, rest=self.downstream_file.tell())
            else:
                self.ncbi_ftp.retrbinary('RETR {}'.format(self.upstream_file), callback=self.download_block, blocksize=1024)
        self.downstream_file.close()
        sys.stdout.write("\nWritten: {}\n".format(self.sra))

    def download_block(self, block):
        percentage = int(math.ceil(100*float(self.downstream_file.tell())/self.upstream_file_size))
        #sys.stdout.write("\r{0}: [{1}{2}] {3}%".format(self.sra, "#"*percentage, ' '*(100-percentage), percentage))
        #sys.stdout.flush()
        self.downstream_file.write(block)


if __name__ == '__main__':
    if __ROOT_DOWNLOAD_LOCATION__ == '':
        stop_err('__ROOT_DOWNLOAD_LOCATION__ not set')
    parser = argparse.ArgumentParser(
        description="Argument parser for downloading GEO files")
    parser.add_argument(
        '--dfile',
        type=str,
        help='Absolute path to markdown dataset file')
    parser.add_argument(
        '--gid',
        type=str,
        nargs='+',
        help="Space separated list of GEO Project IDs")
    parser.add_argument('--path', type=str, help="Root download location")
    args = parser.parse_args(sys.argv[1:])
    temp_dir = tempfile.mkdtemp()
    if args.dfile:
        fs = MarkdownParser(args.dfile)
        geo_links = fs.get_geo_links()
        project_names = geo_links.keys()
        ids_to_download = [geo_links[key] for key in project_names]
    elif args.gid:
        geoQ = GEOQuery()
        ids_to_download = args.gid
        _set_root_download_path(args.path)
    else:
        geoQ = GEOQuery("\"Methylation profiling by high throughput sequencing\"[DataSet Type]")
        records = geoQ.submit_query()
        id_list = records['IdList']
        records = geoQ.get_titles(id_list)
        for i, record in enumerate(records):
            print "{}: {}: {}".format(i+1, record['Accession'].encode('utf-8'), record['title'].encode('utf-8'))
        print("Enter 0 to download all(Warning!)\n"
            "Enter 1-3 to download datasets from 1-3\n"
            "Enter GSE25836, GSE47752 to download GSE25836 and GSE47752")
        to_download = raw_input("")
        to_download = to_download.replace(" ", "")
        if "-" in to_download:
            split = to_download.split('-')
            start = int(split[0])-1
            if split[1]=="":
                end = len(records)-1
            else:
                end = int(split[1])-1
            range_to_download = range(start, end+1)
            ids_to_download = [records['Accession'].encode('utf-8') for i, record in enumerate(records) if i in range_to_download]
        elif "GSE" in to_download:
            ids_to_download = to_download.split(',')
        elif to_download == "0":
            ids_to_download = [records['Accession'].encode('utf-8') for i, record in enumerate(records)]
    qp = QueryProcessor(ids_to_download)
    qp.process()
