import urllib.request
import os
import gzip
import shutil

class FileDownloader:

    def single_file(self, source):

        if source == 'string_edges':

            # check to see if file already exists
            # TODO: check file version too
            file = os.path.join('data', 'input_files', 'string_edges.txt.gz')
            if os.path.exists(file):
                return None

            url = 'http://viruses.string-db.org/download/protein.links.full.v10.5.txt.gz'
            print('Downloading file from: ' + url)

            # download file
            urllib.request.urlretrieve(url, os.path.join('data', 'input_files', 'string_edges.txt.gz'))

            # unzip
            with gzip.open(os.path.join('data', 'input_files', 'string_edges.txt.gz'), 'rb') as f_in:
                with open(os.path.join('data', 'input_files', 'string_edges.txt'), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            print('Done!')