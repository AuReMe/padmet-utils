import argparse
import os
import queue
import sys
import time
from threading import Thread

from ceterach import exceptions as exc
from ceterach.api import MediaWiki

# Loads Aureme results into a mediawiki instanace.
# Uses the ceterach python module: https://github.com/abretaud/ceterach/ (fork from https://github.com/Riamse/ceterach/ with a fix for uploads)
# Requires python 3
# Tested with ceterach commit c890ad91fa0aba845791958fe47a12672a6f0076 and MediaWiki 1.28.1
#
# Install & usage:
#  virtualenv .venv
#  . .venv/bin/activate
#  pip install git+https://github.com/abretaud/ceterach.git
#  python load_page.py


class MWLoader(Thread):

    def __init__(self, mw, queue):

        Thread.__init__(self, daemon=True)

        self.mw = mw
        self.q = queue

    def run(self):
        while True:
            item = self.q.get()
            if item is None:
                break
            if item[0] == 'page':
                self.load_page(item[1])
                print('    %s done (%s remaining)' % (item[1], self.q.qsize()))
            else:
                self.upload_file(item[1])
                print('    %s file uploaded (%s remaining)' % (item[1], self.q.qsize()))
            self.q.task_done()

    def load_page(self, filepath):

        title = os.path.basename(filepath).replace("__47__", "/")
        retries = 0
        success = False
        last_exc = None
        while not success and retries < 5:
            try:
                self._load_page(title, filepath)
                success = True
            except (exc.PermissionsError, exc.NonexistentPageError, exc.InvalidPageError, exc.EditError, exc.EditFilterError, exc.SpamFilterError, exc.EditConflictError) as e:
                print(e)
                retries += 1
                print('Failed, retrying in %s seconds...' % retries)
                last_exc = e
                time.sleep(retries)

                # It can be caused by a token too old, generate a new one just in case
                del self.mw.tokens['edit']
                self.mw.set_token("edit")

        if not success:
            if last_exc:
                print('Failed to load page "%s" (%s).' % (title, filepath), file=sys.stderr)
                raise last_exc
            else:
                # Not supposed to happen, but who knows...
                raise Exception('Failed to load page "%s" (%s).' % (title, filepath))

    def _load_page(self, title, filepath, summary=''):

        p = self.mw.page(title)
        with open(filepath, 'r') as content_file:
            if p.exists:
                p.edit(content_file.read(), summary, minor=True, force=True)
            else:
                p.create(content_file.read(), summary, minor=True, force=True)

    def upload_file(self, filepath):

        title = os.path.basename(filepath).replace("__47__", "/")
        retries = 0
        success = False
        last_exc = None
        while not success and retries < 5:
            try:
                self._upload_file(title, filepath)
                success = True
            except (exc.PermissionsError, exc.NonexistentPageError, exc.InvalidPageError, exc.EditError, exc.EditFilterError, exc.SpamFilterError, exc.EditConflictError, exc.CeterachError) as e:
                print(e)
                retries += 1
                print('Failed, retrying in %s seconds...' % retries)
                last_exc = e
                time.sleep(retries)

                # It can be caused by a token too old, generate a new one just in case
                del self.mw.tokens['edit']
                self.mw.set_token("edit")

        if not success:
            if last_exc:
                print('Failed to upload file "%s" (%s).' % (title, filepath), file=sys.stderr)
                raise last_exc
            else:
                # Not supposed to happen, but who knows...
                raise Exception('Failed to upload file "%s" (%s).' % (title, filepath))

    def _upload_file(self, title, filepath, text='', summary=''):

        p = self.mw.file(title)

        with open(filepath, 'rb') as content_file:
            p.upload(content_file, text, summary)


def list_pages(dirname):

    aureme_dirs_to_load = ['genes', 'metabolites', 'pathways', 'reactions', 'navigation']

    for d in aureme_dirs_to_load:
        path_to_files = os.path.join(dirname, d)

        if os.path.exists(path_to_files):
            list_fns = os.listdir(path_to_files)
            for filename in list_fns:

                full_path = os.path.join(path_to_files, filename)
                yield full_path


def list_files(dirname):

    aureme_dirs_to_load = ['files', 'forms']

    for d in aureme_dirs_to_load:
        path_to_files = os.path.join(dirname, d)

        if os.path.exists(path_to_files):
            list_fns = os.listdir(path_to_files)
            for filename in list_fns:

                full_path = os.path.join(path_to_files, filename)
                yield full_path


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to Aureme output dir")
    parser.add_argument("url", help="Url to the MediaWiki api (should have api.php at the end for a standard install)")
    parser.add_argument("login", help="MediaWiki user login")
    parser.add_argument("password", help="MediaWiki user pasword")
    parser.add_argument("-t", "--threads", help="Number of threads (default: 1)", default="1", type=int)
    args = parser.parse_args()

    path = args.path
    url = args.url
    login = args.login
    password = args.password
    num_worker_threads = args.threads

    mw = MediaWiki(url)
    mw.login(login, password)

    q = queue.Queue()

    # Start all the threads
    threads = []
    for i in range(num_worker_threads):
        t = MWLoader(mw, q)
        t.start()
        threads.append(t)

    # Fill the queues
    for item in list_files('./'):  # Load forms in local dir
        q.put(('file', item))
    for item in list_files(path):
        q.put(('file', item))
    for item in list_pages(path):
        q.put(('page', item))

    # Wait for queues to be empty
    q.join()
