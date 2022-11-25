#!/usr/bin/env python3
from __future__ import print_function
import sys
min_python = (3,6)
try:
    assert(sys.version_info >= min_python)
except:
    from platform import python_version
    print("Python version", python_version(), "is too old.")
    print("Please use Python", ".".join(map(str,min_python)), "or later.")
    sys.exit()

import argparse
import atexit
import glob
import json
import os
import platform
import queue
import re
import shutil
import subprocess
import tarfile
import threading
import time
import tempfile
import contextlib

from io import open
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, urlretrieve, Request
from urllib.error import HTTPError

class urlopen_progress:
    timeout = 60
    retries = 10

    def __init__(self, url, quiet, teamcity):
        self.url = url
        self.bytes_so_far = 0
        self.urlopen()

        self.quiet = quiet
        self.teamcity = teamcity
        if teamcity:
            self.EOL = '\n'
        else:
            self.EOL = '\r'
        self.cur_row = -1
        total_size = 0
        total_size = self.remote_file.getheader('Content-Length', 0) # More modern method

        self.total_size = int(total_size)

    def urlopen(self):
        headers = dict()
        if self.bytes_so_far > 0:
            headers['Range'] = 'bytes={}-'.format(self.bytes_so_far)
        request = Request(self.url, headers=headers)
        self.remote_file = urlopen(request, timeout=self.timeout)

    def read(self, n=8388608):
        delay = 1
        for attempt in range(self.retries):
            try:
                if self.remote_file is None:
                    self.urlopen()
                buffer = self.remote_file.read(n)
                if not buffer:
                    if not self.quiet:
                        sys.stdout.write('\n')
                    return ''
                break
            except Exception as ex:
                self.remote_file = None
                time.sleep(delay)
                delay += delay

        self.bytes_so_far += len(buffer)
        percent = float(self.bytes_so_far) / self.total_size
        percent = round(percent*100, 2)

        do_print = True
        if self.teamcity:
            do_print = False
            row = int(percent)
            if row > self.cur_row:
                self.cur_row = row
                do_print = True

        if do_print and not self.quiet:
            sys.stderr.write("Downloaded %d of %d bytes (%0.2f%%)%s" % (self.bytes_so_far, self.total_size, percent, self.EOL))

        return buffer

def install_url(url, path, quiet, teamcity):
    basename = os.path.basename(urlparse(url).path)
    try:
        local_file =  os.path.join(path, basename)
        if os.path.exists(local_file):
            if not quiet:
                print('Extracting local tarball: {}'.format(local_file))
            fileobj = open(local_file, 'rb')
        else:
            if not quiet:
                print('Downloading and extracting tarball: {}'.format(url))
            fileobj = urlopen_progress(url, quiet, teamcity)
        with tarfile.open(mode='r|*', fileobj=fileobj) as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar, path=path)
    except:
        sys.stderr.write('''
ERROR: Failed to extract tarball; to install manually, try something like:
    curl -OLC - {}
    tar xvf {}
'''.format(url, basename))
        raise

def quiet_remove(filename):
    with contextlib.suppress(FileNotFoundError):
        os.remove(filename)

def find_failed_step(filename):
    r = "^\[(?P<time>[^\]]+)\] (?P<level>[^ ]+) \[(?P<source>[^ ]*) (?P<name>[^\]]*)\] (?P<status>.*)"
    search = re.compile(r)
    lines = open(filename, "r").readlines()
    nameStarts = {}
    start = -1
    for num, line in enumerate(lines):
        r = search.match(line)
        if r:
            name = r.group("name")
            if name not in nameStarts:
                nameStarts[name] = num

            if r.group("status") == "completed permanentFail":
                start = nameStarts[name]
                break
    if start > -1:
        print("Printing log starting from failed job:\n")
        for i in range(start, len(lines)):
            print(lines[i], end="")
    else:
        print("Unable to find error in log file.")


class Setup:

    def __init__(self, args):
        self.args = args
        self.outputdir = self.get_output_dir()
        #self.rundir = '.'

        # Create a work directory.
        if args.input:
            print("Output will be placed in:", self.outputdir)
            os.mkdir(self.outputdir)


    def get_output_dir(self):
        outputdir = os.path.abspath(self.args.output)
        if os.path.exists(outputdir):
            parent, base = os.path.split(outputdir)
            counter = 0
            for sibling in os.listdir(parent):
                if sibling.startswith(base + '.'):
                    ext = sibling[len(base)+1:]
                    if ext.isdecimal():
                       counter = max(counter, int(ext))
            outputdir = os.path.join(parent, base+'.'+str(counter+1))
        return outputdir

    def install_tgz(self):
      guard_file = "vecscreen_x.tgz"
      remote_path = f"https://s3.amazonaws.com/sapojnik-dev/vecscreen_x.tgz"
      if not os.path.isfile(guard_file):
          #install_url(remote_path, self.rundir, self.args.quiet, False)
          install_url(remote_path, '.', self.args.quiet, False)
          open(guard_file, 'a').close()
      else:
          print(f"Skipping already installed tarball: {remote_path}")

class Pipeline:

    def __init__(self, params, local_input):
        self.params = params
        self.input_file = local_input
        if local_input=='': self.input_file='samples/GCA_000173135.1_reduced.1.fna'
        self.cmd = ['bin/blastn']
        self.cmd.extend(['-db',
                        'CommonContaminants/gcontam1',
                        '-query',
                        self.input_file,
                        '-perc_identity',
                        '90.0',
                        '-outfmt',
                        '6',
                        '-out', params.outputdir + '/vecscreen.out'
                        ])

    def launch(self):
        cwllog = self.params.outputdir + '/vecscreen.log'
        with open(cwllog, 'a', encoding="utf-8") as f:
            # Show original command line in log
            cmdline = "Original command: " + " ".join(sys.argv)
            f.write(cmdline)
            f.write("\n\n")
            # Show docker command line in log
            cmdline = "Executing: " + " ".join(self.cmd)
            f.write(cmdline)
            f.write("\n\n")
            f.flush()
            try:
                proc = subprocess.Popen(self.cmd, stdout=f, stderr=subprocess.STDOUT)
                proc.wait()
            finally:
                if proc.returncode == None:
                    print('\nAbnormal termination, stopping all processes.')
                    proc.terminate()
                elif proc.returncode == 0:
                    print(f'Completed successfully.')
                else:
                    print(f'Failed, sub-command exited with rc =', proc.returncode)
                    #find_failed_step(cwllog)
        return proc.returncode

def main():
    parser = argparse.ArgumentParser(description='Run vecscreen_x.')
    parser.add_argument('input', nargs='?',
                        help='Input FASTA file.')
    parser.add_argument('-o', '--output', metavar='path', default='output',
                        help='Output directory to be created, which may include a full path')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Quiet mode, for scripts')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Debug mode')

    args = parser.parse_args()

    retcode = 0
    try:
        params = Setup(args)

        params.install_tgz()

        if args.input:
          print("Executing vecscreen")
          p = Pipeline(params, args.input)
          # $BLASTN -db $BLAST_DB_DIR/gcontam1 -query $SAMPLES -perc_identity 90.0 -outfmt 6 -out $OUT
          p.launch()


    except (Exception, KeyboardInterrupt) as exc:
        if args.debug:
            raise
        retcode = 1
        print(exc)

    sys.exit(retcode)

if __name__== "__main__":
    main()
