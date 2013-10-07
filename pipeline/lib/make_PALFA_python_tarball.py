#!/usr/bin/env python
import shutil, os, glob
import tarfile
import config
import distutils.dir_util

presto_dir = os.getenv('PRESTO')

# binaries to put in tarball bin dir (if folder will copy contents of it)
bins = [ os.path.join(config.pipelinedir,'bin/search.py'),
         os.path.join(presto_dir,'bin/single_pulse_search.py') ]

# individual files to put in tarball site packages dir
site_packages_file = ['/sb/project/bgf-180-aa/lib/python2.6/site-packages/memory.py',
                      '/sb/project/bgf-180-aa/lib/python2.6/site-packages/argparse.py']

# dirs from which to copy the contents of into tarball site packages dir
site_packages_contents = [ os.path.join(config.pipelinedir,'lib/python/'),
                           os.path.join(presto_dir,'lib/python/') ]

# dirs to copy to the tarball site packages dir
# no trailing slash because copytree is silly
site_packages_folder = ['/sb/project/bgf-180-aa/python_tarballs/for_tarball/ubc_AI',
                        '/sb/project/bgf-180-aa/lib/python2.6/site-packages/pyfits',
                        '/sb/project/bgf-180-aa/python_tarballs/for_tarball/ratings2' ]

# filenames of input (base) and output tarballs
base_tarball = "/sb/project/bgf-180-aa/python_tarballs/python-2.6.8-packages-base.tar.gz"
new_tarball = "/sb/project/bgf-180-aa/python_tarballs/python-2.6.8-packages-GBNCC.tar.gz"

def make_tarball():
    # extract base python tarball
    py_tar = tarfile.open(base_tarball)
    py_tar.extractall()
    py_tar.close()
    
    bin_dir = 'python-2.6.8/bin/'
    packages_dir = 'python-2.6.8/lib64/python2.6/site-packages/'
    
    # copy binaries to python tarball bin dir
    for bin in bins:
        if os.path.isdir(bin):
            distutils.dir_util.copy_tree(bin,bin_dir)
        else:
            shutil.copy(bin,bin_dir)
            
    # copy packages to python tarball site packages bin dir
    for package in site_packages_file:
            shutil.copy(package,packages_dir)
    for package in site_packages_folder:
            new_dir = os.path.join(packages_dir,os.path.basename(package))
            shutil.copytree(package,new_dir)
    for package in site_packages_contents:
            distutils.dir_util.copy_tree(package,packages_dir,preserve_symlinks=1)
    
    # remove all .pyc's
    pyc_glob = glob.glob(os.path.join(packages_dir,'*.pyc'))
    for pyc in pyc_glob:
        os.remove(pyc)
    
    # tar up the new tarball
    tar = tarfile.open(new_tarball, "w:gz")
    tar.add('python-2.6.8/')
    tar.close()
    
    # remove extracted tarball dir
    shutil.rmtree('python-2.6.8/')

if __name__ == '__main__':
    make_tarball()
