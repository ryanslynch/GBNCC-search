#!/usr/bin/env python
import sys
import argparse
import copy
import multiprocessing
import traceback
import textwrap
import warnings

import candidate
import sp_raters
import utils

def warn_to_stdout(message, category, filename, lineno, file=None, line=None):
    """A function to replace warnings.showwarning so that warnings are
         printed to STDOUT instead of STDERR.
         Usage: warnings.showwarning = warn_to_stdout
    """
    sys.stdout.write(warnings.formatwarning(message,category,filename,lineno))

def rate_spd(spdfn, rater_instances):
    """Given the name of an *.spd file and a list of Rater instances
        compute the ratings.
        Inputs:
            spdfn: Name of the *.spd file.
            rater_instances: A list of Rater instances to compute ratings.
        Outputs:
            cand: The resulting (rated) Candidate object.
        ***NOTE: RatingValues are added directly to the Candidate object.
    """
    cand = candidate.read_spd_file(spdfn)
    for rater in rater_instances:
        ratval = rater.rate(cand)
        cand.add_rating(ratval)
    return cand


def main():
    if not args.raters:
        print "No raters are loaded."
        args.list_raters = True

    if args.list_raters:
        utils.print_sp_raters_list(args.verbosity)
        sys.exit(0)

    if args.ignore_warnings:
        warnings.simplefilter('ignore', utils.RatingWarning)

    if args.redirect_warnings:
        warnings.showwarning = warn_to_stdout

    rater_instances = []
    for rater_name in args.raters:
        rater_module = getattr(sp_raters, rater_name)
        rater_instances.append(rater_module.Rater())
    
    cands = []
    if args.num_procs > 1:
        print "Using %d rater threads" % args.num_procs
        rater_pool = multiprocessing.Pool(processes=args.num_procs)
        apply_async = lambda spdfn: rater_pool.apply_async(rate_spd, \
                                                (spdfn, rater_instances))
        inprogress = [apply_async(spdfn) for spdfn in args.infiles]
        failed = []
        while inprogress:
            for ii in range(len(inprogress))[::-1]:
                result = inprogress[ii]
                if result.ready():
                    if result.successful():
                        cands.append(result.get())
                    else:
                        failed.append(result)
                    inprogress.pop(ii)
        print "Number of failures detected: %d" % len(failed)
        for fail in failed:
            try:
                print type(fail), fail
                fail.get()
            except:
                traceback.print_exc()
    else:
        for spdfn in args.infiles:
            cands.append(rate_spd(spdfn, rater_instances))
    for cand in cands:
        if args.write_to_file:
            cand.write_ratings_to_file()
        if args.write_to_screen:
            print cand.spdfn
            print cand.get_ratings_overview()
            print '-'*25

class RemoveAllRatersAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, [])


class AddAllRatersAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        ratlist = copy.deepcopy(sp_raters.registered_raters)
        setattr(namespace, self.dest, ratlist)


class RemoveOneRaterAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        rater_name = values[0]
        if rater_name not in sp_raters.registered_raters:
            sys.stderr.write("Unrecognized rater: %s\nThe following " \
                             "raters are registered:\n    %s\n" % \
                             (rater_name, "\n    ".join(sp_raters.registered_raters)))
            sys.exit(1)
        curr_raters = getattr(namespace, self.dest)
        # Remove any instances of 'values' from curr_raters
        while rater_name in curr_raters:
            curr_raters.remove(rater_name)
        setattr(namespace, self.dest, curr_raters)


class AddOneRaterAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        rater_name = values[0]
        if rater_name not in sp_raters.registered_raters:
            sys.stderr.write("Unrecognized rater: %s\nThe following " \
                             "raters are registered:\n    %s\n" % \
                             (rater_name, "\n    ".join(sp_raters.registered_raters)))
            sys.exit(1)
        curr_raters = getattr(namespace, self.dest)
        # Add 'rater_name' to curr_raters
        if rater_name not in curr_raters:
            curr_raters.append(rater_name)
        setattr(namespace, self.dest, curr_raters)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rate SPD files.")
    parser.add_argument('infiles', metavar="INFILES", \
                         nargs='+', type=str, \
                         help="Single Pulse *.spd files to rate.")
    parser.add_argument('-v', '--more-verbose', dest='verbosity', \
                         default=0, action='count', \
                         help="Turn up verbosity by one notch. " \
                                "(Default: Don't be verbose (verbosity=0).)")
    parser.add_argument('-d', '--debug', dest='debug', \
                         default=False, action='store_true', \
                         help="Turn on debugging output. " \
                                "(Default: Don't print debugging info.)")
    parser.add_argument('-L', '--list-raters', dest='list_raters', \
                        default=False, action='store_true', \
                        help="List registered raters and exit.")
    parser.add_argument('-x', '--exclude', dest='raters', \
                         type=str, default=[], nargs=1, \
                         action=RemoveOneRaterAction, \
                         help="Remove rater from list of ratings to apply.")
    parser.add_argument('--exclude-all', dest='raters', \
                         default=[], nargs=0, \
                         action=RemoveAllRatersAction, \
                         help="Clear list of ratings to apply.")
    parser.add_argument('-i', '--include',  dest='raters', \
                         type=str, default=[], nargs=1, \
                         action=AddOneRaterAction, \
                         help="Include rater from list of ratings to apply.")
    parser.add_argument('--include-all', dest='raters', \
                         default=[], nargs=0, \
                         action=AddAllRatersAction, \
                         help="Include all registered ratings in list " \
                                "of ratings to apply.")
    parser.add_argument('-P', '--num-procs', dest="num_procs", \
                        type=int, default=1, \
                        help="The number of rater processes to use. " \
                                "Each thread rates a separate candidate. " \
                                "(Default: use one rater thread.)")
    parser.add_argument('--no-write-to-file', dest='write_to_file', \
                        default=True, action='store_false', \
                        help="Do not write out ratings to file. " \
                                "(Default: Write out file.)")
    parser.add_argument('--no-write-to-screen', dest='write_to_screen', \
                        default=True, action='store_false', \
                        help="Do not write out ratings to screen. " \
                                "(Default: Write ratings to screen.)")
    parser.add_argument('--ignore-warnings', dest='ignore_warnings', \
                        default=False, action='store_true', \
                        help="Do not display RatingWarning warnings. " \
                                "(Default: display warnings.)")
    parser.add_argument('--redirect-warnings', dest='redirect_warnings', \
                        default=False, action='store_true', \
                        help="Redirect warnings to stdout. " \
                                "(Default: warnings print to stderr.)")
    args = parser.parse_args()
    main()
