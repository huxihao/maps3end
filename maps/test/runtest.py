#! /usr/bin/env python
"""
Test runner:
Collects all files ending in _test.py and executes them.
"""

import os, sys, re, unittest

def walk(rootdir):
    filenames = []
    for root, dirs, files in os.walk(rootdir):
        for name in files:       
            filename = os.path.join(root, name)
            filenames.append(filename)
    return filenames

def all_tests():
    "Returns all file names that end in _test.py"
    patt = re.compile("_test.py$")
    mods = os.listdir(os.path.normpath(os.path.dirname(__file__)))
    mods = filter(patt.search, mods)

    # some predictable order...
    mods.sort() 
    return mods

def run(targets):
    "Imports and runs the modules names that are contained in the 'targets'"
    
    success = errors = 0

    # run the tests by importing the module and getting its test suite
    for name in targets:
        try:
            print >> sys.stderr, 'testing in %s.py ...' % name
            l = unittest.TestLoader()
            suite = l.loadTestsFromName(name)
            runner = unittest.TextTestRunner()
            results = runner.run(suite)
            
            # count tests and errors
            success += results.testsRun - \
                       len(results.errors) - \
                       len(results.failures)
            
            errors  += len(results.errors) + len(results.failures)

        except ImportError:
            print >> sys.stderr, "unable to import module '%s'" % name

    # summarize the run
    print >> sys.stderr, '=' * 59
    print >> sys.stderr, '''\
%s tests passed, %s tests failed, %d total''' % \
                  (success, errors, success + errors)

    return (success, errors)

if __name__ == '__main__':
    # modules: from command line args or all modules
    targets = all_tests()
    print targets
    # get rid of the .py ending in case full module names were passed in
    # the command line
    stripped_targets = []
    for t in targets:
        if t.endswith('.py'): t = t[:-3]
        stripped_targets.append(t)
    targets = stripped_targets

    good, bad = run(targets)
    #l = unittest.TestLoader()
    #suite = l.discover('./', pattern='_test_.py$')
    #runner = unittest.TextTestRunner()
    #results = runner.run(suite)



    if bad:
        sys.exit(-1)
        
    sys.exit(0)
