#!/usr/bin/env python


import argparse
import subprocess
import os
import sys
import shutil


def RunAndLog(cmd, logger):
    logger.write( "\nCommand:\n{0}\n\n".format(" ".join(cmd)) )
    logger.flush()
    retcode = subprocess.call(cmd, stdout=logger, stderr=logger)

    if retcode != 0:
        logger.flush()
        logger.seek(0)
        history = logger.read()
        sys.stderr.write("Error!\nCommand history:\n\n{}".format(history))
        logger.close()
        sys.exit(retcode)


def AppendArgs(args, cmd):
    for con in  args:
        configs = con.split()
        for config in configs:
            cmd.append(config)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--wget", help="URL to wget", default=None)
    parser.add_argument("-g", "--git", help="URL to git clone", default=None)

    parser.add_argument("-c", "--config", help='Configure install and settings. Must be given like --config or --config="--prefix=/path".', default=None, nargs="*")
    parser.add_argument("-p", "--python", help='Python install and setting. Must be given like --python or --python="--prefix=/path".', default=None, nargs="*")

    parser.add_argument("-l", "--logdir", help="Directory where to write logfile", default="./")
    parser.add_argument("-k", "--keep", help="Keep source", action="store_true")
    args = parser.parse_args()


    if args.wget is not None:
        name = args.wget
    if args.git is not None:
        name = args.git
    fname = name.split('/')[-1]
    logfile = os.path.join(args.logdir, "{0}.log".format(fname))
    logger = open(logfile, "w+")

    if args.wget is not None:

        cmd = ["wget", args.wget]
        RunAndLog(cmd, logger)

        p = subprocess.Popen( ["tar", "tzf",  fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        out, err = p.communicate()
        codedir = out.strip().split()[0]

        cmd = ["tar", "xvzf", fname]
        RunAndLog(cmd, logger)
        os.remove(fname)

    if args.git is not None:
        cmd = ["git", "clone", args.git]
        RunAndLog(cmd, logger)
        codedir = fname.rstrip('.git')


    os.chdir(codedir)

    if args.python is not None:
        cmd = ["python", "setup.py", "build"]
        AppendArgs(args.python, cmd)
        RunAndLog(cmd, logger)

        cmd = ["python", "setup.py", "install"]
        RunAndLog(cmd, logger)

    if args.config is not None:
        cmd = ["./configure"]
        AppendArgs(args.config, cmd)
        RunAndLog(cmd, logger)

        cmd = ["make"]
        RunAndLog(cmd, logger)
        cmd = ["make", "install"]
        RunAndLog(cmd, logger)

    os.chdir("../")
    if not args.keep:
        shutil.rmtree(codedir)

    logger.seek(0)
    history = logger.read()
    sys.stdout.write(history)
