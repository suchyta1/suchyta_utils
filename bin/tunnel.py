#!/usr/bin/env python


import argparse
import subprocess
import os
import getpass
import storm


def GetHostInfo(match, info):
    if 'options' not in match.keys():
        return

    opts = match['options']
    opts_keys = opts.keys()

    if 'hostname' in opts_keys:
        info['hostname'] = opts['hostname']
    if 'user' in opts_keys:
        info['user'] = opts['user']


def Overwrite(info, key, value):
    if value is not None:
        info[key] = value


def LookupHost(host, user, config):
    info = {'hostname': host,
            'user': None}

    matches = config.search_host(host)
    for match in matches:
        if match['host'] == host:
            GetHostInfo(match, info)
            break

    Overwrite(info, "user", user)
    return info





if __name__ == "__main__":

    duser = getpass.getuser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-hu", "--huser", help="Username on remote host you want to connect to", default=None)
    parser.add_argument("-gu", "--guser", help="Username on the gateway", default=None)

    parser.add_argument("-p", "--port", help="Local port to forward",  default=3000, type=int)
    parser.add_argument("-H", "--host", help="Remote host you want to connect to", required=True)
    parser.add_argument("-g", "--gateway", help="Gateway you want to connect through", required=True)

    parser.add_argument("-e", "--env", help="Environment variable name", default=None)
    parser.add_argument("-t", "--tempfile", help="Temp file name", default="temp-tunnel-var.sh")
    args = parser.parse_args()

    home = os.path.expanduser("~")
    ssh_config = os.path.realpath(os.path.join(home, ".ssh", "config"))
    config = storm.ConfigParser(ssh_config)
    config.load()

    hinfo = LookupHost(args.host, args.huser, config)
    ginfo = LookupHost(args.gateway, args.guser, config)

    cmd = ["ssh", "-f", "-N",
           "-L", "localhost:{0}:{1}:22".format(args.port, hinfo['hostname']),
           "{0}@{1}".format(ginfo['user'], ginfo['hostname'])]

    print ' '.join(cmd)
    subprocess.call(cmd)
    print "Tunnel opened to {0}".format(hinfo['hostname'])

    shortname = args.env
    if shortname is None:
        shortname = hinfo['hostname']
        ind = shortname.find('.')
        if ind != -1:
            shortname = hinfo['hostname'][0:ind]
    shortname = shortname.upper()

    out = open(args.tempfile, "w")
    out.write('export {0}="-p {1} {2}@localhost"'.format(shortname, args.port, hinfo['user']))
    out.close()

    #print "To connect: ssh -p {0} {1}@localhost ".format(args.port, hinfo['user'])
    print "To connect: ssh ${0}".format(shortname)
