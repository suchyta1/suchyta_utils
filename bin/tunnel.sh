#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
temp="temp-tunnel-var.sh"

$DIR/tunnel.py $@ --tempfile $temp
source $temp
rm $temp
