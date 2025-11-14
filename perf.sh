#!/bin/bash

sudo rm -f perf.data
sudo rm -f out.perf-folded
sudo rm -f perf.svg
sudo perf record -F 100 -g --call-graph dwarf,65528 ./ecq
sudo perf script | ./FlameGraph/stackcollapse-perf.pl > out.perf-folded
./FlameGraph/flamegraph.pl out.perf-folded > perf.svg