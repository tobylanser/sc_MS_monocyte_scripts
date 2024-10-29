#!/bin/bash

compass \
--data compass_input/myeloid_counts.txt \
--species homo_sapiens \
--num-threads 28 \
--output-dir COMPASS_outs/ \
--microcluster-size 15000
