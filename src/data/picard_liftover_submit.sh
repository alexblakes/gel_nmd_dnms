#!/usr/bin/env bash
# Submit LSF job to run Picard liftOver

bsub \
-q inter \
-P re_gecip_enhanced_interpretation \
-J "liftover" \
-o data/logs/%J_liftover_out \
-e data/logs/%J_liftover_err \
-R "rusage[mem=64000]" \
-M 64000 \
bash src/data/picard_liftover.sh $1 $2

bwait -w "done(liftover)"