#!/usr/bin/env bash
# Submit LSF job to run Picard liftOver

bsub \
-q inter \
-P re_gecip_enhanced_interpretation \
-J "liftover" \
-o data/logs/%J_liftover_out \
-e data/logs/%J_liftover_err \
bash src/process_dnms/picard_liftover.sh $1 $2

bwait -w "done(liftover)"

# Compress and index the output
bgzip --keep -f "data/interim/${2}"
tabix "data/interim/${2}.gz"