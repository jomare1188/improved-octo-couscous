#!/bin/bash
echo "host=$(hostname) NSLOTS=$NSLOTS TMPDIR=$TMPDIR"
source /etc/profile.d/modules.sh
which qsub qstat qdel 2>&1
ls -d /Storage/data1/jorge.munoz/DOLORES/nf-rnaseq/conf_cluster && echo shared_fs_ok
module load singularity-ce/4.3.3 2>&1; singularity --version 2>&1
timeout 15 curl -sI https://depot.galaxyproject.org/singularity/ 2>&1 | head -1 || echo no_internet_galaxy
timeout 15 curl -sI https://community-cr-prod.seqera.io/v2/ 2>&1 | head -1 || echo no_internet_seqera
df -h /tmp | tail -1
ulimit -a | grep -E 'open files|max user processes'
