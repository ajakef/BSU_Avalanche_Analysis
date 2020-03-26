#!/bin/bash
while :
do
    echo "starting job"
    python $1 $2 $3
    echo "error: hit ctrl-c twice to stop job"
    sleep 2
done
