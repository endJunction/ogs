#!/usr/bin/env bash

HASH=$1

echo -n `git log --format=%h -n 1 $HASH`

if [[ -f test_revs/${HASH} ]]; then
    if head -n 1 test_revs/${HASH} | grep -q ^sys_user_ ; then
        echo -n " OK "
        cat test_revs/${HASH} | sed 's/_/ /g' | awk '{print $3 " " $4 " " $6 " " $8 " "$9}'
    else
        echo " ERROR"
    fi
else
    echo " BAD"
fi

