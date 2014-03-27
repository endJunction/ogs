#!/usr/bin/env bash

#set -x

#
# Setup and check the executable.
#
BINARY_DIR="./bin"
EXECUTABLE=${BINARY_DIR}/"SimpleGroundWaterFlow"
if [[ ! -f ${EXECUTABLE} ]]; then
    echo "${EXECUTABLE} does not exist."
    exit 1
fi

if [[ ! -e ${EXECUTABLE} ]]; then
    echo "${EXECUTABLE} is not executable."
    exit 1
fi

#
# Check required input files.
#
CASE=$1
if [[ -z ${CASE} ]]; then
    echo "No CASE given."
    exit 2
fi

CASE_DIR=`dirname ${CASE}`
if [[ ! -d ${CASE_DIR} ]]; then
    echo "${CASE_DIR} is not accessable."
    exit 2
fi

if [[ ! -f "${CASE}.cnd" ]]; then
    echo "Could not access ${CASE}.cnd file."
    exit 2
fi
if [[ ! -f "${CASE}.gml" ]]; then
    echo "Could not access ${CASE}.gml file."
    exit 2
fi
if [[ ! -f "${CASE}.vtu" ]]; then
    echo "Could not access ${CASE}.vtu file."
    exit 2
fi

ARGS="--boundary_condition ${CASE}.cnd -g ${CASE}.gml -m ${CASE}.vtu"

#
# Wrapper definitions.
#
WRAPPER_MEMCHECK="valgrind --tool=memcheck --log-file=${CASE}_valgrind.log -v \
                --suppressions=${BINARY_DIR}/../../valgrind.suppressions \
                --leak-check=full \
                --show-reachable=yes \
                --track-origins=yes \
                --malloc-fill=0xff \
                --free-fill=0xff"

WRAPPER_CALLGRIND="valgrind --tool=callgrind --branch-sim=yes --cache-sim=yes --dump-instr=yes --collect-jumps=yes"
WRAPPER_TIME="time --format=sys_user_%S_%U__mem_%M__io_%I_%O"

#
# Choose wrapper depending on cli option, or pass forward if not an option.
#
WRAPPER="$2"
if [[ ! -z "$2" ]]; then
    case "$2" in
        (-m)
            WRAPPER=${WRAPPER_MEMCHECK}
            ;;
        (-c)
            WRAPPER=${WRAPPER_CALLGRIND}
            ;;
        (-t)
            WRAPPER=${WRAPPER_TIME}
            ;;
    esac
fi

#
# Execute program in wrapper
#
${WRAPPER} ${EXECUTABLE} ${ARGS}


#
# Different tester for result/output comparison/checking.
#
COMPARE_VTU_OUTPUT="diff -s ${CASE}_expected_result.vtu ${CASE}_with_results.vtu"
TESTER=${COMPARE_VTU_OUTPUT}

#
# Check results
#
${TESTER}
