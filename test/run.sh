#!/bin/sh

VALGRIND="valgrind --leak-check=full --track-origins=yes"

cd ../build/test/
${VALGRIND} ./check
for exe in ./ut*.exe; do
    echo
    echo "Running test ${exe}..."
    ${VALGRIND} ${exe} || exit 1
    echo
done
