#! /bin/sh

cat > test.exp.1.tmp <<EOF
# mean = 2.75
# sigma = 1.08972
1 2 1
2 3 2
3 4 0
4 5 1
EOF

echo 1 2 2.5 4 | ./gsl-histogram 1 5 4 > test.obs.1.tmp

cmp test.exp.1.tmp test.obs.1.tmp
STATUS=$?
rm test.exp.1.tmp test.obs.1.tmp

exit $STATUS