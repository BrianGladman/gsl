#!/bin/sh -e -x
NEW=$1
OLD=$2
NEWTAR="$NEW.tar.gz"
OLDTAR="$OLD.tar.gz"
BASE=`echo $NEW | sed -e 's/-.*//'`
NEWVER=`echo $NEW | sed -e 's/.*-//; s/\.tar.*//' `
OLDVER=`echo $OLD | sed -e 's/.*-//; s/\.tar.*//' `

if [ ! -e $NEWTAR ] ; then echo $NEWTAR not found ; exit ; fi
if [ ! -e $OLDTAR ] ; then echo $OLDTAR not found ; exit ; fi

rm -rf $OLD $NEW
tar xfz $OLDTAR
tar xfz $NEWTAR
PATCH="$BASE-$OLDVER-$NEWVER.patch.gz"
diff -C 2 -rcP -x '*.info' -x '*.info-*' $OLD $NEW | gzip -9 > $PATCH
(cd $OLD; gunzip -c -d ../$PATCH | patch -p1 ; diff -r -q . ../$NEW; )
rm -rf $OLD $NEW

# make directives
cat > $NEWTAR.directive <<EOF
version: 1.1
directory: gsl
filename: $NEWTAR
comment: new release of GNU Scientific Library
EOF

cat > $PATCH.directive <<EOF 
version: 1.1
directory: gsl
filename: $PATCH
comment: diffs for new release of GNU Scientific Library
EOF

# do signing
for f in  $NEWTAR $PATCH ; do
    gpg --sign --detach --yes $f
    gpg --clearsign $f.directive
    FILES="$FILES $f $f.sig $f.directive.asc"
done

read -p "Ready to upload? $FILES"

for n in $FILES; do
 echo curl -T -v $n ftp://ftp-upload.gnu.org/incoming/ftp/
done
