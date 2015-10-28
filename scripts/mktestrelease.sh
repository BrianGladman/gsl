#!/bin/sh -e -x

if [ -z "$1" ]; then
  echo "Usage: $0 <tarfile>"
  exit
fi

NEWTAR="$1"
if [ ! -e $NEWTAR ] ; then echo $NEWTAR not found ; exit ; fi

NEWVER=$(echo $NEWTAR | sed -e 's/.*-//; s/\.tar.*//')

# make directives
DIRFILE="$NEWTAR.directive"
echo "generating $DIRFILE..."
cat > $DIRFILE <<EOF
version: 1.2
directory: gsl
filename: $NEWTAR
comment: test release of GNU Scientific Library
EOF

# do signing
echo "signing $NEWTAR..."
gpg -b $NEWTAR
echo "signing $DIRFILE..."
gpg --clearsign $DIRFILE

FILES="$NEWTAR $NEWTAR.sig $NEWTAR.directive.asc"
echo "Execute the following commands:"

for n in $FILES; do
 echo ncftpput ftp-upload.gnu.org /incoming/alpha $n
done
