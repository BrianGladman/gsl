#!/usr/bin/perl

print <<'EOF';
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<HTML>
<HEAD>
<meta name="GENERATOR" content="toc2htmlhelp">
<!-- Sitemap 1.0 -->
</HEAD><BODY>
<OBJECT type="text/site properties">
        <param name="Type" value="Reference">
        <param name="TypeDesc" value="Reference topics">
</OBJECT>
EOF

while (<>) {
    #print;
s#<A NAME="(\S+)" HREF="(\S+)(\#\S+)?">(.*)</A>#
    <OBJECT type="text/sitemap">
    <param name="Name" value="$4">
    <param name="Type" value="Reference">
    <param name="Local" value="$2">
    </OBJECT>#;
        print;

}

print <<'EOF';
</BODY></HTML>
EOF

