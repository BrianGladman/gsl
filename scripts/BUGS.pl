use Data::Dumper;
use XML::Parser;

# Get the bug database from Savannah.gnu.org as BUGS.xml
# Then run
#
#   wget --no-check-certificate https://savannah.gnu.org/export/gsl/bjg/28.xml -O BUGS.savannah.xml
#   recode latin1:
#     a) utf8 < BUGS.savannah.xml  > BUGS.xml
#     b) iconv -t utf-8 BUGS.savannah.xml > BUGS.xml
#   perl scripts/BUGS.pl | perl -p -e 's/^[ \t]+$//' | cat -s > BUGS
#
# to generate the BUGS file
binmode STDOUT, ":utf8";

my $date = scalar(localtime);
print <<"EOF";
The GSL Bugs Database is at http://savannah.gnu.org/bugs/?group=gsl

This file was generated from it at $date

EOF

my $p = XML::Parser->new(Style => 'Stream', Pkg => 'MySubs', ProtocolEncoding => 'UTF-8');
my $t = $p->parsefile('BUGS.xml');


print "-" x 72, "\n";

{
    package MySubs;
    my $item;
    
    sub StartTag {
        my ($e, $name) = @_;
        $item = {} if ($name eq 'item') ;
        #print "name : $name\n";
        $key = $name;
    }
  
    sub EndTag {
        my ($e, $name) = @_;
        # do something with end tags
        Format($item) if $name eq 'item';
    }
    
    sub Characters {
         my ($e, $data) = @_;
         # do something with text nodes
         #print "key = $key\n";
         $item->{$key} .= $data;
    }
    
    sub Text {
        # do something with text nodes
        s/^\s+//;
        s/&gt;/>/g;
        s/&lt;/</g;
        s/&amp;/&/g;
        s/&quot;/\"/g;
        #print "key = $key\n";
        $item->{$key} .= $_;
    }
    
    sub PI {
        return undef;
    }
    
    sub Format {
        print "-" x 72, "\n";

        $status = $item->{'status'};

        $item->{'status'} = $item->{'open_closed'} . 
               ($status eq 'None' ? "" : "/$status");

        printf("BUG-ID:   %-8s\nSTATUS:   %-16s\nCATEGORY: %s\n", 
               $item->{'item_id'}, $item->{'status'}, $item->{'category'});
        printf("SUMMARY:  %s\n", $item->{'summary'});
        print "\n";


        print $item->{'original_submission'}, "\n";

        print $item->{'old_value'}, "\n";


        print "\n";
    }
  
}
