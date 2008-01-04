use Data::Dumper;
use XML::Parser;

# Get the bug database from Savannah.gnu.org as BUGS.xml
# Then run
#
#    perl BUGS.pl | cat -s > BUGS
#
# to generate the BUGS file

my $p = XML::Parser->new(Style => 'Stream', Pkg => 'MySubs');
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

        @keys = ('item_id', 'open_closed', 'status', 'category', 'summary');
        for $k (@keys) {
             printf("%s: %s\n", uc($k), $item->{$k});
        }
        print "\n";


        print $item->{'original_submission'}, "\n";

        print $item->{'old_value'}, "\n";


        print "\n";
    }
  
}
