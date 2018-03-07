my $tag=$ARGV[0];

my $lineCounter=1;
my $lastCut=0;
while(my $line = <STDIN>){
    if(($lineCounter%2)==1) {
        print $line;
    } else {
        if(($lineCounter%4)==2) {
            my $tail = substr($line, -5, 4);
            if($tail eq $tag) {
                $lastCut=4;
                print substr($line, 0, -5) . "\n";
            } else {
                print $line;
            }
        } else {
            if($lastCut==4) {
                print substr($line, 0, -5) . "\n";
            } else {
                print $line;
            }
            $lastCut=0;
        }
    }
    ++$lineCounter;
}
