my $lineCount=0;

$sum=0;
$bases=0;
while(my $line=<STDIN>) {
    ++$lineCount;
    if(($lineCount%4)==0) {
        chomp($line);
        @array=unpack("C*", $line);
        foreach my $elem (@array) {
            $sum += $elem - 33;
        }
        $bases+=scalar(@array)
    }
}

printf("%2.2f", $sum/$bases);
