use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# this parser is for the permutation data output files

my$directory = $ARGV[0] ;

opendir(FILES, $directory) ;
my@files = nsort (readdir(FILES)) ;

my%data ;

foreach my$i (0..$#files) {
  open IN, "<$directory/$files[$i]" or die "cannot open $directory/$files[$i]\n" ;

  my$next_line_id = "" ;

  while (<IN>) {
    chomp ;

    if ($_ =~ m/^perc_hmv_1evlvra(\d+)/) {$data{"perc_hmv_1evlv"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^aquaticra(\d+)/) {$data{"aquatic"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^terrestrialra(\d+)/) {$data{"terrestrial"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^routera(\d+)/) {$data{"route"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^nutritionra(\d+)/) {$data{"nutrition"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^defensera(\d+)/) {$data{"defense"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^multicomponentra(\d+)/) {$data{"multi"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^unknownra(\d+)/) {$data{"unknown"}{"$1_$i"} = $_ ;}
    if ($_ =~ m/^manipulativera(\d+)/) {$data{"manipulative"}{"$1_$i"} = $_ ;}
  }

  close IN ;
}

my$basename = basename($directory) ;
open OUT, ">consolidated_${basename}.txt" or die "cannot open consolidated_${basename}.txt\n" ;

foreach my$i (nsort keys %data) {
  print OUT "##", $i, "\n" ;

  foreach my$k (nsort keys %{$data{$i}}) {
    print OUT "n$k", $data{$i}{$k}, "\n" ;
  }
}

close OUT ;
