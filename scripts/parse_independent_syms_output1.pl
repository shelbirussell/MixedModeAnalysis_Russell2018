use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

## this parser is for the transmission mode value output files

my$directory = $ARGV[0] ;

opendir(FILES, $directory) ;
my@files = nsort (readdir(FILES)) ;

my%data ;

foreach my$i (0..$#files) {
  open IN, "<$directory/$files[$i]" or die "cannot open $directory/$files[$i]\n" ;

  while (<IN>) {
    chomp ;

    ##
    if ($_ =~ m/^##/) {next; ;}

    if ($_ =~ m/^perc_hmv_1evlv = c\(/) {
      $data{"rates"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^aquatic = c\(/) {
      $data{"environments_aquatic"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^terrestrial = c\(/) {
      $data{"environments_terrestrial"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^route = c\(/) {
      $data{"routes"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^aquatic_route = c\(/) {
      $data{"routes_aquatic"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^terrestrial_route = c\(/) {
      $data{"routes_terrestrial"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^nutrition = c\(/) {
      $data{"functions_nutrition"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^defense = c\(/) {
      $data{"functions_defense"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^multicomponent = c\(/) {
      $data{"functions_multi"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^unknown = c\(/) {
      $data{"functions_unknown"}{$i} = $_ ;
    }

    elsif ($_ =~ m/^manipulative = c\(/) {
      $data{"functions_manipulative"}{$i} = $_ ;
    }

    else {
      if (length($_) == 0) {
          next ;
      }
      else {
        print "what doesn't match?\n", $_, "\n" ;
      }
    }
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
