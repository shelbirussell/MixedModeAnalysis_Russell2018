use strict ;
use warnings ;
use Sort::Naturally ;
use List::Util 1.33 'any';

## Read in data
my%data ;
my$count = 0 ;

##modes: 0=horizontal, 1=mixed, 2=vertical
##routes: 0=na, 1=external vertical 2=internal vertical
##environments: 0=aquatic, 1=terresterial, 2=freshwater
##functions: 0=nutrition, 1=chemosynthesis, 2=defense, 3,4=multi_component, 5=unknown, 6=manipulation
##functions: 0=nutrition, 1=nutrition+energy, 2=defense, 3=nutrition+defense, 4=energy+nutrition+defense, 5=unknown, 6=manipulation

my%variables ;

while (<STDIN>) {
  chomp ;
  if ($_ =~ m/^#/) { next ;  }
  if ($_ =~ m/^\s+$/) { next ;}
  my@split = split(/\t/, $_) ;

  $count ++ ;

  @{$data{$count}} = () ;

  push @{$data{$count}}, $split[0], $split[1], $split[2], $split[3]  ;

  ##host class
  $variables{"host_class:".$split[0]} = 1 ;
  ##host family
  $variables{"host_family:".$split[1]} = 1 ;
  ##symbiont phylum
  $variables{"sym_phylum:".$split[2]} = 1 ;
  ## symbiont family
  $variables{"sym_family:".$split[3]} = 1 ;
  if ($split[4] == 0) { $variables{"horizontal"} = 1 ; push @{$data{$count}}, "horizontal" ; }
  if ($split[4] == 1) { $variables{"mixed"} = 1 ; push @{$data{$count}}, "mixed" ; }
  if ($split[4] == 2) { $variables{"vertical"} = 1 ; push @{$data{$count}}, "vertical"; }

  if ($split[5] == 0) { $variables{"na"} = 1 ; push @{$data{$count}}, "na" ; }
  if ($split[5] == 1) { $variables{"external"} = 1 ; push @{$data{$count}}, "external" ; }
  if ($split[5] == 2) { $variables{"internal"} = 1 ; push @{$data{$count}}, "internal" ; }

  if ($split[6] == 0 || $split[6] == 2) { $variables{"aquatic"} = 1 ; push @{$data{$count}}, "aquatic" ; }
  if ($split[6] == 1) { $variables{"terrestrial"} = 1 ; push @{$data{$count}}, "terrestrial" ; }

  if ($split[7] == 0) { $variables{"nutrition"} = 1 ; push @{$data{$count}}, "nutrition" ; }
  if ($split[7] == 1) { $variables{"defense"} = 1 ; push @{$data{$count}}, "defense" ; }
  if ($split[7] == 2) { $variables{"multi_component"} = 1 ; push @{$data{$count}}, "multi_component" ; }
  if ($split[7] == 3) { $variables{"unknown"} = 1 ; push @{$data{$count}}, "unknown" ; }
  if ($split[7] == 4) { $variables{"manipulation"} = 1 ; push @{$data{$count}}, "manipulation" ; }
}


## Bin data
## Things to bin...
## % mixed
## % vertical
## % horizontal
## %s for environment, function, mode, route, host taxon, sym taxon
## permutation %s for subsamples of all values
## correlation coefficients
my@variables = (nsort keys %variables) ;
my@characters = ("horizontal", "mixed", "vertical", "na", "external", "internal", "aquatic", "terrestrial", "nutrition", "defense", "multi_component","unknown", "manipulation") ;

##Make lists of taxa
my%HOL ;
@{$HOL{"sym_phyla"}} = () ;
@{$HOL{"sym_families"}} = () ;
@{$HOL{"host_classes"}} = () ;
@{$HOL{"host_families"}} = () ;

foreach my$i (@variables) {
  if (any {/$i/} @characters) {next ;}
  if ($i =~ m/^host_class:(.+)/) {push @{$HOL{"host_classes"}}, $1 ;}
  if ($i =~ m/^host_family:(.+)/) {push @{$HOL{"host_families"}}, $1 ;}
  if ($i =~ m/^sym_phylum:(.+)/) {push @{$HOL{"sym_phyla"}}, $1 ;}
  if ($i =~ m/^sym_family:(.+)/) {push @{$HOL{"sym_families"}}, $1 ;}
}

#### Make contingency table for each taxonomic level
foreach my$division (keys %HOL) {
  open OUT, ">${division}_table.txt" or die "cannot open ${division}_table.txt\n" ;
  my@new_variable_list = () ;
  push @new_variable_list, @characters ;
  push @new_variable_list, @{$HOL{$division}} ;
  print OUT join(",", @new_variable_list), "\n" ;

  foreach my$new_variable_list1 (@new_variable_list) {
    foreach my$new_variable_list2 (@new_variable_list) {
      if ($new_variable_list1 eq $new_variable_list2 ) {

        ##count total number with this value
        my$c = 0 ;
        foreach my$i (keys %data) {
          my$matches = 0 ;
          foreach my$j (@{$data{$i}}) {
            if ($j eq $new_variable_list1) {
              $matches ++ ;
            }
          }
          ## if both characters were in list, add to count
          if ($matches == 1) {
            $c ++ ;
          }
        }
        print OUT "n=", $c, "\t" ;
      }

      else {
        my$c = 0 ;

        ##iterate through data to count number with new_variable_list1 and new_variable_list2
        foreach my$i (keys %data) {
          my$matches = 0 ;
          foreach my$j (@{$data{$i}}) {
            if ($j eq $new_variable_list1 || $j eq $new_variable_list2) {
              $matches ++ ;
            }
          }

          ## if both characters were in list, add to count
          if ($matches == 2) {
            $c ++ ;
          }
        }
        print OUT $c, "\t" ;
      }
    }
    print OUT "\n" ;
  }
  close OUT ;
}



##### Make contingency table - ALL BY ALL
#print join("\t", @variables), "\n" ;

#foreach my$variable1 (@variables) {
#  foreach my$variable2 (@variables) {
#    if ($variable1 eq $variable2 ) {
#
#      ##count total number with this value
#      my$c = 0 ;
#      foreach my$i (keys %data) {
#        my$matches = 0 ;
#        foreach my$j (@{$data{$i}}) {
#          if ($j eq $variable1) {
#            $matches ++ ;
#          }
#        }
#        ## if both variable was in list, add to count
#        if ($matches == 1) {
#          $c ++ ;
#        }
#      }
#      print "n=", $c, "\t" ;
#    }
#
#    else {
#      my$c = 0 ;
#
#      ##iterate through data to count number with variable1 and variable2
#      foreach my$i (keys %data) {
#        my$matches = 0 ;
#        foreach my$j (@{$data{$i}}) {
#          if ($j eq $variable1 || $j eq $variable2) {
#            $matches ++ ;
#          }
#        }
#
#        ## if both variables were in list, add to count
#        if ($matches == 2) {
#          $c ++ ;
#        }
#      }
#      print $c, "\t" ;
#
#    }
#  }
#  print "\n" ;
#}
