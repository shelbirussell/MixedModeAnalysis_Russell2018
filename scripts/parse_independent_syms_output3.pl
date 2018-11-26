use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# this parser is for the permutation pvalue output files

my$directory = $ARGV[0] ;

opendir(FILES, $directory) ;
my@files = nsort (readdir(FILES)) ;

my%data ;

foreach my$i (0..$#files) {
  open IN, "<$directory/$files[$i]" or die "cannot open $directory/$files[$i]\n" ;
  my$id = $files[$i] ;
  $id =~ s/^independent(\d+)_permuted_transmission_pvalues.txt/$1/ ;

  while (<IN>) {
    chomp ;
    my$value = $_ ;
    my$empirical = $_ ;
    $empirical =~ s/^.+of (\d+)$/$1/ ;

    ##Testing: abundance of horizontal transmission in aquatic environments relative to chance?
    if ($_ =~ m/^aquatic horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^aquatic horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_horizontal"}{$id}{"above"} = $value ;
        $data{"aquatic_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^aquatic horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^aquatic horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_horizontal"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of vertical transmission in aquatic environments relative to chance?
    elsif ($_ =~ m/^aquatic vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^aquatic vertical mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_vertical"}{$id}{"above"} = $value ;
        $data{"aquatic_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^aquatic vertical mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^aquatic vertical mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_vertical"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: is there a bias in mixed mode transmission in aquatic environments relative to chance?
    elsif ($_ =~ m/^aquatic mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^aquatic mixed mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_mixed"}{$id}{"above"} = $value ;
        $data{"aquatic_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^aquatic mixed mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^aquatic mixed mode:	(\d+)\(.+/$1/ ;
        $data{"aquatic_mixed"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of horizontal transmission in terrestrial environments relative to chance?
    elsif ($_ =~ m/^terrestrial horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^terrestrial horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_horizontal"}{$id}{"above"} = $value ;
        $data{"terrestrial_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^terrestrial horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^terrestrial horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_horizontal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of vertical transmission in terrestrial environments relative to chance?
    elsif ($_ =~ m/^terrestrial vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^terrestrial vertical mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_vertical"}{$id}{"above"} = $value ;
        $data{"terrestrial_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^terrestrial vertical mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^terrestrial vertical mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_vertical"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: bias in mixed mode transmission in terrestrial environments relative to chance?
    elsif ($_ =~ m/^terrestrial mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^terrestrial mixed mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_mixed"}{$id}{"above"} = $value ;
        $data{"terrestrial_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^terrestrial mixed mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^terrestrial mixed mode:	(\d+)\(.+/$1/ ;
        $data{"terrestrial_mixed"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of internal routes in mixed mode transmission relative to chance?
    elsif ($_ =~ m/^mixed internal route:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^mixed internal route:	(\d+)\(.+/$1/ ;
        $data{"mixed_internal"}{$id}{"above"} = $value ;
        $data{"mixed_internal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^mixed internal route:	(\d+)\(.+/$1/ ;
        $data{"mixed_internal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^mixed internal route:	(\d+)\(.+/$1/ ;
        $data{"mixed_internal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of vertical transmission in terrestrial environments relative to chance?
    elsif ($_ =~ m/^vertical internal route:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^vertical internal route:	(\d+)\(.+/$1/ ;
        $data{"vertical_internal"}{$id}{"above"} = $value ;
        $data{"vertical_internal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^vertical internal route:	(\d+)\(.+/$1/ ;
        $data{"vertical_internal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^vertical internal route:	(\d+)\(.+/$1/ ;
        $data{"vertical_internal"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of horizontal transmission in nutritive symbioses relative to chance?
    elsif ($_ =~ m/^nutrition functions horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^nutrition functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_horizontal"}{$id}{"above"} = $value ;
        $data{"nutrition_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^nutrition functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^nutrition functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_horizontal"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of mixed mode transmission in nutritive symbioses relative to chance?
    elsif ($_ =~ m/^nutrition functions mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^nutrition functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_mixed"}{$id}{"above"} = $value ;
        $data{"nutrition_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^nutrition functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^nutrition functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_mixed"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of vertical transmission in nutritive symbioses relative to chance?
    elsif ($_ =~ m/^nutrition functions vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^nutrition functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_vertical"}{$id}{"above"} = $value ;
        $data{"nutrition_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^nutrition functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^nutrition functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"nutrition_vertical"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of horizontal transmission in defensive symbioses relative to chance?
    elsif ($_ =~ m/^defense functions horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^defense functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"defense_horizontal"}{$id}{"above"} = $value ;
        $data{"defense_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^defense functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"defense_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^defense functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"defense_horizontal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of mixed mode transmission in defensive symbioses relative to chance?
    elsif ($_ =~ m/^defense functions mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^defense functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"defense_mixed"}{$id}{"above"} = $value ;
        $data{"defense_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^defense functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"defense_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^defense functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"defense_mixed"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: bias in vertical transmission in defensive symbioses relative to chance?
    elsif ($_ =~ m/^defense functions vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^defense functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"defense_vertical"}{$id}{"above"} = $value ;
        $data{"defense_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^defense functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"defense_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^defense functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"defense_vertical"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of horizontal transmission in manipulative symbioses relative to chance?
    elsif ($_ =~ m/^manipulative functions horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^manipulative functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_horizontal"}{$id}{"above"} = $value ;
        $data{"manipulative_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^manipulative functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^manipulative functions horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_horizontal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of mixed transmission in manipulative symbioses relative to chance?
    elsif ($_ =~ m/^manipulative functions mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^manipulative functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_mixed"}{$id}{"above"} = $value ;
        $data{"manipulative_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^manipulative functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^manipulative functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_mixed"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of vertical transmission in manipulative symbioses relative to chance?
    elsif ($_ =~ m/^manipulative functions vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^manipulative functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_vertical"}{$id}{"above"} = $value ;
        $data{"manipulative_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^manipulative functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^manipulative functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"manipulative_vertical"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of horizontal transmission in multiple function symbioses relative to chance?
    elsif ($_ =~ m/^multiple-function horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^multiple-function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"multi_horizontal"}{$id}{"above"} = $value ;
        $data{"multi_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^multiple-function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"multi_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^multiple-function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"multi_horizontal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of mixed transmission in multiple function symbioses relative to chance?
    elsif ($_ =~ m/^multiple-functions mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^multiple-functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"multi_mixed"}{$id}{"above"} = $value ;
        $data{"multi_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^multiple-functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"multi_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^multiple-functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"multi_mixed"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of vertical transmission in multiple function symbioses relative to chance?
    elsif ($_ =~ m/^multiple-functions vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^multiple-functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"multi_vertical"}{$id}{"above"} = $value ;
        $data{"multi_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^multiple-functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"multi_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^multiple-functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"multi_vertical"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: lack of horizontal transmission in unknown function symbioses relative to chance?
    elsif ($_ =~ m/^unknown function horizontal mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^unknown function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_horizontal"}{$id}{"above"} = $value ;
        $data{"unknown_horizontal"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^unknown function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_horizontal"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^unknown function horizontal mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_horizontal"}{$id}{"below"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: bias in mixed transmission unknown function symbioses relative to chance?
    elsif ($_ =~ m/^unknown functions mixed mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^unknown functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_mixed"}{$id}{"above"} = $value ;
        $data{"unknown_mixed"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^unknown functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_mixed"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^unknown functions mixed mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_mixed"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }

    ##Testing: abundance of vertical transmission in unknown function symbioses relative to chance?
    elsif ($_ =~ m/^unknown functions vertical mode:	(.+)/) {
      if ($1 =~ m/above/) {
        $value =~ s/^unknown functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_vertical"}{$id}{"above"} = $value ;
        $data{"unknown_vertical"}{$id}{"empirical"} = $empirical ;
      }
      elsif ($1 =~ m/below/) {
        $value =~ s/^unknown functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_vertical"}{$id}{"below"} += $value ;
      }
      elsif ($1 =~ m/equal/) {
        $value =~ s/^unknown functions vertical mode:	(\d+)\(.+/$1/ ;
        $data{"unknown_vertical"}{$id}{"above"} += $value ;
      }
      else {print "what doesn't match?\n", $_, "\n" ;}
    }



  }

  close IN ;
}

my$basename = basename($directory) ;
open OUT, ">consolidated_${basename}.txt" or die "cannot open consolidated_${basename}.txt\n" ;

## print out set of pvalues for each bin and compute average
## foreach character bin
foreach my$i (nsort keys %data) {
  print OUT "##", $i, "\n" ;
  my$above_sum = 0 ;
  my$below_sum = 0 ;
  my$permuted = 0 ;
  ## tally each permuted set of independent associations (k)
  foreach my$k (nsort keys %{$data{$i}}) {
    if (($data{$i}{$k}{"above"} + $data{$i}{$k}{"below"}) > 0) {
      $permuted += 100 ;
      print OUT $k, "\tempirical:", $data{$i}{$k}{"empirical"}, "\t" ;
      print OUT "above:", $data{$i}{$k}{"above"}, "\t" ;
      $above_sum += $data{$i}{$k}{"above"} ;
      print OUT "below:", $data{$i}{$k}{"below"}, "\n" ;
      $below_sum += $data{$i}{$k}{"below"} ;
    }
  }

  print OUT "proportion permutations above empirical value: ", $above_sum/$permuted, "(${above_sum} / ${permuted})", "\n" ;
  print OUT "proportion permutations below empirical value: ", $below_sum/$permuted, "(${below_sum} / ${permuted})", "\n\n" ;
}

close OUT ;
