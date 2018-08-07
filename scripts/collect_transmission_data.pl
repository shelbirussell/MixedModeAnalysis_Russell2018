use strict ;
use warnings ;
use Sort::Naturally ;

##usage: perl permute_collect_transmission_data.pl dataset.txt

my$symbioses = $ARGV[0] ;

open IN, "<$symbioses" or die "cannot open $symbioses\n" ;

my@column_names ;
my%symbioses ;
my%symbioses2lev ;
my%symbioses3lev ;


while (<IN>) {
  chomp ;
  ## get column labels
  ## #0=host_class	1=host_family	2=symbiont_phylum	3=symbiont_family	4=mode	5=route	6=environment	7=function  8=evidence
  if ($_ =~ /^#/) {
    $_ =~ s/^#// ;

    @column_names = split(/\t/, $_) ;

    ## make hash of columns, each with an array for its rows
    ## therefore, the relationships between data points in a row will be maintained by their order across arrays in the hash
    foreach my$i (0..$#column_names) {
      @{$symbioses{$i}} = () ;
      @{$symbioses2lev{$i}} = () ;
      @{$symbioses3lev{$i}} = () ;
    }
  }

  else {
    my@split = split(/\t/, $_) ;
    foreach my$i (0..$#split) {push @{$symbioses{$i}}, $split[$i] ;}

    ## only keep data with min 2 evidence levels
    if ($split[8] =~ m/[\d,]{2}/) {
      foreach my$i (0..$#split) {push @{$symbioses2lev{$i}}, $split[$i] ;}
    }

    ## only keep data with all 3 evidence levels
    if ($split[8] =~ m/\d,\d,\d/) {
      foreach my$i (0..$#split) {push @{$symbioses3lev{$i}}, $split[$i] ;}
    }
  }
}

close IN ;

open MAINOUT, ">transmission_analysis_output.txt" or die "cannot open transmission_analysis_output.txt\n" ;

## all data
my($horizontal, $mixed, $vertical) = get_rates(\%symbioses) ;
my($hterrestrial, $mterrestrial, $vterrestrial, $haquatic, $maquatic, $vaquatic) = get_env(\%symbioses) ;
my@env_variables = (4,6) ;
print_contingency_table(\@env_variables, \%symbioses) ;

my($externalvertical, $internalvertical, $externalmixed, $internalmixed) = get_routes(\%symbioses) ;
my@route_variables = (4,5) ;
print_contingency_table(\@route_variables, \%symbioses) ;

my($horiznut, $mixnut, $vertnut, $horizdef, $mixdef, $vertdef, $horizmult, $mixmult, $vertmult, $horizunk, $mixunk, $vertunk, $horizmanip, $mixmanip, $vertmanip) = get_function(\%symbioses) ;
my@func_variables = (4,7) ;
print_contingency_table(\@func_variables, \%symbioses) ;

my@host_class_variables = (4,0) ;
print_contingency_table(\@host_class_variables, \%symbioses) ;
my@host_family_variables = (4,1) ;
print_contingency_table(\@host_family_variables, \%symbioses) ;
my@symbiont_phylum_variables = (4,2) ;
print_contingency_table(\@symbiont_phylum_variables, \%symbioses) ;
my@symbiont_family_variables = (4,3) ;
print_contingency_table(\@symbiont_family_variables, \%symbioses) ;


## subsample by evidence levels
my($horizontal2, $mixed2, $vertical2) = get_rates(\%symbioses2lev) ;
my($horizontal3, $mixed3, $vertical3) = get_rates(\%symbioses3lev) ;
my($hterrestrial2, $mterrestrial2, $vterrestrial2, $haquatic2, $maquatic2, $vaquatic2) = get_env(\%symbioses2lev) ;
my($hterrestrial3, $mterrestrial3, $vterrestrial3, $haquatic3, $maquatic3, $vaquatic3) = get_env(\%symbioses3lev) ;
my($externalvertical2, $internalvertical2, $externalmixed2, $internalmixed2) = get_routes(\%symbioses2lev) ;
my($externalvertical3, $internalvertical3, $externalmixed3, $internalmixed3) = get_routes(\%symbioses3lev) ;
my($horiznut2, $mixnut2, $vertnut2, $horizdef2, $mixdef2, $vertdef2, $horizmult2, $mixmult2, $vertmult2, $horizunk2, $mixunk2, $vertunk2, $horizmanip2, $mixmanip2, $vertmanip2) = get_function(\%symbioses2lev) ;
my($horiznut3, $mixnut3, $vertnut3, $horizdef3, $mixdef3, $vertdef3, $horizmult3, $mixmult3, $vertmult3, $horizunk3, $mixunk3, $vertunk3, $horizmanip3, $mixmanip3, $vertmanip3) = get_function(\%symbioses3lev) ;


##### RATES
#print "# associations horizontally transmitted: ", $horizontal, " (", ($horizontal/($vertical+$horizontal+$mixed)*100), "%)\n" ;
#print "# associations with mixed transmission: ", $mixed, " (", ($mixed/($vertical+$horizontal+$mixed)*100), "%)\n" ;
#print "# associations vertically transmitted: ", $vertical, " (", ($vertical/($vertical+$horizontal+$mixed)*100), "%)\n" ;
print "total: ", ($vertical+$horizontal+$mixed),"\n" ;
print "total min 2 evidence levels: ", ($vertical2+$horizontal2+$mixed2),"\n" ;
print "total min 3 evidence levels: ", ($vertical3+$horizontal3+$mixed3),"\n" ;

## Print stats to file for plotting in R
print MAINOUT "##RATES\n" ;
my@variables = ($horizontal, $mixed, $vertical) ;
print_to_open_file (*MAINOUT, "perc_hmv_1evlv", \@variables) ;
@variables = ($horizontal2, $mixed2, $vertical2) ;
print_to_open_file (*MAINOUT, "perc_hmv_2evlv", \@variables) ;
@variables = ($horizontal3, $mixed3, $vertical3) ;
print_to_open_file (*MAINOUT, "perc_hmv_3evlv", \@variables) ;

print MAINOUT "\n##ENVIRONMENTS\n" ;
@variables = ($haquatic,$maquatic,$vaquatic) ;
print_to_open_file (*MAINOUT, "aquatic", \@variables) ;
@variables = ($haquatic2,$maquatic2,$vaquatic2) ;
print_to_open_file (*MAINOUT, "aquatic2", \@variables) ;
@variables = ($haquatic3,$maquatic3,$vaquatic3) ;
print_to_open_file (*MAINOUT, "aquatic3", \@variables) ;
@variables = ($hterrestrial,$mterrestrial,$vterrestrial) ;
print_to_open_file (*MAINOUT, "terrestrial", \@variables) ;
@variables = ($hterrestrial2,$mterrestrial2,$vterrestrial2) ;
print_to_open_file (*MAINOUT, "terrestrial2", \@variables) ;
@variables = ($hterrestrial3,$mterrestrial3,$vterrestrial3) ;
print_to_open_file (*MAINOUT, "terrestrial3", \@variables) ;

print MAINOUT "\n##ROUTES (internal)\n" ;
print $internalmixed, "\t", $internalvertical, "\n" ;
@variables = ($internalmixed,$externalmixed,$internalvertical,$externalvertical) ;
print_to_open_file2 (*MAINOUT, "route", \@variables) ;
@variables = ($internalmixed2,$externalmixed2,$internalvertical2,$externalvertical2) ;
print_to_open_file2 (*MAINOUT, "route2", \@variables) ;
@variables = ($internalmixed3,$externalmixed3,$internalvertical3,$externalvertical3) ;
print_to_open_file2 (*MAINOUT, "route3", \@variables) ;

print MAINOUT "\n##FUNCTIONS\n" ;
@variables = ($horiznut,$mixnut,$vertnut) ;
print_to_open_file (*MAINOUT, "nutrition", \@variables) ;
@variables = ($horizdef,$mixdef,$vertdef) ;
print_to_open_file (*MAINOUT, "defense", \@variables) ;
@variables = ($horizmult,$mixmult,$vertmult) ;
print_to_open_file (*MAINOUT, "multicomponent", \@variables) ;
@variables = ($horizunk,$mixunk,$vertunk) ;
print_to_open_file (*MAINOUT, "unknown", \@variables) ;
@variables = ($horizmanip,$mixmanip,$vertmanip) ;
print_to_open_file (*MAINOUT, "manipulative", \@variables) ;

@variables = ($horiznut2,$mixnut2,$vertnut2) ;
print_to_open_file (*MAINOUT, "nutrition2", \@variables) ;
@variables = ($horizdef2,$mixdef2,$vertdef2) ;
print_to_open_file (*MAINOUT, "defense2", \@variables) ;
@variables = ($horizmult2,$mixmult2,$vertmult2) ;
print_to_open_file (*MAINOUT, "multicomponent2", \@variables) ;
@variables = ($horizunk2,$mixunk2,$vertunk2) ;
print_to_open_file (*MAINOUT, "unknown2", \@variables) ;
@variables = ($horizmanip2,$mixmanip2,$vertmanip2) ;
print_to_open_file (*MAINOUT, "manipulative2", \@variables) ;

@variables = ($horiznut3,$mixnut3,$vertnut3) ;
print_to_open_file (*MAINOUT, "nutrition3", \@variables) ;
@variables = ($horizdef3,$mixdef3,$vertdef3) ;
print_to_open_file (*MAINOUT, "defense3", \@variables) ;
@variables = ($horizmult3,$mixmult3,$vertmult3) ;
print_to_open_file (*MAINOUT, "multicomponent3", \@variables) ;
@variables = ($horizunk3,$mixunk3,$vertunk3) ;
print_to_open_file (*MAINOUT, "unknown3", \@variables) ;
@variables = ($horizmanip3,$mixmanip3,$vertmanip3) ;
print_to_open_file (*MAINOUT, "manipulative3", \@variables) ;

close MAINOUT ;


sub print_to_open_file {
  my$handle = $_[0] ;
  my$arrayname = $_[1] ;
  my$variables = $_[2] ;
  my@variables = @{$variables} ;
  my$sum = 0 ;
  foreach my$v (@variables) {$sum += $v ;}
  if ($sum > 0) {
    print $handle $arrayname, " = c(" ;
    foreach my$i (0..($#variables-1)) {
      print $handle ($variables[$i]/$sum*100), "," ;
    }
    print $handle ($variables[$#variables]/$sum*100), ")\n" ;
  }
  else {#print $arrayname, " won't print\n" ;
  }
}

sub print_to_open_file2 {
  my$handle = $_[0] ;
  my$arrayname = $_[1] ;
  my$variables = $_[2] ;
  my@variables = @{$variables} ;
  my$sum1 = 0 ;
  my$sum2 = 0 ;
  foreach my$v (0..1) {$sum1 += $variables[$v] ;}
  foreach my$v (2..3) {$sum2 += $variables[$v] ;}
  if ($sum1 > 0 && $sum2 >0) {
    print $handle $arrayname, " = c(" ;
    print $handle ($variables[0]/$sum1*100), "," ;
    print $handle ($variables[2]/$sum2*100), ")\n" ;
  }
}

sub get_rates {
  my%data = %{$_[0]} ;

  ## RATES: get mode %
  my$H = 0 ;
  my$M = 0 ;
  my$V = 0 ;
  foreach my$i (@{$data{4}}) {
    if ($i == 0) {$H ++ ;}
    if ($i == 1) {$M ++ ;}
    if ($i == 2) {$V ++;}
  }

  return ($H, $M, $V) ;
}

sub get_env {
  my%data = %{$_[0]} ;

  ## ENV TYPE: get % of mode in each ENV
  ## unfiltered
  my$HT = 0 ;
  my$HA = 0 ;
  my$MT = 0 ;
  my$MA = 0 ;
  my$VT = 0 ;
  my$VA = 0 ;
  foreach my$i (0..$#{$data{4}}) {
    if ($data{6}[$i] == 0 || $data{6}[$i] == 2) {
      if ($data{4}[$i] == 0) {$HA ++ ;}
      if ($data{4}[$i] == 1) {$MA ++ ;}
      if ($data{4}[$i] == 2) {$VA ++;}
    }
    if ($data{6}[$i] == 1) {
      if ($data{4}[$i] == 0) {$HT ++ ;}
      if ($data{4}[$i] == 1) {$MT ++ ;}
      if ($data{4}[$i] == 2) {$VT ++;}
    }
  }

    return ($HT, $MT, $VT, $HA, $MA, $VA) ;
}

sub get_routes {
  my%data = %{$_[0]} ;

  ## ROUTE TYPE: get % of mode in each route
  ## unfiltered
  my$VE = 0 ;
  my$VI = 0 ;
  my$ME = 0 ;
  my$MI = 0 ;
  foreach my$i (0..$#{$data{4}}) {
    # if vertically transmitted or have mixed modes
    if ($data{4}[$i] == 1) {
      if ($data{5}[$i] == 1) {$ME ++ ;}
      if ($data{5}[$i] == 2) {$MI ++ ;}
    }

    elsif ($data{4}[$i] == 2) {
      if ($data{5}[$i] == 1) {$VE ++ ;}
      if ($data{5}[$i] == 2) {$VI ++ ;}
    }

    else {
      next ;
    }
  }

    return ($VE, $VI, $ME, $MI) ;
}

sub get_function {
  my%data = %{$_[0]} ;

  ## ENV TYPE: get % of mode in each ENV
  ## unfiltered
  my$HN = 0 ;
  my$HD = 0 ;
  my$HMu = 0 ;
  my$HUk = 0 ;
  my$HMan = 0 ;

  my$MN = 0 ;
  my$MD = 0 ;
  my$MMu = 0 ;
  my$MUk = 0 ;
  my$MMan = 0 ;

  my$VN = 0 ;
  my$VD = 0 ;
  my$VMu = 0 ;
  my$VUk = 0 ;
  my$VMan = 0 ;
  foreach my$i (0..$#{$data{4}}) {
    ##nutrition
    if ($data{7}[$i] == 0) {
      if ($data{4}[$i] == 0) {$HN ++ ;}
      if ($data{4}[$i] == 1) {$MN ++ ;}
      if ($data{4}[$i] == 2) {$VN ++;}
    }
    ##defense
    if ($data{7}[$i] == 1) {
      if ($data{4}[$i] == 0) {$HD ++ ;}
      if ($data{4}[$i] == 1) {$MD ++ ;}
      if ($data{4}[$i] == 2) {$VD ++;}
    }
    ##multi
    if ($data{7}[$i] == 2) {
      if ($data{4}[$i] == 0) {$HMu ++ ;}
      if ($data{4}[$i] == 1) {$MMu ++ ;}
      if ($data{4}[$i] == 2) {$VMu ++;}
    }
    ##unknown
    if ($data{7}[$i] == 3) {
      if ($data{4}[$i] == 0) {$HUk ++ ;}
      if ($data{4}[$i] == 1) {$MUk ++ ;}
      if ($data{4}[$i] == 2) {$VUk ++;}
    }
    ##manipulation
    if ($data{7}[$i] == 4) {
      if ($data{4}[$i] == 0) {$HMan ++ ;}
      if ($data{4}[$i] == 1) {$MMan ++ ;}
      if ($data{4}[$i] == 2) {$VMan ++;}
    }
  }

    return ($HN, $MN, $VN, $HD, $MD, $VD, $HMu, $MMu, $VMu, $HUk, $MUk, $VUk, $HMan, $MMan, $VMan) ;
}

sub print_contingency_table {
  ## variables are the column names
  ## 0=host_class	1=host_family	2=symbiont_phylum	3=symbiont_family	4=mode	5=route	6=environment	7=function  8=evidence
  my@variables = @{$_[0]} ;
  my%data = %{$_[1]} ;

  my%table ;

  ##iterate through all rows, tallying the # occurrances of variable_class1_variable_value1 + variable_class2_variable_value2
  foreach my$i (0..$#{$data{$variables[1]}}) {
    if (exists $table{$variables[1]."_".$data{$variables[1]}[$i]."\t".$variables[0]."_".$data{$variables[0]}[$i]}) {
      $table{$variables[1]."_".$data{$variables[1]}[$i]."\t".$variables[0]."_".$data{$variables[0]}[$i]} ++ ;
    }
    else {
      $table{$variables[1]."_".$data{$variables[1]}[$i]."\t".$variables[0]."_".$data{$variables[0]}[$i]} = 1 ;
    }
  }

  ## assign labels based on values
  ## environment
  if (grep { $_ == 6} @variables) {
    ##options are 6_0, 6_1, 6_2 paired with 4_0, 4_1, 4_2
    ## check for zeros
    foreach my$n1 (0..2) {
      foreach my$n2 (0..2) {
        if (exists $table{"6_".$n1."\t"."4_".$n2}) {next ;}
        else {$table{"6_".$n1."\t"."4_".$n2} = 0 ;}}}
    open OUT, "> environment_vs_mode_contingency_table.txt" or die "cannot open environment_vs_mode_contingency_table.txt" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;
    print OUT "aquatic\t", ($table{"6_0\t4_0"} + $table{"6_2\t4_0"}), "\t", ($table{"6_0\t4_1"} + $table{"6_2\t4_1"}), "\t", ($table{"6_0\t4_2"} + $table{"6_2\t4_2"}), "\n" ;
    print OUT "terrestrial\t", $table{"6_1\t4_0"}, "\t", $table{"6_1\t4_1"}, "\t", $table{"6_1\t4_2"} ;
    close OUT ;
  }
  ## route
  if (grep { $_ == 5} @variables) {
    ##options are 5_1, 5_2 paired with 4_1, 4_2
    ##skip 5_0 and 4_0, as these are horizontal and irrelevant to this measure
    open OUT, "> route_vs_mode_contingency_table.txt" or die "cannot open route_vs_mode_contingency_table.txt" ;
    print OUT "\tmixed\tvertical\n" ;
    print OUT "external\t", $table{"5_1\t4_1"}, "\t", $table{"5_1\t4_2"}, "\n" ;
    print OUT "internal\t", $table{"5_2\t4_1"}, "\t", $table{"5_2\t4_2"} ;
    close OUT ;
  }
  ## function
  if (grep { $_ == 7} @variables) {
    ##options are 7_0, 7_1, 7_2, 7_3, 7_4 paired with 4_0, 4_1, 4_2
    ## check for zeros
    foreach my$n1 (0..4) {
      foreach my$n2 (0..2) {
        if (exists $table{"7_".$n1."\t"."4_".$n2}) {next ;}
        else {$table{"7_".$n1."\t"."4_".$n2} = 0 ;}}}
    open OUT, "> function_vs_mode_contingency_table.txt" or die "cannot open function_vs_mode_contingency_table.txt" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;
    print OUT "nutrition\t", $table{"7_0\t4_0"}, "\t", $table{"7_0\t4_1"}, "\t", $table{"7_0\t4_2"}, "\n" ;
    print OUT "defense\t", $table{"7_1\t4_0"}, "\t", $table{"7_1\t4_1"}, "\t", $table{"7_1\t4_2"}, "\n" ;
    print OUT "multicomponent\t", $table{"7_2\t4_0"}, "\t", $table{"7_2\t4_1"}, "\t", $table{"7_2\t4_2"}, "\n" ;
    print OUT "unknown\t", $table{"7_3\t4_0"}, "\t", $table{"7_3\t4_1"}, "\t", $table{"7_3\t4_2"}, "\n" ;
    print OUT "manipulation\t", $table{"7_4\t4_0"}, "\t", $table{"7_4\t4_1"}, "\t", $table{"7_4\t4_2"} ;
    close OUT ;
  }
  ## taxa
  if ((grep { $_ == 0 } @variables) || (grep { $_ == 1 } @variables) || (grep { $_ == 2 } @variables) || (grep { $_ == 3 } @variables)) {
    ##check for emptys, and add zeros
    my%taxa ;
    foreach my$k (nsort keys %table) {
      my@split = split(/\t/, $k) ;
      my$taxon = $split[0] ;
      $taxon =~ s/\d_// ;
      $taxa{$taxon} = "" ;
    }

    ### sort and print counts
    my$outfile = "" ;
    if ($variables[1] == 0) {$outfile = "hostClass_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 1) {$outfile = "hostFamily_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 2) {$outfile = "symPhylum_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 3) {$outfile = "symFamily_vs_mode_contingency_table.txt" ;}
    open OUT, "> $outfile" or die "cannot open $outfile" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;

    foreach my$t (keys %taxa) {
      if (exists $table{$variables[1]."_".$t."\t4_0"}) {
        print OUT $t, "\t", $table{$variables[1]."_".$t."\t4_0"}, "\t" ;}
      else {print OUT $t, "\t", 0, "\t" ;}
      if (exists $table{$variables[1]."_".$t."\t4_1"}) {
        print OUT $table{$variables[1]."_".$t."\t4_1"}, "\t" ;}
      else {print OUT 0, "\t" ;}
      if (exists $table{$variables[1]."_".$t."\t4_2"}) {
        print OUT $table{$variables[1]."_".$t."\t4_2"}, "\n" ;      }
      else {print OUT 0, "\n" ;}
    }
    close OUT ;
  }
}
