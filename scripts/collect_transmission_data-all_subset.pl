use strict ;
use warnings ;
use Sort::Naturally ;
use List::Util qw/shuffle/;

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
    if ($split[12] =~ m/[\d,]{2}/) {
      foreach my$i (0..$#split) {push @{$symbioses2lev{$i}}, $split[$i] ;}
    }

    ## only keep data with all 3 evidence levels
    if ($split[12] =~ m/\d,\d,\d/) {
      foreach my$i (0..$#split) {push @{$symbioses3lev{$i}}, $split[$i] ;}
    }
  }
}

close IN ;

open MAINOUT, ">transmission_analysis_output.txt" or die "cannot open transmission_analysis_output.txt\n" ;


my%empdata ;

## all data
my($horizontal, $mixed, $vertical) = get_rates(\%symbioses) ;
$empdata{"H"} = $horizontal ;
$empdata{"M"} = $mixed ;
$empdata{"V"} = $vertical;

my($hterrestrial, $mterrestrial, $vterrestrial, $haquatic, $maquatic, $vaquatic) = get_env(\%symbioses) ;
my@env_variables = (6,8) ;
print_contingency_table(\@env_variables, \%symbioses) ;
$empdata{"HA"} = $haquatic ;
$empdata{"MA"} = $maquatic ;
$empdata{"VA"} = $vaquatic ;
$empdata{"HT"} = $hterrestrial ;
$empdata{"MT"} = $mterrestrial ;
$empdata{"VT"} = $vterrestrial ;

my($externalvertical, $internalvertical, $externalmixed, $internalmixed) = get_routes(\%symbioses) ;
my@route_variables = (6,7) ;
print_contingency_table(\@route_variables, \%symbioses) ;
$empdata{"MI"} = $internalmixed ;
$empdata{"VI"} = $internalvertical;

my($horiznut, $mixnut, $vertnut, $horizdef, $mixdef, $vertdef, $horizmult, $mixmult, $vertmult, $horizunk, $mixunk, $vertunk, $horizmanip, $mixmanip, $vertmanip) = get_function(\%symbioses) ;
my@func_variables = (6,9) ;
print_contingency_table(\@func_variables, \%symbioses) ;
$empdata{"HN"} = $horiznut ;
$empdata{"HD"} = $horizdef ;
$empdata{"HMu"} = $horizmult ;
$empdata{"HUk"} = $horizunk ;
$empdata{"HMa"} = $horizmanip ;
$empdata{"MD"} = $mixdef ;
$empdata{"MN"} = $mixnut ;
$empdata{"MMu"} = $mixmult ;
$empdata{"MUk"} = $mixunk ;
$empdata{"MMa"} = $mixmanip ;
$empdata{"VN"} = $vertnut ;
$empdata{"VD"} = $vertdef ;
$empdata{"VMu"} = $vertmult ;
$empdata{"VUk"} = $vertunk ;
$empdata{"VMa"} = $vertmanip ;

my@host_class_variables = (6,2) ;
print_contingency_table(\@host_class_variables, \%symbioses) ;
my@host_family_variables = (6,3) ;
print_contingency_table(\@host_family_variables, \%symbioses) ;
my@symbiont_phylum_variables = (6,4) ;
print_contingency_table(\@symbiont_phylum_variables, \%symbioses) ;
my@symbiont_family_variables = (6,5) ;
print_contingency_table(\@symbiont_family_variables, \%symbioses) ;

### Permutations and subsampling
resample_data(\%symbioses, \%empdata, "all_permuted", "permute") ;
resample_data(\%symbioses, \%empdata, "all_subsampleEnv", "environment") ;
resample_data(\%symbioses, \%empdata, "all_subsamplePhyla", "hitaxa") ;
resample_data(\%symbioses, \%empdata, "all_subsampleFamilies", "lotaxa") ;


### subsample by evidence levels
my($horizontal2, $mixed2, $vertical2) = get_rates(\%symbioses2lev) ;
my($horizontal3, $mixed3, $vertical3) = get_rates(\%symbioses3lev) ;
my($hterrestrial2, $mterrestrial2, $vterrestrial2, $haquatic2, $maquatic2, $vaquatic2) = get_env(\%symbioses2lev) ;
my($hterrestrial3, $mterrestrial3, $vterrestrial3, $haquatic3, $maquatic3, $vaquatic3) = get_env(\%symbioses3lev) ;
my($externalvertical2, $internalvertical2, $externalmixed2, $internalmixed2) = get_routes(\%symbioses2lev) ;
my($externalvertical3, $internalvertical3, $externalmixed3, $internalmixed3) = get_routes(\%symbioses3lev) ;
my($horiznut2, $mixnut2, $vertnut2, $horizdef2, $mixdef2, $vertdef2, $horizmult2, $mixmult2, $vertmult2, $horizunk2, $mixunk2, $vertunk2, $horizmanip2, $mixmanip2, $vertmanip2) = get_function(\%symbioses2lev) ;
my($horiznut3, $mixnut3, $vertnut3, $horizdef3, $mixdef3, $vertdef3, $horizmult3, $mixmult3, $vertmult3, $horizunk3, $mixunk3, $vertunk3, $horizmanip3, $mixmanip3, $vertmanip3) = get_function(\%symbioses3lev) ;


##### RATES - print to file for plotting
#print "# associations horizontally transmitted: ", $horizontal, " (", ($horizontal/($vertical+$horizontal+$mixed)*100), "%)\n" ;
#print "# associations with mixed transmission: ", $mixed, " (", ($mixed/($vertical+$horizontal+$mixed)*100), "%)\n" ;
#print "# associations vertically transmitted: ", $vertical, " (", ($vertical/($vertical+$horizontal+$mixed)*100), "%)\n" ;
print "total: ", ($vertical+$horizontal+$mixed),"\n" ;
print "total min 2 evidence levels: ", ($vertical2+$horizontal2+$mixed2),"\n" ;
print "total min 3 evidence levels: ", ($vertical3+$horizontal3+$mixed3),"\n" ;

### Print stats to file for plotting in R
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
##print $internalmixed, "\t", $internalvertical, "\n" ;
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
  my@new_variables ;
  foreach my$v (@variables) {$sum += $v ;}
  if ($sum > 0) {
    print $handle $arrayname, " = c(" ;
    foreach my$i (0..($#variables-1)) {
      print $handle ($variables[$i]/$sum*100), "," ;
      push @new_variables, ($variables[$i]/$sum*100) ;
    }
    print $handle ($variables[$#variables]/$sum*100), ")\n" ;
    push @new_variables, ($variables[$#variables]/$sum*100) ;
  }
  else {#print $arrayname, " won't print\n" ;
  }
  return \@new_variables ;
}

sub print_to_open_file2 {
  my$handle = $_[0] ;
  my$arrayname = $_[1] ;
  my$variables = $_[2] ;
  my@variables = @{$variables} ;
  my$sum1 = 0 ;
  my$sum2 = 0 ;
  my@new_variables = () ;
  foreach my$v (0..1) {$sum1 += $variables[$v] ;}
  foreach my$v (2..3) {$sum2 += $variables[$v] ;}
  if ($sum1 > 0 && $sum2 >0) {
    print $handle $arrayname, " = c(" ;
    print $handle ($variables[0]/$sum1*100), "," ;
    print $handle ($variables[2]/$sum2*100), ")\n" ;
    push @new_variables, ($variables[0]/$sum1*100) ;
    push @new_variables, ($variables[2]/$sum2*100) ;
  }
  return \@new_variables ;
}

sub get_rates {
  my%data = %{$_[0]} ;

  ## RATES: get mode %
  my$H = 0 ;
  my$M = 0 ;
  my$V = 0 ;
  foreach my$i (@{$data{6}}) {
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
  foreach my$i (0..$#{$data{6}}) {
    if ($data{8}[$i] == 0 || $data{8}[$i] == 2) {
      if ($data{6}[$i] == 0) {$HA ++ ;}
      if ($data{6}[$i] == 1) {$MA ++ ;}
      if ($data{6}[$i] == 2) {$VA ++;}
    }
    if ($data{8}[$i] == 1) {
      if ($data{6}[$i] == 0) {$HT ++ ;}
      if ($data{6}[$i] == 1) {$MT ++ ;}
      if ($data{6}[$i] == 2) {$VT ++;}
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
  foreach my$i (0..$#{$data{6}}) {
    # if vertically transmitted or have mixed modes
    if ($data{6}[$i] == 1) {
      if ($data{7}[$i] == 1) {$ME ++ ;}
      if ($data{7}[$i] == 2) {$MI ++ ;}
    }

    elsif ($data{6}[$i] == 2) {
      if ($data{7}[$i] == 1) {$VE ++ ;}
      if ($data{7}[$i] == 2) {$VI ++ ;}
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
  foreach my$i (0..$#{$data{6}}) {
    ##nutrition
    if ($data{9}[$i] == 0) {
      if ($data{6}[$i] == 0) {$HN ++ ;}
      if ($data{6}[$i] == 1) {$MN ++ ;}
      if ($data{6}[$i] == 2) {$VN ++;}
    }
    ##defense
    if ($data{9}[$i] == 1) {
      if ($data{6}[$i] == 0) {$HD ++ ;}
      if ($data{6}[$i] == 1) {$MD ++ ;}
      if ($data{6}[$i] == 2) {$VD ++;}
    }
    ##multi
    if ($data{9}[$i] == 2) {
      if ($data{6}[$i] == 0) {$HMu ++ ;}
      if ($data{6}[$i] == 1) {$MMu ++ ;}
      if ($data{6}[$i] == 2) {$VMu ++;}
    }
    ##unknown
    if ($data{9}[$i] == 3) {
      if ($data{6}[$i] == 0) {$HUk ++ ;}
      if ($data{6}[$i] == 1) {$MUk ++ ;}
      if ($data{6}[$i] == 2) {$VUk ++;}
    }
    ##manipulation
    if ($data{9}[$i] == 4) {
      if ($data{6}[$i] == 0) {$HMan ++ ;}
      if ($data{6}[$i] == 1) {$MMan ++ ;}
      if ($data{6}[$i] == 2) {$VMan ++;}
    }
  }

    return ($HN, $MN, $VN, $HD, $MD, $VD, $HMu, $MMu, $VMu, $HUk, $MUk, $VUk, $HMan, $MMan, $VMan) ;
}

sub print_contingency_table {
  ## variables are the column names
  ## 2=host_class	3=host_family	4=symbiont_phylum	5=symbiont_family	6=mode	7=route	8=environment	9=function  12=evidence
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
  if (grep { $_ == 8} @variables) {
    ##options are 8_0, 8_1, 8_2 paired with 6_0, 6_1, 6_2
    ## check for zeros
    foreach my$n1 (0..2) {
      foreach my$n2 (0..2) {
        if (exists $table{"8_".$n1."\t"."6_".$n2}) {next ;}
        else {$table{"8_".$n1."\t"."6_".$n2} = 0 ;}}}
    open OUT, "> environment_vs_mode_contingency_table.txt" or die "cannot open environment_vs_mode_contingency_table.txt" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;
    print OUT "aquatic\t", ($table{"8_0\t6_0"} + $table{"8_2\t6_0"}), "\t", ($table{"8_0\t6_1"} + $table{"8_2\t6_1"}), "\t", ($table{"8_0\t6_2"} + $table{"8_2\t6_2"}), "\n" ;
    print OUT "terrestrial\t", $table{"8_1\t6_0"}, "\t", $table{"8_1\t6_1"}, "\t", $table{"8_1\t6_2"} ;
    close OUT ;

    open OUT, "> environment_vs_mode_contingency_table_fw.txt" or die "cannot open environment_vs_mode_contingency_table_fw.txt" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;
    print OUT "marine\t", $table{"8_0\t6_0"}, "\t", $table{"8_0\t6_1"}, "\t", $table{"8_0\t6_2"}, "\n" ;
    print OUT "freshwater\t", $table{"8_2\t6_0"}, "\t", $table{"8_2\t6_1"}, "\t", $table{"8_2\t6_2"}, "\n" ;
    print OUT "terrestrial\t", $table{"8_1\t6_0"}, "\t", $table{"8_1\t6_1"}, "\t", $table{"8_1\t6_2"} ;
    close OUT ;

  }
  ## route
  if (grep { $_ == 7} @variables) {
    ##options are 7_1, 7_2 paired with 6_1, 6_2
    ##skip 7_0 and 6_0, as these are horizontal and irrelevant to this measure
    open OUT, "> route_vs_mode_contingency_table.txt" or die "cannot open route_vs_mode_contingency_table.txt" ;
    print OUT "\tmixed\tvertical\n" ;
    print OUT "external\t", $table{"7_1\t6_1"}, "\t", $table{"7_1\t6_2"}, "\n" ;
    print OUT "internal\t", $table{"7_2\t6_1"}, "\t", $table{"7_2\t6_2"} ;
    close OUT ;
  }
  ## function
  if (grep { $_ == 9} @variables) {
    ##options are 9_0, 9_1, 9_2, 9_3, 9_4 paired with 6_0, 6_1, 6_2
    ## check for zeros
    foreach my$n1 (0..4) {
      foreach my$n2 (0..2) {
        if (exists $table{"9_".$n1."\t"."6_".$n2}) {next ;}
        else {$table{"9_".$n1."\t"."6_".$n2} = 0 ;}}}
    open OUT, "> function_vs_mode_contingency_table.txt" or die "cannot open function_vs_mode_contingency_table.txt" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;
    print OUT "nutrition\t", $table{"9_0\t6_0"}, "\t", $table{"9_0\t6_1"}, "\t", $table{"9_0\t6_2"}, "\n" ;
    print OUT "defense\t", $table{"9_1\t6_0"}, "\t", $table{"9_1\t6_1"}, "\t", $table{"9_1\t6_2"}, "\n" ;
    print OUT "multicomponent\t", $table{"9_2\t6_0"}, "\t", $table{"9_2\t6_1"}, "\t", $table{"9_2\t6_2"}, "\n" ;
    print OUT "unknown\t", $table{"9_3\t6_0"}, "\t", $table{"9_3\t6_1"}, "\t", $table{"9_3\t6_2"}, "\n" ;
    print OUT "manipulation\t", $table{"9_4\t6_0"}, "\t", $table{"9_4\t6_1"}, "\t", $table{"9_4\t6_2"} ;
    close OUT ;
  }
  ## taxa
  if ((grep { $_ == 2 } @variables) || (grep { $_ == 3 } @variables) || (grep { $_ == 4 } @variables) || (grep { $_ == 5 } @variables)) {
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
    if ($variables[1] == 2) {$outfile = "hostClass_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 3) {$outfile = "hostFamily_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 4) {$outfile = "symPhylum_vs_mode_contingency_table.txt" ;}
    if ($variables[1] == 5) {$outfile = "symFamily_vs_mode_contingency_table.txt" ;}
    open OUT, "> $outfile" or die "cannot open $outfile" ;
    print OUT "\thorizontal\tmixed\tvertical\n" ;

    foreach my$t (keys %taxa) {
      if (exists $table{$variables[1]."_".$t."\t6_0"}) {
        print OUT $t, "\t", $table{$variables[1]."_".$t."\t6_0"}, "\t" ;}
      else {print OUT $t, "\t", 0, "\t" ;}
      if (exists $table{$variables[1]."_".$t."\t6_1"}) {
        print OUT $table{$variables[1]."_".$t."\t6_1"}, "\t" ;}
      else {print OUT 0, "\t" ;}
      if (exists $table{$variables[1]."_".$t."\t6_2"}) {
        print OUT $table{$variables[1]."_".$t."\t6_2"}, "\n" ;      }
      else {print OUT 0, "\n" ;}
    }
    close OUT ;
  }
}

sub randomize_data {
  ## Shuffle the values of every characteristic (column), randomizing the association between characteristics
  my%data = %{$_[0]} ;
  my%new_dataset ;

  foreach my$col (0..12) {
    my@shuffled_data = shuffle(@{$data{$col}}) ;
    $new_dataset{$col} =  \@shuffled_data ;
  }

  return \%new_dataset ;
}

sub subsample_envt {
  ## Find minimum number of samples in an environment type (aquatic), and repeatedly subsample the other to evaluate the effect of unequal sample size
  my%data = %{$_[0]} ;
  my@symbioses = () ;
  my%new_dataset ;
  foreach my$col (0..12) {$new_dataset{$col} = () ;}
  my$colnum ;

  ## get environment counts
  my$aquatic = 0 ;
  my$terrestrial = 0 ;
  foreach my$i (@{$data{8}}) {
    if ($i == 0 || $i == 2) {$aquatic ++ ;}
    if ($i == 1) {$terrestrial ++ ;}
  }

  ## put all samples from min habitat count in new hash, and randomly subsample larger habitat count to same value
  my$min_count ;
  if ($aquatic < $terrestrial) {$min_count = $aquatic ;
    #print "min count aquatic: ", $min_count, "\n";
  }
  else {$min_count = $terrestrial ; print "min count terrestrial: ", $min_count, "\n";}

  my$Tcount = 0 ;
  my$Acount = 0 ;

  my@shuffled_indices = shuffle(0..$#{$data{8}}) ;

  foreach my$i (@shuffled_indices) {
    if ($data{8}[$i] == 0 || $data{8}[$i] == 2) {
      if ($Acount < $min_count) {
        foreach my$column (0..12) {push @{$new_dataset{$column}}, $data{$column}[$i] ;}
        $Acount ++ ;
      }
    }
    if ($data{8}[$i] == 1) {
      if ($Tcount < $min_count) {
        foreach my$column (0..12) {push @{$new_dataset{$column}}, $data{$column}[$i] ;}
        $Tcount ++ ;
      }
    }
  }

  return \%new_dataset ;
}

sub subsample_taxa {
  ## Find minimum number of samples per taxon and subsample all to that number
  my$taxon_level = $_[0] ;
  my%data = %{$_[1]} ;
  my@htaxa ;
  my@staxa ;
  my%new_dataset ;
  foreach my$col (0..12) {$new_dataset{$col} = () ;}
  my$colnum ;

  ## get taxa
  if ($taxon_level eq "high") {
    @htaxa = @{$data{2}} ;
    @staxa = @{$data{4}} ;
    $colnum = 2 ;
  }

  if ($taxon_level eq "low") {
    @htaxa = @{$data{3}} ;
    @staxa = @{$data{5}} ;
    $colnum = 3 ;
  }

  ## remove duplicates
  my%htaxa ;
  my%staxa ;
  foreach my$i (@htaxa) {$htaxa{$i} = "" ;}
  foreach my$i (@staxa) {$staxa{$i} = "" ;}
  ## save non-redundant taxa in new arrays
  my@nr_htaxa = keys %htaxa ;
  my@nr_staxa = keys %staxa ;

  ## make new dataset by randomly selecting one symbiosis from each taxon
  my@indices = () ;
  ## shuffle symbiont taxa indices
  ## Approach: This approach picks the first row that matches each taxon in the List
  ## by shuffling the order the list of rows is accessed in accomplishes random taxon selection
  my@shuffled_syms = shuffle(0..$#staxa) ;
  ## iterate through shuffled indices, recording one of each taxon to %new_dataset, and deleting the taxon from %staxa to avoid duplicates
  foreach my$index (@shuffled_syms) {
    ## if taxon is still in the %staxa hash, count this as the representative for that taxon
    if (exists $staxa{$data{$colnum}[$index]} && exists $htaxa{$data{($colnum-2)}[$index]}) {
      ## record values for every column for that row into the new, subsampled hash
      foreach my$i (0..12) {push @{$new_dataset{$i}}, $data{$i}[$index] ;}
      delete $staxa{$data{$colnum}[$index]} ;
      delete $htaxa{$data{($colnum-2)}[$index]} ;
    }
    else {
    }
  }

  return \%new_dataset ;
}

sub resample_data {
  my$dataset = $_[0] ;
  my%empirical_data = %{$_[1]} ;
  my$string = $_[2] ;
  my$function = $_[3] ;

  ## randomize characters 1000x
  my@Hp = () ;
  my@Mp = () ;
  my@Vp = () ;
  my@HAp = () ;
  my@MAp = () ;
  my@VAp = () ;
  my@HTp = () ;
  my@MTp = () ;
  my@VTp = () ;
  my@MIp = () ;
  my@VIp = () ;
  my@HNp = () ;
  my@HDp = () ;
  my@HMup = () ;
  my@HUkp = () ;
  my@HMap = () ;
  my@MNp = () ;
  my@MDp = () ;
  my@MMup = () ;
  my@MUkp = () ;
  my@MMap = () ;
  my@VNp = () ;
  my@VDp = () ;
  my@VMup = () ;
  my@VUkp = () ;
  my@VMap = () ;

  open SUBSAMPLEOUT, ">${string}_transmission_data.txt" or die "cannot open ${string}_transmission_data.txt\n" ;
  open OUT1, ">${string}_transmission_pvalues.txt" or die "cannot open ${string}_transmission_pvalues.txt\n" ;

  foreach my$count (0..999) {
    my$new_data ;

    if ($function =~ m/permute/) {$new_data = randomize_data($dataset) ;}
    elsif ($function =~ m/environment/) {$new_data = subsample_envt($dataset) ;}
    elsif ($function =~ m/hitaxa/) {$new_data = subsample_taxa("high", $dataset) ;}
    elsif ($function =~ m/lotaxa/) {$new_data = subsample_taxa("low", $dataset) ;}
    else {print "\"", $function, "\" term not working\n" ;}

    my($horizontalra, $mixedra, $verticalra) = get_rates($new_data) ;
    my($hterrestrialra, $mterrestrialra, $vterrestrialra, $haquaticra, $maquaticra, $vaquaticra) = get_env($new_data) ;
    my($externalverticalra, $internalverticalra, $externalmixedra, $internalmixedra) = get_routes($new_data) ;
    my($horiznutra, $mixnutra, $vertnutra, $horizdefra, $mixdefra, $vertdefra, $horizmultra, $mixmultra, $vertmultra, $horizunkra, $mixunkra, $vertunkra, $horizmanipra, $mixmanipra, $vertmanipra) = get_function($new_data) ;

    # in and out arrays are in the same order (horizontal,mixed,vertical)
    my@variables = ($horizontalra, $mixedra, $verticalra) ;
    my$percmode = print_to_open_file (*SUBSAMPLEOUT, "perc_hmv_1evlvra".$count, \@variables) ;
    my@percmode = @{$percmode} ;
    push @Hp, $percmode[0] ;
    push @Mp, $percmode[1] ;
    push @Vp, $percmode[2] ;

    @variables = ($haquaticra,$maquaticra,$vaquaticra) ;
    my$percaq = print_to_open_file (*SUBSAMPLEOUT, "aquaticra".$count, \@variables) ;
    my@percaq = @{$percaq} ;
    @variables = ($hterrestrialra,$mterrestrialra,$vterrestrialra) ;
    my$percter = print_to_open_file (*SUBSAMPLEOUT, "terrestrialra".$count, \@variables) ;
    my@percter = @{$percter} ;
    push @HAp, $percaq[0] ;
    push @MAp, $percaq[1] ;
    push @VAp, $percaq[2] ;
    push @HTp, $percter[0] ;
    push @MTp, $percter[1] ;
    push @VTp, $percter[2] ;

    ## Since this one is calculated differently - % internal transmission for mixed vs vertical - it is calculated with at different subroutine
    ## The array output only has % internal output (mixed, vertical)
    @variables = ($internalmixedra,$externalmixedra, $internalverticalra, $externalverticalra) ;
    my$percroute = print_to_open_file2 (*SUBSAMPLEOUT, "routera".$count, \@variables) ;
    my@percroute = @{$percroute} ;
    push @MIp, $percroute[0] ;
    push @VIp, $percroute[1] ;

    @variables = ($horiznutra,$mixnutra,$vertnutra) ;
    my$percnut = print_to_open_file (*SUBSAMPLEOUT, "nutritionra".$count, \@variables) ;
    my@percnut = @{$percnut} ;
    @variables = ($horizdefra,$mixdefra,$vertdefra) ;
    my$percdef = print_to_open_file (*SUBSAMPLEOUT, "defensera".$count, \@variables) ;
    my@percdef = @{$percdef} ;
    @variables = ($horizmultra,$mixmultra,$vertmultra) ;
    my$percmult = print_to_open_file (*SUBSAMPLEOUT, "multicomponentra".$count, \@variables) ;
    my@percmult = @{$percmult} ;
    @variables = ($horizunkra,$mixunkra,$vertunkra) ;
    my$percunk = print_to_open_file (*SUBSAMPLEOUT, "unknownra".$count, \@variables) ;
    my@percunk = @{$percunk} ;
    @variables = ($horizmanipra,$mixmanipra,$vertmanipra) ;
    my$percmanip = print_to_open_file (*SUBSAMPLEOUT, "manipulativera".$count, \@variables) ;
    my@percmanip = @{$percmanip} ;
    push @HNp, $percnut[0] ;
    push @HDp, $percdef[0] ;
    push @HMup, $percmult[0] ;
    push @HUkp, $percunk[0] ;
    push @HMap, $percmanip[0] ;
    push @MNp, $percnut[1] ;
    push @MDp, $percdef[1] ;
    push @MMup, $percmult[1] ;
    push @MUkp, $percunk[1] ;
    push @MMap, $percmanip[1] ;
    push @VNp, $percnut[2] ;
    push @VDp, $percdef[2] ;
    push @VMup, $percmult[2] ;
    push @VUkp, $percunk[2] ;
    push @VMap, $percmanip[2] ;

    $count ++ ;
  }
  close SUBSAMPLEOUT ;

  ########################
  my$below = 0 ;
  my$above = 0 ;
  my$equal = 0 ;

  foreach my$p (@Hp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"H"}) {$above ++ ;}
    elsif ($p < $empirical_data{"H"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"H"}, "\n" ;
  print OUT1 "horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"H"}, "\n" ;
  print OUT1 "horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"H"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@Mp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"M"}) {$above ++ ;}
    elsif ($p < $empirical_data{"M"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"M"}, "\n" ;
  print OUT1 "mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"M"}, "\n" ;
  print OUT1 "mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"M"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@Vp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"V"}) {$above ++ ;}
    elsif ($p < $empirical_data{"V"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"V"}, "\n" ;
  print OUT1 "vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"V"}, "\n" ;
  print OUT1 "vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"V"}, "\n\n" ;


  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HAp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HA"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HA"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "aquatic horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HA"}, "\n" ;
  print OUT1 "aquatic horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HA"}, "\n" ;
  print OUT1 "aquatic horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HA"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MAp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MA"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MA"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "aquatic mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MA"}, "\n" ;
  print OUT1 "aquatic mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MA"}, "\n" ;
  print OUT1 "aquatic mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MA"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VAp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VA"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VA"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "aquatic vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VA"}, "\n" ;
  print OUT1 "aquatic vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VA"}, "\n" ;
  print OUT1 "aquatic vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VA"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HTp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HT"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HT"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "terrestrial horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HT"}, "\n" ;
  print OUT1 "terrestrial horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HT"}, "\n" ;
  print OUT1 "terrestrial horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HT"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MTp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MT"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MT"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "terrestrial mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MT"}, "\n" ;
  print OUT1 "terrestrial mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MT"}, "\n" ;
  print OUT1 "terrestrial mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MT"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VTp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VT"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VT"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "terrestrial vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VT"}, "\n" ;
  print OUT1 "terrestrial vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VT"}, "\n" ;
  print OUT1 "terrestrial vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VT"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MIp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MI"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MI"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "mixed internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MI"}, "\n" ;
  print OUT1 "mixed internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MI"}, "\n" ;
  print OUT1 "mixed internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MI"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VIp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VI"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VI"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "vertical internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VI"}, "\n" ;
  print OUT1 "vertical internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VI"}, "\n" ;
  print OUT1 "vertical internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VI"}, "\n\n" ;

  ##############################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HNp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HN"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HN"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "nutrition functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HN"}, "\n" ;
  print OUT1 "nutrition functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HN"}, "\n" ;
  print OUT1 "nutrition functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HN"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HDp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "defense functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HD"}, "\n" ;
  print OUT1 "defense functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HD"}, "\n" ;
  print OUT1 "defense functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HD"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HMup) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "multiple-function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMu"}, "\n" ;
  print OUT1 "multiple-function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMu"}, "\n" ;
  print OUT1 "multiple-function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMu"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HUkp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HUk"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HUk"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "unknown function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HUk"}, "\n" ;
  print OUT1 "unknown function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HUk"}, "\n" ;
  print OUT1 "unknown function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HUk"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@HMap) {
		if (!$p) {next;}
    if ($p > $empirical_data{"HMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "manipulative functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMa"}, "\n" ;
  print OUT1 "manipulative functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMa"}, "\n" ;
  print OUT1 "manipulative functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMa"}, "\n\n" ;

  ##########################
  ########################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MNp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MN"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MN"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "nutrition functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MN"}, "\n" ;
  print OUT1 "nutrition functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MN"}, "\n" ;
  print OUT1 "nutrition functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MN"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MDp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "defense functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MD"}, "\n" ;
  print OUT1 "defense functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MD"}, "\n" ;
  print OUT1 "defense functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MD"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MMup) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "multiple-functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMu"}, "\n" ;
  print OUT1 "multiple-functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMu"}, "\n" ;
  print OUT1 "multiple-functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMu"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MUkp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MUk"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MUk"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "unknown functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MUk"}, "\n" ;
  print OUT1 "unknown functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MUk"}, "\n" ;
  print OUT1 "unknown functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MUk"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@MMap) {
		if (!$p) {next;}
    if ($p > $empirical_data{"MMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "manipulative functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMa"}, "\n" ;
  print OUT1 "manipulative functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMa"}, "\n" ;
  print OUT1 "manipulative functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMa"}, "\n\n" ;

  ###############################
  ########################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VNp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VN"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VN"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "nutrition functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VN"}, "\n" ;
  print OUT1 "nutrition functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VN"}, "\n" ;
  print OUT1 "nutrition functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VN"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VDp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "defense functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VD"}, "\n" ;
  print OUT1 "defense functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VD"}, "\n" ;
  print OUT1 "defense functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VD"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VMup) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "multiple-functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMu"}, "\n" ;
  print OUT1 "multiple-functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMu"}, "\n" ;
  print OUT1 "multiple-functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMu"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VUkp) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VUk"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VUk"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "unknown functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VUk"}, "\n" ;
  print OUT1 "unknown functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VUk"}, "\n" ;
  print OUT1 "unknown functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VUk"}, "\n\n" ;

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  foreach my$p (@VMap) {
		if (!$p) {next;}
    if ($p > $empirical_data{"VMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
  print OUT1 "manipulative functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMa"}, "\n" ;
  print OUT1 "manipulative functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMa"}, "\n" ;
  print OUT1 "manipulative functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMa"}, "\n\n" ;

  close OUT1 ;
}
