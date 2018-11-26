use strict ;
use warnings ;
use Sort::Naturally ;
use List::Util qw/shuffle/;

##usage: perl permute_collect_transmission_data.pl dataset.txt

my$unfiltered_symbioses = $ARGV[0] ;

open IN, "<$unfiltered_symbioses" or die "cannot open $unfiltered_symbioses\n" ;

my@column_names ;
my%unfiltered_symbioses ;
my$symbiosis_count = 0 ;
my%total_count ;

my%total_counts ;
my%total_aqcounts ;
my%total_tercounts ;

while (<IN>) {
  chomp ;
  ## get column labels
  ## #0=host_class	1=host_family	2=symbiont_phylum	3=symbiont_family	4=mode	5=route	6=environment	7=function  8=evidence
  if ($_ =~ /^#/) {
    $_ =~ s/^#// ;
    @column_names = split(/\t/, $_) ;

    ## make hash of columns, each with an array for its rows
    ## therefore, the relationships between data points in a row will be maintained by their order across arrays in the hash
    foreach my$i (0..$#column_names) {@{$unfiltered_symbioses{$i}} = () ;}
  }

  else {
    my@split = split(/\t/, $_) ;
    foreach my$i (0..$#split) {push @{$unfiltered_symbioses{$i}}, $split[$i] ;}
    $symbiosis_count ++ ;
  }
}

close IN ;

foreach my$i (nsort keys %unfiltered_symbioses) {
#  print $i, "\n", join(",", @{$unfiltered_symbioses{$i}}), "\n\n" ;
}

### Identify phylogenetically independent symbiotic host+sym combinations for independent symbiosis events
## cycle through symbioses and tally how many symbiont+host pairings there are
my%unique_symbioses ;
foreach my$i (0..$#{$unfiltered_symbioses{0}}) {
  #use families for deliniation if there isn't data on symbiotic events
  if ($unfiltered_symbioses{0}[$i] =~ m/\?/ && $unfiltered_symbioses{1}[$i] =~ m/\?/) {
    $unique_symbioses{"$unfiltered_symbioses{5}[$i]\t$unfiltered_symbioses{3}[$i]"} += 1 ;
  }
  elsif ($unfiltered_symbioses{0}[$i] =~ m/\?/) {
    $unique_symbioses{"$unfiltered_symbioses{5}[$i]\t$unfiltered_symbioses{1}[$i]"} += 1 ;
  }
  elsif ($unfiltered_symbioses{1}[$i] =~ m/\?/) {
    $unique_symbioses{"$unfiltered_symbioses{0}[$i]\t$unfiltered_symbioses{3}[$i]"} += 1 ;
  }
  else {
    $unique_symbioses{"$unfiltered_symbioses{0}[$i]\t$unfiltered_symbioses{1}[$i]"} += 1 ;
  }
}

my$highest = 0 ;
my$highest_symbiosis = "" ;
foreach my$i (keys %unique_symbioses) {
  if ($unique_symbioses{$i} > $highest) {
    $highest = $unique_symbioses{$i} ;
    $highest_symbiosis = $i ;
  }
}

print "unique symbiotic groups: ", scalar(keys %unique_symbioses), "\n" ;
print "largest sampled symbiotic group: ", $highest_symbiosis, "\t", $highest, "\n" ;

#resample independent symbiotic groups # times equal to number of taxa in most oversampled group
foreach my$i (0 ..($highest-1)) {
  print "resample: ", $i, "\n" ;
  my%unique_subsample ;

  ## iterate through list of symbioses to resample the unique symbiosis events with different constituting taxa
  foreach my$u (nsort keys %unique_symbioses) {
    ## split sym and host names obtained from unique list of associations
    my@split = split(/\t/, $u) ;
    my$sym = $split[0] ;
    my$host = $split[1] ;

    ## set flag to indicate that the association has been fulfilled
    my$found = "no" ;

    while ($found eq "no") {
      ## pick a random integer between 0 and the total number of associations in $unfiltered_symbioses
      my$x = int rand($symbiosis_count) ;

      if ($unfiltered_symbioses{0}[$x] eq $sym || $unfiltered_symbioses{5}[$x] eq $sym) {
        if ($unfiltered_symbioses{1}[$x] eq $host || $unfiltered_symbioses{3}[$x] eq $host) {
          foreach my$k (0..12) {
            push @{$unique_subsample{$k}}, $unfiltered_symbioses{$k}[$x] ;
          }

          $found = "yes" ;
        }
      }
    }
  }

  ##check length of subsample_envt
  if (scalar(@{$unique_subsample{0}}) != scalar(keys %unique_symbioses)) {
    print "check why ", scalar(@{$unique_subsample{0}}), " does not equal ", scalar(keys %unique_symbioses), "\n" ;
  }

  process_data(\%unique_subsample, "independent$i") ;

}



sub print_to_open_file {
  my$handle = $_[0] ;
  my$arrayname = $_[1] ;
  my$variables = $_[2] ;
  my@variables = @{$variables} ;
  my$sum = 0 ;
  my@new_variables ;
  #add up counts for horizontal, mixed, and vertical transmission under given category (e.g., environment type, functional type, etc.)
  foreach my$v (@variables) {$sum += $v ;}
  if ($sum > 0) {
    print $handle $arrayname, " = c(" ;
    #calculate and print percent of horizontal, and vertical transmission
    foreach my$i (0..($#variables-1)) {
      my$percent = $variables[$i]/$sum*100 ;
      print $handle $percent, "," ;
      push @new_variables, $percent ;
    }
    #calculate and print percent of vertical transmission
    my$percent = $variables[$#variables]/$sum*100 ;
    print $handle $percent, ")\n" ;
    push @new_variables, $percent ;
  }
  #some permutations have sums of zero because none of the given category were sampled
  else {
    #  print $arrayname, " won't print\n" ;
    #add placeholders so the arrays aren't empty
    push @new_variables, "NA", "NA","NA" ;
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
  my%aquatic ;
  my%terrestrial ;

  ## ENV TYPE: get % of mode in each ENV
  ## unfiltered
  my$HT = 0 ;
  my$HA = 0 ;
  my$MT = 0 ;
  my$MA = 0 ;
  my$VT = 0 ;
  my$VA = 0 ;
  foreach my$i (0..$#{$data{6}}) {
    #if environment type is marine (0) or freshwater (2), add to aquatic lists
    if ($data{8}[$i] == 0 || $data{8}[$i] == 2) {
      if ($data{6}[$i] == 0) {$HA ++ ;}
      if ($data{6}[$i] == 1) {$MA ++ ;}
      if ($data{6}[$i] == 2) {$VA ++;}
      #add aquatic symbiosis to separate dataset
      foreach my$j (nsort keys %data) {push @{$aquatic{$j}}, $data{$j}[$i] ;}
    }
    #if environment type is terrestrial (1), add to terrestrial lists
    if ($data{8}[$i] == 1) {
      if ($data{6}[$i] == 0) {$HT ++ ;}
      if ($data{6}[$i] == 1) {$MT ++ ;}
      if ($data{6}[$i] == 2) {$VT ++;}
      #add aquatic symbiosis to separate dataset
      foreach my$j (nsort keys %data) {push @{$terrestrial{$j}}, $data{$j}[$i] ;}
    }
  }

    return ($HT, $MT, $VT, $HA, $MA, $VA, \%aquatic, \%terrestrial) ;
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
  if (exists $data{2} && exists $data{4}) {
    if ($taxon_level eq "high") {
      @htaxa = @{$data{2}} ;
      @staxa = @{$data{4}} ;
      $colnum = 2 ;
    }
  }

  if (exists $data{3} && exists $data{5}) {
    if ($taxon_level eq "low" && @{$data{3}} && @{$data{5}}) {
      @htaxa = @{$data{3}} ;
      @staxa = @{$data{5}} ;
      $colnum = 3 ;
    }
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

  foreach my$count (0..99) {
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
    push @HNp, $percnut[0] ;
    push @MNp, $percnut[1] ;
    push @VNp, $percnut[2] ;
    @variables = ($horizdefra,$mixdefra,$vertdefra) ;
    my$percdef = print_to_open_file (*SUBSAMPLEOUT, "defensera".$count, \@variables) ;
    my@percdef = @{$percdef} ;
    push @HDp, $percdef[0] ;
    push @MDp, $percdef[1] ;
    push @VDp, $percdef[2] ;
    @variables = ($horizmultra,$mixmultra,$vertmultra) ;
    my$percmult = print_to_open_file (*SUBSAMPLEOUT, "multicomponentra".$count, \@variables) ;
    my@percmult = @{$percmult} ;
    push @HMup, $percmult[0] ;
    push @MMup, $percmult[1] ;
    push @VMup, $percmult[2] ;
    @variables = ($horizunkra,$mixunkra,$vertunkra) ;
    my$percunk = print_to_open_file (*SUBSAMPLEOUT, "unknownra".$count, \@variables) ;
    my@percunk = @{$percunk} ;
    push @HUkp, $percunk[0] ;
    push @MUkp, $percunk[1] ;
    push @VUkp, $percunk[2] ;
    @variables = ($horizmanipra,$mixmanipra,$vertmanipra) ;
    my$percmanip = print_to_open_file (*SUBSAMPLEOUT, "manipulativera".$count, \@variables) ;
    my@percmanip = @{$percmanip} ;
    push @HMap, $percmanip[0] ;
    push @MMap, $percmanip[1] ;
    push @VMap, $percmanip[2] ;

    $count ++ ;
  }
  close SUBSAMPLEOUT ;

  ########################
  my$below = 0 ;
  my$above = 0 ;
  my$equal = 0 ;
  my$sum = 0 ;
  my$total ;
  foreach my$i (@Hp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@Hp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"H"}) {$above ++ ;}
      elsif ($p < $empirical_data{"H"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"H"}, "\n" ;
    print OUT1 "horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"H"}, "\n" ;
    print OUT1 "horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"H"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@Mp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@Mp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"M"}) {$above ++ ;}
      elsif ($p < $empirical_data{"M"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"M"}, "\n" ;
    print OUT1 "mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"M"}, "\n" ;
    print OUT1 "mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"M"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@Vp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@Vp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"V"}) {$above ++ ;}
      elsif ($p < $empirical_data{"V"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"V"}, "\n" ;
    print OUT1 "vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"V"}, "\n" ;
    print OUT1 "vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"V"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HAp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HAp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HA"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HA"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "aquatic horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HA"}, "\n" ;
    print OUT1 "aquatic horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HA"}, "\n" ;
    print OUT1 "aquatic horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HA"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MAp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MAp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MA"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MA"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "aquatic mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MA"}, "\n" ;
    print OUT1 "aquatic mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MA"}, "\n" ;
    print OUT1 "aquatic mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MA"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VAp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VAp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VA"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VA"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "aquatic vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VA"}, "\n" ;
    print OUT1 "aquatic vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VA"}, "\n" ;
    print OUT1 "aquatic vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VA"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HTp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HTp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HT"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HT"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "terrestrial horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HT"}, "\n" ;
    print OUT1 "terrestrial horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HT"}, "\n" ;
    print OUT1 "terrestrial horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HT"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MTp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MTp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MT"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MT"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "terrestrial mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MT"}, "\n" ;
    print OUT1 "terrestrial mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MT"}, "\n" ;
    print OUT1 "terrestrial mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MT"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VTp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VTp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VT"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VT"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "terrestrial vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VT"}, "\n" ;
    print OUT1 "terrestrial vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VT"}, "\n" ;
    print OUT1 "terrestrial vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VT"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MIp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MIp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MI"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MI"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "mixed internal route:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MI"}, "\n" ;
    print OUT1 "mixed internal route:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MI"}, "\n" ;
    print OUT1 "mixed internal route:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MI"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VIp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VIp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VI"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VI"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "vertical internal route:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VI"}, "\n" ;
    print OUT1 "vertical internal route:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VI"}, "\n" ;
    print OUT1 "vertical internal route:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VI"}, "\n\n" ;
  }

  ##############################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HNp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HNp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HN"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HN"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "nutrition functions horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HN"}, "\n" ;
    print OUT1 "nutrition functions horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HN"}, "\n" ;
    print OUT1 "nutrition functions horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HN"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HDp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HDp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HD"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HD"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "defense functions horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HD"}, "\n" ;
    print OUT1 "defense functions horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HD"}, "\n" ;
    print OUT1 "defense functions horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HD"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HMup) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HMup) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HMu"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HMu"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "multiple-function horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HMu"}, "\n" ;
    print OUT1 "multiple-function horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HMu"}, "\n" ;
    print OUT1 "multiple-function horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HMu"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HUkp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HUkp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HUk"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HUk"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "unknown function horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HUk"}, "\n" ;
    print OUT1 "unknown function horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HUk"}, "\n" ;
    print OUT1 "unknown function horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HUk"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@HMap) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@HMap) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"HMa"}) {$above ++ ;}
      elsif ($p < $empirical_data{"HMa"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "manipulative functions horizontal mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"HMa"}, "\n" ;
    print OUT1 "manipulative functions horizontal mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"HMa"}, "\n" ;
    print OUT1 "manipulative functions horizontal mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"HMa"}, "\n\n" ;
  }

  ##########################
  ########################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MNp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MNp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MN"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MN"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "nutrition functions mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MN"}, "\n" ;
    print OUT1 "nutrition functions mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MN"}, "\n" ;
    print OUT1 "nutrition functions mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MN"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MDp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MDp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MD"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MD"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "defense functions mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MD"}, "\n" ;
    print OUT1 "defense functions mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MD"}, "\n" ;
    print OUT1 "defense functions mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MD"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MMup) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MMup) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MMu"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MMu"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "multiple-functions mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MMu"}, "\n" ;
    print OUT1 "multiple-functions mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MMu"}, "\n" ;
    print OUT1 "multiple-functions mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MMu"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MUkp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MUkp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MUk"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MUk"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "unknown functions mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MUk"}, "\n" ;
    print OUT1 "unknown functions mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MUk"}, "\n" ;
    print OUT1 "unknown functions mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MUk"}, "\n\n" ;
  }
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@MMap) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@MMap) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"MMa"}) {$above ++ ;}
      elsif ($p < $empirical_data{"MMa"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "manipulative functions mixed mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"MMa"}, "\n" ;
    print OUT1 "manipulative functions mixed mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"MMa"}, "\n" ;
    print OUT1 "manipulative functions mixed mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"MMa"}, "\n\n" ;
  }

  ###############################
  ########################
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VNp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VNp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VN"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VN"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "nutrition functions vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VN"}, "\n" ;
    print OUT1 "nutrition functions vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VN"}, "\n" ;
    print OUT1 "nutrition functions vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VN"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VDp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VDp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VD"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VD"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "defense functions vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VD"}, "\n" ;
    print OUT1 "defense functions vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VD"}, "\n" ;
    print OUT1 "defense functions vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VD"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VMup) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VMup) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VMu"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VMu"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "multiple-functions vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VMu"}, "\n" ;
    print OUT1 "multiple-functions vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VMu"}, "\n" ;
    print OUT1 "multiple-functions vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VMu"}, "\n\n" ;
  }

  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VUkp) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VUkp) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VUk"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VUk"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "unknown functions vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VUk"}, "\n" ;
    print OUT1 "unknown functions vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VUk"}, "\n" ;
    print OUT1 "unknown functions vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VUk"}, "\n\n" ;
  }
  ########################
  $below = 0 ;
  $above = 0 ;
  $equal = 0 ;
  $sum = 0 ;
  foreach my$i (@VMap) {
    if (!defined($i)) {next;}
		if($i eq "NA") {next ;}
    $sum += $i ;
  }
  if ($sum > 0) {
		foreach my$p (@VMap) {
			if ($p eq "NA") {next ;}
      if ($p > $empirical_data{"VMa"}) {$above ++ ;}
      elsif ($p < $empirical_data{"VMa"}) {$below ++ ;}
      else {$equal ++ ;}
    }
  }
  $total = $above+$below+$equal ;
  if ($total > 0) {
    print OUT1 "manipulative functions vertical mode:\t", $above, "(", $above/$total*100, "%)", " above the empirical value of ", $empirical_data{"VMa"}, "\n" ;
    print OUT1 "manipulative functions vertical mode:\t", $below, "(", $below/$total*100, "%)", " below the empirical value of ", $empirical_data{"VMa"}, "\n" ;
    print OUT1 "manipulative functions vertical mode:\t", $equal, "(", $equal/$total*100, "%)", " equal to the empirical value of ", $empirical_data{"VMa"}, "\n\n" ;
  }

  close OUT1 ;
}

sub process_data {
  my%symbioses = %{$_[0]} ;
  my$sample = $_[1] ;

  open OUTPUT, ">independent${sample}_transmission_analysis_output.txt" or die "cannot open independent${sample}_transmission_analysis_output.txt\n" ;

  my%empdata ;
  my%tally ;

  #### Print stats to file for plotting in R
  print OUTPUT "##RATES\n" ;
  ## all data
  my($horizontal, $mixed, $vertical) = get_rates(\%symbioses) ;
  my@variables = ($horizontal, $mixed, $vertical) ;
  my$variablesperc = print_to_open_file (*OUTPUT, "perc_hmv_1evlv", \@variables) ;
  my@variablesperc = @{$variablesperc} ;
  #@variables = ($horizontal2, $mixed2, $vertical2) ;
  #print_to_open_file (*OUTPUT, "perc_hmv_2evlv", \@variables) ;
  #@variables = ($horizontal3, $mixed3, $vertical3) ;
  #print_to_open_file (*OUTPUT, "perc_hmv_3evlv", \@variables) ;

  $empdata{"H"} = $variablesperc[0] ;
  $empdata{"M"} = $variablesperc[1] ;
  $empdata{"V"} = $variablesperc[2] ;


  ##### RATES - print to file for plotting proportions
  print "total: ", ($vertical+$horizontal+$mixed),"\n" ;

  print OUTPUT "\n##ENVIRONMENTS\n" ;
  my($hterrestrial, $mterrestrial, $vterrestrial, $haquatic, $maquatic, $vaquatic, $aquaticdata, $terrestrialdata) = get_env(\%symbioses) ;
  @variables = ($haquatic,$maquatic,$vaquatic) ;
  $variablesperc  = print_to_open_file (*OUTPUT, "aquatic", \@variables) ;
  @variablesperc = @{$variablesperc} ;

  $empdata{"HA"} = $variablesperc[0] ;
  $empdata{"MA"} = $variablesperc[1] ;
  $empdata{"VA"} = $variablesperc[2] ;

  @variables = ($hterrestrial,$mterrestrial,$vterrestrial) ;
  $variablesperc = print_to_open_file (*OUTPUT, "terrestrial", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HT"} = $variablesperc[0] ;
  $empdata{"MT"} = $variablesperc[1] ;
  $empdata{"VT"} = $variablesperc[2] ;
  #
  print OUTPUT "\n##ROUTES (internal)\n" ;
  ##print $internalmixed, "\t", $internalvertical, "\n" ;
  my@route_variables = (6,7) ;
  my($externalvertical, $internalvertical, $externalmixed, $internalmixed) = get_routes(\%symbioses) ;
  @variables = ($internalmixed,$externalmixed,$internalvertical,$externalvertical) ;
  $variablesperc = print_to_open_file2 (*OUTPUT, "route", \@variables) ;
  @variablesperc = @{$variablesperc} ;

  $empdata{"MI"} = $variablesperc[0] ;
  $empdata{"VI"} = $variablesperc[1] ;


  ##Get routes for aquatic and terrestrial %symbioses
  my%aquaticdata = %{$aquaticdata} ;
  my%terrestrialdata = %{$terrestrialdata} ;

  ($externalvertical, $internalvertical, $externalmixed, $internalmixed) = get_routes(\%aquaticdata) ;
  @variables = ($internalmixed,$externalmixed,$internalvertical,$externalvertical) ;
  $variablesperc = print_to_open_file2 (*OUTPUT, "aquatic_route", \@variables) ;

  ($externalvertical, $internalvertical, $externalmixed, $internalmixed) = get_routes(\%terrestrialdata) ;
  @variables = ($internalmixed,$externalmixed,$internalvertical,$externalvertical) ;
  $variablesperc = print_to_open_file2 (*OUTPUT, "terrestrial_route", \@variables) ;

  #
  print OUTPUT "\n##FUNCTIONS\n" ;
  my($horiznut, $mixnut, $vertnut, $horizdef, $mixdef, $vertdef, $horizmult, $mixmult, $vertmult, $horizunk, $mixunk, $vertunk, $horizmanip, $mixmanip, $vertmanip) = get_function(\%symbioses) ;
  @variables = ($horiznut,$mixnut,$vertnut) ;
  $variablesperc = print_to_open_file (*OUTPUT, "nutrition", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HN"} = $variablesperc[0] ;
  $empdata{"MN"} = $variablesperc[1] ;
  $empdata{"VN"} = $variablesperc[2] ;
  @variables = ($horizdef,$mixdef,$vertdef) ;
  $variablesperc = print_to_open_file (*OUTPUT, "defense", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HD"} = $variablesperc[0] ;
  $empdata{"MD"} = $variablesperc[1] ;
  $empdata{"VD"} = $variablesperc[2] ;
  @variables = ($horizmult,$mixmult,$vertmult) ;
  $variablesperc = print_to_open_file (*OUTPUT, "multicomponent", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HMu"} = $variablesperc[0] ;
  $empdata{"MMu"} = $variablesperc[1] ;
  $empdata{"VMu"} = $variablesperc[2] ;
  @variables = ($horizunk,$mixunk,$vertunk) ;
  $variablesperc = print_to_open_file (*OUTPUT, "unknown", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HUk"} = $variablesperc[0] ;
  $empdata{"MUk"} = $variablesperc[1] ;
  $empdata{"VUk"} = $variablesperc[2] ;
  @variables = ($horizmanip,$mixmanip,$vertmanip) ;
  $variablesperc = print_to_open_file (*OUTPUT, "manipulative", \@variables) ;
  @variablesperc = @{$variablesperc} ;
  $empdata{"HMa"} = $variablesperc[0] ;
  $empdata{"MMa"} = $variablesperc[1] ;
  $empdata{"VMa"} = $variablesperc[2]  ;


  ### Permutations and subsampling
  resample_data(\%symbioses, \%empdata, "${sample}_permuted", "permute") ;
  resample_data(\%symbioses, \%empdata, "${sample}_subsampleEnv", "environment") ;
  resample_data(\%symbioses, \%empdata, "${sample}_subsamplePhyla", "hitaxa") ;
  resample_data(\%symbioses, \%empdata, "${sample}_subsampleFamilies", "lotaxa") ;



  close OUTPUT ;

}
