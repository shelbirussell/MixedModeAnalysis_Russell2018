use strict ;
use warnings ;
use Sort::Naturally ;
use List::Util qw/shuffle/;

##usage: perl subsample_collect_transmission_data.pl dataset.txt

my$symbioses = $ARGV[0] ;

open IN, "<$symbioses" or die "cannot open $symbioses\n" ;

my@column_names ;
my%symbioses ;

while (<IN>) {
  chomp ;
  ## get column labels
  ## #0=host_class	1=host_family	2=symbiont_phylum	3=symbiont_family	4=mode	5=route	6=environment	7=function  8=evidence
  if ($_ =~ /^#/) {
    $_ =~ s/^#// ;

    @column_names = split(/\t/, $_) ;

    ## make hash of columns, each with an array for its rows
    ## therefore, the relationships between data points in a row will be maintained by their order across arrays in the hash
    foreach my$i (0..$#column_names) {@{$symbioses{$i}} = () ;}
  }

  else {
    my@split = split(/\t/, $_) ;
    foreach my$i (0..$#split) {push @{$symbioses{$i}}, $split[$i] ;}
  }
}

close IN ;

## subsample taxa 1000x
my@Hphi = () ;
my@Mphi = () ;
my@Vphi = () ;
my@HAphi = () ;
my@MAphi = () ;
my@VAphi = () ;
my@HTphi = () ;
my@MTphi = () ;
my@VTphi = () ;
my@MIphi = () ;
my@VIphi = () ;
my@HNphi = () ;
my@HDphi = () ;
my@HMuphi = () ;
my@HUkphi = () ;
my@HMaphi = () ;
my@MNphi = () ;
my@MDphi = () ;
my@MMuphi = () ;
my@MUkphi = () ;
my@MMaphi = () ;
my@VNphi = () ;
my@VDphi = () ;
my@VMuphi = () ;
my@VUkphi = () ;
my@VMaphi = () ;

my@Hplo = () ;
my@Mplo = () ;
my@Vplo = () ;
my@HAplo = () ;
my@MAplo = () ;
my@VAplo = () ;
my@HTplo = () ;
my@MTplo = () ;
my@VTplo = () ;
my@MIplo = () ;
my@VIplo = () ;
my@HNplo = () ;
my@HDplo = () ;
my@HMuplo = () ;
my@HUkplo = () ;
my@HMaplo = () ;
my@MNplo = () ;
my@MDplo = () ;
my@MMuplo = () ;
my@MUkplo = () ;
my@MMaplo = () ;
my@VNplo = () ;
my@VDplo = () ;
my@VMuplo = () ;
my@VUkplo = () ;
my@VMaplo = () ;

open SUBSAMPLEOUT, "> subsampledtaxa_transmission_output.txt" or die "cannot open subsampledtaxa_transmission_output.txt\n" ;

foreach my$count (0..999) {
  my$subsampled_highorder_taxa = subsample_taxa("high", \%symbioses) ;
  my%subsampled_highorder_taxa = %{$subsampled_highorder_taxa} ;
  my($horizontalhi, $mixedhi, $verticalhi) = get_rates(\%subsampled_highorder_taxa) ;
  #print $horizontalhi, "\t", $mixedhi, "\t", $verticalhi, "\n" ;
  my($hterrestrialhi, $mterrestrialhi, $vterrestrialhi, $haquatichi, $maquatichi, $vaquatichi) = get_env(\%subsampled_highorder_taxa) ;
  my($externalverticalhi, $internalverticalhi, $externalmixedhi, $internalmixedhi) = get_routes(\%subsampled_highorder_taxa) ;
  my($horiznuthi, $mixnuthi, $vertnuthi, $horizdefhi, $mixdefhi, $vertdefhi, $horizmulthi, $mixmulthi, $vertmulthi, $horizunkhi, $mixunkhi, $vertunkhi, $horizmaniphi, $mixmaniphi, $vertmaniphi) = get_function(\%subsampled_highorder_taxa) ;

  my$subsampled_loworder_taxa = subsample_taxa("low", \%symbioses) ;
  my%subsampled_loworder_taxa = %{$subsampled_loworder_taxa} ;
  my($horizontallo, $mixedlo, $verticallo) = get_rates(\%subsampled_loworder_taxa) ;
  #print $horizontallo, "\t", $mixedlo, "\t", $verticallo, "\n" ;
  my($hterrestriallo, $mterrestriallo, $vterrestriallo, $haquaticlo, $maquaticlo, $vaquaticlo) = get_env(\%subsampled_loworder_taxa) ;
  my($externalverticallo, $internalverticallo, $externalmixedlo, $internalmixedlo) = get_routes(\%subsampled_loworder_taxa) ;
  my($horiznutlo, $mixnutlo, $vertnutlo, $horizdeflo, $mixdeflo, $vertdeflo, $horizmultlo, $mixmultlo, $vertmultlo, $horizunklo, $mixunklo, $vertunklo, $horizmaniplo, $mixmaniplo, $vertmaniplo) = get_function(\%subsampled_loworder_taxa) ;

  my@variables = ($horizontalhi, $mixedhi, $verticalhi) ;
  my$percmodehi = print_to_open_file (*SUBSAMPLEOUT, "perc_hmv_1evlvhi".$count, \@variables) ;
  my@percmodehi = @{$percmodehi} ;
  @variables = ($horizontallo, $mixedlo, $verticallo) ;
  my$percmodelo = print_to_open_file (*SUBSAMPLEOUT, "perc_hmv_1evlvlo".$count, \@variables) ;
  my@percmodelo = @{$percmodelo} ;
  push @Hphi, $percmodehi[0] ;
  push @Mphi, $percmodehi[1] ;
  push @Vphi, $percmodehi[2] ;
  push @Hplo, $percmodelo[0] ;
  push @Mplo, $percmodelo[1] ;
  push @Vplo, $percmodelo[2] ;

  @variables = ($haquatichi,$maquatichi,$vaquatichi) ;
  my$percaqhi = print_to_open_file (*SUBSAMPLEOUT, "aquatichi".$count, \@variables) ;
  my@percaqhi = @{$percaqhi} ;
  @variables = ($haquaticlo,$maquaticlo,$vaquaticlo) ;
  my$percaqlo = print_to_open_file (*SUBSAMPLEOUT, "aquaticlo".$count, \@variables) ;
  my@percaqlo = @{$percaqlo} ;

  @variables = ($hterrestrialhi,$mterrestrialhi,$vterrestrialhi) ;
  my$percterhi = print_to_open_file (*SUBSAMPLEOUT, "terrestrialhi".$count, \@variables) ;
  my@percterhi = @{$percterhi} ;
  @variables = ($hterrestriallo,$mterrestriallo,$vterrestriallo) ;
  my$percterlo = print_to_open_file (*SUBSAMPLEOUT, "terrestriallo".$count, \@variables) ;
  my@percterlo = @{$percterlo} ;
  push @HAplo, $percaqlo[0] ;
  push @MAplo, $percaqlo[1] ;
  push @VAplo, $percaqlo[2] ;
  push @HTplo, $percterlo[0] ;
  push @MTplo, $percterlo[1] ;
  push @VTplo, $percterlo[2] ;
  push @HAphi, $percaqhi[0] ;
  push @MAphi, $percaqhi[1] ;
  push @VAphi, $percaqhi[2] ;
  push @HTphi, $percterhi[0] ;
  push @MTphi, $percterhi[1] ;
  push @VTphi, $percterhi[2] ;

  @variables = ($internalmixedhi,$externalmixedhi, $internalverticalhi, $externalverticalhi) ;
  my$percroutehi = print_to_open_file2 (*SUBSAMPLEOUT, "routehi".$count, \@variables) ;
  my@percroutehi = @{$percroutehi} ;
  push @MIphi, $percroutehi[0] ;
  push @VIphi, $percroutehi[1] ;
  @variables = ($internalmixedlo,$externalmixedlo, $internalverticallo, $externalverticallo) ;
  my$percroutelo = print_to_open_file2 (*SUBSAMPLEOUT, "routelo".$count, \@variables) ;
  my@percroutelo = @{$percroutelo} ;
  push @MIplo, $percroutelo[0] ;
  push @VIplo, $percroutelo[1] ;

  @variables = ($horiznuthi,$mixnuthi,$vertnuthi) ;
  my$percnuthi = print_to_open_file (*SUBSAMPLEOUT, "nutritionhi".$count, \@variables) ;
  my@percnuthi = @{$percnuthi} ;
  @variables = ($horizdefhi,$mixdefhi,$vertdefhi) ;
  my$percdefhi = print_to_open_file (*SUBSAMPLEOUT, "defensehi".$count, \@variables) ;
  my@percdefhi = @{$percdefhi} ;
  @variables = ($horizmulthi,$mixmulthi,$vertmulthi) ;
  my$percmulthi = print_to_open_file (*SUBSAMPLEOUT, "multicomponenthi".$count, \@variables) ;
  my@percmulthi = @{$percmulthi} ;
  @variables = ($horizunkhi,$mixunkhi,$vertunkhi) ;
  my$percunkhi = print_to_open_file (*SUBSAMPLEOUT, "unknownhi".$count, \@variables) ;
  my@percunkhi = @{$percunkhi} ;
  @variables = ($horizmaniphi,$mixmaniphi,$vertmaniphi) ;
  my$percmaniphi = print_to_open_file (*SUBSAMPLEOUT, "manipulativehi".$count, \@variables) ;
  my@percmaniphi = @{$percmaniphi} ;

  @variables = ($horiznutlo,$mixnutlo,$vertnutlo) ;
  my$percnutlo = print_to_open_file (*SUBSAMPLEOUT, "nutritionlo".$count, \@variables) ;
  my@percnutlo = @{$percnutlo} ;
  @variables = ($horizdeflo,$mixdeflo,$vertdeflo) ;
  my$percdeflo = print_to_open_file (*SUBSAMPLEOUT, "defenselo".$count, \@variables) ;
  my@percdeflo = @{$percdeflo} ;
  @variables = ($horizmultlo,$mixmultlo,$vertmultlo) ;
  my$percmultlo = print_to_open_file (*SUBSAMPLEOUT, "multicomponentlo".$count, \@variables) ;
  my@percmultlo = @{$percmultlo} ;
  @variables = ($horizunklo,$mixunklo,$vertunklo) ;
  my$percunklo = print_to_open_file (*SUBSAMPLEOUT, "unknownlo".$count, \@variables) ;
  my@percunklo = @{$percunklo} ;
  @variables = ($horizmaniplo,$mixmaniplo,$vertmaniplo) ;
  my$percmaniplo = print_to_open_file (*SUBSAMPLEOUT, "manipulativelo".$count, \@variables) ;
  my@percmaniplo = @{$percmaniplo} ;

  push @HNphi, $percnuthi[0] ;
  push @HDphi, $percdefhi[0] ;
  push @HMuphi, $percmulthi[0] ;
  push @HUkphi, $percunkhi[0] ;
  push @HMaphi, $percmaniphi[0] ;
  push @MNphi, $percnuthi[1] ;
  push @MDphi, $percdefhi[1] ;
  push @MMuphi, $percmulthi[1] ;
  push @MUkphi, $percunkhi[1] ;
  push @MMaphi, $percmaniphi[1] ;
  push @VNphi, $percnuthi[2] ;
  push @VDphi, $percdefhi[2] ;
  push @VMuphi, $percmulthi[2] ;
  push @VUkphi, $percunkhi[2] ;
  push @VMaphi, $percmaniphi[2] ;

  push @HNplo, $percnutlo[0] ;
  push @HDplo, $percdeflo[0] ;
  push @HMuplo, $percmultlo[0] ;
  push @HUkplo, $percunklo[0] ;
  push @HMaplo, $percmaniplo[0] ;
  push @MNplo, $percnutlo[1] ;
  push @MDplo, $percdeflo[1] ;
  push @MMuplo, $percmultlo[1] ;
  push @MUkplo, $percunklo[1] ;
  push @MMaplo, $percmaniplo[1] ;
  push @VNplo, $percnutlo[2] ;
  push @VDplo, $percdeflo[2] ;
  push @VMuplo, $percmultlo[2] ;
  push @VUkplo, $percunklo[2] ;
  push @VMaplo, $percmaniplo[2] ;
  $count ++ ;
}
close SUBSAMPLEOUT ;

######
## Test if empirical value is within permuted values
##input measured values to test against permuted dataset
my%empirical_data ;
$empirical_data{"H"} = 15.9453302961276 ;
$empirical_data{"M"} = 33.7129840546697 ;
$empirical_data{"V"} = 50.3416856492027 ;
$empirical_data{"HA"} = 51.5151515151515 ;
$empirical_data{"MA"} = 39.3939393939394 ;
$empirical_data{"VA"} = 9.09090909090909 ;
$empirical_data{"HT"} = 9.6514745308311 ;
$empirical_data{"MT"} = 32.7077747989276 ;
$empirical_data{"VT"} = 57.6407506702413 ;
$empirical_data{"MI"} = 64.8648648648649 ;
$empirical_data{"VI"} = 81.9004524886878 ;
$empirical_data{"HN"} = 25.1046025104602 ;
$empirical_data{"HD"} = 2.04081632653061 ;
$empirical_data{"HMu"} = 0 ;
$empirical_data{"HUk"} = 8.10810810810811 ;
$empirical_data{"HMa"} = 0 ;
$empirical_data{"MN"} = 21.3389121338912 ;
$empirical_data{"MD"} = 67.3469387755102 ;
$empirical_data{"MMu"} = 52.9411764705882 ;
$empirical_data{"MUk"} = 28.8288288288288 ;
$empirical_data{"MMa"} = 100 ;
$empirical_data{"VN"} = 53.5564853556485 ;
$empirical_data{"VD"} = 30.6122448979592 ;
$empirical_data{"VMu"} = 47.0588235294118 ;
$empirical_data{"VUk"} = 63.0630630630631 ;
$empirical_data{"VMa"} = 0 ;


########################
print "##high taxonomic-level filtering: \n" ;
########################
my$below = 0 ;
my$above = 0 ;
my$equal = 0 ;

foreach my$p (@Hphi) {
  if ($p > $empirical_data{"H"}) {$above ++ ;}
  elsif ($p < $empirical_data{"H"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"H"}, "\n" ;
print "horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"H"}, "\n" ;
print "horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"H"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@Mphi) {
  if ($p > $empirical_data{"M"}) {$above ++ ;}
  elsif ($p < $empirical_data{"M"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"M"}, "\n" ;
print "mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"M"}, "\n" ;
print "mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"M"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@Vphi) {
  if ($p > $empirical_data{"V"}) {$above ++ ;}
  elsif ($p < $empirical_data{"V"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"V"}, "\n" ;
print "vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"V"}, "\n" ;
print "vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"V"}, "\n\n" ;


########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HAphi) {
  if ($p > $empirical_data{"HA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HA"}, "\n" ;
print "aquatic horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HA"}, "\n" ;
print "aquatic horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MAphi) {
  if ($p > $empirical_data{"MA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MA"}, "\n" ;
print "aquatic mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MA"}, "\n" ;
print "aquatic mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VAphi) {
  if ($p > $empirical_data{"VA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VA"}, "\n" ;
print "aquatic vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VA"}, "\n" ;
print "aquatic vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HTphi) {
  if ($p > $empirical_data{"HT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HT"}, "\n" ;
print "terrestrial horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HT"}, "\n" ;
print "terrestrial horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MTphi) {
  if ($p > $empirical_data{"MT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MT"}, "\n" ;
print "terrestrial mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MT"}, "\n" ;
print "terrestrial mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VTphi) {
  if ($p > $empirical_data{"VT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VT"}, "\n" ;
print "terrestrial vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VT"}, "\n" ;
print "terrestrial vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MIphi) {
  if ($p) {
    if ($p > $empirical_data{"MI"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MI"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "mixed internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MI"}, "\n" ;
print "mixed internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MI"}, "\n" ;
print "mixed internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MI"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VIphi) {
  if ($p) {
    if ($p > $empirical_data{"VI"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VI"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "vertical internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VI"}, "\n" ;
print "vertical internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VI"}, "\n" ;
print "vertical internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VI"}, "\n\n" ;

##############################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HNphi) {
  if ($p > $empirical_data{"HN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HN"}, "\n" ;
print "nutrition functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HN"}, "\n" ;
print "nutrition functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HDphi) {
  if ($p) {
    if ($p > $empirical_data{"HD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "defense functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HD"}, "\n" ;
print "defense functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HD"}, "\n" ;
print "defense functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HMuphi) {
  if ($p) {
    if ($p > $empirical_data{"HMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HUkphi) {
  if ($p > $empirical_data{"HUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HUk"}, "\n" ;
print "unknown function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HUk"}, "\n" ;
print "unknown function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HMaphi) {
  if ($p) {
    if ($p > $empirical_data{"HMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMa"}, "\n" ;
print "manipulative functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMa"}, "\n" ;
print "manipulative functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMa"}, "\n\n" ;

##########################
########################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MNphi) {
  if ($p > $empirical_data{"MN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MN"}, "\n" ;
print "nutrition functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MN"}, "\n" ;
print "nutrition functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MDphi) {
  if ($p) {
    if ($p > $empirical_data{"MD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "defense functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MD"}, "\n" ;
print "defense functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MD"}, "\n" ;
print "defense functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MMuphi) {
  if ($p) {
    if ($p > $empirical_data{"MMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MUkphi) {
  if ($p > $empirical_data{"MUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MUk"}, "\n" ;
print "unknown functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MUk"}, "\n" ;
print "unknown functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MMaphi) {
  if ($p) {
    if ($p > $empirical_data{"MMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMa"}, "\n" ;
print "manipulative functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMa"}, "\n" ;
print "manipulative functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMa"}, "\n\n" ;

###############################
########################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VNphi) {
  if ($p > $empirical_data{"VN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VN"}, "\n" ;
print "nutrition functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VN"}, "\n" ;
print "nutrition functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VDphi) {
  if ($p) {
    if ($p > $empirical_data{"VD"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VD"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "defense functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VD"}, "\n" ;
print "defense functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VD"}, "\n" ;
print "defense functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VMuphi) {
  if ($p) {
    if ($p > $empirical_data{"VMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VUkphi) {
  if ($p > $empirical_data{"VUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VUk"}, "\n" ;
print "unknown functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VUk"}, "\n" ;
print "unknown functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VMaphi) {
  if ($p) {
    if ($p > $empirical_data{"VMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMa"}, "\n" ;
print "manipulative functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMa"}, "\n" ;
print "manipulative functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMa"}, "\n\n" ;

########################
#######################
print "\n\nlow taxonomic-level filtering \n" ;
#######################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;

foreach my$p (@Hplo) {
  if ($p > $empirical_data{"H"}) {$above ++ ;}
  elsif ($p < $empirical_data{"H"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"H"}, "\n" ;
print "horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"H"}, "\n" ;
print "horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"H"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@Mplo) {
  if ($p > $empirical_data{"M"}) {$above ++ ;}
  elsif ($p < $empirical_data{"M"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"M"}, "\n" ;
print "mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"M"}, "\n" ;
print "mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"M"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@Vplo) {
  if ($p > $empirical_data{"V"}) {$above ++ ;}
  elsif ($p < $empirical_data{"V"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"V"}, "\n" ;
print "vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"V"}, "\n" ;
print "vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"V"}, "\n\n" ;


########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HAplo) {
  if ($p > $empirical_data{"HA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HA"}, "\n" ;
print "aquatic horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HA"}, "\n" ;
print "aquatic horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MAplo) {
  if ($p > $empirical_data{"MA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MA"}, "\n" ;
print "aquatic mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MA"}, "\n" ;
print "aquatic mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VAplo) {
  if ($p > $empirical_data{"VA"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VA"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "aquatic vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VA"}, "\n" ;
print "aquatic vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VA"}, "\n" ;
print "aquatic vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VA"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HTplo) {
  if ($p > $empirical_data{"HT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HT"}, "\n" ;
print "terrestrial horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HT"}, "\n" ;
print "terrestrial horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MTplo) {
  if ($p > $empirical_data{"MT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MT"}, "\n" ;
print "terrestrial mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MT"}, "\n" ;
print "terrestrial mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VTplo) {
  if ($p > $empirical_data{"VT"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VT"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "terrestrial vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VT"}, "\n" ;
print "terrestrial vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VT"}, "\n" ;
print "terrestrial vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VT"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MIplo) {
  if ($p > $empirical_data{"MI"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MI"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "mixed internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MI"}, "\n" ;
print "mixed internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MI"}, "\n" ;
print "mixed internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MI"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VIplo) {
  if ($p > $empirical_data{"VI"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VI"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "vertical internal route:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VI"}, "\n" ;
print "vertical internal route:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VI"}, "\n" ;
print "vertical internal route:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VI"}, "\n\n" ;

##############################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HNplo) {
  if ($p > $empirical_data{"HN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HN"}, "\n" ;
print "nutrition functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HN"}, "\n" ;
print "nutrition functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HDplo) {
  if ($p > $empirical_data{"HD"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HD"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "defense functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HD"}, "\n" ;
print "defense functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HD"}, "\n" ;
print "defense functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HMuplo) {
  if ($p) {
    if ($p > $empirical_data{"HMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HUkplo) {
  if ($p > $empirical_data{"HUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HUk"}, "\n" ;
print "unknown function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HUk"}, "\n" ;
print "unknown function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HMaplo) {
  if ($p) {
    if ($p > $empirical_data{"HMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"HMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMa"}, "\n" ;
print "manipulative functions horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMa"}, "\n" ;
print "manipulative functions horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMa"}, "\n\n" ;

##########################
########################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MNplo) {
  if ($p > $empirical_data{"MN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MN"}, "\n" ;
print "nutrition functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MN"}, "\n" ;
print "nutrition functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MDplo) {
  if ($p > $empirical_data{"MD"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MD"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "defense functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MD"}, "\n" ;
print "defense functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MD"}, "\n" ;
print "defense functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MMuplo) {
  if ($p) {
    if ($p > $empirical_data{"MMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MUkplo) {
  if ($p > $empirical_data{"MUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MUk"}, "\n" ;
print "unknown functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MUk"}, "\n" ;
print "unknown functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MMaplo) {
  if ($p) {
    if ($p > $empirical_data{"MMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"MMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMa"}, "\n" ;
print "manipulative functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMa"}, "\n" ;
print "manipulative functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMa"}, "\n\n" ;

###############################
########################
########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VNplo) {
  if ($p > $empirical_data{"VN"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VN"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "nutrition functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VN"}, "\n" ;
print "nutrition functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VN"}, "\n" ;
print "nutrition functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VN"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VDplo) {
  if ($p > $empirical_data{"VD"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VD"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "defense functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VD"}, "\n" ;
print "defense functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VD"}, "\n" ;
print "defense functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VD"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VMuplo) {
  if ($p) {
    if ($p > $empirical_data{"VMu"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMu"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "multiple-functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VUkplo) {
  if ($p > $empirical_data{"VUk"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VUk"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "unknown functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VUk"}, "\n" ;
print "unknown functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VUk"}, "\n" ;
print "unknown functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VUk"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VMaplo) {
  if ($p) {
    if ($p > $empirical_data{"VMa"}) {$above ++ ;}
    elsif ($p < $empirical_data{"VMa"}) {$below ++ ;}
    else {$equal ++ ;}
  }
}
print "manipulative functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMa"}, "\n" ;
print "manipulative functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMa"}, "\n" ;
print "manipulative functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMa"}, "\n\n" ;



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
    print OUT "\tmixed\tvertical\n" ;
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

sub subsample_taxa {
  ## Find minimum number of samples per taxon and subsample all to that number
  my$taxon_level = $_[0] ;
  my%data = %{$_[1]} ;
  my@htaxa ;
  my@staxa ;
  my%new_dataset ;
  foreach my$col (0..8) {$new_dataset{$col} = () ;}
  my$colnum ;

  ## get taxa
  if ($taxon_level eq "high") {
    @htaxa = @{$data{0}} ;
    @staxa = @{$data{2}} ;
    $colnum = 2 ;
  }

  if ($taxon_level eq "low") {
    @htaxa = @{$data{1}} ;
    @staxa = @{$data{3}} ;
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
      foreach my$i (0..8) {push @{$new_dataset{$i}}, $data{$i}[$index] ;}
      delete $staxa{$data{$colnum}[$index]} ;
      delete $htaxa{$data{($colnum-2)}[$index]} ;
    }
    else {
    }
  }

  ##Unmask the following line to report the number of samples retained in subsamples
  print "taxon filter level: ", $taxon_level, "\t", "taxa included: ", scalar @{$new_dataset{1}}, "\n" ;
  return \%new_dataset ;
}
