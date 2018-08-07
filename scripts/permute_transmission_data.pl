use strict ;
use warnings ;
use Sort::Naturally ;
use List::Util qw/shuffle/;

##usage: perl permute_collect_transmission_data.pl dataset.txt > output_permutation_pvalues.txt

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

open SUBSAMPLEOUT, "> permuted_transmission_output.txt" or die "cannot open subsampledtaxa_transmission_output.txt\n" ;

foreach my$count (0..999) {
  my$randomized_data = randomize_data(\%symbioses) ;
  my%randomized_data = %{$randomized_data} ;
  my($horizontalra, $mixedra, $verticalra) = get_rates(\%randomized_data) ;
  my($hterrestrialra, $mterrestrialra, $vterrestrialra, $haquaticra, $maquaticra, $vaquaticra) = get_env(\%randomized_data) ;
  my($externalverticalra, $internalverticalra, $externalmixedra, $internalmixedra) = get_routes(\%randomized_data) ;
  my($horiznutra, $mixnutra, $vertnutra, $horizdefra, $mixdefra, $vertdefra, $horizmultra, $mixmultra, $vertmultra, $horizunkra, $mixunkra, $vertunkra, $horizmanipra, $mixmanipra, $vertmanipra) = get_function(\%randomized_data) ;

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
my$below = 0 ;
my$above = 0 ;
my$equal = 0 ;

foreach my$p (@Hp) {
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
foreach my$p (@Mp) {
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
foreach my$p (@Vp) {
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
foreach my$p (@HAp) {
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
foreach my$p (@MAp) {
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
foreach my$p (@VAp) {
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
foreach my$p (@HTp) {
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
foreach my$p (@MTp) {
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
foreach my$p (@VTp) {
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
foreach my$p (@MIp) {
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
foreach my$p (@VIp) {
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
foreach my$p (@HNp) {
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
foreach my$p (@HDp) {
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
foreach my$p (@HMup) {
  if ($p > $empirical_data{"HMu"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HMu"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "multiple-function horizontal mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"HMu"}, "\n" ;
print "multiple-function horizontal mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"HMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@HUkp) {
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
foreach my$p (@HMap) {
  if ($p > $empirical_data{"HMa"}) {$above ++ ;}
  elsif ($p < $empirical_data{"HMa"}) {$below ++ ;}
  else {$equal ++ ;}
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
foreach my$p (@MNp) {
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
foreach my$p (@MDp) {
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
foreach my$p (@MMup) {
  if ($p > $empirical_data{"MMu"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MMu"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "multiple-functions mixed mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"MMu"}, "\n" ;
print "multiple-functions mixed mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"MMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@MUkp) {
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
foreach my$p (@MMap) {
  if ($p > $empirical_data{"MMa"}) {$above ++ ;}
  elsif ($p < $empirical_data{"MMa"}) {$below ++ ;}
  else {$equal ++ ;}
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
foreach my$p (@VNp) {
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
foreach my$p (@VDp) {
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
foreach my$p (@VMup) {
  if ($p > $empirical_data{"VMu"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VMu"}) {$below ++ ;}
  else {$equal ++ ;}
}
print "multiple-functions vertical mode:\t", $above, "(", $above/1000*100, "%)", " above the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $below, "(", $below/1000*100, "%)", " below the empirical value of ", $empirical_data{"VMu"}, "\n" ;
print "multiple-functions vertical mode:\t", $equal, "(", $equal/1000*100, "%)", " equal to the empirical value of ", $empirical_data{"VMu"}, "\n\n" ;

########################
$below = 0 ;
$above = 0 ;
$equal = 0 ;
foreach my$p (@VUkp) {
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
foreach my$p (@VMap) {
  if ($p > $empirical_data{"VMa"}) {$above ++ ;}
  elsif ($p < $empirical_data{"VMa"}) {$below ++ ;}
  else {$equal ++ ;}
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

sub randomize_data {
  ## Shuffle the values of every characteristic (column), randomizing the association between characteristics
  my%data = %{$_[0]} ;
  my%new_dataset ;

  foreach my$col (0..8) {
    my@shuffled_data = shuffle(@{$data{$col}}) ;
    $new_dataset{$col} =  \@shuffled_data ;
  }

  return \%new_dataset ;
}
