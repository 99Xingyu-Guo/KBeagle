use strict;
use warnings;

# 
# Conver vcf to hpm
# demo usage:
# perl vcf2hpm2.pl 100k.vcf demo
#			

my $fichier_in;
my $dt ; 

  # part one: check argument
  if ($#ARGV == -1) {
   print " <Name  of input VCF>    <Name  of output VCF>  conditional value Ancestral Allele field in VCF <0: n oAA field, 1 AA field>  ? "; 
   exit;
  }
  else {
   $fichier_in = $ARGV[0];
   $dt = $ARGV[1];
  }

  # part two: lots of variables
  my ($i, $j);
  my $k;
  my $ref; #reference allele
  my $alt; #alternative allele
  my $o;
  my $p; 

  # part three: store vcf into AoA1
  open (IN,"<$fichier_in")  or die "Impossible to open file $fichier_in because : $!\n";  
  
  my $Chrant;
  my $status = 0;
  my ($h1,$h2); #chr_name position
  my @header;
  my $line=1;

while (<IN>) { 
  next if m/^##/;
  chomp $_;
  my @AoA1 = split(/\t/,$_);

  if(/#CHR/){
    #1...9: annotation
    #10...end: individuals
    my $nsouche = $#AoA1-8;
    print "Number of individuals  $nsouche\n";  

    # prepare header line
    @header = (@AoA1);

    # For Chrant's initial value
    $status = 1;

    next;
  }

######
######

  #We should give $Chrant a initial value;
  if($status == 1){
    $Chrant = $AoA1[0];

    # print header line; 
    my $dt_out = $dt.$AoA1[0].".txt";
    open(OUT,">$dt_out") or die "Erreur de creation de DT format hapmap";
    print OUT "rs# SNPalleles $header[0] $header[1] Strand Genome_build Center ProtLSID assayLSID panelLSID QC_code S288C ";
    for $j ( 9..$#header) {  print OUT "$header[$j] ";}
    print OUT "\n";


    $status = 2;
  }

  # part five: deal with vcf file
  #Check for  chromosome change , if yes makes a new file 
  # 0                      3                        6      7
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  qh100
     
  #different chromosome to different files
  if ($Chrant ne $AoA1[0]) {
    close OUT;
    my $dt_out=   $dt.$AoA1[0].".txt"; 

    open(OUT,">$dt_out") or die "Error creation output DT  $dt_out";
    print OUT "rs# SNPalleles $header[0] $header[1] Strand Genome_build Center ProtLSID assayLSID panelLSID QC_code S288C ";
    for $j (9..$#header) {  print OUT "$header[$j] ";}
    print OUT "\n";
    $Chrant = $AoA1[0];
  }

  # obtain ref and alt alleles
  $ref  = $AoA1[3];
  $alt  = $AoA1[4];

  print OUT "rs$line $ref\/$alt $AoA1[0] $AoA1[1] $AoA1[5] $AoA1[6] unknown urn\:unknown urn\:unknown urn\:unknown QC+ $ref$ref ";
  $line += 1;

  # deal with each individuals
  for $k (9..$#AoA1){  
    my %mapping;
    &prepare_mapping(\%mapping,\@AoA1);
    if($AoA1[$k] =~ /(.)\|(.)/){
      die if (not exists($mapping{$1})) or (not exists($mapping{$2}));
      print OUT "$mapping{$1}$mapping{$2} ";
    }elsif($AoA1[$k] =~ /(.)\/(.)/){
      die if (not exists($mapping{$1})) or (not exists($mapping{$2}));
      print OUT "$mapping{$1}$mapping{$2} ";
    }else{die "strange genotype!\n@AoA1\n"}
  } 

  print OUT "\n";
  #####
  #####
}   
close IN;

    
close OUT; 
print  "That's all!\n";	         
exit;

sub prepare_mapping{
  my ($mapping,$AoA1) = @_;
  ${$mapping}{'.'} = 'N';
  ${$mapping}{0} = @{$AoA1}[3];
  my @alt_arr = split /,/,@{$AoA1}[4];
  for(my $i = 0;$i <= $#alt_arr;$i ++){
    my $index = $i + 1;
    ${$mapping}{$index} = $alt_arr[$i];
  }
}

