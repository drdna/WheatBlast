#!/usr/bin/perl

##############################
#
# Generate_STRUCTURE.pl
#
# written by Mark L. Farman
#
# Purpose: Read SNPcaller outfiles, score all SNP positions, and record haplotypes for all individuals in a specified list
#
# Note: Must be re-run every time a new strain is added (unavoidable)
#
##############################

#use strict;

#use warnings;

use FetchGenome;

die "Usage: perl Generate_STRUCTURE.pl <SAMPLES_LIST> <SNPFILE_DIR> <REF_GENOME_FASTA> <SKIP>\n" if @ARGV != 4;

# SNP outfiles must be named according to the format: Ref_v_Subject_out 

use MasterList;

### declare global variables

my ($Strains, $snpsDir, $refGenome, $nth) = @ARGV;

$ChrLensref = FetchGenome::getLens($refGenome);

%ChrLens = %{$ChrLensref};

@KeysList = keys %ChrLens;

print "@KeysList\n";

print "$ChrLens{$KeysList[0]}\n";

foreach my $key (keys %ChrLens) {

  $variants{$key} = "0" x $ChrLens{$key}

}

$StrainsHashRef = MasterList::strains($Strains);

%StrainsHash = %$StrainsHashRef;

& SNP_ALLELES;

& REF_GENOME_SNPs;

& MAKE_SNP_LISTs;

& MAKE_HEADERs;

& READ_WRITE_SNPs;


## SUBROUTINES


sub SNP_ALLELES {

  # Loop through all SNP files and record positions of all SNPs relative to B71 reference genome

  $snpsDir =~ s/\/$//;

  opendir(OUTFILES, $snpsDir);

  @snpFilesList = readdir(OUTFILES);

  print "Identifying all relevant SNP loci\n";

  foreach $snpFile (@snpFilesList) {

    next unless $snpFile =~ /out$/;

    my($Q, $S) = $snpFile =~ /(.+)_v_(.+)_out/;

    next unless exists($StrainsHash{$S});

    print "$S\n";

    open(FILE, "$snpsDir/$snpFile") || die "Problem\n";

    while($L = <FILE>) {

      next if $L =~ /scaf/;

      chomp($L);

      @Data = split(/\t/, $L);

      if(@Data == 7) {

        ($B71Ref, $Other, $B71Pos, $OtherPos, $B71Nucl, $OtherNucl, $dir) = @Data;

        next if $B71Nucl !~ /^[AGTC]$/ || $OtherNucl !~ /^[AGTC]$/;

        $B71Ref =~ s/.+(\d)$/$1/;	# strip off everything except contig identifier (at end)

        ## record number of alleles at each variant position

        $diffAllelesHash{$B71Ref}{$B71Pos}{$OtherNucl}++;

      }

    }

    close FILE  

  }

  close OUTFILES;

}



sub REF_GENOME_SNPs {

  ## grab nucleotides at variant positions of the reference genome

  $GenomeHashRef = FetchGenome::getSeqs($refGenome);

  %GenomeHash = %$GenomeHashRef;

}



sub MAKE_SNP_LISTs {

  # if SNP is informative (i.e. > 1 SNP at that position), record number-based allele at corresponding position of string in variants hash

  # and create a hash of snp positions by chromosome

  %baseNums = (A, "1", a, "1", T, "2", t, "2", C, "3", c, "3", G, "4", g, "4");

  foreach my $chr (sort {$a cmp $b} keys %diffAllelesHash) {

    foreach my $pos (sort {$a <=> $b} keys %{$diffAllelesHash{$chr}}) {

      $multiAllelic = 'no';

      foreach my $allele (sort {$a cmp $b} keys %{$diffAllelesHash{$chr}{$pos}}) {

        $numAlleles = $diffAllelesHash{$chr}{$pos}{$allele};

        $multiAllelic = 'yes' if $numAlleles > 1 && $numAlleles < 94

      }

      if ($multiAllelic eq 'yes') {

        push @{$snpsByChromosome{$chr}}, $pos;

        $refBase = substr($GenomeHash{$chr}, $pos-1, 1);

        $refNum = $baseNums{$refBase};

        $thisChrLen = length($variants{$chr});

        print "$chr, $pos, $refNum, $thisChrLen\n";

        substr($variants{$chr}, $pos-1, 1, $refNum);

      }

    }

  }

}


sub MAKE_HEADERs {

  $outfile = $snpsDir."_structure$nth";

  open(OUT, '>', $outfile) || die "Can't create outfile: $outfile\n";

  foreach my $chr (sort {$a cmp $b} keys %snpsByChromosome) {

  $lastIndex = @{$snpsByChromosome{$chr}} - $nth - 1;

   for($i = 0; $i <= $lastIndex; $i += $nth) {

      push @tempHeader1, ${$snpsByChromosome{$chr}}[$i];

      push @tempHeader2, ${$snpsByChromosome{$chr}}[$i]-${$snpsByChromosome{$chr}}[$i-$nth];

    }

    $tempHeader2[0] = -1;

    push @Header1, @tempHeader1;

    push @Header2, @tempHeader2;

    @tempHeader1 = ();

    @tempHeader2 = ();

  }

  unshift @Header1, '';

  unshift @Header2, '';

  print OUT join ("\t", @Header1), "\n";

  print OUT join ("\t", @Header2), "\n";

  close OUT

}


sub READ_WRITE_SNPs {

  # read SNP reports again and record genotypes at each possible SNP position:

  foreach my $File (@snpFilesList) {

    %subjSNPsHash = undef;

    @genotypeArray = ();

    next if($File !~ /out$/);

    ($Q, $S, $outsffx) = split(/_v_|_/, $File);			        # capture genome identifiers

    next unless exists($StrainsHash{$S});

    open(SNPs, "$snpsDir/$File") || die "Can't open SNPs file\n";

    print "Query: $Q\tSubject: $S\n";

    while(my $L = <SNPs>) {

      chomp($L);

      my @SNPs = split(/\t/, $L);

      ($qid, $sid, $qpos, $qend, $qnucl, $snucl, $dir) = @SNPs;

      next if @SNPs != 7;

      next if $qid =~ /scaf/;

      ($ChromoNum = $qid) =~ s/.+(\d)$/$1/;

      next if $ChromoNum !~ /^[1-7]$/;

      next if $qnucl !~ /^[AGTC]$/ || $snucl !~ /^[AGTC]$/;

      $subjSNPsHash{$ChromoNum}{$qpos} = $snucl;

    }

    close SNPs;

    & PRINT_GENOTYPE

  }

}

sub PRINT_GENOTYPE {

  open(OUT, '>>', $outfile) || die "Can't create outfile: $outfile\n";

  $alignFile = "B71v2sh.$S".'_alignments';

  my $chrAlignRef = READ_ALIGN_FILE($alignFile);

  my %chrAlign = %{$chrAlignRef};

  ## read all positions in the snpsByChromosome lists and check if they occur in subjectSNPsHash

  ## if they do, assign alt base

  ## if not and they align uniquely, assign ref base

  ## if there is not a unique alignment (0 or >1) assign -9
  
  foreach my $chr (sort {$a cmp $b} keys %snpsByChromosome) {

    

    $alignString = $chrAlign{$chr};

#   $subAlign = substr($alignString, 150, 50);

#   print "$subAlign\n";

    print "warning\n" if length($alignString) == 0;

    $lastIndex = @{$snpsByChromosome{$chr}} - $nth - 1;

    for($i = 0; $i <= $lastIndex; $i += $nth) {

      $snpPosition = ${$snpsByChromosome{$chr}}[$i];

      if(exists($subjSNPsHash{$chr}{$snpPosition})) {

       push @genotypeArray, $baseNums{$subjSNPsHash{$chr}{$snpPosition}};	# for numerals

#      push @genotypeArray, $subjSNPsHash{$chr}{$snpPosition};		 	# for bases

      }

      elsif(substr($alignString, $snpPosition-1, 1) == 1) {

        my $refBase = substr($variants{$chr}, $snpPosition-1, 1);		# for numerals

#       my $refBase = substr($GenomeHash{$chr}, $snpPosition-1, 1); 		# for bases

#       $refRegion = substr($variants{$chr}, $snpPosition-3, 5);

        push @genotypeArray, $refBase;
  
      } 

      else {

        push @genotypeArray, '-9';						# for numerals

#        push @genotypeArray, 'N';						# for bases

      }

    }

  }

  unshift @genotypeArray, $S;

  print OUT join ("\t", @genotypeArray), "\n";

  @genotypeArray = ();

  close OUT;

}


sub READ_ALIGN_FILE {
    
  $alignFile = $_[0];

  open(ALIGN, "ALIGN_STRINGS/$alignFile") || die "Can't open alignment file: $alignFile\n";

  %chrAlign = undef;

  while($alignRecord = <ALIGN>) {

    if($alignRecord =~ /(Chr\d)/) {

      chomp($alignRecord);

      my ($chr, $alignString) = split(/\t/, $alignRecord);

      $chr =~ s/Chr(\d)/$1/;

      $chrAlign{$chr} = $alignString

    }

  }

  close ALIGN;

  return(\%chrAlign)

}

