#!/usr/bin/perl


##############################
#
# Generate_SNP_dataset.pl
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

use MasterList;

use FetchGenome;


die "Usage: perl Generate_SNP_dataset.pl <SAMPLES_LIST> <SNPFILE_DIR> <MASKED_REF_GENOME_FSA> <SKIP> <OUTFILE>\n" if @ARGV != 4;

# SNP outfiles must be named according to the format: Ref_v_Subject_out 


### declare global variables

my ($Strains, $snpsDir, $blockSize, $outFile) = @ARGV;

$StrainsHashRef = MasterList::strains($Strains);

%StrainsHash = %{$StrainsHashRef};

& SNP_ALLELES;

& REF_GENOME_SNPs;

& MASK_INDELS;

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

    #print "$S\n";

    open(FILE, "$snpsDir/$snpFile") || die "Problem\n";

    print "Reading $snpFile...\n";

    while($L = <FILE>) {

      next if $L =~ /scaf/;

      chomp($L);

      next if $L =~ /repeat/;

      @Data = split(/\t/, $L);

      if(@Data == 7) {

        ($B71Chr, $Other, $B71Pos, $OtherPos, $B71Nucl, $OtherNucl, $dir) = @Data;

        $B71Chr =~ s/.+?(\d)$/$1/;      # strip off everything except contig identifier (at end)

        $B71Chr = "Chr"."$B71Chr";      # add back a prefix to signal a B71 chromosome ID

        $indelHash{$B71Chr}{$B71Pos} = 1 if $OtherNucl eq '-';

        next if $B71Nucl !~ /^[AGTC]$/ || $OtherNucl !~ /^[AGTC]$/;

        ## record the presence of variant allele at each chromosome position

        $diffAllelesHash{$B71Ref}{$B71Pos} = 1

      }

    }

    close FILE  

  }

  close OUTFILES;

}


sub REF_GENOME_SNPs {

  ## Hash masked reference to allow grabbing of nucleotides at variant positions

  # NB: genome masked for all repeats and positions not aligned in any one isolate


  $maskedGenome = "Fully_filtered_alignments";

  $maskedGenomeHashRef = FetchGenome::getAlign($maskedGenome);

  $Seq =~ s/N//g;

  my $Length = length($Seq);

  %maskedGenomeHash = %$maskedGenomeHashRef;

## following lines for testing masking

#  my $Seq = $maskedGenomeHash{Chr2};

#  $Seq =~ s/N//g;

#  my $Length = length($Seq);  

#  print "PreLength: $Length\n";  
 
}

sub MASK_INDELS {

  print "Masking indel sites...\n";

  foreach my $chr (sort {$a cmp $b} keys %indelHash) {

    foreach my $pos (sort {$a <=> $b} keys %{$indelHash{$chr}}) {

      substr($maskedGenomeHash{$chr}, $pos-1, 1, 'N')

    }

  }

## following lines for testing masking

#  my $Seq = $maskedGenomeHash{Chr2};

#  $Seq =~ s/N//g;

#  my $Length = length($Seq);

#  print "PostLength: $Length\n";

}



sub READ_WRITE_SNPs {

  # read SNP reports again and examine genotypes at each possible SNP position:

  foreach my $File (@snpFilesList) {

    %WorkingHash = %maskedGenomeHash;

    next if($File !~ /out$/);

    ($Q, $S, $outsffx) = split(/_v_|_/, $File);			        # capture genome identifiers

    next unless exists($StrainsHash{$S});

    open(SNPs, "$snpsDir/$File") || die "Can't open SNPs file\n";

    #print "Query: $Q\tSubject: $S\n";

    while(my $L = <SNPs>) {

      chomp($L);

      my @SNPs = split(/\t/, $L);

      next if @SNPs != 7;

      ($qid, $sid, $qpos, $qend, $qnucl, $snucl, $dir) = @SNPs;

      next if $qid =~ /scaf/;

      ($ChromoNum = $qid) =~ s/.+?(\d+)$/$1/;

      next if $ChromoNum !~ /^[1-7]$/;					# inactivate for other projects

      $ChromoNum = 'Chr'.$ChromoNum;

      next if $qnucl !~ /^[AGTC]$/ || $snucl !~ /^[AGTC]$/;

      $maskedGenomeNucl = substr($WorkingHash{$ChromoNum}, $qpos-1, 1);

      substr($WorkingHash{$ChromoNum}, $qpos-1, 1, "$snucl") if $maskedGenomeNucl ne 'N';

    }

    close SNPs;
    
    & PRINT_GENOTYPE if $blockSize < 1;

    & PRINT_BLOCKS if $blockSize > 0;

  }

}

sub PRINT_GENOTYPE {

  print "Writing genotype for $S...\n";

  my $Genotype = undef;

  my $TotalLen = 0;

  open(OUT, '>>', "$outFile");

  print OUT ">$S\n";

  foreach my $chr (sort {$a cmp $b} keys %WorkingHash) {

    my $chromosomeBlock = $WorkingHash{$chr};

    $chromosomeBlock =~ s/N//g;

    $Genotype .= $chromosomeBlock;

  }

  print OUT "$Genotype\n";

  close OUT

}
    

sub PRINT_BLOCKS {

  # Retrieve genotype blocks of select length from masked genome

  # remove sites with non-aligned/repeated/missing data (Ns)

  # print blocks to separate files named by Chr#_block-start#

  foreach my $chr (sort {$a cmp $b} keys %WorkingHash) {

    for($i = 0; $i <= length($WorkingHash{$chr})-$blockSize; $i += $blockSize) {
     
      $blockStart = $i;

      open(OUT, '>>', "$chr"."_"."$blockStart");

      $block = substr($WorkingHash{$chr}, $blockStart, $blockSize);

      $block =~ s/N//g;

      $blockLength = length($block);

      print OUT "$S\n$block\n" if $blockLength >= 500;

      close OUT      

    }

  }

}

