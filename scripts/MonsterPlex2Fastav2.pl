#!/usr/bin/perl

# MonsterPlex2Fasta.pl

# written by Mark L. Farman

# i) Reads raw vcf files produced by GATK and identify genotypes at specific targeted sites (one site per locus)

# ii) Records positions and read coverage for new variants that surface upon analysis of additional strains

# iii) Updates bedCoverage file for all variants (targeted and new)



die "Usage: perl MonsterPlex2Fasta_noMGG.pl <ref-genome> <VCF-directory> <bed-file> <outdir>\n" if @ARGV < 4;

# load modules

use bedCoverage;
use FetchGenome;


# ARGUMENTS

($ref, $vcf, $bedfile, $outdir) = @ARGV;


# MAIN PROGRAM

mkdir $outdir;

& READ_REF_GENOME;

& COVERAGE;

& READ_SNP_SITES;

& READ_VCFs;

& PRINT_GENOTYPES;

& FIND_NEW_SNP_SITES;

& PRINT_NEW_SNP_SITES;

& PRINT_DEPTHS;

& CHECK_PRINT_ORDER;

exit();


## SUBROUTINES

sub READ_REF_GENOME {

# Read reference genome to a hash

  $GenomeRef = FetchGenome::getSeqs($ref);
  %Genome = %{$GenomeRef};
  #print "$Genome{1}\n";
}


sub COVERAGE {

# Keep track of coverage at variant positions that were not originally targeted by assay

  $bedCovRef = bedCoverage::bedCov($bedfile);
  %bedCoverage = %{$bedCovRef};
}


sub READ_SNP_SITES {

# Read list of targeted SNP sites (this can be updated as new, confident sites are identified)

  open(MPS, "MonsterPlexSites") || die "Can't read targeted sites file\n";
  while($L = <MPS>) {
    chomp($L);
    ($Chr, $refID, $pos, $MGG) = split(/\t/, $L);
  #  print "($Chr, $refID, $pos, $MGG)\n";
    $Chr =~ s/chr0//; 
    $SitesHash{$Chr}{$pos} = 1;
  }
  close MPS;
}


sub READ_VCFs {

  opendir(VCFDIR, $vcf);
  @VCFDIR = readdir(VCFDIR);
  foreach $VCFfile (@VCFDIR) {
    next if $VCFfile =~ /^\./;
    & PROCESS_VCF
  }
  closedir VCFDIR
}


sub PROCESS_VCF {

# read VCF files, record coverage across all sites and score genotypes at variant sites

  open(VCF, "$vcf/$VCFfile") || die "can't open VCF file: $VCFfile\n";
  ($SeqID = $VCFfile) =~ s/\.vcf$//g;
  while($V = <VCF>) {
    chomp($V);
    @VCF = split(/\t/, $V) if $V =~ /^supercont/;
    @GT = split(/:/, $VCF[9]) if $V =~ /^supercont/;

    ## Record sites with >5x coverage that contain ref alleles

    if($GT[1] =~ /^(\d+)$/ && $1 > 5) {		# only look at 5x covered sites with just ref alleles
      $depth = $1;
      ($Chr = $VCF[0]) =~ s/supercont8\.//;
      next if $Chr > 7;
      $start = $VCF[1];
      ($end = $VCF[7]) =~ s/END=//;
      for($i = $start; $i <= $end; $i++) {		# iterate through alignment block and add deeply covered sites to %RefHash
        $refAllele = substr($Genome{$Chr}, $i-1, 1);	# retrieve ref allele from reference genome
        $RefHash{$SeqID}{$Chr}{$i} = $refAllele;		# store ref allele for those sites
        $Depth{$Chr}{$i}{$SeqID} = $depth if exists ($SitesHash{$Chr}{$i})	# Record depth of coverage (this is primarily for non-specifically amplified loci)
      }
    }

    ## check variant genotype blocks

    elsif($GT[1] =~ /,/) {   
      ($a, $b, $c) = split(/,/, $GT[1]);
      if($b > 7) {							# only consider sites with coverage > 7 alternate allele reads
        if($a == 0) {
          ($VariantHash{$SeqID}{$Chr}{$VCF[1]} = $VCF[4]) =~ s/,.+//;	# record alt allele if all reads are alt
          $Depth{$Chr}{$VCF[1]}{$SeqID} = $b				# record coverage over variant site
        }
        elsif($b/$a > 7) {						# evaluate sites with some ref alleles present; only capture if #alt/ #ref > 7
          ($alt = $VCF[4]) =~ s/,.+//; 
          if(length($alt) > 1) {					# score indels as missing site
            $VariantHash{$SeqID}{$Chr}{$VCF[1]} = 'N';
          }
          else {
            $VariantHash{$SeqID}{$Chr}{$VCF[1]} = $alt;	                # record alt if	#alt / #ref > 7
            $Depth{$Chr}{$VCF[1]}{$SeqID} = $b	                        # record coverage over site
          }
        }
      }
    } 									# end of checking variant genotype block
  }
  close VCF;
}


sub PRINT_GENOTYPES {

# For every strain, loop through variant hash looking for presence of each targeted site
# if site is listed in variant hash, call the relevant alt allele
# elsif site was interrogated but recorded as ref, call ref allele
# else call as 'N'

  my $i = 0;
  while(-f "$outdir/MPgenotypes".$i.".fasta") {
    $i += 1;
  }
  open(GENOTYPES, '>', "$outdir/MPgenotypes".$i.".fasta");
  foreach $SeqID (sort {$a cmp $b} keys %VariantHash ) {
    foreach my $chr (sort {$a <=> $b} keys %SitesHash) {
      foreach my $pos (sort {$a <=> $b} keys %{$SitesHash{$chr}}) {
        ## add genotypes to sequence string
        if(exists($VariantHash{$SeqID}{$chr}{$pos})) {
          $Seq .= $VariantHash{$SeqID}{$chr}{$pos}
        }
        elsif(exists($RefHash{$SeqID}{$chr}{$pos})) {
          $Seq .= $RefHash{$SeqID}{$chr}{$pos};
        }
        else {
          $Seq .= 'N'
        }
      }
    }
    $Seq =~ s/\*/N/g;		# needed for a few genotypes that GATK calls as '*'
    print GENOTYPES ">$SeqID\n$Seq\n";
    $Seq = ''; 
  }
}


sub FIND_NEW_SNP_SITES {

# Iterate through Depth hash to find sites that are adequately covered for new variant discovery
# Then iterate through Variant hash to find the new variants

  my $Depth = 0;
  foreach my $chr (sort {$a <=> $b} keys %Depth) {
    foreach my $pos (sort {$a <=> $b} keys %{$Depth{$chr}}) {
      foreach my $strainID (sort {$a <=> $b} keys %{$Depth{$chr}{$pos}}) {
        $Depth += $Depth{$chr}{$pos}{$strainID};
      }
      SITES($Depth)
    }
  }
}
 
sub SITES {
  my $Depth = $_[0];
  foreach my $strainID (sort {$a cmp $b} keys %VariantHash) {
    foreach my $chr (sort {$a <=> $b} keys %{$VariantHash{$strainID}}) {
      foreach my $pos (sort {$a <=> $b} keys %{$VariantHash{$strainID}{$chr}}) {       
        next if exists $SitesHash{$chr}{$pos};
#        print join ("\t", ($chr, $pos)), "\t", join(" ", @strainIDs), "\t", $Depth, "\n" if $Depth > 10000;
        push @{$NewSNPsHash{$chr}{$pos}}, $strainID if $Depth > 10000;		# Add site to a hash recording new SNPs if depth of coverage across all samples is > 50,000
      }
    }
  }
}    


sub COUNT_STRAINS {

# Count_strains that exhibit variants at previously untargeted sites

  foreach my $chr (sort {$a <=> $b} keys %NewSNPsHash) {
    foreach my $pos (sort {$a <=> $b} keys %{$NewSNPsHash{$chr}}) {
      @numStrains = @{$NewSNPsHash{$chr}{$pos}};
      $numStrains = @numStrains;
#      print "Chr$chr\t$pos\t$numStrains\t";
#      print join (", ", @{$NewSNPsHash{$chr}{$pos}}), "\n";
    }
  }

}


sub PRINT_NEW_SNP_SITES {

# print extra sites with deep coverage

  my $i = 0;
  while(-f "$outdir/MP_new_sites".$i.".txt") {
    $i += 1;
  }
  open(NEW_SITES, '>', "$outdir/MP_new_sites".$i.".txt");
  print NEW_SITES "chromosome\treference\tposition\tMGG_identifier\n";
    foreach my $chr (sort {$a <=> $b} keys %NewSNPsHash) {
    foreach my $pos (sort {$a <=> $b} keys %{$NewSNPsHash{$chr}}) {
      next if exists($SitesHash{$chr}{$pos});
      unless(exists($redundant{$chr}{$pos})) {
        print NEW_SITES "chr0$chr\tNC_017850.1\t$pos\tMGG\n" if exists($bedCoverage{$chr}{$pos});
        $redundant{$chr}{$pos} = 1
      }
    }
  }
}


sub PRINT_DEPTHS {
 
# print depths for plotting using R script
# this only works for SNPs in the current SNPs list

  my $i = 0;
  while(-f "$outdir/MP_target_depthsDF".$i.".txt") {
    $i += 1;
  }
  open(TDDF, '>', "$outdir/MP_target_depthsDF".$i.".txt");
  print TDDF "Site\tPos\tStrain\n";
  foreach $chr (sort {$a <=> $b} keys %Depth) {
    foreach my $pos (sort {$a <=> $b} keys %{$Depth{$chr}}) {
      $Site += 1;
      foreach my $strainID (sort {$a <=> $b} keys %{$Depth{$chr}{$pos}}) {
        if(exists($SitesHash{$chr}{$pos})) {
          print TDDF "$Site\t$Depth{$chr}{$pos}{$strainID}\t$strainID\n"
        }
      }
    }
  }
} 


    
sub CHECK_PRINT_ORDER {

# checks that sites get printed in Chr,pos order

  foreach my $chr (sort {$a <=> $b} keys %SitesHash) {
    foreach my $pos (sort {$a <=> $b} keys %{$SitesHash{$chr}}) {
      print "$chr\t $pos\n";
    }
  }
}
