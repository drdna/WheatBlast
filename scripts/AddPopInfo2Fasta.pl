#!/usr/bin/perl

open(L, $ARGV[0]);

while($L = <L>) {

  chomp($L);

  ($id, $pop) = split(/_/, $L);

  ($newid = $id) =~ s/-//g;

  $idHash{$newid} = $id;
 
  $Hash{$newid} = $pop;

}

close L;

open(F,  $ARGV[1]);

while($F = <F>) {

  chomp($F);

  if($F =~ />(.+)/ ) {

    $old = $1;

    unless (exists $Hash{$old}) {

      print "$F\n";

      next

    }

    $new = $idHash{$old}.'_'.$Hash{$old};

    next if $F =~ /_/;

    $F =~ s/>$old/>$new/;   

    print "$F\n"

  }


  else {

    print "$F\n";

  }

} 

close F;
