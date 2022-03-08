#!/usr/local/bin/perl -w
#

=head1 description

/*****************************************************************************
**
**  FILE:      structC2for.pl
**
**  PURPOSE:   Convert C struct definition to Fortran 2003
**             fortran-iso-c-binding compatible common block.
**
*****************************************************************************/
=cut

    $, = " ";  $\ = "\n"; #set print defaults;
    select (STDERR); $|=1;      # Make unbuffered
    select (STDOUT); $|=1;      # Make unbuffered

open(In, "<", $ARGV[0]) || die "can not open $ARGV[0]";
$ARGV[0] =~ s/\.h/.inc/;
open(Out, ">", $ARGV[0]) || die "can not open $ARGV[0]";

$comnt = 0;
$snam  = '';
@struct = ();
while($_=<In>) {$_=~s/\s*$//;
         if (!m/\S/) {print Out; next;} # print blank lines
         last if (m/^\s*EXIT/);  # terminate count on EXIT flag
#         last if (m/^\s*EXIT/);  # terminate count on EXIT flag
#         s/[#!].*//;             # remove everything after a comment character
         s/[#!].*//;             # remove preprocessor commands
         next if (!m/\S/);       # skip blank (preprocessor) lines
#         $_=~ s/['"]//g;         # remove string quotes


    $_=~s|\]\[|,|g; #change multiple dimensions
    $_=~s|\[|\(|; #change array open parenthesis
    $_=~s|\]|\)|; #change array close parenthesis

    my $cc = '';

  if($comnt) {          # in multi-line comment
    $comnt =0 if(s|\*\/||); # clear the comment flag at the end of the comment
    if ($_=~m/\S/) {
      if(@struct) {push(@struct,'!'.$_);}
      else {print Out '!'.$_;} # if ($_=!m/\S/);
    }
    next;
  }
  elsif(s|(\s*)\/\*(.*)||) { # start a comment line
    $comnt =1; # set comment flag
    $cc = join('',$1,'!',$2);
    if($cc =~ s|\s*\*\/||) {
       $comnt =0; # clear comment flag at the end of the comment
    }
#    if((length($_) .eq. 0) or ($_ =~ m/^\s*/$)) {
    if(length($_) == 0) {
      if(@struct) {push @struct, $cc;} else {print Out $cc;}
      next;
    }
  }

  if(s|\s*}\s+(\w+?),(.*)||) { # } structname
#  if(s|\s*}\s+(\w+?),\s+.*?([A-Za-z]+).*||) { # } structname
  print "structure type",$1,"pointer type",$2;
    $snam  = lc($1);
    if(!($snam =~ s/_.*/struct/)) {$snam .= 'struct';}
  print "global name",$snam;
    $struct[0] .= $snam;
    {local $,="\n"; print Out @struct; @struct = ();}
    print Out  '!      END TYPE',$snam,"\n";

    my $line = join('','       common /', $snam, "/");
    for $n (@names) {
      if(length($line) + length($n) > 69) {
        print Out join(" ",$line, ' 'x(71-length($line)),"&");
        $line = join('','     &   ',$n,',');
      }
      else {
        $line = join('', $line," ",$n,',');
      }
    }
    $line = substr($line,0,-1);
    print Out $line;
    print Out  '       SAVE /'.$snam.'/';
    print Out  '       BIND(C) :: /'.$snam.'/',"\n\n";
  }

  elsif(m|typedef +struct|) {
    if($#struct) {
      local $, = "\n";
      print Out @struct;
    }
    @struct = ();
    @names  = ();
    push @struct, '!      TYPE, BIND(C) :: '.$cc; # comment the fortran type
    next;
  }

  elsif(s/\s*(double|int|float) +//) { # double swc[MAXLYR];
    my $ctyp = $1;
    m/(\w+)/; push (@names, $1);
    s|(\s*=\s*){(.*)}|$1(/$2/)|; # remove array initialization
    if(substr($ctyp,0,3) =~ /int/i) {
      push @struct, "         integer(c_$ctyp) :: ".$_.$cc;
    } else {
      push @struct, "         real(c_$ctyp) :: ".$_.$cc;
    }
  }

#  elsif(s|\s*float +||) { # float  dpthmn[MAXLYR];
#    m/(\w+)/; push (@names, $1);
#    push @struct, '         real(c_FLOAT) :: '.$_.$cc;
#  }
#
#  elsif(s|\s*int +||) { # int    numlyrs;
#    m/(\w+)/; push (@names, $1);
#    push @struct, '         integer(c_int) :: '.$_.$cc;
#  }
}
exit;
