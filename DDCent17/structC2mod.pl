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

   use File::Basename;
   use File::Copy ("cp", "mv");
   use File::Path;
   use File::Spec::Functions (":ALL");
   use feature qw(switch say);


    $, = " ";  $\ = "\n"; #set print defaults;
    select (STDERR); $|=1;      # Make unbuffered
    select (STDOUT); $|=1;      # Make unbuffered

my $sfx = '.inc';
open(In, "<", $ARGV[0]) || die "can not open $ARGV[0]";
 my($filename, $dirs, $suffix) = fileparse($ARGV[0],".h");
open(Out, ">", catfile($dirs,$filename.$sfx)) || die "can not open ".catfile($dirs,$filename.$sfx);

my $openmod  = "   MODULE $filename\n         USE ISO_C_BINDING";
my $closemod = "   END MODULE $filename";

$comnt = 0;
$typnam  = '';
$snam     = '';
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
    $_=~s|\[|\(|;   #change array open parenthesis
    $_=~s|\]|\)|;   #change array close parenthesis

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

# print the header before the first executable line
if($openmod) {print Out $openmod; undef($openmod);}
#  print $_;

  if(s|\s*}\s+(\w+?),(.*)||) { # } structname
#  if(s|\s*}\s+(\w+?),\s+.*?([A-Za-z]+).*||) { # } structname
  print "structure type",$1,"pointer type",$2;
    $typnam  = lc($1);
    ($pnam   = lc($2)) =~ s/\*//;
    if(!($snam =~ s/_.*/struct/)) {$snam .= 'struct';}
  print "global name",$snam;

    $struct[0] .= $typnam;
    {local $,="\n"; print Out @struct; @struct = ();}
    print Out  '         END TYPE',$typnam,"\n";

#    my $line = join('','       common /', $snam, "/");
    my $line = join('','         type (', $typnam, "), save,", $pnam);
#    for $n (@names) {
#      if(length($line) + length($n) > 69) {
#        print Out join(" ",$line, ' 'x(71-length($line)),"&");
#        $line = join('','     &   ',$n,',');
#      }
#      else {
#        $line = join('', $line," ",$n,',');
#      }
#    }
    $line = substr($line,0,-1);
    print Out $line;
#    print Out  '       SAVE /'.$snam.'/';
#    print Out  '       BIND(C) :: /'.$snam.'/',"\n\n";
  }

  elsif(m|typedef +struct|) {
    if($#struct) {
      local $, = "\n";
      print Out @struct;
    }
    @struct = ();
    @names  = ();
    push @struct, '      TYPE, BIND(C) :: '.$cc; # comment the fortran type
    next;
  }

  elsif(s/extern\s+(double|int|float) +//) {
    my $ctyp = $1;
    (my $vn = $_) =~ s/\s*;\s*$//;
        $vn =~ s/[()].*//;
        $vn =~ s/\s*//;

    s|(\s*=\s*){(.*)}|$1(/$2/)|; # remove array initialization
    if   ($ctyp eq 'int')    {print Out "     integer          (c_$ctyp), save, bind(C, name=\"$vn\") :: ".$_.$cc;}
    elsif($ctyp eq 'float')  {print Out "     real             (c_$ctyp), save, bind(C, name=\"$vn\") :: ".$_.$cc;}
    elsif($ctyp eq 'double') {print Out "     double precision (c_$ctyp), save, bind(C, name=\"$vn\") :: ".$_.$cc;}
  }

  elsif(s/\s*(double|int|float) +//) { # double swc[MAXLYR];
    my $ctyp = $1;
    m/(\w+)/; push (@names, $1);
    s|(\s*=\s*){(.*)}|$1(/$2/)|; # remove array initialization
    if   ($ctyp eq 'int')    {push @struct, "     integer          (c_$ctyp) :: ".$_.$cc;}
    elsif($ctyp eq 'float')  {push @struct, "     real             (c_$ctyp) :: ".$_.$cc;}
    elsif($ctyp eq 'double') {push @struct, "     double precision (c_$ctyp) :: ".$_.$cc;}
  }

  elsif(m/^\s*[{}]\s*$/) {} # do nothing

  else { print "oops:", $_;}
}

print Out $closemod;
close (Out);
exit;


# module flowstak
#
#   implicit none
#!  private                      ! hide everything except the function calls
#!  public flow, flowup, floclr, nflows
#   public flowdat, nflows, stackhead, ierr
#
#
#   type flowdat
#!      private
#!      integer                   :: when           !/* Time to flow */
#      real                      :: when            !/* Time to flow */
#!      character(LEN=16)         :: labl           !/* Flow frame label */
#      real,   pointer           :: donr => null()  !/* Source */
#      double precision, pointer :: ddonr => null() !/* Source */
#      real,   pointer           :: rcvr => null()  !/* Destination */
#      double precision, pointer :: drcvr => null() !/* Destination */
#      double precision          :: amt             !/* Amount */
#      type (flowdat), pointer   :: next => null()  !next entry
#   end type flowdat
#
#   integer                      :: ierr
#   integer                      :: nflows = 0      !/* Number of flows */
#   type (flowdat), POINTER      :: stackhead => null()
#   save nflows, stackhead, ierr
#
# contains
