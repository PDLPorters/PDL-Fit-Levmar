use strict;
use warnings;
use PDL;
use PDL::Fit::Levmar;
use PDL::Fit::Levmar::Func;
use PDL::NiceSlice;
use Test::More;
use Test::PDL;

#  @g is global options to levmar
my @g = ( NOCOVAR => undef );

# used to check some return types to make sure computaton was float
sub check_type {
    my ($Type, @d) = @_;
    my $i=0;
    foreach ( @d )  {
	is $_->type, $Type, "type of var $i ".$_->info;
	$i++;
    }   
}

sub dimst {
    my $x = shift;
    return  "(" . join(',',$x->dims) . ")";
}

sub deb  { print STDERR $_[0],"\n" }
sub cpr  { print $_[0],"\n" }

cpr "# Test implicit threading over levmar()";
cpr "# Compiling fit function...";

# Need to use jacobian so fitting is more robust
my $Gf = '
       function
       x = p0 * exp( -t*t * p1);
       jacobian
       FLOAT ex, arg;
       loop
       arg = -t*t * p1;
       ex = exp(arg);
       d0 = ex;
       d1 = -p0 * t*t * ex ;
      ';

=pod
$Gf = '
       function
       x = p0 * exp( -t*t * p1);
      ';
=cut


my $Gh = levmar_func(FUNC=>$Gf);

cpr "# Done compiling fit function.";

sub keep_work_space {
    my ($Type) = @_;
    my $n = 100;
    my $t = zeroes($n)->xlinvals(-5,4.9)->convert($Type);
    my $x = zeroes($Type,$n);
    my $p = pdl($Type, 1,2);
    my $ip = pdl($Type, 3,4);
    $x .= $p((0)) * exp(-$t*$t * $p((1)) );
    my $h = levmar($ip,$x,$t,$Gh,@g);
    is_pdl $h->{P}, $p;
    check_type($Type, $h->{COVAR});
}

keep_work_space(double);
keep_work_space(float);

done_testing;
