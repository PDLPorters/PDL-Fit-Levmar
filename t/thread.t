use strict;
use warnings;
use PDL;
use PDL::Fit::Levmar;
use PDL::Fit::Levmar::Func;
use PDL::NiceSlice;
use Test::More;
use Test::PDL -atol => 1e-4;

# Check pdl 'threading'. That is, automatically looping over
# extra dimensions in pdls

#  @g is global options to levmar
my @g = ( NOCOVAR => undef );

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

note "Test implicit threading over levmar()";
note "Compiling fit function...";

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

# there is a big difference in speed here!

my $Gh = levmar_func(FUNC=>$Gf);

note "Done compiling fit function.";

# Thread x. Try the same parameters on different sets of data.
# Also test workspace allocation.
sub thread1 {
    my ($Type) = @_;
    my $n = 10000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
    my $params =  [ [3,.2], [ 9, .1] , [2,.01], [3,.3] ];
    my $x = zeroes($Type,$n,scalar @$params);
    my $i = 0;
    map {  $x(:,$i++)  .= $_->[0] * exp(-$t*$t * $_->[1]  ) }  @$params;
    my $p = pdl $Type, [ 5, 1]; # starting guess
    check_type($Type, $p,$x,$t);
    my $h = levmar(  $p, $x, $t, $Gh, @g, DERIVATIVE => 'numeric');
    check_type($Type, $h->{INFO});
    is_pdl $h->{P}, pdl($Type, $params), "Thread x, 1 thread dim ($Type)"
        or diag "report=", levmar_report($h);
    my $m = 2;
    my $s = 4*$n+4*$m + $n*$m + $m*$m;
    $h = levmar(  $p, $x, $t, $Gh, @g);
    $h = levmar(  $p, $x, $t, $Gh, @g);
    $s = 2*$n+4*$m + $n*$m + $m*$m;
}

# Change the following routines to use map the same way

# Thread p. Not the right expression, I think.
# ie, try multiple parameters a single data set.
sub thread2 {
    my ($Type) = @_;
    my $n = 10000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
    my $x = zeroes($Type, $n);
    my $params =   [[0,3,.2]]; # only 1 dimension
    map {  $x(:,$_->[0])  .= $_->[1] * exp(-$t*$t * $_->[2]  ) }  @$params;
    my $p = pdl $Type, [ [ 5, 1], [ 2,4] ]; # starting guess
    my $outp = pdl ($Type, $params);
    my $correct =  pdl $Type, [$outp(1:2,(0)), $outp(1:2,(0))]; #x Ugly
    check_type($Type, $p,$x,$t);
    my $h = levmar($p, $x, $t, $Gh, @g);
    check_type($Type, $h->{INFO});
    is_pdl $h->{P}, $correct, "Thread p, 1 thread dim";
}

# This one threads over both p and x, with one
# extra dimension
sub thread3 {
    my ($Type) = @_;
    my $n = 10000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
    my $params =  [ [0,3,.2], [1, 2, .1] ];
    my $x = zeroes($Type, $n,scalar(@$params));
    my $res =  pdl $Type, $params;
    map {  $x(:,$_->[0])  .= $_->[1] * exp(-$t*$t * $_->[2]  ) }  @$params;
    my $p = pdl $Type, [ [ 5, 1], [2,4]] ; # starting guess
    check_type($Type, $p,$x,$t);
    my $h = levmar($Gh , $p,$x,$t,  @g );
    check_type($Type, $h->{INFO});
    is_pdl $h->{P}, $res(1:2,:), "Thread both x and p, 1 thread dim";
}

sub thread4 {
    my ($Type) = @_;
    my $n = 1000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
# Put any number of triples of actual parameters here.
    my $params =  [ [0,3,.2], [1, 28, .1] , [2,2,.01], [3,3,.3] ];
    my $nx = scalar(@$params);
    my $x = zeroes($Type, $n,$nx);
    my $res =  pdl $Type, $params;
    map {  $x(:,$_->[0])  .= $_->[1] * exp(-$t*$t * $_->[2]  ) }  @$params;
# put any number of initial parameter pairs here
    my $p = pdl $Type, [ [ 5, 1], [2,4], [2,3], [1,1], [1.5, 3] ] ; # starting guess
    my $np = $p->dim(1);
    note "Trying x" . dimst $x->dummy(-1,$np);
    note "input  p" . dimst $p->dummy(1,$nx);
    check_type($Type, $p,$x,$t);
    my $h = levmar($p->dummy(1,$nx), $x->dummy(-1,$np), $t, $Gh , @g );
    note "check that output p has correct shape and values";

# Disabled for levmar-2.6 , not working (broken in 2.5)
#    is_pdl $h->{P}, $res(1:,:)->dummy(-1,$np),
#	"Thread both x and p, 2 thread dims";

    note "returned  p" . dimst $h->{P};
    note "and  covar" . dimst $h->{COVAR};

    my $covar = PDL->null;
    my $save_covar = $covar;
    check_type($Type, $p,$x,$t);
    $h = levmar($p->dummy(1,$nx), $x->dummy(-1,$np), $t, $Gh , @g,
		   COVAR => $covar);
    check_type($Type, $h->{INFO});
    my $count = $h->{COVAR}->nelem;
    $h = levmar($p->dummy(1,$nx), $x->dummy(-1,$np), $t, $Gh , @g,
		   COVAR => $covar);
    check_type($Type, $h->{INFO});
    $save_covar .= 1;
    my $sum = $h->{COVAR}->sum;
    ok( $sum == $count, "Test passing null COVAR pdl");
}

sub thread5 {
    my ($Type) = @_;
    my $n = 10000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
    my $x = zeroes($Type, $n,4);
    my $params =  [ [3,.2], [ 28, .1] , [2,.01], [3,.3] ];
    my $i = 0;
    map {  $x(:,$i++)  .= $_->[0] * exp(-$t*$t * $_->[1]  ) }  @$params;
    my $p = pdl $Type, [ 5, 1]; # starting guess
    check_type($Type, $p,$x,$t);
    my $h = levmar(  $p, $x, $t, $Gh, FIX=> [1,0], @g);
    check_type($Type, $h->{INFO});
    my $outp = pdl $Type, [[ 5, 0.4730849], [5, 1 ],
		   [5, 0.16286478],  [5, 0.70962698], ];

   is_pdl $h->{P}, $outp,
	{atol=>1e-3, test_name=>"Thread x, 1 thread dim, FIX=>[1,0] (linear constr.)"};

#   is_pdl $h->{P}, $outp,
#	{atol=>1e-3, test_name=>"Thread x, 1 thread dim, FIX=>[1,0] (linear constr.)"};

}

# same but easier to read
sub thread6 {
    my ($Type) = @_;
    my $n = 1000;
    my $t = zeroes($n)->xlinvals(-5,4.999)->convert($Type);
# Put any number of pairs of actual parameters here.
    my $params =  [ [500,.01], [3, .1] , [2,.01], [50,.3] ];
    my $nx = scalar(@$params);
    my $x = zeroes($n,$nx);
    my $i = 0;
    foreach( @$params ) {
	$x(:,$i++) .= $_->[0] * exp(-$t*$t * $_->[1]  );
    }
# put any number of initial parameter pairs here
   my $p = pdl [ [ 5, 1], [2,1], [2,3], [40,1], [1.5, 3] ] ; # starting guess
    my $np = $p->dim(1);
    note "Trying x" . dimst $x->dummy(-1,$np);
    note "input  p" . dimst $p->dummy(1,$nx);
    my $pd = $p->dummy(1,$nx);
    my $xd = $x->dummy(-1,$np);
    my $h = levmar($p->dummy(1,$nx), $x->dummy(-1,$np), $t, $Gh , @g );
    note "check that output p has correct shape and values";
    is_pdl $h->{P}, pdl($params)->dummy(-1,$np),
	"Thread both x and p, 2 thread dims";
    note "returned  p" . dimst $h->{P};
    note "and  covar" . dimst $h->{COVAR};
    note "and  info " . dimst $h->{INFO}->slice('(0),:,:');
    diag $h->{INFO}->slice('(0),:,:');
    my $inf = $h->{INFO};
    note "and  info " . dimst $h->{INFO}->slice('(0)');
    note "finally ".  dimst $inf->((0));
    note "finally ".  dimst $h->{REASON};
    my $r = $h->{REASON};
    diag $r;
    my $inds = which($r != 6);
    diag  pdl( [ $inds % $nx, $inds / $nx])->transpose;
    diag  pdl( [ $inds % $nx, $inds / $nx])->transpose;
#    diag $h->{RET};
}

#thread6(double);
thread1(double);
thread2(double);
thread3(double);
thread4(double);

if ($PDL::Fit::Levmar::HAVE_LAPACK) {
 thread5(double);
}

#thread1(float);
#thread2(float);
#thread3(float);
#thread4(float);
#thread5(float);

done_testing;
