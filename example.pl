use strict;
use warnings;
use PDL;
use PDL::Fit::Levmar;
use PDL::Fit::Levmar::Func;
use PDL::NiceSlice;

sub example1 {
    my $n = 10;
    my $t = 10.0*(sequence($n)/$n -1/2);
    my $x = 3 * exp(-$t*$t * .3  );
    my $p = pdl [ 1, 1 ]; # initial guesses
    print levmar($p,$x,$t, FUNC =>
          '   function gaussian
              x = p0 * exp( -t*t * p1);
           ')->{REPORT};
}

sub myexp {
    my ($p,$x,$t) = @_;
    my $p0 = $p->at(0);
    my $p1 = $p->at(1);
    $x .= $p0 * exp(-$t*$t * $p1);
}

sub example2 {
    my $n = 100000;
    my $t = 10.0*(sequence($n)/$n -1/2);
    my $x = 3 * exp(-$t*$t * .3  );
    my $p = pdl [ 1, 1 ]; # initial guesses
    my $h = levmar($p,$x,$t, FUNC => \&myexp);
    print $h->{REPORT};
}

example1();
example2();
