use strict;
use warnings;
use PDL;
use PDL::Fit::Levmar;
use PDL::Fit::Levmar::Func;
use PDL::NiceSlice;
use Test::More;

my @t = ( TESTSYNTAX => 1);

sub test1 {
  my    $str1 = '
function gaussian1
x = p0 * exp(-(t-p1)*(t-p1)*p2);
';
  my $f = levmar_func(FUNC=>$str1, @t);
  like $f->{SYNTAXRESULTS},qr/:LOOP:/;
} # end test1

sub test2 {
    my  $MODROSLAM  = 1e2;

    my $str1 = "

    function modros
    x0 = 10 * (p1 -p0*p0);
    x1 = 1.0 - p0;
    x2 = $MODROSLAM;
    loop

    jacobian jacmodros
    d0[0] = -20 * p0;
    d1[0] = 10;
    d0[1] = -1;
    d1[1] = 0;
    d0[2] = 0;
    d1[2] = 0;
    loop
    
";

    my $str2 = "
    function modros
    noloop
    x0 = 10 * (p1 -p0*p0);
    x1 = 1.0 - p0;
    x2 = $MODROSLAM;

    jacobian jacmodros
    noloop
    d0[0] = -20 * p0;
    d1[0] = 10;
    d0[1] = -1;
    d1[1] = 0;
    d0[2] = 0;
    d1[2] = 0;
    
";
    my $f = levmar_func(FUNC=>$str1,@t);
    like $f->{SYNTAXRESULTS},qr/:LOOP:/, "explicit empty loop via loop directive";
    # I guess I'll let that go for now, because it's harmless, even
    # though there will be an empty loop

    $f = levmar_func(FUNC=>$str2, @t);
    unlike $f->{SYNTAXRESULTS},qr/:LOOP:/, "no loop via noloop directive";
}

sub test3 {
    my $str1 = '
#define ROSD 105.0
 function mros
     x =((1.0-p0)*(1.0-p0) + ROSD*(p1-p0*p0)*(p1-p0*p0));

 jacobian jacmros
    d1=(-2 + 2*p0-4*ROSD*(p1-p0*p0)*p0);
    d2=(2*ROSD*(p1-p0*p0));

';
    my $f = levmar_func(FUNC=>$str1, @t);
    like $f->{SYNTAXRESULTS},qr/:PREFUNC:/, "PREFUNC present";
} 

test1();
test2();
test3();

done_testing;
