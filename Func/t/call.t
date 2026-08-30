use PDL;
use PDL::Fit::Levmar::Func;
use PDL::LiteF;
use Test::More;
use Test::PDL -atol => 1e-4;
use strict;
use warnings;

my $Gf = levmar_func( FUNC => '
   function
   x = p0 * t * t;
 ');
is_pdl $Gf->call([2],sequence(10)),
  pdl(0, 2, 8, 18, 32, 50, 72, 98, 128, 162), " call func from lpp";

$Gf = levmar_func( FUNC => '
#include<string.h>
   function
   memset( x, 0, n );
   loop
   x = p0 * t * t;
 ');
is_pdl $Gf->call([2],sequence(10)),
  pdl(0, 2, 8, 18, 32, 50, 72, 98, 128, 162), "lpp func with #include";

done_testing;
