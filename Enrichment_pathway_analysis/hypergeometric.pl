use strict;
use warnings;
use File::Slurp;
use Math::GSL::Randist;

my $population_size = (562);
my $population_successes = (38);
my $cluster_size = (4);
my $cluster_successes = (2);

my $pvalue = hypergeom_pvalue($cluster_successes, $cluster_size, $population_successes, $population_size);

print "P-valor de hipergeometrica para la ruta con genes diferencialmente expresados: $pvalue\n";

sub hypergeom_pvalue {
    my ($sample_successes, $sample_size, $population_successes, $population_size) = @_;
    my $accum_pvalue = 0;
    my $population_non_successes = $population_size - $population_successes;

    for my $x ($sample_successes .. $population_successes) {
        $accum_pvalue += Math::GSL::Randist::gsl_ran_hypergeometric_pdf($x, $population_successes, $population_non_successes, $sample_size);
    }

    return $accum_pvalue;
}
