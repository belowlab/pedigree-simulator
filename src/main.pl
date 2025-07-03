#! /usr/bin/perl

use strict;
use Carp;
use List::Util 'shuffle';

# Set up relative directories:
use Cwd ();
use File::Basename ();
use File::Spec ();

# This gets relative directory to this script, and loads our PRIMUS
use lib File::Spec->catdir(File::Basename::dirname(Cwd::abs_path __FILE__), '../dependencies/PRIMUS_v1.9.0/lib/perl_modules/');

use PRIMUS::node_v7;
use threads;

my $project_dir = File::Basename::dirname(Cwd::abs_path __FILE__); # top level for this script
my $output_root = "$project_dir/../output";	# Where to write output
my $output_directory = "";
my $simulation_dir = "$project_dir/morrison"; # python simulation scripts
my $cranefoot_binary = "$project_dir/../dependencies/cranefoot/cranefoot"; # check in the cranefoot binary

my $plink_binary = "plink";
my $plink2_binary = "plink2";

my $data_dir = "$project_dir/../data";
my $bgl2ped = "$project_dir/../dependencies/bgl_to_ped/bgl_to_ped";
my $chr_lengths = "$project_dir/chr_lengths.txt";

# Global flags
my $ECHO_SYSTEM = 0;

## Define variables
my $mean_children = 3;
my $sample_ctr = 0;
my $MAX_SAMPLES = 20;
my $ONEKG_pop = "EUR";
my $MAX_GENERATIONS = 15;
my $NUM_SIMULATIONS = 100;
my $num_complex_relationships = 1;
my $HALF_SIB_RATE = .3;
my $NUM_HS = 2;
my $parallel_status = "false";

my @PC_k0s;
my @FS_k0s;
my @HAG_k0s;
my @CGH_k0s;
my @DR_k0s;
my @UN_k0s;

my @PC_k1s;
my @FS_k1s;
my @HAG_k1s;
my @CGH_k1s;
my @DR_k1s;
my @UN_k1s;

my $Missing_ctr = 1;

my $nb_process = 30;
my $nb_compute;
my $i=0;
my @running = ();
my @Threads;

my %children; ## $children{sample} = [sample's children]
my %parents; ## $parents{sample} = [dad,mom]
my %samples; ## $samples{sample} = gender
my %generations; ## $generations{sample} = generation
my %married_into_family;
my @samples_to_visit;
my %pairwise_rels;

my $sim = shift; # Number (identifier of the simulation)
my $type = shift; # e.g., uniform3 

if($type eq "uniform3")
{
	$mean_children = 3;
	$sample_ctr = 0;
	$MAX_SAMPLES = shift;
	my $ONEKG_pop = shift;
	$parallel_status = shift // $parallel_status; # Set to "parallel" if provided correctly

	print("\nPopulation: $ONEKG_pop\n\n");
	$MAX_GENERATIONS = 20;
	$NUM_SIMULATIONS = 100;
	
	$output_directory = "$output_root/$ONEKG_pop";
	if(!-d $output_directory){run_system("mkdir $output_directory")}

	make_simulated_pedigree("$output_directory",$ONEKG_pop,"$type\_size$MAX_SAMPLES\_sim$sim");

	# missing ID steps
	# print ("\nmoving on to missing versions\n");
	# print "output_directory: $output_directory\n";
	# print "type: $type\n";
	# print "MAX_SAMPLES: $MAX_SAMPLES\n";
	# print "sim: $sim\n";
	# print "Missing_ctr: $Missing_ctr\n\n";
	
	make_incrementally_increasing_missing_versions($output_directory,$type,$MAX_SAMPLES,$sim,$Missing_ctr);
}

if($type eq "uniform2")
{
	$mean_children = 2;
	$sample_ctr = 0;
	$MAX_SAMPLES = shift;
	my $ONEKG_pop = shift;
	print("\npop: $ONEKG_pop\n\n");
	$MAX_GENERATIONS = 20;
	$NUM_SIMULATIONS = 100;
	
	$output_directory = "$output_root/$ONEKG_pop";
	if(!-d $output_directory){run_system("mkdir $output_directory")}

	make_simulated_pedigree("$output_directory",$ONEKG_pop,"$type\_size$MAX_SAMPLES\_sim$sim");

	# missing ID steps
	print ("\nmoving on to missing versions\n");
	print "output_directory: $output_directory\n";
	print "type: $type\n";
	print "MAX_SAMPLES: $MAX_SAMPLES\n";
	print "sim: $sim\n";
	print "Missing_ctr: $Missing_ctr\n\n";
	
	make_incrementally_increasing_missing_versions($output_directory,$type,$MAX_SAMPLES,$sim,$Missing_ctr);
}

# Added 1/30/25

elsif($type eq "halfsib3")
{
	$mean_children = 3;
	$sample_ctr = 0;
	$MAX_SAMPLES = shift;
	my $ONEKG_pop = shift;
	print("\npop: $ONEKG_pop\n\n");
	$MAX_GENERATIONS = 20;
	$NUM_SIMULATIONS = 100;
	
	if($MAX_SAMPLES > 3)
	{
		# $output_directory = "../data/simulations/$type/";
		# if(!-d $output_directory){system("mkdir $output_directory")}
		$output_directory = "$output_root/$ONEKG_pop";
		#if(!-d $output_directory){run_system("mkdir $output_directory")}

		#make_simulated_pedigree("$output_directory","CEU","$type\_size$MAX_SAMPLES\_sim$sim","halfsib");
		make_simulated_pedigree("$output_directory",$ONEKG_pop,"$type\_size$MAX_SAMPLES\_sim$sim", "halfsib"); # difference here is the inclusion of "halfsib"

		make_incrementally_increasing_missing_versions($output_directory,$type,$MAX_SAMPLES,$sim,$Missing_ctr);
	}
}

else
{
	print "\nUsage:\n $0 [sim_number] [sim_type] [sim_size] [sim_pop]\n\n";
	exit;
}

exit;



######################
## Subroutines
######################

sub check_if_reconstructable
{
	my $fam_file = shift;
	#print "fam_file: $fam_file\n";
	reset_variables();
	
	## Read in file
	open(IN,$fam_file) or die "can't open $fam_file\n";
	my $network_ref;
	
	## Build pedigree
	while(my $line = <IN>)
	{
		chomp($line);
		my ($FID,$IID,$PID,$MID,$SEX,$PHENOTYPE) = split(/\s+/,$line);
		my $child = "$IID";
		my $mom = "$MID";
		my $dad = "$PID";
		
		if(!exists $$network_ref{$child})
		{
			my $node = new PRIMUS::node_v7($child);
			$$network_ref{$child} = $node;
		}
		if(!exists $$network_ref{$dad} && $PID ne 0)
		{
			my $node = new PRIMUS::node_v7($dad);
			$$network_ref{$dad} = $node;
		}
		if(!exists $$network_ref{$mom} && $MID ne 0)
		{
			my $node = new PRIMUS::node_v7($mom);
			$$network_ref{$mom} = $node;
		}
		if($PID ne 0)
		{
			$$network_ref{$child}->add_parent($dad);
			$$network_ref{$dad}->add_child($child);
		}
		
		if($MID ne 0)
		{
			$$network_ref{$child}->add_parent($mom);
			$$network_ref{$mom}->add_child($child);
		}
	}
	close(IN);

	## Load network
	foreach my $node_name (keys %$network_ref)
	{
		make_relative_network_from_pedigree($$network_ref{$node_name},$network_ref);
		#my %rels = $$network_ref{$node_name}->relatives();
		my %rels = $$network_ref{$node_name}->all_possible_relationships();
		foreach(keys %rels)
		{
			#print "$node_name -> $_  = @{ $rels{$_} } \n";
		}
	}

	## Check that each sample is connected
	foreach my $node1 (keys %$network_ref)
	{
		#print "node1: $node1\n";
		if($node1 =~ /Missing/i){next}
		foreach my $node2 (keys %$network_ref)
		{
			if($node2 =~ /Missing/i){next}
			#print "node2: $node2\n";
			if($node1 eq $node2)
			{
				next;
			}
			my %visited;
			my $val = are_samples_connected($network_ref,$node1,$node2,\%visited);
			if($val == 0)
			{
				#print "NO $node1 <-> $node2\n";
				return 0;
			}
			else
			{
				#print "YES $node1 <-> $node2\n";
			}
		}
	}
	#print "YES\n";
	return 1;
}

sub are_samples_connected
{
	my $network_ref = shift;
	my $node1 = shift;
	my $node2 = shift;
	my $visited_ref = shift;
	
	#print "connected?   $node1 <-> $node2\n";
	
	$$visited_ref{$node1} = 1;

	my %node1_rels = $$network_ref{$node1}->all_possible_relationships();
	
	foreach(keys %node1_rels)
	{
		#print "$node1 -> $_  = @{ $node1_rels{$_} } \n";
	}
	
	if(@{ $node1_rels{$node2} }[0] =~ /(PC|FS|HAG|CGH)/)
	{
		#print "yes $node1 <-> $node2\n";
		return 1;
	}
	foreach my $new_node1 (keys %node1_rels)
	{
		if(exists $$visited_ref{$new_node1}){next}
		if(@{ $node1_rels{$new_node1} }[0] =~ /(DR|UN)/){next}
		my $val = are_samples_connected($network_ref,$new_node1,$node2,$visited_ref);
		if($val == 1)
		{
			return 1;
		}
	}
	#print "no $node1 <-> $node2\n";
	return 0;
}

sub make_incrementally_increasing_missing_versions {

	my $output_directory = shift;
	my $type = shift;
	my $MAX_SAMPLES =shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	my $sim_name = "$type\_size$MAX_SAMPLES\_sim$sim";
	my $sim_dir = "$output_directory/$sim_name";

	## testing to see what the inputs are 
	#print "running make_simple_gap\n";
	#print "sim_dir: $sim_dir\n";
	#print "sim_name: $sim_name\n";
	print "Missing counter: $Missing_ctr\n\n";
	
	#################################################################################################
	# OLD 

	$Missing_ctr = make_simple_gap($sim_dir,$sim_name,$Missing_ctr);

	for my $num_to_remove (1..($MAX_SAMPLES / 5))  # UPDATE 6/3/24 so that we can adjust for size 40 pedigrees 
	{
		my $ctr_before = $Missing_ctr;

		#print "for loop: running make_simple_gap again\n";
		#print "sim_dir: $sim_dir\n";
		#print "sim_dir: $sim_dir\n";
		#print "sim_name: $type\_size$MAX_SAMPLES\_sim$sim-$num_to_remove\n";
		print "Missing counter: $Missing_ctr\n\n";

		$Missing_ctr = make_simple_gap($sim_dir,"$type\_size$MAX_SAMPLES\_sim$sim-$num_to_remove",$Missing_ctr);

		# run plink to generate filesets 
		print "Running PLINK to generate plink-binary for $num_to_remove missing samples\n";
		run_system("$plink_binary --file $sim_dir/${sim_name}-$num_to_remove --make-bed --out $sim_dir/${sim_name}-$num_to_remove > /dev/null 2>&1");
		
		if($ctr_before eq $Missing_ctr)
		{
			print "Max removed of $num_to_remove\n";
			last;
		}
	}

	# Once done, exit the program
	print "\nSimple gap loop complete.\n"
	#exit;

	#################################################################################################
	# Jewett et al. approach (nodes and leaves, at random)

	# Initial removal of a small number of leaves and parents to create the first -X suffix .ped file
	# $Missing_ctr = remove_leaves_and_parents($sim_dir, $sim_name, $Missing_ctr, 1);

	# # Recalculate the number of removable individuals after the first removal
	# my $total_leaves_and_parents = count_leaves_and_parents($sim_dir, $sim_name);


	# for my $num_to_remove (1..int($total_leaves_and_parents * 0.9))  # Allow up to 90% of leaves + parents
	# {
	# 	my $ctr_before = $Missing_ctr;

	# 	print "for loop: running remove_leaves_and_parents again\n";
	# 	print "sim_dir: $sim_dir\n";
	# 	print "sim_name: $type\_size$MAX_SAMPLES\_sim$sim-$num_to_remove\n";
	# 	print "Missing_ctr: $Missing_ctr\n\n";

	# 	# New version that only removes leaves and parents
	# 	$Missing_ctr = remove_leaves_and_parents($sim_dir, "$type\_size$MAX_SAMPLES\_sim$sim-$num_to_remove", $Missing_ctr, $num_to_remove);

	# 	# Run PLINK to generate filesets
	# 	print "running plink to generate plink-binary for $num_to_remove\n\n";
	# 	run_system("$plink_binary --file $sim_dir/${sim_name}-$num_to_remove --make-bed --out $sim_dir/${sim_name}-$num_to_remove");
		
	# 	# Stop if no more individuals are removed
	# 	if ($ctr_before eq $Missing_ctr)
	# 	{
	# 		print "Max removed of $num_to_remove\n";
	# 		last;
	# 	}

	# 	# Recalculate the number of eligible individuals after removal
	# 	$total_leaves_and_parents = count_leaves_and_parents($sim_dir, $sim_name);
	# }


}

#########################
# New - testing a new subroutine that removes leaves and parents of leaves instead of nodes every time
#########################

sub count_leaves_and_parents {
    my ($simulation_dir, $sim) = @_;

    my $ped_file = "$simulation_dir/$sim.ped";
    open(my $PED, '<', $ped_file) or die "Can't open $ped_file: $!";
    my %parents;
    my %children;
    my %all_individuals;

    while (my $line = <$PED>) {
        chomp($line);
        my ($FID, $IID, $PID, $MID) = split(/\s+/, $line);
        next if $IID =~ /Missing/;
        $all_individuals{$IID} = 1;

        $parents{$IID} = [$PID, $MID];
        push @{$children{$PID}}, $IID if $PID ne '0';
        push @{$children{$MID}}, $IID if $MID ne '0';
    }
    close($PED);

    my @leaves = grep { !$children{$_} } keys %parents;
    my %parent_candidates;
    
    foreach my $leaf (@leaves) {
        my ($dad, $mom) = @{$parents{$leaf}};
        $parent_candidates{$dad} = 1 if $dad ne '0';
        $parent_candidates{$mom} = 1 if $mom ne '0';
    }

    return scalar(@leaves) + scalar(keys %parent_candidates);
}

sub remove_leaves_and_parents {
    my ($simulation_dir, $sim, $Missing_ctr, $num_to_remove) = @_;

    my $ped_file = "$simulation_dir/$sim.ped";
    #print "ped file: $simulation_dir/$sim.ped\n\n";
    my $fam_file = "$simulation_dir/$sim.fam";
    my $map_file = "$simulation_dir/$sim.map";
    my $sex_file = "$simulation_dir/$sim.sex";
    my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
    my $cranefoot_config_file = "$simulation_dir/$sim.config";
    my $genome_file = "$simulation_dir/$sim.genome";

    # Create new simulation name with incremented suffix
    my $new_sim_name = "$sim-1";
    if($sim =~ /(\w+)-(\d+)$/) {
        my $num_missing = $2;
        $num_missing++;
        $new_sim_name = "$1-$num_missing";
    }

    # First copy all files to new names
    #print "copying ped file: cp $ped_file $simulation_dir/$new_sim_name.ped\n\n";
    run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
    run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
    run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
    run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
    run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
    run_system("cp $map_file $simulation_dir/$new_sim_name.map");
    run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

    # Read pedigree information from the NEW files
    open(my $PED, '<', "$simulation_dir/$new_sim_name.ped") or die "Can't open $ped_file: $!";
    my %parents;
    my %children;
    my %all_individuals;

    # First pass: collect all relationships
    while (my $line = <$PED>) {
        chomp($line);
        my ($FID, $IID, $PID, $MID) = split(/\s+/, $line);
        next if $IID =~ /Missing/;
        $all_individuals{$IID} = 1;

        # Store parent-child relationships
        $parents{$IID} = [$PID, $MID];
        push @{$children{$PID}}, $IID if $PID ne '0';
        push @{$children{$MID}}, $IID if $MID ne '0';
    }
    close($PED);

    # Identify ALL initial leaves and their parents at once
    my %valid_candidates;
    
    # First identify all leaves
    my @initial_leaves = grep { !$children{$_} } keys %parents;
    foreach my $leaf (@initial_leaves) {
        $valid_candidates{$leaf} = 'leaf';
    }

    # Then identify all parents of initial leaves
    foreach my $leaf (@initial_leaves) {
        my ($dad, $mom) = @{$parents{$leaf}};
        if ($dad ne '0') {
            # Only add parent if ALL their children are leaves
            my $all_children_are_leaves = 1;
            foreach my $child (@{$children{$dad}}) {
                if (!exists $valid_candidates{$child}) {
                    $all_children_are_leaves = 0;
                    last;
                }
            }
            if ($all_children_are_leaves) {
                $valid_candidates{$dad} = 'parent';
            }
        }
        if ($mom ne '0') {
            my $all_children_are_leaves = 1;
            foreach my $child (@{$children{$mom}}) {
                if (!exists $valid_candidates{$child}) {
                    $all_children_are_leaves = 0;
                    last;
                }
            }
            if ($all_children_are_leaves) {
                $valid_candidates{$mom} = 'parent';
            }
        }
    }

    # Convert to array and shuffle
    my @removal_candidates = keys %valid_candidates;
    @removal_candidates = shuffle(@removal_candidates);

    # Limit to requested number
    if (scalar(@removal_candidates) > $num_to_remove) {
        @removal_candidates = @removal_candidates[0 .. $num_to_remove - 1];
    }

    # Remove selected individuals from the NEW files
    foreach my $ind (@removal_candidates) {
        remove_sample("$simulation_dir/$new_sim_name", $ind, "Missing$Missing_ctr", 0);
        $Missing_ctr++;
    }

    # Update config and generate new cranefoot plot
    rewrite_config("$simulation_dir/$new_sim_name.config", 
                  $simulation_dir, 
                  $simulation_dir,
                  $sim,
                  $new_sim_name);
    run_cranefoot($simulation_dir, $new_sim_name);

    return $Missing_ctr;
}

sub remove_leaves_and_parents_works {
    my ($simulation_dir, $sim, $Missing_ctr, $num_to_remove) = @_;

    my $ped_file = "$simulation_dir/$sim.ped";
    #print "ped file: $simulation_dir/$sim.ped\n\n";
    my $fam_file = "$simulation_dir/$sim.fam";
    my $map_file = "$simulation_dir/$sim.map";
    my $sex_file = "$simulation_dir/$sim.sex";
    my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
    my $cranefoot_config_file = "$simulation_dir/$sim.config";
    my $genome_file = "$simulation_dir/$sim.genome";

    # Create new simulation name with incremented suffix
    my $new_sim_name = "$sim-1";
    if($sim =~ /(\w+)-(\d+)$/) {
        my $num_missing = $2;
        $num_missing++;
        $new_sim_name = "$1-$num_missing";
    }

    # First copy all files to new names
    #print "copying ped file: cp $ped_file $simulation_dir/$new_sim_name.ped\n\n";
    run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
    run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
    run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
    run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
    run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
    run_system("cp $map_file $simulation_dir/$new_sim_name.map");
    run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

    # Read pedigree information from the NEW files
    open(my $PED, '<', "$simulation_dir/$new_sim_name.ped") or die "Can't open $ped_file: $!";
    my %parents;
    my %children;
    my %all_individuals;

    while (my $line = <$PED>) {
        chomp($line);
        my ($FID, $IID, $PID, $MID) = split(/\s+/, $line);
        next if $IID =~ /Missing/;
        $all_individuals{$IID} = 1;

        # Store parent-child relationships
        $parents{$IID} = [$PID, $MID];
        push @{$children{$PID}}, $IID if $PID ne '0';
        push @{$children{$MID}}, $IID if $MID ne '0';
    }
    close($PED);

    # Identify leaves (individuals without children)
    my @leaves = grep { !$children{$_} } keys %parents;

    # Collect removal candidates (leaves first, then parents)
    my @removal_candidates;
    my %removed;

    # First, add leaves to removal candidates
    foreach my $leaf (@leaves) {
        push @removal_candidates, $leaf;
        last if scalar(@removal_candidates) >= $num_to_remove;
    }

    # Then, add parents if more removal is needed
    foreach my $leaf (@leaves) {
        last if scalar(@removal_candidates) >= $num_to_remove;
        my ($dad, $mom) = @{$parents{$leaf}};
        foreach my $parent ($dad, $mom) {
            next if $parent eq '0' || exists $removed{$parent};
            my @siblings = grep { $_ ne $leaf } @{$children{$parent} || []};
            unless (@siblings) {
                push @removal_candidates, $parent;
                last if scalar(@removal_candidates) >= $num_to_remove;
            }
        }
    }

    # Shuffle and pick exactly the required number to remove
    @removal_candidates = shuffle(@removal_candidates);
    @removal_candidates = @removal_candidates[0 .. $num_to_remove - 1] 
        if scalar(@removal_candidates) > $num_to_remove;

    # Remove selected individuals from the NEW files
    # Key change: Always set is_leaf to 0 to keep individuals in structure
    foreach my $ind (@removal_candidates) {
        # Always pass is_leaf as 0 to keep them in structure but marked as Missing
        remove_sample("$simulation_dir/$new_sim_name", $ind, "Missing$Missing_ctr", 0);
        $removed{$ind} = 1;
        $Missing_ctr++;
    }

    # Update config and generate new cranefoot plot
    rewrite_config("$simulation_dir/$new_sim_name.config", 
                  $simulation_dir, 
                  $simulation_dir,
                  $sim,
                  $new_sim_name);
    run_cranefoot($simulation_dir, $new_sim_name);

    return $Missing_ctr;
}

sub remove_leaves_and_parents_old {
    my ($simulation_dir, $sim, $Missing_ctr, $num_to_remove) = @_;

    my $ped_file = "$simulation_dir/$sim.ped";
    print "ped file: $simulation_dir/$sim.ped\n\n";
    my $fam_file = "$simulation_dir/$sim.fam";
    my $map_file = "$simulation_dir/$sim.map";
    my $sex_file = "$simulation_dir/$sim.sex";
    my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
    my $cranefoot_config_file = "$simulation_dir/$sim.config";
    my $genome_file = "$simulation_dir/$sim.genome";

    # Create new simulation name with incremented suffix
    my $new_sim_name = "$sim-1";
    if($sim =~ /(\w+)-(\d+)$/) {
        my $num_missing = $2;
        $num_missing++;
        $new_sim_name = "$1-$num_missing";
    }

    # First copy all files to new names
    #print "copying ped file: cp $ped_file $simulation_dir/$new_sim_name.ped\n\n";
    run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
    run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
    run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
    run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
    run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
    run_system("cp $map_file $simulation_dir/$new_sim_name.map");
    run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

    # Read pedigree information from the NEW files
    open(my $PED, '<', "$simulation_dir/$new_sim_name.ped") or die "Can't open $ped_file: $!";
    my %parents;
    my %children;
    my %all_individuals;

    while (my $line = <$PED>) {
        chomp($line);
        my ($FID, $IID, $PID, $MID) = split(/\s+/, $line);
        next if $IID =~ /Missing/;
        $all_individuals{$IID} = 1;

        # Store parent-child relationships
        $parents{$IID} = [$PID, $MID];
        push @{$children{$PID}}, $IID if $PID ne '0';
        push @{$children{$MID}}, $IID if $MID ne '0';
    }
    close($PED);

    # Identify leaves (individuals without children)
    my @leaves = grep { !$children{$_} } keys %parents;

    # Collect removal candidates (leaves first, then parents)
    my @removal_candidates;
    my %removed;

    # First, add leaves to removal candidates
    foreach my $leaf (@leaves) {
        push @removal_candidates, $leaf;
        last if scalar(@removal_candidates) >= $num_to_remove;
    }

    # Then, add parents if more removal is needed
    foreach my $leaf (@leaves) {
        last if scalar(@removal_candidates) >= $num_to_remove;
        my ($dad, $mom) = @{$parents{$leaf}};
        foreach my $parent ($dad, $mom) {
            next if $parent eq '0' || exists $removed{$parent};
            my @siblings = grep { $_ ne $leaf } @{$children{$parent} || []};
            unless (@siblings) {
                push @removal_candidates, $parent;
                last if scalar(@removal_candidates) >= $num_to_remove;
            }
        }
    }

    # Shuffle and pick exactly the required number to remove
    @removal_candidates = shuffle(@removal_candidates);
    @removal_candidates = @removal_candidates[0 .. $num_to_remove - 1] 
        if scalar(@removal_candidates) > $num_to_remove;

    # Remove selected individuals from the NEW files
    foreach my $ind (@removal_candidates) {
        my $is_leaf = exists $children{$ind} ? 0 : 1;
        remove_sample("$simulation_dir/$new_sim_name", $ind, "Missing$Missing_ctr", $is_leaf);
        $removed{$ind} = 1;
        $Missing_ctr++;
    }

    # Update config and generate new cranefoot plot
    rewrite_config("$simulation_dir/$new_sim_name.config", 
                  $simulation_dir, 
                  $simulation_dir,
                  $sim,
                  $new_sim_name);
    run_cranefoot($simulation_dir, $new_sim_name);

    return $Missing_ctr;
}



sub add_genotypes {

	my $fam_file_root = shift; ## Must have 6 columns
	my $ONEKG_pop = shift;
	my $fam_file = "$fam_file_root.fam";
	if(!-e $fam_file){die "$fam_file does not exists. Must be a .fam file\n";}

	#print "Generating ped and map files for $fam_file\n";
	
	## Diagram pedigree
	run_system("python $simulation_dir/recomb_sim.py --fam=$fam_file --chr-lengths=$chr_lengths --out=$fam_file_root\_diag.txt");

	## Get IBD from diagrams
	run_system("python $simulation_dir/diagram_ibd.py --out=$fam_file_root\_diag.IBD --chrom_num=23 $fam_file_root\_diag.txt");
	
	# string to concatenate with all chromosomal commands to execute with parallel
	my $command_str = '';

	for my $chr(1..22)
	{

		my $genetic_map_file = "$data_dir/map.gz";

		# Set reference file for dropping in genotypes 

		## UPDATE 11/22/24
		my $reference_sample_file = "$data_dir/reference/$ONEKG_pop/1KG.$ONEKG_pop.GRCH38.rsID.chr$chr.vcf.gz";   
		

		## build command string for each chromosome
		my $full_command = "python $simulation_dir/sim_to_genotypes.py --mode=IBD --chr=$chr --simulation=$fam_file_root\_diag.txt --map $genetic_map_file --out=$fam_file_root.chr$chr\_diag.vcf --chunk=10000 --write-names=$fam_file_root.chr$chr\_diag.names --gzip $reference_sample_file\n";

		## concatenate the 22-chr string to be read into parallel
		$command_str = $command_str . $full_command;
	}

	# Run genotype dropping commands in parallel if notated at runtime

	if ($parallel_status eq "parallel") {
		print "\nRunning add_genotypes commands in parallel\n";
		run_system("echo \"$command_str\" | parallel -j22");
	}
	else {
		print "\nRunning add_genotypes commands with single thread\n";
		run_system("$command_str");
	}

	my $bcftools_command = 'bcftools concat ';
	for my $chr(1..22)
	{
		$bcftools_command = $bcftools_command . "$fam_file_root.chr$chr\_diag.vcf ";
	}
	$bcftools_command = $bcftools_command . "-o $fam_file_root\_all_chr.vcf.gz -O z > /dev/null 2>&1";  
	#print "\n\n$bcftools_command\n\n";

	print "\nConcatenating output with bcftools ...\n";

	run_system($bcftools_command);

	run_system("rm $fam_file_root.chr*");

}


sub rewrite_config
{
	my $file = shift;
	my $old_path = shift;
	my $new_path = shift;
	my $old_name = shift;
	my $new_name = shift;
	
	#print "old: $old_path\n";
	#print "new: $new_path\n";
	#print "old_name: $old_name\n";
	#print "new_name: $new_name\n";

	## Rewrite .config file
	open(C,$file) or die "can't open $file\n";
	my @lines = <C>;
	close(C);

	open(NEW,">$file") or die "can't open $file\n";
	foreach my $line (@lines)
	{
		if($line =~ /\s+$old_name/)
		{
			$line =~ s/$old_name/$new_name/;
		}
		$line =~ s/$old_path/$new_path/;
		$line =~ s/$old_name.cranefoot.ped/$new_name.cranefoot.ped/;
		print NEW $line;
	}
	close(NEW);
}

sub remove_sample
{
	my $file_root = shift;
	my $id = shift;
	my $Missing_name = shift;
	my $is_leaf = shift;
	#print "id:$id; Missing_name: $Missing_name\n";
	
	my $ped_file = "$file_root.ped"; 
	my $fam_file = "$file_root.fam";
	my $cranefoot_ped_file = "$file_root.cranefoot.ped";
	my $cranefoot_config_file = "$file_root.config";
	my $genome_file = "$file_root.genome";

	## rewrite fam file
	open(FAM,$fam_file) or die "can't open $fam_file\n";
	my @lines = <FAM>;
	close(FAM);
	open(NEW_FAM,">$fam_file") or die "can't open $fam_file\n";
	foreach my $line (@lines)
	{
		if($line =~ /^Fam\d+\s+$id\s+/)
		{
			if($Missing_name == -1){next}
		}
		if(($line =~ /^Fam\d+\s+$id\s+/) && $is_leaf)
		{
			next;
		}
		$line =~ s/\s+$id\s+/\t$Missing_name\t/g; # removed ` here
		$line =~ s/\s+$id\n/\t$Missing_name\n/g;
		print NEW_FAM "$line";
	}
	close(NEW_FAM);

	## rewrite ped file
	my $temp = `wc -l $file_root.map`;
	my ($num_SNPs,$map_file_path) = split(/\s+/,$temp);
	my $missing_SNPs;
	for(1..$num_SNPs*2) # the number of SNPs = 2 times the .map file length
	{
		$missing_SNPs .= "0 ";
	}
	chop($missing_SNPs); #remove that extra space at the end of the line

	run_system("mv $ped_file $ped_file.temp");
	open(PED,"$ped_file.temp") or die "can't open line 615 $ped_file.temp\n";
	open(NEW_PED,">$ped_file") or die "can't open line 616 $ped_file\n";
	while( my $line = <PED>)
	{
		if($line =~ /^\w+\d+\s$id\s/)
		{
			if($Missing_name == -1 || $is_leaf){next}
			my ($FID,$IID,$PID,$MID,$SEX,$A_STATUS) = split(/\s+/,$line);
			$line = "$FID $IID $PID $MID $SEX $A_STATUS $missing_SNPs\n";
			#print "$line";
		}
		$line =~ s/^(\w+)\s+$id\s+/$1 $Missing_name /g;	# IID
		$line =~ s/^(\w+\s+\w+)\s+$id\s+/$1 $Missing_name /g; # PID
		$line =~ s/^(\w+\s+\w+\s+\w+)\s+$id\s+/$1 $Missing_name /g; # MID
		#$line =~ s/\s$id\n/ $Missing_name\n/g;
		print NEW_PED "$line";
	}
	close(NEW_PED);
	close(PED);
	run_system("rm $ped_file.temp");
	
	## Rewrite cranefoot file
	open(C,$cranefoot_ped_file) or die "can't open $cranefoot_ped_file\n";
	my @lines = <C>;
	close(C);
	
	open(NEW_C,">$cranefoot_ped_file") or die "can't open $cranefoot_ped_file\n";
	foreach my $line (@lines)
	{
		if($line =~ /\t$id\t/)
		{
			if($line =~ /^Fam\d+\t$id\t/)
			{
				$line =~ s/1\n/11\n/;
				if($Missing_name == -1 || $is_leaf){
				next;}
			}
			$line =~ s/\t$id\t/\t$Missing_name\t/g;
		}
		print NEW_C "$line";
	}
	close(NEW_C);
	
	## Rewrite .genome file
	open(G,$genome_file) or die "can't open $genome_file\n";
	my @lines = <G>;
	close(G);
	open(NEW_G,">$genome_file") or die "can't open $genome_file\n";
	foreach my $line (@lines)
	{
		if($line !~ /\s$id\s/)
		{
			print NEW_G $line;
		}
	}
	close(NEW_G);
	
	return 1;
}

sub remove_samples
{
	my $complete_dir = shift;
	my $simulation_dir = shift;
	my $sim = shift;
	
	my $ped_file = "$complete_dir/$sim.ped";
	my $fam_file = "$complete_dir/$sim.fam";
	my $sex_file = "$complete_dir/$sim.sex";
	my $cranefoot_ped_file = "$complete_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$complete_dir/$sim.config";
	my $genome_file = "$complete_dir/$sim.genome";
	

	## Make a simple 1 gap removed
	$Missing_ctr = 1;
	my $new_sim_name = "$sim.1S";
	my $output_dir = "$simulation_dir/one_small_gap";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	make_simple_gap($output_dir,$new_sim_name,$Missing_ctr);
	run_cranefoot("$output_dir","$new_sim_name");


	## Put multiple 1 person gaps in the pedigree
	$Missing_ctr = 1;

	my $new_sim_name = "$sim.MS";
	my $output_dir = "$simulation_dir/multiple_small_gap";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	my ($sample_size, $name) = `wc -l $ped_file $output_dir/$new_sim_name.ped`;
	for(my $i = 0; $i < $sample_size/15; $i++)
	{
		make_simple_gap($output_dir,$new_sim_name,$Missing_ctr);
		$Missing_ctr++;
	}
	
	run_cranefoot("$output_dir","$new_sim_name");
	
	
	## Make simulations with single 2 person gap
	$Missing_ctr = 1;
	my $new_sim_name = "$sim.1D";
	my $output_dir = "$simulation_dir/one_double_gap";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	my $pass = make_double_gap($output_dir,$new_sim_name,$Missing_ctr);
	if($pass == 0)
	{
		run_system("rm $output_dir/$new_sim_name.ped");
		run_system("rm $output_dir/$new_sim_name.sex");
		run_system("rm $output_dir/$new_sim_name.cranefoot.ped");
		run_system("rm $output_dir/$new_sim_name.config");
		run_system("rm $output_dir/$new_sim_name.genome");
	}
	else
	{
		run_cranefoot("$output_dir","$new_sim_name");
	}

	## REMOVE Parents and children
	$Missing_ctr = 1;
	my $new_sim_name = "$sim.1DG";
	my $output_dir = "$simulation_dir/one_double_generation_gap";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	my $pass = remove_parents_and_children($output_dir,$new_sim_name,$Missing_ctr);
	if($pass == 0)
	{
		run_system("rm $output_dir/$new_sim_name.ped");
		run_system("rm $output_dir/$new_sim_name.sex");
		run_system("rm $output_dir/$new_sim_name.cranefoot.ped");
		run_system("rm $output_dir/$new_sim_name.config");
		run_system("rm $output_dir/$new_sim_name.genome");
	}
	else
	{
		run_cranefoot("$output_dir","$new_sim_name");
	}
	

	## Remove parent pair
	$Missing_ctr = 1;
	my $new_sim_name = "$sim.1PP";
	my $output_dir = "$simulation_dir/one_parent_pair_gap";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	my $pass = remove_parent_pair($output_dir,$new_sim_name,$Missing_ctr);
	if($pass == 0)
	{
		run_system("rm $output_dir/$new_sim_name.ped");
		run_system("rm $output_dir/$new_sim_name.sex");
		run_system("rm $output_dir/$new_sim_name.cranefoot.ped");
		run_system("rm $output_dir/$new_sim_name.config");
		run_system("rm $output_dir/$new_sim_name.genome");
	}
	else
	{
		run_cranefoot("$output_dir","$new_sim_name");
	}


	## REMOVE A COMBINATINO
	$Missing_ctr = 1;
	my $new_sim_name = "$sim.COMBO";
	my $output_dir = "$simulation_dir/COMBO";
	
	if(!-d $output_dir)
	{
		run_system("mkdir $output_dir");
	}
	
	run_system("cp $ped_file $output_dir/$new_sim_name.ped");
	run_system("cp $sex_file $output_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $output_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $output_dir/$new_sim_name.config");
	run_system("cp $genome_file $output_dir/$new_sim_name.genome");
	
	rewrite_config("$output_dir/$new_sim_name.config",$complete_dir,$output_dir,$sim,$new_sim_name);
	$Missing_ctr = make_simple_gap($output_dir,$new_sim_name,$Missing_ctr);
	$Missing_ctr = make_double_gap($output_dir,$new_sim_name,$Missing_ctr);
	$Missing_ctr = remove_parents_and_children($output_dir,$new_sim_name,$Missing_ctr);
	$Missing_ctr = remove_parent_pair($output_dir,$new_sim_name,$Missing_ctr);
	

	if($pass == 0)
	{
		run_system("rm $output_dir/$new_sim_name.ped");
		run_system("rm $output_dir/$new_sim_name.sex");
		run_system("rm $output_dir/$new_sim_name.cranefoot.ped");
		run_system("rm $output_dir/$new_sim_name.config");
		run_system("rm $output_dir/$new_sim_name.genome");
	}
	else
	{
		run_cranefoot("$output_dir","$new_sim_name");
	}
}

sub remove_parent_pair
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	
	## Get ped info
	open(PED,$ped_file) or die "can't open line 879 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;

	foreach my $line (@lines)
	{
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
	}

	## Select sample to remove
	my $continue = 1;
	my ($FID,$IID,$PID,$MID);
	my $parent;
	while($continue != 0 && $continue < 1000)
	{
		$continue++;
		my $size = @lines;
		my $val = int(rand($size));
		my $line = @lines[$val];
		($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		#print "\n\n\n\nIIDs to remove: $IID ($PID $MID) and $parent\n";
		if($IID =~ /Missing/){next;}
		if($PID =~ /Missing/){next;}
		if($MID =~ /Missing/){next;}
		if(exists $has_a_parent{$IID}){$continue = 0}
	}
	if($continue >= 1000){return $Missing_ctr}
	print "\n\nIIDs parents to remove: $PID $MID\n";
	
	## Remove samples
	remove_sample("$simulation_dir/$sim",$PID,"Missing$Missing_ctr");
	$Missing_ctr++;
	remove_sample("$simulation_dir/$sim",$MID,"Missing$Missing_ctr");
	$Missing_ctr++;

	return $Missing_ctr;
}

## This makes a double gap, but more eliminates the option of connected with hag (AV) relationships
sub remove_parents_and_children
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	
	## Get ped info
	open(PED,$ped_file) or die "can't open line 938 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;
	my %parents;

	foreach my $line (@lines)
	{
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
			push(@{ $parents{$PID}{$MID} }, $IID);
		}
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
	}

	## Select sample to remove
	my $continue = 1;
	my ($FID,$IID,$PID,$MID);
	my $parent;
	while($continue != 0 && $continue < 1000)
	{
		$continue++;
		my $size = @lines;
		my $val = int(rand($size));
		my $line = @lines[$val];
		($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		#print "\n\n\n\nIIDs to remove: $IID ($PID $MID) and $parent\n";
		if($IID =~ /Missing/){next;}
		if($PID =~ /Missing/){next;}
		if($MID =~ /Missing/){next;}
		
		if((exists $has_a_parent{$PID} || exists $has_a_parent{$PID}) && exists $is_a_parent{$IID}){$continue = 0;}
	}
	if($continue >= 1000){return $Missing_ctr}
	print "\n\nIIDs to remove: $IID ($PID $MID)\n";
	
	## Remove samples
	my @children = @{ $parents{$PID}{$MID} };
	print "$PID $MID @children\n";
	
	remove_sample("$simulation_dir/$sim",$PID,"Missing$Missing_ctr");
	$Missing_ctr++;
	remove_sample("$simulation_dir/$sim",$MID,"Missing$Missing_ctr");
	$Missing_ctr++;

	foreach my $child (@children)
	{
		if(exists $is_a_parent{$child})
		{
			remove_sample("$simulation_dir/$sim",$child,"Missing$Missing_ctr");
			$Missing_ctr++;
		}
		else
		{
			remove_sample("$simulation_dir/$sim",$child,-1);
		}
	}
	return $Missing_ctr;
	
}

## Create a double person gap (remove one individual with parents and child)
sub make_double_gap
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$simulation_dir/$sim.config";
	my $genome_file = "$simulation_dir/$sim.genome";
	
	## Get ped info
	open(PED,$ped_file) or die "can't open line 1019 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;

	foreach my $line (@lines)
	{
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
	}

	## Select sample to remove
	my $continue = 1;
	my ($FID,$IID,$PID,$MID);
	my $parent;
	while($continue != 0 && $continue < 1000)
	{
		$continue++;
		my $size = @lines;
		#print "size: $size\n";
		my $val = int(rand($size));
		#print "$val\n";
		#exit;
		my $line = @lines[$val];
		($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		$parent = $PID;
		#print "\n\n\n\nIIDs to remove: $IID ($PID $MID) and $parent\n";
		if($IID =~ /Missing/){next;}
		if($PID =~ /Missing/){$parent=$MID}
		if($parent =~ /Missing/){next;}
		
		if(exists $has_a_parent{$parent} && exists $is_a_parent{$IID}){$continue = 0;}

	}
	if($continue >= 1000){return $Missing_ctr}
	print "\n\nIIDs to remove: $IID ($PID $MID) and $parent\n";
	
	## Remove samples
	remove_sample("$simulation_dir/$sim",$IID,"Missing$Missing_ctr");
	$Missing_ctr++;
	remove_sample("$simulation_dir/$sim",$parent,"Missing$Missing_ctr");
	$Missing_ctr++;
	
	return $Missing_ctr;
}

## Create a version with only the leaf samples remaining
sub remove_all_samples_except_leafs
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	my $fam_file = "$simulation_dir/$sim.fam";
	my $map_file = "$simulation_dir/$sim.map";
	my $sex_file = "$simulation_dir/$sim.sex";
	my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$simulation_dir/$sim.config";
	my $genome_file = "$simulation_dir/$sim.genome";
	
	my $new_sim_name = "$sim-leafs";

	run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
	run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
	run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
	run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
	run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
	run_system("cp $map_file $simulation_dir/$new_sim_name.map");
	run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

	
	## Get ped info
	open(PED,$ped_file) or die "can't open 1353 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;


	foreach my $line (@lines)
	{
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
	}

	## Select sample to remove
	my $reconstructable = 0;
	my %invalid_vals;
	
	my $continue = 1;
	my $continue_MAX = 2000;
	while($reconstructable == 0 && $continue < $continue_MAX)
	{
		$continue = 1;
		my ($FID,$IID,$PID,$MID);
		my $val;
		while($continue != 0 && $continue < $continue_MAX)
		{
			$continue++;
			my $size = @lines;
			#print "size: $size\n";
			$val = int(rand($size));
			#print "$val\n";
			#exit;
			my $line = @lines[$val];
			($FID,$IID,$PID,$MID) = split(/\s+/,$line);
			#print "IID to remove: $IID ($PID $MID)\n";
			if($IID =~ /Missing/){next;}
			#if($PID =~ /id/ && $MID =~ /id/ && exists $is_a_parent{$IID}){$continue = 0;}
			if(exists $is_a_parent{$IID} && !exists $invalid_vals{$val}){$continue = 0;}
		}
		
		$invalid_vals{$val} = 1; ## Add all $vals because if it is valid then it won't come up again.
		
		if($continue >= 2000){return $Missing_ctr}
		print "IID to remove: $IID ($PID $MID)\n";
		
		if($sim =~ /(\w+)-(\d+)$/)
		{
			my $num_missing = $2;
			$num_missing++;
			$new_sim_name = "$1-$num_missing";
		}

		## Rewrite ped file
		remove_sample("$simulation_dir/$new_sim_name",$IID,"Missing$Missing_ctr");
		$Missing_ctr++;
	}
	
	## Fix cranefoot plot
	rewrite_config("$simulation_dir/$new_sim_name.config",$simulation_dir,$simulation_dir,$sim,$new_sim_name);
	run_cranefoot($simulation_dir,$new_sim_name);
}

sub remove_unnecessary
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	print "Remove Unnecessary\n";

	my $ped_file = "$simulation_dir/$sim.ped";
	my $fam_file = "$simulation_dir/$sim.fam";
	my $map_file = "$simulation_dir/$sim.map";
	my $sex_file = "$simulation_dir/$sim.sex";
	my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$simulation_dir/$sim.config";
	my $genome_file = "$simulation_dir/$sim.genome";
	
	## Get ped info
	open(PED,$ped_file) or die "can't open 1440 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;
	my %IIDs;

	foreach my $line (@lines)
	{
		#print "line: $line\n";
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		$IIDs{$IID} = 1;
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
	}
	
	foreach my $IID (keys %IIDs)
	{
		print "IID: $IID\n";
		if(!exists $is_a_parent{$IID} && $IID =~ /Missing/)
		{
			print "removing: $IID\n";
			my $is_leaf = 1;
			remove_sample("$simulation_dir/$sim",$IID,"Missing$Missing_ctr",$is_leaf);
		}
	}
}

## Remove any samples
sub remove_any_sample
{
	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	my $fam_file = "$simulation_dir/$sim.fam";
	my $map_file = "$simulation_dir/$sim.map";
	my $sex_file = "$simulation_dir/$sim.sex";
	my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$simulation_dir/$sim.config";
	my $genome_file = "$simulation_dir/$sim.genome";
	
	## Get ped info
	#print "$fam_file\n";
	open(PED,$fam_file) or die "can't open $fam_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;


	foreach my $line (@lines)
	{
		#print "line: $line\n";
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
	}

	## Select sample to remove
	my $reconstructable = 0;
	my %invalid_vals;
	
	my $continue_MAX = 2000;
	my $continue = 1;
	while($reconstructable == 0 && $continue < $continue_MAX)
	{
		$continue = 1;
		my ($FID,$IID,$PID,$MID);
		my $val;
		while($continue != 0 && $continue < $continue_MAX)
		{
			$continue++;
			my $size = @lines;
			#print "size: $size\n";
			$val = int(rand($size));
			#print "$val\n";
			#exit;
			my $line = @lines[$val];
			($FID,$IID,$PID,$MID) = split(/\s+/,$line);
			#print "IID to remove: $IID ($PID $MID)\n";
			if($IID =~ /Missing/){next;}
			#if($PID =~ /id/ && $MID =~ /id/ && exists $is_a_parent{$IID}){$continue = 0;}
			if(!exists $invalid_vals{$val}){$continue = 0;}
		}
		

		$invalid_vals{$val} = 1; ## Add all $vals because if it is valid then it won't come up again.
		
		if($continue >= 2000){return $Missing_ctr}
		#print "IID to remove: $IID ($PID $MID)\n";
		#my $new_sim_name = "temp";
		my $sim_name = $sim;
		my $new_sim_name = "$sim";

		#run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
		#run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
		#run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
		#run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
		#run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
		#run_system("cp $map_file $simulation_dir/$new_sim_name.map");
		#run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

		## Rewrite ped file
		my $is_leaf = 0;
		if(!exists $is_a_parent{$IID}){$is_leaf = 1}
		remove_sample("$simulation_dir/$new_sim_name",$IID,"Missing$Missing_ctr",$is_leaf);

		#$reconstructable = check_if_reconstructable("$simulation_dir/$new_sim_name.fam");
		#print "reconstructable: $reconstructable\n";
		$reconstructable = 1;
		
		## Fix cranefoot plot
		rewrite_config("$simulation_dir/$new_sim_name.config",$simulation_dir,$simulation_dir,$sim,$new_sim_name);
		#run_cranefoot($simulation_dir,$new_sim_name);
		
		#if($reconstructable == 1)
		if(-1 == 1)
		{
			run_system("cp $simulation_dir/$new_sim_name.ped $simulation_dir/$sim_name.ped");
			run_system("cp $simulation_dir/$new_sim_name.sex $simulation_dir/$sim_name.sex");
			run_system("cp $simulation_dir/$new_sim_name.cranefoot.ped $simulation_dir/$sim_name.cranefoot.ped");
			run_system("cp $simulation_dir/$new_sim_name.config $simulation_dir/$sim_name.config");
			run_system("cp $simulation_dir/$new_sim_name.genome $simulation_dir/$sim_name.genome");
			run_system("cp $simulation_dir/$new_sim_name.map $simulation_dir/$sim_name.map");
			run_system("cp $simulation_dir/$new_sim_name.fam $simulation_dir/$sim_name.fam");
			#print "invalid $simulation_dir/$new_sim_name.fam\n";
			#exit;
		}
	}

	$Missing_ctr++;
	return $Missing_ctr;
}

## Create a single person gap (remove one individual with at least one child)
sub make_simple_gap {

	my $simulation_dir = shift;
	my $sim = shift;
	my $Missing_ctr = shift;
	
	my $ped_file = "$simulation_dir/$sim.ped";
	#print "ped file: $simulation_dir/$sim.ped\n\n";
	my $fam_file = "$simulation_dir/$sim.fam";
	my $map_file = "$simulation_dir/$sim.map";
	my $sex_file = "$simulation_dir/$sim.sex";
	my $cranefoot_ped_file = "$simulation_dir/$sim.cranefoot.ped";
	my $cranefoot_config_file = "$simulation_dir/$sim.config";
	my $genome_file = "$simulation_dir/$sim.genome";
	
	## Get ped info
	open(PED,$ped_file) or die "can't open line 1355 $ped_file\n";
	my @lines = <PED>;
	close(PED);
	
	## Get parent/child info
	my %has_a_parent;
	my %is_a_parent;

	foreach my $line (@lines)
	{
		my ($FID,$IID,$PID,$MID) = split(/\s+/,$line);
		$is_a_parent{$PID} = $IID;
		$is_a_parent{$MID} = $IID;
		if($IID =~ /Missing/){next;}
		if($PID =~ /id/ && $MID =~ /id/)
		{
			$has_a_parent{$IID} = 1;
		}
	}

	## Select sample to remove
	my $reconstructable = 0;
	my %invalid_vals;
	
	my $continue = 1;
	my $continue_MAX = 25;
	while($reconstructable == 0 && $continue < $continue_MAX)
	{
		$continue = 1;
		my ($FID,$IID,$PID,$MID);
		my $val;
		while($continue != 0 && $continue < $continue_MAX)
		{
			$continue++;
			my $size = @lines;
			$val = int(rand($size));
			my $line = @lines[$val];
			($FID,$IID,$PID,$MID) = split(/\s+/,$line);
			if($IID =~ /Missing/){next;}
			if(exists $is_a_parent{$IID} && !exists $invalid_vals{$val}){$continue = 0;}
		}
		
		$invalid_vals{$val} = 1; ## Add all $vals because if it is valid then it won't come up again.
		
		if($continue >= 2000){return $Missing_ctr} # script is dying here 

		print "IID to remove: $IID\n";
		
		my $new_sim_name = "$sim-1";
		if($sim =~ /(\w+)-(\d+)$/)
		{
			my $num_missing = $2;
			$num_missing++;
			$new_sim_name = "$1-$num_missing";
		}

		# we are somehow not getting to this step where we copy the ped files + everything else 
		#print "copying ped file: cp $ped_file $simulation_dir/$new_sim_name.ped\n\n";
		run_system("cp $ped_file $simulation_dir/$new_sim_name.ped");
		run_system("cp $sex_file $simulation_dir/$new_sim_name.sex");
		run_system("cp $cranefoot_ped_file $simulation_dir/$new_sim_name.cranefoot.ped");
		run_system("cp $cranefoot_config_file $simulation_dir/$new_sim_name.config");
		run_system("cp $genome_file $simulation_dir/$new_sim_name.genome");
		run_system("cp $map_file $simulation_dir/$new_sim_name.map");
		run_system("cp $fam_file $simulation_dir/$new_sim_name.fam");

		## Rewrite ped file

		#print "running remove_sample with arguments: $simulation_dir/$new_sim_name, $IID, and Missing$Missing_ctr\n\n";

		remove_sample("$simulation_dir/$new_sim_name",$IID,"Missing$Missing_ctr");

		$reconstructable = check_if_reconstructable("$simulation_dir/$new_sim_name.fam");
		
		## Fix cranefoot plot
		rewrite_config("$simulation_dir/$new_sim_name.config",$simulation_dir,$simulation_dir,$sim,$new_sim_name);
		run_cranefoot($simulation_dir,$new_sim_name);
		
		if($reconstructable == 0) # is this meant to happen?
		{
			#print "Test\n\n";
			run_system("rm $simulation_dir/$new_sim_name.ped");
			run_system("rm $simulation_dir/$new_sim_name.sex");
			run_system("rm $simulation_dir/$new_sim_name.cranefoot.ped");
			run_system("rm $simulation_dir/$new_sim_name.config");
			# run_system("rm $simulation_dir/$new_sim_name.genome");
			run_system("rm $simulation_dir/$new_sim_name.map");
			run_system("rm $simulation_dir/$new_sim_name.fam");
			#print "invalid $simulation_dir/$new_sim_name.fam\n";
			#exit;
		}

	}
	$Missing_ctr++;
	return $Missing_ctr;
}

sub run_cranefoot
{
	my $dir = shift;
	my $sim = shift;

    if ("$^O" eq "darwin") {
        run_system("$cranefoot_binary.mac $dir/$sim.config . > /dev/null 2>&1");
    } else {
        # I hope we're on 64-bit linux!
        run_system("$cranefoot_binary $dir/$sim.config . > /dev/null 2>&1");
    }

	run_system("mv $sim.ps $dir");
	run_system("rm $sim.topology.txt");

	## Implement ps2pdf 
	run_system("ps2pdf -dFirstPage=2 $dir/$sim.ps $dir/$sim.pdf");
	run_system("rm $dir/$sim.ps");

}

sub build_network_from_ped_file
{
	my $ped_file = shift;
	my $output_dir = shift;
	my $sim = shift;
	open(G,">$output_dir/$sim.genome");
	print G "RELATIONSHIP	FID1	IID1	FID2	IID2	RT	EZ	Z0	Z1	Z2	PI_HAT	PHE	DST	PPC	RATIO\n";

	## Read in file
	open(IN,$ped_file);
	my %all_nodes_network;
	my $network_ref = \%all_nodes_network;
	my @network_refs;
	
	## Build pedigree
	while(my $line = <IN>)
	{
		chomp($line);
		my ($FID,$IID,$PID,$MID,$SEX,$PHENOTYPE) = split(/\s+/,$line);
		my $child = "$IID";
		my $mom = "$MID";
		my $dad = "$PID";
		
		if(!exists $$network_ref{$child})
		{
			my $node = new PRIMUS::node_v7($child);
			$$network_ref{$child} = $node;
		}
		if(!exists $$network_ref{$dad} && $PID ne 0)
		{
			my $node = new PRIMUS::node_v7($dad);
			$$network_ref{$dad} = $node;
		}
		if(!exists $$network_ref{$mom} && $MID ne 0)
		{
			my $node = new PRIMUS::node_v7($mom);
			$$network_ref{$mom} = $node;
		}
		if($PID ne 0)
		{
			$$network_ref{$child}->add_parent($dad);
			$$network_ref{$dad}->add_child($child);
		}
		
		if($MID ne 0)
		{
			$$network_ref{$child}->add_parent($mom);
			$$network_ref{$mom}->add_child($child);
		}
	}
	close(IN);


	## Break network into individual pedigrees and
	## Build relationships
	my @keys = keys %$network_ref;
	print "# keys = " . @keys . "\n";
	while (@keys > 0)
	{
		my $node_name = @keys[0];
		
		my %names = $$network_ref{$node_name}->get_subpedigree_names($network_ref);
		my %pedigree;
		
		foreach my $node_name (keys %names)
		{
			$pedigree{$node_name} = $$network_ref{$node_name};
			delete $$network_ref{$node_name};
		}
		push(@network_refs,\%pedigree);
		@keys = keys %$network_ref;
		print "# keys = " . @keys . "\n";
		print "# network_refs = " . @network_refs . "\n";
	}
	

	
	## Write out relationships
	foreach my $network_ref (@network_refs)
	{
		foreach my $node_name (keys %$network_ref)
		{
			make_relative_network_from_pedigree($$network_ref{$node_name},$network_ref);
			my %rels = $$network_ref{$node_name}->relatives();
			foreach(keys %rels)
			{
				#print "$node_name -> $_  = ". $rels{$_}."\n";
			}
		}
	}
	
	return @network_refs;
}

sub load_ibd_estimates
{
	my $file = shift;

	## Loading data
	open(IN,$file) or die "$!\n";
	my $header = <IN>;
	while(my $line = <IN>)
	{
		chomp($line);
		my @temp = split(/\s+/,$line);
		
		my $k0 = @temp[6];
		my $k1 = @temp[7];
		my $k2 = @temp[8];
		
		my $relationship = @temp[14];
		
		if($relationship eq "PC"){push(@PC_k0s,$k0);push(@PC_k1s,$k1);}
		if($relationship eq "FS"){push(@FS_k0s,$k0);push(@FS_k1s,$k1);}
		if($relationship eq "HAG"){push(@HAG_k0s,$k0);push(@HAG_k1s,$k1);}
		if($relationship eq "CGH"){push(@CGH_k0s,$k0);push(@CGH_k1s,$k1);}
		if($relationship eq "DR"){push(@DR_k0s,$k0);push(@DR_k1s,$k1);}
		if($relationship eq "UN" && @UN_k1s <= 6000 && $k0 >= .95){push(@UN_k0s,$k0);push(@UN_k1s,$k1);}
	}
}

sub get_PC_IBD
{
	my $val = rand(@PC_k1s);
	my $k0 = @PC_k0s[$val];
	my $k1 = @PC_k1s[$val];
	my $k2 = 1 - @PC_k0s[$val] - @PC_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_FS_IBD
{
	my $val = rand(@FS_k1s);
	my $k0 = @FS_k0s[$val];
	my $k1 = @FS_k1s[$val];
	my $k2 = 1 - @FS_k0s[$val] - @FS_k1s[$val];
	return($k0,$k1,$k2);
}

sub get_HS_IBD
{
	my $val = rand(@HAG_k1s);
	my $k0 = @HAG_k0s[$val];
	my $k1 = @HAG_k1s[$val];
	my $k2 = 1 - @HAG_k0s[$val] - @HAG_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_GG_IBD
{
	my $val = rand(@HAG_k1s);
	my $k0 = @HAG_k0s[$val];
	my $k1 = @HAG_k1s[$val];
	my $k2 = 1 - @HAG_k0s[$val] - @HAG_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_AV_IBD
{
	my $val = rand(@HAG_k1s);
	my $k0 = @HAG_k0s[$val];
	my $k1 = @HAG_k1s[$val];
	my $k2 = 1 - @HAG_k0s[$val] - @HAG_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_GGG_IBD
{
	my $val = rand(@CGH_k1s);
	my $k0 = @CGH_k0s[$val];
	my $k1 = @CGH_k1s[$val];
	my $k2 = 1 - @CGH_k0s[$val] - @CGH_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_GAV_IBD
{
	my $val = rand(@CGH_k1s);
	my $k0 = @CGH_k0s[$val];
	my $k1 = @CGH_k1s[$val];
	my $k2 = 1 - @CGH_k0s[$val] - @CGH_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_HAV_IBD
{
	my $val = rand(@CGH_k1s);
	my $k0 = @CGH_k0s[$val];
	my $k1 = @CGH_k1s[$val];
	my $k2 = 1 - @CGH_k0s[$val] - @CGH_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_1C_IBD
{
	my $val = rand(@CGH_k1s);
	my $k0 = @CGH_k0s[$val];
	my $k1 = @CGH_k1s[$val];
	my $k2 = 1 - @CGH_k0s[$val] - @CGH_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_DIS_IBD
{
	my $val = rand(@DR_k1s);
	my $k0 = @DR_k0s[$val];
	my $k1 = @DR_k1s[$val];
	my $k2 = 1 - @DR_k0s[$val] - @DR_k1s[$val];
	return($k0,$k1,$k2);
}
sub get_UN_IBD
{
	my $val = rand(@UN_k1s);
	my $k0 = @UN_k0s[$val];
	my $k1 = @UN_k1s[$val];
	my $k2 = 1 - @UN_k0s[$val] - @UN_k1s[$val];
	return($k0,$k1,$k2);
}


sub make_relative_network_from_pedigree {	
	my $self = shift;
	my $network_ref = shift;
    	my $phase = 10;
    	my $self_name = $self->name();
   	if($self_name =~/Missing/i){next;}
    	my @parents = @{$self->{PARENTS} };
   	my @children = @{$self->{CHILDREN} };
    	my %rels = % {$self->{RELATIVES} };
    
    	my @PC_probs = (1,0,0,0,0,0);
	my @FS_probs = (0,1,0,0,0,0);
	my @HAG_probs = (0,0,1,0,0,0);
	my @CGH_probs = (0,0,0,1,0,0);
	my @D1C_probs = (0,0,0,2,0,0);
	my @HS1C_probs = (0,0,1,1,0,0);
	my @DR_probs = (0,0,0,0,1,0);
	my @UN_probs = (0,0,0,0,0,1);


	## Check if rel is child of self
    	foreach(@children)
    	{
		if($_ =~ /Missing/i){next;}
		$$network_ref{$self_name}->add_relative($_,@PC_probs);
    	}    	
   		
	if ($phase >= 2)
	{
		## Check if rel is grandchild of self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
				$$network_ref{$self_name}->add_relative($_,@HAG_probs);
			}
		}
	}
	
	if ($phase >= 3)
	{
		## Check if rel is great grandchild self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				my @GG_children = $$network_ref{$_}->children();
				foreach(@GG_children) # These are neices/nephews
				{
					if($_ =~ /Missing/i){next;}
					$$network_ref{$self_name}->add_relative($_,@CGH_probs);
				}
			}
		}
	}
	
	#print "parents: @parents\n";
    	if(@parents ne 0)
    	{
    		## Check parents
 		foreach(@parents)
    		{
    			if($_ =~/Missing/i){next;}
    			$$network_ref{$self_name}->add_relative($_,@PC_probs);
    		}
    		
    		## Check full siblings
    		my @P1_children = $$network_ref{@parents[0]}->children();
    		my @P2_children = $$network_ref{@parents[1]}->children();
    		#print "@parents[0] children: @P1_children\n";
    		#print "@parents[1] children: @P2_children\n";
		my @FS;
		my @HS;
		my @UNCLES_AUNTS;
 
		foreach(@P1_children)
    		{
    			my $child = $_;
    			if(grep($_ eq $child,@P2_children)) #full sibling
    			{
				push(@FS,$child);
    				if($child =~/Missing/i || $child eq $self_name){next;}
				$$network_ref{$self_name}->add_relative($child,@FS_probs);
    			}
    		}
    		

    		## Check half siblings
    		foreach(@P1_children)
    		{
    			my $child = $_;
    			if(!grep($_ eq $child, @P2_children)) # Half sib
    			{
				push(@HS,$child);
    				if($child =~ /Missing/i){next;}
				$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    			}
    		}
    		foreach(@P2_children)
    		{
    			my $child = $_;
    			if(!grep($_ eq $child, @P1_children)) # Half sib
    			{
				push(@HS,$child);
    				if($child =~ /Missing/i){next;}
				$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    			}
    		}

		## Check Nieces/nephews
		foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
    				$$network_ref{$self_name}->add_relative($_,@HAG_probs);
			}
		}


    		## Check grandparents
    		my @P1_parents = $$network_ref{@parents[0]}->parents();
    		my @P2_parents = $$network_ref{@parents[1]}->parents();
    		#print "@parents[0] parents: @P1_parents\n";
    		#print "@parents[1] parents: @P2_parents\n";
		foreach(@P1_parents)
    		{
    			if($_ =~/Missing/i){next;}
    			$$network_ref{$self_name}->add_relative($_,@HAG_probs);
    		}
		foreach(@P2_parents)
		{
    			if($_ =~/Missing/i){next;}
    			$$network_ref{$self_name}->add_relative($_,@HAG_probs);
    		}
    		
    		
    		## Check if $rel is uncle/aunt of $self
    		if(@P1_parents > 0)
    		{
    			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
    			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
    			#print "P1 @P1_parents[0] children: @G1_1_children\n";
    			#print "P1 @P1_parents[1] children: @G1_2_children\n";
		 
 			foreach(@G1_1_children)
    			{
    				my $child = $_;
    				if(grep($_ eq $child,@G1_2_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
					push(@UNCLES_AUNTS,$child);
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
					$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    				}
    			}
    		}
    		if(@P2_parents > 0)
    		{
    			#print "$self_name parent @parents[1] parents @P2_parents\n";
    			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
    			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
    			#print "P2 @P2_parents[0] children: @G2_1_children\n";
    			#print "P2 @P2_parents[1] children: @G2_2_children\n";
		 
 			foreach(@G2_1_children)
    			{
    				my $child = $_;
    				if(grep($_ eq $child,@G2_2_children)) #Unclde/aunt
    				{
					push(@UNCLES_AUNTS,$child);    				
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
					$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    				}
    			}
    		}
    		
    		## Check half uncle/aunts
    		if(@P1_parents > 0)
    		{
    			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
    			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
    			#print "P1 @P1_parents[0] children: @G1_1_children\n";
    			#print "P1 @P1_parents[1] children: @G1_2_children\n";
		 
 			foreach(@G1_1_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G1_2_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
    				}
    			}
    			foreach(@G1_2_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G1_1_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
    				}
    			}
    		}
    		if(@P2_parents > 0)
    		{
    			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
    			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
    			#print "P2 @P2_parents[0] children: @G2_1_children\n";
    			#print "P2 @P2_parents[1] children: @G2_2_children\n";
    			
 			foreach(@G2_1_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G2_2_children)) #Unclde/aunt
    				{
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
    				}
    			}
			foreach(@G2_2_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G2_1_children)) #Unclde/aunt
    				{
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
    				}
    			}
    		}
    			
    		## Check if $rel is half neice/nephew of $self
		foreach(@HS)
		{
			my @HS_children = $$network_ref{$_}->children();
			foreach(@HS_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
			}
		}
				
				
		## Check if $rel is GREAT-grandparent of $self
    		my @grandparents;
    		push(@grandparents, @P1_parents);
    		push(@grandparents, @P2_parents);
    		#my @P1_parents = $$network_ref{@parents[0]}->parents();
    		#my @P2_parents = $$network_ref{@parents[1]}->parents();
    		#print "@parents[0] parents: @P1_parents\n";
    		#print "@parents[1] parents: @P2_parents\n";
    		
    		#print "$self_name grandparents: @grandparents\n";
		foreach(@grandparents)
		{
			my @GG_parents = $$network_ref{$_}->parents();
			foreach(@GG_parents)
    			{
    				if($_ =~/Missing/i){next;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
    			}
		}
   	 		
   	 	## CHECK if COUSINS
		if($self_name eq "E" && @UNCLES_AUNTS[0] eq "F")
   	 	{
			#print "$self_name\'s uncles_aunts: @UNCLES_AUNTS#######################################\n";
		}
		foreach(@UNCLES_AUNTS)
		{
			my @cousins = $$network_ref{$_}->children();
			foreach(@cousins)
			{
				if($_ =~ /Missing/i){next;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
			}
		}
    	    
    	    ## CHECK IF REL IS GREAT UNCLE/AUNT of self
    	    my @G_parents;
    	    push(@G_parents, @P1_parents);
    	    push(@G_parents, @P2_parents);
    	    
    	    foreach my $G_parent (@G_parents)
    	    {
    	    	my @GG_parents = $$network_ref{$G_parent}->parents();
    	    	if(@GG_parents eq 0){next;}
    	    	
    	    	my @GG1_children = $$network_ref{@GG_parents[0]}->children();
    			my @GG2_children = $$network_ref{@GG_parents[1]}->children();
    			#print "GG1 @GG_parents[0] children: @GG1_children\n";
    			#print "GG2 @GG_parents[1] children: @GG2_children\n";
		 
 			foreach my $child (@GG1_children)
    			{
    				if(grep($_ eq $child,@GG2_children)) #great uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq $G_parent){next;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
    				}
    			}
    	    }
    	      	    
    	    ## Check if rel is great neice or nephew of self
	foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				my @FS_Gchildren = $$network_ref{$_}->children();
				foreach(@FS_Gchildren) # These are neices/nephews
				{
					if($_ =~ /Missing/i){next;}
    					$$network_ref{$self_name}->add_relative($_,@CGH_probs);
				}
			}
		}
    		
    	}
	
	## Check for more distant relationships
	my %all_relatives = $self->get_all_relatives($network_ref);
	my %current_rels = $self->relatives();
	
	foreach(keys %all_relatives)
	{
		if(!exists $current_rels{$_} && $_ ne $self_name)
		{
			if($_ =~ /Missing/i){next;}
			$$network_ref{$self_name}->add_relative($_,@DR_probs);
		}
	}

	## get all unrelateds
	my %all_pedigree = $self->get_subpedigree_names($network_ref);
	foreach(keys %all_pedigree)
	{
		#print "$self_name->node $_\n";
		if(!exists $all_relatives{$_})
		{
			if($_ =~ /Missing/i){next;}
			$$network_ref{$self_name}->add_relative($_,@UN_probs);
		}
	}
}



sub make_relative_network_from_pedigree_OLD
{	
	my $self = shift;
	my $network_ref = shift;
    	my $phase = 10;
    	my $self_name = $self->name();
   	if($self_name =~/Missing/i){next;}
    	my @parents = @{$self->{PARENTS} };
   	my @children = @{$self->{CHILDREN} };
    	my %rels = % {$self->{RELATIVES} };
    	my $num_fails = 0;
    
    	my @PC_probs = (1,0,0,0,0,0);
	my @FS_probs = (0,1,0,0,0,0);
	my @HAG_probs = (0,0,1,0,0,0);
	my @CGH_probs = (0,0,0,1,0,0);
	my @D1C_probs = (0,0,0,2,0,0);
	my @HS1C_probs = (0,0,1,1,0,0);
	my @DR_probs = (0,0,0,0,1,0);
	my @UN_probs = (0,0,0,0,0,1);
	
	## Check if rel is child of self
    	foreach(@children)
    	{
		if($_ =~ /Missing/i){next;}
		$$network_ref{$self_name}->add_relative($_,@PC_probs);
		if(exists $pairwise_rels{$self_name}{$_}){next;}
		$pairwise_rels{$self_name}{$_} = 1;
		$pairwise_rels{$_}{$self_name} = 1;
		my @ibd = get_PC_IBD();
		print G "fam1 $self_name fam1 $_ PC NA @ibd\n";
    	}    	
   		
	if ($phase >= 2)
	{
		## Check if rel is grandchild of self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
				#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not grandchild of $self_name\n";$num_fails++;}
				my @ibd = get_GG_IBD();
				$$network_ref{$self_name}->add_relative($_,@HAG_probs);
				if(exists $pairwise_rels{$self_name}{$_}){next;}
				$pairwise_rels{$self_name}{$_} = 1;
				$pairwise_rels{$_}{$self_name} = 1;
				print G "fam1 $self_name fam1 $_ GG NA @ibd\n";
			}
		}
	}
	
	if ($phase >= 3)
	{
		## Check if rel is great grandchild self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				my @GG_children = $$network_ref{$_}->children();
				foreach(@GG_children) # These are neices/nephews
				{
					if($_ =~ /Missing/i){next;}
					#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not great-grandchild of $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($_,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$_}){next;}
					$pairwise_rels{$self_name}{$_} = 1;
					$pairwise_rels{$_}{$self_name} = 1;
					my @ibd = get_GGG_IBD();
					print G "fam1 $self_name fam1 $_ GGG NA @ibd\n";
				}
			}
		}
	}
	
	#print "parents: @parents\n";
    	if(@parents ne 0)
    	{
    		## Check parents
 		foreach(@parents)
    		{
    			if($_ =~/Missing/i){next;}
    			#if($rels{$_} ne "PC"){print "FAIL!!! $_ not PC with $self_name. it is " . $rels{$_} ."\n";$num_fails++;}
    			$$network_ref{$self_name}->add_relative($_,@PC_probs);
			
			if(exists $pairwise_rels{$self_name}{$_}){next;}
			$pairwise_rels{$self_name}{$_} = 1;
			$pairwise_rels{$_}{$self_name} = 1;
			my @ibd = get_PC_IBD();
			print G "fam1 $self_name fam1 $_ PC NA @ibd\n";
    		}
    		
    		## Check full siblings
    		my @P1_children = $$network_ref{@parents[0]}->children();
    		my @P2_children = $$network_ref{@parents[1]}->children();
    		#print "@parents[0] children: @P1_children\n";
    		#print "@parents[1] children: @P2_children\n";
		my @FS;
		my @HS;
		my @UNCLES_AUNTS;
 
		foreach(@P1_children)
    		{
    			my $child = $_;
    			if(grep($_ eq $child,@P2_children)) #full sibling
    			{
    				if($child =~/Missing/i || $child eq $self_name){next;}
    				#if($rels{$child} ne "FS"){print "FAIL!!! $child not FS with $self_name\n";$num_fails++;}
    				#else
    				#{
    					$$network_ref{$self_name}->add_relative($child,@FS_probs);
    					push(@FS,$child);
					
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_FS_IBD();
					print G "fam1 $self_name fam1 $child FS NA @ibd\n";
    				#}
    			}
    		}
    		
    		## Don't do the phase 2 checks if we are not in Phase 2
    		if ($phase < 2){
    			return $num_fails;
    		}

    		
    		## Check half siblings
    		foreach(@P1_children)
    		{
    			my $child = $_;
    			if(!grep($_ eq $child, @P2_children)) # Half sib
    			{
    				if($child =~ /Missing/i){next;}
    				#if($rels{$child} ne "HAG"){print "FAIL!!! $child not HS with $self_name\n";$num_fails++;}
    				#else
    				#{
    					$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    					push(@HS,$child);
					
					if(exists $pairwise_rels{$self_name}{$_}){next;}
					$pairwise_rels{$self_name}{$_} = 1;
					$pairwise_rels{$_}{$self_name} = 1;
					my @ibd = get_HS_IBD();
					print G "fam1 $self_name fam1 $_ HS NA @ibd\n";
    				#}
    			}
    		}
    		foreach(@P2_children)
    		{
    			my $child = $_;
    			if(!grep($_ eq $child, @P1_children)) # Half sib
    			{
    				if($child =~ /Missing/i){next;}
    				#if($rels{$child} ne "HAG"){print "FAIL!!! $child not HS with $self_name\n";$num_fails++;}
    				#else
    				#{
    					$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    					push(@HS,$child);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_HS_IBD();
					print G "fam1 $self_name fam1 $_ HS NA @ibd\n";
    				#}
    			}
    		}

		## Check Nieces/nephews
		foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
    				#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not neices/nephew with $self_name\n";$num_fails++;}
    				$$network_ref{$self_name}->add_relative($_,@HAG_probs);
				if(exists $pairwise_rels{$self_name}{$_}){next;}
				$pairwise_rels{$self_name}{$_} = 1;
				$pairwise_rels{$_}{$self_name} = 1;
				my @ibd = get_AV_IBD();
				print G "fam1 $self_name fam1 $_ AV NA @ibd\n";
			}
		}


    		## Check grandparents
    		my @P1_parents = $$network_ref{@parents[0]}->parents();
    		my @P2_parents = $$network_ref{@parents[1]}->parents();
    		#print "@parents[0] parents: @P1_parents\n";
    		#print "@parents[1] parents: @P2_parents\n";
		foreach(@P1_parents)
    		{
    			if($_ =~/Missing/i){next;}
    			#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not GG with $self_name\n";$num_fails++;}
			if(exists $pairwise_rels{$self_name}{$_}){next;}
			$pairwise_rels{$self_name}{$_} = 1;
			$pairwise_rels{$_}{$self_name} = 1;
			my @ibd = get_GG_IBD();
			print G "fam1 $self_name fam1 $_ GG NA @ibd\n";
    			$$network_ref{$self_name}->add_relative($_,@HAG_probs);
    		}
		foreach(@P2_parents)
		{
    			if($_ =~/Missing/i){next;}
    			#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not GG with $self_name\n";$num_fails++;}
    			$$network_ref{$self_name}->add_relative($_,@HAG_probs);
			if(exists $pairwise_rels{$self_name}{$_}){next;}
			$pairwise_rels{$self_name}{$_} = 1;
			$pairwise_rels{$_}{$self_name} = 1;
			my @ibd = get_GG_IBD();
			print G "fam1 $self_name fam1 $_ GG NA @ibd\n";
    		}
    		
    		
    		## Check if $rel is uncle/aunt of $self
    		if(@P1_parents > 0)
    		{
    			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
    			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
    			#print "P1 @P1_parents[0] children: @G1_1_children\n";
    			#print "P1 @P1_parents[1] children: @G1_2_children\n";
		 
 			foreach(@G1_1_children)
    			{
    				my $child = $_;
    				if(grep($_ eq $child,@G1_2_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
    					#if($rels{$child} ne "HAG"){print "FAIL!!! $self_name not neice or nephew1 of $child\n";$num_fails++;}
    					#else
    					#{
    						$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    						push(@UNCLES_AUNTS,$child);
						
						if(exists $pairwise_rels{$self_name}{$child}){next;}
						$pairwise_rels{$self_name}{$child} = 1;
						$pairwise_rels{$child}{$self_name} = 1;
						my @ibd = get_AV_IBD();
						print G "fam1 $self_name fam1 $_ AV NA @ibd\n";
    					#}
    				}
    			}
    		}
    		if(@P2_parents > 0)
    		{
    			#print "$self_name parent @parents[1] parents @P2_parents\n";
    			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
    			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
    			#print "P2 @P2_parents[0] children: @G2_1_children\n";
    			#print "P2 @P2_parents[1] children: @G2_2_children\n";
		 
 			foreach(@G2_1_children)
    			{
    				my $child = $_;
    				if(grep($_ eq $child,@G2_2_children)) #Unclde/aunt
    				{
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
    					#if($rels{$child} ne "HAG"){print "FAIL!!! $self_name not neice or nephew2 of $child\n";$num_fails++;}
    					#else
    					#{
    						$$network_ref{$self_name}->add_relative($child,@HAG_probs);
    						push(@UNCLES_AUNTS,$child);    				
						
						if(exists $pairwise_rels{$self_name}{$child}){next;}
						$pairwise_rels{$self_name}{$child} = 1;
						$pairwise_rels{$child}{$self_name} = 1;
						my @ibd = get_AV_IBD();
						print G "fam1 $self_name fam1 $_ AV NA @ibd\n";
    					#}
    				}
    			}
    		}
    		
    		
    		## Don't do the phase 3 checks if we are not in Phase 3
    		if ($phase < 3){
    			return $num_fails;
    		}
    		
    		## Check half uncle/aunts
    		if(@P1_parents > 0)
    		{
    			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
    			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
    			#print "P1 @P1_parents[0] children: @G1_1_children\n";
    			#print "P1 @P1_parents[1] children: @G1_2_children\n";
		 
 			foreach(@G1_1_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G1_2_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
    					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_HAV_IBD();
					print G "fam1 $self_name fam1 $_ HAV NA @ibd\n";
    				}
    			}
    			foreach(@G1_2_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G1_1_children)) #uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq @parents[0]){next;}
    					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_HAV_IBD();
					print G "fam1 $self_name fam1 $_ HAV NA @ibd\n";
    				}
    			}
    		}
    		if(@P2_parents > 0)
    		{
    			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
    			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
    			#print "P2 @P2_parents[0] children: @G2_1_children\n";
    			#print "P2 @P2_parents[1] children: @G2_2_children\n";
    			
 			foreach(@G2_1_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G2_2_children)) #Unclde/aunt
    				{
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
    					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_HAV_IBD();
					print G "fam1 $self_name fam1 $_ HAV NA @ibd\n";
    				}
    			}
			foreach(@G2_2_children)
    			{
    				my $child = $_;
    				if(!grep($_ eq $child,@G2_1_children)) #Unclde/aunt
    				{
    					if($child =~/Missing/i || $_ eq @parents[1]){next;}
    					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_HAV_IBD();
					print G "fam1 $self_name fam1 $_ HAV NA @ibd\n";
    				}
    			}
    		}
    			
    		## Check if $rel is half neice/nephew of $self
		foreach(@HS)
		{
			my @HS_children = $$network_ref{$_}->children();
			foreach(@HS_children) # These are neices/nephews
			{
				if($_ =~ /Missing/i){next;}
    				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not half neices/nephew with $self_name\n";$num_fails++;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
				if(exists $pairwise_rels{$self_name}{$_}){next;}
				$pairwise_rels{$self_name}{$_} = 1;
				$pairwise_rels{$_}{$self_name} = 1;
				my @ibd = get_HAV_IBD();
				print G "fam1 $self_name fam1 $_ HAV NA @ibd\n";
			}
		}
				
				
		## Check if $rel is GREAT-grandparent of $self
    		my @grandparents;
    		push(@grandparents, @P1_parents);
    		push(@grandparents, @P2_parents);
    		#my @P1_parents = $$network_ref{@parents[0]}->parents();
    		#my @P2_parents = $$network_ref{@parents[1]}->parents();
    		#print "@parents[0] parents: @P1_parents\n";
    		#print "@parents[1] parents: @P2_parents\n";
    		
    		#print "$self_name grandparents: @grandparents\n";
		foreach(@grandparents)
		{
			my @GG_parents = $$network_ref{$_}->parents();
			foreach(@GG_parents)
    			{
    				if($_ =~/Missing/i){next;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
    				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not GG-parent of $self_name\n";$num_fails++;}
				if(exists $pairwise_rels{$self_name}{$_}){next;}
				$pairwise_rels{$self_name}{$_} = 1;
				$pairwise_rels{$_}{$self_name} = 1;
				my @ibd = get_GGG_IBD();
				print G "fam1 $self_name fam1 $_ GGG NA @ibd\n";
    			}
		}
   	 		
   	 		## CHECK if COUSINS
		if($self_name eq "E" && @UNCLES_AUNTS[0] eq "F")
   	 	{
			#print "$self_name\'s uncles_aunts: @UNCLES_AUNTS#######################################\n";
		}
		foreach(@UNCLES_AUNTS)
		{
			my @cousins = $$network_ref{$_}->children();
			foreach(@cousins)
			{
				if($_ =~ /Missing/i){next;}
    				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not cousin with $self_name\n";$num_fails++;}
    				$$network_ref{$self_name}->add_relative($_,@CGH_probs);
				
				if(exists $pairwise_rels{$self_name}{$_}){next;}
				$pairwise_rels{$self_name}{$_} = 1;
				$pairwise_rels{$_}{$self_name} = 1;
				my @ibd = get_1C_IBD();
				print G "fam1 $self_name fam1 $_ 1C NA @ibd\n";
			}
		}
    	    
    	    ## CHECK IF REL IS GREAT UNCLE/AUNT of self
    	    my @G_parents;
    	    push(@G_parents, @P1_parents);
    	    push(@G_parents, @P2_parents);
    	    
    	    foreach my $G_parent (@G_parents)
    	    {
    	    	my @GG_parents = $$network_ref{$G_parent}->parents();
    	    	if(@GG_parents eq 0){next;}
    	    	
    	    	my @GG1_children = $$network_ref{@GG_parents[0]}->children();
    			my @GG2_children = $$network_ref{@GG_parents[1]}->children();
    			#print "GG1 @GG_parents[0] children: @GG1_children\n";
    			#print "GG2 @GG_parents[1] children: @GG2_children\n";
		 
 				foreach my $child (@GG1_children)
    			{
    				if(grep($_ eq $child,@GG2_children)) #great uncle/aunt
    				{
    					#print "$child : $self_name\n";
    					if($child =~ /Missing/i || $child eq $G_parent){next;}
    					#if($rels{$child} ne "CGH"){print "FAIL!!! $self_name not great neice or nephew of $child\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($child,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$child}){next;}
					$pairwise_rels{$self_name}{$child} = 1;
					$pairwise_rels{$child}{$self_name} = 1;
					my @ibd = get_GAV_IBD();
					print G "fam1 $self_name fam1 $child GAV NA @ibd\n";
    				}
    			}
    	    }
    	      	    
    	    ## Check if rel is great neice or nephew of self
	foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				my @FS_Gchildren = $$network_ref{$_}->children();
				foreach(@FS_Gchildren) # These are neices/nephews
				{
					if($_ =~ /Missing/i){next;}
    					#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not great neices/nephew with $self_name\n";$num_fails++;}
    					$$network_ref{$self_name}->add_relative($_,@CGH_probs);
					if(exists $pairwise_rels{$self_name}{$_}){next;}
					$pairwise_rels{$self_name}{$_} = 1;
					$pairwise_rels{$_}{$self_name} = 1;
					my @ibd = get_GAV_IBD();
					print G "fam1 $self_name fam1 $_ GAV NA @ibd\n";
				}
			}
		}
    		
    	}
	
	## Check for more distant relationships
	my %all_relatives = $self->get_all_relatives($network_ref);
	my %current_rels = $self->relatives();
	
	foreach(keys %all_relatives)
	{
		if(!exists $current_rels{$_} && $_ ne $self_name)
		{
			if($_ =~ /Missing/i){next;}
			$$network_ref{$self_name}->add_relative($_,@DR_probs);
			if(exists $pairwise_rels{$self_name}{$_}){next;}
			$pairwise_rels{$self_name}{$_} = 1;
			$pairwise_rels{$_}{$self_name} = 1;
			my @ibd = get_DIS_IBD();
			print G "fam1 $self_name fam1 $_ DIS NA @ibd\n";
		}
	}

	## get all unrelateds
	my %all_pedigree = $self->get_subpedigree_names($network_ref);
	foreach(keys %all_pedigree)
	{
		print "$self_name->node $_\n";
		if(!exists $all_relatives{$_})
		{
			print "here\n";
			if($_ =~ /Missing/i){next;}
			$$network_ref{$self_name}->add_relative($_,@UN_probs);
			if(exists $pairwise_rels{$self_name}{$_}){next;}
			$pairwise_rels{$self_name}{$_} = 1;
			$pairwise_rels{$_}{$self_name} = 1;
			my @ibd = get_UN_IBD();
			print G "fam1 $self_name fam1 $_ UN NA @ibd\n";
		}
	}
	return $num_fails;
}


sub make_simulated_pedigree {

	my $out_dir = shift;
	my $ONEKG_pop = shift;
	my $sim_name = shift;
	my $type = shift;

	my $sim_dir = "$out_dir/$sim_name";
	if(!-d $sim_dir){run_system("mkdir $sim_dir")}

	print "Generating simulation: $sim_dir/$sim_name\n\n";
	## Reset all variables
	reset_variables();
	
	## Prime the array
	my $f1 = add_sample(0,0,1,1);
	push(@samples_to_visit,$f1);
	
	#my $ORIGINAL_MAX_SAMPLES = $MAX_SAMPLES;

	## Alter max_samples to allow for added halfsibs
	if($type eq "halfsib")
	{
		my $num_parents = $MAX_SAMPLES
	}

	## Build pedigree
	my $curr_ID = shift(@samples_to_visit);
	my $complex_added = 0;
	while($curr_ID && $generations{$curr_ID} < $MAX_GENERATIONS && $sample_ctr < $MAX_SAMPLES)
	{
		print "$curr_ID\n";
		my @children = add_children($curr_ID);
		push(@samples_to_visit,@children);
		
		## Add half sibs
		if($type =~ /halfsib/ || $type =~ /complex/)
		{
			my $val = rand();
			if($val <= $HALF_SIB_RATE)
			{
				if($sample_ctr + $NUM_HS + 1 < $MAX_SAMPLES)
				{
					my @children = add_children($curr_ID,-1,$NUM_HS);
					push(@samples_to_visit,@children);
				}
				elsif($sample_ctr + $NUM_HS < $MAX_SAMPLES)
				{
					my @children = add_children($curr_ID,-1,($NUM_HS-1));
					push(@samples_to_visit,@children);
				}
			}
		}

		if($type =~ /complex/ && $sample_ctr > $MAX_SAMPLES - 15 && $complex_added == 0)
		{
			my @children = add_complex_relationship();
			print "ADDED COMPLEX\n";
			push(@samples_to_visit,@children);
			$complex_added = 1;
		}


		print "samples to visit: @samples_to_visit\n";
		$curr_ID = shift(@samples_to_visit);
		print "sample: $curr_ID = $generations{$curr_ID}\n";
		#exit;
	}

	## OLD ADD HALD SIB
	# Break up homes (add half siblings)
	if($type eq "halfsib")
	{
		my @samples_w_children = keys %children;
		my $num_halves = @samples_w_children/4;
	
		for(my $i = 0; $i < $num_halves; $i++)
		{
			my $val = rand(@samples_w_children);
			my $sample = @samples_w_children[$val];
			if($samples{$sample} != 2)
			{
				add_children($sample)
			}
			else
			{
				add_children("",$sample)
			}
	
		}
	}

	## Add complex relationships
	#add_complex_relationship();

	## Print simulation data - can probably stop this from being printed to the console
	print_fam_file($sim_name,$sim_dir,$sim_name);
	write_sex_file($sim_name,$sim_dir,$sim_name);
	write_cranefoot_files($sim_name,$sim_dir);

	add_genotypes("$sim_dir/$sim_name", $ONEKG_pop, $parallel_status);

	my $freq_file = "$data_dir/all.frq";
	my $mega_bim = "$data_dir/extract_mega.bim";
	
	### NEW: filter vcf file ahead of time to stop all the downstream files from being massive
	run_system("$plink2_binary --vcf $sim_dir/${sim_name}_all_chr.vcf.gz --geno 0.1 --maf 0.05 --mind 0.1 --extract ${mega_bim} --double-id --allow-extra-chr --set-missing-var-ids \@:\# --export vcf bgz --out $sim_dir/${sim_name}_all_chr_qced > /dev/null 2>&1"); 	
	print("\nQCed output VCF file\n\n");

	## get a list of snps to the location you're running the script 
	run_system("$plink_binary --genome --vcf $sim_dir/${sim_name}_all_chr_qced.vcf.gz --const-fid --write-snplist --out $sim_dir/${sim_name}_plink > /dev/null 2>&1"); 	

	## get duplicates from this list and write to new file
	run_system("cat $sim_dir/${sim_name}_plink.snplist | sort | uniq -d > $sim_dir/${sim_name}_dupsnps.txt > /dev/null 2>&1");

	## run main command with dup exclusion list -- this builds ped and map 
	run_system("$plink_binary --genome --vcf $sim_dir/${sim_name}_all_chr_qced.vcf.gz --out $sim_dir/$sim_name --biallelic-only --snps-only \"just-acgt\" --read-freq $freq_file --const-fid --exclude $sim_dir/${sim_name}_dupsnps.txt --set-missing-var-ids \@:\# > /dev/null 2>&1");

	## Generate plink-binary files -- ADD FILTERING STEPS INTO HERE TOO 
	run_system("$plink_binary --vcf $sim_dir/${sim_name}_all_chr_qced.vcf.gz --make-bed --out $sim_dir/${sim_name} > /dev/null 2>&1");
	print "\nBuilt PLINK-binary output files\n";

	# overwrite the new .fam file with the .diag.txt.fam file since it has parental IDs and sex in it 
	run_system("cp $sim_dir/${sim_name}_diag.txt.fam $sim_dir/${sim_name}.fam");

	## Generate ped and map files for good measure
	run_system("$plink_binary --bfile $sim_dir/${sim_name} --recode --out $sim_dir/${sim_name} > /dev/null 2>&1"); 

	print("\nDone with all PLINK steps\n\n");

	# print "\nSamples w/ children: " .@samples_w_children ."\n";
}

sub add_complex_relationship
{
	my @unrooted_samples = keys %married_into_family;
	my @all_samples = keys %samples;

	## Select sample that married into the family
	my $p1;
	my $p1_generation;
	my $p1_gender;
	
	while(!$p1 || $p1_generation == 1)
	{
		my $val = rand(@unrooted_samples);
		$p1 = @unrooted_samples[$val];
		$p1_generation = $generations{$p1};
		$p1_gender = $samples{$p1};
		print "$p1\n";
	}

	my $p2;
	my $p2_generation;
	my $p2_gender;

	## Both must be in the same generation and an opposite gender
	while(!$p2 || exists $married_into_family{$p2} || $p2_generation != $p1_generation || $p2_gender eq $p1_gender)
	{
		my $val2 = rand(@all_samples);
		$p2 = @all_samples[$val2];
		$p2_generation = $generations{$p2};
		$p2_gender = $samples{$p2};
		print "$p2\n";
	}
	if($p1_gender != 2)
	{
		my @children = add_children($p1,$p2,$NUM_HS);
		return @children;
	}
	elsif($p2_gender != 2)
	{
		my @children = add_children($p2,$p1,$NUM_HS);
		return @children;
	}
	else
	{
		die "ERROR!!!\n";
	}
}

sub add_sample
{
	my $dad = shift;
	my $mom = shift;
	my $gender = shift;
	my $generation = shift;
	$sample_ctr++;
	my $sample = "id$sample_ctr";
	
	$generations{$sample} = $generation;
	$samples{$sample} = $gender;
	@{ $parents{$sample} } = ($dad,$mom);
	return $sample;
}

sub add_children
{
	my $dad = @_[0];
	my $mom = @_[1];
	my $num_children = @_[2];
	my $generation;
	
	if($dad){$generation = 1 + $generations{$dad};}
	else{$generation = 1 + $generations{$mom};}
	
	if(!$num_children)
	{
		#$num_children = random_poisson(1,$mean_children);
		$num_children = $mean_children;
	}
	
	## Make the founders have three kids
	if($generation <= 2 && $num_children < 2){$num_children = $mean_children;}
	
	#print "num children: $num_children\n";
	
	## Don't proceed if there are no kids to be added
	if($num_children == 0){return};

	################################################################
	
	## Don't allow making families larger than $MAX_SAMPLES
	# adding in another round of kids requires at least one kid and one mom, so if you are already one away from max, just at 1 extra kid
	# if(($sample_ctr + $num_children + 2) == $MAX_SAMPLES)
	# {
	# 	print "HERE1\n";
	# 	print "sample_ctr: $sample_ctr; $num_children\n";
	# 	$num_children++;
	# }
	# while(($sample_ctr + $num_children + 1) > $MAX_SAMPLES)
	# {
	# 	print "HERE2\n";
	# 	print "sample_ctr: $sample_ctr; $num_children\n";
	# 	if($num_children == 0){die "$num_children is 0; didn't do add 4 chilren in previous add\n";}
	# 	$num_children--;
	# }
	# print "sample_ctr: $sample_ctr; $num_children\n";

	################################################################

	### UPDATE 1/30/25

	## Don't allow making families larger than $MAX_SAMPLES
	# If we're near the max sample cutoff, adjust number of children accordingly
	if ($sample_ctr >= $MAX_SAMPLES) {
		return (); # Return empty array if we've hit max samples
	}

	# Calculate how many more samples we can add
	my $remaining_samples = $MAX_SAMPLES - $sample_ctr;

	# We need at least 1 spot for a potential new parent
	$remaining_samples--;

	# Adjust number of children if needed
	if ($num_children > $remaining_samples) {
		$num_children = $remaining_samples;
	}

	# If we can't add any children, return empty array
	if ($num_children <= 0) {
		return ();
	}

	print "sample_ctr: $sample_ctr; $num_children\n";

	################################################################

	## If mom does not exist, make her
	if(!$mom || $mom == -1){$mom = add_sample(0,0,2,$generations{$dad});$married_into_family{$mom}=1;}
	if(!$dad || $mom == -1){$dad = add_sample(0,0,1,$generations{$mom});$married_into_family{$dad}=1;}

	if($samples{$dad} == 0){$samples{$dad} = 1;} 
	if($samples{$mom} == 0){$samples{$mom} = 2;} 
	if($samples{$dad} != 1 || $samples{$mom} != 2){die "ERROR!!! dad $dad($samples{$dad}) and/or mom $mom($samples{$mom}) genders invalid\n";}
	#print "$generations{$mom}: $generation\n";
	my @children;
	for(my $i = 0; $i < $num_children; $i++)
	{
		my $child = add_sample($dad,$mom,0,$generation);
		push(@{ $children{$dad} }, $child);
		push(@{ $children{$mom} }, $child);
		push(@children,$child);
	}
	return(@children);
}

sub print_fam_file
{
	my $file = shift;
	my $out_dir = shift;
	my $FID = shift;
	
	open(PED,">$out_dir/$file.fam");
	foreach my $sample (sort {substr($a,2) <=> substr($b,2)} keys %parents)
	{
		my ($dad,$mom) = @{$parents{$sample} };
		my $sex = $samples{$sample};
		print PED "$FID\t$sample\t$dad\t$mom\t$sex\t0\n";
	}
	close(PED);
}

sub write_sex_file
{
	my $sim_name = shift;
	my $out_dir = shift;

	open(SEX,">$out_dir/$sim_name.sex");

	foreach my $sample (keys %samples)
	{
		my $sex = $samples{$sample};
		if($sex == 0){$sex = 2;}
		print SEX "$sim_name\t$sample\t$sex\n";
	}
	close(SEX);
}

sub write_cranefoot_files
{
	my $file = shift;
	my $out_dir = shift;
	
	open(OUT,">$out_dir/$file.config");
	print OUT "PedigreeFile\t$out_dir/$file.cranefoot.ped\n";
	print OUT "PedigreeName\t$file\n";
	print OUT "SubgraphVariable FID\n";
	print OUT "NameVariable IID\n";
	print OUT "FatherVariable PID\n";
	print OUT "MotherVariable MID\n";
	print OUT "GenderVariable GENDER\n";
	print OUT "PatternVariable AFFECTED_STATUS\n";
	print OUT "TextVariable IID\n";
	close(OUT);
	
	## Write the cranefoot ped file
	open(OUT, ">$out_dir/$file.cranefoot.ped");
	print OUT "FID\tIID\tPID\tMID\tGENDER\tAFFECTED_STATUS\n";
	#foreach my $ID (sort keys %$network_ref)
	foreach my $ID (sort {substr($a,2) <=> substr($b,2)} keys %parents)
	{
		my ($dad,$mom) = @{$parents{$ID} };
		my $FID = "fam1";
		my $node_name = $ID;
			
		my @parent_IDs = @{$parents{$ID} };
		my @parents = @{$parents{$ID} }; ## without the FID
		
		my $g = "0";
		if (exists $samples{$ID})
		{
			$g = $samples{$ID};
			#print "node: $node_name\n";
			#print "gender $g\n";
			if($g eq 1){$g = "M";}
			else{$g = "F";}
		}
		my $a = 1;
		if(@parents eq 0){@parents = (0,0);}
		if(@parents eq 1){push @parents, 0;}
		if(@parents ne 2)
		{
			die "ERROR!!! Incorrect number of parents ". @parents .": @parents\n";
		}
		my $pat;
		my $mat;
		if(exists $samples{@parent_IDs[0]} && $samples{@parent_IDs[0]} eq $samples{@parent_IDs[1]} && $samples{@parent_IDs[0]} ne 0)
		{
			print "WARNING!!! Parents @parents[0] and @parents[1] are the same gender " . $samples{@parents[1]}." !!!\n";
		}
		elsif($samples{@parent_IDs[0]} eq 1)
		{
			$pat = @parents[0];
			$mat = @parents[1];
		}
		elsif($samples{@parent_IDs[0]} eq 0 && $samples{@parent_IDs[1]} eq 2)
		{
			#print "2. $node_name\t$pat\t$mat\n";
			$pat = @parents[0]; 
			$mat = @parents[1]; 
		}
		else
		{
			$pat = @parents[1];
			$mat = @parents[0]; 
		}
		
		## Set the Missing fill
		if($node_name =~ /Missing/i){$a = 11}
		
		my $line = "Fam1\t$node_name\t$pat\t$mat\t$g\t$a";
		print OUT "$line\n";
	}
	close(OUT);

    run_cranefoot($out_dir, $file);
	#print "file $file\n";
	#exit;

}

sub reset_variables
{
	%children = (); ## $children{sample} = [sample's children]
	%parents = (); ## $parents{sample} = [dad,mom]
	%samples = (); ## $samples{sample} = gender
	%generations = (); ## $generations{sample} = generation
	%married_into_family = ();
	@samples_to_visit = ();
	#$mean_children = 3;
	$sample_ctr = 0;
	%pairwise_rels = ();
}

sub run_system
{
	my $command_to_run = shift;
	if ($ECHO_SYSTEM) {
		print("$command_to_run\n");
	}
	system($command_to_run);
	if ($? != 0) {
		my $exit_code = $? >> 8;
		confess "Exit code was $exit_code Failed to run $command_to_run";
	}
}
