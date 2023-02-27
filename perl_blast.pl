#!/usr/local/bin/perl 
# 1. Whenever a 4-mer in S is determined to be in Q, extract the location of the first occurrence of that 4-mer in Q.
# 2. Then put the characters of Q and S in arrays 
# 3. Scan left from the k-mer in Q and in S, as long as you find matching characters. Repeat to the right. Let L denote the length of the whole match obtained in this way. 
# 4. If L is greater than 10, then print a message that a good HSP has been found between Q and S, and print S. 
#!/usr/local/bin/perl 

print("k: ");
$k = int(<STDIN>);

print("threshold: ");
$threshold = int(<STDIN>);


print("query filename \n");
$handle = "query.txt";
open(MYFILE, $handle);
$query = <MYFILE>; 

# find location of each different kmer in query
%qkmer = ();                      
$i = 1;
$q = $query;
while (length($q) >= $k) {
  $q =~ m/(.{$k})/; 
  #print "$1, $i \n";
   if (! defined $kmer{$1}) {
    $qkmer{$1} = $i;       
   }
 $i++;
  $q = substr($q, 1, length($q) -1);
}

print("sequence filename \n");
$handle = "sequences.txt";
open(MYFILE, $handle);
while (<MYFILE>) {
	chomp;
	$seq .= $_;
}



# find location of each different kmer in sequence
%skmer = ();                      
$i = 1;
$s = $seq;
while (length($s) >= $k) {
  $s =~ m/(.{$k})/; 
 #print "$1, $i \n";
   if (! defined $kmer{$1}) {
    $skmer{$1} = $i;       
   }
 $i++;
  $s = substr($s, 1, length($s) -1);
}
# get substrings of query
@queries = ();
for $i (0..length($query)-$k){
    push(@queries, substr($query, $i, $k));
}

# get substrings of the sequence
@sequences = ();
for $j (0..length($seq)-$k){
    push(@sequences, substr($seq, $j, $k));
}

# split query / seq into list of characters
@q_char = split("", $query);
@s_char = split("", $seq);

# create database of location of kmers
# key is kmer, value is array of locations
%s_database = {};
for $s_index (0..$#sequences){
	$s = $sequences[$s_index];
	if (exists $s_database{$s}){
		push @{$s_database{$s}}, $s_index;
	}
	else{
		$s_database{$s} = [$s_index];
	}
}
# print everything in %s_database
=begin
while (($key, $value) = each %s_database) {
	print "$key => @$value \n";
}

=cut

my %string_hash = ();
# iterate through kmers in query and check if in database
for my $q_index (0 .. $#queries-$k+1) {
	$l = $queries[$q_index];
	$l_length = 0;
	$r_length = 0;
	@ls_l_lengths = ();
	@ls_r_lengths = ();
	if (exists $s_database{$l}){
		# $len_l_s_database = @{$s_database{$l}};
		# iterate thru all loc of kmer in seq
		# use loc with longest match
		if (scalar @{$s_database{$l}} > 1){ # checks if length of matches is > 1
			# check left
			for $s (@{$s_database{$l}}){
				$s_index = $s;
				# iterate to the left of the kmer in query and seq
				$q_index = $q_index - 1;
				$s_index = $s_index - 1;
				while ($q_index >= 0 and $s_index >= 0){
					if ($q_char[$q_index] eq $s_char[$s_index]){
						$l_length += 1;
						$q_index -= 1;
						$s_index -= 1;}
					else{
						push (@ls_l_lengths, $l_length);
						$l_length = 0;
						$q_index = index($queries, $l);
						last;}
				}
			}
			# check right
			$q_index = index($query, $l) + $k;
			for $s (@{$s_database{$l}}){
				$s_index = $s;
				$s_index = $s_index + $k;
				while ($q_index <= ($#queries+1-$k) and $s_index <= ($#sequences+1-$k)){
					if ($q_char[$q_index] eq $s_char[$s_index]){
						$r_length += 1;
						$q_index += 1;
						$s_index += 1;}
					else{
						push @ls_r_lengths, $r_length;
						$r_length = 0;
						$q_index = index($queries, $l)+ $k;
						last;}
				}
			}		
			# get longest match
			$s_index = @{$s_database{$l}}[0];
			$longest_match = 0;
			$len_l = @ls_l_lengths;
			while (my ($i, $v) = each @ls_l_lengths){
				#print("$l, $i, $v, \n");
				$match = $ls_l_lengths[$i] + $ls_r_lengths[$i];
				if ($match > $longest_match){
					$longest_match = $match;
					$s_index = $s_database{$l}[$i];
					$l_length = $ls_l_lengths[$i];
					$r_length = $ls_r_lengths[$i];}

			}
		}
		else{
			$s_index = @{$s_database{$l}}[0]; # dereference array
			# check left
			while ($q_index >= 0 and $s_index >= 0){
				if ($q_char[$q_index] eq $s_char[$s_index]){
					$l_length += 1;
					$q_index -= 1;
					$s_index -= 1;}
				else{
					last;
				}
			}
			# check right
			# get index of query kmer using List utils

			# dereference index of s from database
			$s_index = @{$s_database{$l}}[0] + $k;
			$q_index = index($query, $l) + $k;
			while ($q_index <= $#queries and $s_index <= $#sequences){ # $# basically the same as length - 1
				if ($q_char[$q_index] eq $s_char[$s_index]){
					$r_length += 1;
					$q_index += 1;
					$s_index += 1;
				}
				else{
					last;
				}
			}
			$s_index = @{$s_database{$l}}[0]; 
		}
	$length = $l_length + $r_length + $k;
	if ($length >= $threshold){
		$output = substr($seq, $s_index-$l_length+1, $length-1);
		if (! defined $string_hash{$output}){
			$string_hash{$output} = 1;
			print("$l, $output, $s_index, $l_length, $r_length \n");
		}
	}	
	}
}
=begin

=cut

=begin
while (my ($key, $value) = each(%s_database)) {
	$x = scalar @{$s_database{$key}};
		print("$key, @$value, $x \n")
	}
=cut