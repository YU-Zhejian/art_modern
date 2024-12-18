#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20120309
=cut
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

my $SAMTOOLSBIN="samtools";
#$SAMTOOLSBIN="/ifs1/ST_ASMB/USER/yuanjy/huxuesong/tmp/group/rev/test/samtools";
my $MAXREADStoCHECK=10000;
my $MAXINDELEN=3;

die "Usage: $0 <single_sam_bam_file> <output_prefix> [max_running_minutes]\n" if @ARGV<2;
my $name=shift;
my $out=shift;
my $timelimit=shift;
my ($READLEN,$lines)=(0,0);
if ($name =~ /\.bam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam\.gz$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.indel\.gz$/) {
	open IN,'-|',"zcat $name" or die "Error opening $name : $!\n";
	$READLEN=100;
	goto START;
} else {
	die "[x]Unsupport file type.";
}

while (<IN>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read1=split /\t/;
	next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
	next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
	next unless $read1[5] =~ /^(\d+)M$/;
	$READLEN = $1 if $READLEN < $1;
	++$lines;
	last if $lines > $MAXREADStoCHECK;
}
close IN;
$lines=0;
if ($name =~ /\.bam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam\.gz$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} else {
	die "[x]Unsupport file type.";
}

START:
my $start_time = [gettimeofday];
my (%Cnt,%InDel,%LenInDel,%DistInDelMatrix,%DistAll,%RL);
my ($Read12,$PosShift,$cigar,$DelShift);
while (<IN>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read1=split /\t/;
	#next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
	#next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
	next if $read1[5] =~ /S/;	# only IDM, no S.
	my $OPT = join("\t",@read1[11 .. $#read1]);
	next if $OPT =~ /\bXT:A:R\b/;
	if ($read1[1] & 64) {
		$Read12 = 1;
	} elsif ($read1[1] & 128) {
		$Read12 = 2;
	} else {
		warn "[w]",join("\t",@read1),"\n";
		next;
	}
	$Cnt{$Read12}{'All'} += $READLEN;
	++$InDel{$Read12}{'All'};
	if ($read1[1] & 16) {	#  | r  | 0x0010 | strand of the query (1 for reverse)   |
		$PosShift = -$READLEN - 1;
		$DelShift = 1;
		#next;
	} else {	# +
		$PosShift = 0;
		$DelShift = -1;
	}
	$cigar=$read1[5];
	if ($read1[5] =~ /(\d+)[ID]/) {	# I/D in reads
		if ($read1[5] =~ /(\d+)I/) {
			++$InDel{$Read12}{'Ins'};
		}
		if ($read1[5] =~ /(\d+)D/) {
			++$InDel{$Read12}{'Del'};
		}
      # http://davetang.org/muse/2011/01/28/perl-and-sam/
      my $position = '1';
	  my $flag=0;
      while ($cigar !~ /^$/){
         if ($cigar =~ /^([0-9]+[MIDS])/){
            my $cigar_part = $1;
            if ($cigar_part =~ /(\d+)M/){
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)I/){
				if ($1 <= $MAXINDELEN) {
					$flag=1;
					$Cnt{$Read12}{'Ins'} += $1;
					my $p=abs($position+$PosShift);
					if ($DelShift==1) {	# minus strand
						$p -= $1;
					} else {	# + strand
						$p -= 1;
					}
					$DistInDelMatrix{$1}{$p}{$Read12}++;
					++$LenInDel{$1}{$Read12};
#if ($cigar eq "3I2M2D94M" or $p == 98 and $1> 1) {
#	warn "$position -> ",$position+$PosShift,"\t$1,$p,$Read12\t$cigar\t$cigar_part\n$_\n";
#}
					#for my $p ($position .. ($position + $1 -1)) {
						#$DistIns{abs($p+$PosShift)}{$Read12}++;
#warn "$position -> ",$position+$PosShift,"\t$1\t$cigar\t$cigar_part\n$_\n" if abs($position+$PosShift)<=1 or abs($position+$PosShift)>=$READLEN;
					#}
#warn "$position -> ",$position+$PosShift,"\t$cigar\t$cigar_part\n$_\n" if abs($position+$PosShift)<1 or abs($position+$PosShift)>$READLEN;
				}
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)D/){
				if ($1 <= $MAXINDELEN) {
					$flag=1;
					$Cnt{$Read12}{'Del'} += $1;
					my $p=abs($position+$DelShift+$PosShift);
					$DistInDelMatrix{-$1}{$p}{$Read12}++;	# 99M1D1M: D@99, not 100.
					++$LenInDel{-$1}{$Read12};
					#$DistInDelMatrix{-1}{$Read12}++;
#warn "$position -> ",$p,"\t$1\t$cigar\t$cigar_part\n$_\n" if $p<1 or $p>=$READLEN;
				}
               #$position += $1;
            } elsif ($cigar_part =~ /(\d+)S/){
               die "[!]Not ready for this!\n";
               #my $insertion = 'x' x $1;
               #substr($new_ref,$position,0,$insertion);
               #$position += $1;
            }
            $cigar =~ s/$cigar_part//;
         } else {
            die "Unexpected cigar: $cigar\n";
         }
      }
	  if ($flag == 1) {
		  #++$DistAll{abs($_+$PosShift)}{$Read12} for (1 .. $position);
		  #++$RL{$position};
		  die "[x]Read Length Error !\n" if $position != $READLEN + 1;
	  }
	}
	#++$lines;
	#last if $lines > 100*$MAXREADStoCHECK;
	if (defined $timelimit) {
		last if tv_interval( $start_time, [gettimeofday] ) > $timelimit * 60;
	}
}
close IN;

sub getValue($) {
	if (defined $_[0]) {
		return $_[0];
	} else {
		return '0';
	}
}
sub getRatio($$) {
	if (defined $_[0]) {
		if ((defined $_[1]) && ($_[1] != 0)) {
			return $_[0]/$_[1];
		} else {
			return '0';
		}
	} else {
		return '0';
	}
}

open O,'>',$out.'.InDel.matrix' or die "Error opening ${out}.InDel.matrix : $!\n";
print O "[Info]
File = $name
Read_Length = $READLEN
Read_1_Count = $InDel{1}{'All'}
Read_2_Count = $InDel{2}{'All'}
MaxInDel_Length = $MAXINDELEN
<<END
\n[Overall]
Read\tType\tBaseCount\tBaseRatio\tReadCnt\tReadCntRatio\n";
for my $Read12 (sort keys %Cnt) {
	for (sort {$b cmp $a} keys %{$Cnt{$Read12}}) {
		#next if $_ eq 'All';
		print O join("\t",$Read12,$_,$Cnt{$Read12}{$_},$Cnt{$Read12}{$_}/$Cnt{$Read12}{'All'},$InDel{$Read12}{$_},$InDel{$Read12}{$_}/$InDel{$Read12}{'All'}),"\n";
	}
}
print O "\nIndel\tCount1\tCount2\tBaseRatio1\tBaseRatio2\n";
for (sort {$a<=>$b} keys %LenInDel) {
	print O join("\t",$_,getValue($LenInDel{$_}{1}),getValue($LenInDel{$_}{2}),getRatio($LenInDel{$_}{1},$Cnt{1}{'All'})*abs($_),getRatio($LenInDel{$_}{2},$Cnt{2}{'All'})*abs($_)),"\n";
}

my @Ins = (-$MAXINDELEN .. -1, 1 .. $MAXINDELEN);
print O "<<END\n\n[InDel]\nCycle\t",join("\t",@Ins),"\n";

open D,'>',$out.'.InDel.dat' or die "Error opening ${out}.InDel.dat : $!\n";
print D join("\t",'#Cycle',@Ins),"\n";
for my $cyc (1 .. $READLEN) {
	my (@Counts,@Ratios);
	for my $ins (@Ins) {
		push @Counts,getValue($DistInDelMatrix{$ins}{$cyc}{1});
		push @Ratios,getRatio($DistInDelMatrix{$ins}{$cyc}{1},$InDel{1}{'All'})*abs($ins);
	}
	print O join("\t",$cyc,@Counts),"\n";
	print D join("\t",$cyc,@Ratios),"\n";
}

for my $cyc (1 .. $READLEN) {
	my (@Counts,@Ratios);
	for my $ins (@Ins) {
		push @Counts,getValue($DistInDelMatrix{$ins}{$cyc}{2});
		push @Ratios,getRatio($DistInDelMatrix{$ins}{$cyc}{2},$InDel{1}{'All'})*abs($ins);
	}
	print O join("\t",$cyc+$READLEN,@Counts),"\n";
	print D join("\t",$cyc+$READLEN,@Ratios),"\n";
}

close D;
print O "<<END\n";
close O;
