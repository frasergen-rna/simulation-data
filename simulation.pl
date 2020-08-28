#! /usr/bin/perl

my $hseq;	
my %h;	
my %hc;
my %hl;
my %hs;	
my %base2;
open FILE,"$ARGV[0]";
my %hr;
$/=">";
while(<FILE>){
        chomp;
        while(<FILE>){
                chomp;
                my ($qid,$seq)=split/\n/,$_,2;
                $seq=~s/\s+//g;
                @sca=split/\s/,$qid;
                $hr{$sca[0]}=$seq;
        }
}
$/="\n";
print "[loading] genome.fa\n";
open INN,"$ARGV[1]";
while(<INN>){
	chomp $_;
	my($chr,$site,$a,$col,$ler,$b,$c,$d,$e)=split/\t/,$_;
	$h{$chr}.="$site;";
	$hc{"$chr:$site"}=$col;
	$hl{"$chr:$site"}=$ler;
}
print "[loading] homozygote info\n";
open IN,"$ARGV[2]";
open $o1,">R1.fa";
open $o2,">R2.fa";
open $b1,">r1.tmp";
open $b2,">r2.tmp";
open $l,">log";
while(<IN>){
	chomp $_;
	my($chr,$gene,$s1,$e1,$trans,$s2,$e2,$strand,$number,$exon,$exon_len)=split/\t/,$_;
	@site=split/;/,$h{$chr};
	foreach $sit(@site){
		if($sit>=$s2 && $sit<=$e2){
			$hs{$trans}.="$sit;";
		}
	}
	$deep=sprintf("%01d",$exon_len*100/250);
	if($exon_len>=125){
		for($i=1;$i<=$deep;$i++){
			print $l "$trans:DEEP"."$i\t$exon\t$exon_len\n";
			@jun=split/,/,$exon;
			pop @jun;shift @jun;
			if($exon_len>=300){
				my $range=$exon_len-300+1;
				my $random=int(rand($range));
				print $l "[random] $random\n";
				$r1s=$random;				
				$r2s=$r1s+125+50;
			}else{
				$r1s=0;
				$r2s=$exon_len-125;
			}
			print $l "[r1] start $r1s\n";
			print $l "[r2] start $r2s\n";
			my $parent_random=int(rand 2)+1;
			print $l "[parent_random] $parent_random\n";
			print $l "[homozygote-site] ";
			if($parent_random==1){
				$type="COL";
			}else{
				$type="LER";
			}
			print $o1 ">".$trans.":DEEP".$i.":$type 1 \n";
			print $o2 ">".$trans.":DEEP".$i.":$type 2 \n";
			print $b1 ">".$trans.":DEEP".$i.":$type 1 \n";
			print $b2 ">".$trans.":DEEP".$i.":$type 2 \n";
			for($j=0;$j<=124;$j++){
				$site1=$s2+$r1s+$j;
				$site2=$s2+$r2s+$j;
				foreach $junction(@jun){
					my($left,$right)=split/;/,$junction;
					if($site1>$left){
						$site1=$site1+$right-$left-1;	
					}
					if($site2>$left){
						$site2=$site2+$right-$left-1;	
				}
				print $b1 "$site1-";
				print $b2 "$site2-";
				if($hs{$trans}=~/$site1/){
					print $l "$chr:$site1 ";
					if($parent_random==1){
						$base1=$hc{"$chr:$site1"};
						print $b1 "col-";
					}else{
						$base1=$hl{"$chr:$site1"};
						print $b1 "ler-";
					}
				}else{
					my $base_random=int(rand 100)+1;
					if($base_random==100){
						$base_raw=substr($hr{$chr},$site1-1,1);
						$allbase="ATCG";
						$allbase=~s/$base_raw//g;
						@base=split//,$allbase;
						$base_count=@base;
						my $base_new=int(rand $base_count);
						$base1=@base[$base_new];
						print $b1 "error-";
					}else{
						$base1=substr($hr{$chr},$site1-1,1);
						print $b1 "raw-";
					}
				}
				print $o1 "$base1";
				print $b1 "$base1,";
				if($hs{$trans}=~/$site2/){
					if($parent_random==1){
						$base2=$hc{"$chr:$site2"};
						print $b2 "col-";
					}else{
						$base2=$hl{"$chr:$site2"};
						print $b2 "ler-";
					}
				}else{
					my $base_random=int(rand 100)+1;
					if($base_random==100){
						$base_raw=substr($hr{$chr},$site2-1,1);
						$allbase="ATCG";
						$allbase=~s/$base_raw//g;
						@base=split//,$allbase;
						$base_count=@base;
						my $base_new=int(rand $base_count);
						$base2=@base[$base_new];
						print $b2 "error-";
					}else{
						$base2=substr($hr{$chr},$site2-1,1);
						print $b2 "raw-";
					}
				}
				print $b2 "$base2,";
				$base2{"$trans.$i"}.=$base2;
			}	
			print $o1 "\n";
			$read2 = reverse $base2{"$trans.$i"};
			$read2 =~ tr/GATC/CTAG/;
			print $o2 "$read2\n";
			print $b1 "\n";
			print $b2 "\n";
			print $l "\n";
		}
	}else{
	}
	
}

