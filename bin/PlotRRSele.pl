#!/usr/bin/perl-w
use strict;

use Data::Dumper;
use Getopt::Long;

#############Befor  Start  , open the files ####################

sub usage
{
	print STDERR <<USAGE;
	2019-01-08       hewm\@genomics.cn

	Usage:	  perl   $0  -inFile  In.winRR.gz  -output OutPrefix

	Options
	-inFile   <s>  :  Input RRSelection OutPut Stat File(path/xxx.winRR.gz)
	-inList   <s>  :  Input FileList if multi-File of RRSelection OutPut Stat
	-output   <s>  :  Output Figure File Prefix

	-Pvalue   <f>  :  the cutoff Pvalue to pick out the seletion Region [0.005]
	-keepR         :  keep the R script to modify and RePlot OUT Fig

	-help          :  show this help

USAGE
}

my ($help,$inFile,$inList,$output,$keepR,$Pvalue);

GetOptions(
	"help"=>\$help,
	"inFile:s"=>\$inFile,
	"inList:s"=>\$inList,
	"output:s"=>\$output,
	"Pvalue:s"=>\$Pvalue,
	"keepR"=>\$keepR,
);

if( defined($help) || !defined($output))
{
	usage;
	exit ;
}
if( (!defined($inFile)) && (!defined($inList)))
{
	usage;
	exit ;
}

$Pvalue||=0.005 ;


################ Do what you want to do #######################


my $Format=0;
my  @header =();


my %hash1=();
my %ChrLenCount=();
my %hash2=();
my $PlotPointbin=1000; ## 100 1000 10000 ;

my %RRNum_1=();
my %RRNum_1P=();
my %RRNum_2=();
my %RRNum_2P=();


###########################   input File  or File List ###################



########################   input File  ###################
if ( defined($inFile) )
{

	if  ($inFile =~s/\.gz$/\.gz/)
	{
		open AA,"gzip -cd  $inFile | "  || die "input file can't open $!" ;
	}
	else
	{
		open AA,"$inFile"  || die "input file can't open $!" ;
	}
	@header=();
	$_=<AA> ; chomp  $_ ; @header=split ;
	if ($#header > 9)
	{
		if ( $Format ==0 )
		{
			$Format = 2 ;
		}
		else
		{
			if (  $Format  !=2 )
			{
				print "Two File format Diff ,something wrong in the file format\n";
				exit ;
			}
		}
		<AA> ;  <AA> ;  <AA> ;
	}
	else
	{
		if ( $Format ==0 )
		{
			$Format = 1 ;
		}
		else
		{
			if (   $Format  != 1 )
			{
				print "Two File format Diff ,something wrong in the file format\n";
				exit ;
			}
		}
		<AA> ;
	}



	if ( $Format  == 1 ) 
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if  ($_=~s/#/#/);
			if (!exists $ChrLenCount{$inf[0]} )
			{
				$ChrLenCount{$inf[0]}=1;
			}
			else
			{
				$ChrLenCount{$inf[0]}++;
			}
			my $CCCC=$ChrLenCount{$inf[0]} ;
			$hash1{$inf[0]}{$CCCC}=$inf[3] ;


			next if  ($_=~s/NA/NA/);
			$inf[3]=int($inf[3]*$PlotPointbin)/$PlotPointbin;
			$RRNum_1{$inf[3]}++;
			if  ( $inf[-1]  <  $Pvalue )
			{
				if  ( $inf[3] <0.2 )
				{
					$RRNum_1P{$inf[3]}=2;
				}
				else
				{
					$RRNum_1P{$inf[3]} =3;
				}
			}
		}
	}
	else
	{

		while(<AA>)
		{
			chomp;
			my @inf=split /\s+/;
			next if  ($_=~s/#/#/);
			if (!exists $ChrLenCount{$inf[0]} )
			{
				$ChrLenCount{$inf[0]}=1;
			}
			else
			{
				$ChrLenCount{$inf[0]}++;
			}
			my $CCCC=$ChrLenCount{$inf[0]} ;
			$hash1{$inf[0]}{$CCCC}=$inf[3] ;
			$hash2{$inf[0]}{$CCCC}=$inf[6] ;

			next if  ($_=~s/NA/NA/);
			$inf[3]=int($inf[3]*$PlotPointbin)/$PlotPointbin;
			$inf[6]=int($inf[6]*$PlotPointbin)/$PlotPointbin;
			$RRNum_2{$inf[3]}{$inf[6]}++;
			if  ( $inf[-1] < $Pvalue)
			{
				if ($inf[3] > $inf[6])
				{
					$RRNum_2P{$inf[3]}{$inf[6]}=2;
				}
				else
				{
					$RRNum_2P{$inf[3]}{$inf[6]}=3;
				}
			}
		}
	}
	close AA;

}






########################   input File  List ###################



if ( defined($inList) )
{

	if  ($inList =~s/\.gz$/\.gz/)
	{
		open LIST,"gzip -cd  $inList | "  || die "input file can't open $!" ;
	}
	else
	{
		open LIST,"$inList"  || die "input file can't open $!" ;
	}

	while($_=<LIST>)
	{
		my $FileThis=$_; chomp $FileThis;

		if  ($FileThis =~s/\.gz$/\.gz/)
		{
			open BBC,"gzip -cd  $FileThis | "  || die "input file can't open $!" ;
		}
		else
		{
			open BBC,"$FileThis"  || die "input file can't open $!" ;
		}



		@header=();
		$_=<BBC> ; chomp  $_ ; @header=split ;
		if ($#header > 9)
		{
			if ( $Format ==0 )
			{
				$Format = 2 ;
			}
			else
			{
				if (  $Format  !=2 )
				{
					print "Two File format Diff ,something wrong in the file format\n";
					print "$FileThis\n";
					exit ;
				}
			}
			<BBC> ;  <BBC> ;  <BBC> ;
		}
		else
		{
			if ( $Format ==0 )
			{
				$Format = 1 ;
			}
			else
			{
				if (   $Format  != 1 )
				{
					print "Two File format Diff ,something wrong in the file format\n";
					print "$FileThis\n";
					exit ;
				}
			}
			<BBC> ;
		}



		if ( $Format  == 1 ) 
		{
			while (<BBC>)
			{
				chomp  ;
				my @inf=split ;
				next if  ($_=~s/#/#/);
				if (!exists $ChrLenCount{$inf[0]} )
				{
					$ChrLenCount{$inf[0]}=1;
				}
				else
				{
					$ChrLenCount{$inf[0]}++;
				}
				my $CCCC=$ChrLenCount{$inf[0]} ;
				$hash1{$inf[0]}{$CCCC}=$inf[3] ;


				next if  ($_=~s/NA/NA/);
				$inf[3]=int($inf[3]*$PlotPointbin)/$PlotPointbin;
				$RRNum_1{$inf[3]}++;
				if  ( $inf[-1]  <  $Pvalue )
				{
					if  ( $inf[3] <0.2 )
					{
						$RRNum_1P{$inf[3]}=2;
					}
					else
					{
						$RRNum_1P{$inf[3]} =3;
					}
				}
			}
		}
		else
		{

			while(<BBC>)
			{
				chomp;
				my @inf=split /\s+/;
				next if  ($_=~s/#/#/);
				if (!exists $ChrLenCount{$inf[0]} )
				{
					$ChrLenCount{$inf[0]}=1;
				}
				else
				{
					$ChrLenCount{$inf[0]}++;
				}
				my $CCCC=$ChrLenCount{$inf[0]} ;
				$hash1{$inf[0]}{$CCCC}=$inf[3] ;
				$hash2{$inf[0]}{$CCCC}=$inf[6] ;

				next if  ($_=~s/NA/NA/);
				$inf[3]=int($inf[3]*$PlotPointbin)/$PlotPointbin;
				$inf[6]=int($inf[6]*$PlotPointbin)/$PlotPointbin;
				$RRNum_2{$inf[3]}{$inf[6]}++;
				if  ( $inf[-1] < $Pvalue)
				{
					if ($inf[3] > $inf[6])
					{
						$RRNum_2P{$inf[3]}{$inf[6]}=2;
					}
					else
					{
						$RRNum_2P{$inf[3]}{$inf[6]}=3;
					}
				}
			}
		}
		close BBC;


	}


}






####################################################################



if  ( $Format ==2 )
{

	open OA,">$output.tmp" || die "output file can't open $!" ;
	open OB,">$output.tmp2" || die "output file can't open $!" ;


	print OA "Site";
	my $Max=0;
	my @hhC=split /\_/,$header[3];
	my @hhW=split /\_/,$header[6];
	foreach my $k(sort  keys  %ChrLenCount)
	{
		print OA "\t$hhC[-1]_$k\t$hhW[-1]_$k";
		if ($ChrLenCount{$k}> $Max)
		{
			$Max=$ChrLenCount{$k};
		}
	}
	print OA "\n";

	foreach my $site (1..$Max)
	{
		print OA  $site;
		foreach my $k(sort  keys  %ChrLenCount)
		{
			my $BBBBB=$hash1{$k}{$site} ;	$BBBBB||="NA"; 
			my $CCCCC=$hash2{$k}{$site} ;	$CCCCC||="NA";

			if  ($BBBBB eq "NA") 
			{
				if ( $site > $ChrLenCount{$k})
				{
					$BBBBB=0;
				}
			}
			if  ($CCCCC eq "NA") 
			{
				if ( $site > $ChrLenCount{$k})
				{
					$CCCCC=0;
				}
			}

			print OA  "\t$BBBBB\t$CCCCC";
		}
		print OA "\n";
	}

	close OA ;

	print OB "MeanRRA\tMeanRRB\tNum\tType\n";
	foreach my $key1  ( sort {$a<=>$b} keys %RRNum_2)
	{
		my $BBB=$RRNum_2{$key1};
		foreach my $key2  (sort {$a<=>$b}keys %$BBB)
		{
			if (!exists $RRNum_2P{$key1}{$key2})
			{
				print OB "$key1\t$key2\t$RRNum_2{$key1}{$key2}\tNoSele\n";
			}
			elsif ( $RRNum_2P{$key1}{$key2} ==2 )
			{
				print OB "$key1\t$key2\t$RRNum_2{$key1}{$key2}\tSele$hhC[-1]\n";
			}
			else
			{
				print OB "$key1\t$key2\t$RRNum_2{$key1}{$key2}\tSele$hhW[-1]\n";
			}
		}
	}
	close OB;



	my $Rshell=<<LOVE;

library(ggplot2);
library(reshape2);
library(gridExtra);
library(ggExtra);

read.table("$output.tmp",header=T)->data
data2<-melt(data,id=1);
p <- ggplot(data = data2)+geom_tile(aes(x = Site,y = variable, fill=value))+scale_fill_gradient(limit=c(0,1),low ='white',high ='red',na.value="grey50")
pdf("$output.sliding.pdf",h=6,w=25)
p
dev.off()
png("$output.sliding.png",h=6,w=25)
p
dev.off()


read.table("$output.tmp2",header=T)->data;
MeanRR$hhC[-1]<-data[,1];
MeanRR$hhW[-1]<-data[,2];
Num<-data[,3];
Type<-data[,4];
df <- data.frame(MeanRR$hhC[-1],MeanRR$hhW[-1]);
p <- ggplot(df, aes(MeanRR$hhC[-1], MeanRR$hhW[-1],colour=Num,shape=Type)) + geom_point(size=1.5) +theme(legend.direction="horizontal",legend.position= "bottom",legend.title = element_text(face ="bold",size=7));
pdf("$output.Dis.pdf");
ggExtra::ggMarginal(p, type = "histogram");
dev.off()
png("$output.Dis.png");
ggExtra::ggMarginal(p, type = "histogram");
dev.off()

LOVE

	open TMP,">$output.r" || die "output file can't open $!" ;

	print TMP $Rshell ;
	close TMP;

}
elsif  ( $Format ==1 )
{


	open OA,">$output.tmp" || die "output file can't open $!" ;
	open OB,">$output.tmp2" || die "output file can't open $!" ;


	print OA "Site";
	my $Max=0;
	my @hhC=split /\_/,$header[3];
	foreach my $k(sort  keys  %ChrLenCount)
	{
		print OA "\t$k";
		if ($ChrLenCount{$k}> $Max)
		{
			$Max=$ChrLenCount{$k};
		}
	}
	print OA "\n";

	foreach my $site (1..$Max)
	{
		print OA  $site;
		foreach my $k(sort  keys  %ChrLenCount)
		{
			my $BBBBB=$hash1{$k}{$site} ;	$BBBBB||="NA"; 

			if  ($BBBBB eq "NA") 
			{
				if ( $site > $ChrLenCount{$k})
				{
					$BBBBB=0;
				}
			}

			print OA  "\t$BBBBB";
		}
		print OA "\n";
	}

	close OA ;

	print OB "MeanRRA\tNum\tType\n";
	foreach my $key1  ( sort {$a<=>$b} keys %RRNum_1)
	{
		my $BBB=$RRNum_1{$key1};
			if (!exists $RRNum_1P{$key1})
			{
				print OB "$key1\t$BBB\tNoSele\n";
			}
			elsif ( $RRNum_1P{$key1} ==2 )
			{
				print OB "$key1\t$BBB\tLow\n";
			}
			else
			{
				print OB "$key1\t$BBB\tHighSele\n";
			}
	}
	close OB;



	my $Rshell=<<LOVE;

library(ggplot2);
library(reshape2);
library(gridExtra);
library(ggExtra);
library(ggthemes);

read.table("$output.tmp",header=T)->data
data2<-melt(data,id=1);
p <- ggplot(data = data2)+geom_tile(aes(x = Site,y = variable, fill=value))+scale_fill_gradient(limit=c(0,1),low ='white',high ='red',na.value="grey50")
pdf("$output.sliding.pdf",h=6,w=25)
p
dev.off()
png("$output.sliding.png",h=6,w=25)
p
dev.off()


read.table("$output.tmp2",header=T)->data;
MeanRR$hhC[-1]<-data[,1];
Num<-data[,2];
sele<-data[,3];
df <- data.frame(MeanRR$hhC[-1],Num);
p <- ggplot(df, aes(x = MeanRR$hhC[-1], y = Num, col = sele,fill=factor(sele) ) ) + geom_bar(stat = "identity") + theme_economist();

pdf("$output.Dis.pdf");
p
dev.off()
png("$output.Dis.png");
p
dev.off()

LOVE

	open TMP,">$output.r" || die "output file can't open $!" ;

	print TMP $Rshell ;
	close TMP;




}



##################   Run R  and plot the result ############

my $R="/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/R-3.3.1/bin/Rscript";
if  ( !(-e $R) )
{
	$R=`which Rscript`;chomp $R;
}


if  ( !(-e $R) )
{
	print "Can't find the [ Rscript ] bin in your [\$PATH] path, You shoud install the R First,then:\n";
	print " $R  $output.r  \n";
	exit(1);
}



system (" $R  $output.r  ");

if  (  defined($keepR) )
{
	system ("echo  $R  $output.r  ");
}
else
{
	system ("rm -rf   $output.r  $output.tmp $output.tmp2 ");
}


######################swiming in the sky and flying in the sea #############################


