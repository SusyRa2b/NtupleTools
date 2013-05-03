#!/usr/bin/perl -w
use strict;

######
#This is a "prepare_run_reduced.C" file for the T3
#The only thing that needs user modification is the inputdir
#Then do: > ./prepare_filelist.pl
######

# Declare the subroutines
sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

######################################################
#change this!
my $inputdir = 'joshmt/cfA/v66/';
#my $inputdir = 'joshmt/cfA/v66_noskim/';
######################################################


#no need to change these
my $prefix = 'root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/'; 
my $CUbase = 'srm://osg-se.cac.cornell.edu:8443/srm/v2/server?SFN=/xrootdfs/cms/store/user/';


my $directory = '"' . $CUbase . $inputdir . '"';

#my $command = "srmls $directory | grep T1tttt";
my $command = "srmls $directory";


print $command;
print "\n";
#my @files = `$command`;
my @samples = `$command`;

my @samplenames;

#parse the subdirectories, each one corresponding to a sample
foreach(@samples) {
    print $_;
    #my $len = length($_);
    my $pos_lastslash = rindex($_, "/");
    my $result = rindex($_,"/", $pos_lastslash-1);
    #print "string length = $len, pos_lastslash = $pos_lastslash , result = $result\n";
    my $samplename = trim(substr($_, $result+1));
    $samplename =~ s/\/+$//;
    #print "sample name = $samplename\n";
    if( length($samplename) gt 1){
	push( @samplenames, $samplename );
    }
}


print "Samples found : @samplenames\n";


#my $outfilename = "filelist.txt";
#open(my $out, ">",  $outfilename) or die "Can't open filelist.txt: $!";

my $sampledirectory = $CUbase . $inputdir ;
#create one output file list per sample
foreach(@samplenames) {

    my $samplecommand = "srmls " . '"' .  "$sampledirectory$_/" . '"';
    print "$samplecommand\n";

    my @files = `$samplecommand`;

    my $iname = $_;

    #create file list
    my $outfilename = $iname . '.txt';
    open(my $out, ">",  $outfilename) or die "Can't open txt file: $!";

    foreach(@files) {
	#print "$_";
	my $filename = substr($_,rindex($_, "/")+1) ;
	if ($filename ne "\n" ){ #remove blank lines
	    #add the necessary prepending
	    $filename = $prefix . $inputdir . $iname . '/' . $filename; 
	    #print $filename;
	    print $out $filename;
	}
    }
    close $out;
}

print "All done.\n";
