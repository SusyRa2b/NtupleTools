#!/usr/bin/perl
use strict;
use warnings;

#Use to fix the automatic .C output if you use SaveAs("myFile.C").
#Run by doing "perl fixC.pl myFile.C".
#Note: doesn't work if you save a root file as a C file using the TBrowser.
#--Ben

print "Fixing $ARGV[0]...\n";

open(my $fin,  "<", "$ARGV[0]") or die "Or wait, file does not exist.";
open(my $fout, ">", "fixed_$ARGV[0]") or die "Can't create output file.";

while(<$fin>){
    my $line = $_;
    
    if(($ARGV[0] =~ /mindpPassOverFail/) && !($ARGV[0] =~ /Data/)){ #jumbo
	if($line =~ /SetBottomMargin/){
	    print "   - Fixing right margin.\n";
	    print $fout $line;
	    print $fout "thecanvas->SetRightMargin(0.1);\n";
	    next;
	}
    }
    
    if(($ARGV[0] =~ /mindpPassOverFail/) && ($ARGV[0] =~ /Data/)){ #pfdata
	if(($line =~ /SetBottomMargin/) && ($line =~ /thecanvas_1/)){
	    print "   - Fixing right margin.\n";
	    print $fout $line;
	    print $fout "thecanvas_1->SetRightMargin(0.1);\n";
	    next;
	}
	
	if($line =~ /totalsm->SetLineWidth/){
	    print "   - Fixing marker style of SM.\n";
	    print $fout $line;
	    print $fout "totalsm->SetMarkerStyle(1);\n";
	    next
	    }
    }
    
    if($ARGV[0] =~ /METshape/){
	if($line =~ /SIGplot->GetYaxis\(\)->SetTitleSize/){
	    print "   - Fixing Y axis offset.\n";
	    print $fout $line;
	    print $fout "   SIGplot->GetYaxis()->SetTitleOffset(1.);\n";
	    next;
	}
    }
    
    if($ARGV[0] =~ /effAtBestUL/){
	if($line =~ /ROOT/){
	    print "   - Fixing text precision.\n";
	    print "   - Fixing error bar style.\n";
	    print $fout $line;
	    print $fout "   gStyle->SetErrorX(0.);\n";
	    print $fout "   gStyle->SetPaintTextFormat(\"2.0f\");\n";
	    next;
	}
    }
    
    if($ARGV[0] =~ /bestUL_/){
	if($line =~ /ROOT/){
	    print "   - Fixing text precision.\n";
	    print "   - Fixing error bar style.\n";
	    print $fout $line;
	    print $fout "   gStyle->SetErrorX(0.);\n";
	    print $fout "   gStyle->SetPaintTextFormat(\"3.2f\");\n";
	    next;
	}
    }
    
    
    if($line =~ /(.*)(cloned_)(\d)(.)(.*)(cloned_)(\d)(.)(.*)/){
	print $fout "   $1$2$3$5$6$7$9\n";
	next;	
    }
    elsif( !($line =~ /SetHistogram/) && ($line =~ /(.*)(cloned_)(\d)(.)(.*)/)){ 
	print $fout "   $1$2$3$5\n";
	next;
    }
    elsif(($line =~ /SetHistogram/) && ($line =~ /(.*)(cloned_)(\d)(\.)(.*)(\d)(\)\;)/)){
	print $fout "$1$2$3$5$6__$6$7\n";
	next;
    }
    

    if($line =~ /ROOT/){
	print "   - Fixing error bar style.\n";
	print $fout $line;
	print $fout "   gStyle->SetErrorX(0.);\n";
	next;
    }

    if($line =~ /(thestack->SetHistogram\(thestack)(\d*)/){
	print "   - Fixing thestack->SetHistogram().\n";
	print $fout "   thestack->SetHistogram(thestack$2__$2);\n";
	next;
    }
    
    if($line =~ /(TLegend)(\s)(\*)(\S*)(\s)(=)/){
	print "   - Fixing the TLegend: $4\n";
	print $fout $line;
	print $fout "   $4->SetTextFont(42);\n";
	next;
    }
    
    print $fout $line;

}#end loop over fin 
