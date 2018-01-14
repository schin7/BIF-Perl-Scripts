#!/usr/bin/perl

use strict;
use warnings;
use LWP::Simple;
#use Mail::Sendmail;
use CGI ":standard";              # Include CGI functions
#use CGI::Carp "fatalsToBrowser"; # Send error messages to browser

# Using perl CGI module and reading in parameters
my $cgi = CGI->new;
#print "Content-Type: text/html\n\n";

my @attributes = $cgi->param('attributes');
my $checkall = $cgi->param('CHECKALL');
my $viruses = $cgi->param('viruses');
my $mailto = $cgi->param('mailto');

# HTML Webpage to display results
print $cgi->header( );  # print the HTML header
print $cgi->start_html( );
print "<html>\n";
print "<head>\n";
print "<title>Genbank Results...</title>\n";
print "</head>\n";
print "<body>\n";
print "<pre>\n";

# Subroutine to count bases

sub countBases($){
   my $seq = shift;
   my $char;
   my $countA = 0;
   my $countT = 0;
   my $countC = 0;
   my $countG = 0;
   my $inval = 0;
 
   $seq =~ s/\bORIGIN\b//;
   $seq =~ s/[0-9]//g;
   $seq =~ s/\s+//g;

   for(my $i = 0; $i < length($seq); $i++){ 

      $char = uc(substr($seq, $i, 1));
      #print "$char\n"; #display character

      if($char eq "A"){
         $countA = $countA + 1;
      }elsif($char eq "T"){
         $countT = $countT + 1;
      }elsif($char eq "C"){
         $countC = $countC + 1;
      }elsif($char eq "G"){
         $countG = $countG + 1;
      }else{
         $inval = $inval + 1;
         print "Invalid Character: $char\n";
      }
   }

   print "BASE COUNT = $countA A, $countT T, $countC C, $countG G, $inval Invalid\n";

}
# Subroutine to enter subdirectory of URL
sub getURL {
   my $virus = $_[1];
   my $baseURL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/";
   my $virusURL = $baseURL . $viruses; 
}

my @tempArray = split('/', $viruses);  # capture the accession number from the string
my $ncbiFile = $tempArray[1];          # the 2nd element of the array after the split '/'
my $rawData;
my $ncbiURL = getURL ($viruses);

# Get URL using LWP::Simple
unless(-e $ncbiFile) {
   $rawData = get($ncbiURL); # Download genbank file and store it in current working directory
   open(FD, "> $ncbiFile") || die("Unable to get page... $ncbiFile $!\n");
   print FD $rawData;
   close(FD);
}

$/ = undef;
open(FD, "< $ncbiFile") || die("Unable to get page... $ncbiFile $!\n");
$rawData = <FD>;
close(FD);

# Use regular expression to process, extract and create webpage
my $storedResult = "";
my $i = 1;
my $tempAttr;
my @lines =();
my $baseRef;

foreach $tempAttr (@attributes) {
   if($tempAttr =~ /LOCUS/) {
      $rawData =~ /(LOCUS.*)DEFINITION/s;
      print "$1";
      $storedResult .= $1;  # $1 stores last regex match, to allows for the data to be sent by mail
   }
   elsif($tempAttr =~ /DEFINITION/) {
      $rawData =~ /(DEFINITION.*?)ACCESSION/s;
      print "$1";
      $storedResult .= $1;  
   }
   elsif($tempAttr =~ /ACCESSION/) {
      $rawData =~ /(ACCESSION.*?)VERSION/s;
      print "$1";
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /VERSION/) {
      $rawData =~ /(VERSION.*?)KEYWORDS/s;
      print "$1";
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /KEYWORDS/) {
      $rawData =~ /(KEYWORDS.*?)SOURCE/s;
      print "$1";
      $storedResult .= $1;  
   }
   elsif($tempAttr =~ /SOURCE/) {
      $rawData =~ /(SOURCE.*?)ORGANISM/s;
      print "$1";
      $storedResult .= $1;
   }
   elsif($tempAttr =~ /ORGANISM/) {
      $rawData =~ /(ORGANISM.*?)REFERENCE/s;
      print "$1";
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /REFERENCE/) {
      @lines = $rawData =~ /(REFERENCE.*?)(?=AUTHORS|CONSRTM)/gs;
      print "$1";
      print @lines;
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /AUTHORS/) {
      @lines = $rawData =~ /(AUTHORS.*?)TITLE/gs;
      print "$1";
      print @lines;
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /TITLE/) {
      @lines = $rawData =~ /(TITLE.*?)JOURNAL/gs;
      print "$1";
      print @lines;
      $storedResult .= $1;
   }
   elsif($tempAttr =~ /JOURNAL/) { 
      @lines = $rawData =~ /(JOURNAL.*?)(?=PUBMED|REMARK|REFERENCE|COMMENT)/gs;
      print "$1";
      print @lines;
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /MEDLINE/) { 
      @lines = $rawData =~ /(PUBMED.*?)REFERENCE/gs;
      print "$1";
      print @lines;
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /FEATURES/) {
      $rawData =~ /(FEATURES.*?)3'UTR/s;
      print "$1";
      $storedResult .= $1; 
   }
   elsif($tempAttr =~ /BASECOUNT/) { #Count number of nucleotides in sequence
      $rawData =~ /(ORIGIN.*?)\/\//gs;
      $baseRef = $1;
      countBases($baseRef);
      #print "$1";
      $storedResult .= $1;  
   }
   elsif($tempAttr =~ /ORIGIN/) {
      $rawData =~ /(ORIGIN.*$)/s;
      print "$1";
      $storedResult .= $1; 
   }
}

my %mail = ( 
      From    => 'schin8@myseneca.ca',
      To      => '$mailto',
      Subject => 'Genbank File',
      Message => "$storedResult"
      );
sendmail(%mail) or die $Mail::Sendmail::error;


print "</pre>";
print "</body>";
print "</html>\n";
print $cgi->end_html( );

