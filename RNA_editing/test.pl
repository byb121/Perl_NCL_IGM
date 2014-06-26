#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
#use LWP::Simple;
#use Getopt::Long;

my $string = "asdfhalksdjfhlqkjwehflkajsdhflkjahwe;kfaposidufypqwoieualsdkjf;alwejfpoasidfj;alkj;lkhoiuadsghqwkjeghalskdjghklasdhgklajshdgklahsdlkhfaklsjhdf";
print length($string)."\n";
print substr($string, 3,1)."\n";
 print substr($string, 0, 80, '')."\n" while (length($string));
exit;
