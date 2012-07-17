#!/usr/bin/perl -w
use strict;

my($my_cmake_version);
my($ccmake);
my($cmake);
my($ctest);
my($cpack);

$my_cmake_version = `cmake --version`;
$ccmake = `which ccmake`;
$cmake = `which cmake`;
$ctest = `which ctest`;
$cpack = `which cpack`;
print "$ccmake";
print "$cmake";
print "$ctest";
print "$cpack";

chomp $ccmake;
chomp $cmake;
chomp $ctest;
chomp $cpack;

sub compile_cmake {
    my $major = shift;
    my $minor = shift;
    my $rel   = shift;
    system "wget http://www.cmake.org/files/v$major.$minor/cmake-$major.$minor.$rel.tar.gz";
    system "tar zxvf cmake-$major.$minor.$rel.tar.gz";
    chdir("./cmake-$major.$minor.$rel");
    system "cmake .";
    system "make";
    if ( -e "/usr/bin/apt-get" ) {
        system "sudo apt-get -y remove git-core";
    } else {
        if ( -e "/usr/bin/yum" ) {
            system "sudo yum -y remove git-core";
        }
    }
    system "sudo make install";
}

if ( -f "$cmake" ) {
	if ($my_cmake_version =~ /.*2.[8-9].[6-9].*$/) {
		print "cmake version : $my_cmake_version";
	} else {
		if ($my_cmake_version =~ /.*2.6.[1-9].*$/) {
			print "version > 2.6.0\n";
            compile_cmake(2,8,8);
		} else {
			print "version < 2.6.1\n";
            compile_cmake(2,6,4);
			chdir("..");
            compile_cmake(2,8,8);
		}

        $my_cmake_version = `cmake --version`;
        print "cmake version : $my_cmake_version";
        $cmake = `which cmake`;
        $ctest = `which ctest`;
        $cpack = `which cpack`;
        print "$cmake";
        print "$ctest";
        print "$cpack";
	}
}
