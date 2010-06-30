
solution "MySolution"
	configurations {'Debug', 'Release'}

        includedirs { '.', 'vendor/include' }
        libdirs { '.', 'vendor/lib' }

        configuration 'Debug'
        	targetdir 'build/debug'
        	defines { "_DEBUG" }
                flags { 
                   'ExtraWarnings', 
                   'FatalWarnings',
                   'Symbols'
                }

        configuration 'Release'
                targetdir 'build/release'
                flags {
                   'FloatStrict',
                   'OptimizeSpeed'
                }

project 'lmcmc'
	kind 'SharedLib'
	language 'C++'

        files {
           "classes/*.cc"
        }

        links { 'csv_parser', 'gsl', 'gslcblas' }

project 'test'
        kind 'ConsoleApp'
        language 'C++'
        files 'test/test.cc'

	targetname 'test_me'
	flags { "ExtraWarnings", "FatalWarnings", "Symbols" }
	links { 
           'lmcmc',
           'gsl'
        }
        targetdir 'build/test'