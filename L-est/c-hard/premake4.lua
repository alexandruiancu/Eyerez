dofile "p4defs.lua"

newoption {
   trigger     = "valgrind",
   description = "Check tests using Valgrind."
}

solution "MySolution"
	configurations {'Debug', 'Release'}

        includedirs { 'vendor/include', 'src/' }
        libdirs { 'vendor/lib' }
        links { 'csv_parser', 'gsl', 'gslcblas' }

        configuration 'Debug'
        	targetdir 'build/debug'
        	defines { "_DEBUG" }
                flags { 
                   'ExtraWarnings', 
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
           "src/classes/*.cc"
        }

test('foo')

if _ACTION == 'clean' then
   for k, name in pairs(tests) do
      os.execute ('rm -rf build')
      os.remove ('tests/' .. name .. '.cc')
   end
end