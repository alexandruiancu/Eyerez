alltests = {}

function tests( tests )
   for k, name in pairs(tests) do
      test (name)
   end
end

function test(test_name)
   table.insert(alltests, test_name)

   project ('test_' .. test_name)
      kind 'ConsoleApp'
      language 'C++'

      prebuildcommands {          
         ('@cd tests; python ../cxxtestgen.py ' ..
          '--error-printer ' ..
          '-o ' .. test_name .. '.cc ' ..
          test_name .. '.hpp')
      }

      files { 
         ('tests/' .. test_name .. '.hpp'), 
         ('tests/' .. test_name .. '.cc')
      }
      flags { 'ExtraWarnings', 'FatalWarnings', 'Symbols' }

      includedirs { 'vendor' }
      libdirs { 'vendor/cxxtest' }
      links { 'lmcmc', 'csv_parser' }
      buildoptions { "-pg" }

      targetdir 'build/tests'
      postbuildcommands { 
         ('@echo "----> Running test (' .. test_name .. ')"'),
         "@build/tests/test_" .. test_name 
      }

      configuration 'valgrind'
         postbuildcommands {
            ('@mkdir -p build/tests/log'),
            ('@echo "====> Valgrinding test (' .. test_name .. ')"'),
            ('@valgrind --leak-check=full ' .. 
             '--log-file=build/tests/log/' .. test_name .. '_valgrind.log ' ..
             'build/tests/test_' .. test_name),
            ('@egrep "(ERROR SUMMARY|lost)" build/tests/log/' .. test_name .. '_valgrind.log')
         }
end