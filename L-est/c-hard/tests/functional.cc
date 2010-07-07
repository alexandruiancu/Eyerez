/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "functional.hpp"

static MetropolisSetup suite_MetropolisSetup;

static CxxTest::List Tests_MetropolisSetup = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_MetropolisSetup( "functional.hpp", 6, "MetropolisSetup", suite_MetropolisSetup, Tests_MetropolisSetup );

static class TestDescription_MetropolisSetup_testMakeMetropolis : public CxxTest::RealTestDescription {
public:
 TestDescription_MetropolisSetup_testMakeMetropolis() : CxxTest::RealTestDescription( Tests_MetropolisSetup, suiteDescription_MetropolisSetup, 9, "testMakeMetropolis" ) {}
 void runTest() { suite_MetropolisSetup.testMakeMetropolis(); }
} testDescription_MetropolisSetup_testMakeMetropolis;

static MetropolisFunction suite_MetropolisFunction;

static CxxTest::List Tests_MetropolisFunction = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_MetropolisFunction( "functional.hpp", 29, "MetropolisFunction", suite_MetropolisFunction, Tests_MetropolisFunction );

static class TestDescription_MetropolisFunction_testDidJump : public CxxTest::RealTestDescription {
public:
 TestDescription_MetropolisFunction_testDidJump() : CxxTest::RealTestDescription( Tests_MetropolisFunction, suiteDescription_MetropolisFunction, 52, "testDidJump" ) {}
 void runTest() { suite_MetropolisFunction.testDidJump(); }
} testDescription_MetropolisFunction_testDidJump;

#include <cxxtest/Root.cpp>
