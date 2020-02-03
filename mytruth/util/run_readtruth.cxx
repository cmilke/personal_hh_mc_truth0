#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "mytruth/readtruth_azheavyh.h"
#include "mytruth/readtruth_hhbbzzllvv.h"
#include "mytruth/readtruth_hhbbzz4l.h"
#include "mytruth/readtruth_hhbbwwlvlv.h"
#include "mytruth/readtruth_hhbbwwlvlv_chkweight.h"
#include "mytruth/readtruth_hhbbttlvlv.h"
#include "mytruth/readtruth_hh4b.h"
#include "mytruth/readtruth_vbfhh.h"

int main( int argc, char* argv[] ){

  // Take the submit directory from the input if provided:
  if( argc < 5 ){
    std::cout << "run_readtruth submitDir inputDir inputPattern channel [number_of_events] [others]" << std::endl;
	exit(1);
  }
  std::string submitDir = argv[ 1 ];
  std::string inputDir = argv[ 2 ];
  std::string inputPattern = argv[ 3 ];
  std::string channel = argv[4];

  double nEvtMax = -1;
  if( argc == 6 ){
    nEvtMax = atof(argv[5]);
  }

  std::string option7 = "";
  if( argc == 7 ){
    option7 = argv[6];
  }

  // channel specific

  // azheavyh
  bool isbbA = false;
  int mH = 300; // only used when isbbA true
  if( argc == 8 ){
    isbbA = (std::string( argv[6] ) == "bbA");
    mH = atoi( argv[7] );
  }

  // Set up the job for xAOD access
  xAOD::Init().ignore();
  gErrorIgnoreLevel = kError;

  // Construct the samples to run on
  SH::SampleHandler sh;

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  const char* inputFilePath = gSystem->ExpandPathName( inputDir.c_str() );
  SH::ScanDir().filePattern(inputPattern).scan(sh,inputFilePath);

  // set the name of tree
  sh.setMetaString("nc_tree","CollectionTree");

  // sample information
  sh.print();

  // basic description of job
  EL::Job job;
  job.sampleHandler( sh );
  if( nEvtMax > 0 ){
    job.options()->setDouble( EL::Job::optMaxEvents, nEvtMax );
  }

  // add algo
  if( channel == "azheavyh" ){
    readtruth_azheavyh* rt = new readtruth_azheavyh( isbbA, mH );
    job.algsAdd( rt );
  }else if( channel == "hhbbzzllvv" ){
    readtruth_hhbbzzllvv* rt = new readtruth_hhbbzzllvv();
	rt->setOption( option7 );
    job.algsAdd( rt );
  }
  else if( channel == "hhbbwwlvlv" ){
    readtruth_hhbbwwlvlv* rt = new readtruth_hhbbwwlvlv();
	rt->setOption( option7 );
	job.algsAdd( rt );
  }
  else if( channel == "hhbbwwlvlv_chkweight" ){
    readtruth_hhbbwwlvlv_chkweight* rt = new readtruth_hhbbwwlvlv_chkweight();
	job.algsAdd( rt );
  }
  else if( channel == "hhbbttlvlv" ){
    readtruth_hhbbttlvlv* rt = new readtruth_hhbbttlvlv();
	rt->setOption( option7 );
	job.algsAdd( rt );
  }
  else if( channel == "hhbbzz4l" ){
    readtruth_hhbbzz4l* rt = new readtruth_hhbbzz4l();
    job.algsAdd( rt );
  }
  else if( channel == "vbfhh" ){
    readtruth_vbfhh* rt = new readtruth_vbfhh();
    job.algsAdd( rt );
  }
  else if( channel == "hh4b" ){
    readtruth_hh4b* rt = new readtruth_hh4b();
    job.algsAdd( rt );
  }
  else{
    Error("main()", "Channel %s not recognised", channel.c_str() );
  }

  // make the driver: local Grid PROOF etc.
  EL::DirectDriver driver;

  // process job via driver
  driver.submit( job, submitDir );

  return 0;
}


