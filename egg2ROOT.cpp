#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
extern "C" { 
  #include "monarch.h"
  #include "cFFTransform.h"
}
#include "cROOTFile.h"
using namespace std;

int main(int argc, char *argv[]) {

  int c;
  char *fileName = "sample.egg";
  char *outputFilePrefix;
  int num_events = 0;

  while((c = getopt(argc, argv, "f:p:n:")) != -1)
    switch(c)
    {
    case 'f':
      fileName = optarg;
      break;
    case 'p':
      outputFilePrefix = optarg;
      break;
    case 'n':
      num_events = atoi(optarg);
      break;
    }

  printf("%s\t%s\t%d\n",fileName,outputFilePrefix,num_events);

  struct egg *current = (egg*)malloc(sizeof(struct egg));
  //open file, read header curtesy of monarch
  mBreakEgg(fileName,current);
  mParseEggHeader(current);
  printf("%s\t%d %s\n","Sample Length :",current->data->sample_length," sec");
  printf("%s\t%d %s\n","Sample Rate:",current->data->sample_rate," MHz");
  printf("%s\t%f %s\n","Actual length of time series:",(float)(current->data->record_size/current->data->sample_rate)," usec");
  
  //initpxlize FFT structure curtesy of cFFTransform
  struct FFTinput *input = (FFTinput*)malloc(sizeof(struct FFTinput));
  setUp(current, input);//set up from header values
  createPlan(input);//create plan before initpxlizing input (b/c overwritten)
  
  //initpxlize Root file  
  cROOTFile rFile(fileName);

  //loop through events in the file
  int i;
  int flag;
  int count = 0;
  int nPeaks;
  const int maxP = 15;
  double px[maxP];
    //"hatch" or read first event
  if (num_events == 0) {
    //default is process all events
    while((flag = mHatchNextEvent(current)) != 1) {
      count++;
      inputData(current, input);//store record into array for fft
      //rFile.inputTestData(input);//put fake data into array for fft
      executePlan(input);//execute fft
      rFile.writeTimeTrees(input, count);
      rFile.fitTimeNoise(count);
      rFile.writeFreqTrees(input, count);
      rFile.createPS(count);
      rFile.findFreqPeaks(count, px, nPeaks, maxP);
      rFile.fitFreqNoise(count, px, nPeaks);
      rFile.fitFreqPeak(count, px, nPeaks);
    } 
  } else {
    //otherwise only process events specified 
    for(i = 0; i < num_events; i++) {
      if ((flag = mHatchNextEvent(current)) != 1) {
        count++;
        inputData(current, input);//store record into array for fft
        executePlan(input);//execute fft
        rFile.writeTimeTrees(input, i);
        rFile.fitTimeNoise(i);
        rFile.writeFreqTrees(input, i);
        rFile.createPS(i);
        rFile.findFreqPeaks(i, px, nPeaks, maxP);
        rFile.fitFreqNoise(i, px, nPeaks);
        rFile.fitFreqPeak(i, px, nPeaks);
      }
    }
  }
  printf("%s%d%s\n","Processed ", count, "events");
  mCleanUp(current);
  cleanUp(input);
  rFile.closeFile();
  return 0;
}
