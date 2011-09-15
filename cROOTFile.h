#ifndef CROOTFILE_H 
#define  CROOTFILE_H

class TFile;

class cROOTFile {
  public:
    cROOTFile(char*);
    ~cROOTFile();
    void writeFreqTrees(struct FFTinput *,int );
    void writeTimeTrees(struct FFTinput *,int );
    void inputTestData(struct FFTinput *);
    void fitTimeNoise(int );
    void createPS(int);
    void findFreqPeaks(int,double*, int &, const int);
    void fitFreqNoise(int, double*, int );
    void fitFreqPeak(int, double*, int );
    void closeFile();

    TFile *fFile;
};
#endif
