// This script is part of HPG-Hmapper
//
// HPG-Hmapper is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// HPG-Hmapper is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HPG-Hmapper. If not, see <http://www.gnu.org/licenses/>.

// Creator: César González
// Date: March 2017

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <stdint.h>

using namespace std;

const string gen_coverage_p = 
"# GNU Plot script to generate the coverage signal for 5mC and 5hmC\n \
# HPG-Mapper, 2017 César González\n \
\n \
set   autoscale                        # scale axes automatically\n \
unset log                              # remove any log-scaling\n \
unset label                            # remove any previous labels\n \
set xtic auto                          # set xtics automatically\n \
set ytic auto                          # set ytics automatically\n \
set terminal postscript eps enhanced color size 3.5,2.5 font 'Helvetica,10'\n \
set output \"result/coverage.eps\"\n \
\n \
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5   # --- blue\n \
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.5   # --- red\n \
\n \
set title \"Sequence coverage for chr%CHR%:%START%-%END%\"\n \
set xlabel \"Position\"\n \
set ylabel \"Sequencing coverage\"\n \
set xtics auto\n \
set format x \"%.0f\"\n \
set xrange [%XTICS_START%:%XTICS_END%] \n \
set key left top box 3\n \
\n \
plot \"result/coverage_mc.dlm\" using 1:2 title '5mC Cov'  with lines ls 1, \\\n \
     \"result/coverage_hmc.dlm\" using 1:2 title '5hmC Cov' with lines ls 2\n";

const string gen_meth_signal_p = 
"# GNU Plot script to generate the methylantion signal for 5mC and 5hmC\n \
# HPG-Mapper, 2017 César González\n \
\n \
set   autoscale                        # scale axes automatically\n \
unset log                              # remove any log-scaling\n \
unset label                            # remove any previous labels\n \
set xtic auto                          # set xtics automatically\n \
set ytic auto                          # set ytics automatically\n \
set terminal postscript eps enhanced color size 3.5,2.5 font 'Helvetica,10'\n \
set output \"result/methylation.eps\"\n \
\n \
set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 7 ps 1.5   # --- blue\n \
set style line 2 lc rgb '#dd181f' lt 1 lw 1 pt 5 ps 1.5   # --- red\n \
\n \
set title \"5mC and 5hmC signal for chr%CHR%:%START%-%END%\"\n \
set xlabel \"Position\"\n \
set ylabel \"Signal intensity\"\n \
set xtics auto\n \
set format x \"%.0f\"\n \
set xrange [%XTICS_START%:%XTICS_END%] \n \
set key left top box 3\n \
\n \
plot \"result/signal_mc.dlm\" using 1:2 title '5mC %'  with lines ls 1, \\\n \
     \"result/signal_hmc.dlm\" using 1:2 title '5hmC %' with lines ls 2\n";

typedef struct meth {
  uint64_t position;
  unsigned int c;
  unsigned int nc;
  unsigned int mc;
  unsigned int hmc;
} meth_t;

int main(int argc, char* argv[]) {
  // Get start and end positions from argv
  uint64_t start, end;
  unsigned int chromosome;
  double samples;

  if (argc != 6) {
    cout << "Usage: hpg-hmapper-graph-tool [file] [start_position] [end_position] [samples] [chr]" << endl;
    exit(1);
  }

  start = atoi(argv[2]);
  end = atoi(argv[3]);
  samples = atof(argv[4]);
  chromosome = atoi(argv[5]);

  // Read until the start position has been reached or surpassed, and
  // store the methylation data until the end position has been reached
  vector<meth_t> methData;
  ifstream file(argv[1], ifstream::in);

  string csv, temp;
  uint64_t position;
  bool completed = false;
  bool first = true;

  while (!file.eof() && !completed) {
    getline(file, csv);
    stringstream csvDelim(csv);

    // Check the position
    getline(csvDelim, temp, ' ');
    position = atoi(temp.c_str());

    if (position >= start && position <= end) {
      if (first) {
        cout << "Found first position: " << position << endl;
        first = false;
      }

      // Start storing data if its in the interval
      meth_t m;
      m.position = position;

      getline(csvDelim, temp, ' ');
      m.c = atoi(temp.c_str());

      getline(csvDelim, temp, ' ');
      m.nc = atoi(temp.c_str());
      
      getline(csvDelim, temp, ' ');
      m.mc = atoi(temp.c_str());

      getline(csvDelim, temp, '\n');
      m.hmc = atoi(temp.c_str());

      methData.push_back(m);
    } else if (position > end) {
      completed = true;
      cout << "Reached last position: " << position << endl;
    }
  }

  // Close the file and start processing
  file.close();

  vector<double> mcCoverage, hmcCoverage, mcSignal, hmcSignal, sigPositions;
  uint64_t steps = methData.size() / samples;
  uint64_t inputStep = (end - start) / samples;

  for (int i = 0; i < samples; ++i) {
    double tempCoverageMc = 0.0;
    double tempCoverageHmc = 0.0;

    double tempMcSignal = 0.0;
    double tempHmcSignal = 0.0;
    double maxMcSignal = 1.0;
    double maxHmcSignal = 1.0;

    for (int j = 0; j  < steps; ++j) {
      meth_t m = methData[j+i];
      tempCoverageMc += m.c + m.nc + m.mc; 
      tempCoverageHmc += m.c + m.nc + m.hmc;

      tempMcSignal += m.mc;
      tempHmcSignal += m.hmc;

      if (m.mc > maxMcSignal) {
        maxMcSignal = m.mc;
      }

      if (m.hmc > maxHmcSignal) {
        maxHmcSignal = m.hmc;
      }
    }

    tempCoverageMc /= (double)steps;
    tempCoverageHmc /= (double)steps;

    tempMcSignal /= (double)steps;
    tempHmcSignal /= (double)steps;

    tempMcSignal /= maxMcSignal;
    tempHmcSignal /= maxHmcSignal;

    mcCoverage.push_back(tempCoverageMc);
    hmcCoverage.push_back(tempCoverageHmc);
    mcSignal.push_back(tempMcSignal);
    hmcSignal.push_back(tempHmcSignal);
    sigPositions.push_back(i*inputStep + start);
  }

  // Print the data to files
  mkdir("./result/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  cout << endl ;
  cout << "Saving 5mC and 5hmC signal and coverage from " << start << " to " << end << endl;

  //----------------------------------

  char outTemp[256];

  cout << "File: coverage_mc.dlm" << endl;
  ofstream mcCov("result/coverage_mc.dlm", ofstream::out);

  for (int i = 0; i < mcCoverage.size(); ++i) {
    sprintf(outTemp, "%.0f %.3f\n", sigPositions[i], mcCoverage[i]);
    mcCov << outTemp;
  }

  mcCov.close();

  //----------------------------------

  cout << "File: coverage_hmc.dlm" << endl;
  ofstream hmcCov("result/coverage_hmc.dlm", ofstream::out);

  for (int i = 0; i < hmcCoverage.size(); ++i) {
    sprintf(outTemp, "%.0f %.3f\n", sigPositions[i], hmcCoverage[i]);
    hmcCov << outTemp;
  }

  hmcCov.close();

  //----------------------------------

  cout << "File: signal_mc.dlm" << endl;
  ofstream mcSig("result/signal_mc.dlm", ofstream::out);

  for (int i = 0; i < mcSignal.size(); ++i) {
    sprintf(outTemp, "%.0f %.3f\n", sigPositions[i], mcSignal[i]);
    mcSig << outTemp;
  }

  mcSig.close();

  //----------------------------------

  cout << "File: signal_hmc.dlm" << endl;
  ofstream hmcSig("result/signal_hmc.dlm", ofstream::out);

  for (int i = 0; i < hmcSignal.size(); ++i) {
    sprintf(outTemp, "%.0f %.3f\n", sigPositions[i], hmcSignal[i]);
    hmcSig << outTemp;
  }

  hmcSig.close();

  // Generate the plots using GNUplot
  char sedStr[512];
  memset(sedStr, 0, 256);
  int sysPid = -1;

  ofstream gplotFile("gen_coverage_temp.p", ofstream::out);
  gplotFile << gen_coverage_p;
  gplotFile.close();
  
  sprintf(sedStr, "sed -i -e 's/%%START%%/%i/g' gen_coverage_temp.p", start);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%END%%/%i/g' gen_coverage_temp.p", end);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%CHR%%/%i/g' gen_coverage_temp.p", chromosome);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%XTICS_START%%/%i/g' gen_coverage_temp.p", start);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%XTICS_END%%/%i/g' gen_coverage_temp.p", end);
  sysPid = system(sedStr);

  sysPid = system("gnuplot gen_coverage_temp.p");
  sysPid = system("rm gen_coverage_temp.p");

  //----------------------------------

  ofstream gPlotMeth("gen_meth_sig.p", ofstream::out);
  gPlotMeth << gen_meth_signal_p;
  gPlotMeth.close();
  
  sprintf(sedStr, "sed -i -e 's/%%START%%/%i/g' gen_meth_sig.p", start);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%END%%/%i/g' gen_meth_sig.p", end);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%CHR%%/%i/g' gen_meth_sig.p", chromosome);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%XTICS_START%%/%i/g' gen_meth_sig.p", start);
  sysPid = system(sedStr);

  sprintf(sedStr, "sed -i -e 's/%%XTICS_END%%/%i/g' gen_meth_sig.p", end);
  sysPid = system(sedStr);

  sysPid = system("gnuplot gen_meth_sig.p");
  sysPid = system("rm gen_meth_sig.p");

  // Open both graphs
  sysPid = system("xdg-open result/coverage.eps &");
  sysPid = system("xdg-open result/methylation.eps");
  sysPid = -1;
  wait(&sysPid);

  return 0;
}
